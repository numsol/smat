#include "smat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <libconfig.h>
#include <mpi.h>

/* getopt definitions */
static const struct option long_opts [] =
{
  {"input_file", optional_argument, NULL, 'i'},
  {NULL, no_argument, NULL, 0}
};

/* libconfig structure */
struct config_t input_cfg;
config_setting_t *setting_cfg;

int
main (int argc, char **argv)
{
  /* initialize mpi */
  int MPI_STARTUP = MPI_Init (&argc, &argv);
  if (MPI_STARTUP != MPI_SUCCESS)
    MPI_Abort (MPI_COMM_WORLD, MPI_STARTUP);

  MPI_Comm_size (MPI_COMM_WORLD, &NUM_TASKS);
  MPI_Comm_rank (MPI_COMM_WORLD, &RANK);

  /* print header */
  output_print_header ();

  /* get input from flags and input file */
  get_input (argc, argv);

  /* can't have more fabric blocks than available mpi tasks */
  if (FABRIC_ROW * FABRIC_COL > NUM_TASKS)
    exit_failure ("CAN NOT INITIALIZE MORE FABRIC BLOCKS THAN NUMBER OF TASKS");

  /* warn user if the number of fabric blocks is less than the number of mpi tasks */
  if (FABRIC_ROW * FABRIC_COL < NUM_TASKS)
    exit_failure ("SOME MPI TASKS ARE NOT BEING USED IN BLACS FABRIC");

  /* initialize blacs fabric */
  MONE = -1; ZERO = 0; ONE = 1;
  blacs_pinfo_ (&RANK, &NUM_TASKS);
  blacs_get_ (&MONE, &ZERO, &CONTEXT);
  blacs_gridinit_ (&CONTEXT, "Row", &FABRIC_ROW, &FABRIC_COL);
  blacs_gridinfo_ (&CONTEXT, &NPROW, &NPCOL, &MYPROW, &MYPCOL);

  LOCAL_ROW = numroc_ (&SIZE_DIM, &BLOCK_SIZE, &MYPROW, &ZERO, &NPROW);
  LOCAL_COL = numroc_ (&SIZE_DIM, &BLOCK_SIZE, &MYPCOL, &ZERO, &NPCOL);

  MAX_ROW = MAX (1, LOCAL_ROW);
  MAX_COL = MAX (1, LOCAL_COL);

  /* local and total matrix sizes */
  LOCAL_SIZE = MAX_ROW * MAX_COL;
  TOTAL_SIZE = SIZE_DIM * SIZE_DIM;

  /* blacs descriptor */
  descinit_ (DESC, &SIZE_DIM, &SIZE_DIM, &BLOCK_SIZE, &BLOCK_SIZE, \
    &ZERO, &ZERO, &CONTEXT, &MAX_ROW, &INFO);

  /* assign scalapack function pointers */
  #ifdef SINGLE_PRECISION
    matrix_value_set = pselset_;
    matrix_matrix_multiply = psgemm_;
    qr_pack = psgeqrf_;
    qr_unpack = psorgqr_;
  #else
    matrix_value_set = pdelset_;
    matrix_matrix_multiply = pdgemm_;
    qr_pack = pdgeqrf_;
    qr_unpack = pdorgqr_;
  #endif /* SINGLE_PRECISION */

  /* start up the random number subsystem with a different
     seed on each processor */
  srand ((unsigned long int) (time (NULL) + RANK));

  /* determine the run type and execute */
  if (string_match (TEST_TYPE, "multiply"))
    run_multiply ();
  else if (string_match (TEST_TYPE, "qr"))
    run_qr ();
  else
    exit_failure ("DO NOT UNDERSTAND INPUT OPTION TEST TYPE");

  output_finalize ();
  blacs_exit_ (&ONE);
  MPI_Finalize ();
  return (EXIT_SUCCESS);
}

void
run_multiply (void)
{
  output_start_section ("MATRIX MULTIPLY TEST");

  /* counters */
  int i, timing;
  double counter_a, counter_b;

  /* allocate memory */
  output_simple_line ("ALLOCATING AND FILLING MEMORY");
  counter_a = MPI_Wtime ();

  decimal *MATRIX_A;
  decimal *MATRIX_B;
  decimal *MATRIX_C;

  MATRIX_A = malloc (LOCAL_SIZE * sizeof (decimal));
  if (MATRIX_A == NULL)
    exit_failure ("RAN OUT OF MEMORY");

  MATRIX_B = malloc (LOCAL_SIZE * sizeof (decimal));
  if (MATRIX_B == NULL)
    exit_failure ("RAN OUT OF MEMORY");

  MATRIX_C = malloc (LOCAL_SIZE * sizeof (decimal));
  if (MATRIX_C == NULL)
    exit_failure ("RAN OUT OF MEMORY");

  /* fill the matricies on each process with random data */
  for (i = 0; i < LOCAL_SIZE; i++)
  {
    MATRIX_A[i] = random_return ();
    MATRIX_B[i] = random_return ();
  }
  
  counter_b = MPI_Wtime ();
  if (VERBOSITY && NUM_TASKS > 1)
    for (timing = 0; timing < NUM_TASKS; timing++)
      output_timing (timing, counter_b - counter_a, "MEMORY ALLOCATION AND FILL - ALL TASKS");
  else
    output_timing (ROOT, counter_b - counter_a, "MEMORY ALLOCATION AND FILL - ROOT TASK");

  /* mpi barrier to wait for output */
  output_barrier ();

  /* scaling parameters */
  decimal alpha = 1.0;
  decimal beta = 0.0;

  /* loop and perform matrix multiplication */
  for (i = 0; i < NUM_LOOPS; i++)
  {
    if (VERBOSITY && NUM_TASKS > 1)
    {
      sprintf (OUTPUT_STRING, "%s%d", "STARTING ITERATION ", i + 1);
      output_simple_line (OUTPUT_STRING);
    }

    counter_a = MPI_Wtime ();

    matrix_matrix_multiply ("N", "N", &SIZE_DIM, &SIZE_DIM, &SIZE_DIM, \
      &alpha, MATRIX_A, &ONE, &ONE, DESC, MATRIX_B, &ONE, &ONE, DESC, &beta, \
      MATRIX_C, &ONE, &ONE, DESC);

    counter_b = MPI_Wtime ();
    if (VERBOSITY && NUM_TASKS > 1)
      for (timing = 0; timing < NUM_TASKS; timing++)
        output_timing (timing, counter_b - counter_a, "ITERATION - ALL TASKS");
    else
      output_timing (ROOT, counter_b - counter_a, "ITERATION - ROOT TASK");

    /* mpi barrier to wait for output */
    output_barrier ();
  }

  /* free memory */
  free (MATRIX_A);
  free (MATRIX_B);
  free (MATRIX_C);
}

void
run_qr (void)
{
  output_start_section ("QR PACK/UNPACK TEST");

  /* counters */
  int i, timing;
  double counter_a, counter_b;

  /* allocate memory */
  output_simple_line ("ALLOCATING AND FILLING MEMORY");
  counter_a = MPI_Wtime ();

  decimal *MATRIX_A;
  decimal *TAU, *WORK;

  MATRIX_A = malloc (LOCAL_SIZE * sizeof (decimal));
  if (MATRIX_A == NULL)
    exit_failure ("RAN OUT OF MEMORY");

  WORK = malloc (LOCAL_SIZE * sizeof (decimal));
  if (WORK == NULL)
    exit_failure ("RAN OUT OF MEMORY");

  TAU = malloc (SIZE_DIM * sizeof (decimal));
  if (TAU == NULL)
    exit_failure ("RAN OUT OF MEMORY");

  /* fill the matrix on each process with random data */
  for (i = 0; i < LOCAL_SIZE; i++)
    MATRIX_A[i] = random_return ();

  counter_b = MPI_Wtime ();
  if (VERBOSITY && NUM_TASKS > 1)
    for (timing = 0; timing < NUM_TASKS; timing++)
      output_timing (timing, counter_b - counter_a, "MEMORY ALLOCATION AND FILL - ALL TASKS");
  else
    output_timing (ROOT, counter_b - counter_a, "MEMORY ALLOCATION AND FILL - ROOT TASK");

  /* mpi barrier to wait for output */
  output_barrier ();

  /* loop and perform matrix multiplication */
  for (i = 0; i < NUM_LOOPS; i++)
  {
    if (VERBOSITY && NUM_TASKS > 1)
    {
      sprintf (OUTPUT_STRING, "%s%d", "STARTING ITERATION ", i + 1);
      output_simple_line (OUTPUT_STRING);
    }

    counter_a = MPI_Wtime ();

    qr_pack (&SIZE_DIM, &SIZE_DIM, MATRIX_A, &ONE, &ONE, DESC, TAU, \
      WORK, &LOCAL_SIZE, &INFO);

    counter_b = MPI_Wtime ();
    if (VERBOSITY && NUM_TASKS > 1)
      for (timing = 0; timing < NUM_TASKS; timing++)
        output_timing (timing, counter_b - counter_a, "QR PACK ITERATION - ALL TASKS");
    else
      output_timing (ROOT, counter_b - counter_a, "QR PACK ITERATION - ROOT TASK");

    /* mpi barrier to wait for output */
    output_barrier ();

    counter_a = MPI_Wtime ();

    qr_unpack (&SIZE_DIM, &SIZE_DIM, &SIZE_DIM, MATRIX_A, &ONE, &ONE, \
      DESC, TAU, WORK, &LOCAL_SIZE, &INFO);

    counter_b = MPI_Wtime ();
    if (VERBOSITY && NUM_TASKS > 1)
      for (timing = 0; timing < NUM_TASKS; timing++)
        output_timing (timing, counter_b - counter_a, "QR UNPACK ITERATION - ALL TASKS");
    else
      output_timing (ROOT, counter_b - counter_a, "QR UNPACK ITERATION - ROOT TASK");

    /* mpi barrier to wait for output */
    output_barrier ();
  }

  /* free memory */
  free (MATRIX_A);
  free (TAU);
  free (WORK);
}

decimal
random_return (void)
{
  decimal value;

  value = rand () / (RAND_MAX + 1.0);

  return value;
}

void
get_input (int argc, char **argv)
{
  /* default zero verbosity and input file at input.dat */
  VERBOSITY = 0;
  strcpy (INPUT_FILE, "input.dat");

  /* parse the command line flags */
  int opt;
  char *opt_string = "i::v";
  while (1)
  {
    int long_index = 0;
    opt = getopt_long (argc, argv, opt_string, long_opts, &long_index);
    if (opt == -1) break;

    switch (opt)
    {
      case 'i':
        if (optarg)
        {
          strcpy (INPUT_FILE, optarg);
          break;
        }
        break;
      case 'v':
        VERBOSITY++;
        break;
      default:
        break;
    }
  }

  /* read and verify format of input file */
  config_init (&input_cfg);

  if (!config_read_file (&input_cfg, INPUT_FILE))
    exit_failure ("INPUT FILE FORMAT IS BAD");

  output_start_section ("INPUT OPTIONS");

  /* read data from input file */
  input_get_str ("test_type", TEST_TYPE, "multiply");
  input_get_int ("block_size", &BLOCK_SIZE, 16);
  input_get_int ("fabric_row", &FABRIC_ROW, 1);
  input_get_int ("fabric_col", &FABRIC_COL, 1);
  input_get_int ("num_loops", &NUM_LOOPS, 10);
  input_get_int ("size_dim", &SIZE_DIM, 10);
}

void
input_get_int (char *id, int *var, int def)
{
  setting_cfg = config_lookup (&input_cfg, id);
  if (setting_cfg)
    *var = config_setting_get_int (setting_cfg);
  else
    *var = def;

  /* variable must be positive */
  if (*var < 1)
    exit_failure ("ALL INPUT VARIABLES MUST BE POSITIVE INTEGERS");

  /* tell the user about what we've read from file */
  if (RANK == ROOT)
    fprintf (stdout, "## %-46s: %-31d\n", id, *var);
}

void
input_get_str (char *id, char *var, char *def)
{
  setting_cfg = config_lookup (&input_cfg, id);
  if (setting_cfg)
    strcpy (var, config_setting_get_string (setting_cfg));
  else
    strcpy (var, def);

  if (RANK == ROOT)
    fprintf (stdout, "## %-46s: %-31s\n", id, var);
}

int
string_match (char *stringa, char *stringb)
{
  if (strcmp (stringa, stringb))
    return (0);
  else
    return (1);
}

void
output_barrier (void)
{
  /* hack, mpi doesn't respect scalapack very well */
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
}

void
output_print_header (void)
{
  if (RANK == ROOT)
  {
    /* decide whether this is single or double precision binary */
    #ifdef SINGLE_PRECISION
      char *precision = "Single";
    #else
      char *precision = "Double";
    #endif

    /* date and time */
    const char *date = __DATE__;
    const char *time = __TIME__;

    fprintf (stdout, "\n");
    output_print_bar ();
    fprintf (stdout, "##                      SCALAPCK MPI Applications Tester                      ##\n");
    fprintf (stdout, "##                     Copyright 2011 by Anthony B. Costa                     ##\n");
    fprintf (stdout, "##               Version %s Compiled in %s Precision               ##\n", VERSION, precision);
    output_print_bar ();
    output_start_section ("BOOTING UP");
    fprintf (stdout, "## %-46s: %-31d\n", "PROCESSOR COUNT", NUM_TASKS);
    fprintf (stdout, "## %-46s: %-31s\n", "BUILT ON DATE", date);
    fprintf (stdout, "## %-46s: %-31s\n", "BUILT AT TIME", time);
  }
}

void
output_print_bar (void)
{
  if (RANK == ROOT)
    fprintf (stdout, "################################################################################\n");
}

void
output_start_section (char *string)
{
  if (RANK == ROOT)
  {
    fprintf (stdout, "\n");
    fprintf (stdout, "%80s\n", string);
    output_print_bar ();
  }
}

void
output_simple_line (char *message)
{
  if (RANK == ROOT)
    fprintf (stdout, "## %-77s\n", message);
}

void
output_message (char *name, char *message)
{
  if (RANK == ROOT)
  {
    output_start_section (name);
    fprintf (stdout, "## %-77s\n", message);
  }
}

void
output_finalize (void)
{
  output_message ("COMPLETE", "SMAT EXIT");
  if (RANK == ROOT)
    fprintf (stdout, "\n");
}

void
output_timing (int proc, double var, char *name)
{
  if (RANK == proc)
  {
    if (var < 60)
      fprintf (stdout, "## %-46s: %-21.4f SECONDS\n", name, var);
    else if (var < 3600)
      fprintf (stdout, "## %-46s: %-21.4f MINUTES\n", name, var / 60);
    else if (var < 86400)
      fprintf (stdout, "## %-46s: %-23.4f HOURS\n", name, var / 3600);
    else if (var < 604800)
      fprintf (stdout, "## %-46s: %-24.4f DAYS\n", name, var / 86400);
    else
      fprintf (stdout, "## %-46s: %-23.4f WEEKS\n", name, var / 604800);
  }
}

void
exit_failure (char *message)
{
  if (RANK == ROOT)
  {
    output_message ("EXITING ON ERROR", message);
    fprintf (stdout, "\n");
  }

  MPI_Finalize ();
  exit (EXIT_FAILURE);
}
