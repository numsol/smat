/* ----------------------------------------------------------------------
   SMAT - Copyright (2014) Anthony B. Costa.
   anthony.costa@numericalsolutions.org, Numerical Solutions, Inc.

   A simple C program with a configuration file-driven interface 
   supporting matrix multiplication and QR decomposition implemented in 
   Scalapack for the purposes of "burning in" an arbitrary number of 
   cluster nodes while testing performance.
------------------------------------------------------------------------- */

#ifndef SMAT_H
#define SMAT_H

#define VERSION "0.0.2-RLS"
#define MAX_LEN 128

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/* mpi standard variables */
int RANK, ROOT, NUM_TASKS;

/* blacs standard variables */
int DESC[9];
int INFO, MONE, ZERO, ONE;
int CONTEXT, BLOCK_SIZE;
int FABRIC_ROW, FABRIC_COL;
int MAX_ROW, MAX_COL;
int LOCAL_ROW, LOCAL_COL;
int LOCAL_SIZE, TOTAL_SIZE;
int NPROW, NPCOL, MYPROW, MYPCOL;

/* run parameters */
int VERBOSITY;
char TEST_TYPE[MAX_LEN];
char INPUT_FILE[MAX_LEN];
int NUM_LOOPS, SIZE_DIM;

/* output string concatenator */
char OUTPUT_STRING[MAX_LEN];

/* select precision */
#ifdef SINGLE_PRECISION
  typedef float decimal;
  #define MPI_DECIMAL MPI_FLOAT
#else
  typedef double decimal;
  #define MPI_DECIMAL MPI_DOUBLE
#endif /* SINGLE_PRECISION */

/* helper functions */
void run_multiply (void);
void run_qr (void);
decimal random_return (void);
void get_input (int, char **);
void input_get_int (char *, int *, int);
void input_get_str (char *, char *, char *);
int string_match (char *, char *);
void output_barrier (void);
void output_print_header (void);
void output_print_bar (void);
void output_start_section (char *);
void output_simple_line (char *);
void output_message (char *, char *);
void output_timing (int, double, char *);
void output_finalize (void);
void exit_failure (char *);

/* external blacs and scalapack functions from fortran interface */
extern void blacs_pinfo_ (int *, int *);
extern void blacs_get_ (int *, int *, int *);
extern void blacs_gridinit_ (int *, char *, int *, int *);
extern void blacs_gridinfo_ (int *, int *, int *, int *, int *);
extern void blacs_exit_ (int *);
extern int numroc_ (int *, int *, int *, int *, int *);
extern void descinit_ (int *, int *, int *, int *, int *, \
  int *, int *, int *, int *, int *);

void (*matrix_value_set)(decimal *, int *, int *, int *, decimal *);
void (*matrix_matrix_multiply)(char *, char *, int *, int *, int *, decimal *, \
  decimal *, int *, int *, int *, decimal *, int *, int *, int *, \
  decimal *, decimal *, int *, int *, int *);
void (*qr_pack)(int *, int *, decimal *, int *, int *, int *, \
  decimal *, decimal *, int *, int *);
void (*qr_unpack)(int *, int *, int *, decimal *, int *, int *, \
  int *, decimal *, decimal *, int *, int *);

#ifdef SINGLE_PRECISION
  extern void pselset_ (decimal *, int *, int *, int *, decimal *);
  extern void psgemm_ (char *, char *, int *, int *, int *, decimal *, \
    decimal *, int *, int *, int *, decimal *, int *, int *, int *, \
    decimal *, decimal *, int *, int *, int *);
  extern void psgeqrf_ (int *, int *, decimal *, int *, int *, int *, \
    decimal *, decimal *, int *, int *);
  extern void psorgqr_ (int *, int *, int *, decimal *, int *, int *, \
    int *, decimal *, decimal *, int *, int *);
#else
  extern void pdelset_ (decimal *, int *, int *, int *, decimal *);
  extern void pdgemm_ (char *, char *, int *, int *, int *, decimal *, \
    decimal *, int *, int *, int *, decimal *, int *, int *, int *, \
    decimal *, decimal *, int *, int *, int *);
  extern void pdgeqrf_ (int *, int *, decimal *, int *, int *, int *, \
    decimal *, decimal *, int *, int *);
  extern void pdorgqr_ (int *, int *, int *, decimal *, int *, int *, \
    int *, decimal *, decimal *, int *, int *);
#endif /* SINGLE_PRECISION */

#endif /* SMAT_H */
