TARGET = smat

CC = mpicc

CFLAGS = -std=gnu99 -Wall -O3
LFLAGS = -lm -L/opt/local/lapack-3.5.0/lib -lblas -llapack 
LFLAGS += -L/opt/local/scalapack-2.0.2/lib -lscalapack -lconfig -lmpi
LFLAGS += -L/opt/local/macports-latest/lib/gcc48 -lgfortran

FEATURES =
# FEATURES += -DSINGLE_PRECISION

CFLAGS += $(FEATURES)

SRC = smat.c

all:			$(TARGET)

$(TARGET):		$(SRC:.c=.o)
			$(CC) -o $(TARGET) $(SRC:.c=.o) $(LFLAGS)

clean:
			@rm -f $(SRC:.c=.o) $(SRC:.c=.d) *~
