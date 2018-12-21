CC=gcc
OMP_PARALLEL=TRUE
LIB_DIR=/cm/shared/apps/gromacs/4.6.7/lib
INC_DIRS=/cm/shared/apps/gromacs/4.6.7/include/gromacs/


ifeq ($(OMP_PARALLEL), TRUE) 
CFLAGS=-O3 -Wall -fopenmp -DOMP_PARALLEL=TRUE
LIBFLAGS=-lgmx -ldl -lm -pthread -lgomp
else
CFLAGS=-O3 -Wall
LIBFLAGS=-lgmx -ldl -lm 
endif

OBJECTS=AmideImaps.o

all: $(OBJECTS)
	$(CC) $(OBJECTS) -L $(LIB_DIR) $(LIBFLAGS) -o AmideImaps

AmideImaps.o: AmideImaps_backbone_sidechain.h AmideImaps.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) AmideImaps.c

clean:
	rm *.o
