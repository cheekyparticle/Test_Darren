USR := $(shell whoami)

#NESTLIBDIR = ../
ifndef NESTLIBDIR
  NESTLIBDIR = ../
endif

LAPACKLIB = -llapack
ifneq (, $(findstring andrewcheek, $(USR)))
	GSLLIB = /sw/lib
else 
	GSLLIB =/usr/local/lib
endif

LIBS =  -L$(NESTLIBDIR) -L$(GSLLIB) -lnest3 $(LAPACKLIB)
#LIBS =  -L$(NESTLIBDIR) -lmultinest_mpi -L$(GSLLIB) $(LAPACKLIB) 
 

EXPS = ../../DirectDetectionCalc/XENON/*.c ../../DirectDetectionCalc/GERMANIUM/*.c ../../DirectDetectionCalc/Argon/*.c ../../DirectDetectionCalc/FLUORINE/*.c 
OBJFILES = *.o
FC = gfortran -fopenmp 
CC = gcc -fopenmp
FFLAGS += -w -O3 -ffree-line-length-none -fopenmp
CFLAGS += -O3 -fopenmp

SOURCES= eggbox.c ../../DirectDetectionCalc/hcubature.c ../../DirectDetectionCalc/source/*.c

all: eggbox 

%.o: %.c
	$(CC) $(CFLAGS) -c $(SOURCES) $(EXPS)
 
eggbox: $(OBJFILES)
	$(FC) $(FFLAGS) -o ../eggboxC $(OBJFILES) $(LIBS) -lgsl -lgslcblas -DHCUBATURE  -lm 

clean:
	rm -f *.o *.mod ../eggboxC
