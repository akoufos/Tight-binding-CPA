# Makefile for TB-CPA program
#--------------------------------------------------------------------------
# Defaults
#--------------------------------------------------------------------------
MAKE = make
F90 = ifort
F90_OPTS = -O3
OMP = -openmp
MKL = -I$(MKLROOT)/include/intel64/lp64 -mkl=parallel
FFLAGS = $(OMP) $(MKL) $(F90_OPTS)

SRC_modules = global.f90

SRC_cpa = ccrlu.f90 clubk.f90 cmplxInv.f90 cmplxLin.f90 cpaDOS.f90 \
	cpagen.f90 crete.f90 kpts.f90 readin.f90 readSec.f90 simp.f90 \
	mainn.f90

SRC = $(SRC_modules) $(SRC_cpa)

OBJ = $(SRC:.f90=.o)
EXE = cpaFeSe

#--------------------------------------------------------------------------
# Suffix rules
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#--------------------------------------------------------------------------
# Targets
#--------------------------------------------------------------------------
cpaTB:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ)

all:	cpaTB

clean:
	rm -f *.o *.mod $(EXE)
