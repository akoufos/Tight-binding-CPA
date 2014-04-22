# Makefile for TB-CPA program
#--------------------------------------------------------------------------
# Defaults
#--------------------------------------------------------------------------
MAKE = make
G90 = gfortran
G90_OPTS = -ffree-form -O3 -ftree-vectorize
F90 = ifort
F90_OPTS = -O3
OMP = -openmp
GOMP = -fopenmp
MKL = -I$(MKLROOT)/include/intel64/lp64 -mkl=parallel
FFLAGS = $(OMP) $(MKL) $(F90_OPTS)

SRC_modules = global.f90

SRC_cpa = bisect.f90 calcSig.f90 ccrlu.f90 clubk.f90 cpaDOS.f90 \
	cpagen.f90 falsi.f90 fixpt.f90 kpts.f90 greens.f90 mainn.f90 \
	newton.f90 readin.f90 readSec.f90 setInitHam.f90 setHam.f90 \
	setOnsites.f90 simp.f90
#greens_omp.f90
#	kpts.f90 mainn.f90 readin.f90 readSec.f90 setHam.f90 \

SRC_omp = cmplxInv.f90 #greens_omp.f90

SRC = $(SRC_modules) $(SRC_cpa) $(SRC_omp)

OBJ_omp = $(SRC_omp:.f90=.o)
$(OBJ_omp): F90_OPTS+=$(OMP)
#$(OBJ_omp): G90_OPTS+=$(GOMP)
OBJ = $(SRC:.f90=.o)
EXE = cpaFeSe

#--------------------------------------------------------------------------
# Suffix rules
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -warn all -c $<
#	$(G90) $(G90_OPTS) -c $<


#--------------------------------------------------------------------------
# Targets
#--------------------------------------------------------------------------
cpaTB:	$(OBJ)
#	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ)
	$(F90) $(F90_OPTS) $(OMP) -o $(EXE) $(OBJ)
#	$(G90) $(G90_OPTS) $(GOMP) -o $(EXE) $(OBJ)

all:	cpaTB

clean:
	rm -f *.o *.mod *__genmod.f90 $(EXE)
