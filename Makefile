#CAMB Makefile for macOS with gfortran and GSL

# GSL installation path (adjust if needed for your system)
# On macOS with Homebrew: brew install gsl
GSL_DIR = /opt/homebrew/Cellar/gsl/2.8

#Set FISHER=Y to compile bispectrum fisher matrix code
FISHER=Y

#Gfortran compiler configuration for macOS
F90C     = gfortran
FFLAGS = -O3 -fopenmp -ffast-math -march=native -funroll-loops -fPIC

# GSL flags
FFLAGS += -I$(GSL_DIR)/include
LIBS = -L$(GSL_DIR)/lib -lgsl -lgslcblas

# CFITSIO library
CFITSIO_DIR = /opt/homebrew/Cellar/cfitsio/4.6.2
FFLAGS += -I$(CFITSIO_DIR)/include
LIBS += -L$(CFITSIO_DIR)/lib -lcfitsio

# LAPACK/BLAS libraries (using Accelerate framework on macOS)
ifneq ($(FISHER),)
FFLAGS += -framework Accelerate
LIBS += -framework Accelerate
endif

#Settings for building camb_fits
#Location of FITSIO and name of library
FITSDIR       = /opt/homebrew/lib
FITSLIB       = cfitsio
#Location of HEALPIX for building camb_fits
HEALPIXDIR    = /opt/homebrew/healpix

ifneq ($(FISHER),)
FFLAGS += -DFISHER
EXTCAMBFILES = Matrix_utils.o
else
EXTCAMBFILES =
endif

include ./Makefile_main
