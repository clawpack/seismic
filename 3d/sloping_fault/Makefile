# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran
GFORTRAN_UNBUFFERED_PRECONNECTED=y

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?=  -fopenmp

# ---------------------------------
# List of sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/3d
SEISLIB = $(CLAW)/seismic/src/3d
include $(AMRLIB)/Makefile.amr_3d

MODULES = \
  $(SEISLIB)/regions_module.f90 \
  $(AMRLIB)/slices_module.f90 

SOURCES = \
  ./b4step3.f90 \
  ./bc3amr.f \
  ./mapc2p.f90 \
  ./qinit.f90 \
  ./setaux.f90 \
  ./setprob.f90 \
  ../rp/rpn3_vc_elasticity_mapped_Yinvariant_fault.f90 \
  ../rp/rpt3_vc_elasticity_mapped_Yinvariant.f90 \
  ../rp/rptt3_vc_elasticity_mapped_Yinvariant.f90 \
  $(SEISLIB)/amr3.f90 \
  $(SEISLIB)/flag2refine.f90 \


#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

