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
FFLAGS ?=  -fopenmp -fdefault-integer-8

# ---------------------------------
# List of sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/3d
SEISLIB = $(CLAW)/seismic/src/3d
#include $(AMRLIB)/Makefile.amr_3d

MODULES = \
  $(CLAW)/classic/src/utility_module.f90 \
  $(SEISLIB)/amr_module.f90 \
  $(SEISLIB)/regions_module.f90 \
  $(SEISLIB)/slices_module.f90 \
  $(AMRLIB)/gauges_module.f90

SOURCES = \
  ./advanc_debug.f \
  ./bc3amr.f \
  ./mapc2p.f90 \
  ./qinit.f90 \
  ./setaux.f90 \
  ./setaux1d.f90 \
  ./setprob.f90 \
  $(CLAW)/seismic/src/rp/rpn3_vc_elasticity_mapped_Yinvariant_fault.f90 \
  $(CLAW)/seismic/src/rp/rpt3_vc_elasticity_mapped_Yinvariant.f90 \
  $(CLAW)/seismic/src/rp/rptt3_vc_elasticity_mapped_Yinvariant.f90 \
  $(SEISLIB)/amr3.f90 \
  $(SEISLIB)/flag2refine.f90 \
  $(SEISLIB)/init_alloc.f90 \
  $(SEISLIB)/step3.f \
  $(SEISLIB)/stepgrid.f \
  $(SEISLIB)/tick.f \
  $(AMRLIB)/b4step3.f90 \
  $(AMRLIB)/opendatafile.f \
  $(AMRLIB)/qad.f \
  $(AMRLIB)/src3.f90 \
  $(AMRLIB)/src1d.f90 \
  $(AMRLIB)/bound.f \
  $(AMRLIB)/stepgrid_dimSplit.f \
  $(AMRLIB)/auxcoarsen.f \
  $(AMRLIB)/fixcapaq.f \
  $(AMRLIB)/estdt.f \
  $(AMRLIB)/igetsp.f \
  $(AMRLIB)/reclam.f \
  $(AMRLIB)/birect.f \
  $(AMRLIB)/cleanup.f \
  $(AMRLIB)/colate.f \
  $(AMRLIB)/errest.f \
  $(AMRLIB)/allowflag.f \
  $(AMRLIB)/bufnst.f \
  $(AMRLIB)/spest.f \
  $(AMRLIB)/errf1.f \
  $(AMRLIB)/gfixup.f \
  $(AMRLIB)/filval.f \
  $(AMRLIB)/filpatch.f \
  $(AMRLIB)/prefilp.f \
  $(AMRLIB)/flglvl.f \
  $(AMRLIB)/fluxad.f \
  $(AMRLIB)/fluxsv.f \
  $(AMRLIB)/ginit.f \
  $(AMRLIB)/grdfit.f \
  $(AMRLIB)/intfil.f \
  $(AMRLIB)/moment.f \
  $(AMRLIB)/nestck.f \
  $(AMRLIB)/prepc.f \
  $(AMRLIB)/prepf.f \
  $(AMRLIB)/projec.f \
  $(AMRLIB)/signs.f \
  $(AMRLIB)/findcut.f \
  $(AMRLIB)/smartbis.f \
  $(AMRLIB)/putnod.f \
  $(AMRLIB)/putsp.f \
  $(AMRLIB)/regrid.f \
  $(AMRLIB)/setgrd.f \
  $(AMRLIB)/setuse.f \
  $(AMRLIB)/stst1.f \
  $(AMRLIB)/trimbd.f \
  $(AMRLIB)/update.f \
  $(AMRLIB)/nodget.f \
  $(AMRLIB)/upbnd.f  \
  $(AMRLIB)/basic.f \
  $(AMRLIB)/outval.f \
  $(AMRLIB)/copysol.f \
  $(AMRLIB)/outvar.f \
  $(AMRLIB)/outmsh.f \
  $(AMRLIB)/outtre.f \
  $(AMRLIB)/domain.f  \
  $(AMRLIB)/setflags.f \
  $(AMRLIB)/shiftset.f \
  $(AMRLIB)/conck.f \
  $(AMRLIB)/domshrink.f \
  $(AMRLIB)/domprep.f \
  $(AMRLIB)/domup.f \
  $(AMRLIB)/domcopy.f \
  $(AMRLIB)/coarsen.f \
  $(AMRLIB)/intcopy.f \
  $(AMRLIB)/preintcopy.f \
  $(AMRLIB)/icall.f \
  $(AMRLIB)/preicall.f \
  $(AMRLIB)/step3x.f \
  $(AMRLIB)/step3y.f \
  $(AMRLIB)/step3z.f \
  $(AMRLIB)/flux3.f \
  $(AMRLIB)/flux3_dimSplit.f \
  $(AMRLIB)/cstore.f \
  $(AMRLIB)/saveqc.f \
  $(AMRLIB)/valout.f \
  $(AMRLIB)/check.f \
  $(AMRLIB)/restrt.f \
  $(AMRLIB)/quick_sort1.f \
  $(AMRLIB)/resize_alloc.f90 \
  $(AMRLIB)/restrt_alloc.f90 \
  $(AMRLIB)/inlinelimiter.f





#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

