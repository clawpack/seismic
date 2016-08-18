"""
Module to set up run time parameters for Clawpack -- classic code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np

#------------------------------
def setrun(claw_pkg='amrclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "amrclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data


    assert claw_pkg.lower() == 'amrclaw',  "Expected claw_pkg = 'amrclaw'"

    num_dim = 3
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    # Sample setup to write one line to setprob.data ...
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('domain_depth', 200e3, 'depth of domain')
    probdata.add_param('domain_width', 400e3, 'width of domain (in x direction)')
    probdata.add_param('domain_length', 400e3, 'length of domain (in y direction)')
    probdata.add_param('fault_xcenter', 25e3, 'centeroid of fault (x)')
    probdata.add_param('fault_ycenter', 0.0, 'centeroid of fault (y)')
    probdata.add_param('fault_width', 50735, 'width of fault (in x direction)')
    probdata.add_param('fault_length', 50735, 'width of fault (in y direction)')
    probdata.add_param('fault_dip', 0.17, 'angle of fault dip')
    probdata.add_param('fault_depth', 19.3e3, 'depth of fault centroid')

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # ---------------
    # Spatial domain
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Number of grid cells
    num_cells_fault_width = 10
    num_cells_fault_length = 10
    dx = probdata.fault_width/num_cells_fault_width
    dy = probdata.fault_length/num_cells_fault_length
    ## specify dz using dx,dy
    num_cells_above = np.rint(probdata.fault_depth/min(dx,dy))
    dz = probdata.fault_depth/num_cells_above
    clawdata.num_cells[0] = int(np.ceil(probdata.domain_width/dx))
    clawdata.num_cells[1] = int(np.ceil(probdata.domain_length/dy))
    clawdata.num_cells[2] = int(np.ceil(probdata.domain_depth/dz))

    # Lower and upper edgse of computational domain:
    ## note the size of domain is likely expanded here
    num_cells_remain = clawdata.num_cells[0] - num_cells_fault_width
    clawdata.lower[0] = probdata.fault_xcenter-0.5*probdata.fault_width - np.floor(num_cells_remain/2.0)*dx
    clawdata.upper[0] = probdata.fault_xcenter+0.5*probdata.fault_width + np.floor(num_cells_remain/2.0)*dx
    num_cells_remain = clawdata.num_cells[1] - num_cells_fault_length
    clawdata.lower[1] = probdata.fault_ycenter-0.5*probdata.fault_length - np.floor(num_cells_remain/2.0)*dy
    clawdata.upper[1] = probdata.fault_ycenter+0.5*probdata.fault_length + np.floor(num_cells_remain/2.0)*dy
    clawdata.lower[2] = -clawdata.num_cells[2]*dz
    clawdata.upper[2] = 0.0
    probdata.domain_width = clawdata.upper[0] - clawdata.lower[0]
    probdata.domain_length = clawdata.upper[1] - clawdata.lower[1]
    probdata.domain_depth = clawdata.upper[2] - clawdata.lower[2]

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 9

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 12

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 12


    # -------------
    # Initial time:
    # -------------

    # should be set to final time in reload file if reloading
    clawdata.t0 = 0.000000

    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data


    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 50 #80
        clawdata.tfinal = 100.0 #6e-5 #0.0001800000 #0.00003
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  [48.0, 65.0, 90.0, 135.0, 180.0]

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 10
        clawdata.output_t0 = True  # output at initial (or restart) time?


    clawdata.output_format = 'ascii'      # 'ascii' or 'binary'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    clawdata.verbosity = 1



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==False: fixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 0.25

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1.000000e+99

    # Desired Courant number if variable dt used
    clawdata.cfl_desired = 0.900000
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.000000

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 1000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Use dimensional splitting?
    clawdata.dimensional_split = False

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 22


    # Number of waves in the Riemann solution:
    clawdata.num_waves = 6

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = ['vanleer', 'vanleer', 'vanleer', 'vanleer', 'vanleer','vanleer']

    clawdata.use_fwaves = False    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 0  # use axi-symmetric source initially


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'extrap'   # at yupper

    clawdata.bc_lower[2] = 'extrap'   # at zlower
    clawdata.bc_upper[2] = 'user'   # at zupper

    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.
      clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5


    # ---------------
    # Gauges:
    # ---------------
    gauges = rundata.gaugedata.gauges

    # Gauges across transducer face
    #top edgs
    xgauges = np.linspace(clawdata.lower[0]+1, clawdata.upper[0]-1,100)
    ygauges = np.linspace(clawdata.lower[1]+1, clawdata.upper[1]-1,100)
    gcount = 0
    for x in xgauges:
        for y in ygauges:
            gauges.append([gcount, x, y, clawdata.upper[2]-1, 0.0, 1e9])
            gcount = gcount + 1
    # above fault place
    xgauges = np.linspace(probdata.fault_xcenter-0.5*probdata.fault_width,probdata.fault_xcenter+0.5*probdata.fault_width,100)
    ygauges = np.linspace(probdata.fault_ycenter-0.5*probdata.fault_length,probdata.fault_ycenter+0.5*probdata.fault_length,100)
    for x in xgauges:
        for y in ygauges:
            gauges.append([gcount, x, y, -probdata.fault_depth+1, 0.0, 1e9])
            gcount = gcount + 1


    # -----------------
    # Slice parameters:
    # -----------------

#    rundata.slicedata.slices_xy = [0.0]
#    rundata.slicedata.slices_xz = [0.0]

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 4

    # List of refinement ratios at each level (length at least amr_level_max-1)
    amrdata.refinement_ratios_x = [4,4,2]
    amrdata.refinement_ratios_y = [4,4,2]
    amrdata.refinement_ratios_z = [4,4,2]
    amrdata.refinement_ratios_t = [4,4,2]

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one
    # of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center', 'center', 'center', 'center', 'xleft', 'xleft', 'xleft', 'yleft', 'yleft', 'yleft', 'zleft','capacity']

    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.000000e+00  # Richardson to
    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 1.0e-4  # tolerance used in this routine
    # Flags regions to be refined where tol is roughly the smallest amplitude of normal velocity to capture

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 2

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells
    # refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged
    # cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 4


    # -------------------
    # Refinement Regions:
    # -------------------
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2,z1,z2]
    #regions.append([1,1, 0,1e9, -1.e9, 1e9, -1e9, 1e9, -1e9, 0.])
    #regions.append([1,2, 0,1e9, -90e3, 140e3, -90e3, 140e3, -150e3, 0.])
    #regions.append([1,3, 0,1e9, -75e3, 125e3, -75e3, 125e3, -80e3, 0.])
    #regions.append([1,4, 0,1e9, -50e3, 105e3, -50e3, 105e3, -70e3, 0.])
    #regions.append([4,4, 0,1, -10.e3, 60e3, 0.85, 0.95])
    #regions.append([amrdata.amr_levels_max,amrdata.amr_levels_max,
    #                0,1, probdata.fault_xcenter-0.5*probdata.fault_width-2*dx,
    #                probdata.fault_xcenter+0.5*probdata.fault_width+2*dx,
    #                probdata.fault_ycenter-0.5*probdata.fault_length-2*dy,
    #                probdata.fault_ycenter+0.5*probdata.fault_length+2*dy,
    #                -probdata.fault_depth-2*dz, -probdata.fault_depth+2*dz])



    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting


    return rundata

    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
