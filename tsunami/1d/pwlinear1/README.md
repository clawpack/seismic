
Test problem with piecewise linear topography profile consisting of:
flat ocean, continental slope, flat or sloping shelf, linear beach.

First adjust parameters in `make_topo_and_grid.py` and run to create grid.data
file specifying cell edges and topo values at these edges.

Note that currently `mx` set in `make_topo_and_grid.py` should agree with
`num_cells[0]` in setrun.py`.

Then run geoclaw code, e.g. via `make .plots`.  
This code depends on some 1d geoclaw routines not yet in
geoclaw, but in the directory `seismic/tsunami/geoclaw_src/1d`.

This initial configuration just sets a square pulse in `qinit.f90`.

The maximum depth at any time is recorded and written to `fort.hmax`.  This is
plotted as a red line by `setplot.py`.

