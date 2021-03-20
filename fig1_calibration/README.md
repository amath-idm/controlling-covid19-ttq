# Fig. 1: Calibration of the model to data from Seattle-King County, Washington, from 271 JanuaryFebruary to 9 June 2020

This folder contains the code for reproducing Fig. 1 of the manuscript. The files are as follows:

- `fig1_plot.py` -- the script file that generates the figure.
- `cache_populations.py` -- generates the populations used by `create_sim.py`.
- `create_sim.py` -- configures the simulation for `run_fig1.py` and `run_optimization.py`.
- `run_optimization.py` -- runs the optimization process calibrating the model to data, for simulations in `run_fig1a.py`.
- `run_fig1a.py` -- uses the calibrated model to run the simulations.
- `fig1a_sims_sep20_sg1.sims` -- the cached simulations (output of `run_fig1a.py`).

