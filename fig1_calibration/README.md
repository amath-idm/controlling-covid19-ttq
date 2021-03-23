# Fig. 1: Calibration of the model to data from Seattle-King County, Washington, from 27 January to 9 June 2020

This folder contains the code for reproducing Fig. 1 of the paper. The files are as follows:

- `plot_fig1.py` -- the script file that generates the figure.
- `cache_populations.py` -- generates the populations used by `create_sim.py` (both in this and other folders).
- `create_sim.py` -- configures the simulation for `run_fig1.py` and `run_optimization.py`.
- `run_optimization.py` -- runs the optimization process calibrating the model to data, for simulations in `run_fig1a.py`.
- `analyze_optimization.py`-- processes the calibrations and converts the outputs to JSON files with the calibrated parameters.
- `run_fig1a.py` -- uses the JSON files for the calibrated model to run the simulations.
- `fig1a.sims` -- the cached simulations (output of `run_fig1a.py`).

