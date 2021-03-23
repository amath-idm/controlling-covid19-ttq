# Code for "Controlling COVID-19 via testing, tracing, and quarantine"

This repository includes the code for reproducing the results in the manuscript "Controlling COVID-19 via testing, tracing, and quarantine". The citation for the manuscript is:

> **Controlling COVID-19 via test-trace-quarantine**. Kerr CC, Mistry D, Stuart RM, Rosenfeld R, Hart G, Núñez RC, Selvaraj P, Cohen JA, Abeysuriya RG, George L, Hagedorn B, Jastrzębski M, Fagalde M, Duchin J, Famulare M, Klein DJ (under review; posted 2020-07-16). *medRxiv* 2020.07.15.20154765; doi: https://doi.org/10.1101/2020.07.15.20154765.

A webapp with an interactive version of the figures is available at http://ttq-app.covasim.org.


## Organization

The repository is organized as follows:

- `fig1_calibration` through `fig5_projections` are the main folders containing the code for reproducing each figure of the manuscript. In each of these folders, there is a file called `plot_figX.py` that will generate the corresponding figure from the manuscript. These scripts can be run directly, without the need to re-run the scripts for generating the underlying data (some of which take considerable CPU time). These are the scripts that are being run by the webapp.
- `inputs` and `outputs` are folders containing the input data and the model-based outputs, respectively.
- `webapp` includes the code for serving the figures via the interactive [webapp](http://ttq-app.covasim.org) or a local Jupyter notebook.


## Installation

Use `pip install -e .` to install dependencies, then run the scripts in the figure folders. Alternatively, run `pip install -e .[web]`, which will also let you run the notebook/webapp in the `webapp` folder. 

Note: installation is intended for Ubuntu 18.04/20.04 using a Python 3.8 [conda](https://www.anaconda.com/products/individual) virtual environment. Other platforms or environments may work with minor modifications to the scripts, but they have not been extensively tested.


## Further information

Further information on Covasim is available [here](http://docs.covasim.org). If you have further questions or would like technical assistance, please reach out to us at covasim@idmod.org.