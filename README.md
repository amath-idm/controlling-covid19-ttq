# Code for "Controlling COVID-19 via testing, tracing, and quarantine"

This repository includes the code for reproducing the results in the manuscript "Controlling COVID-19 via testing, tracing, and quarantine". The citation for the manuscript is:

> **Controlling COVID-19 via test-trace-quarantine**. Kerr CC, Mistry D, Stuart RM, Rosenfeld R, Hart G, Núñez RC, Selvaraj P, Cohen JA, Abeysuriya RG, George L, Hagedorn B, Jastrzębski M, Fagalde M, Duchin J, Famulare M, Klein DJ (under review; posted 2020-07-16). *medRxiv* 2020.07.15.20154765; doi: https://doi.org/10.1101/2020.07.15.20154765.

A webapp with an interactive version of the figures is available at http://ttq-app.covasim.org.


## Organization

The repository is organized as follows:

- `fig1` through `fig5` are the main folders containing the code for reproducing each figure of the manuscript. In each of these folders, there is a file called `figX_plot.py` which will generate the corresponding figure from the manuscript.
- `inputs` and `outputs` are folders containing the input data and the model-based outputs, respctively.
- `webapp` includes the code for serving the figures via an interactive webapp or Jupyter notebook.


## Installation

To install dependencies, use `pip install -e .`, then run the scripts in the figure folders. Alternatively, run `pip install -e .[web]` and run the notebook/webapp in the `webapp` folder. 

Note: installation is intended for Ubuntu 18.04/20.04 using a Python 3.8 [conda](https://www.anaconda.com/products/individual) virtual environment. Other platforms or environments may work, but they have not been extensively tested.


## Further information

If you have further questions or would like technical assistance, please reach out to us at covasim@idmod.org.