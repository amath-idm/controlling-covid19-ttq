'''
Analyze the output of run_scenarios.py and convert to dataframes. This script
requires the output of run_scenarios.py, and the output of this script is used
by plot_fig4.py.
'''

import pandas as pd
import numpy as np
import sciris as sc
import pylab as pl
import seaborn as sns
import covasim as cv
import run_scenarios as rs

T = sc.tic()

which = 'a' # Switch between A and B here

regenerate = 1
if which == 'a':
    msimfile   = 'sensitivity_scenarios_sep28_final.msim'
    dffile     = 'sensitivity_scenarios_sep28_final.df'
else:
    msimfile   = 'reopening_sweeps_sep28_final.msim'
    dffile     = 'reopening_sweeps_sep28_final.df'
sim_start  = '2020-01-27'
start_date = '2020-06-01'
start_day  = cv.day(start_date, start_day=sim_start)
reff_window = 30 # Don't go for too long since hit herd immunity; take only this many days
reff2_window = 14 # Take this many days before, one infectious period
if which == 'a':
    f = rs.get_f_a()
else:
    f = rs.get_f_b()
fkeys = list(f.keys())


def splitlabel(lb):
    ''' Split a string label into pieces for a key '''
    out = sc.objdict()
    ief = lb.split(' ')
    out.ind = int(ief[0].split('=')[1])
    out.eind  = int(ief[1].split('=')[1])

    fv1 = ief[2].split('=')[1]
    fv2 = fv1.split(',') # New version
    fvals = np.array([float(v) for v in fv2])

    for k,v in zip(fkeys, fvals):
        out[k] = v
    return out


def get_params(pars):
    ''' Load parameter values '''
    out = sc.objdict()
    out.beta = pars['beta']
    for interv in pars['interventions']:
        if hasattr(interv, 'label'):
            if   interv.label == 'beta_c': out.last_bc_wc = interv.changes[-1]
            elif interv.label == 'beta_l': out.last_bc_l  = interv.changes[-1]
            elif interv.label == 'tn2':    out.last_tn    = interv.symp_test
    return out


def get_results(res):
    ''' Compute derived results '''
    out = sc.objdict()
    out.inf_start       = res['cum_infections'][start_day]
    out.inf_end         = res['cum_infections'][-1]
    out.diag_start      = res['cum_diagnoses'][start_day]
    out.diag_end        = res['cum_diagnoses'][-1]
    out.quar_start      = res['cum_quarantined'][start_day]
    out.quar_end        = res['cum_quarantined'][-1]
    out.cum_infections  = out.inf_end - out.inf_start
    out.cum_deaths      = res['cum_deaths'][-1] - res['cum_deaths'][start_day]
    out.cum_tests       = res['cum_tests'][-1] - res['cum_tests'][start_day]
    out.cum_diagnoses   = out.diag_end - out.diag_start
    out.cum_quarantined = out.quar_end - out.quar_start
    out.ave_infections  = res['new_infections'][start_day:].mean()
    out.ave_infectious  = res['n_infectious'][start_day:].mean()
    out.ave_tests       = res['new_tests'][start_day:].mean()
    out.ave_quarantined = res['new_quarantined'][start_day:].mean()
    out.log_infections  = np.log10(out.cum_infections)
    out.start_active    = res['n_infectious'][start_day]
    out.r_eff           = res['r_eff'][start_day:start_day+reff_window].mean()
    out.r_eff2          = out.ave_infectious/res['n_infectious'][start_day-reff2_window:start_day].mean()
    out.n_days          = len(res['r_eff'][start_day:]) # Calculate the number of days being averaged
    return out


def pairplotpars(df, inds=None, keys=None, color_column=None, bounds=None, cmap='parula', bins=None, edgecolor='w', facecolor='#F8A493', figsize=(20,16)):
    ''' Plot scatterplots, histograms, and kernel densities '''

    if inds is not None:
        df = df.iloc[inds,:].copy()

    # Choose the colors
    if color_column:
        colors = sc.vectocolor(df[color_column].values, cmap=cmap)
    else:
        colors = [facecolor for i in range(len(df))]
    df['color_column'] = [sc.rgb2hex(rgba[:-1]) for rgba in colors]

    if keys is not None:
        df = df.loc[:,keys+['color_column']].copy()

    # Make the plot
    grid = sns.PairGrid(df)
    grid = grid.map_lower(pl.scatter, **{'facecolors':df['color_column']})
    grid = grid.map_diag(pl.hist, bins=bins, edgecolor=edgecolor, facecolor=facecolor)
    grid = grid.map_upper(pl.hexbin, cmap="Blues", edgecolor="none" , gridsize=25)
    grid.fig.set_size_inches(figsize)
    grid.fig.tight_layout()

    # Set bounds
    if bounds:
        for ax in grid.axes.flatten():
            xlabel = ax.get_xlabel()
            ylabel = ax.get_ylabel()
            if xlabel in bounds:
                ax.set_xlim(bounds[xlabel])
            if ylabel in bounds:
                ax.set_ylim(bounds[ylabel])

    return grid


def findinds(df, *args):
    ''' Find matching indices in a large dataframe '''
    inds = np.arange(len(df))
    for arg in args:
        filterinds = sc.findinds(arg)
        inds = np.intersect1d(inds, filterinds)
    return inds


if __name__ == '__main__':

    if regenerate:
        print('Loading...')
        msim = cv.load(msimfile)

        print('Processing...')
        data = []
        for sim in msim.sims:
            info = splitlabel(sim.label)
            pars = get_params(sim.pars)
            results = get_results(sim.results)
            data.append(sc.mergedicts(info, pars, results))

        df = pd.DataFrame(data)
        cv.save(dffile, df)
    else:
        print('Reloading...')
        df = cv.load(dffile)

    print('Done.')
    sc.toc(T)
