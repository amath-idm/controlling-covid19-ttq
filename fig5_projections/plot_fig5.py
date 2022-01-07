'''
Create Fig. 5 from "Controlling COVID-19 via test-trace-quarantine". Relies on
output from run_fig5.py (this output is already saved to the repository).
'''

import numpy as np
import pandas as pd
import sciris as sc
import datetime as dt
import pylab as pl
import matplotlib.ticker as ticker
import covasim as cv

# General settings
do_save = 0
do_show = 1
fig_path = 'ttq_fig5.png'
pl.rcParams['font.size'] = 20 # Set general figure options
T = sc.tic()

# Load data
print('Loading data...')
msimsfile = 'fig5.msims'
tsfn = './kc_data/20200614chop5_KingCounty_Covasim_extended.xlsx' # Time series data
ctfn = './kc_data/contact_tracing.xlsx' # Contact tracing
tsdf = pd.read_excel(tsfn, engine='openpyxl')
ctdf = pd.read_excel(ctfn, sheet_name='final', engine='openpyxl')
ctdf['date'] = pd.to_datetime(ctdf['date']).dt.date
tsdf['date'] = pd.to_datetime(tsdf['date']).dt.date

msims = cv.load(msimsfile)
refsim = msims[0].sims[0] # A reference simulation, any will do


#%% Plot setup
scargs = dict(color='k', marker='d', linestyle='none', alpha=0.5, markersize=10, lw=3)
plargs = dict(lw=4)
fillargs = dict(alpha=0.2)
xlims1 = refsim.day('2020-06-01', '2020-08-31')
xlims2 = refsim.day('2020-01-27', '2020-08-31')

intervcolor = [1, 0.5, 0.05] # Was [0.1, 0.4, 0.1]
scencolors = dict(actual=[0.12, 0.47, 0.71], high=[0.1, 0.4, 0.0])
alpha = dict(actual=1.0, high=0.5)


def sub(data, wind=2):
    ''' Subsample data '''
    if hasattr(data, 'values'):
        data = data.values
    return data[::wind]


def roll(data):
    ''' Rolling window smoothing '''
    if hasattr(data, 'values'):
        data = data.values
    output = pd.Series(data).rolling(window=7).mean()
    return output


def simplotter(ax, key, label=None, each_sim=False):
    ''' Plot a simulation '''
    zorders = {'actual':10, 'high':0}
    linelabels = {'actual': 'Status quo\n(modeled)', 'high': 'High test + trace\n(modeled)'}
    for which in ['actual', 'high']:
        res = msims[which].results
        date = refsim.day(res['date']) # Convert to days
        best = res[key].values
        low  = res[key].low
        high = res[key].high
        ax.fill_between(date, roll(low), roll(high), facecolor=scencolors[which], **fillargs, zorder=zorders[which])
        ax.plot(date, roll(best), **plargs, c=scencolors[which], alpha=alpha[which], label=linelabels[which], zorder=zorders[which])
        ax.set_ylabel(label)
        if each_sim:
            for sim in msims[which].sims:
                qres = sim.results
                best = qres[key].values
                ax.plot(date, roll(best), **plargs, label=label + sim.label)
    return


def dataplotter(ax, date, data):
    ''' Plot data '''
    if hasattr(data, 'values'):
        data = data.values
    date = refsim.day(date.tolist())
    ax.plot(date, data, **scargs, label='Data', zorder=20)
    return


def format_axs(axs, key=None):
    ''' Format axes nicely '''
    @ticker.FuncFormatter
    def date_formatter(x, pos):
        # print(x)
        return (refsim['start_day'] + dt.timedelta(days=x)).strftime('%b-%d')
    for i,ax in enumerate(axs):
        bbox = None if i!=1 else (1.05,1.05) # Move legend up a bit
        day_stride = 21
        xmin,xmax = ax.get_xlim()
        ax.set_xticks(np.arange(xmin, xmax, day_stride))
        ax.xaxis.set_major_formatter(date_formatter)
        ax.legend(frameon=False, bbox_to_anchor=bbox)
        sc.boxoff(ax=ax)
        sc.setylim(ax=ax)
        sc.commaticks(ax=ax)
    return


def plot_intervs(ax):
    ''' Plot interventions '''
    color = intervcolor
    may31 = refsim.day('2020-05-31')
    ax.axvline(may31, c=color, linestyle='--', alpha=0.4, lw=3)
    yl = ax.get_ylim()
    labely = yl[1]*1.03
    ax.text(may31, labely, 'Stay-at-home lifted', color=color, alpha=0.9, style='italic', horizontalalignment='center')
    return


#%% Plotting
def plot():
    fig = pl.figure(num='Fig. 5: Projections and validation', figsize=(22, 14))

    x1 = 0.07 # Panel and text locations
    x2 = 0.48
    y1 = 0.05
    y2 = 0.55
    dx1 = 0.32
    dx2 = 0.51
    dy = 0.40
    ax1 = fig.add_axes([x1, y2, dx1, dy])
    ax2 = fig.add_axes([x1, y1, dx1, dy])
    ax3 = fig.add_axes([x2, y2, dx2, dy])
    ax4 = fig.add_axes([x2, y1, dx2, dy])
    axs = [ax1, ax2, ax3, ax4]

    fsize = 40
    tx1 = -0.06
    tx2 = -0.04
    ty1 = 0.40
    ty2 = 0.41
    pl.figtext(x1+tx1, y2+ty2, 'a', fontsize=fsize)
    pl.figtext(x1+tx1, y1+ty1, 'b', fontsize=fsize)
    pl.figtext(x2+tx2, y2+ty2, 'c', fontsize=fsize)
    pl.figtext(x2+tx2, y1+ty1, 'd', fontsize=fsize)

    # Panel 1
    dataplotter(ax1, sub(tsdf['date']), sub(roll(tsdf['new_tests'])))
    simplotter(ax1, 'new_tests', label='Tests conducted per day')
    ax1.set_xlim(xlims1)

    # Panel 2
    dataplotter(ax2, ctdf['date'], ctdf['estimated_daily'])
    simplotter(ax2, 'new_quarantined', label='Contacts traced per day')
    ax2.set_xlim(xlims1)

    # Panel 3
    simplotter(ax3, 'new_infections', label='New infections per day')
    plot_intervs(ax3)
    ax3.set_xlim(xlims2) # Starts later because of rolling window

    # Panel 4
    dataplotter(ax4, sub(tsdf['date']), sub(roll(tsdf['new_diagnoses'])))
    simplotter(ax4, 'new_diagnoses', label='New diagnoses per day')
    plot_intervs(ax4)
    ax4.set_xlim(xlims2)

    # Tidy up
    format_axs(axs)

    return fig

# Actually plot
fig = plot()


#%% Tidy up

if do_save:
    cv.savefig(fig_path, dpi=150)

if do_show:
    pl.show()

print('Done.')
sc.toc(T)