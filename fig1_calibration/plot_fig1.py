'''
Create Fig. 1 from "Controlling COVID-19 via test-trace-quarantine". Relies on
output from run_fig1a.py (this output is already saved to the repository).
'''

import numpy as np
import pandas as pd
import sciris as sc
import datetime as dt
import pylab as pl
import seaborn as sns
import matplotlib.ticker as ticker
import covasim as cv


#%% Settings

# General settings
do_save = 0
do_show = 1
fig_path = 'ttq_fig1.png'

# Figure data settings
simsfile = 'fig1a.sims'
datafiles = sc.objdict()
datafiles.sg  = '../outputs/opt_merged_sep20_sg1.json'
datafiles.nsg = '../outputs/opt_merged_sep20_sg0.json'
safegraph_file = '../inputs/KC_weeklyinteractions_20200811_trim.csv'
scan_file = '../inputs/scanprev_5_21.csv'
mismatch_cutoff = 30 # Exclude simulations with a worse mismatch value than this
sims_cutoff = 9999 # Exclude sims with a higherindex than this
day_stride = 21 # Plot ticks in 3-week intervals
pl.rcParams['font.size'] = 20 # Set general figure options


#%% Helper functions

# Store plotted data
plotres = sc.objdict()

def format_ax(ax, sim, key=None):
    ''' Format the axes nicely '''
    @ticker.FuncFormatter
    def date_formatter(x, pos):
        return (sim['start_day'] + dt.timedelta(days=x)).strftime('%b-%d')
    ax.xaxis.set_major_formatter(date_formatter)
    if key != 'r_eff':
        sc.commaticks()
    pl.xlim([0, sim['n_days']])
    sc.boxoff()
    return


def plotter(key, sims, ax, ys=None, calib=False, label='', ylabel='', low_q=0.025, high_q=0.975, flabel=True, subsample=2):
    ''' Plot a single time series with uncertainty '''

    which = key.split('_')[1]
    try:
        color = cv.get_colors()[which]
    except:
        color = [0.5,0.5,0.5]
    if which == 'deaths':
        color = [0.5,0.0,0.0]

    if ys is None:
        ys = []
        for i,s in enumerate(sims):
            if i<sims_cutoff:
                ys.append(s.results[key].values)

    yarr = np.array(ys)
    best = pl.median(yarr, axis=0) # Changed from median to mean for smoother plots
    low  = pl.quantile(yarr, q=low_q, axis=0)
    high = pl.quantile(yarr, q=high_q, axis=0)

    sim = sims[0] # For having a sim to refer to

    tvec = np.arange(len(best))
    data, data_t = None, None
    if key in sim.data:
        data_t = np.array((sim.data.index-sim['start_day'])/np.timedelta64(1,'D'))
        inds = np.arange(0, len(data_t), subsample)
        data = sim.data[key][inds]
        pl.plot(data_t[inds], data, 'd', c=color, markersize=10, alpha=0.5, label='Data')

    end = None
    if flabel:
        if which == 'infections':
            fill_label = '95% predic-\ntion interval'
        else:
            fill_label = '95% prediction\ninterval'
    else:
        fill_label = None

    # Trim the beginning for r_eff and actually plot
    start = 2 if key == 'r_eff' else 0
    pl.fill_between(tvec[start:end], low[start:end], high[start:end], facecolor=color, alpha=0.2, label=fill_label)
    pl.plot(tvec[start:end], best[start:end], c=color, label=label, lw=4, alpha=1.0)

    sc.setylim()
    xmin,xmax = ax.get_xlim()
    ax.set_xticks(np.arange(xmin, xmax, day_stride))
    pl.ylabel(ylabel)

    plotres[key] = sc.objdict(dict(tvec=tvec, best=best, low=low, high=high, data=data, data_t=data_t))

    return


def plot_intervs(sim, labels=True):
    ''' Plot interventions as vertical lines '''

    color = [0.1, 0.4, 0.1]
    mar12 = sim.day('2020-03-12')
    mar23 = sim.day('2020-03-23')
    for day in [mar12, mar23]:
        pl.axvline(day, c=color, linestyle='--', alpha=0.4, lw=3)

    if labels:
        yl = pl.ylim()
        labely = yl[1]*1.03
        pl.text(mar12-25, labely, 'Schools close', color=color, alpha=0.9, style='italic')
        pl.text(mar23,    labely, 'Stay-at-home', color=color, alpha=0.9, style='italic')
    return


#%% Handle Fig. 1E

# Process data
data = sc.objdict()
dfs = sc.objdict()
for k,fn in datafiles.items():
    data[k] = sc.loadjson(datafiles[k])
    flat_data = []
    for entry in data[k]:
        flat = {}
        for k2 in ['index', 'mismatch']:
            flat[k2] = entry[k2]
        for k3,val3 in entry['pars'].items():
            flat[k3] = val3
            if k3 == 'beta':
                flat[k3] *= 100*3 # Household transmission probability and percentage
            elif k3.startswith('bc'):
                flat[k3] = 100*(1-flat[k3]) # To percentage, convert relative beta to beta reductions
        flat_data.append(flat)
    dfs[k] = pd.DataFrame(flat_data)


# Convert to dataframes
bestdfs = sc.objdict()
for k,df in dfs.items():
    df = df.loc[df['mismatch']<mismatch_cutoff,:]
    bestdfs[k] = df.drop(columns=['index', 'mismatch'])
    print('Samples remaining: ', len(df))
df1 = bestdfs[0]
df2 = bestdfs[1]
keys = list(bestdfs[0].columns)


# Define mapping to labels
mapping = sc.objdict({
    'beta'  : r'Overall $\beta$ (%)',
    'bc_wc1': r'Work/community $\beta$ reduction (%)',
    'bc_lf' : r'LTCF $\beta$ reduction (%)',
    'tn'   : 'Symptomatic testing OR',
    })
keys = mapping.keys()

for k in keys:
    q1 = 0.025
    q2 = 0.975
    sgb = np.median(df1[k])
    sgl = np.quantile(df1[k], q1)
    sgh = np.quantile(df1[k], q2)
    nsgb = np.median(df2[k])
    nsgl = np.quantile(df2[k], q1)
    nsgh = np.quantile(df2[k], q2)
    print(f'{k:8s}: SG={sgb:5.2f} ({sgl:5.2f}, {sgh:5.2f}) NSG={nsgb:5.2f} ({nsgl:5.2f}, {nsgh:5.2f})')


#%% Do the actual plotting

sims = cv.load(simsfile)
base_sim = sims[0] # For having a sim to refer to

def plot():

    # Create the figure
    fig = pl.figure(num='Fig. 1: Calibration', figsize=(24,20))
    tx1, ty1 = 0.005, 0.97
    tx2, ty2 = 0.545, 0.66
    ty3      = 0.34
    fsize = 40
    pl.figtext(tx1, ty1, 'a', fontsize=fsize)
    pl.figtext(tx1, ty2, 'b', fontsize=fsize)
    pl.figtext(tx2, ty1, 'c', fontsize=fsize)
    pl.figtext(tx2, ty2, 'd', fontsize=fsize)
    pl.figtext(tx1, ty3, 'e', fontsize=fsize)
    pl.figtext(tx2, ty3, 'f', fontsize=fsize)


    #%% Fig. 1A: diagnoses
    x0, y0, dx, dy = 0.055, 0.73, 0.47, 0.24
    ax1 = pl.axes([x0, y0, dx, dy])
    format_ax(ax1, base_sim)
    plotter('cum_diagnoses', sims, ax1, calib=True, label='Model', ylabel='Cumulative diagnoses')
    pl.legend(loc='lower right', frameon=False)


    #%% Fig. 1B: deaths
    y0b = 0.42
    ax2 = pl.axes([x0, y0b, dx, dy])
    format_ax(ax2, base_sim)
    plotter('cum_deaths', sims, ax2, calib=True, label='Model', ylabel='Cumulative deaths')
    pl.legend(loc='lower right', frameon=False)


    #%% Fig. 1A-B inserts (histograms)

    agehists = []

    for s,sim in enumerate(sims):
        agehist = sim['analyzers'][0]
        if s == 0:
            age_data = agehist.data
        agehists.append(agehist.hists[-1])

    # Observed data
    x = age_data['age'].values
    pos = age_data['cum_diagnoses'].values
    death = age_data['cum_deaths'].values

    # Model outputs
    mposlist = []
    mdeathlist = []
    for hists in agehists:
        mposlist.append(hists['diagnosed'])
        mdeathlist.append(hists['dead'])
    mposarr = np.array(mposlist)
    mdeatharr = np.array(mdeathlist)

    low_q = 0.025
    high_q = 0.975
    mpbest = pl.median(mposarr, axis=0)
    mplow  = pl.quantile(mposarr, q=low_q, axis=0)
    mphigh = pl.quantile(mposarr, q=high_q, axis=0)
    mdbest = pl.median(mdeatharr, axis=0)
    mdlow  = pl.quantile(mdeatharr, q=low_q, axis=0)
    mdhigh = pl.quantile(mdeatharr, q=high_q, axis=0)

    w = 4
    off = 2

    # Insets
    x0s, y0s, dxs, dys = 0.11, 0.84, 0.17, 0.13
    ax1s = pl.axes([x0s, y0s, dxs, dys])
    c1 = [0.3,0.3,0.6]
    c2 = [0.6,0.7,0.9]
    xx = x+w-off
    pl.bar(x-off,pos, width=w, label='Data', facecolor=c1)
    pl.bar(xx, mpbest, width=w, label='Model', facecolor=c2)
    for i,ix in enumerate(xx):
        pl.plot([ix,ix], [mplow[i], mphigh[i]], c='k')
    ax1s.set_xticks(np.arange(0,81,20))
    pl.xlabel('Age')
    pl.ylabel('Cases')
    sc.boxoff(ax1s)
    pl.legend(frameon=False, bbox_to_anchor=(0.7,1.1))

    y0sb = 0.53
    ax2s = pl.axes([x0s, y0sb, dxs, dys])
    c1 = [0.5,0.0,0.0]
    c2 = [0.9,0.4,0.3]
    pl.bar(x-off,death, width=w, label='Data', facecolor=c1)
    pl.bar(x+w-off, mdbest, width=w, label='Model', facecolor=c2)
    for i,ix in enumerate(xx):
        pl.plot([ix,ix], [mdlow[i], mdhigh[i]], c='k')
    ax2s.set_xticks(np.arange(0,81,20))
    pl.xlabel('Age')
    pl.ylabel('Deaths')
    sc.boxoff(ax2s)
    pl.legend(frameon=False)
    sc.boxoff(ax1s)


    #%% Fig. 1C: infections
    x0, dx = 0.60, 0.38
    ax3 = pl.axes([x0, y0, dx, dy])
    format_ax(ax3, sim)

    # Plot SCAN data
    pop_size = 2.25e6
    scan = pd.read_csv(scan_file)
    for i,r in scan.iterrows():
        label = "Data" if i==0 else None
        ts = np.mean(sim.day(r['since'], r['to']))
        low  = r['lower']*pop_size
        high = r['upper']*pop_size
        mean = r['mean']*pop_size
        ax3.plot([ts]*2, [low, high], alpha=1.0, color='k', zorder=1000)
        ax3.plot(ts, mean, 'o', markersize=7, color='k', alpha=0.5, label=label, zorder=1000)

    # Plot simulation
    plotter('cum_infections', sims, ax3, calib=True, label='Cumulative\ninfections\n(modeled)', ylabel='Infections')
    plotter('n_infectious', sims, ax3, calib=True, label='Active\ninfections\n(modeled)', ylabel='Infections', flabel=False)
    pl.legend(loc='upper left', frameon=False)
    pl.ylim([0, 130e3])
    plot_intervs(sim)



    #%% Fig. 1C: R_eff
    ax4 = pl.axes([x0, y0b, dx, dy])
    format_ax(ax4, sim, key='r_eff')
    plotter('r_eff', sims, ax4, calib=True, label='$R_{eff}$ (modeled)', ylabel=r'Effective reproduction number')
    pl.axhline(1, linestyle='--', lw=3, c='k', alpha=0.5)
    pl.legend(loc='upper right', frameon=False)
    plot_intervs(sim)



    #%% Fig. 1E

    # Do the plotting
    pl.subplots_adjust(left=0.04, right=0.52, bottom=0.03, top=0.35, wspace=0.12, hspace=0.50)

    for i,k in enumerate(keys):
        eax = pl.subplot(2,2,i+1)

        c1 = [0.2, 0.5, 0.8]
        c2 = [1.0, 0.5, 0.0]
        c3 = [0.1, 0.6, 0.1]
        sns.kdeplot(df1[k], shade=1, linewidth=3, label='', color=c1, alpha=0.5)
        sns.kdeplot(df2[k], shade=0, linewidth=3, label='', color=c2, alpha=0.5)

        pl.title(mapping[k])
        pl.xlabel('')
        pl.yticks([])
        if not i%4:
            pl.ylabel('Density')

        yfactor = 1.3
        yl = pl.ylim()
        pl.ylim([yl[0], yl[1]*yfactor])

        m1 = np.median(df1[k])
        m2 = np.median(df2[k])
        m1std = df1[k].std()
        m2std = df2[k].std()
        pl.axvline(m1, c=c1, ymax=0.9, lw=3, linestyle='--')
        pl.axvline(m2, c=c2, ymax=0.9, lw=3, linestyle='--')

        def fmt(lab, val, std=-1):
            if val<0.1:
                valstr = f'{lab} = {val:0.4f}'
            elif val<1.0:
                valstr = f'{lab} = {val:0.2f}±{std:0.2f}'
            else:
                valstr = f'{lab} = {val:0.1f}±{std:0.1f}'
            if std<0:
                valstr = valstr.split('±')[0] # Discard STD if not used
            return valstr

        if k.startswith('bc'):
            pl.xlim([0,100])
        elif k == 'beta':
            pl.xlim([3,5])
        elif k.startswith('tn'):
            pl.xlim([0,50])
        else:
            raise Exception(f'Please assign key {k}')

        xl = pl.xlim()
        xfmap = dict(
            beta   = 0.15,
            bc_wc1 = 0.30,
            bc_lf  = 0.35,
            tn     = 0.55,
        )

        xf = xfmap[k]
        x0 = xl[0] + xf*(xl[1]-xl[0])

        ypos1 = yl[1]*0.97
        ypos2 = yl[1]*0.77
        ypos3 = yl[1]*0.57

        if k == 'beta': # Use 2 s.f. instead of 1
            pl.text(x0, ypos1, f'M: {m1:0.2f} ± {m1std:0.2f}', c=c1)
            pl.text(x0, ypos2, f'N: {m2:0.2f} ± {m2std:0.2f}', c=c2)
            pl.text(x0, ypos3, rf'$\Delta$: {(m2-m1):0.2f} ± {(m1std+m2std):0.2f}', c=c3)
        else:
            pl.text(x0, ypos1, f'M: {m1:0.1f} ± {m1std:0.1f}', c=c1)
            pl.text(x0, ypos2, f'N: {m2:0.1f} ± {m2std:0.1f}', c=c2)
            pl.text(x0, ypos3, rf'$\Delta$: {(m2-m1):0.1f} ± {(m1std+m2std):0.1f}', c=c3)

        sc.boxoff(ax=eax)


    #%% Fig. 1F: SafeGraph
    x0, y0c, dyc = 0.60, 0.03, 0.30
    ax5 = pl.axes([x0, y0c, dx, dyc])
    format_ax(ax5, sim, key='r_eff')
    fn = safegraph_file
    df = pd.read_csv(fn)
    week = df['week']
    days = sim.day(week.values.tolist())
    s = df['p.tot.schools'].values*100
    w = df['p.tot.no.schools'].values*100

    # From Fig. 2
    colors = sc.gridcolors(5)
    wcolor = colors[3] # Work color/community
    scolor = colors[1] # School color

    pl.plot(days, w, 'd-', c=wcolor, markersize=15, lw=3, alpha=0.9, label='Workplace and\ncommunity mobility data')
    pl.plot(days, s, 'd-', c=scolor, markersize=15, lw=3, alpha=0.9, label='School mobility data')
    sc.setylim()
    xmin,xmax = ax5.get_xlim()
    ax5.set_xticks(np.arange(xmin, xmax, day_stride))
    pl.ylabel('Relative mobility (%)')
    pl.legend(loc='upper right', frameon=False)
    plot_intervs(sim)

    return fig

# Actually plot
fig = plot()


#%% Tidy up

if do_save:
    cv.savefig(fig_path, dpi=150)

if do_show:
    pl.show()


# Calculate stats
tot_inf = sc.objdict()
r_e = sc.objdict()
r_e_start = 10 # Days to start and end "initial" R_e calculation
r_e_end = 24
for key in ['best', 'low', 'high']:
    tot_inf[key] = plotres.cum_infections[key][-1]
    r_e[key] = plotres.r_eff[key][r_e_start:r_e_end+1].mean()

tot_diagnoses = plotres.cum_diagnoses.data[-1]
peak_active = plotres.n_infectious.best.max()

print('Stats for paper:')
print(f'Total infections: {tot_inf.best:0.1f} ({tot_inf.low:0.1f}, {tot_inf.high:0.1f})')
print(f'Total diagnoses: {tot_diagnoses:0.1f}')
print(f'Diagnosis rate: {tot_diagnoses/tot_inf.best*100:0.1f} ({tot_diagnoses/tot_inf.low*100:0.1f}, {tot_diagnoses/tot_inf.high*100:0.1f})')
print(f'R_e: {r_e.best:0.3f} ({r_e.low:0.3f}, {r_e.high:0.3f})')
print(f'Peak active: {peak_active:0.1f}')


print('Done.')
