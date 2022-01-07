'''
Create Fig. 2 from "Controlling COVID-19 via test-trace-quarantine". Relies on
output from analyze_scenarios.py (this output is already saved to the repository).
'''

import numpy as np
import sciris as sc
import pylab as pl
import matplotlib as mpl
import statsmodels.api as sm
import covasim as cv
import create_sim as cs
import run_scenarios as rs

#%% Setup

# Basic setup
do_save = 0
do_show = 1
verbose = 0
fig_path = 'ttq_fig4.png'
dffile1 = 'sensitivity_scenarios_sep28_final.df'
dffile2 = 'reopening_sweeps_sep28_final.df'

T = sc.tic()

# Figure configuration
figw, figh = 24, 20
pl.rcParams['font.size'] = 20 # Set general figure options


#%% Top row: scatter plots
print('Calculating slopes...')
df1 = cv.load(dffile1)

ykey = ['cum_infections', 'r_eff'][0] # r_eff just for debugging
logy = ykey in ['cum_infections'] # Boolean for what to do with y
kcpop = 2.25e6


# Process baseline sims
f = rs.get_f_a()
fkeys = list(f.keys())
baseinds = np.arange(len(df1))
for k in fkeys:
    naninds = sc.findinds(np.isnan(df1[k].values))
    baseinds = np.intersect1d(baseinds, naninds)


df1map = sc.objdict(
    iqfactor  = 'Isolation/quarantine effectiveness',
    testprob  = 'Routine testing probability',
    testdelay = 'Swab-to-result delay',
    trprob    = 'Contact tracing probability',
    testqprob = 'Quarantine testing probability',
    trtime    = 'Contact tracing delay',
    )

xlabelmap = sc.objdict(
    iqfactor  = 'Reduction in transmission (%)',
    trprob    = 'Number of contacts traced per index case',
    testprob  = 'Routine tests per 1000 people per day',
    testqprob = 'Number of contacts tested per index case',
    testdelay = 'Average time for a test result to be returned (days)',
    trtime    = 'Average time for a contact to be traced (days)',
    )

slopelabels = sc.objdict(
    iqfactor  = 'infections averted\nper person fully isolated',
    trprob    = 'infections averted\nper index case',
    testprob  = 'infections averted\nper routine test',
    testqprob = 'infections averted\nper index case',
    testdelay = 'additional infections\nper day delay per positive test',
    trtime    = 'additional infections\nper day delay per index case',
    )

# Map the actual value onto the x coordinate
n_days = df1['n_days'][0] # All the same
max_delay = 7
cum_quar = df1['cum_quarantined']
cum_diag = df1['cum_diagnoses']
all_tests = df1['cum_tests']
quar_tests = cs.d_calcs.quar_test*cum_quar
non_quar_tests = all_tests - quar_tests

xvals = sc.objdict( # Can't use ones outside the one being used since nan
    iqfactor  = (1-df1['iqfactor'])*100, # Invert this so higher is more isolated
    trprob    = cum_quar/cum_diag, # Per case
    testprob  = non_quar_tests*1000/kcpop/n_days,
    testqprob = df1['testqprob']*cum_quar.values[baseinds].mean()/cum_diag.values[baseinds].mean(),
    testdelay = np.round(df1['testdelay']*max_delay),
    trtime    = np.round(df1['trtime']*max_delay),
)

n = len(df1)
ones = np.ones(n)
default_xvals = sc.objdict( # Can't use ones outside the one being used since nan
    iqfactor  = (1-cs.d_calcs.iso_factor*ones)*100, # Invert this so higher is more isolated
    trprob    = cum_quar/cum_diag, # Per case
    testprob  = non_quar_tests*1000/kcpop/n_days,
    testqprob = cs.d_calcs.quar_test*ones*cum_quar.values[baseinds].mean()/cum_diag.values[baseinds].mean(),
    testdelay = np.round(cs.d_calcs.t_delay*ones),
    trtime    = np.round(cs.d_calcs.tr_delay*ones),
)


def dvbm(key):
    ''' df1 values baseinds mean = dvbm '''
    return df1[key].values[baseinds].mean()

# Pull out the default values
base_cum_quar = dvbm('cum_quarantined')
base_cum_diag = dvbm('cum_diagnoses')

# Default positions of the lines
slopepoint = sc.objdict()
for k in fkeys:
    slopepoint[k] = default_xvals[k][baseinds].mean()

# Denominator for the slopes -- number of people at status quo
slopedenom = sc.objdict(
    iqfactor  = (1-cs.d_calcs.iso_factor)*base_cum_diag + (1-min(1,cs.d_calcs.iq_ratio*cs.d_calcs.iso_factor))*base_cum_quar,
    trprob    = base_cum_diag,
    testprob  = dvbm('cum_tests'),
    testqprob = (1*cs.d_calcs.quar_test)*base_cum_diag,
    testdelay = base_cum_diag,
    trtime    = base_cum_diag, # Denominator is people you start contact tracing from, not those traced
)

pointcolor = [0.1, 0.4, 0.0] # High mobility, high intervention
pointcolor2 = 'k' # Status quo, darker
c1 = [0.4, 0.4, 0.4]
c2 = [0.2, 0.2, 0.2]

# Plot each seed (eind) separately
sepinds = 0
if sepinds:
    einds = df1['eind'].values
    eis = np.unique(einds)
    cols = sc.gridcolors(len(eis))
if not sepinds:
    eis = [0]
    cols = [c1]


#%% Bottom row: surface plots
print('Calculating surfaces...')

df2 = cv.load(dffile2)

bottom = pl.cm.get_cmap('Oranges', 128)
top = pl.cm.get_cmap('Blues_r', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
newcmp = mpl.colors.ListedColormap(newcolors, name='OrangeBlue')

colormap = newcmp


def gauss2d(x, y, z, xi, yi, eps=1.0, xscale=1.0, yscale=1.0):
    ''' Gaussian 2D kernel '''

    # Doubles the speed to convert to 32 bit, functions faster than lambdas
    def arr32(arr): return np.array(arr, dtype=np.float32)
    def f32(x):     return np.float32(x)
    x, y, z, xi, yi = arr32(x), arr32(y), arr32(z), arr32(xi), arr32(yi)
    eps, xscale, yscale = f32(eps), f32(xscale), f32(yscale)

    # Actual computation
    nx = len(xi)
    ny = len(yi)
    zz = np.zeros((ny, nx))
    for i in range(nx):
        for j in range(ny):
            dist = np.sqrt(((x - xi[i])/xscale)**2 + ((y - yi[j])/yscale)**2)
            weights = np.exp(-(dist/eps)**2)
            weights = weights/np.sum(weights)
            val = np.sum(weights*z)
            zz[j,i] = val

    return np.array(zz, dtype=np.float64) # Convert back


def plot_surface(ax, dfr, col=0, colval=0):
    ''' Plot one of the surfaces '''

    all_tests = dfr['cum_tests']
    quar_tests = cs.d_calcs.quar_test*dfr['cum_quarantined']
    non_quar_tests = all_tests - quar_tests
    scaled_nq_tests = non_quar_tests*1000/kcpop/n_days

    x = scaled_nq_tests
    y = np.array(dfr['trprob'])
    z = np.array(dfr['r_eff'])
    min_x = 0
    min_y = 0
    max_x = 6
    max_y = 1.0
    min_z = 0.25
    max_z = 1.75

    eps = 0.08
    npts = 100
    xi = np.linspace(min_x, max_x, npts)
    yi = np.linspace(min_y, max_y, npts)
    xx, yy = np.meshgrid(xi, yi)
    zz = gauss2d(x, y, z, xi, yi, eps=eps, xscale=max_x-min_x, yscale=max_y-min_y)

    im = ax.contourf(xx, yy, zz, cmap=colormap, levels=np.linspace(min_z, max_z, 300))
    if col == 0:
        ax.set_ylabel('Contact tracing probability at home and work', labelpad=20)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    if col == 1:
        ax.set_xlabel('Number of routine tests per 1,000 people per day')

    axhandle = ax.contour(xx, yy, zz, levels=[1], colors='k', linestyles=':', linewidths=2)
    if   col == 0: ax.clabel(axhandle, fmt=r'$R_{e}=%2.1d$', colors='k', fontsize=18, manual=[(1, 0.3)])
    elif col == 1: ax.clabel(axhandle, fmt=r'$R_{e}=%2.1d$', colors='k', fontsize=18, manual=[(2, 0.4)])
    else:          ax.clabel(axhandle, fmt=r'$R_{e}=%2.1d$', colors='k', fontsize=18, manual=[(3, 0.5)])

    ax.set_title('%i%% mobility' % (colval*100))

    scolors = sc.vectocolor(z, cmap=colormap, minval=min_z, maxval=max_z)
    ax.scatter(x, y, marker='o', c=scolors, edgecolor=[0.3]*3, s=15, linewidth=0.1, alpha=0.5) # Add the measles plots

    ax.set_xlim([min_x, max_x])
    ax.set_ylim([min_y, max_y])
    if verbose:
        print(f'Plot: {col}, min zz: {zz.min():0.2f}; max zz: {zz.max():0.2f}; min z: {z.min():0.2f}; max z: {z.max():0.2f}')

    return im


#%% Actually plot

def plot():
    fig = pl.figure(num='Fig. 4: Suppression scenarios', figsize=(figw, figh))

    rx   = 0.07
    r1y  = 0.74
    rdx  = 0.26
    r1dy = 0.20
    rδ   = 0.30
    r1δy = 0.29
    r1ax = {}
    for i in range(6):
        xi = i%3
        yi = i//3
        r1ax[i] = pl.axes([rx+rδ*xi, r1y-r1δy*yi, rdx, r1dy])

    r2y = 0.05
    r2dy = rdx*figw/figh # To ensure square
    r2ax = {}
    for i in range(3):
        r2ax[i] = pl.axes([rx+rδ*i, r2y, rdx, r2dy])

    cax = pl.axes([0.96, r2y, 0.01, r2dy])

    # Labels
    lx = 0.015
    pl.figtext(lx, r1y+r1dy+0.02, 'a', fontsize=40)
    pl.figtext(lx, r2y+r2dy+0.02, 'b', fontsize=40)

    slopes = sc.objdict().make(keys=df1map.keys(), vals=[])
    slopes2 = sc.objdict().make(keys=df1map.keys(), vals=[])
    for plotnum,key,label in df1map.enumitems():
        cv.set_seed(plotnum)
        ax = r1ax[plotnum]
        for ei in eis:
            theseinds = sc.findinds(~np.isnan(df1[key].values))
            if sepinds:
                ei_ok = sc.findinds(df1['eind'].values == ei)
                theseinds = np.intersect1d(theseinds, ei_ok)
            x = xvals[key][theseinds]
            rawy = df1[ykey].values[theseinds]
            y = rawy/kcpop*100 if logy else rawy
            xm = x.max()
            if key in ['testdelay', 'trtime']:
                xnoise = 0.01
                ynoise = 0.05
            else:
                xnoise = 0
                ynoise = 0
            rndx = (np.random.randn(len(x)))*xm*xnoise
            rndy = (np.random.randn(len(y)))*ynoise
            ax.scatter(x+rndx, y*(1+rndy), alpha=0.2, c=[cols[ei]], edgecolor='w')

            # Calculate slopes
            slopey = np.log(rawy) if logy else rawy
            slopex = xvals[key].values[theseinds]
            tmp, res = np.polyfit(slopex, slopey, 1, cov=True)
            fitm, fitb = tmp
            factor = np.exp(fitm*slopepoint[key]+fitb) / slopedenom[key] * slopex.max()
            slope = fitm * factor
            slopes[key].append(slope)

            # Calculate slopes, method 2 -- used for std calculation
            X = sm.add_constant(slopex)
            mod = sm.OLS(slopey, X)
            res = mod.fit()
            conf = res.conf_int(alpha=0.05, cols=None)
            best = res.params[1]*factor
            high = conf[1,1]*factor
            slopes2[key].append(res)

            # Plot fit line
            bflx = np.array([x.min(), x.max()])
            ploty = np.log10(y) if logy else y
            plotm, plotb = np.polyfit(x, ploty, 1)
            bfly = plotm*bflx + plotb
            if logy:
                ax.semilogy(bflx, 10**(bfly), lw=3, c=c2)
            else:
                ax.plot(bflx, bfly, lw=3, c=c2)

        plot_baseinds = False
        if plot_baseinds:
            default_x = default_xvals[key][baseinds]
            default_rawy = df1[ykey].values[baseinds]
            default_y = default_rawy/kcpop*100 if logy else default_rawy
            ax.scatter(default_x, default_y, marker='x', alpha=1.0, c=[cols[i]])
        if verbose:
            print(f'Slope for {key:10s}: {np.mean(slopes[key]):0.3f} ± {high-best:0.3f}')

        sc.boxoff(ax=ax)
        ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        if plotnum in [0, 3]:
            if logy:
                ax.set_ylabel('Attack rate (%)')
                ax.set_yticks((0.01, 0.1, 1.0, 10, 100))
                ax.set_yticklabels(('0.01', '0.1', '1.0', '10', '100'))
            else:
                ax.set_ylabel(r'$R_{e}$')
        ax.set_xlabel(xlabelmap[key])

        if key in ['iqfactor']:
            ax.set_xticks(np.linspace(0,100,6))
        elif key in ['trprob']:
            ax.set_xticks(np.linspace(0,10,6))
        elif key in ['testprob']:
            ax.set_xticks(np.arange(7))
        elif key in ['testqprob']:
            ax.set_xticks(np.arange(7))
        else:
            ax.set_xticks(np.arange(8))

        if logy:
            ax.set_ylim([0.1,100])
        else:
            ax.set_ylim([0.7,1.7])

        xl = ax.get_xlim()
        if logy:
            ypos1 = 150
            ypos2 = 40
        else:
            ypos1 = 1.9
            ypos2 = 1.7

        xlpos = dict(
            iqfactor  = 0.86,
            testprob  = 0.13,
            trprob    = 0.83,
            testqprob = 0.86,
            testdelay = 0.13,
            trtime    = 0.00,
        )

        if key in ['iqfactor', 'testqprob', 'trprob']:
            align = 'right'
        else:
            align = 'left'

        ax.text((xl[0]+xl[1])/2, ypos1, label, fontsize=26, horizontalalignment='center')
        ax.text(xlpos[key]*xl[1], ypos2, f'{abs(best):0.2f} ± {high-best:0.2f} {slopelabels[key]}', color=pointcolor, horizontalalignment=align)
        ax.axvline(slopepoint[key], ymax=0.83, linestyle='--', c=pointcolor, alpha=0.5, lw=2)

    reop = [0.6, 0.8, 1.0]

    for ri, r in enumerate(reop):
        dfr = df2[df2['reopen'] == r]
        im = plot_surface(r2ax[ri], dfr, col=ri, colval=r)

        bbox = dict(facecolor='w', alpha=0.0, edgecolor='none')

        pointsize = 150
        if ri == 0:
            dotx = 1900*1000/kcpop
            doty = 0.06
            r2ax[ri].scatter([dotx], [doty], c='k', s=pointsize, zorder=10, marker='d')
            r2ax[ri].text(dotx*1.20, doty*1.50, 'Estimated\nconditions\non June 1', bbox=bbox)
        if ri == 1:
            dotx = 3000*1000/kcpop
            doty = 0.20
            r2ax[ri].scatter([dotx], [doty], c=[pointcolor2], s=pointsize, zorder=10, marker='d')
            r2ax[ri].text(dotx*1.10, doty*0.20, 'Estimated\nconditions\non July 15', color=pointcolor2, bbox=bbox)
        if ri == 2:
            dotx = 2.80# 7200*1000/kcpop
            doty = 0.70# 0.66
            r2ax[ri].scatter([dotx], [doty], c=[pointcolor], s=pointsize, zorder=10, marker='d')
            r2ax[ri].text(dotx*1.05, doty*1.05, 'High mobility,\nhigh test + trace\nscenario', color=pointcolor, bbox=bbox)

    cbar = pl.colorbar(im, ticks=np.linspace(0.4, 1.6, 7), cax=cax)
    cbar.ax.set_title('$R_{e}$', rotation=0, pad=20, fontsize=24)

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