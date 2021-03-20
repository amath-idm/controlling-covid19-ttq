'''
Run the simulations used in Fig. 3E-F. The file saved by this script is used by
fig3_plot.py to generate the figure.

The simulations take about 5 min to run.
'''

import numpy as np
import pylab as pl
import sciris as sc
import covasim as cv

# Basic setup
do_run = 1
do_save = 0
simsfile = 'fig3.sims'
T = sc.tic()
s = sc.objdict()
a = sc.objdict()

# Figure configuration
font_size = 22
font_family = 'Proxima Nova'
pl.rcParams['font.size'] = font_size
pl.rcParams['font.family'] = font_family
fig = pl.figure(figsize=(24,14))

# Placeholder axes
ay0 = 0.06
adx = 0.2
ady = 0.91
axa, axb, axc = 0.03, 0.26, 0.49
a['none']  = pl.axes([axa, ay0, adx, ady])
a['test']  = pl.axes([axb, ay0, adx, ady])
a['trace'] = pl.axes([axc, ay0, adx, ady])
xoff = 0.02
yoff = 0.05
pl.figtext(axa-xoff, ady+yoff, 'A', fontsize=40)
pl.figtext(axb-xoff, ady+yoff, 'B', fontsize=40)
pl.figtext(axc-xoff, ady+yoff, 'C', fontsize=40)

# For this figure
adx = 0.21
ady = 0.25
hspace = 0.32
aya, ayb, ayc = ay0+2*hspace, ay0+hspace, ay0
axx = 0.77
a['best'] = pl.axes([axx, aya, adx, ady])
a['low']  = pl.axes([axx, ayb, adx, ady])
a['high'] = pl.axes([axx, ayc, adx, ady])
xoff = 0.05
yoff = 0.26
pl.figtext(axx-xoff, aya+yoff, 'D', fontsize=40)
pl.figtext(axx-xoff, ayb+yoff, 'E', fontsize=40)
pl.figtext(axx-xoff, ayc+yoff, 'F', fontsize=40)

max_n = 0; ζ = {'↓':-2, '↑':10}


#%% Configure for running

p = sc.objdict(
    pop_type     = 'hybrid',
    pop_size     = 30e3,
    pop_infected = 100,
    start_day    = '2020-02-01',
    end_day      = '2020-07-01',
    rand_seed    = 1,
    verbose      = 0.1,
)

# Overall settings
betadict = dict(best=0.018, low=0.014, high=0.022)
n_seeds = 10

# Define scenarios
sd = 10
ed = 15
bc = 1.0
baseline_testing = cv.test_num(start_day=0, end_day=sd-1, daily_tests=p.pop_size/200)

scenarios = dict(
    beta = [
        baseline_testing,
        cv.change_beta(days=[sd, ed], changes=[bc, 0.40]),
        cv.test_prob(start_day=sd, symp_prob=0.075, asymp_prob=0.0075, test_delay=4),
        ],
    test = [
        baseline_testing,
        cv.change_beta(days=[sd, ed], changes=[bc, 1.0]),
        cv.test_prob(start_day=sd, symp_prob=0.75, asymp_prob=0.075, test_delay=0),
        ],
    trace = [
        baseline_testing,
        cv.change_beta(days=[sd, ed], changes=[bc, 1.0]),
        cv.test_prob(start_day=sd, symp_prob=0.08, asymp_prob=0.008, symp_quar_prob=0.75, asymp_quar_prob=0.5, test_delay=0, quar_policy='both'),
        cv.contact_tracing(start_day=sd, trace_probs=0.9, trace_time=0),
        ],
)

# Configure and run sims
if do_run:

    # Configure sims
    sims = []
    for betakey,beta in betadict.items():
        for scenkey,interventions in scenarios.items():
            for seed in np.arange(n_seeds):
                sim = cv.Sim(p, beta=beta, rand_seed=seed, interventions=interventions, label=f'{betakey}_{scenkey}_{seed}')
                sims += [sim]

    # Run sims
    msim = cv.MultiSim(sims)
    msim.run(par_args={'ncpus':5})
    sims = msim.sims
    if do_save:
        cv.save(simsfile, sims)

# Loading from disk
else:
    print(f'Loading {simsfile}...')
    sims = cv.load(simsfile)

for sim in sims:
    s[sim.label] = sim


#%% Plotting
colors = dict(beta='#b22222', test='#80d380', trace='#6495ed')
betamap = dict(best='Medium', low='Low', high='High')
scenmap = dict(beta='Distancing', test='Testing', trace='TTQ')

for betakey in betadict.keys():
    for scenkey in scenarios.keys():
        tvec = s[0].results['t']
        average = np.zeros(s[0].npts)
        for seed in np.arange(n_seeds):
            label = f'{betakey}_{scenkey}_{seed}'
            res = s[label].results
            cum_infections = res['cum_infections'].values
            a[betakey].plot(tvec, cum_infections, c=colors[scenkey], lw=3, alpha=0.2)
            average += cum_infections/n_seeds
        a[betakey].plot(tvec, average, '--', c=colors[scenkey], lw=3, alpha=1.0, zorder=10, label=scenmap[scenkey])

for betakey in betadict.keys():
    ax = a[betakey]
    sc.commaticks(ax=ax)
    sc.setylim(ax=ax)
    ax.set_xlim([0,150])
    ax.set_ylabel('Cumulative infections')
    betanorm = betadict[betakey]*7*100 # Convert from a relative to absolute beta
    ax.set_title(rf'{betamap[betakey]} transmission, $\beta$ = {betanorm:0.0f}%')
    sc.boxoff(ax=ax)
    if betakey == 'best':
        ax.legend(frameon=False)
    if betakey == 'high':
        ax.set_yticks(np.arange(8)*2e3)
        ax.set_xlabel('Days since seed infections')


pl.show()

print('Done.')
sc.toc(T)