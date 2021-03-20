'''
Create Fig. 3 from "Controlling COVID-19 via test-trace-quarantine". Relies on
output from run_fig3.py (this output is already saved to the repository).
'''

import numpy as np
import pylab as pl
import sciris as sc
import covasim as cv

do_save = 0
do_show = 1
fig_path = 'ttq_fig3_sep25.png'

font_size = 22
font_family = 'Proxima Nova'
pl.rcParams['font.size'] = font_size
pl.rcParams['font.family'] = font_family
T = sc.tic()

#%% Left hand side: transmission trees

iday = 20 # Day interventions begin
p = sc.objdict(
    pop_size = 100,
    n_days = 125,
    pop_type = 'hybrid',
    pop_infected = 1,
    rand_seed = 8, # For the tree, select a random seed where the epidemic doesn't immediately die out
    verbose = 0,
    iso_factor  = dict(h=0.3, s=0.0, w=0.0, c=0.1),
    quar_factor = dict(h=0.8, s=0.0, w=0.0, c=0.3),
)

s = sc.objdict()
t = sc.objdict()
a = sc.objdict()

fig = pl.figure(num='Fig. 3: Theoretical TTQ', figsize=(24,14))

euclid = False
ay0 = 0.06
adx = 0.2
ady = 0.91
axa, axb, axc = 0.04, 0.27, 0.49
a['none']  = pl.axes([axa, ay0, adx, ady])
a['test']  = pl.axes([axb, ay0, adx, ady])
a['trace'] = pl.axes([axc, ay0, adx, ady])

xoff = 0.02
yoff = 0.05
pl.figtext(axa-xoff, ady+yoff, 'A', fontsize=40)
pl.figtext(axb-xoff, ady+yoff, 'B', fontsize=40)
pl.figtext(axc-xoff, ady+yoff, 'C', fontsize=40)

max_n = 0; ζ = {'↓':-2, '↑':10} # Configure plot zorder

for k in ['none', 'test', 'trace']:
    interventions = []
    if k in ['test', 'trace']:
        interventions += [cv.test_prob(start_day=iday, symp_prob=0.15, symp_quar_prob=1.0, asymp_quar_prob=1.0, quar_policy='start')]
    if k in ['trace']:
        interventions += [cv.contact_tracing(start_day=iday, trace_probs=dict(h=0.7, s=0.1, w=0.1, c=0.0))]

    s[k] = cv.Sim(p, interventions=interventions)
    s[k].run()
    t[k] = s[k].make_transtree()
    tt = t[k]
    max_n = max(max_n, len(tt.infection_log))

    # Get the history of each chain
    hist = [None]*p.pop_size
    for i,entry in enumerate(tt.detailed):
        if entry:
            hist[i] = [i]
            source = entry['source']
            while source is not None:
                hist[i].insert(0, source)
                source = tt.detailed[source]['source']

    # Figure out order
    histstrs = []
    for h in hist:
        if h is not None:
            string = ''.join([str(n) for n in h])
            histstrs.append(string)
    iorder = np.arange(p.pop_size)
    inf_inds = cv.false(s[k].people.susceptible)
    iorder[inf_inds] = inf_inds[np.argsort(np.argsort(histstrs))] # 🤔🙈
    min_inf_ind = inf_inds.min()
    min_ind_ind = sc.findinds(iorder==min_inf_ind)[0]
    orig_0_ind = sc.findinds(iorder==0)[0]
    iorder[min_ind_ind] = 0
    iorder[orig_0_ind] = min_inf_ind

    # Initialization
    n = p.n_days + 1
    frames = [list() for i in range(n)]
    quar_color  = '#6495ed'
    diag_color  = '#80d380'
    inf_color   = '#b22222'
    quar_color2 = '#b3ceff'

    # Construct each frame of the animation
    for ddict in tt.detailed:  # Loop over every person
        if ddict is None:
            continue # Skip the 'None' node corresponding to seeded infections

        frame = {}
        tdq = {}  # Short for "tested, diagnosed, or quarantined"
        target_ind = ddict['target']

        if not np.isnan(ddict['date']): # If this person was infected
            source_ind = ddict['source'] # Index of the person who infected the target
            target_date = ddict['date']
            if source_ind is not None:  # Seed infections and importations won't have a source
                source_date = tt.detailed[source_ind]['date']
            else:
                source_ind = target_ind
                source_date = 0

            # Construct this frame
            frame['x'] = [source_date, target_date]
            frame['y'] = [source_ind, target_ind]
            tdict = ddict['t']
            is_quar = not np.isnan(tdict['date_quarantined'])
            is_diag = not np.isnan(tdict['date_diagnosed'])
            frame['c'] = inf_color
            frames[int(target_date)].append(frame)

    # Plot
    pl.sca(a[k])
    print('Working on', k)
    for frame in frames:
        for entry in frame:
            x0 = entry['x'][0]
            x1 = entry['x'][1]
            y0 = iorder[entry['y'][0]]
            y1 = iorder[entry['y'][1]]
            color = entry['c']
            if euclid:
                pl.plot([x0, x1], [y0, y1], lw=3, c=[0.7]*3)
            else:
                pl.plot([x0, x0, x1], [y0, y1, y1], lw=3, c=[0.7]*3)
            pl.scatter([x1], [y1], s=120, zorder=ζ['↑'], c=[color])

    for i in range(p.pop_size):
        ii = iorder[i]
        dq = s[k].people.date_quarantined[i]
        dd = s[k].people.date_diagnosed[i]
        de = s[k].people.date_exposed[i]
        qdy = 0.2
        α = 0.9
        mark = ['>', 'x'][1]
        if not np.isnan(dq) and not np.isnan(de) and np.isnan(dd):
            print(f'Person {i} is infected and quarantined but not diagnosed, due to presymptomatic period')
        if not np.isnan(dq):
            pl.plot([de, dq], [ii, ii], c=quar_color2, linestyle='--', zorder=ζ['↓']-4)
            pl.fill_between([dq, dq+14], [ii-qdy]*2, [ii+qdy]*2, facecolor=quar_color2, zorder=ζ['↓']-2, alpha=α)
            if np.isnan(de):
                pl.scatter([dq], [ii], s=120, c='w', edgecolor=quar_color2, lw=2, alpha=1.0, zorder=ζ['↓'])
        if not np.isnan(dd):
            if not np.isnan(dq):
                color = quar_color2
            else:
                color = diag_color
            pl.plot([de, dd], [ii, ii], c=color, linestyle='--', zorder=ζ['↓']-4)
            pl.scatter([dd], [ii], marker=mark, linewidth=3, s=120, c=[color], zorder=ζ['↓']+1)
            pl.fill_between([dd, dd+14], [ii-qdy]*2, [ii+qdy]*2, facecolor=diag_color, zorder=ζ['↓']-1, alpha=α)
    sc.boxoff()

    a[k].spines['left'].set_visible(False)
    if k != 'none':
        a[k].set_yticks([])
    else:
        a[k].set_yticks(np.arange(0,p.pop_size+1,10))
        a[k].set_ylabel('Person')
    pl.xlabel('Days since seed infection')
    ylimmap = dict(none=1.16, test=1.16, trace=1.16)
    pl.ylim([-1, p.pop_size*ylimmap[k]])

    # Labels
    n_inf  = len(tt.infection_log)
    n_diag_dir = len(np.intersect1d(cv.defined(s[k].people.date_diagnosed), cv.undefined(s[k].people.date_quarantined)))
    n_diag_quar = len(np.intersect1d(cv.defined(s[k].people.date_diagnosed), cv.defined(s[k].people.date_quarantined)))
    n_quar = len(np.intersect1d(cv.defined(s[k].people.date_quarantined), cv.true(s[k].people.susceptible)))
    txtargs = dict(horizontalalignment='right')

    ps = p.pop_size
    ty1, ty2, ty25, ty3, ty4, ty5 = ps*1.02, ps*1.05, ps*1.065, ps*1.08, ps*1.11, ps*1.14
    δy = 0.8

    xl = pl.xlim()
    tx0map = {
        'none': xl[1]*0.6,
        'test': xl[1]*0.6,
        'trace': xl[1]*0.8,
        }

    tx0 = tx0map[k]
    tx1 = tx0 + xl[1]*0.04
    tx2 = tx0 + xl[1]*0.11
    tx3 = tx0 + xl[1]*0.25

    if k == 'none':
        pl.text(tx0, ty25, f'Infected: {n_inf:2d}', **txtargs)
        pl.scatter([tx1], [ty25+δy], s=120, c=inf_color)

    elif k == 'test':
        pl.text(tx0, ty3, f'Infected: {n_inf:2d}', **txtargs)
        pl.scatter([tx1], [ty3+δy], s=120, c=inf_color)

        pl.text(tx0, ty2, f'Diagnosed: {n_diag_dir:2d}', **txtargs)
        pl.plot([tx1, tx2], [ty2+δy]*2, c=diag_color, linestyle='--')
        pl.scatter([tx2], [ty2+δy], marker=mark, linewidth=3, s=120, c=diag_color)
        pl.fill_between([tx2, tx3], [ty2+δy-qdy]*2, [ty2+δy+qdy]*2, facecolor=diag_color)

    elif k == 'trace':
        pl.text(tx0, ty5, f'Infected: {n_inf:2d}', **txtargs)
        pl.scatter([tx1], [ty5+δy], s=120, c=inf_color)

        pl.text(tx0, ty4, f'Diagnosed directly: {n_diag_dir:2d}', **txtargs)
        pl.plot([tx1, tx2], [ty4+δy]*2, c=diag_color, linestyle='--')
        pl.scatter([tx2], [ty4+δy], marker=mark, linewidth=3, s=120, c=diag_color)

        pl.text(tx0, ty3, f'Diagnosed via tracing: {n_diag_quar:2d}', **txtargs)
        pl.plot([tx1, tx2], [ty3+δy]*2, c=quar_color2, linestyle='--')
        pl.scatter([tx2], [ty3+δy], marker=mark, linewidth=3, s=120, c=quar_color2)

        pl.text(tx0, ty2, f'Isolated after diagnosis: {n_diag_dir+n_diag_quar:2d}', **txtargs)
        pl.fill_between([tx1, tx3], [ty2+δy-qdy]*2, [ty2+δy+qdy]*2, facecolor=diag_color)

        pl.text(tx0, ty1, f'Quarantined: {n_quar:2d}', **txtargs)
        pl.scatter([tx1], [ty1+δy], s=120, c=quar_color)
        pl.scatter([tx1], [ty1+δy], s=120, c='w', edgecolor=quar_color2, lw=2, zorder=10)
        pl.fill_between([tx1, tx3], [ty1+δy-qdy]*2, [ty1+δy+qdy]*2, facecolor=quar_color2)



#%% Right hand side: theoretical TTQ

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


# Configure

betadict = dict(best=0.018, low=0.014, high=0.022)
scenkeys = ['beta', 'test', 'trace']
n_seeds = 10
simsfile = 'fig3.sims'

sims = sc.objdict()
simsobj = cv.load(simsfile)
for sim in simsobj:
    sims[sim.label] = sim


# Plotting

colors = dict(beta='#b22292', test='#80d380', trace='#6495ed')
betamap = dict(best='Medium', low='Low', high='High')
scenmap = dict(beta='Distancing', test='Testing', trace='TTQ')

for betakey in betadict.keys():
    for scenkey in scenkeys:
        tvec = sims[0].results['t']
        average = np.zeros(sims[0].npts)
        for seed in np.arange(n_seeds):
            label = f'{betakey}_{scenkey}_{seed}'
            res = sims[label].results
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
    betanorm = betadict[betakey]*100*3*0.78 # Normalized to match Fig. 1 beta
    ax.set_title(rf'{betamap[betakey]} transmission, $\beta$ = {betanorm:0.1f}%')
    sc.boxoff(ax=ax)
    if betakey == 'best':
        ax.legend(frameon=False)
    if betakey == 'high':
        ax.set_yticks(np.arange(8)*2e3)
        ax.set_xlabel('Days since seed infections')


#%% Tidy up

if do_save:
    cv.savefig(fig_path, dpi=150)

if do_show:
    pl.show()

print('Done.')
sc.toc(T)