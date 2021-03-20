'''
Create Fig. 2 from "Controlling COVID-19 via test-trace-quarantine". Relies on
output from run_fig2.py (this output is already saved to the repository).
'''

import matplotlib.dates as mdates
import sciris as sc
import pylab as pl
import numpy as np
import covasim as cv

#%% Settings

# General settings
do_save = 0
do_show = 1
fig_path = 'ttq_fig2_sep25.png'
simfile = 'fig2.sim'  # File to load -- produced by fig2_run.py
T = sc.tic()

# Set general figure options
font_size = 22
font_family = 'Proxima Nova'
pl.rcParams['font.size'] = font_size
pl.rcParams['font.family'] = font_family
pieargs = dict(startangle=90, counterclock=False, labeldistance=1.25)

# Load the sim
print('Loading... (takes ~5 s)')
sim = cv.load(simfile)
tt = sim.results.transtree

# Make the plots
print('Plotting...')
fig = pl.figure(num='Fig. 2: Transmission dynamics', figsize=(20,14))
piey, tsy, r3y = 0.68, 0.50, 0.07
piedx, tsdx, r3dx = 0.2, 0.9, 0.25
piedy, tsdy, r3dy = 0.2, 0.47, 0.35
pie1x, pie2x = 0.12, 0.65
tsx = 0.07
dispx, cumx, sympx = tsx, 0.33+tsx, 0.66+tsx
ts_ax   = pl.axes([tsx, tsy, tsdx, tsdy])
pie_ax1 = pl.axes([pie1x, piey, piedx, piedy])
pie_ax2 = pl.axes([pie2x, piey, piedx, piedy])
symp_ax = pl.axes([sympx, r3y, r3dx, r3dy])
disp_ax = pl.axes([dispx, r3y, r3dx, r3dy])
cum_ax  = pl.axes([cumx, r3y, r3dx, r3dy])

off = 0.06
txtdispx, txtcumx, txtsympx = dispx-off, cumx-off, sympx-off+0.02
tsytxt = tsy+tsdy
r3ytxt = r3y+r3dy
pl.figtext(txtdispx, tsytxt, 'A', fontsize=40)
pl.figtext(txtdispx, r3ytxt, 'B', fontsize=40)
pl.figtext(txtcumx, r3ytxt, 'C', fontsize=40)
pl.figtext(txtsympx,  r3ytxt, 'D', fontsize=40)


#%% Fig. 2A -- Time series plot

layer_keys = list(sim.layer_keys())
layer_mapping = {k:i for i,k in enumerate(layer_keys)}
n_layers = len(layer_keys)
colors = sc.gridcolors(n_layers)

layer_counts = np.zeros((sim.npts, n_layers))
for source_ind, target_ind in tt.transmissions:
    dd = tt.detailed[target_ind]
    date = dd['date']
    layer_num = layer_mapping[dd['layer']]
    layer_counts[date, layer_num] += sim.rescale_vec[date]

mar12 = cv.date('2020-03-12')
mar23 = cv.date('2020-03-23')
mar12d = sim.day(mar12)
mar23d = sim.day(mar23)

labels = ['Household', 'School', 'Workplace', 'Community', 'LTCF']
for l in range(n_layers):
    ts_ax.plot(sim.datevec, layer_counts[:,l], c=colors[l], lw=3, label=labels[l])
sc.setylim(ax=ts_ax)
sc.boxoff(ax=ts_ax)
ts_ax.set_ylabel('Transmissions per day')
ts_ax.set_xlim([sc.readdate('2020-01-18'), sc.readdate('2020-06-09')])
ts_ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
ts_ax.set_xticks([sim.date(d, as_date=True) for d in np.arange(0, sim.day('2020-06-09'), 14)])
ts_ax.legend(frameon=False, bbox_to_anchor=(0.85,0.1))

color = [0.2, 0.2, 0.2]
ts_ax.axvline(mar12, c=color, linestyle='--', alpha=0.4, lw=3)
ts_ax.axvline(mar23, c=color, linestyle='--', alpha=0.4, lw=3)
yl = ts_ax.get_ylim()
labely = yl[1]*1.015
ts_ax.text(mar12, labely, 'Schools close                     ', color=color, alpha=0.9, style='italic', horizontalalignment='center')
ts_ax.text(mar23, labely, '                   Stay-at-home', color=color, alpha=0.9, style='italic', horizontalalignment='center')


#%% Fig. 2A inset -- Pie charts

pre_counts = layer_counts[0:mar12d, :].sum(axis=0)
post_counts = layer_counts[mar23d:, :].sum(axis=0)
pre_counts = pre_counts/pre_counts.sum()*100
post_counts = post_counts/post_counts.sum()*100

lpre = [
    f'Household\n{pre_counts[0]:0.1f}%',
    f'School\n{pre_counts[1]:0.1f}% ',
    f'Workplace\n{pre_counts[2]:0.1f}%    ',
    f'Community\n{pre_counts[3]:0.1f}%',
    f'LTCF\n{pre_counts[4]:0.1f}%',
]

lpost = [
    f'Household\n{post_counts[0]:0.1f}%',
    f'School\n{post_counts[1]:0.1f}%',
    f'Workplace\n{post_counts[2]:0.1f}%',
    f'Community\n{post_counts[3]:0.1f}%',
    f'LTCF\n{post_counts[4]:0.1f}%',
]

pie_ax1.pie(pre_counts, colors=colors, labels=lpre, **pieargs)
pie_ax2.pie(post_counts, colors=colors, labels=lpost, **pieargs)

pie_ax1.text(0, 1.75, 'Transmissions by layer\nbefore schools closed', style='italic', horizontalalignment='center')
pie_ax2.text(0, 1.75, 'Transmissions by layer\nafter stay-at-home', style='italic', horizontalalignment='center')


#%% Fig. 2B -- histogram by overdispersion

# Process targets
n_targets = tt.count_targets()

# Handle bins
max_infections = n_targets.max()
bins = np.arange(0, max_infections+2)

# Analysis
counts = np.histogram(n_targets, bins)[0]
bins = bins[:-1] # Remove last bin since it's an edge
total_counts = counts*bins
total_counts = total_counts/total_counts.sum()*100
n_bins = len(bins)
index = np.linspace(0, 100, len(n_targets))
sorted_arr = np.sort(n_targets)
sorted_sum = np.cumsum(sorted_arr)
sorted_sum = sorted_sum/sorted_sum.max()*100
change_inds = sc.findinds(np.diff(sorted_arr) != 0)

pl.set_cmap('Spectral_r')
sscolors = sc.vectocolor(n_bins)

width = 1.0
for i in range(n_bins):
    disp_ax.bar(bins[i], total_counts[i], width=width, facecolor=sscolors[i])
disp_ax.set_xlabel('Number of transmissions per case')
disp_ax.set_ylabel('Proportion of all transmissions (%)')
sc.boxoff()
disp_ax.set_xlim([0.5, 32.5])
disp_ax.set_xticks(np.arange(0, 32.5, 4))
sc.boxoff(ax=disp_ax)

dpie_ax = pl.axes([dispx+0.05, 0.20, 0.2, 0.2])
trans1 = total_counts[1:3].sum()
trans2 = total_counts[3:6].sum()
trans3 = total_counts[6:].sum()
labels = [
    f'1-2: {trans1:0.0f}%',
    f'3-5: {trans2:0.0f}%',
    f'>5: {trans3:0.0f}%',
    ]
dpie_args = sc.mergedicts(pieargs, {'labeldistance':1.1}) # Slightly smaller label distance
dpie_ax.pie([trans1, trans2, trans3], labels=labels, colors=sscolors[[0,3,6]], **dpie_args)


#%% Fig. 2C -- cumulative distribution function

index = 100 - index
sorted_sum = 100 - sorted_sum
change_inds = change_inds[::-1]
n_change_inds = len(change_inds)
for i in range(n_change_inds):
    ci = change_inds[i]
    ici = index[ci]
    sci = sorted_sum[ci]
    if i>0:
        cim1 = change_inds[i-1]
        icim1 = index[cim1]
        scim1 = sorted_sum[cim1]
        color = sscolors[n_change_inds-i-1]
    else:
        icim1 = 0
        scim1 = 0
        color = [0.2,0.2,0.2]
    cum_ax.plot([icim1, ici], [scim1, sci], lw=4, c=color)
    cum_ax.scatter([ici], [sci], s=150, zorder=50-i, c=[sscolors[n_change_inds-i-1]], edgecolor='w', linewidth=0.2)
cum_ax.set_xlabel('Proportion of primary infections (%)')
cum_ax.set_ylabel('Proportion of secondary infections (%)    ')
xmin = -2
ymin = -2
cum_ax.set_xlim([xmin, 100])
cum_ax.set_ylim([ymin, 102])
sc.boxoff(ax=cum_ax)

# Draw horizontal lines and annotations
ancol1 = [0.2, 0.2, 0.2]
ancol2 = sscolors[0]
ancol3 = sscolors[6]

i01 = sc.findinds(sorted_sum==100)[-1]
i20 = sc.findinds(sorted_sum>=80)[-1]
i50 = sc.findinds(sorted_sum>=50)[-1]
cum_ax.plot([index[i01], index[i01], 100], [ymin, 100, 100], '--', lw=2, c=ancol1)
cum_ax.plot([xmin, index[i20], index[i20]], [80, 80, ymin], '--', lw=2, c=ancol2)
cum_ax.plot([xmin, index[i50], index[i50]], [50, 50, ymin], '--', lw=2, c=ancol3)

bbox = dict(facecolor='w', alpha=0.9, lw=0)
cum_ax.text(44, 85, f'{100-index[i01]:0.0f}% of infections\ndo not transmit', c=ancol1, bbox=bbox)
cum_ax.text(25, 50, f'{index[i20]:0.0f}% of infections cause\n80% of transmissions', c=ancol2, bbox=bbox)
cum_ax.text(14, 15, f'{index[i50]:0.0f}% of infections cause\n50% of transmissions', c=ancol3, bbox=bbox)


#%% Fig. 2D -- histogram by date of symptom onset

# Calculate
asymp_count = 0
symp_counts = {}
minind = -5
maxind = 15
for _, target_ind in tt.transmissions:
    dd = tt.detailed[target_ind]
    date = dd['date']
    delta = sim.rescale_vec[date] # Increment counts by this much
    if dd['s']:
        if tt.detailed[dd['source']]['date'] <= date: # Skip dynamical scaling reinfections
            sdate = dd['s']['date_symptomatic']
            if np.isnan(sdate):
                asymp_count += delta
            else:
                ind = int(date - sdate)
                if ind not in symp_counts:
                    symp_counts[ind] = 0
                symp_counts[ind] += delta

# Convert to an array
days = list(symp_counts.keys())
xax = np.arange(minind-1, maxind+1)
counts = np.zeros(len(xax))
for i,val in symp_counts.items():
    if i<minind:
        ind = 0
    elif i>maxind:
        ind = -1
    else:
        ind = sc.findinds(xax==i)[0]
    counts[ind] += val

# Plot
total_count = asymp_count + counts.sum()
counts = counts/total_count*100
presymp = sc.findinds(xax<=0)[-1]
colors = ['#eed15b', '#ee943a', '#c3211a']

asymp_frac = asymp_count/total_count*100
pre_frac = counts[:presymp].sum()
symp_frac = counts[presymp:].sum()
symp_ax.bar(xax[0]-2, asymp_frac, label='Asymptomatic', color=colors[0])
symp_ax.bar(xax[:presymp], counts[:presymp], label='Presymptomatic', color=colors[1])
symp_ax.bar(xax[presymp:], counts[presymp:], label='Symptomatic', color=colors[2])
symp_ax.set_xlabel('Days since symptom onset')
symp_ax.set_ylabel('Proportion of transmissions (%)')
symp_ax.set_xticks([minind-3, 0, 5, 10, maxind])
symp_ax.set_xticklabels(['Asymp.', '0', '5', '10', f'>{maxind}'])
sc.boxoff(ax=symp_ax)

spie_ax = pl.axes([sympx+0.05, 0.20, 0.2, 0.2])
labels = [f'Asymp-\ntomatic\n{asymp_frac:0.0f}%', f'Presymp-\ntomatic\n{pre_frac:0.0f}%', f'Symp-\ntomatic\n{symp_frac:0.0f}%']
spie_ax.pie([asymp_frac, pre_frac, symp_frac], labels=labels, colors=colors, **pieargs)


#%% Tidy up

if do_save:
    cv.savefig(fig_path, dpi=150)

if do_show:
    pl.show()

sc.toc(T)

print('Done.')
