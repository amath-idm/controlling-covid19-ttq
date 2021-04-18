'''
Run the simulations used in Fig. 3A-C. The file saved by this script is used by
plot_fig3.py to generate the figure.

The simulations only take a few seconds to run.
'''

import covasim as cv
import sciris as sc

# cv.check_save_version('2.0.0', filename='fig3a_run.gitinfo', die=True) # Must be this version or else transtrees are wrong

do_save = True
tt_filename = 'fig3_transtrees.obj'

iday = 20
p = sc.objdict(
    pop_size = 100,
    n_days = 125,
    pop_type = 'hybrid',
    pop_infected = 1,
    rand_seed = 8, # Clearest illustration of dynamics (in many cases the epidemic dies out immmediately)
    verbose = 0,
    iso_factor = dict(h=0.3, s=0.0, w=0.0, c=0.1),
    quar_factor = dict(h=0.8, s=0.0, w=0.0, c=0.3),
)

# Run sims
s = sc.objdict()
t = sc.objdict()
for k in ['none', 'test', 'trace']:
    interventions = []
    if k in ['test', 'trace']:
        interventions += [cv.test_prob(start_day=iday, symp_prob=0.15, symp_quar_prob=1.0, asymp_quar_prob=1.0, quar_policy='start')]
    if k in ['trace']:
        interventions += [cv.contact_tracing(start_day=iday, trace_probs=dict(h=0.7, s=0.1, w=0.1, c=0.0))]

    s[k] = cv.Sim(p, interventions=interventions, version='2.0.0') # Version 2.0 parameters are what were used in the manuscript
    s[k].run()
    t[k] = s[k].make_transtree()

# Save
cv.save(tt_filename, [p, s, t])