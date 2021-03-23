'''
Run the simulations used in Fig. 3E-F. The file saved by this script is used by
plot_fig3.py to generate the figure.

The simulations take about 5 min to run.
'''

import numpy as np
import sciris as sc
import covasim as cv

# Basic setup
do_save = 0
simsfile = 'fig3.sims'
T = sc.tic()
s = sc.objdict()
a = sc.objdict()


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
n_seeds = 10 # Number of replicates

# Define scenarios
sd = 10 # Start day
ed = 15 # End day
bc = 1.0 # Beta change
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


if __name__ == '__main__':

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

    print('Done.')
    sc.toc(T)