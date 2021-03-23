'''
Run the simulations used in Fig. 5. The file saved by this script is used by
plot_fig5.py to generate the figure.

The simulations take about 5 min to run.
'''

import sciris as sc
import covasim as cv
import create_sim as cs

# Run configuration
indices    = range(5) # Number of replicates
ncpus      = min(len(indices), 5)
do_plot    = 0
do_save    = 1
from_cache = 1
msimsfile  = 'fig5.msims'
T = sc.tic()

# Set base parameters
base = {
        'trace_time'      : {k:1   for k in 'hswcl'},
        'i_factor'        : {k:0.1 for k in 'hswcl'}, # Isolation factor
        'q_factor'        : {k:0.2 for k in 'hswcl'}, # Quarantine factor
        'reopen'          : [0.80, 0.60], # Levels of reopening
}

# Set the three scenarios -- low TTQ, actual levels of TTQ, and high levels of TTQ
scens = dict(
    low = {
        'n_tests'         : 1500,
        'test_delay'      : 3,
        'trace_probs'     : sc.mergedicts({k:0.0 for k in 'scl'}, {'h': 0.2}, {'w':0.05}),
        'which'           : 'low',
    },
    actual = {
        'test_delay'      : 1,
        'trace_probs'     : sc.mergedicts({k:0.6 for k in 'hsl'}, {'w':0.05}, {'c': 0.0}),
        'which'           : 'actual',
    },
    high = {
        'symp_prob'       : 2.5*0.08,
        'asymp_prob'      : 2.5*0.001,
        'symp_quar_prob'  : 1.0,
        'asymp_quar_prob' : 1.0,
        'test_delay'      : 1,
        'trace_probs'     : sc.mergedicts({k:0.66 for k in 'hswl'}, {'c': 0.0}),
        'which'           : 'high',
    },
)


if __name__ == '__main__':

    msims = sc.objdict()

    for which in ['actual', 'low', 'high']:

        scenpars = sc.objdict(sc.mergedicts(base, scens[which]))
        kwargs = dict(scenpars=scenpars, from_cache=from_cache, do_shrink=0)

        # Run the simulations
        if ncpus>1:
            sims = sc.parallelize(cs.run_sim, iterarg=indices, ncpus=ncpus, kwargs=kwargs)
        else:
            sims = [cs.run_sim(index, **kwargs) for index in indices]

        # Merge into a multisim and optionally plot
        msim = cv.MultiSim(sims=sims)
        msim.reduce()
        msims[which] = msim
        if do_plot:
            msim.plot(to_plot='overview', fig_args={'figsize':(38,19)})

    if do_save:
        print('Saving...')
        for msim in msims.values():
            for sim in msim.sims:
                sim.shrink()
        cv.save(msimsfile, msims)

    sc.toc(T)
    print('Done.')