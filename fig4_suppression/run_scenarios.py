'''
Run the simulations for Fig. 4A (sensitivity scenarios) and Fig. 4B (reopening
sweeps). After running this script, run analyze_scenarios.py to parse into the
format used by plot_fig4.py.

To choose, change 'which'. You might also want to update the filenames with dates
in them.

Requires ~20,000 runs, each taking about 30 s, so running on an HPC or cluster
is required.
'''

import os
import sys
import shutil as sh
import numpy as np
import sciris as sc
import covasim as cv
import create_sim as cs

# Change this to make the runs for Fig. 4A or Fig. 4B, or call from the command line
which = 'a'
if len(sys.argv)>1:
    which = sys.argv[1]
randomize = True # False is not implemented yet

jsonfile = '../outputs/opt_merged_sep20_sg1.json'
doplot = 0
ncpus = 36
start_seed = 0
entry_offset = 0
npars_a = 6
npars_b = 3
reopenings = [0.6, 0.8, 1.0]

n_top_entries = 10 # Number of top entries to keep
resample = True # If true, resample seeds to only the top entries

if which == 'a':
    n_sweep_points = 51 # for a -- 51 points is ~3k sims
    msimfile = 'sensitivity_scenarios_sep28_final.msim'
else:
    ntrials = 12000 # For b
    msimfile = 'reopening_sweeps_sep28_final.msim'
json = sc.loadjson(jsonfile)


def get_f_a(v=None):
    ''' Get factors for panel A -- 1D sweeps '''
    if v is None:
        v = np.random.random(npars_a) # Must match length of dict
    f = sc.objdict(
        iqfactor  = v[0],
        testprob  = v[1],
        testqprob = v[2],
        trprob    = v[3],
        testdelay = v[4],
        trtime    = v[5],
        )
    return f


def get_f_b(v=None, index=None):
    ''' Get factor for panel B -- reopening '''
    if v is None:
        v = np.zeros(npars_b)
        if index is None:
            v[0] = 1.0
        else:
            v[0] = reopenings[index%len(reopenings)]
        v[1] = np.random.random() # Testing
        v[2] = np.random.random() # Tracing
    f = sc.objdict(
        reopen    = v[0],
        testprob  = v[1],
        trprob    = v[2],
        )
    return f


def rmnan(raw_val, calc_val, default_val):
    ''' If the raw input is is nan, return the default value; else, return the calculated value '''
    if np.isnan(raw_val):
        return default_val
    else:
        return calc_val


def scenconvert_a(f=None):
    ''' Convert scenarios for panel A '''
    if f is None:
        print('Using default f')
        f = get_f_a()
    f = sc.objdict(f)
    p = sc.objdict(
        symp_prob       = rmnan(f.testprob  , f.testprob*cs.d_calcs.max_test                               , cs.d_pars.symp_prob)  ,
        asymp_prob      = rmnan(f.testprob  , f.testprob*cs.d_calcs.max_test/cs.d_calcs.symp_asymp_ratio   , cs.d_pars.asymp_prob) ,
        symp_quar_prob  = rmnan(f.testqprob , f.testqprob                                                  , cs.d_pars.symp_quar_prob)  ,
        asymp_quar_prob = rmnan(f.testqprob , f.testqprob                                                  , cs.d_pars.asymp_quar_prob)  ,
        test_delay      = rmnan(f.testdelay , np.round(f.testdelay*cs.d_calcs.max_delay)                   , cs.d_pars.test_delay) ,
        trace_probs     = rmnan(f.trprob    , sc.mergedicts({k:f.trprob for k in 'hswl'}, {'c':0})         , cs.d_pars.trace_probs) ,
        trace_time      = rmnan(f.trtime    , {k:np.round(f.trtime*cs.d_calcs.max_delay) for k in 'hswcl'} , cs.d_pars.trace_time) ,
        i_factor        = rmnan(f.iqfactor  , {k:f.iqfactor for k in 'hswcl'}                              , cs.d_pars.i_factor) ,
        q_factor        = rmnan(f.iqfactor  , {k:min(1, f.iqfactor*2) for k in 'hswcl'}                    , cs.d_pars.q_factor) ,
        reopen          = cs.d_pars.reopen, # Static
    )
    return p


def scenconvert_b(f=None):
    ''' Convert scenarios for panel B '''
    if f is None:
        print('Using default f')
        f = get_f_b()
    f = sc.objdict(f)
    p = sc.objdict(
        symp_prob       = f.testprob*cs.d_calcs.max_test, # As above
        asymp_prob      = f.testprob*cs.d_calcs.max_test/cs.d_calcs.symp_asymp_ratio, # As above
        symp_quar_prob  = cs.d_pars.symp_quar_prob, # Static
        asymp_quar_prob = cs.d_pars.asymp_quar_prob, # Static
        test_delay      = cs.d_pars.test_delay, # Static
        trace_probs     = sc.mergedicts({k:f.trprob for k in 'hswl'}, {'c':0}),  # As above
        trace_time      = cs.d_pars.trace_time, # Static
        i_factor        = cs.d_pars.i_factor,  # Static
        q_factor        = cs.d_pars.q_factor,  # Static
        reopen          = f.reopen, # Different
    )
    return p


def run_scenario(index, v=None, which=None, doshrink=True, verbose=True):
    if verbose:
        print(f'Now running {index}:')
    cv.set_seed(index+start_seed)

    if which == 'a':
        if v is None:
            sweeps = construct1dsweeps()
            v = sweeps[index]
        f = get_f_a(v)
        scenpars = scenconvert_a(f)
    elif which == 'b':
        f = get_f_b(v, index=index)
        scenpars = scenconvert_b(f)

    if resample:
        eind = (index + entry_offset) % n_top_entries
    else:
        eind = index + entry_offset

    labelvals = [f'{v:0.2f}' for v in f.values()]
    labelstr = f'i={index} e={eind} f=' + ','.join(labelvals)

    runinfo = sc.objdict()
    runinfo.f = f
    runinfo.scenpars = scenpars
    runinfo.ind = index
    runinfo.eind = eind

    sim = cs.run_sim(index=eind, scenpars=scenpars, label=labelstr, runinfo=runinfo, verbose=verbose)

    # Try to save progress, but don't worry about it
    try:
        sc.savetext(f'progress/{index}.ind', str(f))
    except Exception as E:
        print(f'Could not save file {index}: {str(E)}')

    return sim


def construct1dsweeps():
    kernel = np.append(np.nan, np.linspace(0,1.0,n_sweep_points)) # np.sort(np.random.random(n_sweep_points))
    sweeps = []
    for p in range(npars_a):
        for k in kernel:
            for e in range(n_top_entries):
                sweep = np.nan*np.ones(npars_a)
                sweep[p] = k
                sweeps.append(sweep)

    return sweeps


if __name__ == '__main__':

    T = sc.tic()

    try:
        sh.rmtree('./progress/', ignore_errors=True)
        os.makedirs('./progress/', exist_ok=True)
    except Exception as E:
        print(f'Could not make progress folder: {E}')

    if which == 'a':
        sweeps = construct1dsweeps()
        ntrials = len(sweeps)

    sc.heading(f'Beginning run for type "{which}" for {ntrials} trials...')
    sc.timedsleep(3)

    sims = sc.parallelize(run_scenario, iterarg=np.arange(ntrials), ncpus=ncpus, kwargs={'which':which})

    msim = cv.MultiSim(sims)
    msim.save(msimfile)
    if doplot:
        msim.plot(max_sims=8, to_plot='overview', fig_args={'figsize':(38,19)})

    print('Done.')
    sc.toc(T)