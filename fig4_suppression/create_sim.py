'''
Create the calibrated sim for the King County results. Used by run_scenarios.py. Very
similar to create_sim.py in the other folders, with the exception of the scenario definitions.

Note: this script relies the population files created by cache_population.py in the
Fig. 1 folder, as well as several CSV data files.
'''

import os
import numpy as np
import pandas as pd
import sciris as sc
import covasim as cv


# Define the input files
inputs         = '../inputs'
epi_data_file  = f'{inputs}/20200614chop5_KingCounty_Covasim.csv'
age_data_file  = f'{inputs}/20200614chop5_KingCounty_AgeHist.csv'
safegraph_file = f'{inputs}/KC_weeklyinteractions_20200811_trim.csv'
popfile_stem   = f'{inputs}/kc_rnr_seed'
jsonfile       = '../outputs/opt_merged_sep20_sg1.json' # Calibrated runs
json = sc.loadjson(jsonfile)

# Define default values for calculating the scenarios
d_calcs = sc.objdict()
d_calcs.routine_test     = 0.08 # Baseline daily symptomatic testing
d_calcs.testing_factor   = 2.7 # Scale-up compared to current testing
d_calcs.symp_asymp_ratio = 80 # How much less asymptomatic vs symptomatic testing there is
d_calcs.max_test         = 0.5 # Maximum (routine) testing probability to consider
d_calcs.max_delay        = 7 # Maximum testing or tracing delays to consider
d_calcs.quar_test        = 0.90 # Baseline quarantine testing
d_calcs.iso_factor       = 0.1 # How much people isolate
d_calcs.iq_ratio         = 2 # How much less people quarantine than isolate
d_calcs.tr_prob          = 0.70 # Tracing probability
d_calcs.t_delay          = 1 # Time for test
d_calcs.tr_delay         = 0 # Time for trace
d_calcs.reopen           = 1.0 # Amount by which to reopen to

# Define actual scenario parameters based on d_calcs values
d_pars = sc.objdict()
d_pars.symp_prob        = d_calcs.routine_test*d_calcs.testing_factor # 0.2 by default
d_pars.asymp_prob       = d_calcs.routine_test*d_calcs.testing_factor/d_calcs.symp_asymp_ratio # 0.0025 by default
d_pars.symp_quar_prob   = d_calcs.quar_test
d_pars.asymp_quar_prob  = d_calcs.quar_test
d_pars.test_delay       = d_calcs.t_delay
d_pars.trace_probs      = sc.mergedicts({k:d_calcs.tr_prob for k in 'hswl'}, {'c':0})
d_pars.trace_time       = {k:d_calcs.tr_delay for k in 'hswcl'}
d_pars.i_factor         = {k:d_calcs.iso_factor for k in 'hswcl'} # Isolation
d_pars.q_factor         = {k:min(1, d_calcs.iso_factor*d_calcs.iq_ratio) for k in 'hswcl'} # Quarantine
d_pars.reopen           = d_calcs.reopen


def load_pars(index):
    ''' Load the parameters from JSON '''
    entry = json[index]
    pars = sc.objdict(entry['pars'])
    pars.rand_seed = int(entry['index'])
    print(f'Loading parameters from trial {index}, mismatch {entry["mismatch"]}...')
    sc.pp(pars)
    return pars


# Generate the population filename
def get_popfile(pars):
    ''' Load the population files '''
    n_popfiles = 5
    popfile = popfile_stem + str(pars['rand_seed']%n_popfiles) + '.ppl'

    # Check that the population file exists
    if not os.path.exists(popfile):
        errormsg = f'WARNING: could not find population file {popfile}! Please regenerate first'
        raise FileNotFoundError(errormsg)

    return popfile


def check_contacts(sim, check=False, verbose=True):
    ''' Store the number of contacts in the sim '''
    if check:
        contacts = {}
        for lkey in ['h','w','s','c']:
            contacts[lkey] = len(sim.people.contacts[lkey])
        if not hasattr(sim, 'n_contacts'):
            sim.n_contacts = sc.odict()
        sim.n_contacts[sim.date(sim.t)] = contacts
        if verbose:
            print(f'>>> On day {sim.t}, there were {contacts} contacts')
    return


def make_safegraph(sim):
    ''' Create interventions representing SafeGraph data '''

    # Load data.values
    fn = safegraph_file
    df = pd.read_csv(fn)
    week = df['week']
    s = df['p.tot.schools'].values
    w = df['p.tot.no.schools'].values
    c = sc.dcp(w) # Not different enough to warrant different values

    # Do processing
    days = sim.day(week.values.tolist())
    last_day = days[-1]+1
    i_days = np.arange(days[0], last_day)
    s = np.interp(i_days, days, s)
    w = np.interp(i_days, days, w)
    c = np.interp(i_days, days, c)
    days = i_days

    # Create interventions
    interventions = [
        cv.clip_edges(days=days, changes=s, layers='s', label='clip_s'),
        cv.clip_edges(days=days, changes=w, layers='w', label='clip_w'),
        cv.clip_edges(days=days, changes=c, layers='c', label='clip_c'),
        ]

    return interventions


def remove_ltcf_community(sim, debug=False):
    ''' Ensure LTCF residents don't have community contacts '''
    over_65 = sc.findinds(sim.people.age>65)
    llayer = sim.people.contacts['l']
    clayer = sim.people.contacts['c']
    in_ltcf = np.union1d(llayer['p1'], llayer['p2'])
    over_65_in_ltcf = np.intersect1d(over_65, in_ltcf)
    p1inds = sc.findinds(np.isin(clayer['p1'], over_65_in_ltcf))
    p2inds = sc.findinds(np.isin(clayer['p2'], over_65_in_ltcf))
    popinds = np.union1d(p1inds, p2inds)
    clayer.pop_inds(popinds)
    return clayer


def test_num_subtarg(sim, sev=100.0, u20=0.5):
    ''' Subtarget severe people with more testing, and young people with less '''
    sev_inds = sim.people.true('severe')
    u20_inds = sc.findinds((sim.people.age<20) * (~sim.people.severe)) # People who are under 20 and severe test as if they're severe; * is element-wise "and"
    u20_vals = u20*np.ones(len(u20_inds))
    sev_vals = sev*np.ones(len(sev_inds))
    inds = np.concatenate([u20_inds, sev_inds])
    vals = np.concatenate([u20_vals, sev_vals])
    return {'inds':inds, 'vals':vals}


def create_sim(eind, verbose=True):
    ''' Create a single simulation for further use '''

    p = load_pars(eind)

    # Basic parameters and sim creation
    pars = {  'pop_size'      : 225e3,
              'pop_scale'     : 10,
              'pop_type'      : 'synthpops',
              'pop_infected'  : 300,
              'beta'          : p.beta,
              'start_day'     : '2020-01-27',
              'end_day'       : '2020-08-31',
              'rescale'       : True,
              'rescale_factor': 1.1,
              'verbose'       : p.get('verbose', 0.01*verbose),
              'rand_seed'     : int(p.rand_seed),
              'analyzers'     : cv.daily_stats(days=[]),
              'beta_layer'    : dict(h=3.0, s=0.6, w=0.6, c=0.3, l=1.5),
            }

    # Create and initialize the sim
    if pars['verbose']:
        print(f'Creating sim! seed={p.rand_seed}')
    sim = cv.Sim(pars, label=f'base_sim_{p.rand_seed}', popfile=get_popfile(pars), load_pop=True, datafile=epi_data_file) # Create this here so can be used for test numbers etc.

    sim.sceninfo = sc.objdict()
    sim.sceninfo.calib_end  = '2020-05-31'
    sim.sceninfo.scen_start = '2020-06-01'

    # Define testing interventions -- 97% sensitivity from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7177629/
    test_kwargs = dict(daily_tests=sim.data['new_tests'], test_delay=2, sensitivity=0.97, subtarget=test_num_subtarg)
    tn = cv.test_num(symp_test=p.tn, start_day='2020-01-27', end_day=sim.sceninfo.calib_end, **test_kwargs, label='tn')
    interventions = [tn]

    # Define beta interventions
    hwc_days  = ['2020-02-24', '2020-03-23'] # Change date here, 04-27 or 05-04
    hwc_days = sim.day(hwc_days)
    b_wc_ch   = [1.0, p.bc_wc1] # Simple interpoation
    b_h_ch    = [1.0, 1.1] # Optional household

    all_b_days = np.arange(hwc_days[0], hwc_days[-1]+1) # Full time series
    all_ch_wc = np.interp(all_b_days, hwc_days, b_wc_ch) # Linearly interpolate
    all_ch_h = np.interp(all_b_days, hwc_days, b_h_ch) # Linearly interpolate
    interventions += [cv.change_beta(days=all_b_days, changes=all_ch_h, layers='h', label='beta_h')]
    lkeys = ['w','c','s']
    for lkey in lkeys: # Assume schools also use masks, etc. so have same non-movement beta change
        cb = cv.change_beta(days=all_b_days, changes=all_ch_wc, layers=lkey, label=f'beta_{lkey}')
        interventions += [cb]

    # LTCF beta change
    b_l_days = ['2020-02-24', '2020-03-23']
    b_l_days = np.arange(sim.day(b_l_days[0]), sim.day(b_l_days[1]))
    b_l_ch   = np.linspace(1.0, p.bc_lf, len(b_l_days))
    interventions += [cv.change_beta(days=b_l_days, changes=b_l_ch, layers='l', label='beta_l')]
    sim.people.contacts['c'] = remove_ltcf_community(sim) # Remove community contacts from LTCF

    # SafeGraph intervention & tidy up
    interventions += make_safegraph(sim)
    sim['interventions'] = interventions

    # Don't show interventions in plots, there are too many
    for interv in sim['interventions']:
        interv.do_plot = False

    # These are copied from parameters.py -- modified to capture change in work status at age 65
    sim['prognoses']['age_cutoffs'] = np.array([0,      10,      20,      30,      40,      50,      65,      70,      80,      90]) # Age cutoffs (upper limits)

    return sim


def get_interv(sim, label, ):
    ''' Return an intervention based on a label '''
    match = []
    nomatch = []
    for interv in sim['interventions']:
        if interv.label == label:
            match += [interv]
        else:
            nomatch += [interv.label]
    if len(match) == 1:
        return match[0]
    elif len(match) > 1:
        raise Exception(f'Multiple matches for "{label}" found ({len(match)}!')
    else:
        raise Exception(f'Label "{label}" not found amongst {nomatch}')
    return


def modify_sim(sim, scenpars, label=None, runinfo=None):
    ''' Modify the simulation for the scenarios '''

    print(f'  Note: modifying simulation {label} at day={sim.t}, date={sim.date(sim.t)}, scenpars:\n{scenpars}')

    # Do reopening: modify clip_edges
    last_calib_day = sim.day(sim.sceninfo.calib_end)
    first_scen_day = sim.day(sim.sceninfo.scen_start)
    for ilabel in ['clip_w', 'clip_c']:
        interv = get_interv(sim,ilabel)
        valid_days = sc.findinds(interv.days<=last_calib_day)
        interv.days = np.append(interv.days[valid_days], first_scen_day) # NB, repr of intervention will be wrong with direct modification!
        interv.changes = np.append(interv.changes[valid_days], scenpars['reopen'])

    # Change iso_factor and quar_factor
    sim.pars['iso_factor']  = sc.dcp(scenpars['i_factor'])
    sim.pars['quar_factor'] = sc.dcp(scenpars['q_factor'])

    # Implement testing & tracing interventions
    tppars = {k:scenpars[k] for k in ['symp_prob', 'asymp_prob', 'symp_quar_prob', 'asymp_quar_prob', 'test_delay']}
    ctpars = {k:scenpars[k] for k in ['trace_probs', 'trace_time']}
    tp = cv.test_prob(start_day=first_scen_day, quar_policy=[0], **tppars)
    ct = cv.contact_tracing(start_day=first_scen_day, label='contact_tracing', **ctpars)
    sim['interventions'] += [tp, ct]

    # Final tidying
    sim.label = label
    sim.runinfo = runinfo
    sim.sceninfo.scenpars = scenpars

    return sim


def run_sim(index, scenpars=None, label=None, runinfo=None, do_shrink=True, from_cache=True, verbose=True):
    ''' Load, modify, and run the simulation '''

    if scenpars is None:
        print('WARNING, loading default parameters!')
        scenpars = sc.dcp(d_pars)

    if label is None:
        label = f'full_sim_trial{index}_nolabel'

    # Run sim or load up until scenarios start
    sim_loaded = 0
    filename = f'./simcache/sim_trial{index}.sim'
    if from_cache and os.path.exists(filename):
        tries = 0
        while not sim_loaded and tries<3:
            sc.timedsleep(0.5)
            try:
                print(f'run_sim(): Loading sim from cache ({filename})...')
                sim = cv.load(filename)
                assert isinstance(sim, cv.Sim) # Sometimes unpickling fails
                sim_loaded = 1
            except Exception as E:
                print(f'Loading failed on try {tries}! {E}')
                string = '\n\n'.join([f'Index {index}', filename, f'Try {tries}', str(E), sc.traceback()])
                fn = f'simcache_exception_index{index}_try{tries}.err'
                sc.savetext(fn, string)
                tries += 1

    if not sim_loaded:
        print(f'run_sim(): Creating new sim ({filename})...')
        sim = create_sim(index, verbose=verbose)
        sim.run(until=sim.sceninfo.calib_end)
        sim.save(filename, keep_people=True)

    # Modify the sim with these scenario values
    sim = modify_sim(sim, scenpars=scenpars, label=label, runinfo=runinfo)

    # Rerun till the end and optionally shrink
    sim.run(restore_pars=False) # So we can see the changes made
    if do_shrink:
        sim.shrink()

    return sim


if __name__ == '__main__':

    T = sc.tic()

    indices  = [0]
    ncpus    = 1
    do_plot  = 1
    scenpars = None

    # Settings
    if ncpus>1:
        sims = sc.parallelize(run_sim, iterarg=indices, ncpus=ncpus, kwargs={'scenpars':scenpars, 'do_shrink':0})
    else:
        sims = [run_sim(index, scenpars=scenpars) for index in indices]

    msim = cv.MultiSim(sims=sims)
    if do_plot:
        get_interv(msim.sims[0], 'contact_tracing').do_plot = True # Turn on plotting for this intervention
        msim.plot(to_plot='overview', fig_args={'figsize':(38,19)}, plot_sims=True)

    sc.toc(T)
