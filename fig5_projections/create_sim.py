'''
Create the calibrated sim for the King County results. Used by run_scenarios.py. Very
similar to create_sim.py in the other folders, but projects further into the future.

Note: this script relies the population files created by cache_population.py in the
Fig. 1 folder, as well as several CSV data files.
'''

import os
import numpy as np
import pandas as pd
import sciris as sc
import covasim as cv


# Define the input files
epi_data_file  = './kc_data/20200614chop5_KingCounty_Covasim_extended.xlsx'
safegraph_file = './kc_data/KC_weeklyinteractions_20200811_extended.xlsx'
popfile_stem   = '../inputs/kc_rnr_seed'
jsonfile       = '../outputs/opt_merged_sep20_sg1.json' # Calibrated runs
json = sc.loadjson(jsonfile)

# Define default values for calculating the scenarios
d_calcs = sc.objdict()
d_calcs.routine_test     = 0.08 # Baseline daily symptomatic testing
d_calcs.testing_factor   = 2.5 # Scale-up compared to current testing
d_calcs.symp_asymp_ratio = 80 # How much less asymptomatic vs symptomatic testing there is
d_calcs.max_test         = 0.5 # Maximum (routine) testing probability to consider
d_calcs.max_delay        = 7 # Maximum testing or tracing delays to consider
d_calcs.quar_test        = 0.90 # Baseline quarantine testing
d_calcs.iso_factor       = 0.1 # How much people isolate
d_calcs.iq_ratio         = 2 # How much less people quarantine than isolate
d_calcs.tr_prob          = 0.66 # Tracing probability
d_calcs.t_delay          = 1 # Time for test
d_calcs.tr_delay         = 2 # Time for trace
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
    df = pd.read_excel(fn)
    week = df['week']
    s = df['p.tot.schools'].values
    w = df['merged'].values
    c = sc.dcp(w) # Not different enough to warrant different values

    # Do processing
    days = sim.day([pd.Timestamp(v) for v in week.values])
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


def test_num_subtarg(sim, sev=100.0, u20=0.5, quar=1.0):
    ''' Subtarget severe people with more testing, and young people with less '''
    inds = np.arange(len(sim.people))
    vals = np.ones(len(sim.people))
    vals[sc.findinds((sim.people.age<20) * (~sim.people.severe))] *= u20 # People who are under 20 and severe test as if they're severe; * is element-wise "and"
    vals[sim.people.true('severe')] *= sev
    vals[sim.people.true('quarantined')] *= quar
    return {'inds':inds, 'vals':vals}


def create_sim(eind, verbose=True): # Difference
    ''' Create a single simulation for further use '''

    p = load_pars(eind) # Difference

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

    # Define testing interventions -- 97% sensitivity from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7177629/
    test_kwargs = dict(daily_tests=sim.data['new_tests'], test_delay=2, sensitivity=0.97, subtarget=test_num_subtarg)
    tn = cv.test_num(symp_test=p.tn, start_day='2020-01-27', end_day=None, **test_kwargs, label='tn')
    interventions = [tn]

    # Define beta interventions
    hwc_days  = ['2020-02-24', '2020-03-23'] # Change date here, 04-27 or 05-04; 'hwc' = household, work, community
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


def modify_sim(sim, scenpars, label=None, runinfo=None):
    ''' Modify the simulation for the scenarios '''

    print(f'  Note: modifying simulation {label} at day={sim.t}, date={sim.date(sim.t)}, scenpars:\n{scenpars}')

    sim.sceninfo.scen_start = '2020-06-01'
    sim.sceninfo.scale_down = '2020-07-01'
    last_calib_day = sim.day(sim.sceninfo.calib_end)
    first_scen_day = sim.day(sim.sceninfo.scen_start)
    scale_down_day = sim.day(sim.sceninfo.scale_down)
    n_days = 5
    scale_range = np.arange(n_days)
    ramp_up = np.linspace(0, 1, n_days)
    ramp_down = 1 - ramp_up

    for ilabel in ['clip_w', 'clip_c']:
        interv = sim.get_interventions(ilabel)
        valid_days = sc.findinds(interv.days<=last_calib_day)

        # Do reopening: modify clip_edges
        v0 = interv.changes[valid_days[-1]]
        v1 = scenpars['reopen'][0]
        v2 = scenpars['reopen'][1]
        scale_up   = ramp_down*v0 + ramp_up*v1
        scale_down = ramp_down*v1 + ramp_up*v2
        scale_days = np.append(scale_range+first_scen_day, scale_range+scale_down_day)
        scale_vals = np.append(scale_up, scale_down)

        interv.days = np.append(interv.days[valid_days], scale_days)
        interv.changes = np.append(interv.changes[valid_days], scale_vals)

    # Change iso_factor and quar_factor
    sim.pars['iso_factor']  = sc.dcp(scenpars['i_factor'])
    sim.pars['quar_factor'] = sc.dcp(scenpars['q_factor'])

    # Implement testing & tracing interventions
    ctpars = {k:scenpars[k] for k in ['trace_probs', 'trace_time']}
    tn = sim.get_interventions('tn')
    if scenpars['which'] == 'actual':
        tn.test_delay = scenpars.test_delay
        sim['interventions'] += [cv.contact_tracing(start_day=first_scen_day, label='contact_tracing', **ctpars)]
    elif scenpars['which'] == 'low':
        offset = 7 # So it's constant
        tn.daily_tests[first_scen_day-tn.start_day-offset:] = scenpars['n_tests']
        sim['interventions'] += [cv.contact_tracing(start_day=first_scen_day, label='contact_tracing', **ctpars)]
    elif scenpars['which'] == 'high':
        cv.set_seed(sim['rand_seed']) # Just in case since we use a random number
        tn.end_day = last_calib_day # End the test_num intervention
        n_weeks = 4 # Ramp up the test_prob and contact_tracing interventions
        days = first_scen_day + np.arange(n_weeks)*7
        tpramp = np.maximum(0, np.linspace(0.2, 0.7, n_weeks+1)[1:]+np.random.randn(n_weeks)*0.1) # Don't start from 0 and don't make it a perfectly smooth ramp
        ctramp = np.linspace(0.0, 0.7, n_weeks+1)[1:] # Don't start from 0
        tppars = {k:scenpars[k] for k in ['symp_prob', 'asymp_prob', 'symp_quar_prob', 'asymp_quar_prob', 'test_delay']}
        tplist = []
        ctlist = []
        for w in range(n_weeks):
            tp = cv.test_prob(quar_policy=[0], **tppars, label=f'tp_{w}')
            ct = cv.contact_tracing(label=f'ct_{w}', **ctpars)
            tp.symp_prob *= tpramp[w]
            tp.asymp_prob *= tpramp[w]
            ct.trace_probs = sc.dcp(ct.trace_probs)
            for key in ct.trace_probs.keys(): ct.trace_probs[key] *= ctramp[w]
            tplist.append(tp)
            ctlist.append(ct)
            print(f'I AM {w} {ct.trace_probs}')
        tpseq = cv.sequence(days=days, interventions=tplist)
        ctseq = cv.sequence(days=days, interventions=ctlist)
        tpseq.initialize(sim)
        ctseq.initialize(sim)
        sim['interventions'] += [tpseq, ctseq]

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
    filename = f'./simcache/projection_trial{index}.sim'
    if from_cache and os.path.exists(filename):
        tries = 0
        while not sim_loaded and tries<3:
            try:
                print(f'run_sim(): Loading sim from cache ({filename})...')
                sim = cv.load(filename)
                assert isinstance(sim, cv.Sim) # Sometimes unpickling fails, but can reload locally
                tn = sim.get_interventions('tn')
                tn.subtarget = test_num_subtarg
                sim_loaded = 1
            except Exception as E:
                print(f'Loading failed on try {tries}! {E}')
                string = '\n\n'.join([f'Index {index}', filename, f'Try {tries}', str(E), sc.traceback()])
                fn = f'simcache_exception_index{index}_try{tries}.err'
                sc.savetext(fn, string)
                tries += 1
                sc.timedsleep(0.5)

    if not sim_loaded:
        print(f'run_sim(): Creating new sim (and saving to {filename})...')
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
    sc.toc(T)
