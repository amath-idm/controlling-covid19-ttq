'''
Create the calibrated sim for the King County results. Used by run_fig1a.py as
well as run_optimization.py.

Note: this script relies the population files created by cache_population.py, as
well as several CSV data files.
'''

# Standard packages
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


# Generate the population filename
def get_popfile(pars):
    n_popfiles = 5
    popfile = popfile_stem + str(pars['rand_seed']%n_popfiles) + '.ppl'

    # Check that the population file exists
    if not os.path.exists(popfile):
        errormsg = f'WARNING: could not find population file {popfile}! Please regenerate first'
        raise FileNotFoundError(errormsg)

    return popfile


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
    sim.intervention_info.ce = sc.objdict({'days':days, 'changes':w})

    return interventions


def define_pars(which='best', use_safegraph=True):
    ''' Define the parameter best guesses and bounds '''

    pardata = dict(
        # name          best       low       high
        beta         = [  0.0146,  0.013  ,  0.016],
        bc_lf        = [  0.20  ,  0.05   ,  0.40],
        tn           = [ 25.0   , 10.0    , 60.0],
    )
    if use_safegraph:
        pardata.update(dict(
            bc_wc1   = [  0.85  ,  0.50   ,   1.00],
        ))
    else:
        pardata.update(dict(
            bc_wc1   = [  0.20  ,  0.10   ,   0.60],
        ))

    output = {}
    for key,arr in pardata.items():
        if which == 'best':
            output[key] = arr[0]
        elif which == 'bounds':
            output[key] = arr[1:3]

    return output


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


def create_sim(pars=None, use_safegraph=True, label=None, show_intervs=False):
    ''' Create a single simulation for further use '''

    p = sc.objdict(sc.mergedicts(define_pars(which='best', use_safegraph=use_safegraph), pars))
    if 'rand_seed' not in p:
        seed = 1
        print(f'Note, could not find random seed in {pars}! Setting to {seed}')
        p['rand_seed'] = seed # Ensure this exists

    # Basic parameters and sim creation
    pars = {  'pop_size'      : 225e3,
              'pop_scale'     : 10,
              'pop_type'      : 'synthpops',
              'pop_infected'  : 300,
              'beta'          : p.beta,
              'start_day'     : '2020-01-27',
              'end_day'       : '2020-06-08',
              'rescale'       : True,
              'rescale_factor': 1.1,
              'verbose'       : p.get('verbose', 0.01),
              'rand_seed'     : int(p.rand_seed),
              'analyzers'     : cv.age_histogram(datafile=age_data_file),
              'beta_layer'    : dict(h=3.0, s=0.6, w=0.6, c=0.3, l=1.5),
            }

    # Create and initialize the sim
    if pars['verbose']:
        print(f'Creating sim! safegraph={use_safegraph}, seed={p.rand_seed}')
    sim = cv.Sim(pars, label=label, popfile=get_popfile(pars), load_pop=True, datafile=epi_data_file) # Create this here so can be used for test numbers etc.

    # Define testing interventions -- 97% sensitivity from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7177629/
    test_kwargs = dict(daily_tests=sim.data['new_tests'], test_delay=2, sensitivity=0.97, subtarget=test_num_subtarg)
    tn = cv.test_num(symp_test=p.tn, start_day='2020-01-27', end_day=None, **test_kwargs, label='tn')
    interventions = [tn]

    # Define beta interventions
    sim.intervention_info = sc.objdict()
    hwc_days  = ['2020-02-24', '2020-03-23', '2020-05-31'] # Change date here, 04-27 or 05-04
    hwc_days = sim.day(hwc_days)
    b_wc_ch   = [1.0, p.bc_wc1, p.get('bc_wc2', p.bc_wc1)] # To allow either one or two beta change parameters
    b_h_ch    = [1.0, 1.1, 1.1] # Optional household

    all_b_days = np.arange(hwc_days[0], hwc_days[-1]+1) # Full time series
    all_ch_wc = np.interp(all_b_days, hwc_days, b_wc_ch) # Linearly interpolate
    all_ch_h = np.interp(all_b_days, hwc_days, b_h_ch) # Linearly interpolate
    interventions += [cv.change_beta(days=all_b_days, changes=all_ch_h, layers='h', label='beta_h')]
    lkeys = ['w','c','s'] if use_safegraph else ['w','c'] # Skip schools if not using SafeGraph
    for lkey in lkeys: # Assume schools also use masks, etc. so have same non-movement beta change
        cb = cv.change_beta(days=all_b_days, changes=all_ch_wc, layers=lkey, label=f'beta_{lkey}')
        interventions += [cb]
        sim.intervention_info.bc = sc.objdict({'days':all_b_days, 'changes':all_ch_wc}) # Store for plotting later

    # LTCF beta change
    b_l_days = ['2020-02-24', '2020-03-23']
    b_l_days = np.arange(sim.day(b_l_days[0]), sim.day(b_l_days[1]))
    b_l_ch   = np.linspace(1.0, p.bc_lf, len(b_l_days))
    interventions += [cv.change_beta(days=b_l_days, changes=b_l_ch, layers='l', label='beta_l')]
    sim.people.contacts['c'] = remove_ltcf_community(sim) # Remove community contacts from LTCF

    # SafeGraph intervention & tidy up
    if use_safegraph:
        interventions += make_safegraph(sim)
    else:
        interventions += [cv.clip_edges(days='2020-03-12', changes=0.1, layers='s', label='clip_s')]
    sim['interventions'] = interventions

    # Don't show interventions in plots, there are too many
    for interv in sim['interventions']:
        interv.do_plot = False

    # These are copied from parameters.py -- modified to capture change in work status at age 65
    sim['prognoses']['age_cutoffs'] = np.array([0,      10,      20,      30,      40,      50,      65,      70,      80,      90]) # Age cutoffs (upper limits)

    return sim


def run_sim(pars=None, interactive=False, sim=None, use_safegraph=True, do_plot=True): # Difference
    ''' Create and run a simulation from a given set of parameters '''

    # Create and run the sim
    if sim is None:
        sim = create_sim(pars=pars, use_safegraph=use_safegraph)
    if not sim.results_ready:
        sim.run()

    # Get the age histogram and compute the fit
    ageh = sim['analyzers'][0]
    agehd = ageh.data
    agehh = ageh.hists[-1]
    custom = sc.objdict(dict(
        age_tests     = dict(data=agehd['cum_tests'][:],     sim=agehh['tested'],    weight=0.0),
        age_diagnoses = dict(data=agehd['cum_diagnoses'][:], sim=agehh['diagnosed'], weight=0.2),
        age_deaths    = dict(data=agehd['cum_deaths'][:],    sim=agehh['dead'],      weight=1),
        age_yield     = dict(data=agehd['cum_diagnoses'][:]/agehd['cum_tests'].values, sim=agehh['diagnosed']/agehh['tested'], weight=1.0),
    ))

    # Add rolling
    def roll(vec):
        rolled = pd.Series(vec).rolling(window=window).mean().to_numpy()
        return rolled[~np.isnan(rolled)]

    offset = 0 # Difference between start data for the data and start date for the sim
    window = 7 # A week
    data_roll_diag  = roll(sim.data['new_diagnoses'])
    data_roll_death = roll(sim.data['new_deaths'])
    sim_roll_diag   = roll(sim.results['new_diagnoses'].values[offset:])
    sim_roll_death  = roll(sim.results['new_deaths'].values[offset:])
    custom.rolling_diag = dict(data=data_roll_diag,  sim=sim_roll_diag,  weight=0.5)
    custom.rolling_death = dict(data=data_roll_death, sim=sim_roll_death, weight=0.5)

    ramp_weight = 1 + np.linspace(0, 1, sim.npts)
    weights = {'cum_diagnoses':ramp_weight, 'cum_deaths':ramp_weight, 'new_diagnoses':0.01, 'new_deaths':0.01}
    fit = sim.compute_fit(custom=custom, keys=['cum_diagnoses', 'cum_deaths', 'new_diagnoses', 'new_deaths'], weights=weights)

    # Handle output
    if interactive:
        if do_plot:
            sim.plot()
            fit.plot()
        return sim, fit
    else:
        return fit.mismatch


if __name__ == '__main__':

    T = sc.tic()

    use_multisim  = 1
    use_safegraph = 1

    # Settings
    reps = 5 # Set multiple runs to average over likelihood
    base_sim = create_sim(use_safegraph=use_safegraph)

    # Plot calibration
    if use_multisim:
        msim = cv.MultiSim(base_sim, n_runs=reps)
        msim.run(reseed=True, noise=0.0, keep_people=True, par_args={'ncpus':5})
        sims = msim.sims
        msim.plot(to_plot='overview', plot_sims=True)
        msim.reduce()
        sim = msim.base_sim
    else:
        sim = base_sim
        sim.run()
        sims = [sim]

    # Do plotting
    sim, fit = run_sim(sim=sim, use_safegraph=use_safegraph, interactive=True)

    sc.toc(T)
