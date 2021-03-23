
'''
Run the simulation used in Fig. 2. The file saved by this script is used by
plot_fig2.py to generate the figure. Requires a JSON file produced by processing
the results of run_optimization.py (included in the repository). Also requires a
population file produced by cache_populations.py.

The simulation takes about 5 min to run.
'''

import sciris as sc
import covasim as cv
import create_sim as cs

# Settings
index    = 0
do_plot  = 1
do_save  = 1
verbose  = 1
sg       = 1 # Use the calibration with SafeGraph data
jsonfile = f'../outputs/opt_merged_sep20_sg{sg}.json'
pplfile  = '../inputs/kc_big_seed0.ppl'
simfile  = 'fig2.sim'


print(f'Original parameters for index {index}:')
json = sc.loadjson(jsonfile)
entry = json[index]
sc.pp(entry)
pars = entry['pars']
pars['rand_seed'] = int(entry['index'])+0 # To run with a different seed

# Adjust for the larger population
pars['pop_size'] = 2.25e6 # Reset for the full population
pars['beta'] *= 0.99 # Adjust the calibration slightly for the larger population
pars['verbose'] = verbose


print('Loading population file...')
with sc.Timer():
    people = cv.load(pplfile)


print('Creating sim...')
with sc.Timer():
    sim = cs.create_sim(pars=pars, people=people, use_safegraph=sg)


print('Running sim...')
with sc.Timer():
    sim, fit = cs.run_sim(sim=sim, use_safegraph=sg, interactive=do_plot)


print('Analyzing sim...')
with sc.Timer():
    sim.make_transtree(output=False)


if do_save:
    print('Saving sim...')
    with sc.Timer():
        sim.save(simfile, keep_people=True)

print('Done.')
