'''
Run the simulations used in Fig. 1. The file saved by this script is used by
plot_fig1.py to generate the figure. Requires a JSON file produced by processing
the results of run_optimization.py (included in the repository).

Each simulation takes about 30 s to run. By default 200 simulations are run, so
it is strongly recommended that you run on an HPC (the default set below is for
36 cores).
'''

import numpy as np
import sciris as sc
import covasim as cv
import create_sim as cs

T = sc.tic()

use_safegraph = 1
do_save = 1
do_plot = 0
parallelize = 1
n_sims = 200
jsonfile = '../outputs/opt_merged_sep20_sg1.json' # Calibration results
simsfile = 'fig1a.sims' # Output file

json = sc.loadjson(jsonfile)

total_sample_size = len(json)-1
indices = np.linspace(0, total_sample_size, n_sims, dtype=int)

def make_sim(index):
    print(f'Now creating {index}:')
    entry = json[index]
    sc.pp(entry)
    pars = entry['pars']
    pars['rand_seed'] = int(entry['index'])+0 # To run with a different seed
    sim = cs.create_sim(pars=pars, use_safegraph=use_safegraph)
    sim.label = f'Trial {index}'
    sim.jsonpars = entry
    return sim


if parallelize:
    sims = sc.parallelize(make_sim, indices)
else:
    sims = []
    for index in indices:
        sims.append(make_sim(index))


msim = cv.MultiSim(sims)
msim.run(par_args={'ncpus':36})
if do_save:
    cv.save(simsfile, msim.sims)
if do_plot:
    msim.plot()

sc.toc(T)
print('Done.')
