'''
Pre-generate the population. Takes about 10 min. Almost identical to
the version for Fig. 1, except 10x larger of a population (all of King
County).

Requires SynthPops, which must be installed from the repository:

https://github.com/institutefordiseasemodeling/synthpops
'''

import psutil
import sciris  as sc
import covasim as cv
import synthpops as sp
sp.config.set_nbrackets(20) # Essential for getting the age distribution right
sp.logger.setLevel('DEBUG') # Show additional information during population creation

# Settings
pop_size = 2.25e6 # 100% of the King County population
inputs   = '../inputs'
popfile_stem = f'{inputs}/kc_big_seed' # Stands for "King County big population random seed"

def cache_populations(seed=0, popfile=None):
    ''' Pre-generate the synthpops population '''

    pars = sc.objdict(
        pop_size = pop_size,
        pop_type = 'synthpops',
        rand_seed = seed,
    )

    if popfile is None:
        popfile = f'{popfile_stem}{pars.rand_seed}.ppl'

    T = sc.tic()
    print(f'Making "{popfile}"...')
    sim = cv.Sim(pars)
    cv.make_people(sim, popfile=popfile, save_pop=True, with_facilities=True, generate=True, layer_mapping={'LTCF':'l'})
    sc.toc(T)

    print('Done')
    return


if __name__ == '__main__':

    seeds = [0] # NB, each one takes at least 8 GB of RAM
    ram = psutil.virtual_memory().available/1e9
    required = pop_size/225e3*len(seeds)
    if required < ram:
        print(f'You have {ram} GB of RAM, and this is estimated to require {required} GB: you should be fine')
    else:
        raise ValueError(f'You have {ram:0.2f} GB of RAM, but this is estimated to require {required} GB')
    sc.parallelize(cache_populations, iterarg=seeds) # Run them in parallel
