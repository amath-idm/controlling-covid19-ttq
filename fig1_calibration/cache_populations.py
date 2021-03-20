'''
Pre-generate the population. Takes about 70 s. 

Requires SynthPops, which must be installed from the repository:

https://github.com/institutefordiseasemodeling/synthpops
'''

import psutil
import sciris  as sc
import covasim as cv
import synthpops as sp
sp.config.set_nbrackets(20) # Essential for getting the age distribution right
sp.logger.setLevel('DEBUG')


pop_size = 225e3

def cache_populations(seed=0, popfile=None):
    ''' Pre-generate the synthpops population '''

    pars = sc.objdict(
        pop_size = pop_size,
        pop_type = 'synthpops',
        rand_seed = seed,
    )

    if popfile is None:
        popfile = f'inputs/kc225_rnr_seed{pars.rand_seed}.ppl'

    T = sc.tic()
    print(f'Making "{popfile}"...')
    sim = cv.Sim(pars)
    cv.make_people(sim, popfile=popfile, save_pop=True, with_facilities=True, generate=True, layer_mapping={'LTCF':'l'})
    sc.toc(T)

    print('Done')
    return


if __name__ == '__main__':

    seeds = [0,1,2,3,4] # NB, each one takes 8 GB of RAM! -- split up 0-4 in pieces
    ram = psutil.virtual_memory().available/1e9
    required = pop_size/225e3*len(seeds)
    if required < ram:
        print(f'You have {ram} GB of RAM, and this is estimated to require {required} GB: you should be fine')
    else:
        raise ValueError(f'You have {ram:0.2f} GB of RAM, but this is estimated to require {required} GB')
    sc.parallelize(cache_populations, iterarg=seeds) # Run them in parallel
