
import numpy as np
import stan
from build import model
import matplotlib
matplotlib.use('TkAgg')


def brsindy(antstring, data_frame, rxn_types=None):

    ansatz_rxns, stan_model, num_func = model(len(data_frame.columns)-1, rxn_types)

    ansatzfile = antstring.split('.')[0] + '.ansatz'
    with open(ansatzfile, 'w') as f:
        f.write(ansatz_rxns)

    stanfile = antstring.split('.')[0] + '.stan'
    with open(stanfile, 'w') as f:
        f.write(stan_model)

    num_species = len(data_frame.columns)-1
    spec_ind = [i+1 for i in range(num_species)]

    data = {'N': len(data_frame) - 1,
            'M': num_species,
            'obs_idx': spec_ind,
            'D': num_func,
            'y0': np.array(data_frame.iloc[0, 1:].T),
            'y': np.array(data_frame.iloc[1:, spec_ind]),
            'ts': np.array(data_frame.iloc[:, 0].T),
            'm0': 3,
            'slab_scale': 4,
            'slab_df': 4,
            }

    posterior = stan.build(stan_model, data=data)
    fit = posterior.sample(num_chains=1, num_warmup=2000, num_samples=1000)

    return fit, ansatz_rxns
