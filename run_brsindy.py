
import brsindy
import numpy as np
import pandas as pd
from copy import deepcopy
import tellurium as te
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

if __name__ == "__main__":

    # todo: algorithm to set inference conditions and generate data
    modelname = 'lv.txt'
    model = te.loada(modelname)
    species = model.getFloatingSpeciesIds() + model.getBoundarySpeciesIds()
    cols = deepcopy(species)
    cols.insert(0, 'time')
    sim = model.simulate(0, 15, 30, selections=cols)
    sim_df = pd.DataFrame(sim, columns=cols)

    # It may be necessary to pass the maximum data set and instructions to subset it.
    fit, ansatz_rxns = brsindy.brsindy(modelname, sim_df)

    sparse_parameters = np.round(np.mean(fit['rates'], axis=1), 3)
    for i, line in enumerate(ansatz_rxns.splitlines()):
        print(sparse_parameters[i], line)
