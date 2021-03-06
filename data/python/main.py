# %% Import libraries

import pandas as pd

from helper import *
from cantera import CanteraError
from SALib.sample import saltelli

# %% Define inputs using SALib with saltelli space filling sampling

# Description of the input variables:
# T_in:  inlet temperature [K]
# eta:   ammonia decomposition rate
# phi:   equivalence ratio
# omega: steam-to-air ratio

input_labels = ['T_in', 'eta', 'phi', 'omega']

problem = {'num_vars': len(input_labels),
           'names': input_labels,
           'bounds': [[300, 1000],
                      [0, 1],
                      [0.5, 2],
                      [0, 0.8]]}

inputs = list(saltelli.sample(problem, 2**10, calc_second_order=False))

# %% Run case and get outputs

# Description of the output variables:
# Tad:     Adiabatic flame temperature [K]
# NO:      NO fraction [ppmvd]
# NO2:     NO2 fraction [ppmvd]
# NH3:     Ammonia fraction [ppmvd]
# SL:      Flame speed [m/s]
# delta:   Thermal thickness [m]
# runtime: Code run time [s]

output_labels = ['Tad', 'NO', 'NO2', 'NH3', 'SL', 'delta', 'runtime']
outputs = []

idx_first = 2736 - 2

for idx in range(idx_first, len(inputs)):
    try:
        outputs.append(burn_ammonia(*inputs[idx]))
    except CanteraError as e:
        print(f'Error: {e}. Inputs used were ' + ', '.join(
            [f'{label}={value}' for label, value in zip(input_labels, inputs[idx])]))
        outputs.append([None] * len(output_labels))

#for input in inputs:
#    try:
#        outputs.append(burn_ammonia(*input))
#    except CanteraError as e:
#        print(f'Error: {e}. Inputs used were ' + ', '.join(
#            [f'{label}={value}' for label, value in zip(input_labels, input)]))
#        outputs.append([None] * len(output_labels))

# %% Make dataFrame and export to CSV file

idx_last = idx + 1

df = pd.concat([pd.DataFrame(data=inputs[idx_first:idx_last], columns=input_labels),
                pd.DataFrame(data=outputs, columns=output_labels)], axis=1)

# Indexes first and last point to orig database (both included)
df.to_csv('database_' + str(idx_first + 2) + '_to_' + str(idx_last + 1) + '.csv', index=False, na_rep='NaN')

# %%
