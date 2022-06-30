#%% Imports
import time

from helper import *
from SALib.sample import saltelli
from UQpy.RunModel import RunModel

#%% Problem definition and sample generation with SALib

# Description of input variables:
# T_in:  inlet temperature
# eta:   ammonia decomposition rate
# phi:   equivalence ratio
# omega: steam-to-air ratio

names = ['T_in', 'eta', 'phi', 'omega']

problem = {
    'num_vars': len(names),
    'names': names,
    'bounds': [[300, 1000],
               [0, 1],
               [0.5, 2],
               [0, 0.8]]
}

params_values = saltelli.sample(problem, 2**0, calc_second_order=False)
params_list = list(params_values)

#%% Run simulations with UQpy

#t = time.time()
#m = RunModel(samples=params_list, ntasks=4, model_script='runner.py',
#             input_template='input_vars.txt', var_names=names, model_object_name="runner",
#             output_script='postprocess.py', output_object_name='postprocess',
#             resume=False, model_dir='runs', verbose=True)
#t_parallel = time.time() - t
#print("\nTime for parallel execution:")
#print(t_parallel)
