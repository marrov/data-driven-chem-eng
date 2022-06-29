# %% Import libraries

from helper import *

# %% Define inputs

T_in = 300
eta = 0.1  # ammonia decomposition rate
phi = 0.5  # equivalence ratio
omega = 0.1  # steam-to-air ratio

# %% Run case and get outputs

(Tad, NO, NO2, NH3, SL, delta) = burn_ammonia(T_in, eta, phi, omega)

# %%
