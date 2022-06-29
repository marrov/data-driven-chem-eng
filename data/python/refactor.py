# %% Import libraries

import numpy as np
import pandas as pd
import cantera as ct

# %% Define all inputs

# Fixed inputs
P_in = 101325
mechanism = './SanDiego_NH3-H2.cti'
#mechanism = 'data/python/SanDiego_NH3-H2.cti'

# Variable inputs
T_in = 743
eta = 0.1  # ammonia decomposition rate
phi = 0.5  # equivalence ratio
omega = 0.1  # steam-to-air ratio

# %% Helper functions for adiabatic flame temp


def get_fuel(eta: float) -> dict:
    N2 = 0.25 * eta
    H2 = 0.75 * eta
    NH3 = 1 - eta
    return {'NH3': NH3, 'H2': H2, 'N2': N2}


def get_species(gas) -> dict:
    "Units of output: mole fraction"
    nsteam = (omega * (gas['N2'].molecular_weights[0] * (gas['O2'].X[0] * 3.762) +
                       gas['O2'].molecular_weights[0] * gas['O2'].X[0])) / gas['H2O'].molecular_weights[0]
    n_wet_sum = gas['N2'].X[0] + gas['O2'].X[0] + \
        gas['H2'].X[0] + gas['NH3'].X[0] + nsteam
    XN2 = gas['N2'].X[0] / n_wet_sum
    XO2 = gas['O2'].X[0] / n_wet_sum
    XH2 = gas['H2'].X[0] / n_wet_sum
    XNH3 = gas['NH3'].X[0] / n_wet_sum
    XH2O = nsteam / n_wet_sum
    return {'NH3': XNH3, 'H2': XH2, 'N2': XN2, 'O2': XO2, 'H2O': XH2O}


def get_products(gas) -> tuple:
    "Units of output: ppmvd (parts per million by volume, dry)"
    XNO = (gas['NO'].X[0] / (1 - gas['H2O'].X[0]) * (20.9 - 15) /
           (20.9 - gas['O2'].X[0] / (1 - gas['H2O'].X[0]))) * 1e6
    XNO2 = (gas['NO2'].X[0] / (1 - gas['H2O'].X[0]) * (20.9 - 15) /
            (20.9 - gas['O2'].X[0] / (1 - gas['H2O'].X[0]))) * 1e6
    XNH3 = (gas['NH3'].X[0] / (1 - gas['H2O'].X[0]) * (20.9 - 15) /
            (20.9 - gas['O2'].X[0] / (1 - gas['H2O'].X[0]))) * 1e6
    return XNO, XNO2, XNH3

# %% Main function for adiabatic flame temp (Tad)

gas = ct.Solution(mechanism)
fuel = get_fuel(eta)
gas.set_equivalence_ratio(phi, fuel, 'O2:0.21, N2:0.79')
X = get_species(gas)
gas.TPX = T_in, P_in, X
gas.equilibrate('HP')
(XNO, XNO2, XNH3) = get_products(gas)

# %% Helper functions for computing flame speed


def get_thermal_thickness(f):
    """ Calculate flame thickness """
    x = f.grid
    T = f.T
    Dx = np.diff(x)
    DT = np.diff(T)
    return (np.amax(T)-np.amin(T))/np.amax(np.abs(DT/Dx))


def NO_index(f):
    """ Get index where emissions are calculated """
    x = f.grid
    u = f.u[0]
    T = f.T
    Dx = np.diff(x)
    DT = np.diff(T)
    array = DT/Dx
    mid_index = np.argmax(np.abs(array))
    length = 0.0
    for i in range(len(np.arange(mid_index+1))):
        length = length + Dx[i]
    # 0.1s after the point 'length' defined by where T gradient is max
    NOposition2 = length + u*0.1
    NO_ind = np.argmin(abs(x - NOposition2))
    return NO_ind


def get_flame(gas):
    Lx = 0.02
    tol_ss = [1.0e-6, 1.0e-14] # [rtol atol] for steady-state problem
    tol_ts = [1.0e-5, 1.0e-13] # [rtol atol] for time stepping
    loglevel = 0 # amount of diagnostic output
    refine_grid = True # True to enable refinement
    f = ct.FreeFlame(gas, width=Lx)
    f.transport_model = 'Multi'
    f.soret_enabled = True
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.set_refine_criteria(ratio=3, slope=0.01, curve=0.01)
    f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=True)
    return gas, f


# %% Main function for flame speed (SL)

gas.TPX = T_in, P_in, X # Reset gas state
(gas, f) = get_flame(gas) 
SL = f.u[0]
delta = get_thermal_thickness(f)

labels = ['eta', 'phi', 'omega', 'T', 'NO', 'NO2', 'NH3', 'SL', 'delta']
data = [eta, phi, omega, gas.T, XNO, XNO2, XNH3, SL, delta]

# %%
