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

# %% Helper functions

def get_fuel(eta: float) -> dict:
    N2 = 0.25 * eta
    H2 = 0.75 * eta
    NH3 = 1 - eta
    return {'NH3': NH3, 'H2': H2, 'N2': N2}


def get_species(gas) -> dict:
    nsteam = (omega * (gas['N2'].molecular_weights[0] * (gas['O2'].X[0] * 3.762) +
                       gas['O2'].molecular_weights[0] * gas['O2'].X[0])) / gas['H2O'].molecular_weights[0]
    n_wet_sum = gas['N2'].X[0] + gas['O2'].X[0] + gas['H2'].X[0] + gas['NH3'].X[0] + nsteam
    XN2 = gas['N2'].X[0] / n_wet_sum
    XO2 = gas['O2'].X[0] / n_wet_sum
    XH2 = gas['H2'].X[0] / n_wet_sum
    XNH3 = gas['NH3'].X[0] / n_wet_sum
    XH2O = nsteam / n_wet_sum
    return {'NH3': XNH3, 'H2': XH2, 'N2': XN2, 'O2': XO2, 'H2O': XH2O}


def get_products(gas) -> tuple:
    XNH3 = (gas['NH3'].X[0] / (1 - gas['H2O'].X[0]) * (20.9 - 15) / (20.9 - gas['O2'].X[0] / (1 - gas['H2O'].X[0]))) * 1e6 
    XNO = (gas['NO'].X[0] / (1 - gas['H2O'].X[0]) * (20.9 - 15) / (20.9 - gas['O2'].X[0] / (1 - gas['H2O'].X[0]))) * 1e6 
    XNO2 = (gas['NO2'].X[0] / (1 - gas['H2O'].X[0]) * (20.9 - 15) / (20.9 - gas['O2'].X[0] / (1 - gas['H2O'].X[0]))) * 1e6 
    return XNH3, XNO, XNO2

# %% Main function for adiabatic flame temp (Tad)

data = []
gas = ct.Solution(mechanism)
fuel = get_fuel(eta)
gas.set_equivalence_ratio(phi, fuel, 'O2:0.21, N2:0.79')
X = get_species(gas)
gas.TPX = T_in, P_in, X
gas.equilibrate('HP')
(XNO, XNO2, XNH3) = get_products(gas)
data.append((eta, phi, omega, gas.T, XNO, XNO2, XNH3))
labels = ['eta', 'phi', 'omega', 'T', 'XNO', 'XNO2', 'XNH3']
df = pd.DataFrame(data, columns=labels)

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

# %% Main function for flame speed (SL)


data = pd.read_csv('eq.csv')

Lx = 0.02
tol_ss = [1.0e-6, 1.0e-14]        # [rtol atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-13]        # [rtol atol] for time stepping
loglevel = 0                        # amount of diagnostic output (0
refine_grid = True                     # True to enable refinement

data2 = []
for i in range(len(data)):
    df = data.iloc[i]
    XN2 = df['XN2']  # wet base reactant
    XO2 = df['XO2']  # wet base reactant
    XH2 = df['XH2']  # wet base reactant
    XNH3 = df['XNH3']  # wet base reactant
    Xsteam = df['Xsteam']  # wet base reactant
    eta = df['eta']  # decomposition
    phi = df['phi']  # wet base equivalence ratio
    omega = df['omega']  # wet base steam/air
    X = 'NH3:{},H2:{},N2:{},O2:{},H2O:{}'.format(XNH3, XH2, XN2, XO2, Xsteam)
    gas.TPX = T_in, P_in, X
    f = ct.FreeFlame(gas, width=Lx)
    f.transport_model = 'Multi'
    f.soret_enabled = True
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.set_refine_criteria(ratio=3, slope=0.01, curve=0.01)
    f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=True)
    SL = f.u[0]
    delta = get_thermal_thickness(f)
    index = NO_index(f)  # return the index at which NO is picked up
    idx_NO = gas.species_index('NO')
    idx_O2 = gas.species_index('O2')
    idx_NH3 = gas.species_index('NH3')
    idx_steam = gas.species_index('H2O')
    idx_NO2 = gas.species_index('NO2')
    XNO_p = f.X[idx_NO][index]  # wet base product
    XNO2_p = f.X[idx_NO2][index]  # wet base product
    XNH3_p = f.X[idx_NH3][index]  # wet base product
    Xsteam_p = f.X[idx_steam][index]  # wet base product
    XO2_p = f.X[idx_O2][index]  # wet base product
    XNH3_ppmvd = (XNH3_p / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06  # ppmvd product
    XNO_ppmvd = (XNO_p / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06  # ppmvd product
    XNO2_ppmvd = (XNO2_p / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06  # ppmvd product
    f.write_csv('flame{}.csv'.format(i+1), species='X')
    data2.append((eta, phi, omega, gas.T, Xsteam, XO2, XN2, XH2, XNH3, XNO_ppmvd, XNO2_ppmvd, XNH3_ppmvd))
df = pd.DataFrame(data, columns=['eta', 'phi', 'omega', 'T', 'Xsteam', 'XO2', 'XN2', 'XH2', 'XNH3', 'XNO_ppmvd', 'XNO2_ppmvd', 'XNH3_ppmvd'])
df.to_csv('SL.csv')
