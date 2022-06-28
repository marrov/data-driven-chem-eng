# %% Import libraries

import os
import numpy as np
import pandas as pd
import cantera as ct

# %% Adiabatic temperature: set variables

P_in = 101325
T_in = 743
mechanism = './SanDiego_NH3-H2.cti'
#mechanism = 'data/python/SanDiego_NH3-H2.cti'

# %% Main function for adiabatic flame temp (Tad)

gas   = ct.Solution(mechanism)
eta   = [0.1]#np.linspace(0.0, 1.0, num=11) #101ammonia decomposition rate
eq    = [0.5]#np.linspace(0.5, 1.3, num=41) #81step = 0.1
omega = [0.1]#np.linspace(0, 1, num=1001) #101step = 0.01
    #print (eta, eq, omega)

data  = []
for i in range(len(eta)):
    #sum of these must be 1        
    nN2  = 0.25 * eta[i]
    nH2  = 0.75 * eta[i]
    nNH3 = 1 - eta [i]
    for j in range(len(eq)):
        fuel  = 'NH3:{},H2:{},N2:{}'.format(nNH3, nH2, nN2)
        gas.set_equivalence_ratio(eq[j], fuel, 'O2:0.21, N2:0.79')
        nN2_in_dry = gas['N2'].X
        nO2_in_dry = gas['O2'].X
        nH2_in_dry = gas['H2'].X
        nNH3_in_dry= gas['NH3'].X
        nN2_in_air = nO2_in_dry * 3.762
        mair_in_dry= gas['N2'].molecular_weights * nN2_in_air + gas['O2'].molecular_weights * nO2_in_dry

        for k in range(len(omega)):
            msteam     = omega[k] * mair_in_dry
            nsteam     = msteam / gas['H2O'].molecular_weights
            n_wet_sum  = nN2_in_dry + nO2_in_dry + nH2_in_dry + nNH3_in_dry + nsteam
            XN2        = nN2_in_dry  / (n_wet_sum) # wet base reactant
            XO2        = nO2_in_dry  / (n_wet_sum) # wet base reactant
            XH2        = nH2_in_dry  / (n_wet_sum) # wet base reactant
            XNH3       = nNH3_in_dry / (n_wet_sum) # wet base reactant
            Xsteam     = nsteam      / (n_wet_sum) # wet base reactant

            X = 'NH3:{},H2:{},N2:{},O2:{},H2O:{}'.format(XNH3[0], XH2[0], XN2[0], XO2[0], Xsteam[0])

            gas.TPX = T_in, P_in, X
            eq_wet  = gas.get_equivalence_ratio()
            gas.equilibrate('HP')

            # check sum of X = 1 during wet
            XNO_p      = gas['NO'].X # wet base product
            XNO2_p     = gas['NO2'].X # wet base product
            XNH3_p     = gas['NH3'].X # wet base product
            Xsteam_p   = gas['H2O'].X # wet base product
            XO2_p      = gas['O2'].X # wet base product
            XNH3_ppmvd = (XNH3_p / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06 ##ppmvd product
            XNO_ppmvd  = (XNO_p  / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06 ##ppmvd product
            XNO2_ppmvd = (XNO2_p  / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06 ##ppmvd product
            data.append((eta[i], eq[j], omega[k], gas.T, Xsteam[0], XO2[0], XN2[0], XH2[0], XNH3[0], XNO_ppmvd[0], XNO2_ppmvd[0], XNH3_ppmvd[0]))

df    = pd.DataFrame(data, columns = ['eta', 'eq', 'omega', 'T', 'Xsteam', 'XO2', 'XN2', 'XH2', 'XNH3', 'XNO_ppmvd', 'XNO2_ppmvd', 'XNH3_ppmvd'])
df.to_csv('eq.csv', index=False)

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
    Dx=np.diff(x)
    DT=np.diff(T)
    array = DT/Dx
    mid_index = np.argmax(np.abs(array))
    length = 0.0
    for i in range(len(np.arange(mid_index+1))):
        length = length + Dx[i]
    NOposition2 = length + u*0.1 ##0.1s after the point 'length' defined by where T gradient is max
    NO_ind = np.argmin(abs(x - NOposition2))
    return NO_ind

# %% Main function for flame speed (SL)

data = pd.read_csv('eq.csv')

Lx=0.02
tol_ss      = [1.0e-6, 1.0e-14]        # [rtol atol] for steady-state problem
tol_ts      = [1.0e-5, 1.0e-13]        # [rtol atol] for time stepping
loglevel    = 0                        # amount of diagnostic output (0
refine_grid = True                     # True to enable refinement

os.makedirs('./flame', exist_ok=True)
data2  = []
for i in range(len(data)):
    df = data.iloc[i]
    XN2        = df['XN2'] # wet base reactant
    XO2        = df['XO2'] # wet base reactant
    XH2        = df['XH2'] # wet base reactant
    XNH3       = df['XNH3'] # wet base reactant
    Xsteam     = df['Xsteam'] # wet base reactant
    eta        = df['eta'] # decomposition
    eq         = df['eq'] # wet base equivalence ratio
    omega      = df['omega'] # wet base steam/air
    X = 'NH3:{},H2:{},N2:{},O2:{},H2O:{}'.format(XNH3, XH2, XN2, XO2, Xsteam)
    gas.TPX = T_in, P_in, X
    f = ct.FreeFlame(gas, width=Lx)
    f.transport_model = 'Multi'
    f.soret_enabled=True
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.set_refine_criteria(ratio=3, slope=0.01, curve=0.01)
    f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=True)
    SL = f.u[0]
    delta = get_thermal_thickness(f)
    index     = NO_index(f) ##return the index at which NO is picked up
    idx_NO    = gas.species_index('NO')
    idx_O2    = gas.species_index('O2')
    idx_NH3   = gas.species_index('NH3')
    idx_steam = gas.species_index('H2O')
    idx_NO2   = gas.species_index('NO2')
    XNO_p      = f.X[idx_NO][index] # wet base product
    XNO2_p     = f.X[idx_NO2][index] # wet base product
    XNH3_p     = f.X[idx_NH3][index] # wet base product
    Xsteam_p   = f.X[idx_steam][index] # wet base product
    XO2_p      = f.X[idx_O2][index] # wet base product
    XNH3_ppmvd = (XNH3_p / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06 ##ppmvd product
    XNO_ppmvd  = (XNO_p  / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06 ##ppmvd product
    XNO2_ppmvd = (XNO2_p  / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06 ##ppmvd product
    f.write_csv('./flame/flame{}.csv'.format(i+1), species='X')
    data2.append((eta, eq, omega, gas.T, Xsteam, XO2, XN2, XH2, XNH3, XNO_ppmvd, XNO2_ppmvd, XNH3_ppmvd))
df    = pd.DataFrame(data, columns = ['eta', 'eq', 'omega', 'T', 'Xsteam', 'XO2', 'XN2', 'XH2', 'XNH3', 'XNO_ppmvd', 'XNO2_ppmvd', 'XNH3_ppmvd'])
df.to_csv('SL.csv')

# %%
