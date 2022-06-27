# %% Import libraries

import cantera as ct
import pandas as pd

from calc import *

# %% Adiabatic temperature: set variables

P_in = 101325
T_in = 743
T_out = 1600
mechanism = 'data/python/SanDiego_NH3-H2.cti'
#X_H2_ref = 0.5

# %% Main function

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

            if abs(gas.T-T_out)<1: # if equilibrium temperature < 1700K
                # check sum of X = 1 during wet
                print ('eq dry = {}, eq wet = {} -> they should be equal!'.format(eq[j], eq_wet))
                print ('Omega = {}, wet mixture sum = {}, T = {}\n'.format(omega[k], (XN2 + XO2 + XH2 + XNH3 + Xsteam), gas.T))
                XNO_p      = gas['NO'].X # wet base product
                XNO2_p     = gas['NO2'].X # wet base product
                XNH3_p     = gas['NH3'].X # wet base product
                Xsteam_p   = gas['H2O'].X # wet base product
                XO2_p      = gas['O2'].X # wet base product

                XNH3_ppmvd = (XNH3_p / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06 ##ppmvd product
                XNO_ppmvd  = (XNO_p  / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06 ##ppmvd product
                XNO2_ppmvd = (XNO2_p  / (1 - Xsteam_p) * (20.9 - 15) / (20.9 - XO2_p / (1 - Xsteam_p))) * 1e06 ##ppmvd product

                data.append((eta[i], eq[j], omega[k], gas.T, Xsteam[0], XO2[0], XN2[0], XH2[0], XNH3[0], XNO_ppmvd[0], XNO2_ppmvd[0], XNH3_ppmvd[0]))
                break
df    = pd.DataFrame(data, columns = ['eta', 'eq', 'omega', 'T', 'Xsteam', 'XO2', 'XN2', 'XH2', 'XNH3', 'XNO_ppmvd', 'XNO2_ppmvd', 'XNH3_ppmvd'])
df.to_csv('eq.csv')

# %%
