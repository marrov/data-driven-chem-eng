# %% Imports

import pandas as pd
import numpy as np 

# %% Read dataset

df = pd.read_csv('database.csv')

# %% Add noise to database and save as csv

# Define function to add noise
def add_noise(x, rel_std):
  return x + np.random.normal(0, rel_std * np.std(x), size = x.shape)

# Define a relarive standard deviation of 5%
rel_std = 0.05

# Apply function to each column individually and save to csv
df.apply(lambda x: add_noise(x, rel_std)).to_csv('noisy_database.csv', index=False)

# %%
