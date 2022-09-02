# %% Imports

import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression

# %% Read dataset

df = pd.read_csv('database.csv')

# %% Clean database

# Keep only input variables and adiabatic flame temperature
df = df[['T_in', 'eta', 'phi', 'omega', 'Tad']]

# Remove rows with outliers based on IQR
def remove_outlier_IQR(df):
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    df = df[~((df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR)))]
    return df

df = remove_outlier_IQR(df)

# Drop NaN values
df.dropna(inplace=True)

# Reindex dataframe
df.reset_index(drop=True, inplace=True)

# %% Split into test and train

X = df[['T_in', 'eta', 'phi', 'omega']]
y = df['Tad']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# %% Scala data with standard scaler

scaler = StandardScaler()

X_train_scaled = pd.DataFrame(scaler.fit_transform(X_train), columns = X.columns)
X_test_scaled = pd.DataFrame(scaler.transform(X_test), columns = X.columns)

# %% Run linear regression

linear_model = LinearRegression().fit(X_train_scaled, y_train)
print(linear_model.score(X_train_scaled, y_train))

# %% Show scatterplots
sns.scatterplot(data=df, x='eta', y='Tad')

# %%
sns.scatterplot(data=df, x='T_in', y='Tad')

# %%
sns.scatterplot(data=df, x='phi', y='Tad')

# %%
sns.scatterplot(data=df, x='omega', y='Tad')

# %%
