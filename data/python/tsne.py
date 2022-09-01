# %% Imports

import keyfi as kf
import pandas as pd

from math import sqrt
from sklearn.preprocessing import PowerTransformer

# %% Read database

df = pd.read_csv('database.csv')

# %% Clean database

# Remove input variables and run time output variable
input_labels = ['T_in', 'eta', 'phi', 'omega']
df.drop(labels=input_labels + ['runtime'], axis=1, inplace=True)

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

# %% Scale data

#df[:] = PowerTransformer(method="box-cox").fit_transform(df)

# %% Run dimensionality reduction with t-SNE or UMAP

# Find rounded sqrt of number of points
N = len(df.index)
n = 2 * round(sqrt(N))

# Use flag to change between algorithms
flag_tsne = 1

if flag_tsne:
    embedding, mapper = kf.embed_data(
        data=df,
        algorithm=kf.dimred.TSNE,
        scale=True,
        perplexity=n,
        init='pca'
    )
else:
    embedding, mapper = kf.embed_data(
        data=df,
        algorithm=kf.dimred.UMAP,
        scale=True,
        n_neighbors=n,
        min_dist=0.01,
        random_state=42,
        n_components=2
    )

# %% Plot mapping

kf.plot_embedding(
    embedding=embedding,
    data=df,
    scale_points=True,
    cmap_var='SL',
    cmap_minmax=[]
)

# %% Run clustering with HDBSCAN and plot it

clusterer = kf.cluster_embedding(
    embedding=embedding,
    algorithm=kf.cluster.HDBSCAN,
    min_cluster_size=n,
    prediction_data=True
)

kf.plot_cluster_membership(
    embedding=embedding,
    clusterer=clusterer,
    soft=True
)

# %% Plot MI scores for both clusters

cluster_mi_scores = kf.get_cluster_mi_scores(
    data=df,
    clusterer=clusterer,
    embedding=embedding,
    cluster_num=0,
    scale=False,
    flag_print=False,
    flag_plot=True
)

cluster_mi_scores = kf.get_cluster_mi_scores(
    data=df,
    clusterer=clusterer,
    embedding=embedding,
    cluster_num=1,
    scale=False,
    flag_print=False,
    flag_plot=True
)

# %%
