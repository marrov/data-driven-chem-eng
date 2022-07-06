# %% Imports

import keyfi as kf
import pandas as pd

from math import sqrt

# %% Read database

df = pd.read_csv('database.csv')

# %% Clean database

df.dropna(inplace=True)

# %% Run dimensionality reduction with t-SNE or UMAP

# Find rounded sqrt of number of points
N = len(df.index)
n = round(sqrt(N))

# Find a reasonable learning rate using default early_exaggeration=12
learning_rate = max(N / 12 / 4, 50)

# Use flag to change between algorithms
flag_tsne = 1

if flag_tsne:
    embedding, mapper = kf.embed_data(
        data=df,
        algorithm=kf.dimred.TSNE,
        scale=True,
        perplexity=n,
        #learning_rate = learning_rate,
        #init='pca'
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
    cmap_var='phi',
    cmap_minmax=[]
)

# %% Run clustering with HDBSCAN and plot it

clusterer = kf.cluster_embedding(
    embedding=embedding,
    algorithm=kf.cluster.HDBSCAN,
    min_cluster_size=n,
    min_samples=15,
    prediction_data=True
)

kf.plot_cluster_membership(
    embedding=embedding,
    clusterer=clusterer,
    soft=True
)

# %%
