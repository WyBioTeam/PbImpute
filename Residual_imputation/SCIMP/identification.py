import os
import h5py
import sklearn
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os.path import join
from functools import partial
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
def process_and_identify_dropouts(file_path, num_components=25, clusters=3, dropout_thr=0.4, cv_thr=0.4):
    # Read the raw data
    rawdata = pd.read_csv(file_path, sep=',', header=0, index_col=0)
    row_index = rawdata.index
    column_index = rawdata.columns
    data_norm = rawdata.values
    data_norm1 = data_norm.copy()
    # Perform PCA
    pca = PCA(n_components=num_components, svd_solver="randomized")
    pca_data = pca.fit_transform(data_norm)
    # Perform K-means clustering
    kmeans = KMeans(n_clusters=clusters, random_state=0, n_init=10).fit(pca_data)
    label_pr = kmeans.labels_
    print(len(label_pr))
    def find_cluster_cell_idx(l, label):
        return label == l
    def identify_dropout(cluster_cell_idxs, X, dropout_thr, cv_thr):
        for idx in cluster_cell_idxs:
            dropout = (X[:, idx] == 0).sum(axis=1) / (X[:, idx].shape[1])
            dropout_upper_thr, dropout_lower_thr = np.nanquantile(dropout, q=dropout_thr), np.nanquantile(dropout, q=0)
            gene_index1 = (dropout <= dropout_upper_thr) & (dropout >= dropout_lower_thr)
            means = X[:, idx].mean(axis=1)
            means[means == 0] = np.nan  # 避免除以零
            cv = X[:, idx].std(axis=1) / means
            cv_upper_thr, cv_lower_thr = np.nanquantile(cv, q=cv_thr), np.nanquantile(cv, q=0)
            gene_index2 = (cv <= cv_upper_thr) & (cv >= cv_lower_thr)

            include_faslezero_gene = np.logical_and(gene_index1, gene_index2)
            tmp = X[:, idx]
            tmp[include_faslezero_gene] = tmp[include_faslezero_gene] + (tmp[include_faslezero_gene] == 0) * -1
            X[:, idx] = tmp
        return X
    label_set = np.unique(label_pr)
    cluster_cell_idxs = list(map(partial(find_cluster_cell_idx, label=label_pr), label_set))
    data_identi = identify_dropout(cluster_cell_idxs, X=data_norm.T, dropout_thr=dropout_thr, cv_thr=cv_thr).T
    data_identi = pd.DataFrame(data_identi, index=row_index, columns=column_index)
    num_minus_ones = np.count_nonzero(data_identi == -1)
    print("Number of -1 values:", num_minus_ones)
    return data_identi
