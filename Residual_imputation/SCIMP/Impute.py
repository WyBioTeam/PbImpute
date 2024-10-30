import numpy as np
import pandas as pd
import networkx as nx
from sklearn.metrics.pairwise import pairwise_distances
from .CellEmbeddings.node2vec import Node2Vec
import networkx as nx

def _calculate_distance(X, num_jobs=5):

    D = pairwise_distances(X.T, n_jobs=num_jobs, metric='euclidean')
    return D

def GraphBuild(input_file,log_rawfile,k=3):
    nor_file=log_rawfile
    distance_matrix=_calculate_distance(nor_file)
    adj_matrix=np.zeros((distance_matrix.shape[0],distance_matrix.shape[0]))
    for i in range(distance_matrix.shape[0]):
        adj_matrix[np.argsort(distance_matrix[:,i])[1:k+1],i]=1
        adj_matrix[i,np.argsort(distance_matrix[:,i])[1:k+1]]=1
    adj_df=pd.DataFrame(adj_matrix)
    adj_df.index=input_file.columns.tolist()
    adj_df.columns=input_file.columns.tolist()
    return adj_df

def trainCellEmbeddings(Graph,Seed=0,Workers=1):
    G = nx.from_pandas_adjacency(Graph)
    Cell_node2vec = Node2Vec(G,dimensions=256,p=4,q=2, walk_length=5, num_walks=20, workers=Workers,seed=Seed)
    tmodel=Cell_node2vec.fit(window=3, epochs=3)
    emb=Cell_node2vec.get_embeddings(tmodel)
    return emb

def imputation(scfile,embeddingfile,AdjGraph):
    prediction_distance=pairwise_distances(embeddingfile.T,n_jobs=5,metric="euclidean")
    adj_matrix_merge=np.array(AdjGraph)
    for j in range(adj_matrix_merge.shape[0]):
        adj_matrix_merge[np.argsort(prediction_distance[:,j])[0:np.sum(adj_matrix_merge[:,j]==1)],j]=1
    np_scfile=np.array(scfile)
    imputation_file=np.array(scfile)
    for k in range(adj_matrix_merge.shape[1]):
        imputation_file[np.where(imputation_file[:,k]==0)[0],k] = np.mean(np_scfile[:,np.where(adj_matrix_merge[:,k]==1)[0]][np.where(imputation_file[:,k]==0)[0],:],axis=1)
    df_imputation=pd.DataFrame(imputation_file)
    df_imputation.columns=scfile.columns.tolist()
    df_imputation.index=scfile.index.tolist()
    return df_imputation
