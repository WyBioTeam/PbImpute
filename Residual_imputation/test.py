#change directory to src
from SCIMP import Impute
from SCIMP import identification
from SCIMP import whole_reconstrtor
import pandas as pd
import numpy as np
import scanpy as sc
import os
import numpy as np
import pandas as pd
##############################graph#########################################
#Step1 load the raw count matrix of scRNA-seq data, where rows are genes and columns are cells.
print("----------------------------------------")
print("Residual imputation is starting, please wait...")
file_name = r'D:\PbImpute\ZINB_imputation\\PbImpute1.csv'
rawfile = pd.read_csv(file_name, sep=',', header=0, index_col=0)

log_file_name = r'D:\PbImpute\ZINB_imputation\\log_PbImpute1.csv'
log_rawfile = pd.read_csv(log_file_name, sep=',', header=0, index_col=0)

#Step2 build adjacent matrix for scRNA-seq data.
graph_adj=Impute.GraphBuild(rawfile,log_rawfile)

#Step3 learn cell embeddings.
cell_emb=Impute.trainCellEmbeddings(graph_adj)

#Step4 scRNA-seq data imputation, the format of output file is genes x cells expression matrix.
data_imp=Impute.imputation(scfile=rawfile,embeddingfile=cell_emb,AdjGraph=graph_adj)
Identification_data=identification.process_and_identify_dropouts(log_file_name)

Identification_data =pd.DataFrame(Identification_data)
data_imp =pd.DataFrame(data_imp)
rawfile =pd.DataFrame(rawfile)
mask = Identification_data == -1
rawfile[mask] = data_imp[mask]
rawfile11=rawfile.copy()
###################repair` 2#################################################
adata = sc.AnnData(rawfile)
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)  # 进行 log-normali    zation
# 现在 `adata` 对象中包含了归一化后的数据
normalized_data = adata.to_df()  # 转换回 pandas DataFrame
aa=whole_reconstrtor.balance_data(normalized_data,rawfile11)
print(aa.shape)
aa.index = rawfile.index
aa.columns = rawfile.columns
# Save the aligned data to a CSV file
print(" Congratulations! Residual imputation has finished!")
print("----------------------------------------")
aa.to_csv(r'D:\PbImpute\Residual_imputation\final_PbImpute.csv', index=True)