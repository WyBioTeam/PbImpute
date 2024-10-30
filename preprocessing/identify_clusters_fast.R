args <- commandArgs()
baseName <- "D:\\PbImpute\\preprocessing\\"
if (!require("Seurat")) {
  install.packages("Seurat", dependencies = TRUE, repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(data.table)
print("666666666666666666666")
#sessionInfo()
# pre
data_file <- as.character(file.path(baseName, "selected_genes_expression_matrix.csv"))
# data <- read.table(data_file, sep = '\t', header = TRUE)
#data <- read.table(data_file,sep = '\t')

counts_matrix <- read.table("D:\\PbImpute\\preprocessing\\selected_genes_expression_matrix.csv", sep = ',')
print("666666666666666666666")
print(dim(counts_matrix))


#counts_matrix <- data[, -1]# 
#counts_matrix <- data[-1, ]
w10x_new <- CreateSeuratObject(counts = counts_matrix, raw.data = data, min.cells = 3, min.genes = 200, project = "MUT")
print(ncol(w10x_new))
print("Seurat Process")
print(length(VariableFeatures(w10x_new)))
# norm
#w10x_new <- NormalizeData(object = w10x_new, normalization.method = "LogNormalize", scale.factor = 10000)
w10x_new <- FindVariableFeatures(object = w10x_new, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.01, x.high.cutoff = 5, y.cutoff = 0.01)
print(length(VariableFeatures(w10x_new)))
w10x_new <- ScaleData(object = w10x_new)
# PCA
w10x_new <- RunPCA(w10x_new, pc.genes = w10x_new@var.genes, pcs.compute = 30, do.print = FALSE) 
w10x_new <- FindNeighbors(w10x_new, dims = 1:30)

# Specify the range of resolution values to explore
res <- c(0.05,0.1,0.15,0.2,0.25)

# Iterate through different resolution values
for (i in 1:length(res)) {
 w10x_new <- FindClusters(w10x_new, reduction.type = "pca",
       print.output = 0,force.recalc = T,
                           algorithm = 1, 
                           n.start = 800,    # Set the number of times for different initialization points
                           save.SNN = TRUE, resolution = res[i])  # Save similarity matrix
 
  
  # Visualize t-SNE results
   #DimPlot(w10x_new, group.by = "seurat_clusters")

  # Obtain clustering results
  cluster_results <- data.frame(Cell = names(w10x_new@active.ident), Cluster = w10x_new@active.ident)
  baseName <- "D:\\PbImpute\\preprocessing\\"
  # Write clustering results
   out <- paste(baseName, paste("identity_clustering_res", res[i], ".txt", sep = ""), sep = "/")
   write.table(w10x_new@active.ident, file = out, sep = '\t')

  # Save t-SNE plot
  #plot_filename <- paste(baseName, paste("tSNE_plot_res", res[i], ".jpg", sep = ""), sep = "/")
  #ggsave(plot_filename, device = "jpg", width = 6, height = 6)  
}

  





