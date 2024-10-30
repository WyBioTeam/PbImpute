#---------------------------------------------------------------1
#calculate the response of EM algorithm
response<- function(expr, phi, alpha, lambda){
  p<- alpha/(alpha+lambda)
  b<- phi/(phi+ (1-phi)*(p^alpha)) 
  results<- rep(0, length(expr))
  results[expr==0]<- b
  return(results)
}
# The response can distinguish dropout zeros.

f<- function(paramt, weights, x){ #Update the two parameters of the Negative Binomial Distribution.
  alpha<- paramt[1]
  lambda<- paramt[2]
  results<- (1-weights)*log(dnbinom(x, alpha, prob= alpha/(alpha+lambda)))
  return(-sum(results))
}


update_nbinom <- function(x, weights, paramt) {
  # print("The initial parameters for updating are:")
  # print(paramt)
  # Use the optim function for parameter optimization.
  pars <- optim(paramt, f, weights = weights, x = x)
  # print(pars$par)
  
  return(pars$par)
}

#对数似然函数：
loglik<- function(zinb_data, weights, phi, alpha, lambda){
  part1<- (zinb_data==0)*log(phi+(1-phi)*((alpha/(alpha+lambda))^alpha))
  part2<- (zinb_data>0)*(log(1-phi)+ dnbinom(zinb_data, alpha, prob= alpha/(alpha+lambda),
                                             log = T))
  return(sum(part1)+sum(part2))
}


zinb<- function(zinb_data){
  # 1.Two initial values
  # print(zinb_data[zinb_data != 0])
  
  lambda<- mean(zinb_data)
  
  
  if(lambda==0){ 
  # After cell clustering, it is possible that a certain gene is not expressed in all cells of a given cluster."
  #   "Later, it may be necessary to save these genes that are completely unexpressed in a specific cluster as true zeros."
  #   "Keep all of them directly
    estimates<- rep(NA, 4)
    names(estimates)<- c('alpha(size)', 'lambda', 'phi', 'prob')
    return(estimates)
  }
  
  alpha<- (lambda^2)/((sd(zinb_data)^2)- lambda)
  # print("alpha1111111111")
  # print(alpha)
  if(alpha<=0||alpha==Inf){
    estimates<- rep(NA, 4)
    names(estimates)<- c('alpha(size)', 'lambda', 'phi', 'prob')
    return(estimates)
  }
  # print("alpha22222")
  # 2. Estimation of π0
  #  Generate samples from a Negative Binomial Distribution.
  est_sample<- rnbinom(length(zinb_data), alpha, alpha/(alpha+lambda))
  if(sum(zinb_data==0)>sum(est_sample==0)){
    # If the number of zeros in the original data zinb_data is greater than the number of zeros in the generated sample est_sample, calculate the zero-inflation parameter
    # 𝜙
    phi<- (sum(zinb_data==0)-sum(est_sample==0))/sum(zinb_data!=0)
    if(phi>=1){
      phi<- sum(zinb_data==0)/length(zinb_data)
    }
    #The previous line may result in phi being greater than 1, which is unreasonable.
  }else{
    phi<- sum(zinb_data==0)/length(zinb_data) 
  }
  
  # 3.Parameter Estimation
  old_lik<- 0
  t<- 1
  paramt<- c(alpha, lambda)
  estimates<- c(alpha, lambda, phi, alpha/(alpha+lambda))
  names(estimates)<- c('alpha(size)', 'lambda', 'phi', 'prob')
  while(t<=20){
    weights<- response(zinb_data, phi, paramt[1], paramt[2])
    # if(t==1){
    #   print("---------------")
    #   print("alpha0")
    #   print(paramt[1])
    #   print("lambda0")
    #   print(paramt[2])
    #   print("phi0")
    #   print(phi)
    #   print("---------------")
    # }
    paramt<- update_nbinom(zinb_data, weights, paramt)
    alpha<- paramt[1]
    lambda<- paramt[2]
    phi<- mean(weights)#
    estimates<- rbind(estimates, c(alpha, lambda, phi, alpha/(alpha+lambda)))
    #Whether it converges
    lik<- loglik(zinb_data, weights, phi, alpha, lambda)
    #print(lik)
    eps<- (lik- old_lik)^2
    if(is.na(lik)){
      print(estimates)
    }
    old_lik<- lik
    if(eps<0.000001){
      break
    }
    t<- t+1
  }
  return(estimates[nrow(estimates),])
}
get_parameters<- function(count){
  options(warn= -1)
  parslist = lapply(1:nrow(count), function(ii) { 
    xdata = count[ii, ] #xdata是某个基因向量
    paramt = zinb(xdata)
  
    if (ii %%1000 == 0 && all(!is.na(paramt))) {
      print(paste0("Parameter estimation of row: ", paste(paramt, collapse = ", ")))
      
    }
   
    return(paramt)
  })
  parslist = Reduce(rbind, parslist)
  colnames(parslist) = c('alpha(size)', 'lambda', 'phi', 'prob')
  return(parslist)
}
#--------------------------------------------------------------------------2
zinb_impute<- function(j, subcount, parlist, dropout_thres, scale_factor){
  gene<- subcount[j,]
  phi <- parlist[j, 'phi']
  prob <- parlist[j, 'prob']
  alpha <- parlist[j, 'alpha(size)']
  freq<- sum(gene==0)/ncol(subcount)
  dropout_prob <- phi / (phi + (1 - phi) * (prob ^ alpha))
  dropout_prob<- ifelse(!is.na(dropout_prob), dropout_prob, 0)
  prior<- ifelse((dropout_prob>=dropout_thres)
                 &&(!is.na(parlist[j, 'lambda'])), parlist[j, 'lambda'], 0)
  impute<- (gene==0)*prior*scale_factor*
    ifelse(is.na(dropout_prob), 0, dropout_prob)#
  gene<- gene+ impute
  return(gene)
}
init_impute<- function(count, cell_cluster, dropout_quantile= 0.5, augment= 1){
  for (i in unique(cell_cluster)){
    print(paste0("Distribution imputation for category ", i, " begins..."))
    parlist<- readRDS(paste0('pars_cluster_',i,'.rds'))
    cell_id<- colnames(count[,cell_cluster==i]) #The current cell 
    subcount<- count[ ,cell_id] 
    scale_factor<- colSums(subcount)/mean(colSums(subcount))*augment #Scaling factor.
    dropout_thres<- quantile(parlist[,'phi'], dropout_quantile, na.rm=TRUE)
    ncores<- detectCores()
    cl<- makeCluster(ncores)#quick
    imputed_genes<- parLapply(cl, 1:nrow(subcount), zinb_impute, subcount, parlist,
                              dropout_thres, scale_factor)
    imputed_genes<- Reduce(rbind, imputed_genes)
    #    print(colnames(imputed_genes))
    count[, cell_id]<- imputed_genes
    stopCluster(cl)
    print(paste0("Distribution imputation for category ", i, " is complete..."))
  }
  return(count)
}
#----------------------------------------------------------------------3
#'
#' @param count A count matrix with genes as rows and cells as columns.
#'
#' @param dropout_quantile Genes with a dropout probability exceeding this quantile will be classified as dropouts.
#'
#' @param augment An augmentation parameter to adjust the imputed values by scaling them up or down after ZINB imputation..
#'
#' @param seed Random seed to use.
#'
#' @return An imputed count matrix with rows as genes and columns as cells.

pbimpute1<- function(count, dropout_quantile= 0.5,
                       augment= 1, seed= 1){

  library(Seurat)
  library(factoextra)
  library(parallel)
  count<- as.matrix(count)
  print(paste0('The sample contains ', ncol(count), ' cells and ',nrow(count),
               ' genes'))
  if (is.null(colnames(count))|is.null(rownames(count))){
    stop('The sample should have both column names and row names.')
  }
  row.names(count)<- gsub('_','-', row.names(count))
  #It is recommended to select highly variable features to prepare for imputation.
  library(Seurat)
  pbmc <- CreateSeuratObject(counts = count)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst")
  # test 2000
  top <- head(VariableFeatures(pbmc),2000)
  topmatrix<- pbmc@assays$RNA@counts[top, ]
  filtered<- topmatrix
  filtered <- as.matrix(filtered)  # 
  print('Reading the Cell partitioning file...')
  print(dim(count))
  
  file_path <- "D:\\PbImpute\\preprocessing\\group.txt"
  cluster_numbers <- as.integer(read.table(file_path, header = FALSE)[, 1])
  # Obtain the column names of the count matrix as cell names.
  cell_names <- colnames(filtered)
  if (length(cluster_numbers) != length(cell_names)) {
    stop("The length of cluster numbers does not match the number of columns in the count matrix.")
  }
  cell_cluster <- cluster_numbers
  names(cell_cluster) <- cell_names
  print('estimating the distribution parameters of the genes...')
  for(i in unique(cell_cluster)){
    print("-------- Category n --------")
    print(i)
    subcount<- filtered[,cell_cluster==i]
    #Estimation parameters are placed here:
    abc<- get_parameters(count = subcount)
    saveRDS(abc, paste0('pars_','cluster_',i,'.rds'))
    # First, mark the positions with zero values in the expression matrix. If the dropout probability of the corresponding gene exceeds the threshold, the position is considered a dropout.
    # #  dropout_thres<- quantile(abc[,'phi'], 0.8, na.rm=TRUE)
  }
  print('pbimpute1 imputation...')
  dropout_quantile<- dropout_quantile
  xxx<- init_impute(filtered, cell_cluster, dropout_quantile = dropout_quantile, augment = augment) #这个就是初次填充后的结果
  # If you want to retain all imputations, you can output<- count
  output<- count
  output[rownames(xxx),]<- xxx
  print('pbimpute1 finished...')
  print(' congratulations!')
  return (output) 
}

imputed_data <- read.csv("D:\\PbImpute\\test_data\\gene_cell.csv", header = TRUE, row.names = 1)
print(dim(imputed_data))
output<- pbimpute1(imputed_data,dropout_quantile= 0.5,seed= 1)
mean_vals <- colMeans(output, na.rm = TRUE)  
std_dev_vals <- apply(output, 2, sd, na.rm = TRUE)  
z_scores <- scale(output, center = mean_vals, scale = std_dev_vals)
mask <- z_scores < -0.3
output[mask] <- 0
# 
write.csv(output, 
          "D:/PbImpute/ZINB_imputation/PbImpute1.csv")

library(Seurat)
output_norm<- NormalizeData(output, normalization.method = "LogNormalize", scale.factor = 10000)
write.csv(output_norm, 
          "D:/PbImpute/ZINB_imputation/log_PbImpute1.csv", 
          row.names = TRUE, )  # 






