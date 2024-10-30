function proData = preprocessing(iniData, minCells, minGenes, libararyflag, logNormalize)
% Preprocess scRNAseq data using the following steps
% (1) Filter out low-quality cells
% (2) Filter out low-expressed genes
% (3) log tranformation
% Inputs:
%   iniData: a table variable, store the loaded data
%   minCells: a value, filter genes that are expressed in less than minCells cells; default = 3
%   minGenes: a value, filter low-quality cells in which the number of expressed genes is less than minGenes; default = 200
%   libararyflag: boolean, use a global-scaling normalization method ot
%   not, default = 1
%   logNormalize: boolean, to do log10 normalization or not, default= 1 (true)

% Outputs:
%   proData: a struct variable, store the data after QC
%   proData.data: a matrix giving the single cell data (rows are cells and columns are genes)
%   proData.genes: a cell array giving the gene names
%   proData.cells: a cell array, each cell giving cell attributes (such as cell type, culture condition, day captured)
if ~exist('minGenes','var') || isempty(minGenes)
    minGenes = 200;
end
if ~exist('minCells','var') || isempty(minCells)
    minCells = 3;
end
if ~exist('libararyflag','var') || isempty(libararyflag)
    libararyflag = 0;
end
if ~exist('logNormalize','var') || isempty(logNormalize)
    logNormalize = 0;
end
% disp('processing data:');
data0 = table2array(iniData); gene0 =  iniData.Properties.RowNames; cell0 = iniData.Properties.VariableNames;
%% filter cells that have expressed genes less than #200
dataTemp = data0;
dataTemp(data0 > 0) = 1;
msum = sum(dataTemp,1);
data0(:,msum < minGenes) = [];
cell0(msum < minGenes) = [];
%% filter genes that only express less than #3 cells
dataTemp = data0;
dataTemp(data0 > 0) = 1;
nsum = sum(dataTemp,2);
data0(nsum < minCells,:) = [];gene0(nsum < minCells) = [];%%

%% normalization:we employ a global-scaling normalization method that normalizes the gene expression measurements for each cell by the total expression 
% multiplies this by a scale factor (10,000 by default)
if libararyflag
    sM = sum(data0);
    data0 = data0./repmat(sM,size(data0,1),1)*10000;
end

if logNormalize
    data0 = log10(data0+1);
end
proData.data = data0;
proData.genes = gene0;
proData.cells = cell0;