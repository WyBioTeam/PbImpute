## PbImpute
### Overall architecture of the PbImpute algorithm pipeline:

![Sample Image](model.png)
## 1. Introduction

PbImpute: Precise Zero Discrimination and Balanced Imputation in Single-Cell RNA Sequencing Data

## 2. Run (Please place the entire project in "D:\PbImpute")

### 2.1:&nbsp;&nbsp;Step 1: preprocessing (Please run `D:\PbImpute\preprocessing\PbImpute.mlx`)

The first step is preprocessing. Please run the MATLAB file PbImpute.mlx, which contains detailed execution steps.
Please note that in the Ledein_SNN function, Rscript = 'D:\R\R-4.3.1\bin\Rscript.exe' specifies the path to your R installation. If the operation is successful, there will be a "grout.txt" file in the D:\PbImpute\preprocessing\.

### 2.2:&nbsp;&nbsp;Step 2: ZINB imputation (Please run `D:\PbImpute\ZINB_imputation\pbimpute1.R`)
The second step is  ZINB imputation, specifically including our a new ZINB model and static repair.

### 2.2:&nbsp;&nbsp;Step 3: Residual imputation (Please run `D:\PbImpute\Residual_imputation\test.py`)
The third  step is  Residual imputation, specifically including our a graphic embedding and dynamic repair.

## 3. Comparative Models
The following models and methods are referenced:
- **scImpute**: [Nature Communications](https://www.nature.com/articles/s41467-018-03405-7)
- **SAVER**: [Nature Methods](https://www.nature.com/articles/s41592-018-0033-z)
- **MAGIC**: [Cell](https://www.cell.com/cell/fulltext/S0092-8674(18)30724-4)
- **DCA**: [Nature Communications](https://www.nature.com/articles/s41467-018-07931-2)
- **DeepImpute**: [Genome Biology](https://link.springer.com/article/10.1186/s13059-019-1837-6)
- **bayNorm**: [Bioinformatics](https://academic.oup.com/bioinformatics/article/36/4/1174/5581401)
- **GE-Impute**: [Bioinformatics](https://academic.oup.com/bib/article/23/5/bbac313/6651303?login=false)
- **ALRA**: [Nature Communications](https://www.nature.com/articles/s41467-021-27729-z)
- **CL-Impute**: [Journal of Computational Biology](https://www.sciencedirect.com/science/article/abs/pii/S001048252300728X)
- **TsImpute**: [Bioinformatics](https://academic.oup.com/bioinformatics/article/39/12/btad731/7457483)

## 4. Data availability
### 4.1 Real dataset
The datasets were derived from publicly available sources: 
- The PBMC datasets from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k),
- The worm neuron cells from [Cole Trapnell Lab](https://cole-trapnell-lab.github.io/worm-rna/docs/),
- The LPS datasets from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17721),
- The mouse bladder cells from [Figshare](https://figshare.com/s/865e694ad06d5857db4b).
### 4.2 Simulate dataset
The simulated datasets come from: [SCRABBLE](https://link.springer.com/article/10.1186/s13059-019-1681-8) and [Splatter](https://link.springer.com/article/10.1186/s13059-017-1305-0).  
Since SCRABBLE is constrained by bulk RNA-seq data and introduces other datasets, we will not compare the SCRABBLE method here.
- **Dataset 1**, **Dataset 2**, **Dataset 3**:
  - Each dataset corresponds to dropout rates of 35%, 43%, 57%, 71%, and 83% dropout dataset.
  - Compute the average of the results with the same dropout rate.
- **Large Dataset**:
  - To test runtime and memory usage (from 10,000 cells to 50,000 cells)
  - To visualize clustering (40,000 cells)





