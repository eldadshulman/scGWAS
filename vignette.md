Example for trait-trajectory analysis
================

This is an example of the analysis performed for a dataset of pancreatic
islet development from [Byrnes et
al](https://doi.org/10.1038/s41467-018-06176-3).

The Monocle object is available to download from figshare.com:
[link](https://figshare.com/articles/dataset/Monocle_Objects_-_V2_Dataset/6783554?backTo=/collections/Lineage_dynamics_of_murine_pancreatic_development_at_single-cell_resolution/4158458).

Load the following required packages:

``` r
library(monocle)
library(Seurat)
library(data.table)
```

Source the functions from this github directory:

``` r
source('functions_GWAS_traj.R')
```

### Preprocessing

First, we extract the expression matrix and metadata, including the
pseudotime of the cells. We Load the Monocole object:

``` r
load('E14_fev_lineage_monocle_ob.Rdata')
```

We extract the count matrix and cell metadata from this object:

``` r
# count matrix
exp <- HSMM_seur_var@assayData$exprs

# metadata
meta <- HSMM_seur_var@phenoData@data
```

Might be a good idea to remove the monocle object and save the new
objects:

``` r
rm('HSMM_seur_var')
saveRDS(object = exp, file = "panc.exp.RDS")
saveRDS(object = meta, file = "panc.meta.RDS")
```

#### Processing the count matrix

We keep genes expressed in at least 10 cells

``` r
genes.keep <- rowSums(as.matrix(exp) > 0) >= 10
exp <- exp[genes.keep,]
```

We normalize the data using Seurat’s log-normalization:

``` r
exp <- Seurat::NormalizeData(exp)
```

The following converts mouse genes to human orthologues:

``` r
exp <- conv_hs(exp)
```

#### Calculating cell-trait association scores

This part involves running MAGMA’s gene property analysis for each cell
separately. This might take long, so consider using multiple cores.
First, we generate “covariate files” for MAGMA for each cell. Each
output file is a table with the columns: gene, normalized expression in
the cell, and average normalized expression in the dataset.

``` r
generate_covs_cells(exp = exp, cors = 20, wd = ".")
```

The next step, it running the gene property analysis, for each cell.
This analysis firs the following the regression model to each cell:

![](https://github.com/eldadshulman/scGWAS/blob/master/eq.PNG)

where Z is the vector of the gene’s Z-score converted from the p-values
obtained from MAGMA’s gene analysis for the trait. B is a matrix of
technical confounders, including gene length and SNPs’ LDs, calculated
by MAGMA. C is the vector of normalized expression of the genes in the
cell, and A is a vector of the average normalized expression for the
gene in the dataset. The t-statistic of \(\beta_c\) is taken as the
score for the association between the cell and trait.

Note that this step takes even longer.

``` r
gene_prop_cells(cor = 30)
```

Last, extract the cell scores from MAGMA’s output, namely the
t-statistics.

``` r
tss <- get_scores(cors = 20)
```

#### Calculating trajectory-trait association
