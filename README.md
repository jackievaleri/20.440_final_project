# Analysis of T Cell Populations in Blood and Cerebrospinal Fluid during Latent HIV-1 Infection 
This repository is for Problem Set 8 to replicate the figures, analyses, and scripts accompanying the final paper for 20.440. The content contained in this repository make up the majority of the supplementary figures, as they deal with imputation method comparisons as well as biomarker analysis.

## Purpose 
Human immunodeficiency virus (HIV)-associated neurological disorders (HAND) such as dementia, myelopathy, and sensory neuropathies are well-documented, but molecular mechanisms underlying their onset is not fully understood. Understanding sources of potential neural degeneration caused by aberrant immune activation could allow for development of therapies/strategies to reduce HAND. To understand the populations of immune cells that are activated in HIV+ vs. HIV- patients, single cell RNA-seq (scRNA-seq) datasets were processed and genes with zero expression values were imputed. Two figures were created surrounding this preprocessing and imputation: the first, a comparison of different imputation methods; the second, a visualization of biomarkers superimposed on the imputed data from one imputation method, called MAGIC.

## To Reproduce the Analysis
### Directory Structure
The structure of this repository is based off of Claire Duvallet's Aerodigestion Aspiration Analysis public [repository](https://github.com/cduvallet/aspiration-analysis-public/blob/master/README.md). 

``` 
|-- README.md
|
|-- data
|    |-- raw_data: raw data from GEO, Farhadian 2018
|    |-- imputation_intermediate
|    |    |-- different_imputations 
|    |    |     |-- SAVER and scImpute final files 
|    |    |     |-- scimpute_intermediates: intermediate files created running scImpute
|    |    |-- magic_all_samples: files after MAGIC imputation was run on them
|    |-- filtered_data: data for all samples after preprocessing completed
|
|-- output
|    |-- figures: png files used to assemble final figures
|    
|-- src
|    |-- processing_and_imputation: notebooks and scripts required to process and impute the files
|    |-- figure_making: jupyter notebooks to create supplementary figures 1 and 2
```

### Scripts to Reproduce
The scripts and notebooks can be run in any order as they are all self-contained. 

#### Reproducing Processing and Imputation
To rerun MAGIC imputation and filtering of the data, the `magic_imputation.ipynb` notebook in the `src/processing_and_imputation` folder can be run. This notebook will preprocess the data according to the preprocessing steps outlined below in Methods and then perform MAGIC on the samples as outlined below in the imputation sections of the methods. The output to this notebook will appear in `data/filtered_data` for the preprocessing steps and `data/imputation_intermediate/magic_all_samples` for the MAGIC steps, although these folders are already pre-populated for convenience, as the code takes a while to run.

To rerun scImpute and SAVER imputation of the filtered data, the `imputation_scimpute_saver.R` script in the `src/processing_and_imputation` folder can be run. This script will perform scImpute and SAVER on the samples as outlined below in the imputation sections of the methods. The output to this notebook will appear in `data/imputation_intermediate/different_imputations` for the final results of the scImpute and SAVER steps and `data/imputation_intermediate/different_imputations/scimpute_intermediates` for the intermediate scImpute steps, although these folders are already pre-populated for convenience.

#### Reproducing Figures
To regenerate Supplementary Figure 1 to visualize comparisons of the imputation methods, the `supp_figure_1_imputation_comparison.ipynb` in the `src/figure_making` folder can be run. This script will run PCA and t-SNEs on all three imputation methods and the raw filtered data to visualize the effects of different imputations. The output to this notebook will appear in the `src/output/figures` folder, although this folder is prepopulated for convenience.

To regenerate Supplementary Figure 2 to visualize biomarkers on the MAGIC-imputed data, the `supp_figure_2_magic_biomarkers.ipynb` in the `src/figure_making` folder can be run. This script will run a t-SNE on a representative sample imputed by MAGIC and superimpose the t-SNE plots with 9 biomarkers identified in Farhadian 2018. The output to this notebook will appear in the `src/output/figures` folder, although this folder is prepopulated for convenience.

### Installation Requirements
You will need Python v.3.7.1, although some packages may have backwards compatibility. Python analyses require the following packages: magic v.1.5.3, scprep v.0.11.1, matplotlib v.3.0.3, numpy v.1.15.4, pandas v.0.23.4, csv v.1.0, and sklearn v.0.20.1.

You will also need R v.3.5.1, although again some packages may have backwards compatibility. R analyses require the following packages: scImpute v.0.0.9, SAVER v.1.1.1, doSNOW v.1.0.16, and doParallel v.1.0.14. In the R script, the installation lines are commented out but can easily be added depending on the packages your machine already has installed.

# import statements
import magic
import scprep

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
%matplotlib inline

import gzip
import csv

from sklearn.manifold import TSNE
from scprep.io.csv import load_csv

## Information about the Code
### Data Sources
To understand the characteristic gene expression features of HIV+ and HIV- patients, single cell RNA sequencing (scRNA-seq) data of blood and cerebrospinal fluid (CSF) cells were obtained using SeqWell methodology<sup>1</sup>.  Datasets from [Farhadian, et. al](https://insight.jci.org/articles/view/121718). included CSF scRNA-seq data from 3 HIV+ and 2 HIV- patients; for 2 of the HIV+ patients, matched blood scRNA-seq data was included. The patients were males between the ages of 33 and 59, with varied histories of substance use, alcohol use, and smoking, matched between the patient and control populations. This data is also available on GEO associated with the ID [GEO177397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117397).

### Methods
#### Preprocessing
Several preprocessing steps were utilized via scprep before imputing the data. First, empty columns and rows were removed from the expression matrix. Any cells that expressed less than 500 genes were retained for the analysis, as this step is standard to remove low-quality cells. Genes expressed in less than 3 cells were filtered out to remove extremely rare genes that might confound the analysis, again as is standard in the field. The gene expression matrix was normalized by library size via L1 normalization, meaning each cell was rescaled so that the sum of expression values for each cell sums to 1. This rescales the data as if each cell was sampled evenly. Finally, data was transformed by taking the square root of each element in the gene expression matrix.

#### Imputation with SAVER, scImpute, and MAGIC
To decrease the sparsity of the gene expression matrix prior to downstream analyses, several imputation-based techniques have been published that uncover likely values for dropout genes. Three methods, SAVER<sup>2</sup>, scImpute<sup>3</sup>, and MAGIC<sup>3</sup>, were compared. MAGIC was used with the default settings, using the built-in MAGIC principal component analysis (PCA) and t-distributed stochastic neighborhood embedding (t-SNE) functions to visualize the results, again with the default settings (decay = 15, k nearest neighbors = 5, distance = Euclidean, number of principal components = 100, time step t = 10). scImpute was used with an estimation of 20 clusters, as Farhadian 2018 identified 14 clusters but with some clusters rather large and heterogeneous, we felt it might be too restrictive to use only 14. 20 was chosen because this would ensure that we at least capture the original clusters as well as possibly split the one mixed population and the four large T cell clusters. All other default settings were utilized (labeled cells = False, dropout threshold = 0.5). SAVER was run with the default settings and the computation was parallelized with the doParallel function on six cores to match the number on the server used. To visualize scImpute, SAVER, and raw results with PCA and t-SNE, the sklearn functions were used with default settings.

## References
1. Farhadian, SF, Mehta, SS, Zografou, C. Single-cell RNA sequencing reveals microglia-like cells in cerebrospinal fluid during virologically suppressed HIV. JCI Insight. 2018; 3(18): 121718. doi: 10.1172/jci.insight.121718.

2. Huang, M, Wang, J, Torre, E, et al. SAVER: gene expression recovery for single-cell RNA sequencing. Nat Methods. 2018; 15: 539-542. doi: 10.1038/s41592-018-0033-z.

3. Li, WV, Li, JJ. An accurate and robust imputation method scImpute for single-cell RNA-seq data. Nat Communications. 2018; 9(997). doi: 10.1038/s41467-018-03405-7.

4. van Dijk, D, Sharma, R, Nainys, J, et al. Recovering gene interactions from single-cell data using data diffusion. Cell. 2018; 174(3): 716-729. doi: 10.1016/j.cell.2018.05.061.  
