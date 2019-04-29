# Analysis of T Cell Populations in Blood and Cerebrospinal Fluid during Latent HIV-1 Infection 
This repository is for Problem Set 8 to replicate the figures, analyses, and scripts accompanying the final paper for 20.440. The content contained in this repository make up the majority of the supplementary figures, as they deal with imputation method comparisons as well as initial data analysis with Seurat.

## Purpose 
Human immunodeficiency virus (HIV)-associated neurological disorders (HAND) such as dementia, myelopathy, and sensory neuropathies are well-documented, but molecular mechanisms underlying their onset is not fully understood. Understanding sources of potential neural degeneration caused by aberrant immune activation could allow for development of therapies/strategies to reduce HAND. To understand the populations of immune cells that are activated in HIV+ vs. HIV- patients, single cell RNA-seq (scRNA-seq) datasets were processed and clustered. Differential expression analysis was conducted on subclusters of the samples representing T cells to understand differences in gene expression among diseased and healthy individuals, especially how these differences relate to latently infected CD4 T cells and inflammation.

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

## Information about the Code
### Data Sources
To understand the characteristic gene expression features of HIV+ and HIV- patients, single cell RNA sequencing (scRNA-seq) data of blood and cerebrospinal fluid (CSF) cells were obtained using SeqWell methodology<sup>1</sup>.  Datasets from [Farhadian, et. al](https://insight.jci.org/articles/view/121718). included CSF scRNA-seq data from 3 HIV+ and 2 HIV- patients; for 2 of the HIV+ patients, matched blood scRNA-seq data was included. The patients were males between the ages of 33 and 59, with varied histories of substance use, alcohol use, and smoking, matched between the patient and control populations. This data is also available on GEO associated with the ID [177397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117397).

### Methods
#### Imputation with SAVER, scImpute, and MAGIC
To decrease the sparsity of the gene expression matrix prior to downstream analyses, several imputation-based techniques have been published that uncover likely values for dropout genes. Three methods, SAVER<sup>2</sup>, scImpute<sup>3</sup>, and MAGIC<sup>3</sup>, were compared. MAGIC was used with the default settings, using the built-in MAGIC principal component analysis (PCA) and t-distributed stochastic neighborhood embedding (t-SNE) functions to visualize the results, again with the default settings (decay = 15, k nearest neighbors = 5, distance = Euclidean, number of principal components = 100, time step t = 10). scImpute was used with an estimation of 20 clusters, as Farhadian 2018 identified 14 clusters but with some clusters rather large and heterogeneous, we felt it might be too restrictive to use only 14. 20 was chosen because this would ensure that we at least capture the original clusters as well as possibly split the one mixed population and the four large T cell clusters. All other default settings were utilized (labeled cells = False, dropout threshold = 0.5). SAVER was run with the default settings and the computation was parallelized with the doParallel function on six cores to match the number on the server used. To visualize scImpute, SAVER, and raw results with PCA and t-SNE, the sklearn functions were used with default settings.



1. Farhadian, SF, Mehta, SS, Zografou, C. Single-cell RNA sequencing reveals microglia-like cells in cerebrospinal fluid during virologically suppressed HIV. JCI Insight. 2018; 3(18): 121718. doi: 10.1172/jci.insight.121718.

2. Huang, M, Wang, J, Torre, E, et al. SAVER: gene expression recovery for single-cell RNA sequencing. Nat Methods. 2018; 15: 539-542. doi: 10.1038/s41592-018-0033-z.

3. Li, WV, Li, JJ. An accurate and robust imputation method scImpute for single-cell RNA-seq data. Nat Communications. 2018; 9(997). doi: 10.1038/s41467-018-03405-7.

4. van Dijk, D, Sharma, R, Nainys, J, et al. Recovering gene interactions from single-cell data using data diffusion. Cell. 2018; 174(3): 716-729. doi: 10.1016/j.cell.2018.05.061.  
