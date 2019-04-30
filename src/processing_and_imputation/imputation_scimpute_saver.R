## NOTE: We do not recommend running this yourself. The script can take >4 hours.

## Part 1: ScImpute
#install_github("Vivianstats/scImpute") # need to do this if not already installed
library(scImpute)

# just use the first patient as this is representative
# note: if running in R studio, may need to use setwd() command as this is not automatic
data_path = '../../data/imputation_intermediate/filtered_data/HIV1_Bld.csv'
num_clusters = 20
out_path = "../../data/imputation_intermediate/different_imputations/scimpute_intermediates/"

# Run scimpute
scimpute(# full path to raw count matrix
  count_path = data_path, 
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = out_path,           # full path to output directory
  labeled = FALSE,          # cell type labels not available
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = num_clusters,             # 20 cell subpopulations
  ncores = 10)              # number of cores used in parallel computation

## Part 2: Saver
#install.packages("SAVER") # need to do this if not already installed
library(SAVER)
#install.packages("doSNOW") # need to do this if not already installed
library(doSNOW)
#install.packages("parallel") # need to do this if not already installed
library(parallel)

# detect cores with parallel() package
nCores <- detectCores(logical = FALSE)
cat(nCores, "cores detected.")

# parallelize the computation
#install.packages("doParallel") # need to do this if not already installed
library(doParallel)
cl <- makeCluster(nCores) # fill in with your number of cores detected
registerDoParallel(cl)

# read in the filtered data
data = read.csv(data_path)
X <- data[, 2:ncol(data) ]

# perform the SAVER imputation
X.saver <- saver(X, ncores = 6)
str(X.saver$estimate)
out_name = "../../data/imputation_intermediate/different_imputations/saver_hiv1_bld.csv"
write.csv(X.saver$estimate, out_name)
