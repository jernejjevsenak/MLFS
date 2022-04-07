Dear CRAN

this is the first submission of the R package MLFS. MLFS is designed to enable simulations of future forests based on machine learning algorithms. It is fully data-driven and as such does not need to be parameterized. It is applicable to all types of forests and runs on forest inventory data collected in almost every country in the world. Therefore, we believe it can make a significant contribution to forest and environmental modelling.

We are also preparing a publication to present our new model, which will be submitted in the next few days. Please note that some functions are very computationally intensive, so we have packaged them as 'dontrun', but they have all been tested in different environments and all examples have worked. We have also added some sample data that will be used to test different functions.

Best, 
Jernej


##  First submission
* This is the first submission of the package MLFS.

## Test environments
* local OS X install, R 4.1.1

* rhub Windows Server 2022 (https://builder.r-hub.io/status/original/MLFS_1.2.7.tar.gz-9a38a5c04747483dbd87c42b10ea73b5)
* rhub Ubuntu (https://builder.r-hub.io/status/original/MLFS_1.2.7.tar.gz-91812ea7bef541a8aedc9f2906b0d553)
* rhub Fedora Linux (https://builder.r-hub.io/status/original/MLFS_1.2.7.tar.gz-4f2ad668f9bd454ab2aac91736cf1fca)

* win-check oldrelease (https://win-builder.r-project.org/bE14YGPKjP0w/00check.log)
* win-check release (https://win-builder.r-project.org/o4LbYroGIzo0/00check.log)
* win-check devel (https://win-builder.r-project.org/AONlRRPVa49x/00check.log)

## R CMD check results
There were 0 ERRORs, 0 WARNINGs and 0 NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of MLFS
https://github.com/jernejjevsenak/MLFS/blob/master/revdep/checks.rds

All packages that we could install passed. 
