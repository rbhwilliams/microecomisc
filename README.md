# microecomisc
An R package containing a set of functions useful for basic, streamlined analysis of read count matrices encountered in microbial ecology. Input is a matrix that is assumed to have row and column names that link to annotations and metadata, respectively. Makes use of the an internal function from R/phyloseq called 'rarefaction_subsample'. There are functions for counting NA/zeros/non-detects, subsampling columns of a matrix to uniform depth, double centering a marix and perfomring PCA.
