##################################################################################################
# rarefaction subsample function, one sample
# RW  on 27 August: this function taken from transform_filter-methods.R in R/phyloseq and modified
##################################################################################################
#' @keywords internal
rarefaction_subsample<-function(x,sample.size,replace=FALSE)
{
 #--Create replacement species vector
 rarvec <- numeric(length(x))
 
 #--Perform the sub-sampling. Suppress warnings due to old R compat issue.
 #--Also, make sure to avoid errors from x summing to zero,
 #--and there are no observations to sample.
 #--The initialization of rarvec above is already sufficient.
 if(sum(x)<=0)
 {
  #--Protect against, and quickly return an empty vector,
  #--if x is already an empty count vector
  return(rarvec)
 }
 if(replace)
 {
  #--then resample with replacement
  suppressWarnings(subsample<-sample(1:length(x),sample.size,replace=TRUE,prob=x))
 }
 else
 {
  #--resample without replacement
  obsvec<-apply(data.frame(OTUi=1:length(x),times=x),1,function(x){rep_len(x["OTUi"],x["times"])})
  obsvec<-unlist(obsvec,use.names=FALSE)
  #--use `sample` for subsampling. Hope that obsvec doesn't overflow.
  suppressWarnings(subsample<-sample(obsvec,sample.size,replace=FALSE))
 }
 #--Tabulate the results (these are already named by the order in `x`)
 sstab<-table(subsample)
 #--Assign the tabulated random subsample values to the species vector
 rarvec[as(names(sstab),"integer")]<-sstab
 #--Return abundance vector. Let replacement happen elsewhere.
 return(rarvec)
}

