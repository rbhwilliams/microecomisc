#--functions for performing no-frills PCA--
#--written by RW (this file dated 150916)

#--function internal functions for computing variance related measures
#--'svdObj' is an output from 'svd'.
compute.variance.per.pc<-function(svdObj){(svdObj$d^2)/nrow(svdObj$u)}
compute.cum.per.total.var<-function(svdObj){cumsum(svdObj$d^2)/sum(svdObj$d^2)*100}

#--quick PCA for doubly centered matrices...assumes column-means of X  will be zero
quick.pca<-function(X)
{
 res<-list(scores=NULL,var.per.pc=NULL,cum.prop.var.across.pc=NULL)
 matd<-svd(X)
 res$scores<-matd$u%*%diag(matd$d)
 rownames(res$scores)<-rownames(X)
 res$var.per.pc<-compute.variance.per.pc(matd)
 res$cum.prop.var.across.pc<-compute.cum.per.total.var(matd)
 return(res)
}
