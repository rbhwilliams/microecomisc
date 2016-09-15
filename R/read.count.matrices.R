#--some functions for working with read count matrices--
#--written by Rohan Williams (this file from 150916)

#--for a read count matrix 'X', convert NA elements to 0
convert.na2zero<-function(X)
{
 res<-X
 res[is.na(X)]<-0
 return(res)
}

#--and count number of zero elements in each row of a matrix 'X'
count.zeros.per.row<-function(X)
{
 Xnz<-X
 Xnz[X<1e-15]<-NA
 Xnzc<-rowSums(is.na(Xnz))
 return(Xnzc)
}

#--a function to double centre a matrix (useful for visualisation of time effects using PCA)
double.center.data.matrix<-function(X)
{
 #--compute a double centered version of 'mat'
 nr<-nrow(X)
 nc<-ncol(X)
 Jr<-matrix(1,nr,nr)
 Jc<-matrix(1,nc,nc)
 res<-(X-(X%*%Jc)/nc-(Jr%*%X)/nr+(Jr%*%X%*%Jc)/(nr*nc))
 rownames(res)<-rownames(X)
 colnames(res)<-colnames(X)
 return(res)
}

