############################################################################
###
### This file includes the following functions:
### sleone: the proposed supervised Laplacian embedding (SLE) that produces 
###         the coefficient matrix V
### sleone_predict: same as above but produces the reduced predictor for any
###         set of observations (of X) 
### slesec: an alternative of the proposed SLE mentioned at the end of 
###         Section 2, again with output V 
### slesec_predict: same as above but produces the reduced predictor for any
###         set of observations (of X) 
### tune_hx: the tuning procedure for the bandwidth hx used in the Gaussian 
###         kernel of X, with the optimal bandwidth as the output
### tune_hy: the tuning procedure for the bandwidth hy used in the Gaussian 
###         kernel of Y if Y is continuous, with output being the optimal hy
### od: the proposed BIC-type criterion for determining the dimension of the
###         the reduced predictor, with the dimension being the output
###
############################################################################
library(MASS)
library(energy)
library(cdcsis)

#########################################################################
###
### The proposed supervised Laplacian embedding with output V
### Input:
###   x is X, a n times p matrix
###   y is Y, a n-dimensional vector
###   y.type is continuous/categorical, set to continuous by default
###   hx is the bandwidth used in Gaussian kernel of X
###   hy is the bandwidth in Gaussian kernel of Y, if Y is continuous
###   epsilon is the scalar used to approximate the MP-inverse of G_X 
###   d is the dimension of reduced predictor (see "od" below for BIC)
###
### Output:
###   V is the coeffcient matrix for the reduced predictor U = G_X V
###
########################################################################
sleone<-function(x,y,y.type,hx,hy,epsilon,d)
{
 if(missing(y.type)) {y.type="continuous"}   
   ## regard Y as continuous if not prespecified

 n<-length(y)
 W<-diag(0,n)
 K<-W
 for (i in 1:n)
 for (j in 1:n)
  {
   if (y.type=="continuous") {W[i,j]<-exp(-(y[i]-y[j])^2/(2*hy^2))}
   if (y.type=="categorical") {W[i,j]<-I(y[i]==y[j])}
   K[i,j]<-exp(-sum((x[i,]-x[j,])^2)/(2*hx^2))
  }
 Q<-diag(1,n)-matrix(1/n,n,n)
 G<-Q%*%K%*%Q
 W<-Q%*%W%*%Q
   ## G is the gram matrix G_X, W is the same W as in existing SLE methods

 Ginv<-ginv(G+epsilon*diag(1,n))
 V<-Ginv%*%G%*%svd(Ginv%*%G%*%W%*%G%*%Ginv)$u[,1:d]
 return(V)
}

#########################################################################
###
### The proposed supervised Laplacian embedding with output U for any data
### Input:
###   x is X, a n times p matrix
###   y is Y, a n-dimensional vector
###   y.type is continuous/categorical, set to continuous by default
###   xte is X in test set, a m times p matrix, may equal X if needed
###   hx is the bandwidth used in Gaussian kernel of X
###   hy is the bandwidth in Gaussian kernel of Y, if Y is continuous
###   epsilon is the scalar used to approximate the MP-inverse of G_X 
###   d is the dimension of reduced predictor (see "od" below for BIC)
###
### Output:
###   U is the reduced predictor in test set associated with input Xte
###
########################################################################
sleone_predict<-function(x,y,y.type,xte,hx,hy,epsilon,d)
{
 if(missing(y.type)) {y.type="continuous"}
   ## regard Y as continuous if not prespecified

 n<-length(y)
 W<-diag(0,n)
 K<-W
 for (i in 1:n)
 for (j in 1:n)
  {
   if (y.type=="continuous") {W[i,j]<-exp(-(y[i]-y[j])^2/(2*hy^2))}
   if (y.type=="categorical") {W[i,j]<-I(y[i]==y[j])}
   K[i,j]<-exp(-sum((x[i,]-x[j,])^2)/(2*hx^2))
  }
 Q<-diag(1,n)-matrix(1/n,n,n)
 G<-Q%*%K%*%Q
 W<-Q%*%W%*%Q
   ## G is the gram matrix G_X, W is the same W as in existing SLE methods

 Ginv<-ginv(G+epsilon*diag(1,n))
 V<-svd(Ginv%*%G%*%W%*%G%*%Ginv)$u[,1:d]

 m<-nrow(xte)
 Kte<-matrix(0,m,n)
 for (i in 1:m)
 for (j in 1:n)
  {Kte[i,j]<-exp(-sum((xte[i,]-x[j,])^2)/(2*hx^2))}
 Gte<-Kte%*%Q
 U<-Gte%*%Ginv%*%V
 return(U)
}

#########################################################################
###
### An alternative of the proposed supervised Laplacian embedding with 
### output V, using a different normalization constraint 
### Input:
###   x is X, a n times p matrix
###   y is Y, a n-dimensional vector
###   y.type is continuous/categorical, set to continuous by default
###   hx is the bandwidth used in Gaussian kernel of X
###   hy is the bandwidth in Gaussian kernel of Y, if Y is continuous
###   d is the dimension of reduced predictor (see "od" below for BIC)
###
### Output:
###   V is the coeffcient matrix for the reduced predictor U = G_X V
###
########################################################################
slesec<-function(x,y,y.type,hx,hy,d)
{
 if(missing(y.type)) {y.type="continuous"}
   ## regard Y as continuous if not prespecified

 n<-length(y)
 W<-diag(0,n)
 K<-W
 for (i in 1:n)
 for (j in 1:n)
  {
   if (y.type=="continuous") {W[i,j]<-exp(-(y[i]-y[j])^2/(2*hy^2))}
   if (y.type=="categorical") {W[i,j]<-I(y[i]==y[j])}
   K[i,j]<-exp(-sum((x[i,]-x[j,])^2)/(2*hx^2))
  }
 Q<-diag(1,n)-matrix(1/n,n,n)
 G<-Q%*%K%*%Q
 W<-Q%*%W%*%Q
   ## G is the gram matrix G_X, W is the same W as in existing SLE methods

 eG<-eigen(G)
 hG<-eG$vectors%*%diag(sqrt(abs(eG$values)))%*%t(eG$vectors)

 V<-W%*%hG%*%svd(hG%*%W%*%hG)$u[,1:d]
 return(V)
}

#########################################################################
###
### An alternative of the proposed supervised Laplacian embedding with 
### output V, using a different normalization constraint 
### Input:
###   x is X, a n times p matrix
###   y is Y, a n-dimensional vector
###   y.type is continuous/categorical, set to continuous by default
###   xte is X in test set, a m times p matrix, may equal X if needed
###   hx is the bandwidth used in Gaussian kernel of X
###   hy is the bandwidth in Gaussian kernel of Y, if Y is continuous
###   epsilon is the scalar used to approximate the MP-inverse of G_X 
###   d is the dimension of reduced predictor (see "od" below for BIC)
###
### Output:
###   U is the reduced predictor in test set associated with input Xte
###
########################################################################
slesec_predict<-function(x,y,y.type,xte,hx,hy,epsilon,d)
{
 if(missing(y.type)) {y.type="continuous"}
 n<-length(y)
 W<-diag(0,n)
 K<-W
 for (i in 1:n)
 for (j in 1:n)
  {
   if (y.type=="continuous") {W[i,j]<-exp(-(y[i]-y[j])^2/(2*hy^2))}
   if (y.type=="categorical") {W[i,j]<-I(y[i]==y[j])}
   K[i,j]<-exp(-sum((x[i,]-x[j,])^2)/(2*hx^2))
  }
 Q<-diag(1,n)-matrix(1/n,n,n)
 G<-Q%*%K%*%Q
 W<-Q%*%W%*%Q

 eG<-eigen(G)
 hG<-eG$vectors%*%diag(sqrt(abs(eG$values)))%*%t(eG$vectors)

 V<-W%*%hG%*%svd(hG%*%W%*%hG)$u[,1:d]

 m<-nrow(xte)
 Kte<-matrix(0,m,n)
 for (i in 1:m)
 for (j in 1:n)
  {Kte[i,j]<-exp(-sum((xte[i,]-x[j,])^2)/(2*hx^2))}
 Gte<-Kte%*%Q
 Ginv<-ginv(G+epsilon*diag(1,n))

 U<-Gte%*%Ginv%*%V
 return(U)
}

########################################################################
###
###  The tuning procedure of the bandwidth hx for Gaussian kernel of X,
###  using k-fold cross-validation and the distance correlation between
###  the reduced predictor and Y as the measure of goodness of fit.
###  Input:
###       k is the number of folds in CV, default at 10
###       the other inputs are the same as "sleone" above
###
###  Output:
###       the optimal bandwidth hx
###
########################################################################
tune_hx<-function(x,y,y.type,hy,epsilon,d,k)
{ 
 if(missing(k)) {k=10}
 if(missing(y.type)) {y.type="continuous"}

 n<-nrow(x)
 p<-ncol(x)
 hx<-(2^c((-3):4))*sqrt(mean(diag(var(x))))
 res<-matrix(0,k,length(hx))
 for (i in 1:k)
  {
   ind<-c((as.integer(n/k)*(i-1)+1):(as.integer(n/k)*i))
   xtr<-x[-ind,]
   ytr<-y[-ind]
   xte<-x[ind,]
   yte<-y[ind]
   for (j in 1:length(hx))
    {   
     m<-nrow(xte)
     uij<-sletai_predict(xtr,ytr,y.type="continuous",xte,hx[j],hy,epsilon,d)
     res[i,j]<-dcor(uij,c(yte))
    }
   }
 return(hx[which.max(apply(res,2,mean))])
}

########################################################################
###
###  The tuning procedure of the bandwidth hy for Gaussian kernel of Y if
###  Y is continuous, using k-fold cross-validation and the distance 
###  correlation between the reduced predictor and Y as the measure of 
###  goodness of fit.
###  Input:
###       k is the number of folds in CV, default at 10
###       the other inputs are the same as "sleone" above
###
###  Output:
###       the optimal bandwidth hy
###
########################################################################
tune_hy<-function(x,y,y.type,hx,epsilon,d,k)
{ 
 if(missing(k)) {k=10}
 if(missing(y.type)) {y.type="continuous"}

 n<-nrow(x)
 p<-ncol(x)
 hy<-(2^c((-3):4))*sd(y)
 res<-matrix(0,k,length(hy))
 for (i in 1:k)
  {
   ind<-c((as.integer(n/k)*(i-1)+1):(as.integer(n/k)*i))
   xtr<-x[-ind,]
   ytr<-y[-ind]
   xte<-x[ind,]
   yte<-y[ind]
   for (j in 1:length(hy))
    {   
     m<-nrow(xte)
     uij<-sleone_predict(xtr,ytr,y.type="continuous",xte,hx,hy[j],epsilon,d)
     res[i,j]<-dcor(uij,c(yte))
    }
   }
 return(hy[which.max(apply(res,2,mean))])
}

##########################################################################
### 
### The proposed BIC-type criterion for determining the dimension of the
###         the reduced predictor
### Input:
###      dmax is a strict upper bound of the dimension of the reduced 
###         predictor, set to 10 by default
###      The other inputs are the same as for "sleone" above
###
### Output:
###      The dimension of the reduced predictor 
###
##########################################################################
od<-function(x,y,y.type,hx,hy,epsilon,dmax)
{ 
 if (missing(dmax)) {dmax=10}
 if(missing(y.type)) {y.type="continuous"}

 n<-nrow(x)
 p<-ncol(x)
 dmax<-min(p,10)
 res<-rep(0,dmax)

 U<-sleone_predict(x,y,y.type=y.type,x,hx,hy,eps,dmax)
 tau<-rep(0,dmax)
 tau[1]<-dcor(U[,1],c(y))
 for (i in 2:dmax) {tau[i]<-cdcor(U[,i],y,U[,1:(i-1)],width=1)$statistic}
 
 mt<-tau
 for (i in 1:dmax) {mt[i]<-max(tau[i:dmax])}
 return(which.min(mt+0.5*c(1:dmax)/n^(1/3))-1)
}






