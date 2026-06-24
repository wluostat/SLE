library(MASS)
library(energy)
source("supervised Laplacian embedding.R")

N=500
n=400
p=10
eps=0.001
rec<-matrix(0,N,5)
for (i in 1:N)
 {
  x<-matrix((runif(n*p)-0.5)*4,n,p)
  y<-I(x[,1]+rnorm(n)/3>0) + I(x[,2]+rnorm(n)/3>0)

  m<-n*5
  xte<-matrix((runif(m*p)-0.5)*4,m,p)

  d=2
  rxte<-xte[,1:2]

  vm<-od(x,y,y.type="categorical",2,1,eps,dmax=10)
  rec[i,1]<- vm
  v1<-sleone_predict(x,y,y.type="categorical",xte,16,1,eps,vm)
  rec[i,2]<-dcor(v1,rxte)
  v3<-slesec_predict(x,y,y.type="categorical",xte,16,1,eps,vm)
  rec[i,3]<-dcor(v3,rxte)
 }
apply(rec,2,mean)




