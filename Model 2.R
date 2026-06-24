library(MASS)
library(energy)
source("supervised Laplacian embedding.R")

N=500
n=200
p=10
eps=0.001
rec<-matrix(0,N,3)
for (i in 1:N)
 {
  x<-matrix(rnorm(n*p),n,p)
  y<-abs(x[,1])+abs(x[,2]+1)+rnorm(n)/3

  m<-n*5
  xte<-matrix(rnorm(m*p),m,p)
  rxte<-abs(xte[,1])+abs(xte[,2]+1)

  hx<-16
  hy<-1
  vm<-od(x,y,y.type="continuous",hx,hy,eps,10)
  rec[i,1]<-vm
  v2<-sleone_predict(x,y,y.type="continuous",xte,hx,hy,eps,vm)
  rec[i,2]<-dcor(v2,rxte)
  v3<-slesec_predict(x,y,y.type="continuous",xte,hx,hy,eps,vm)
  rec[i,3]<-dcor(v3,rxte)
 }



