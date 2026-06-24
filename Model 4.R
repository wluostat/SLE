library(MASS)
library(energy)
source("supervised Laplacian embedding.R")

N=500
n=400
p=10
eps=0.001
hx=9.6
hy=3.2
rec<-matrix(0,N,5)
for (i in 1:N)
 {
  x<-matrix(runif(n*p)*4,n,p)-2
  sig<-diag(1,p)
  sig[1,2]<-0
  sig[2,1]<-0
  es<-eigen(sig)
  hsig<-es$vectors%*%diag(sqrt(es$values))%*%t(es$vectors)
  x<-x%*%hsig
  y<-x[,1]^2+(2*rbinom(n,size=1,prob=0.5)-1)*x[,2]^2+rnorm(n)/3

  m<-n*5
  xte<-matrix(runif(m*p)*4,m,p)-2
  rxte<-cbind(xte[,1]^2,xte[,2]^2)

  vm<-od(x,y,y.type="continuous",hx,hy,eps,10)
  rec[i,1]<-vm
  v1<-sleone_predict(x,y,y.type="continuous",xte,hx,hy,eps,vm)
  rec[i,2]<-dcor(v1,rxte)
  v3<-slesec_predict(x,y,y.type="continuous",xte,hx,hy,eps,vm)
  rec[i,3]<-dcor(v3,rxte)
 }

apply(rec,2,mean)



