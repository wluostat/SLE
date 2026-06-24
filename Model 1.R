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
  sig<-diag(1,p)
  delta<-0.5
  for (k in 1:(p-1))
   {
    sig[k,k+1]<-delta
    sig[k+1,k]<-delta
   }
  es<-eigen(sig)
  hsig<-es$vectors%*%diag(sqrt(es$values))%*%t(es$vectors)
  x<-x%*%hsig
  px<-exp(3*(x[,1]-x[,2]))/(1+exp(3*(x[,1]-x[,2])))
  y<-rbinom(n,size=1,prob=px)


  m<-n*5
  xte<-matrix(rnorm(m*p),m,p)
  xte<-xte%*%hsig
  pxte<-exp(3*(xte[,1]-xte[,2]))/(1+exp(3*(xte[,1]-xte[,2])))
  yte<-rbinom(m,size=1,prob=pxte)

  d=1
  rxte<-pxte

  vm<-od(x,y,y.type="categorical",16,1,eps,10)
  rec[i,1]<-vm
  v1<-sleone_predict(x,y,y.type="categorical",xte,16,1,eps,vm)
  rec[i,2]<-dcor(v1,pxte)
  v3<-slesec_predict(x,y,y.type="categorical",xte,16,1,eps,vm)
  rec[i,3]<-dcor(v3,pxte)
 }
apply(rec,2,mean)


#par(mfrow=c(2,2))
#vm<-sletai(x,y,hx,hy,eps,1)
#plot(x[,1],Re(vm))
#vg<-gsir(x,y,hx,hy,eps,1)
#plot(x[,1],vg)
#v1<-sleone(x,y,hx,hy,eps,1)
#plot(x[,1],v1)
#v2<-slesec(x,y,hx,hy,1)
#plot(x[,1],v2)


