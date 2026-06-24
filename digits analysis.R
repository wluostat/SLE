dat<-read.table("digits.txt",sep=",",header=FALSE)

### X is 64-dimensional and Y is a vector
x<-as.matrix(dat[,-65])
y<-c(dat[,65])


source("supervised Laplacian embedding.R")

### We first derive a bivariate reduced predictor
vj<-sleone_predict(x,y,y.type="categorical",x,hx=75,hy=1,0.001,2)
plot(vj[,1],vj[,2],col=y+1,
     xlab="",
     ylab="",
     main="2D reduced data for Digits",
     cex.main=2,cex.lab=1.6,tick.marks=FALSE)

### The proposed BIC-type criterion selects d=9
vm<-od(x,y,y.type="categorical",hx=75,hy=1,0.001,10) 


### Generate the curve of test misclassification error as 
###   the dimension of the reduced predictor varies
library(MASS)
library(ggplot2)
ss<-list()
n<-length(y)
for (j in 1:10)
{
 ss[[j]]<-numeric(0)
}
for (i in 1:10)
{
 si<-(1:n)[y==i-1]
 ssi<-sample(si,size=length(si),replace=FALSE)
 ni<-as.integer(length(si)/10)
 for (j in 1:10)
 { 
  if (j<10) {ss[[j]]<-c(ss[[j]],ssi[((j-1)*ni+1):(j*ni)])} else {
   ss[[j]]<-c(ss[[j]],ssi[(9*ni+1):(length(si))])}
 }
}
res<-matrix(0,10,10)
for (i in 1:10)
for (j in 1:10)
{
 xtr<-x[-ss[[j]],]
 ytr<-y[-ss[[j]]]

 rj<-sleone_predict(xtr,ytr,y.type="categorical",x,hx=75,hy=1,0.001,d=i)
 rj<-rj%*%diag(1/sqrt(diag(var(rj))))
 rjtr<-rj[-ss[[j]],]
 rjte<-rj[ss[[j]],]
 yte<-y[ss[[j]]]
 train_df <- data.frame(rjtr, group = ytr)
 lda_model <- lda(group ~ ., data = train_df)
 hyj <- predict(lda_model, data.frame(rjte))
 res[i,j]<-sum(hyj$class==yte)
}
a<- 1 - apply(res,1,sum)/1797
plot(1:10,a,xlab="Dimension of the reduced predictor",
     ylab="Test mis-classification rate",
     main="Performance of classification analysis",
     cex.main=2,cex.lab=1.6,pch=19)
lines(1:10,a)


### Now use d=3 for the reduced predictor and draw the 3-D scatter plot
vj<-sleone_predict(x,y,y.type="categorical",x,hx=75,hy=1,0.001,d=3)

library(scatterplot3d)
ind<-(1:n)[(y==1)|(y==2)|(y==8)|(y==9)]
scatterplot3d(vj[ind,1],vj[ind,3],vj[ind,2],angle=40,color=y[ind]+1,
     xlab="",
     ylab="",
     zlab="",
     main="3D reduced data for Digits (partial)",tick.marks=FALSE,
     cex.main=2,cex.lab=1.6)


