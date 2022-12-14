#Program Spline
data3<-read.table("E:/1/Data & Syntax/Data_IHS_JII.txt",header=TRUE)
y1 <- data3$IHSG
y2 <- data3$JII
x1 <- data3$Inflasi
x2 <- data3$IDJ
plot(y1~x1)
plot(y1~x2)
plot(y2~x1)
plot(y2~x2)

mpi<-function(x,eps=1.47925e-17)
{
  x<-as.matrix(x)
  xsvd<-svd(x)
  diago<-xsvd$d[xsvd$d>eps]
  if (length(diago)==1)
  {
    xplus<-as.matrix(xsvd$v[,1])%*%t(as.matrix(xsvd$u[,1])/diago)
  }else
  {
    xplus<-xsvd$v[,1:length(diago)]%*% diag(1/diago)%*%t(xsvd$u[,1:length(diago)])
  }
  	return(xplus)
}

k21<-2.2
k22<-2.4
k23<-3.0
l <- 0.01 # lambda
n<-length(y2)
  trun<-function(data,knots,power)
  {((data-knots)^power)*(data>=knots)}
  m<-matrix(0,ncol=5,nrow=n)
  m[,1]<-1
  m[,2]<-x1
  m[,3]<-trun(x1,k21,1)
  m[,4]<-trun(x1,k22,1)
  m[,5]<-trun(x1,k23,1)
  D<-matrix(0,ncol=5,nrow=5)
  D[,1]<-c(0,0,0,0,0)
  D[,2]<-c(0,0,0,0,0)
  D[,3]<-c(0,0,1,0,0)
  D[,4]<-c(0,0,0,1,0)
  D[,5]<-c(0,0,0,0,1)
  beta<-solve(t(m)%*%m+(l*D))%*%t(m)%*%y2
  flamda<-m%*%beta
  alamda<- m%*% solve(t(m)%*%m)%*%t(m)
  ia<-diag(n)-alamda
  bgcv<-(sum(diag(ia))/n)^2
  agcv<-t(ia%*%y2)%*%(ia%*%y2)/n^2
  gcv<-agcv/bgcv
  galat<-y2-flamda
  ybar<-sum(y2)/n
  sse<-t(y2-flamda)%*%(y2-flamda)
  syy<- t(y2-ybar)%*%(y2-ybar)
  ssr<- t(flamda-ybar)%*%( flamda-ybar)
  koef.determinasi<-ssr/syy
  mse<-as.vector(sse)/(n-5)
  msr<- as.vector(ssr)/4
  fhitung<-msr/mse
  pvalue=pf(fhitung,4,(n-5),lower.tail = FALSE)
  covbeta<-solve(t(m)%*%m)*mse
  seb1<-sqrt(covbeta[1,1])
  seb2<-sqrt(covbeta[2,2])
  seb3<-sqrt(covbeta[3,3])
  seb4<-sqrt(covbeta[4,4])
  seb5<-sqrt(covbeta[5,5])
  tb1<-beta[1]/seb1
  tb2<-beta[2]/seb2
  tb3<-beta[3]/seb3
  tb4<-beta[4]/seb4
  tb5<-beta[5]/seb5
  i1<-seq(min(x1),max(x1),length=n)
fest1<-beta[1]+beta[2]*i1+beta[3]*trun(i1,k21,1)+beta[4]*trun(i1,k22,1)+beta[5]*trun(i1,k23,1)
win.graph()
plot(x1,y2,xlim=c(1,4),ylim=c(500,800), xlab="Inflasi", ylab="JII")
par(new=T)
plot(i1,fest1,type="l",xlim=c(1,4),ylim=c(500,800),xlab="Inflasi", ylab="JII")

beta
gcv
fhitung
pvalue