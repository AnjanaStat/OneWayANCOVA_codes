rm(list=ls())
library(matrixcalc)
install.packages("Matrix")
library(Matrix)
library(smoothmest)
n1=10;n2=10;n3=10;n4=25;b1=0.1;b2=0.2;b3=0.3;b4=1;p=1;v1=0.5;v2=0.8;v3=1;v4=1

fun1<-function(n1,n2,n3,n4,p,v1,v2,v3,v4)
{
  #X1=rnorm(n1,2,1);X2=rnorm(n2,2,1);X3=rnorm(n3,2,1);X4=rnorm(n4,2,1)
  h1=10/(n1-1);h2=10/(n2-1);h3=5/(n3-1);h4=(2/(n4-1))
  X1=seq(-5,5,h1);X2=seq(-5,5,h2);X3=seq(0,5,h3);X4=seq(-4,-2,h4)
  Y1=rnorm(n1,p*X1,sqrt(v1));Y2=rnorm(n2,p*X2,sqrt(v2))
  Y3=rnorm(n3,p*X3,sqrt(v3));Y4=rnorm(n4,p*X4,sqrt(v4))
  y1=mean(Y1);y2=mean(Y2);y3=mean(Y3);y4=mean(Y4)
  x1=mean(X1);x2=mean(X2);x3=mean(X3);x4=mean(X4)
  p1=sum((X1-x1)*(Y1-y1));p2=sum((X2-x2)*(Y2-y2));p3=sum((X3-x3)*(Y3-y3));p4=sum((X4-x4)*(Y4-y4))
  q1=sum((Y1-y1)^2);q2=sum((Y2-y2)^2);q3=sum((Y3-y3)^2);q4=sum((Y4-y4)^2)
  r1=sum((X1-x1)^2);r2=sum((X2-x2)^2);r3=sum((X3-x3)^2);r4=sum((X4-x4)^2)
  S1=q1/(n1-2)-(p1^2)/((n1-2)*r1);S2=q2/(n2-2)-(p2^2)/((n2-2)*r2)
  S3=q3/(n3-2)-(p3^2)/((n3-2)*r3);S4=q4/(n4-2)-(p4^2)/((n4-2)*r4)
  bnu=p1+p2+p3+p4
  bdno=r1+r2+r3+r4
  b=bnu/bdno
  a1=y1-b*x1;a2=y2-b*x2;a3=y3-b*x3;a4=y4-b*x4
  phinu=S1*r1+S2*r2+S3*r3+S4*r4
  phidn=(r1+r2+r3+r4)^2
  phi=phinu/phidn
  V1=(S1/n1 +S2/n2 +phi*(x1-x2)^2);V2=(S2/n2 +S3/n3 +phi*(x2-x3)^2);V3=(S3/n3 +S4/n4 +phi*(x3-x4)^2)
  T1=(a2-a1)/sqrt(V1);T2=(a3-a2)/sqrt(V2);T3=(a4-a3)/sqrt(V3)
  T=min(T1,T2,T3)
  return(T)
}

fun2<-function(n1,n2,n3,n4,p,v1,v2,v3,v4)
{
  x<-replicate(1000,fun1(n1,n2,n3,n4,p,v1,v2,v3,v4))
  y<-sort(x,decreasing=FALSE)
  c<-y[950]
  return(c)
}

fun3<-function(n1,n2,n3,n4,mu1,mu2,mu3,mu4,p,v1,v2,v3,v4)
{
  #X1=rnorm(n1,2,1);X2=rnorm(n2,2,1);X3=rnorm(n3,2,1);X4=rnorm(n4,2,1)
  h1=10/(n1-1);h2=10/(n2-1);h3=5/(n3-1);h4=(2/(n4-1))
  X1=seq(-5,5,h1);X2=seq(-5,5,h2);X3=seq(0,5,h3);X4=seq(-4,-2,h4)
  #Y1=rdoublex(n1,mu1+p*X1,sqrt(v1)/sqrt(2))
  #Y2=rdoublex(n2,mu2+p*X2,sqrt(v2)/sqrt(2))
  #Y3=rdoublex(n3,mu3+p*X3,sqrt(v3)/sqrt(2))
  #Y4=rdoublex(n4,mu4+p*X4,sqrt(v4)/sqrt(2))
  #g1=rt(n1,5,ncp=0)/sqrt(5/3);g2=rt(n2,5,ncp=0)/sqrt(5/3)
  #g3=rt(n3,5,ncp=0)/sqrt(5/3);g4=rt(n4,5,ncp=0)/sqrt(5/3)
  #Y1=sqrt(v1)*g1+mu1+p*X1;Y2=sqrt(v2)*g2+mu2+p*X2
  #Y3=sqrt(v3)*g3+mu3+p*X3;Y4=sqrt(v4)*g4+mu4+p*X4
  #g1=rexp(n1,1);g2=rexp(n2,1);g3=rexp(n3,1);g4=rexp(n4,1)
  #Y1=mu1+p*X1+sqrt(v1)*(g1-1);Y2=mu2+p*X2+sqrt(v2)*(g2-1)
  #Y3=mu3+p*X3+sqrt(v3)*(g3-1);Y4=mu4+p*X4+sqrt(v4)*(g4-1)
  #g1=(rlnorm(n1,0,1)-exp(.5))/sqrt(exp(2)-exp(1));g2=(rlnorm(n2,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #g3=(rlnorm(n3,0,1)-exp(.5))/sqrt(exp(2)-exp(1));g4=(rlnorm(n4,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y1=mu1+p*X1+sqrt(v1)*g1;Y2=mu2+p*X2+sqrt(v2)*g2;Y3=mu3+p*X3+sqrt(v3)*g3;Y4=mu4+p*X4+sqrt(v4)*g4
  #g1=(rweibull(n1,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);g2=(rweibull(n2,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #g3=(rweibull(n3,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);g4=(rweibull(n4,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y1=mu1+p*X1+sqrt(v1)*g1;Y2=mu2+p*X2+sqrt(v2)*g2;Y3=mu3+p*X3+sqrt(v3)*g3;Y4=mu4+p*X4+sqrt(v4)*g4
  Y1=rnorm(n1,mu1+p*X1,sqrt(v1));Y2=rnorm(n2,mu2+p*X2,sqrt(v2))
  Y3=rnorm(n3,mu3+p*X3,sqrt(v3));Y4=rnorm(n4,mu4+p*X4,sqrt(v4))
  y1=mean(Y1);y2=mean(Y2);y3=mean(Y3);y4=mean(Y4)
  x1=mean(X1);x2=mean(X2);x3=mean(X3);x4=mean(X4)
  p1=sum((X1-x1)*(Y1-y1));p2=sum((X2-x2)*(Y2-y2));p3=sum((X3-x3)*(Y3-y3));p4=sum((X4-x4)*(Y4-y4))
  q1=sum((Y1-y1)^2);q2=sum((Y2-y2)^2);q3=sum((Y3-y3)^2);q4=sum((Y4-y4)^2)
  r1=sum((X1-x1)^2);r2=sum((X2-x2)^2);r3=sum((X3-x3)^2);r4=sum((X4-x4)^2)
  S1=q1/(n1-2)-(p1^2)/((n1-2)*r1);S2=q2/(n2-2)-(p2^2)/((n2-2)*r2)
  S3=q3/(n3-2)-(p3^2)/((n3-2)*r3);S4=q4/(n4-2)-(p4^2)/((n4-2)*r4)
  bnu=p1+p2+p3+p4
  bdno=(r1+r2+r3+r4)
  b=bnu/bdno
  a1=y1-b*x1;a2=y2-b*x2;a3=y3-b*x3;a4=y4-b*x4
  phinu=S1*r1+S2*r2+S3*r3+S4*r4
  phidn=(r1+r2+r3+r4)^2
  phi=phinu/phidn
  V1=(S1/n1 +S2/n2 +phi*(x1-x2)^2);V2=(S2/n2 +S3/n3 +phi*(x2-x3)^2);V3=(S3/n3 +S4/n4 +phi*(x3-x4)^2)
  T1=(a2-a1)/sqrt(V1);T2=(a3-a2)/sqrt(V2);T3=(a4-a3)/sqrt(V3)
  T=min(T1,T2,T3)
  out<-fun2(n1,n2,n3,n4,b,S1,S2,S3,S4)
  a=0
  if(T>out)
    a=a+1
  return(a)
}
fun4<-function(n1,n2,n3,n4,mu1,mu2,mu3,mu4,p,v1,v2,v3,v4)
{
  # find number of times maxt-T value > critical value among 1000 values
  out<-replicate(1000,fun3(n1,n2,n3,n4,mu1,mu2,mu3,mu4,p,v1,v2,v3,v4))
  p<-sum(out)/1000
  return(p)
}


c<-replicate(4,fun4(8,9,7,5,0,0,0,0,4,1,1,1,1))
p=mean(c);p
c<-replicate(4,fun4(8,9,7,5,0,0,0,0,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(8,9,7,5,0,0,0,0,4,4,3,2,1))
p=mean(c);p
c<-replicate(4,fun4(8,9,7,5,0,0,0,0,4,4,4,4,4))
p=mean(c);p
c<-replicate(4,fun4(8,9,7,5,0,0,0,0,4,0.2,0.3,0.8,1))
p=mean(c);p
c<-replicate(4,fun4(8,9,7,5,0,0,0,0,4,2,3,8,1))
p=mean(c);p
c<-replicate(4,fun4(8,9,7,5,0,0,0,0,4,1,2,0.8,1))
p=mean(c);p





c<-replicate(4,fun4(25,9,30,5,0,0,0,0,4,1,1,1,1))
p=mean(c);p
c<-replicate(4,fun4(25,9,30,5,0,0,0,0,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(25,9,30,5,0,0,0,0,4,4,3,2,1))
p=mean(c);p
c<-replicate(4,fun4(25,9,30,5,0,0,0,0,4,4,4,4,4))
p=mean(c);p
c<-replicate(4,fun4(25,9,30,5,0,0,0,0,4,0.2,0.3,0.8,1))
p=mean(c);p
c<-replicate(4,fun4(25,9,30,5,0,0,0,0,4,2,3,8,1))
p=mean(c);p
c<-replicate(4,fun4(25,9,30,5,0,0,0,0,4,1,2,0.8,1))
p=mean(c);p




c<-replicate(4,fun4(10,10,10,10,0,0,0,0,4,1,1,1,1))
p=mean(c);p
c<-replicate(4,fun4(10,10,10,10,0,0,0,0,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(10,10,10,10,0,0,0,0,4,4,3,2,1))
p=mean(c);p
c<-replicate(4,fun4(10,10,10,10,0,0,0,0,4,4,4,4,4))
p=mean(c);p
c<-replicate(4,fun4(10,10,10,10,0,0,0,0,4,0.2,0.3,0.8,1))
p=mean(c);p
c<-replicate(4,fun4(10,10,10,10,0,0,0,0,4,2,3,8,1))
p=mean(c);p
c<-replicate(4,fun4(10,10,10,10,0,0,0,0,4,1,2,0.8,1))
p=mean(c);p





c<-replicate(4,fun4(25,20,15,10,0,0,0,0,4,1,1,1,1))
p=mean(c);p
c<-replicate(4,fun4(25,20,15,10,0,0,0,0,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(25,20,15,10,0,0,0,0,4,4,3,2,1))
p=mean(c);p
c<-replicate(4,fun4(25,20,15,10,0,0,0,0,4,4,4,4,4))
p=mean(c);p
c<-replicate(4,fun4(25,20,15,10,0,0,0,0,4,0.2,0.3,0.8,1))
p=mean(c);p
c<-replicate(4,fun4(25,20,15,10,0,0,0,0,4,2,3,8,1))
p=mean(c);p
c<-replicate(4,fun4(25,20,15,10,0,0,0,0,4,1,2,0.8,1))
p=mean(c);p




c<-replicate(4,fun4(10,15,20,25,0,0,0,0,4,1,1,1,1))
p=mean(c);p
c<-replicate(4,fun4(10,15,20,25,0,0,0,0,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(10,15,20,25,0,0,0,0,4,4,3,2,1))
p=mean(c);p
c<-replicate(4,fun4(10,15,20,25,0,0,0,0,4,4,4,4,4))
p=mean(c);p
c<-replicate(4,fun4(10,15,20,25,0,0,0,0,4,0.2,0.3,0.8,1))
p=mean(c);p
c<-replicate(4,fun4(10,15,20,25,0,0,0,0,4,2,3,8,1))
p=mean(c);p
c<-replicate(4,fun4(10,15,20,25,0,0,0,0,4,1,2,0.8,1))
p=mean(c);p






c<-replicate(1,fun4(20,25,30,35,0,0,0,0,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1,1.2,1.2,1.6,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*1.2,1.2*1.2,1.2*1.2,1.6*1.2,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*1.4,1.2*1.4,1.2*1.4,1.6*1.4,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*1.6,1.2*1.6,1.2*1.6,1.6*1.6,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*1.8,1.2*1.8,1.2*1.8,1.6*1.8,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*2,1.2*2,1.2*2,1.6*2,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*2.2,1.2*2.2,1.2*2.2,1.6*2.2,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*2.4,1.2*2.4,1.2*2.4,1.6*2.4,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*2.6,1.2*2.6,1.2*2.6,1.6*2.6,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*2.8,1.2*2.8,1.2*2.8,1.6*2.8,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*3,1.2*3,1.2*3,1.6*3,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*3.2,1.2*3.2,1.2*3.2,1.6*3.2,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*3.4,1.2*3.4,1.2*3.4,1.6*3.4,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*3.6,1.2*3.6,1.2*3.6,1.6*3.6,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*3.8,1.2*3.8,1.2*3.8,1.6*3.8,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*4,1.2*4,1.2*4,1.6*4,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*4.2,1.2*4.2,1.2*4.2,1.6*4.2,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*4.4,1.2*4.4,1.2*4.4,1.6*4.4,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*4.6,1.2*4.6,1.2*4.6,1.6*4.6,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*4.8,1.2*4.8,1.2*4.8,1.6*4.8,4,1,2,3,4))
p=mean(c);p
c<-replicate(1,fun4(20,25,30,35,1*5,1.2*5,1.2*5,1.6*5,4,1,2,3,4))
p=mean(c);p





c<-replicate(4,fun4(20,25,30,20,0,0,0,0,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(20,25,30,20,1,1.1,1.2,1.3,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(20,25,30,20,1*1.5,1.1*1.5,1.2*1.5,1.3*1.5,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(20,25,30,20,1*2,1.1*2,1.2*2,1.3*2,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(20,25,30,20,1*2.5,1.1*2.5,1.2*2.5,1.3*2.5,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(20,25,30,20,1*3,1.1*3,1.2*3,1.3*3,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(20,25,30,20,1*3.5,1.1*3.5,1.2*3.5,1.3*3.5,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(20,25,30,20,1*4,1.1*4,1.2*4,1.3*4,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(20,25,30,20,1*4.5,1.1*4.5,1.2*4.5,1.3*4.5,4,1,2,3,4))
p=mean(c);p
c<-replicate(4,fun4(20,25,30,20,1*5,1.1*5,1.2*5,1.3*5,4,1,2,3,4))
p=mean(c);p



