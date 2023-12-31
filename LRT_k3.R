rm(list=ls())
library(Iso)
library(smoothmest)
n1=18;n2=18;n3=18;n4=20;b=1;v1=1;v2=1.5;v3=2;v4=4;mu1=mu2=mu3=1

fun1<-function(n1,n2,n3,b,v1,v2,v3)
{
  #Y1=rnorm(n1,2,1);Y2=rnorm(n2,2,1);Y3=rnorm(n3,2,1);Y4=rnorm(n4,2,1)
  #Y1=runif(n1, 0, 5);Y2=runif(n2,0,5);Y3=runif(n3,0,5)
  h1=10/(n1-1);h2=5/(n2-1);h3=2/(n3-1)
  Y1=seq(-5,5,h1);Y2=seq(0,5,h2);Y3=seq(2,4,h3)
  X1=rnorm(n1,b*Y1,sqrt(v1));X2=rnorm(n2,b*Y2,sqrt(v2));X3=rnorm(n3,b*Y3,sqrt(v3))
  X1.=mean(X1);X2.=mean(X2);X3.=mean(X3)
  Y1.=mean(Y1);Y2.=mean(Y2);Y3.=mean(Y3)
  p1=sum((X1-X1.)*(Y1-Y1.));p2=sum((X2-X2.)*(Y2-Y2.));p3=sum((X3-X3.)*(Y3-Y3.))
  q1=sum((Y1-Y1.)^2);q2=sum((Y2-Y2.)^2);q3=sum((Y3-Y3.)^2)
  r1=sum((X1-X1.)^2);r2=sum((X2-X2.)^2);r3=sum((X3-X3.)^2)
  S1=r1/(n1-2)-(p1^2)/((n1-2)*q1);S2=r2/(n2-2)-(p2^2)/((n2-2)*q2)
  S3=r3/(n3-2)-(p3^2)/((n3-2)*q3)
  w1=n1/S1;w2=n2/S2;w3=n3/S3;w=w1+w2+w3
  Xw.=(w1*X1.+w2*X2.+w3*X3.)/w;Yw.=(w1*Y1.+w2*Y2.+w3*Y3.)/w
  b0n=(sum((X1-Xw.)*(Y1-Yw.))/S1+sum((X2-Xw.)*(Y2-Yw.))/S2+sum((X3-Xw.)*(Y3-Yw.))/S3)
  b0d=(sum((Y1-Yw.)^2)/S1)+(sum((Y2-Yw.)^2)/S2)+(sum((Y3-Yw.)^2)/S3)
  B0=b0n/b0d;A0=Xw.-B0*Yw.
  x13=A0
  x15=B0
  S01=sum((X1-A0-B0*Y1)^2)/n1;S02=sum((X2-A0-B0*Y2)^2)/n2;S03=sum((X3-A0-B0*Y3)^2)/n3
  c(S01,S02,S03)
  repeat
  {
    w1=n1/S01;w2=n2/S02;w3=n3/S03;w=w1+w2+w3
    X0w.=(w1*X1.+w2*X2.+w3*X3.)/w;Y0w.=(w1*Y1.+w2*Y2.+w3*Y3.)/w
    b00n=(sum((X1-Xw.)*(Y1-Yw.))/S01+sum((X2-Xw.)*(Y2-Yw.))/S02+sum((X3-Xw.)*(Y3-Yw.))/S03)
    b00d=(sum((Y1-Yw.)^2)/S01)+(sum((Y2-Yw.)^2)/S02)+(sum((Y3-Yw.)^2)/S03)
    B00=b00n/b00d;A00=X0w.-B00*Y0w.
    D=abs(x13-A00);E=abs(x15-B00)
    if(D<=0.000001 & E<=0.00001)
    {
      break
    }
    S001=sum((X1-A00-B00*Y1)^2)/n1;S002=sum((X2-A00-B00*Y2)^2)/n2;S003=sum((X3-A00-B00*Y3)^2)/n3
    x13<-A00
    x15<-B00
    S01=S001;S02=S002;S03=S003
  }
  P1=sum((X1-X1.)*(Y1-Y1.))+sum((Y2-Y2.)*(X2-X2.))+sum((Y3-Y3.)*(X3-X3.))
  Q1=sum((Y1-Y1.)^2)+sum((Y2-Y2.)^2)+sum((Y3-Y3.)^2)
  beta0=P1/Q1
  alpha0=c(X1.-beta0*Y1., X2.-beta0*Y2., X3.-beta0*Y3.)
  s1=sum((X1-alpha0[1]-beta0*Y1)^2)/n1;s2=sum((X2-alpha0[2]-beta0*Y2)^2)/n2
  s3=sum((X3-alpha0[3]-beta0*Y3)^2)/n3
  x0=alpha0
  x8<-alpha0
  x12=beta0
  repeat
  {
    W1=c(n1/s1,n2/s2,n3/s3)
    alphan=pava(x0, W1)
    pn=(sum(Y1*X1)/s1+sum(Y2*X2)/s2+sum(Y3*X3)/s3)-
      (sum(Y1*alphan[1])/s1+sum(Y2*alphan[2])/s2+sum(Y3*alphan[3])/s3)
    qn=sum(Y1^2)/s1+sum(Y2^2)/s2+sum(Y3^2)/s3
    betan=pn/qn
    x9=c(X1.-betan*Y1., X2.-betan*Y2., X3.-betan*Y3.)
    varn1=sum((X1-alphan[1]-betan*Y1)^2)/n1
    varn2=sum((X2-alphan[2]-betan*Y2)^2)/n2
    varn3=sum((X3-alphan[3]-betan*Y3)^2)/n3
    s1<-varn1;s2<-varn2;s3<-varn3
    x10=abs(x8-alphan)
    M=max(x10[1], x10[2], x10[3],na.rm = FALSE);N=abs(x12-betan)
    if(M<=0.0000001& N<=0.0000001)
    {
      break
    }
    x8<-alphan
    x0=x9
    x12=betan
  }
  value=((s1/S01)^(n1/2))*((s2/S02)^(n2/2))*((s3/S03)^(n3/2))
  return(value)
}

X<-replicate(1000,fun1(n1=18,n2=18,n3=18,b=1,v1=10,v2=20,v3=30))
plot(X)

fun2<-function(n1,n2,n3,b,v1,v2,v3)
{
  x<-replicate(500,fun1(n1,n2,n3,b,v1,v2,v3))
  y<-sort(x,decreasing=FALSE)
  c<-y[25]
  return(c)
}

fun3<-function(n1,n2,n3,a1,a2,a3,b,v1,v2,v3)
{
  #Y1=rnorm(n1,2,1);Y2=rnorm(n2,2,1);Y3=rnorm(n3,2,1)
  #Y1=runif(n1, 0, 5);Y2=runif(n2,0,5);Y3=runif(n3,0,5)
  h1=10/(n1-1);h2=5/(n2-1);h3=2/(n3-1)
  Y1=seq(-5,5,h1);Y2=seq(0,5,h2);Y3=seq(2,4,h3)
  X1=rnorm(n1,a1+b*Y1,sqrt(v1));X2=rnorm(n2,a2+b*Y2,sqrt(v2));X3=rnorm(n3,a3+b*Y3,sqrt(v3))
  #X1=rdoublex(n1,a1+b*Y1,sqrt(v1)/sqrt(2))
  #X2=rdoublex(n2,a2+b*Y2,sqrt(v2)/sqrt(2))
  #X3=rdoublex(n3,a3+b*Y3,sqrt(v3)/sqrt(2))
  #g1=rt(n1,5,ncp=0)/sqrt(5/3);g2=rt(n2,5,ncp=0)/sqrt(5/3);g3=rt(n3,5,ncp=0)/sqrt(5/3)
  #X1=sqrt(v1)*g1+a1+b*Y1;X2=sqrt(v2)*g2+a2+b*Y2;X3=sqrt(v3)*g3+a3+b*Y3
  #g1=rexp(n1,1);g2=rexp(n2,1);g3=rexp(n3,1)
  #X1=a1+b*Y1+sqrt(v1)*(g1-1);X2=a2+b*Y2+sqrt(v2)*(g2-1);X3=a3+b*Y3+sqrt(v3)*(g3-1)
  #g1=(rlnorm(n1,0,1)-exp(.5))/sqrt(exp(2)-exp(1));g2=(rlnorm(n2,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #g3=(rlnorm(n3,0,1)-exp(.5))/sqrt(exp(2)-exp(1));g4=(rlnorm(n4,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y1=mu1+p*X1+sqrt(v1)*g1;Y2=mu2+p*X2+sqrt(v2)*g2;Y3=mu3+p*X3+sqrt(v3)*g3;Y4=mu4+p*X4+sqrt(v4)*g4
  #g1=(rweibull(n1,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);g2=(rweibull(n2,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #g3=(rweibull(n3,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #X1=a1+b*Y1+sqrt(v1)*g1;X2=a2+b*Y2+sqrt(v2)*g2;X3=a3+b*Y3+sqrt(v3)*g3
  #Y1=rnorm(n1,b1+p*X1,sqrt(v1));Y2=rnorm(n2,b2+p*X2,sqrt(v2))
  #Y3=rnorm(n3,b3+p*X3,sqrt(v3));Y4=rnorm(n4,b4+p*X4,sqrt(v4))
  #d1=sqrt(v1)*((rchisq(n1,df=8)-8)/4);d2=sqrt(v2)*((rchisq(n2,df=8)-8)/4)
  #d3=sqrt(v3)*((rchisq(n3,df=8)-8)/4)
  #X1=a1+b*Y1+d1; X2=a2+b*Y2+d2; X3=a3+b*Y3+d3
  X1.=mean(X1);X2.=mean(X2);X3.=mean(X3)
  Y1.=mean(Y1);Y2.=mean(Y2);Y3.=mean(Y3)
  p1=sum((X1-X1.)*(Y1-Y1.));p2=sum((X2-X2.)*(Y2-Y2.));p3=sum((X3-X3.)*(Y3-Y3.))
  q1=sum((Y1-Y1.)^2);q2=sum((Y2-Y2.)^2);q3=sum((Y3-Y3.)^2)
  r1=sum((X1-X1.)^2);r2=sum((X2-X2.)^2);r3=sum((X3-X3.)^2)
  S1=r1/(n1-2)-(p1^2)/((n1-2)*q1);S2=r2/(n2-2)-(p2^2)/((n2-2)*q2)
  S3=r3/(n3-2)-(p3^2)/((n3-2)*q3)
  w1=n1/S1;w2=n2/S2;w3=n3/S3;w=w1+w2+w3
  Xw.=(w1*X1.+w2*X2.+w3*X3.)/w;Yw.=(w1*Y1.+w2*Y2.+w3*Y3.)/w
  b0n=(sum((X1-Xw.)*(Y1-Yw.))/S1+sum((X2-Xw.)*(Y2-Yw.))/S2+sum((X3-Xw.)*(Y3-Yw.))/S3)
  b0d=(sum((Y1-Yw.)^2)/S1)+(sum((Y2-Yw.)^2)/S2)+(sum((Y3-Yw.)^2)/S3)
B0=b0n/b0d;A0=Xw.-B0*Yw.
x13=A0
x15=B0
S01=sum((X1-A0-B0*Y1)^2)/n1;S02=sum((X2-A0-B0*Y2)^2)/n2;S03=sum((X3-A0-B0*Y3)^2)/n3
c(S01,S02,S03)
repeat
{
  w1=n1/S01;w2=n2/S02;w3=n3/S03;w=w1+w2+w3
  X0w.=(w1*X1.+w2*X2.+w3*X3.)/w;Y0w.=(w1*Y1.+w2*Y2.+w3*Y3.)/w
  b00n=(sum((X1-Xw.)*(Y1-Yw.))/S01+sum((X2-Xw.)*(Y2-Yw.))/S02+sum((X3-Xw.)*(Y3-Yw.))/S03)
  b00d=(sum((Y1-Yw.)^2)/S01)+(sum((Y2-Yw.)^2)/S02)+(sum((Y3-Yw.)^2)/S03)
  B00=b00n/b00d;A00=X0w.-B00*Y0w.
  D=abs(x13-A00);E=abs(x15-B00)
  if(D<=0.000001 & E<=0.00001)
  {
    break
  }
  S001=sum((X1-A00-B00*Y1)^2)/n1;S002=sum((X2-A00-B00*Y2)^2)/n2;S003=sum((X3-A00-B00*Y3)^2)/n3
  x13<-A00
  x15<-B00
  S01=S001;S02=S002;S03=S003
}
P1=sum((X1-X1.)*(Y1-Y1.))+sum((Y2-Y2.)*(X2-X2.))+sum((Y3-Y3.)*(X3-X3.))
Q1=sum((Y1-Y1.)^2)+sum((Y2-Y2.)^2)+sum((Y3-Y3.)^2)
beta0=P1/Q1
alpha0=c(X1.-beta0*Y1., X2.-beta0*Y2., X3.-beta0*Y3.)
s1=sum((X1-alpha0[1]-beta0*Y1)^2)/n1;s2=sum((X2-alpha0[2]-beta0*Y2)^2)/n2
s3=sum((X3-alpha0[3]-beta0*Y3)^2)/n3
x0=alpha0
x8<-alpha0
x12=beta0
repeat
{
  W1=c(n1/s1,n2/s2,n3/s3)
  alphan=pava(x0, W1)
  pn=(sum(Y1*X1)/s1+sum(Y2*X2)/s2+sum(Y3*X3)/s3)-
    (sum(Y1*alphan[1])/s1+sum(Y2*alphan[2])/s2+sum(Y3*alphan[3])/s3)
  qn=sum(Y1^2)/s1+sum(Y2^2)/s2+sum(Y3^2)/s3
  betan=pn/qn
  x9=c(X1.-betan*Y1., X2.-betan*Y2., X3.-betan*Y3.)
  varn1=sum((X1-alphan[1]-betan*Y1)^2)/n1
  varn2=sum((X2-alphan[2]-betan*Y2)^2)/n2
  varn3=sum((X3-alphan[3]-betan*Y3)^2)/n3
  s1<-varn1;s2<-varn2;s3<-varn3
  x10=abs(x8-alphan)
  M=max(x10[1], x10[2], x10[3],na.rm = FALSE);N=abs(x12-betan)
  if(M<=0.0000001& N<=0.0000001)
  {
    break
  }
  x8<-alphan
  x0=x9
  x12=betan
}
value=((s1/S01)^(n1/2))*((s2/S02)^(n2/2))*((s3/S03)^(n3/2))
out<-fun2(n1,n2,n3,B0,S1,S2,S3)
a=0
if(value<out)
  a=a+1
return(a)
  #value=-2*log(LR)
  #return(value)
}

fun4<-function(n1,n2,n3,b1,b2,b3,p,v1,v2,v3)
{
  # find number of times maxt-T value > critical value among 1000 values
  out<-replicate(1000,fun3(n1,n2,n3,b1,b2,b3,p,v1,v2,v3))
  p<-sum(out)/1000
  return(p)
}

c<-replicate(5,fun4(10,15,20,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1,1.3,1.8,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*1.2,1.3*1.2,1.8*1.2,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*1.4,1.3*1.4,1.8*1.4,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*1.6,1.3*1.6,1.8*1.6,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*1.8,1.3*1.8,1.8*1.8,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2,1.3*2,1.8*2,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2.2,1.3*2.2,1.8*2.2,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2.4,1.3*2.4,1.8*2.4,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2.6,1.3*2.6,1.8*2.6,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2.8,1.3*2.8,1.8*2.8,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3,1.3*3,1.8*3,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3.2,1.3*3.2,1.8*3.2,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3.4,1.3*3.4,1.8*3.4,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3.6,1.3*3.6,1.8*3.6,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3.8,1.3*3.8,1.8*3.8,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4,1.3*4,1.8*4,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4.2,1.3*4.2,1.8*4.2,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4.4,1.3*4.4,1.8*4.4,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4.6,1.3*4.6,1.8*4.6,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4.8,1.3*4.8,1.8*4.8,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*5,1.3*5,1.8*5,4,1,1,1))
c;crit_value=mean(c);crit_value

c<-replicate(5,fun4(10,15,20,0,0,0,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1,1.3,1.8,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*1.2,1.3*1.2,1.8*1.2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*1.4,1.3*1.4,1.8*1.4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*1.6,1.3*1.6,1.8*1.6,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*1.8,1.3*1.8,1.8*1.8,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2,1.3*2,1.8*2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2.2,1.3*2.2,1.8*2.2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2.4,1.3*2.4,1.8*2.4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2.6,1.3*2.6,1.8*2.6,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*2.8,1.3*2.8,1.8*2.8,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3,1.3*3,1.8*3,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3.2,1.3*3.2,1.8*3.2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3.4,1.3*3.4,1.8*3.4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3.6,1.3*3.6,1.8*3.6,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*3.8,1.3*3.8,1.8*3.8,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4,1.3*4,1.8*4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4.2,1.3*4.2,1.8*4.2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4.4,1.3*4.4,1.8*4.4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4.6,1.3*4.6,1.8*4.6,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*4.8,1.3*4.8,1.8*4.8,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,1*5,1.3*5,1.8*5,4,1,2,3))
c;crit_value=mean(c);crit_value


c<-replicate(5,fun4(20,30,25,0,0,0,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,1,1.1,1.2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,1.2,1.1*1.2,1.2*1.2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,1.4,1.1*1.4,1.2*1.4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,1.6,1.1*1.6,1.2*1.6,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,1.8,1.1*1.8,1.2*1.8,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,2,1.1*2,1.2*2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,2.2,1.1*2.2,1.2*2.2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,2.4,1.1*2.4,1.2*2.4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,2.6,1.1*2.6,1.2*2.6,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,2.8,1.1*2.8,1.2*2.8,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,3,1.1*3,1.2*3,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,3.2,1.1*3.2,1.2*3.2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,3.4,1.1*3.4,1.2*3.4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,3.6,1.1*3.6,1.2*3.6,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,3.8,1.1*3.8,1.2*3.8,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,4,1.1*4,1.2*4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,4.2,1.1*4.2,1.2*4.2,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,4.4,1.1*4.4,1.2*4.4,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,4.6,1.1*4.6,1.2*4.6,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,4.8,1.1*4.8,1.2*4.8,4,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,30,25,5,1.1*5,1.2*5,4,1,2,3))
c;crit_value=mean(c);crit_value







c<-replicate(5,fun4(10,15,8,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,35,35,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,10,10,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(15,15,15,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,25,30,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(40,40,40,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(40,25,20,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,15,10,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,8,9,0,0,0,4,1,1,1))
c;crit_value=mean(c);crit_value



c<-replicate(5,fun4(10,15,8,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,35,35,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,10,10,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(15,15,15,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,25,30,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(40,40,40,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(40,25,20,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,15,10,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,8,9,0,0,0,4,4,4,4))
c;crit_value=mean(c);crit_value



c<-replicate(5,fun4(10,15,8,0,0,0,4,4,2,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,20,0,0,0,4,1,3,2))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,35,35,0,0,0,4,8,1,2))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,10,10,0,0,0,4,2,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(15,15,15,0,0,0,4,10,6,15))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(20,25,30,0,0,0,4,2,4,4))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(40,40,40,0,0,0,4,3,2,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(40,25,20,0,0,0,4,1,4,5))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,15,10,0,0,0,4,1,8,9))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,8,9,0,0,0,4,1,2,3))
c;crit_value=mean(c);crit_value



c<-replicate(5,fun4(10,15,8,0,0,0,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1,1.5,2,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*1.2,1.5*1.2,2*1.2,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*1.4,1.5*1.4,2*1.4,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*1.6,1.5*1.6,2*1.6,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*1.8,1.5*1.8,2*1.8,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*2,1.5*2,2*2,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*2.2,1.5*2.2,2*2.2,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*2.4,1.5*2.4,2*2.4,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*2.6,1.5*2.6,2*2.6,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*2.8,1.5*2.8,2*2.8,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*3,1.5*3,2*3,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*3.2,1.5*3.2,2*3.2,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*3.4,1.5*3.4,2*3.4,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*3.6,1.5*3.6,2*3.6,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*3.8,1.5*3.8,2*3.8,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*4,1.5*4,2*4,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*4.2,1.5*4.2,2*4.2,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*4.4,1.5*4.4,2*4.4,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*4.6,1.5*4.6,2*4.6,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*4.8,1.5*4.8,2*4.8,4,8,4,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,15,8,1*5,1.5*5,2*5,4,8,4,1))
c;crit_value=mean(c);crit_value




c<-replicate(5,fun4(7,8,9,0,0,0,2,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,8,9,0,0,0,2,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,8,9,0,0,0,2,4,8,7))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,8,9,0,0,0,2,1,9,9))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,8,9,0,0,0,2,9,9,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,8,9,0,0,0,2,1,7,6))
c;crit_value=mean(c);crit_value

c<-replicate(5,fun4(30,15,10,0,0,0,2,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,15,10,0,0,0,2,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,15,10,0,0,0,2,4,8,7))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,15,10,0,0,0,2,1,9,9))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,15,10,0,0,0,2,9,9,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,15,10,0,0,0,2,1,7,6))
c;crit_value=mean(c);crit_value

c<-replicate(5,fun4(10,14,15,0,0,0,2,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,14,15,0,0,0,2,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,14,15,0,0,0,2,4,8,7))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,14,15,0,0,0,2,1,9,9))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,14,15,0,0,0,2,9,9,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(10,14,15,0,0,0,2,1,7,6))
c;crit_value=mean(c);crit_value

c<-replicate(5,fun4(30,8,8,0,0,0,2,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,8,8,0,0,0,2,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,8,8,0,0,0,2,4,8,7))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,8,8,0,0,0,2,1,9,9))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,8,8,0,0,0,2,9,9,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(30,8,8,0,0,0,2,1,7,6))
c;crit_value=mean(c);crit_value

c<-replicate(5,fun4(8,7,30,0,0,0,2,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(8,7,30,0,0,0,2,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(8,7,30,0,0,0,2,4,8,7))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(8,7,30,0,0,0,2,1,9,9))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(8,7,30,0,0,0,2,9,9,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(8,7,30,0,0,0,2,1,7,6))
c;crit_value=mean(c);crit_value

c<-replicate(5,fun4(7,10,25,0,0,0,2,1,1,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,10,25,0,0,0,2,1,2,3))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,10,25,0,0,0,2,4,8,7))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,10,25,0,0,0,2,1,9,9))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,10,25,0,0,0,2,9,9,1))
c;crit_value=mean(c);crit_value
c<-replicate(5,fun4(7,10,25,0,0,0,2,1,7,6))
c;crit_value=mean(c);crit_value


c<-replicate(4,fun4(80,80,80,0,0,0,2,1,1,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(80,80,80,0,0,0,2,1,2,3));c
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(80,80,80,0,0,0,2,4,8,8))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(80,80,80,0,0,0,2,1,9,9))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(80,80,80,0,0,0,2,9,9,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(80,80,80,0,0,0,2,1,7,6))
crit_value=sum(c)/4;crit_value


c<-replicate(4,fun4(8,10,15,0,0,0,2,1,1,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(8,10,15,0,0,0,2,1,2,3));c
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(8,10,15,0,0,0,2,4,8,8))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(8,10,15,0,0,0,2,1,9,9))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(8,10,15,0,0,0,2,9,9,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(8,10,15,0,0,0,2,1,7,6))
crit_value=sum(c)/4;crit_value


c<-replicate(4,fun4(30,40,45,0,0,0,2,1,1,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,40,45,0,0,0,2,1,2,3));c
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,40,45,0,0,0,2,4,8,8))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,40,45,0,0,0,2,1,9,9))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,40,45,0,0,0,2,9,9,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,40,45,0,0,0,2,1,7,6))
crit_value=sum(c)/4;crit_value

c<-replicate(4,fun4(20,20,25,0,0,0,2,1,1,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(20,20,25,0,0,0,2,1,2,3));c
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(20,20,25,0,0,0,2,4,8,8))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(20,20,25,0,0,0,2,1,9,9))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(20,20,25,0,0,0,2,9,9,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(20,20,25,0,0,0,2,1,7,6))
crit_value=sum(c)/4;crit_value

c<-replicate(4,fun4(30,25,25,0,0,0,2,1,1,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,25,25,0,0,0,2,1,2,3));c
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,25,25,0,0,0,2,4,8,8))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,25,25,0,0,0,2,1,9,9))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,25,25,0,0,0,2,9,9,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,25,25,0,0,0,2,1,7,6))
crit_value=sum(c)/4;crit_value

c<-replicate(4,fun4(8,10,15,0,0,0,2,1,1,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(8,10,15,0,0,0,2,1,2,3))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(8,10,15,0,0,0,2,4,3,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(8,10,15,0,0,0,2,0.1,0.2,0.3))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(8,10,15,0,0,0,2,0.1,1,0.5))
crit_value=sum(c)/4;crit_value

c<-replicate(4,fun4(30,40,45,0,0,0,2,1,1,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,40,45,0,0,0,2,1,2,3))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,40,45,0,0,0,2,4,3,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,40,45,0,0,0,2,0.1,0.2,0.3))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(30,40,45,0,0,0,2,0.1,1,0.5))
crit_value=sum(c)/4;crit_value

c<-replicate(4,fun2(20,20,25,0,0,0,2,1,1,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun2(20,20,25,0,0,0,2,1,2,3))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun2(20,20,25,0,0,0,2,4,3,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun2(20,20,25,0,0,0,2,0.1,0.2,0.3))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun2(20,20,25,0,0,0,2,0.1,1,0.5))
crit_value=sum(c)/4;crit_value

c<-replicate(4,fun2(30,25,25,0,0,0,2,1,1,1))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun2(30,25,25,0,0,0,2,1,2,3))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun2(30,25,25,0,0,0,2,4,3,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun2(30,25,25,0,0,0,2,0.1,0.2,0.3))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun2(30,25,25,0,0,0,2,0.1,1,0.5))
crit_value=sum(c)/4;crit_value



c<-replicate(4,fun4(10,12,14,0,0,0,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3,3.4,3.6,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*1.2,3.4*1.2,3.6*1.2,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*1.4,3.4*1.4,3.6*1.4,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*1.6,3.4*1.6,3.6*1.6,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*1.8,3.4*1.8,3.6*1.8,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*2,3.4*2,3.6*2,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*2.2,3.4*2.2,3.6*2.2,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*2.4,3.4*2.4,3.6*2.4,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*2.6,3.4*2.6,3.6*2.6,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*2.8,3.4*2.8,3.6*2.8,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*3,3.4*3,3.6*3,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*3.2,3.4*3.2,3.6*3.2,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*3.4,3.4*3.4,3.6*3.4,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*3.6,3.4*3.6,3.6*3.6,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*3.8,3.4*3.8,3.6*3.8,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*4,3.4*4,3.6*4,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*4.2,3.4*4.2,3.6*4.2,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*4.4,3.4*4.4,3.6*4.4,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*4.6,3.4*4.6,3.6*4.6,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*4.8,3.4*4.8,3.6*4.8,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*5,3.4*5,3.6*5,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*5.2,3.4*5.2,3.6*5.2,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*5.4,3.4*5.4,3.6*5.4,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*5.6,3.4*5.6,3.6*5.6,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*5.8,3.4*5.8,3.6*5.8,6,3,2,2))
crit_value=sum(c)/4;crit_value
c<-replicate(4,fun4(10,12,14,3*6,3.4*6,3.6*6,6,3,2,2))
crit_value=sum(c)/4;crit_value