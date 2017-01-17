# Step1: Get the data
# Download the Leukaemia data: (Golub et al., 1999) 
# Save the training and test data as well as the labels into one dataset file
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
library("golubEsets")
data(Golub_Train)
golub_train <- data.frame(Golub_Train)[1:7129]
data(Golub_Test)
golub_test <- data.frame(Golub_Test)[1:7129]
golub_train.pheno=pData(Golub_Train)
attach(golub_train.pheno)
golub_train.type = as.numeric(ALL.AML=="AML")
table(golub_train.type)
detach(golub_train.pheno)
golub_test.pheno=pData(Golub_Test)
attach(golub_test.pheno)
golub_test.type = as.numeric(ALL.AML=="AML")
table(golub_test.type)
detach(golub_test.pheno)
save(golub_train,golub_test,golub_train.type,golub_test.type, file = "saveddataset.RData")


#===================================================================
# Step2: Model 1:Use the basic LDA
#===================================================================
p=7129
n=38
library(matrixcalc)
# construct a new function to get miu, delta and Sigma from a dataset X
# and its corresponding type vector y
Get_MDS<-function(X,y,n){
  golub_train1=X[y==0,]
  golub_train2=X[y==1,]
  mean1=colMeans(golub_train1)
  mean2=colMeans(golub_train2)
  # The estimator for miu 
  miu0=(mean1+mean2)/2 
  # The estimator for delta 
  delta0=mean1-mean2
  # The sample covariance matrix 
  X1=golub_train1-mean1%*%t(c(rep(1,7129)))
  X2=golub_train2-mean2%*%t(c(rep(1,7129)))
  # n is the sample size
  Sigma=(t(X1)%*%as.matrix(X1)+t(X2)%*%as.matrix(X2))/n
  list0=list("Sigma"=Sigma,"delta"=delta0,"miu"=miu0)
  return(list0)
}
a=Get_MDS(golub_train,golub_train.type,38)
Sigma=a$Sigma
delta0=a$delta
miu0=a$miu
#===================================================================
# Now we should find the inverse of Sigma
# We use the function eigs_sym() in package rARPACK to find the first 50 eigen
#===================================================================
library(rARPACK)
# Here we construct a function Solve_sigma() to find the estimator of Sigma^(-1)
Solve_sigma<-function(Sigma,tol){
  svd0<-eigs_sym(Sigma,150,which="LA")
  a=svd0$values>tol
  if (sum(as.numeric(a))==150)
  {print("150 is not enough")}
  U=svd0$vectors[,a]
  DIAG=diag(svd0$values[a]^(-1))
  Sigma_plus_inverse=U%*%DIAG%*%t(U)
  return(Sigma_plus_inverse)
}
Sigma_inverse=Solve_sigma(Sigma,tol)
#===================================================================
#Now we have Sigma_inverse, delta0 and miu0, and then we have the classification rule
# we use the classification rule to find the misclassfication rate
#===================================================================
Rate_mis<-function(Sigma_inverse,delta,miu,type){
  type_predict=matrix(0,1,34)
  for (i in 1:34){
    X=golub_test[i,]
    if (t(delta0)%*%Sigma_inverse%*%(t(X)-miu0)>=0) 
      type_predict[i]=0
    else type_predict[i]=1
  }
  rate=sum(as.numeric(type_predict!=type))/34
  return(rate)
}
# Apply the above function to find the misclassification rate
rate_LDA0=Rate_mis(Sigma_inverse,delta0,miu0,golub_test.type)
rate_LDA0 # 0.1470588


#================================================================
# Each following method applies the same thresholding to the estimate of delta
#===================================================================
alpha=0.3
M2=300
delta_thr<-function(delta0,alpha,M2){
  delta1=delta0*as.numeric(abs(delta0)>M2*(log(p)/n)^alpha)
  return(delta1)
}

#=====================================================================
# Step 3: Model 2: a slightly different LDA, only threshold delta
#================================================================
Rate_LDA1=matrix(0,1,10)
for (i in 1:10){
  M2=50*i
  delta1<-delta_thr(a$delta,alpha,M2)
  Rate_LDA1[i]=Rate_mis(Sigma_inverse,delta1,a$miu,golub_test.type)
}
Rate_LDA1 #all are 0.1470588,keep the same

#======================================================================
# Step 4: Model 3: banding estimation
#================================================================
Sigma_B<-function(Sigma,M1){
  SigmaB=Sigma
  SigmaB[abs(SigmaB)<M1]=0
  return(SigmaB)
}
Rate_LDA2=matrix(0,4,10)
M2=300
for(j in 1:4){
  print("j is",j)
  tol=j*5*10^7
  for (i in 1:10){
    print(i)
    M1=i*10^6
    SigmaB=Sigma_B(a$Sigma,M1)
    SigmaB_inverse=Solve_sigma(SigmaB,tol)
    delta1<-delta_thr(a$delta,alpha,M2)
    Rate_LDA2[j,i]=Rate_mis(SigmaB_inverse,delta1,a$miu,golub_test.type)
  }
}
attach(mtcars)
par(mfrow=c(2,2))
M1=10^6*c(1:10)
plot(M1,Rate_LDA2[1,],type="l",xlab="M1",ylab="misclassification rate", main="tol=0.5*10^8")
points(M1,Rate_LDA2[1,],cex=1,col="dark red")
plot(M1,Rate_LDA2[2,],type="l",xlab="M1",ylab="misclassification rate", main="tol=10^8")
points(M1,Rate_LDA2[2,],cex=1,col="dark red")
plot(M1,Rate_LDA2[3,],type="l",xlab="M1",ylab="misclassification rate", main="tol=1.5*10^8")
points(M1,Rate_LDA2[3,],cex=1,col="dark red")
plot(M1,Rate_LDA2[4,],type="l",xlab="M1",ylab="misclassification rate", main="tol=2*10^8")
points(M1,Rate_LDA2[4,],cex=1,col="dark red")
# Check the image of banding matrix B when tol=10^8 and M1=9*10^6
library(SparseM)
# tol=10^8 
M1=9*10^6
B=matrix(1,p,p)
index=abs(a$Sigma)<M1
B[index]=0
B_sparse=Matrix(B,sparse=TRUE)
image(B_sparse)
image1=image(B_sparse[1:100,1:100])
image2=image(B_sparse[6000:6100,6000:6100])
image1
image2
#================================================================
# Step 5: Model 4: tapering estimation
#================================================================
Sigma_T<-function(Sigma,M1){
  T=matrix(1,p,p)
  T[abs(Sigma)<M1/2]=0
  index=abs(Sigma)<M1&abs(Sigma)>=M1/2 
  T[index]=(2-2*abs(Sigma)/M1)[index]
  SigmaT=Sigma*T
  return(SigmaT)
}
Rate_LDA3=matrix(0,4,10)
for (j in 1:6){
  tol=j*0.5*10^8
  for (i in 1:10){
    print(i)
    M1=i*10^6
    SigmaT=Sigma_T(a$Sigma,M1)
    SigmaT_inverse=Solve_sigma(SigmaT,tol)
    delta1<-delta_thr(a$delta,alpha,M2)
    Rate_LDA3[j,i]=Rate_mis(SigmaT_inverse,delta1,a$miu,golub_test.type)
  }
}
attach(mtcars)
par(mfrow=c(3,2))
M1=10^6*c(1:10)
j=1
plot(M1,Rate_LDA3[j,],type="l",xlab="M1",ylab="misclassification rate", main="tol=0.5*10^8")
points(M1,Rate_LDA3[j,],cex=1,col="dark red")
j=2
plot(M1,Rate_LDA3[j,],type="l",xlab="M1",ylab="misclassification rate", main="tol=1*10^8")
points(M1,Rate_LDA3[j,],cex=1,col="dark red")
j=3
plot(M1,Rate_LDA3[j,],type="l",xlab="M1",ylab="misclassification rate", main="tol=1.5*10^8")
points(M1,Rate_LDA3[j,],cex=1,col="dark red")
j=4
plot(M1,Rate_LDA3[j,],type="l",xlab="M1",ylab="misclassification rate", main="tol=2*10^8")
points(M1,Rate_LDA3[j,],cex=1,col="dark red")
j=5
plot(M1,Rate_LDA3[j,],type="l",xlab="M1",ylab="misclassification rate", main="tol=2.5*10^8")
points(M1,Rate_LDA3[j,],cex=1,col="dark red")
j=6
plot(M1,Rate_LDA3[j,],type="l",xlab="M1",ylab="misclassification rate", main="tol=3*10^8")
points(M1,Rate_LDA3[j,],cex=1,col="dark red")
# Check the image of tapering matrix T when tol=2*10^8, M1=9*10^6 
# tol=2*10^8
M1=9*10^6 
T=matrix(1,p,p)
T[abs(Sigma)<M1/2]=0
index=abs(Sigma)<M1&abs(Sigma)>=M1/2 
T[index]=(2-2*abs(Sigma)/M1)[index]
T_sparse=Matrix(T,sparse=TRUE)
image(T_sparse)
image3=image(T_sparse[1:100,1:100])
image4=image(T_sparse[6000:6100,6000:6100])
image3
image4
#================================================================
# Step 6: Model5: Block Thresholding
#================================================================
# We should notice that the assuption is bandable matrice
# We construct the following function to find the best order 
# This is to make the matrix close to the "bandable" definition
p=7129
Index_order<-function(Sigma){
  index_order=matrix(0,1,p)
  Sigma0=abs(Sigma)
  diag(Sigma0)=0
  index_order[1]=which(Sigma0==max(Sigma0),arr.ind=TRUE)[1]
  index_order[2]=which(Sigma0==max(Sigma0),arr.ind=TRUE)[2]
  l=index_order[2]
  for (i in 2:p-1){
    index_a=index_order[1:i]
    vector=Sigma0[l,]
    vector[index_a]=0
    l=which.max(vector)
    index_order[i+1]=l
  }
  return(index_order)
}

# We construct a new Solve_inverse1() wrt the number of eigenvalues n_eig considered
# here, we don't consider tol any more
# This is for the case when Solve_inverse() is used, "not enough" is always poped up.
Solve_sigma1<-function(Sigma,n_eig){
  svd0<-eigs_sym(Sigma,n_eig,which="LA")
  U=svd0$vectors
  DIAG=diag(svd0$values^(-1))
  Sigma_plus_inverse=U%*%DIAG%*%t(U)
  return(Sigma_plus_inverse)
}

# adaptive estimator
# lambda0 is the parameter that can be changed
# We construct a function to get the block thresholding covariance estimator
Sigma_BT<-function(Sigma){
  # construct the block thresholding estimator
  # k0=log(p)=8
  # n/log(n)=10.45
  # So in this case, need to keep the diagonal blocks 
  Sigma_BT=matrix(0,p,p)
  # we find that 7129=891*8+1, so the last block on the diagonal is of size 1
  # spectral norms of the diagonal blocks are stored in a vector
  norm_diag=matrix(0,892)
  for (i in 1:891){
    a=8*(i-1)
    Sigma_BT[(a+1):(a+8),(a+1):(a+8)]=Sigma[(a+1):(a+8),(a+1):(a+8)]
    norm_diag[i]=spectral.norm(Sigma[(a+1):(a+8),(a+1):(a+8)])
  }
  Sigma_BT[7129,7129]=Sigma[7129,7129]
  norm_diag[892]=spectral.norm(Sigma[7129,7129])
  lambda0=6
  # threshold blocks that d(B)=k0 except the diagonal ones
  # Also need to threshold one block lies in the last row and the last column 
  add_matrix=matrix(0,p,p)
  # 1+445*2=891
  for (j in 1:445){
    b=16*(j-1)
    a1=(b+1):(b+8)
    a2=(b+9):(b+16)
    a3=(b+17):(b+24)
    b1=sqrt((8+log(p))/n)
    add_matrix[a1,a2]=Sigma[a1,a2]*(spectral.norm(Sigma[a1,a2])>lambda0*b1*sqrt(norm_diag[2*j-1]*norm_diag[2*j]))
    add_matrix[a1,a3]=Sigma[a1,a3]*(spectral.norm(Sigma[a1,a3])>lambda0*b1*sqrt(norm_diag[2*j-1]*norm_diag[2*j+1]))
    add_matrix[a2,a3]=Sigma[a2,a3]*(spectral.norm(Sigma[a2,a3])>lambda0*b1*sqrt(norm_diag[2*j]*norm_diag[2*j+1]))
  }
  add_matrix[7121:7128,7129]=Sigma[7121:7128,7129]*(spectral.norm(Sigma[7121:7128,7129])>lambda0*b1*sqrt(norm_diag[891]*norm_diag[892]))
# kill other blocks, no need to write code, as Sigma_BT is initiated with all 0's.
# The final covariance estimator is as follows:
  Sigma_BT=Sigma_BT+add_matrix+t(add_matrix)
  return(Sigma_BT)
}

n=38
a=Get_MDS(golub_train,golub_train.type,n)
IO=Index_order(a$Sigma)
SigmaBT=Sigma_BT(a$Sigma[IO,IO])
alpha=0.3
M2=300
Rate_LDA4=matrix(0,1,25)
for ( i in c(21:25)){
  print(i)
  n_eig=10*i+40
  SigmaBT_inverse=Solve_sigma1(SigmaBT,n_eig)
  delta1<-delta_thr(a$delta,alpha,M2)
  Rate_LDA4[i]=Rate_mis(SigmaBT_inverse,delta1[IO],a$miu[IO],golub_test.type)
}

plot.new()
par(mfrow=c(1,1))
n_eig=seq(50,290,by=10)
plot(n_eig,Rate_LDA3,type="l",xlab="the number of eigenvalues n_eig", ylab="misclassification rate")
points(n_eig,Rate_LDA3,cex=1,col="dark red")