rm(list=ls(all=TRUE))
library(FAdist)
library(rootSolve)
library(pso)
nllike_Gamma3=function(dat,x){
  n0=length(dat)
  n0*x[1]*log(x[2])+n0*log(gamma(x[1]))-(x[1]-1)*sum(log(dat-x[3]))+sum(dat-x[3])/x[2]
}
#BL
BL=function(dat,lower,upper,w=c(1.1,0.1),c0=c(1.49,1.49),s0=10){
  n0=length(dat)
  nllike_G3=function(x){
    n0*x[1]*log(x[2])+n0*log(gamma(x[1]))-(x[1]-1)*sum(log(dat-x[3]))+sum(dat-x[3])/x[2]-log(dat[1]-x[3])+log(dat[2]-x[3])}
  BLE=psoptim(par=rep(NA,3),fn=nllike_G3,lower=lower,upper=upper,control=list(w=w,c.p=c0[1],c.g=c0[2],s=s0))$par
  return(BLE)
}
#MLE
ML=function(dat,lower,upper,w=c(1.1,0.1),c0=c(1.49,1.49),s0=10){
  n0=length(dat)
  nllike_G3=function(x){
    n0*x[1]*log(x[2])+n0*log(gamma(x[1]))-(x[1]-1)*sum(log(dat-x[3]))+sum(dat-x[3])/x[2]}
  MLE=psoptim(par=rep(NA,3),fn=nllike_G3,lower=lower,upper=upper,control=list(w=w,c.p=c0[1],c.g=c0[2],s=s0))$par
  return(MLE)
}
#MPS
MPS=function(dat,lower,upper,w=c(1.1,0.1),c0=c(1.49,1.49),s0=10){
  n0=length(dat)+1
  MPSL=function(x){
    y=0
    for (i in 1:(n0-1)) y[i]=pgamma3(dat[i],x[1],x[2],x[3])
    D=c(y[1],diff(y),1-y[n0-1])
    -sum(log(D))/n0
  }
  MPL=psoptim(par=rep(NA,3),fn=MPSL,lower=lower,upper=upper,control=list(w=w,c.p=c0[1],c.g=c0[2],s=s0))$par
  return(MPL)
}
#thres-est
bmle_f <- function(dat){
  a=2;b=-1;n0=length(dat);k=round(1.6*sqrt(n0))
  min_nums <- head(dat,k+1)
  ftheta <- function(theta){
    term1 <- 0
    for (i in 1:k){term1=term1+(min_nums[k+1]-theta)/(min_nums[i]-theta)-1}
    term2 <- 0
    for (i in 1:k){term2=term2+log((min_nums[k+1]-theta)/(min_nums[i]-theta))}
    term3 <- 0
    for (i in 2:k){term3=term3+(min_nums[k+1]-theta)/(min_nums[i]-theta)}
    d <- (2-a)*(min_nums[k+1]-theta)/(min_nums[1]-theta)
    return(term1-term2*(a+b+1+d+term3)/k)
  }
  min(uniroot.all(ftheta,c(-1e5,min_nums[1]-1e-6),tol = 1e-7))
}
#MLE for gamma2(two step)
twostepgamma=function(dat){
  n0=length(dat)
  par3=bmle_f(dat)
  if (par3==Inf) {par3=dat[1]-1e-5}
  nllike_G2=function(x){
    n0*x[1]*log(x[2])+n0*log(gamma(x[1]))-(x[1]-1)*sum(log(dat-par3))+sum(dat-par3)/x[2]}
  par12=optim(c(1,1),nllike_G2)$par
  return(c(par12,par3))
}
n=20
dat=sort(rgamma3(n,4,1,1));
par0=twostepgamma(dat);bootn=500
partwostep=matrix(0,nrow=3,ncol=bootn)
for (B in 1:bootn) {
  datboot=sort(rgamma3(n,shape=par0[1],scale=par0[2],thres=par0[3]))
  partwostep[,B]=twostepgamma(datboot)}

shaperange=quantile(partwostep[1,],c(0.025,0.975))
scalerange=quantile(partwostep[2,],c(0.025,0.975))
locationrange=quantile(partwostep[3,],c(0.025,0.095))
lower1=c(shaperange[1],scalerange[1],locationrange[1]);lower2=c(0.001,0.001,-5)
upper1=c(shaperange[2],scalerange[2],locationrange[2]);upper2=c(20,10,dat[1]-1e-10)
(BLE1=BL(dat,lower=lower1,upper=upper1))
(BLE2=BL(dat,lower=lower2,upper=upper2))
(MPSE1=MPS(dat,lower=lower1,upper=upper1))
(MPSE2=MPS(dat,lower=lower2,upper=upper2))
(MLE1=ML(dat,lower=lower1,upper=upper1))
(MLE2=ML(dat,lower=lower2,upper=upper2))
################################################################################
###################simulation######start#########simulation#####################
################################################################################
start_time=Sys.time()
ite=1000;n=20
shapev=c(0.8,1,1.5,2,2.5,3,3.5,4);lenshape=length(shapev)
reliability0=c(0.99,0.95,0.9,0.5);lenrelia=length(reliability0)
results_list=list()
################################################################################
library(parallel);library(doParallel);library(foreach)
################################################################################
num_cores <- detectCores() - 1  
cl <- makeCluster(num_cores)
registerDoParallel(cl)  
clusterExport(cl, c("bmle_f","twostepgamma","BL","ML","MPS","nllike_Gamma3","ite","n","shapev","lenshape","reliability0","lenrelia"))
results_list <- foreach(shapei=shapev,.packages = c("FAdist","rootSolve","pso")) %do% {
  #pars
  parmle=matrix(0,3,ite);parmle2=matrix(0,3,ite)
  parbl=matrix(0,3,ite);parbl2=matrix(0,3,ite)
  parmps=matrix(0,3,ite);parmps2=matrix(0,3,ite)
  parple=matrix(0,3,ite);
  #qls
  qlmle=matrix(0,lenrelia,ite);qlmle2=matrix(0,lenrelia,ite)
  qlbl=matrix(0,lenrelia,ite);qlbl2=matrix(0,lenrelia,ite)
  qlmps=matrix(0,lenrelia,ite);qlmps2=matrix(0,lenrelia,ite)
  qlple=matrix(0,lenrelia,ite);
  #relias
  reliamle=matrix(0,lenrelia,ite);reliamle2=matrix(0,lenrelia,ite)
  reliabl=matrix(0,lenrelia,ite);reliabl2=matrix(0,lenrelia,ite)
  reliamps=matrix(0,lenrelia,ite);reliamps2=matrix(0,lenrelia,ite)
  reliaple=matrix(0,lenrelia,ite);
  #meanlife
  lifemmle=rep(0,ite);lifemmle2=rep(0,ite)
  lifembl=rep(0,ite);lifembl2=rep(0,ite)
  lifemmps=rep(0,ite);lifemmps2=rep(0,ite)
  lifemple=rep(0,ite)
  ##############################################################################
  trues=c(shapei,1,1)
  set.seed(20250712)
  datm=matrix(rgamma3(n*ite,trues[1],trues[2],trues[3]),n)
  bootn=500;QL0=qgamma3(1-reliability0,shape=trues[1],scale=trues[2],thres=trues[3])
  shaperange=matrix(0,nrow=ite,ncol=2);scalerange=matrix(0,nrow=ite,ncol=2);locationrange=matrix(0,nrow=ite,ncol=2)
  for (i in 1:ite){
    print(c(shapei,i))
    dat=sort(datm[,i])
    par0=twostepgamma(dat);parple[,i]=par0
    
    for (relia in 1:lenrelia){
      qlple[relia,i]=qgamma3(1-reliability0[relia],par0[1],par0[2],par0[3])}
    reliaple[,i]=1-pgamma3(QL0,par0[1],par0[2],par0[3])
    lifemple[i]=par0[3]+par0[2]*par0[1]
    
    partwostep=matrix(0,nrow=3,ncol=bootn)
    for (B in 1:bootn) {
      datboot=sort(rgamma3(n,shape=par0[1],scale=par0[2],thres=par0[3]))
      partwostep[,B]=twostepgamma(datboot)}
    
    shaperange[i,]=quantile(partwostep[1,],c(0.025,0.975))
    scalerange[i,]=quantile(partwostep[2,],c(0.025,0.975))
    locationrange[i,]=c(quantile(partwostep[3,],0.025),dat[1]-1e-5)
    lower1=c(shaperange[i,][1],scalerange[i,][1],locationrange[i,][1]);lower2=c(0.001,0.001,-5)
    upper1=c(shaperange[i,][2],scalerange[i,][2],locationrange[i,][2]);upper2=c(20,10,dat[1]-1e-5)
    
    parmle[,i]=tryCatch({ML(dat,lower=lower1,upper=upper1)},error=function(e) {rep(1e-5,3)})
    for (relia in 1:lenrelia){
      qlmle[relia,i]=qgamma3(1-reliability0[relia],parmle[,i][1],parmle[,i][2],parmle[,i][3])}
    reliamle[,i]=1-pgamma3(QL0,parmle[,i][1],parmle[,i][2],parmle[,i][3])
    lifemmle[i]=parmle[,i][3]+parmle[,i][2]*parmle[,i][1]
    
    parmle2[,i]=tryCatch({ML(dat,lower=lower2,upper=upper2)},error=function(e) {rep(1e-5,3)})
    for (relia in 1:lenrelia){
      qlmle2[relia,i]=qgamma3(1-reliability0[relia],parmle2[,i][1],parmle2[,i][2],parmle2[,i][3])}
    reliamle2[,i]=1-pgamma3(QL0,parmle2[,i][1],parmle2[,i][2],parmle2[,i][3])
    lifemmle2[i]=parmle2[,i][3]+parmle2[,i][2]*parmle2[,i][1]
    
    parbl[,i]=tryCatch({BL(dat,lower=lower1,upper=upper1)},error=function(e) {rep(1e-5,3)})
    for (relia in 1:lenrelia){
      qlbl[relia,i]=qgamma3(1-reliability0[relia],parbl[,i][1],parbl[,i][2],parbl[,i][3])}
    reliabl[,i]=1-pgamma3(QL0,parbl[,i][1],parbl[,i][2],parbl[,i][3])
    lifembl[i]=parbl[,i][3]+parbl[,i][2]*parbl[,i][1]
    
    parbl2[,i]=tryCatch({BL(dat,lower=lower2,upper=upper2)},error=function(e) {rep(1e-5,3)})
    for (relia in 1:lenrelia){
      qlbl2[relia,i]=qgamma3(1-reliability0[relia],parbl2[,i][1],parbl2[,i][2],parbl2[,i][3])}
    reliabl2[,i]=1-pgamma3(QL0,parbl2[,i][1],parbl2[,i][2],parbl2[,i][3])
    lifembl2[i]=parbl2[,i][3]+parbl2[,i][2]*parbl2[,i][1]
    
    parmps[,i]=tryCatch({MPS(dat,lower=lower1,upper=upper1)},error=function(e) {rep(1e-5,3)})
    for (relia in 1:lenrelia){
      qlmps[relia,i]=qgamma3(1-reliability0[relia],parmps[,i][1],parmps[,i][2],parmps[,i][3])}
    reliamps[,i]=1-pgamma3(QL0,parmps[,i][1],parmps[,i][2],parmps[,i][3])
    lifemmps[i]=parmps[,i][3]+parmps[,i][2]*parmps[,i][1]
    
    parmps2[,i]=tryCatch({MPS(dat,lower=lower2,upper=upper2)},error=function(e) {rep(1e-5,3)})
    for (relia in 1:lenrelia){
      qlmps2[relia,i]=qgamma3(1-reliability0[relia],parmps2[,i][1],parmps2[,i][2],parmps2[,i][3])}
    reliamps2[,i]=1-pgamma3(QL0,parmps2[,i][1],parmps2[,i][2],parmps2[,i][3])
    lifemmps2[i]=parmps2[,i][3]+parmps2[,i][2]*parmps2[,i][1]
  }
  list(parmle=parmle,parmle2=parmle2,parbl=parbl,parbl2=parbl2,parmps=parmps,parmps2=parmps2,parple=parple,
       qlmle=qlmle,qlmle2=qlmle2,qlbl=qlbl,qlbl2=qlbl2,qlmps=qlmps,qlmps2=qlmps2,qlple=qlple,
       reliamle=reliamle,reliamle2=reliamle2,reliabl=reliabl,reliabl2=reliabl2,reliamps=reliamps,reliamps2=reliamps2,reliaple=reliaple,
       lifemmle=lifemmle,lifemmle2=lifemmle2,lifembl=lifembl,lifembl2=lifembl2,lifemmps=lifemmps,lifemmps2=lifemmps2,lifemple=lifemple,
       shaperange=shaperange,scalerange=scalerange,locationrange=locationrange)
}
stopCluster(cl)
end_time=Sys.time()
################################################################################
################################################################################
################################################################################
result_df_list=list()
for (shapei in 1:lenshape){
  try=results_list[[shapei]];result_df=data.frame()
  parstrue=c(shapev[shapei],1,1)
  QL0=qgamma3(1-reliability0,shape=parstrue[1],scale=parstrue[2],thres=parstrue[3])
  lifem0=parstrue[3]+parstrue[2]*parstrue[1]
  
  nonconverge_index=which(try$parmle[1,]==1e-5)
  if (length(nonconverge_index) == 0) {
    bias=apply(try$parmle,1,mean)-parstrue
    rmse=sqrt((apply(try$parmle,1,mean)-parstrue)^2+(apply(try$parmle,1,sd))^2)
    QLbias=apply(try$qlmle,1,mean)-QL0
    QLrmse=sqrt((apply(try$qlmle,1,mean)-QL0)^2+(apply(try$qlmle,1,sd))^2)
    reliabias=apply(try$reliamle,1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliamle,1,mean)-QL0)^2+(apply(try$reliamle,1,sd))^2)
    lifembias=mean(try$lifemmle)-lifem0
    lifemrmse=sqrt((mean(try$lifemmle)-lifem0)^2+sd(try$lifemmle)^2)
  } else {
    bias=apply(try$parmle[,-nonconverge_index],1,mean)-parstrue
    rmse=sqrt((apply(try$parmle[,-nonconverge_index],1,mean)-parstrue)^2+(apply(try$parmle[,-nonconverge_index],1,sd))^2)
    QLbias=apply(try$qlmle[,-nonconverge_index],1,mean)-QL0
    QLrmse=sqrt((apply(try$qlmle[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$qlmle[,-nonconverge_index],1,sd))^2)
    reliabias=apply(try$reliamle[,-nonconverge_index],1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliamle[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$reliamle[,-nonconverge_index],1,sd))^2)
    lifembias=mean(try$lifemmle[-nonconverge_index])-lifem0
    lifemrmse=sqrt((mean(try$lifemmle[-nonconverge_index])-lifem0)^2+sd(try$lifemmle[-nonconverge_index])^2)
  }
  result_df=rbind(result_df,data.frame(method="MLE",shape=shapev[shapei],biassum=sum(abs(bias)),rmsesum=sum(rmse),
                                       QLbiassum=sum(abs(QLbias)),QLrmsesum=sum(QLrmse),shapebias=bias[1],scalebias=bias[2],
                                       thresholdbias=bias[3],shapermse=rmse[1],scalermse=rmse[2],thresholdrmse=rmse[3],
                                       noncon_prop=length(nonconverge_index)/ite,
                                       QL0.01bias=QLbias[1],QL0.05bias=QLbias[2],QL0.1bias=QLbias[3],QL0.5bias=QLbias[4],
                                       QL0.01rmse=QLrmse[1],QL0.05rmse=QLrmse[2],QL0.1rmse=QLrmse[3],QL0.5rmse=QLrmse[4],
                                       RL0.99bias=reliabias[1],RL0.95bias=reliabias[2],RL0.9bias=reliabias[3],RL0.5bias=reliabias[4],
                                       RL0.99rmse=reliarmse[1],RL0.95rmse=reliarmse[2],RL0.9rmse=reliarmse[3],RL0.5rmse=reliarmse[4],
                                       lifembias=lifembias,lifemrmse=lifemrmse))
  
  nonconverge_index=which(try$parmle2[1,]==1e-5)
  if (length(nonconverge_index) == 0) {
    bias=apply(try$parmle2,1,mean)-parstrue
    rmse=sqrt((apply(try$parmle2,1,mean)-parstrue)^2+(apply(try$parmle2,1,sd))^2)
    QLbias=apply(try$qlmle2,1,mean)-QL0
    QLrmse=sqrt((apply(try$qlmle2,1,mean)-QL0)^2+(apply(try$qlmle2,1,sd))^2)
    reliabias=apply(try$reliamle2,1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliamle2,1,mean)-QL0)^2+(apply(try$reliamle2,1,sd))^2)
    lifembias=mean(try$lifemmle2)-lifem0
    lifemrmse=sqrt((mean(try$lifemmle2)-lifem0)^2+sd(try$lifemmle2)^2)
  } else {
    bias=apply(try$parmle2[,-nonconverge_index],1,mean)-parstrue
    rmse=sqrt((apply(try$parmle2[,-nonconverge_index],1,mean)-parstrue)^2+(apply(try$parmle2[,-nonconverge_index],1,sd))^2)
    QLbias=apply(try$qlmle2[,-nonconverge_index],1,mean)-QL0
    QLrmse=sqrt((apply(try$qlmle2[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$qlmle2[,-nonconverge_index],1,sd))^2)
    reliabias=apply(try$reliamle2[,-nonconverge_index],1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliamle2[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$reliamle2[,-nonconverge_index],1,sd))^2)
    lifembias=mean(try$lifemmle2[-nonconverge_index])-lifem0
    lifemrmse=sqrt((mean(try$lifemmle2[-nonconverge_index])-lifem0)^2+sd(try$lifemmle2[-nonconverge_index])^2)
  }
  result_df=rbind(result_df,data.frame(method="MLE2",shape=shapev[shapei],biassum=sum(abs(bias)),rmsesum=sum(rmse),
                                       QLbiassum=sum(abs(QLbias)),QLrmsesum=sum(QLrmse),shapebias=bias[1],scalebias=bias[2],
                                       thresholdbias=bias[3],shapermse=rmse[1],scalermse=rmse[2],thresholdrmse=rmse[3],
                                       noncon_prop=length(nonconverge_index)/ite,
                                       QL0.01bias=QLbias[1],QL0.05bias=QLbias[2],QL0.1bias=QLbias[3],QL0.5bias=QLbias[4],
                                       QL0.01rmse=QLrmse[1],QL0.05rmse=QLrmse[2],QL0.1rmse=QLrmse[3],QL0.5rmse=QLrmse[4],
                                       RL0.99bias=reliabias[1],RL0.95bias=reliabias[2],RL0.9bias=reliabias[3],RL0.5bias=reliabias[4],
                                       RL0.99rmse=reliarmse[1],RL0.95rmse=reliarmse[2],RL0.9rmse=reliarmse[3],RL0.5rmse=reliarmse[4],
                                       lifembias=lifembias,lifemrmse=lifemrmse))
  
  nonconverge_index=which(try$parbl[1,]==1e-5)
  if (length(nonconverge_index) == 0) {
    bias=apply(try$parbl,1,mean)-parstrue
    rmse=sqrt((apply(try$parbl,1,mean)-parstrue)^2+(apply(try$parbl,1,sd))^2)
    QLbias=apply(try$qlbl,1,mean)-QL0
    QLrmse=sqrt((apply(try$qlbl,1,mean)-QL0)^2+(apply(try$qlbl,1,sd))^2)
    reliabias=apply(try$reliabl,1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliabl,1,mean)-QL0)^2+(apply(try$reliabl,1,sd))^2)
    lifembias=mean(try$lifembl)-lifem0
    lifemrmse=sqrt((mean(try$lifembl)-lifem0)^2+sd(try$lifembl)^2)
  } else {
    bias=apply(try$parbl[,-nonconverge_index],1,mean)-parstrue
    rmse=sqrt((apply(try$parbl[,-nonconverge_index],1,mean)-parstrue)^2+(apply(try$parbl[,-nonconverge_index],1,sd))^2)
    QLbias=apply(try$qlbl[,-nonconverge_index],1,mean)-QL0
    QLrmse=sqrt((apply(try$qlbl[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$qlbl[,-nonconverge_index],1,sd))^2)
    reliabias=apply(try$reliabl[,-nonconverge_index],1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliabl[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$reliabl[,-nonconverge_index],1,sd))^2)
    lifembias=mean(try$lifembl[-nonconverge_index])-lifem0
    lifemrmse=sqrt((mean(try$lifembl[-nonconverge_index])-lifem0)^2+sd(try$lifembl[-nonconverge_index])^2)
  }
  result_df=rbind(result_df,data.frame(method="BL",shape=shapev[shapei],biassum=sum(abs(bias)),rmsesum=sum(rmse),
                                       QLbiassum=sum(abs(QLbias)),QLrmsesum=sum(QLrmse),shapebias=bias[1],scalebias=bias[2],
                                       thresholdbias=bias[3],shapermse=rmse[1],scalermse=rmse[2],thresholdrmse=rmse[3],
                                       noncon_prop=length(nonconverge_index)/ite,
                                       QL0.01bias=QLbias[1],QL0.05bias=QLbias[2],QL0.1bias=QLbias[3],QL0.5bias=QLbias[4],
                                       QL0.01rmse=QLrmse[1],QL0.05rmse=QLrmse[2],QL0.1rmse=QLrmse[3],QL0.5rmse=QLrmse[4],
                                       RL0.99bias=reliabias[1],RL0.95bias=reliabias[2],RL0.9bias=reliabias[3],RL0.5bias=reliabias[4],
                                       RL0.99rmse=reliarmse[1],RL0.95rmse=reliarmse[2],RL0.9rmse=reliarmse[3],RL0.5rmse=reliarmse[4],
                                       lifembias=lifembias,lifemrmse=lifemrmse))
  
  nonconverge_index=which(try$parbl2[1,]==1e-5)
  if (length(nonconverge_index) == 0) {
    bias=apply(try$parbl2,1,mean)-parstrue
    rmse=sqrt((apply(try$parbl2,1,mean)-parstrue)^2+(apply(try$parbl2,1,sd))^2)
    QLbias=apply(try$qlbl2,1,mean)-QL0
    QLrmse=sqrt((apply(try$qlbl2,1,mean)-QL0)^2+(apply(try$qlbl2,1,sd))^2)
    reliabias=apply(try$reliabl2,1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliabl2,1,mean)-QL0)^2+(apply(try$reliabl2,1,sd))^2)
    lifembias=mean(try$lifembl2)-lifem0
    lifemrmse=sqrt((mean(try$lifembl2)-lifem0)^2+sd(try$lifembl2)^2)
  } else {
    bias=apply(try$parbl2[,-nonconverge_index],1,mean)-parstrue
    rmse=sqrt((apply(try$parbl2[,-nonconverge_index],1,mean)-parstrue)^2+(apply(try$parbl2[,-nonconverge_index],1,sd))^2)
    QLbias=apply(try$qlbl2[,-nonconverge_index],1,mean)-QL0
    QLrmse=sqrt((apply(try$qlbl2[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$qlbl2[,-nonconverge_index],1,sd))^2)
    reliabias=apply(try$reliabl2[,-nonconverge_index],1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliabl2[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$reliabl2[,-nonconverge_index],1,sd))^2)
    lifembias=mean(try$lifembl2[-nonconverge_index])-lifem0
    lifemrmse=sqrt((mean(try$lifembl2[-nonconverge_index])-lifem0)^2+sd(try$lifembl2[-nonconverge_index])^2)
  }
  result_df=rbind(result_df,data.frame(method="BL2",shape=shapev[shapei],biassum=sum(abs(bias)),rmsesum=sum(rmse),
                                       QLbiassum=sum(abs(QLbias)),QLrmsesum=sum(QLrmse),shapebias=bias[1],scalebias=bias[2],
                                       thresholdbias=bias[3],shapermse=rmse[1],scalermse=rmse[2],thresholdrmse=rmse[3],
                                       noncon_prop=length(nonconverge_index)/ite,
                                       QL0.01bias=QLbias[1],QL0.05bias=QLbias[2],QL0.1bias=QLbias[3],QL0.5bias=QLbias[4],
                                       QL0.01rmse=QLrmse[1],QL0.05rmse=QLrmse[2],QL0.1rmse=QLrmse[3],QL0.5rmse=QLrmse[4],
                                       RL0.99bias=reliabias[1],RL0.95bias=reliabias[2],RL0.9bias=reliabias[3],RL0.5bias=reliabias[4],
                                       RL0.99rmse=reliarmse[1],RL0.95rmse=reliarmse[2],RL0.9rmse=reliarmse[3],RL0.5rmse=reliarmse[4],
                                       lifembias=lifembias,lifemrmse=lifemrmse))
  
  nonconverge_index=which(try$parmps[1,]==1e-5)
  if (length(nonconverge_index) == 0) {
    bias=apply(try$parmps,1,mean)-parstrue
    rmse=sqrt((apply(try$parmps,1,mean)-parstrue)^2+(apply(try$parmps,1,sd))^2)
    QLbias=apply(try$qlmps,1,mean)-QL0
    QLrmse=sqrt((apply(try$qlmps,1,mean)-QL0)^2+(apply(try$qlmps,1,sd))^2)
    reliabias=apply(try$reliamps,1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliamps,1,mean)-QL0)^2+(apply(try$reliamps,1,sd))^2)
    lifembias=mean(try$lifemmps)-lifem0
    lifemrmse=sqrt((mean(try$lifemmps)-lifem0)^2+sd(try$lifemmps)^2)
  } else {
    bias=apply(try$parmps[,-nonconverge_index],1,mean)-parstrue
    rmse=sqrt((apply(try$parmps[,-nonconverge_index],1,mean)-parstrue)^2+(apply(try$parmps[,-nonconverge_index],1,sd))^2)
    QLbias=apply(try$qlmps[,-nonconverge_index],1,mean)-QL0
    QLrmse=sqrt((apply(try$qlmps[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$qlmps[,-nonconverge_index],1,sd))^2)
    reliabias=apply(try$reliamps[,-nonconverge_index],1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliamps[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$reliamps[,-nonconverge_index],1,sd))^2)
    lifembias=mean(try$lifemmps[-nonconverge_index])-lifem0
    lifemrmse=sqrt((mean(try$lifemmps[-nonconverge_index])-lifem0)^2+sd(try$lifemmps[-nonconverge_index])^2)
  }
  result_df=rbind(result_df,data.frame(method="MPS",shape=shapev[shapei],biassum=sum(abs(bias)),rmsesum=sum(rmse),
                                       QLbiassum=sum(abs(QLbias)),QLrmsesum=sum(QLrmse),shapebias=bias[1],scalebias=bias[2],
                                       thresholdbias=bias[3],shapermse=rmse[1],scalermse=rmse[2],thresholdrmse=rmse[3],
                                       noncon_prop=length(nonconverge_index)/ite,
                                       QL0.01bias=QLbias[1],QL0.05bias=QLbias[2],QL0.1bias=QLbias[3],QL0.5bias=QLbias[4],
                                       QL0.01rmse=QLrmse[1],QL0.05rmse=QLrmse[2],QL0.1rmse=QLrmse[3],QL0.5rmse=QLrmse[4],
                                       RL0.99bias=reliabias[1],RL0.95bias=reliabias[2],RL0.9bias=reliabias[3],RL0.5bias=reliabias[4],
                                       RL0.99rmse=reliarmse[1],RL0.95rmse=reliarmse[2],RL0.9rmse=reliarmse[3],RL0.5rmse=reliarmse[4],
                                       lifembias=lifembias,lifemrmse=lifemrmse))
  
  nonconverge_index=which(try$parmps2[1,]==1e-5)
  if (length(nonconverge_index) == 0) {
    bias=apply(try$parmps2,1,mean)-parstrue
    rmse=sqrt((apply(try$parmps2,1,mean)-parstrue)^2+(apply(try$parmps2,1,sd))^2)
    QLbias=apply(try$qlmps2,1,mean)-QL0
    QLrmse=sqrt((apply(try$qlmps2,1,mean)-QL0)^2+(apply(try$qlmps2,1,sd))^2)
    reliabias=apply(try$reliamps2,1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliamps2,1,mean)-QL0)^2+(apply(try$reliamps2,1,sd))^2)
    lifembias=mean(try$lifemmps2)-lifem0
    lifemrmse=sqrt((mean(try$lifemmps2)-lifem0)^2+sd(try$lifemmps2)^2)
  } else {
    bias=apply(try$parmps2[,-nonconverge_index],1,mean)-parstrue
    rmse=sqrt((apply(try$parmps2[,-nonconverge_index],1,mean)-parstrue)^2+(apply(try$parmps2[,-nonconverge_index],1,sd))^2)
    QLbias=apply(try$qlmps2[,-nonconverge_index],1,mean)-QL0
    QLrmse=sqrt((apply(try$qlmps2[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$qlmps2[,-nonconverge_index],1,sd))^2)
    reliabias=apply(try$reliamps2[,-nonconverge_index],1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliamps2[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$reliamps2[,-nonconverge_index],1,sd))^2)
    lifembias=mean(try$lifemmps2[-nonconverge_index])-lifem0
    lifemrmse=sqrt((mean(try$lifemmps2[-nonconverge_index])-lifem0)^2+sd(try$lifemmps2[-nonconverge_index])^2)
  }
  result_df=rbind(result_df,data.frame(method="MPS2",shape=shapev[shapei],biassum=sum(abs(bias)),rmsesum=sum(rmse),
                                       QLbiassum=sum(abs(QLbias)),QLrmsesum=sum(QLrmse),shapebias=bias[1],scalebias=bias[2],
                                       thresholdbias=bias[3],shapermse=rmse[1],scalermse=rmse[2],thresholdrmse=rmse[3],
                                       noncon_prop=length(nonconverge_index)/ite,
                                       QL0.01bias=QLbias[1],QL0.05bias=QLbias[2],QL0.1bias=QLbias[3],QL0.5bias=QLbias[4],
                                       QL0.01rmse=QLrmse[1],QL0.05rmse=QLrmse[2],QL0.1rmse=QLrmse[3],QL0.5rmse=QLrmse[4],
                                       RL0.99bias=reliabias[1],RL0.95bias=reliabias[2],RL0.9bias=reliabias[3],RL0.5bias=reliabias[4],
                                       RL0.99rmse=reliarmse[1],RL0.95rmse=reliarmse[2],RL0.9rmse=reliarmse[3],RL0.5rmse=reliarmse[4],
                                       lifembias=lifembias,lifemrmse=lifemrmse))
  
  nonconverge_index=which(try$parple[1,]==1e-5)
  if (length(nonconverge_index) == 0) {
    bias=apply(try$parple,1,mean)-parstrue
    rmse=sqrt((apply(try$parple,1,mean)-parstrue)^2+(apply(try$parple,1,sd))^2)
    QLbias=apply(try$qlple,1,mean)-QL0
    QLrmse=sqrt((apply(try$qlple,1,mean)-QL0)^2+(apply(try$qlple,1,sd))^2)
    reliabias=apply(try$reliaple,1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliaple,1,mean)-QL0)^2+(apply(try$reliaple,1,sd))^2)
    lifembias=mean(try$lifemple)-lifem0
    lifemrmse=sqrt((mean(try$lifemple)-lifem0)^2+sd(try$lifemple)^2)
  } else {
    bias=apply(try$parple[,-nonconverge_index],1,mean)-parstrue
    rmse=sqrt((apply(try$parple[,-nonconverge_index],1,mean)-parstrue)^2+(apply(try$parple[,-nonconverge_index],1,sd))^2)
    QLbias=apply(try$qlple[,-nonconverge_index],1,mean)-QL0
    QLrmse=sqrt((apply(try$qlple[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$qlple[,-nonconverge_index],1,sd))^2)
    reliabias=apply(try$reliaple[,-nonconverge_index],1,mean)-reliability0
    reliarmse=sqrt((apply(try$reliaple[,-nonconverge_index],1,mean)-QL0)^2+(apply(try$reliaple[,-nonconverge_index],1,sd))^2)
    lifembias=mean(try$lifemple[-nonconverge_index])-lifem0
    lifemrmse=sqrt((mean(try$lifemple[-nonconverge_index])-lifem0)^2+sd(try$lifemple[-nonconverge_index])^2)
  }
  result_df=rbind(result_df,data.frame(method="PLE",shape=shapev[shapei],biassum=sum(abs(bias)),rmsesum=sum(rmse),
                                       QLbiassum=sum(abs(QLbias)),QLrmsesum=sum(QLrmse),shapebias=bias[1],scalebias=bias[2],
                                       thresholdbias=bias[3],shapermse=rmse[1],scalermse=rmse[2],thresholdrmse=rmse[3],
                                       noncon_prop=length(nonconverge_index)/ite,
                                       QL0.01bias=QLbias[1],QL0.05bias=QLbias[2],QL0.1bias=QLbias[3],QL0.5bias=QLbias[4],
                                       QL0.01rmse=QLrmse[1],QL0.05rmse=QLrmse[2],QL0.1rmse=QLrmse[3],QL0.5rmse=QLrmse[4],
                                       RL0.99bias=reliabias[1],RL0.95bias=reliabias[2],RL0.9bias=reliabias[3],RL0.5bias=reliabias[4],
                                       RL0.99rmse=reliarmse[1],RL0.95rmse=reliarmse[2],RL0.9rmse=reliarmse[3],RL0.5rmse=reliarmse[4],
                                       lifembias=lifembias,lifemrmse=lifemrmse))
  
  
  result_df_list[[shapei]]=result_df
}
result_df_merge=do.call(rbind,result_df_list)
View(result_df_merge)
library(openxlsx)
write.xlsx(result_df_merge, file = "gamma_n_20_four methods.csv", rowNames = FALSE)
save.image(file="gamma_n_20_four methods.RData")
cat("\ntotal time:", difftime(end_time,start_time, units = "mins"), "mins\n")















