library(tidyverse)
library(dplyr)
library(rsvd)
library(matrixcalc)
library(MASS)
library(foreach)
library(doParallel)
library(doRNG)
library(parallel)
library(tdROC)
library(survival)

#dirP<-"/Users/ma521/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/"
dirP<-"/hpc/home/ma521/hdld_surv_jm/"
source(paste(dirP,"high_dimensional_long_surv_jm.R",sep=""))
source(paste(dirP,"supFTSVD.R",sep=""))
#source(paste(dirP,"supFTSVD_JMRE.R",sep=""))
# Number of subjects
n<-300

# feature dimension 

pdim<-c(500)

# No of rank-1 components
r<-3
model_rank<-r

# Weight for components
lmd_val<-c(5.20,4.80,3.35) 

# noise variance in tensor
Tau2<-c(0.1)

# Parameters for supervised component
Eta2<-c(1,1,1)

#diag(Eta2)<-c(1.85,2.50,1.3,1.60,1.25,3.5)
gam_par<-matrix(c(2.2,1.6,2.5,3.6,-1.50,2.6),ncol=3)
gam_par<-apply(gam_par,2,function(u){u/sqrt(sum(u^2))})

# Grid of Time points
nres<-101
Time<-seq(0,1,length.out=nres)

# Singular Function
PhiFunc<-list(function(x){(8+(6*x^8)-(3*x^2)-(4*x^3))/sqrt(45.61)},
              function(x){(10*x^2/exp(x^5))/(sqrt(10-10*exp(-2)))},
              function(x){sqrt(2)*sin(2.5*pi*x)})

PhiF<-sapply(1:r,function(k){PhiFunc[[k]](Time)})

# Feature loading
set.seed(pdim)
Bval<<-sapply(1:r, function(b){runif(pdim)})
bval<<-Bval*outer(rep(1,pdim),1/apply(Bval,2,norm,type="2"))

# controlling the signal ot noise ratio
cmp_var<-((lmd_val^2)*Eta2)
CmpV<-Reduce(`+`,lapply(1:r,function(k){cmp_var[k]*(outer(bval[,k],PhiF[,k]))^2}))
mean(apply(CmpV,1,mean))/Tau2 
min(apply(CmpV,1,mean))/Tau2


# Survival parameters
alp_par<-matrix(c(3.50,2.60,-0.15,-0.20,0.15),ncol=1) 
lmd<-0.05
## Subject-specific covariates # ADAS13
set.seed(n) 
Vmat<<-cbind(round(runif(n),2),round(rbeta(n,2.5,1.5),2))
EAval<-Vmat%*%gam_par
miv_subE<-sapply(Eta2,function(u){rnorm(n,mean=0,sd=sqrt(u))}) 
ZetaSL<-EAval+miv_subE

# Training data set
# Survival data generation
set.seed(seed_n*24) # for both survival and longitudinal data
surv_time<-(-log(runif(n,0,1)))/(lmd*exp(as.numeric(cbind(Vmat,ZetaSL)%*%(alp_par))))
summary(surv_time)
cen_llim<-0
cen_ulim<-5
censor_time<-runif(n,cen_llim,cen_ulim)
survT<-apply(cbind(surv_time,censor_time),1,min)
#summary(survT)
cenI<-apply(cbind(surv_time,censor_time),1,which.min)-1
cen_indx<-which(survT>1)
survT[cen_indx]<-1
cenI[cen_indx]<-1
cenP<-mean(cenI)
cenP
sum(as.numeric(survT==1))
#summary(survT[cenI==1])
sur_obj<-Surv(survT,1-cenI)


# training m-omics data
m_i<-sample(5:10,n,replace = TRUE)
bl_time<-sapply(survT, function(u){runif(1,0,u)})
#tr_obsTIME<-lapply(1:n, function(i){c(bl_time[i],sort(runif(m_i[i]-2,bl_time[i],survT[i])),survT[i])})
tr_obsTIME<-lapply(1:n, function(i){c(0,sort(runif(m_i[i]-2,0,survT[i])),survT[i])})
gen_dataCFS<-omics_data_gen_surv(m_i = m_i,Zeta= ZetaSL,obsTIME = tr_obsTIME,Xi = bval,
                                 PsiF = PhiFunc,sing_val = lmd_val,Data_Var = Tau2,surv_time=survT)


# plot(gen_dataCFS$data[[1]][1,],gen_dataCFS$data[[1]][2,],type="l",ylim=c(0,3))
# for(i in 2:n){
#   lines(gen_dataCFS$data[[i]][1,],gen_dataCFS$data[[i]][2,])
# }

# True R^2
comps<-do.call(rbind,lapply(1:n, function(i){
  sapply(1:r, function(k){
    as.numeric((lmd_val[k]*gen_dataCFS$zeta[i,k])*outer(bval[,k],gen_dataCFS$singF[[i]][,k]))
  })
}))

genY<-do.call(c,lapply(1:n,function(i){
  as.numeric(gen_dataCFS$data[[i]][-1,])
}))

tr2<-summary(lm(genY~comps-1))$r.squared
summary(lm(genY~comps[,1]-1))$r.squared
summary(lm(genY~comps[,2]-1))$r.squared
summary(lm(genY~comps[,3]-1))$r.squared

pt0<-as.numeric(proc.time())[3]
fit_model<-supFTSVD_JM(datlist = gen_dataCFS$data, 
                       response=Vmat, interval = c(0,1), r = model_rank,resolution=50, CVPhi=TRUE, K=5, cvT=5, smooth=round(exp(seq(-7,1,length.out=20)),3),
                       surv_time=survT,
                       censor_status=cenI,
                       maxiter=100, epsilon=1e-5,KInd=NULL,rsvd_seed=100,conv_criteria = "cond_lik",
                       survX = Vmat,scale = TRUE,constant_hazard = TRUE)
cmp_time<-c(as.numeric(proc.time())[3]-pt0)

fit_jll<-as.numeric(fit_model$cnvrgnt[nrow(fit_model$cnvrgnt),ncol(fit_model$cnvrgnt)])

### combined omics data
pt4<-as.numeric(proc.time())[3]

tfit_model<-supFTSVD(datlist = gen_dataCFS$data, 
                     response=Vmat, interval = c(0,1), r = model_rank,resolution=50, CVPhi=TRUE, K=5, cvT=5, smooth=round(exp(seq(-7,1,length.out=20)),3),
                     maxiter=10, epsilon=1e-5,KInd=NULL)
cmp_time<-c(cmp_time,as.numeric(proc.time())[3]-pt4)

## Survival obj based on two-stage estimation
coxOBJ<-coxph_mle_fllik(surv_time=survT,censor_status=cenI,Xmat=cbind(Vmat,sapply(1:ncol(tfit_model$A.hat),function(k){tfit_model$A.hat[,k]/sqrt(tfit_model$Sigma2R[k])})),
                        max_iter=100,cnv_type = "log_lik",cnv_crit=1e-6,constant_hazard = TRUE)

tfit_jll<-as.numeric(tfit_model$cnvrgnt[nrow(tfit_model$cnvrgnt),ncol(tfit_model$cnvrgnt)])+coxOBJ$cnvrgn[nrow(coxOBJ$cnvrgn),3]

## singular function differences
phiDIF<-t(sapply(1:model_rank, function(k){
  estP<-approx(fit_model$time,fit_model$Phi.hat[,k],xout = Time,rule = 2)$y
  estP1<-approx(tfit_model$time,tfit_model$Phi.hat[,k],xout = Time,rule = 2)$y
  c(min(mean((estP-PhiF[,k])^2),mean((-estP-PhiF[,k])^2)),
    min(mean((estP1-PhiF[,k])^2),mean((-estP1-PhiF[,k])^2)))
}))



featDIF<-t(sapply(1:model_rank, function(k){
  estP<-fit_model$B.hat[,k]
  estP1<-tfit_model$B.hat[,k]
  c(min(mean((estP-bval[,k])^2),mean((-estP-bval[,k])^2)),
    min(mean((estP1-bval[,k])^2),mean((-estP1-bval[,k])^2)))
}))


mseR2<-c(sqrt(sum((fit_model$accum.r.square[model_rank]-tr2)^2)),
         sqrt(sum((tfit_model$accum.r.square[model_rank]-tr2)^2)))

coxll<-c(fit_model$cnvrgnt[nrow(fit_model$cnvrgnt)-1,ncol(fit_model$cnvrgnt)-1],coxOBJ$cnvrgn[nrow(coxOBJ$cnvrgn),3])

tr_res<-data.frame("JM"=c(as.numeric(rbind(phiDIF,featDIF*1000,mseR2,coxll,fit_jll)[,1]),fit_model$CoxPH_obj$coefficients),
                   "Two"=c(as.numeric(rbind(phiDIF,featDIF*1000,mseR2,coxll,tfit_jll)[,2]),coxOBJ$CoxPH_obj$coefficients),
                   "Comps"=c("k=1","k=2","k=3","k=1","k=2","k=3","all","survLL","jointLL",rep("surv_coef",ncol(Vmat)+r)),
                   "Stat"=c("SingFunc","SingFunc","SingFunc","FeatLoad","FeatLoad","FeatLoad","difR2","sLL","jointLL",rep("surv_coef",ncol(Vmat)+r)))

# Integrated AUC and BS

# test data
lm_time<-c(0.2,0.3,0.4,0.5)

pt5<-as.numeric(proc.time())[3]
prd_res<-t(sapply(1:length(lm_time),function(l){
  # test data: dynamic framework
  # test data
  ssurvT<-NULL
  scenI<-NULL
  sZetaSL<-NULL
  svalid_Vmat<-NULL
  set.seed(n+25)
  vN<-100
  repeat{
    tvalid_Vmat<-cbind(round(runif(vN),2),round(rbeta(vN,2.5,1.5),2))
    tEAval<-tvalid_Vmat%*%gam_par
    
    tmiv_subE<-sapply(Eta2,function(u){rnorm(vN,mean=0,sd=sqrt(u))}) 
    tZetaSL<-tEAval+tmiv_subE
    
    tsurv_time<-(-log(runif(vN,0,1)))/(lmd*exp(as.numeric(cbind(tvalid_Vmat,tZetaSL)%*%matrix(alp_par))))
    #summary(surv_time)
    tcensor_time<-runif(vN,cen_llim,cen_ulim)
    tsurvT<-apply(cbind(tsurv_time,tcensor_time),1,min)
    tcenI<-apply(cbind(tsurv_time,tcensor_time),1,which.min)-1
    tcen_indx<-which(tsurvT>1)
    tsurvT[tcen_indx]<-1
    tcenI[tcen_indx]<-1
    f_ind<-which(tsurvT>lm_time[l])
    ssurvT<-c(ssurvT,tsurvT[f_ind])
    scenI<-c(scenI,tcenI[f_ind])
    sZetaSL<-rbind(sZetaSL,tZetaSL[f_ind,])
    svalid_Vmat<-rbind(svalid_Vmat,tvalid_Vmat[f_ind,])
    if(length(ssurvT)>=vN)
      break
  }
  
  vsam<-sample(1:length(ssurvT),vN,replace = FALSE)
  vsurvT<-ssurvT[vsam]
  vcenI<-scenI[vsam]
  vcenP<-mean(vcenI)
  valid_Vmat<-svalid_Vmat[vsam,]
  vZetaSL<-sZetaSL[vsam,]
  # test m-omics data
  vm_i<-rep(ifelse(l==1,10,ifelse(l==2,15,20)),vN)
  #vobsTIME<-lapply(vm_i,function(u){c(0,sort(runif(8,0,lm_time[1])),lm_time[1],sort(runif(4,lm_time[1],lm_time[2])),lm_time[2],sort(runif(4,lm_time[2],lm_time[3])),lm_time[3])})
  vobsTIME<-lapply(vm_i,function(u){c(0,sort(runif(u-2,0,lm_time[l])),lm_time[l])})
  vgen_dataCFS<-omics_data_gen_surv(m_i = vm_i,Zeta= vZetaSL,Xi = bval,
                                    PsiF = PhiFunc,sing_val = lmd_val,obsTIME = vobsTIME,Data_Var = Tau2,surv_time=rep(max(lm_time),100))
  
  
  #list(vsurvT,vcenI,vcenP,vgen_dataCFS)
  
  vld_Fdata<-lapply(1:length(vm_i),function(i){
    Reduce(`+`,lapply(1:r, function(k){
      (vgen_dataCFS$zeta[i,k]*lmd_val[k])*outer(bval[,k],as.numeric(PhiF[,k]))
    }))
  })
  
  #tsvl_dat<-lapply(vgen_dataCFS$data,function(u){as.matrix(u[,which(u[1,]<=lm_time[l])])})
  
  
  prd_msF<-predict.supFTSVD_JM(obj=fit_model,designM=valid_Vmat,mc_ss=1000,surv_time=rep(lm_time[l],vN),
                               censor_status=rep(1,vN),survX = valid_Vmat,new_dat=vgen_dataCFS$data,TimeG=Time)
  
  tprd_msF<-predict.supFTSVD(obj = tfit_model,designM = valid_Vmat,new_dat = vgen_dataCFS$data,TimeG = Time)
  # two-stage conditional probability
  bz_fit<-coxOBJ$CoxPH_obj$cum_hazard
  bz_grid<-approx(x=bz_fit$time,y=bz_fit$hazard,xout=Time,rule=2)$y
  tsg_par<-as.numeric(coxOBJ$CoxPH_obj$coefficients)
  tpred_survP<-sapply(exp(cbind(valid_Vmat,sapply(1:ncol(tprd_msF$subEFF),function(k){tprd_msF$subEFF[,k]/sqrt(tfit_model$Sigma2R[k])}))%*%matrix(tsg_par)),function(u){exp(-u*bz_grid)})
  
  
  prp_pe<-mean(colMeans(sapply(1:length(vm_i),function(i){
    rowMeans((vld_Fdata[[i]]-prd_msF$subTRJ[[i]])^2)
  })))
  
  tprp_pe<-mean(colMeans(sapply(1:length(vm_i),function(i){
    rowMeans((vld_Fdata[[i]]-tprd_msF$subTRJ[[i]])^2)
  })))
  
  del_time<-seq(lm_time[l],lm_time[l]+0.25,length.out=26) 
  del_pos<-sapply(del_time, function(u){which.min((u-Time)^2)})
  est_cond_prob<-sapply(1:vN, function(i){
    a<-prd_msF$Surv_Prob[del_pos,i]/prd_msF$Surv_Prob[del_pos[1],i]
    1-a[-1]
  })
  
  test_cond_prob<-sapply(1:vN, function(i){
    a<-tpred_survP[del_pos,i]/tpred_survP[del_pos[1],i]
    1-a[-1]
  })
  
  trw<-exp(as.numeric(cbind(valid_Vmat,vZetaSL)%*%matrix(alp_par))) #*c(1,1,lmd_val)
  tsurvP<-sapply(trw,function(u){exp(-u*lmd*Time)})
  true_cond_prob<-sapply(1:vN, function(i){
    a<-tsurvP[del_pos,i]/tsurvP[del_pos[1],i]
    1-a[-1]
  })
  
  
  ROC.est <- sapply(seq_len(nrow(est_cond_prob)),function(i){
    fit<-tdROC(X = est_cond_prob[i,], Y = vsurvT,
               delta = 1-vcenI, tau = del_time[i+1],
               span = 0.1, alpha = 0.05,
               n.grid = 1000, cut.off = 0.5)
    c(fit$main_res$AUC.empirical,as.numeric(fit$calibration_res[1]))
  })
  
  twROC.est <- sapply(seq_len(nrow(test_cond_prob)),function(i){
    fit<-tdROC(X = test_cond_prob[i,], Y = vsurvT,
               delta = 1-vcenI, tau = del_time[i+1],
               span = 0.1, alpha = 0.05,
               n.grid = 1000, cut.off = 0.5)
    c(fit$main_res$AUC.empirical,as.numeric(fit$calibration_res[1]))
  })
  
  
  trROC.est <- sapply(seq_len(nrow(true_cond_prob)),function(i){
    fit<-tdROC(X = true_cond_prob[i,], Y = vsurvT,
               delta = 1-vcenI, tau = del_time[i+1],
               span = 0.1, alpha = 0.05,
               n.grid = 1000, cut.off = 0.5)
    c(fit$main_res$AUC.empirical,as.numeric(fit$calibration_res[1]))
  })
  
  c(lm_time[l],rowMeans(rbind(ROC.est[1,],twROC.est[1,],trROC.est[1,],
                              ROC.est[2,],twROC.est[2,],trROC.est[2,]),na.rm=TRUE),prp_pe,tprp_pe,
    mean(apply((est_cond_prob-true_cond_prob)^2,2,mean)),
    mean(apply((test_cond_prob-true_cond_prob)^2,2,mean)))
}))
cmp_time<-c(cmp_time,as.numeric(proc.time())[3]-pt5)

colnames(prd_res)<-c("LM","prp_auc","tw_auc","tr_auc",
                     "prp_bs","tw_bs","tr_bs",
                     "prp_err","sprp_err","prp_scrv","ts_scrv")

#prd_res
#save_dir<-"/Users/ma521/Academic/PostDocWorks/Research_PostDocs/Longitudinal_Multi_Omics/momics_simstudy_codes/momics_ss_august19_24/"
save_dir<-"/cwork/ma521/hdld_surv_jm/SimResults/"  
write.table(tr_res,file = paste(save_dir,"hdld_tr_res","_",n,"_",seed_n,".txt",sep=""))
write.table(prd_res,file = paste(save_dir,"hdld_prd_res","_",n,"_",seed_n,".txt",sep=""))
write.table(cmp_time,file = paste(save_dir,"hdld_cmp_time","_",n,"_",seed_n,".txt",sep=""))

