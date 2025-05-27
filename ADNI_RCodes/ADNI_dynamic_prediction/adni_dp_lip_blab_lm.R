# necessary packages
library(tidyverse)
library(dplyr)
library(rsvd)
library(matrixcalc)
library(survival)
library(MASS)
library(tdROC)
# directory for for codes
#code_dir<-"~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/"
code_dir<-"/hpc/home/ma521/hdld_surv_jm_adni/"
source(paste(code_dir,"high_dimensional_long_surv_jm.R",sep=""))
source(paste(code_dir,"supFTSVD.R",sep = ""))

# data directory
#data_dir<-"~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/adni_data_files_hdld_surv_jm/"
#res_dir<-"~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/adni_da_results/"

data_dir<-"/hpc/home/ma521/hdld_surv_jm_adni/adni_data_files_processed/"
res_dir<-"/cwork/ma521/adni_da_hdld_jm/adni_dp_res1year/"


# loading full survial and bile acid data
file_name<-"lipidBL_surv_full_17_797.csv"
surv_long_data<-read.csv(paste(data_dir,file_name,sep="")) %>%
  dplyr::select_if(is.numeric) %>%
  group_by(RID) %>%
  mutate(n_obs=n()) %>%
  filter(n_obs>1) %>%
  mutate(censor_status=ifelse(surv_time>max(visit_time),1,censor_status)) %>%
  mutate(surv_time=max(visit_time)) %>%
  ungroup()

# censoring status for 51 patients with dementia got altered when we consider observed survival time
# and the last visit time are equal


# construction of meta data (first 16 columns are survival and clinical variables)
meta_data<-surv_long_data %>%
  dplyr::select(1:16) %>%
  distinct() 


# training-test ID
#lm_time<-c(18)
test_id_pop<-meta_data %>%
  group_by(RID) %>%
  filter(surv_time>=max(lm_time) & min(visit_time)<=min(lm_time)) %>%
  ungroup() %>%
  dplyr::select(RID) %>%
  distinct() %>%
  as.matrix() %>%
  as.numeric()



vN<-floor(length(unique(meta_data$RID))*0.20)
set.seed(seed_n)
test_id<-sample(x = test_id_pop,size = vN,replace = FALSE)

# id_vectors for longitudinal data
tr_id<-surv_long_data %>%
  filter(!RID %in% test_id) %>%
  dplyr::select(RID) %>%
  as.matrix() %>%
  as.numeric()

ts_id<-surv_long_data %>%
  filter(RID %in% test_id) %>%
  dplyr::select(RID) %>%
  as.matrix() %>%
  as.numeric()


#### training and test data
tr_meta_data<-meta_data  %>%
  filter(!RID %in% test_id)

ts_meta_data<-meta_data  %>%
  filter(RID %in% test_id)

tr_surv_long_data<-surv_long_data %>%
  filter(!RID %in% test_id)

ts_surv_long_data<-surv_long_data %>%
  filter(RID %in% test_id)


# design matrix for supervised decomposition  
tr_res_vec<-tr_meta_data  %>%
  mutate(Interp=1) %>%
  dplyr::select(RID,Interp,AGE,PTGENDER,PTEDUCAT,APOE4) %>%
  distinct() %>%
  dplyr::select(-RID) %>%
  as.matrix()

ts_res_vec<-ts_meta_data  %>%
  mutate(Interp=1) %>%
  dplyr::select(RID,Interp,AGE,PTGENDER,PTEDUCAT,APOE4) %>%
  distinct() %>%
  dplyr::select(-RID) %>%
  as.matrix()


# survival covariates
tr_surv_pred<-tr_meta_data  %>%
  dplyr::select(RID,AGE,PTGENDER,PTEDUCAT,APOE4) %>%
  distinct() %>%
  dplyr::select(-RID) %>%
  as.matrix()

ts_surv_pred<-ts_meta_data  %>%
  dplyr::select(RID,AGE,PTGENDER,PTEDUCAT,APOE4) %>%
  distinct() %>%
  dplyr::select(-RID) %>%
  as.matrix()


# survival data 
tr_surv_tdat<-tr_meta_data %>%
  dplyr::select(RID,surv_time,censor_status) %>%
  distinct() 

ts_surv_tdat<-ts_meta_data %>%
  dplyr::select(RID,surv_time,censor_status) %>%
  distinct() 


# omics data preparation
tr_omics_datF<-tr_surv_long_data %>%
  dplyr::select(-c(1:16))


tr_ls_omics_dat<-split.data.frame(tr_omics_datF,tr_id)
tr_ls_obs_dat<-split(tr_surv_long_data$visit_time,tr_surv_long_data$RID)

tr_n_sub<-length(tr_ls_obs_dat)

tr_omics_dat<-lapply(1:tr_n_sub, function(i){
  xx<-t(tr_ls_omics_dat[[i]])
  yy<-ifelse(xx==0,0.01,xx)
  rbind(tr_ls_obs_dat[[i]],log(yy))
})

#load(file = paste(res_dir,"cv_res_lipid_blab.RData",sep=""))
load(file = paste(data_dir,"cv_res_lipid_blab.RData",sep=""))

# fitting of the joint model
cpTime<-NULL
cp_time<-proc.time()
adni_blab_fobj_jm<-supFTSVD_JM(datlist = tr_omics_dat, 
                               response=tr_res_vec, 
                               interval = NULL, 
                               r = cv_res_lipid_blab$opt_r,
                               resolution=50, 
                               CVPhi=TRUE, 
                               K=5, 
                               cvT=5, 
                               smooth=exp(seq(-5,2,length.out=20)),
                               maxiter=200, 
                               epsilon=1e-6,
                               KInd=NULL,
                               rsvd_seed=100,
                               conv_criteria = "cond_like",
                               surv_time = tr_surv_tdat$surv_time,
                               censor_status = tr_surv_tdat$censor_status,
                               survX = tr_surv_pred,
                               constant_hazard=FALSE)
cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])
save(adni_blab_fobj_jm,file = paste(res_dir,"adni_blab_fobj_jm_",seed_n,".RData",sep = ""))


#load(paste("~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/adni_da_results/adni_dp_results/adni_blab_fobj_jm_",seed_n+100,".RData",sep=""))

# dynamic prediction
# omics data preparation
ts_omics_datF<-ts_surv_long_data %>%
  dplyr::select(-c(1:16))


ts_ls_omics_dat<-split.data.frame(ts_omics_datF,ts_id)
ts_ls_obs_dat<-split(ts_surv_long_data$visit_time,ts_surv_long_data$RID)

ts_n_sub<-length(ts_ls_obs_dat)

# Prediction grid
time_grid<-round(seq(0,max(tr_meta_data$surv_time),length=101),3)[-1]

# prediction
cp_time<-proc.time()
auc_bs_jm<-sapply(lm_time, function(u){
  set.seed(seed_n+which(lm_time==u))
  ts_omics_dat<-lapply(1:ts_n_sub, function(i){
    lm_ind<-which(ts_ls_obs_dat[[i]]<=u)
    xx<-t(ts_ls_omics_dat[[i]][lm_ind,])
    yy<-ifelse(xx==0,0.01,xx)
    if(length(lm_ind)==1){
      rbind(ts_ls_obs_dat[[i]][lm_ind],as.matrix(log(yy)))
    } else{
      rbind(ts_ls_obs_dat[[i]][lm_ind],log(yy))
    }
  })
  
  prd_msF<-predict.supFTSVD_JM(obj=adni_blab_fobj_jm,
                               designM=ts_res_vec,
                               mc_ss=1000,
                               surv_time=rep(u,vN),
                               censor_status=rep(1,vN),
                               survX = ts_surv_pred,
                               new_dat=ts_omics_dat,
                               TimeG=time_grid)
  
  del_time<-seq(u,u+12,by=0.25)
  del_pos<-sapply(del_time, function(u){which.min((u-time_grid)^2)})
  est_cond_prob<-sapply(1:vN, function(i){
    a<-prd_msF$Surv_Prob[del_pos,i]/prd_msF$Surv_Prob[del_pos[1],i]
    1-a[-1]
  })
  
  ROC.est <- sapply(seq_len(nrow(est_cond_prob)),function(i){
    fit<-tdROC(X = est_cond_prob[i,], Y = ts_surv_tdat$surv_time,
               delta = 1-ts_surv_tdat$censor_status, tau = del_time[i+1],
               span = 0.1, alpha = 0.05,
               n.grid = 1000, cut.off = 0.5)
    c(fit$main_res$AUC.empirical,as.numeric(fit$calibration_res[1]))
  })
  
  rowMeans(ROC.est,na.rm = TRUE)
})
cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])

# fitting of the joint model without direct effect
cp_time<-proc.time()
adni_blab_fobj_jm1<-supFTSVD_JM(datlist = tr_omics_dat, 
                                response=tr_res_vec, 
                                interval = NULL, 
                                r = cv_res_lipid_blab$opt_r,
                                resolution=50, 
                                CVPhi=TRUE, 
                                K=5, 
                                cvT=5, 
                                smooth=exp(seq(-5,2,length.out=20)),
                                maxiter=200, 
                                epsilon=1e-6,
                                KInd=NULL,
                                rsvd_seed=100,
                                conv_criteria = "cond_like",
                                surv_time = tr_surv_tdat$surv_time,
                                censor_status = tr_surv_tdat$censor_status,
                                survX = NULL,
                                constant_hazard=FALSE)
cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])
save(adni_blab_fobj_jm1,file = paste(res_dir,"adni_blab_fobj_jm1_",seed_n,".RData",sep = ""))

# prediction
cp_time<-proc.time()
auc_bs_jm1<-sapply(lm_time, function(u){
  set.seed(seed_n+which(lm_time==u))
  ts_omics_dat<-lapply(1:ts_n_sub, function(i){
    lm_ind<-which(ts_ls_obs_dat[[i]]<=u)
    xx<-t(ts_ls_omics_dat[[i]][lm_ind,])
    yy<-ifelse(xx==0,0.01,xx)
    if(length(lm_ind)==1){
      rbind(ts_ls_obs_dat[[i]][lm_ind],as.matrix(log(yy)))
    } else{
      rbind(ts_ls_obs_dat[[i]][lm_ind],log(yy))
    }
  })
  
  prd_msF<-predict.supFTSVD_JM(obj=adni_blab_fobj_jm1,
                               designM=ts_res_vec,
                               mc_ss=1000,
                               surv_time=rep(u,vN),
                               censor_status=rep(1,vN),
                               survX = NULL,
                               new_dat=ts_omics_dat,
                               TimeG=time_grid)
  
  del_time<-seq(u,u+12,by=0.25)
  del_pos<-sapply(del_time, function(u){which.min((u-time_grid)^2)})
  est_cond_prob<-sapply(1:vN, function(i){
    a<-prd_msF$Surv_Prob[del_pos,i]/prd_msF$Surv_Prob[del_pos[1],i]
    1-a[-1]
  })
  
  ROC.est <- sapply(seq_len(nrow(est_cond_prob)),function(i){
    fit<-tdROC(X = est_cond_prob[i,], Y = ts_surv_tdat$surv_time,
               delta = 1-ts_surv_tdat$censor_status, tau = del_time[i+1],
               span = 0.1, alpha = 0.05,
               n.grid = 1000, cut.off = 0.5)
    c(fit$main_res$AUC.empirical,as.numeric(fit$calibration_res[1]))
  })
  
  rowMeans(ROC.est,na.rm = TRUE)
})
cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])


## fitting of two-stage model
cp_time<-proc.time()
adni_blab_fobj_ts<-supFTSVD(datlist = tr_omics_dat, 
                            response=tr_res_vec,
                            interval = NULL,
                            r = cv_res_lipid_blab$opt_r,
                            resolution=50, 
                            CVPhi=TRUE, 
                            K=5, 
                            cvT=5, 
                            smooth=exp(seq(-5,2,length.out=20)),
                            maxiter=200, 
                            epsilon=1e-6,
                            KInd=NULL,
                            rsvd_seed=100,
                            conv_criteria = "cond_like")

cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])
save(adni_blab_fobj_ts,file = paste(res_dir,"adni_blab_fobj_ts_",seed_n,".RData",sep = ""))

#load(paste("~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/adni_da_results/adni_dp_results/adni_blab_fobj_ts_",seed_n+100,".RData",sep=""))

cp_time<-proc.time()

coxOBJ<-coxph_mle_fllik(surv_time=tr_surv_tdat$surv_time,censor_status=tr_surv_tdat$censor_status,Xmat=cbind(tr_surv_pred,sapply(1:ncol(adni_blab_fobj_ts$A.hat),function(k){adni_blab_fobj_ts$A.hat[,k]/sqrt(adni_blab_fobj_ts$Sigma2R[k])})),max_iter=100,cnv_type = "LL",cnv_crit=1e-6,constant_hazard = FALSE)
save(coxOBJ,file = paste(res_dir,"coxOBJ_ts_",seed_n,".RData",sep = ""))


auc_bs_ts<-sapply(lm_time, function(u){
  ts_omics_dat<-lapply(1:ts_n_sub, function(i){
    lm_ind<-which(ts_ls_obs_dat[[i]]<=u)
    xx<-t(ts_ls_omics_dat[[i]][lm_ind,])
    yy<-ifelse(xx==0,0.01,xx)
    if(length(lm_ind)==1){
      rbind(ts_ls_obs_dat[[i]][lm_ind],as.matrix(log(yy)))
    } else{
      rbind(ts_ls_obs_dat[[i]][lm_ind],log(yy))
    }
  })
  
  tprd_msF<-predict.supFTSVD(obj = adni_blab_fobj_ts,designM = ts_res_vec,new_dat = ts_omics_dat,TimeG = time_grid)
  # two-stage conditional probability
  bz_fit<-coxOBJ$CoxPH_obj$cum_hazard
  bz_grid<-approx(x=bz_fit$time,y=bz_fit$hazard,xout=time_grid,rule=2)$y
  tsg_par<-as.numeric(coxOBJ$CoxPH_obj$coefficients)
  tpred_survP<-sapply(exp(cbind(ts_surv_pred,sapply(1:ncol(adni_blab_fobj_ts$A.hat),function(k){tprd_msF$subEFF[,k]/sqrt(adni_blab_fobj_ts$Sigma2R[k])}))%*%matrix(tsg_par)),function(u){exp(-u*bz_grid)})
  
  del_time<-seq(u,u+12,by=0.25)
  del_pos<-sapply(del_time, function(u){which.min((u-time_grid)^2)})
  
  test_cond_prob<-sapply(1:vN, function(i){
    a<-tpred_survP[del_pos,i]/tpred_survP[del_pos[1],i]
    1-a[-1]
  })
  
  
  twROC.est <- sapply(seq_len(nrow(test_cond_prob)),function(i){
    fit<-tdROC(X = test_cond_prob[i,], Y = ts_surv_tdat$surv_time,
               delta = 1-ts_surv_tdat$censor_status, tau = del_time[i+1],
               span = 0.1, alpha = 0.05,
               n.grid = 1000, cut.off = 0.5)
    c(fit$main_res$AUC.empirical,as.numeric(fit$calibration_res[1]))
  })
  
  rowMeans(twROC.est,na.rm = TRUE)
})

cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])

## fitting of two-stage model with no direct effect
cp_time<-proc.time()
cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])
#save(adni_blab_fobj_ts,file = paste(res_dir,"adni_blab_fobj_ts_",seed_n,".RData",sep = ""))

#load(paste("~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/adni_da_results/adni_dp_results/adni_blab_fobj_ts_",seed_n+100,".RData",sep=""))

cp_time<-proc.time()

coxOBJ1<-coxph_mle_fllik(surv_time=tr_surv_tdat$surv_time,censor_status=tr_surv_tdat$censor_status,Xmat=sapply(1:ncol(adni_blab_fobj_ts$A.hat),function(k){adni_blab_fobj_ts$A.hat[,k]/sqrt(adni_blab_fobj_ts$Sigma2R[k])}),max_iter=100,cnv_type = "LL",cnv_crit=1e-6,constant_hazard = FALSE)
save(coxOBJ1,file = paste(res_dir,"coxOBJ_ts1_",seed_n,".RData",sep = ""))


auc_bs_ts1<-sapply(lm_time, function(u){
  ts_omics_dat<-lapply(1:ts_n_sub, function(i){
    lm_ind<-which(ts_ls_obs_dat[[i]]<=u)
    xx<-t(ts_ls_omics_dat[[i]][lm_ind,])
    yy<-ifelse(xx==0,0.01,xx)
    if(length(lm_ind)==1){
      rbind(ts_ls_obs_dat[[i]][lm_ind],as.matrix(log(yy)))
    } else{
      rbind(ts_ls_obs_dat[[i]][lm_ind],log(yy))
    }
  })
  
  tprd_msF<-predict.supFTSVD(obj = adni_blab_fobj_ts,designM = ts_res_vec,new_dat = ts_omics_dat,TimeG = time_grid)
  # two-stage conditional probability
  bz_fit<-coxOBJ1$CoxPH_obj$cum_hazard
  bz_grid<-approx(x=bz_fit$time,y=bz_fit$hazard,xout=time_grid,rule=2)$y
  tsg_par<-as.numeric(coxOBJ1$CoxPH_obj$coefficients)
  tpred_survP<-sapply(exp(sapply(1:ncol(adni_blab_fobj_ts$A.hat),function(k){tprd_msF$subEFF[,k]/sqrt(adni_blab_fobj_ts$Sigma2R[k])})%*%matrix(tsg_par)),function(u){exp(-u*bz_grid)})
  
  #del_time<-seq(max(lm_time),max(lm_time)+12,length.out=13)
  del_time<-seq(u,u+12,by=0.25)
  del_pos<-sapply(del_time, function(u){which.min((u-time_grid)^2)})
  
  test_cond_prob<-sapply(1:vN, function(i){
    a<-tpred_survP[del_pos,i]/tpred_survP[del_pos[1],i]
    1-a[-1]
  })
  
  
  twROC.est <- sapply(seq_len(nrow(test_cond_prob)),function(i){
    fit<-tdROC(X = test_cond_prob[i,], Y = ts_surv_tdat$surv_time,
               delta = 1-ts_surv_tdat$censor_status, tau = del_time[i+1],
               span = 0.1, alpha = 0.05,
               n.grid = 1000, cut.off = 0.5)
    c(fit$main_res$AUC.empirical,as.numeric(fit$calibration_res[1]))
  })
  
  rowMeans(twROC.est,na.rm = TRUE)
})

cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])


##
## fitting basic cox-ph model
cp_time<-proc.time()
cx_datF<-data.frame(tr_surv_pred)
cx_datF$surv_time<-tr_surv_tdat$surv_time
cx_datF$event_status<-(1-tr_surv_tdat$censor_status)
adni_coxph<-coxph(Surv(surv_time,event_status)~AGE+PTGENDER+PTEDUCAT+APOE4,data=cx_datF)
adni_cxbh<-basehaz(adni_coxph)
cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])
save(adni_coxph,file = paste(res_dir,"adni_coxph_",seed_n,".RData",sep = ""))
#load(paste("~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/adni_da_results/adni_dp_results/adni_blab_fobj_ts_",seed_n+100,".RData",sep=""))
cp_time<-proc.time()

#coxOBJ<-coxph_mle_fllik(surv_time=tr_surv_tdat$surv_time,censor_status=tr_surv_tdat$censor_status,Xmat=adni_blab_fobj_ts$A.hat,max_iter=100,cnv_crit=1e-6)

auc_bs_cx<-sapply(lm_time, function(u){
  # two-stage conditional probability
  bz_fit<-adni_cxbh
  bz_grid<-approx(x=bz_fit$time,y=bz_fit$hazard,xout=time_grid,rule=2)$y
  tsg_par<-as.numeric(coef(adni_coxph))
  tpred_survP<-sapply(exp(ts_surv_pred%*%matrix(tsg_par)),function(u){exp(-u*bz_grid)})
  
  del_time<-seq(u,u+12,by=0.25)
  del_pos<-sapply(del_time, function(u){which.min((u-time_grid)^2)})
  
  test_cond_prob<-sapply(1:vN, function(i){
    a<-tpred_survP[del_pos,i]/tpred_survP[del_pos[1],i]
    1-a[-1]
  })
  
  
  twROC.est <- sapply(seq_len(nrow(test_cond_prob)),function(i){
    fit<-tdROC(X = test_cond_prob[i,], Y = ts_surv_tdat$surv_time,
               delta = 1-ts_surv_tdat$censor_status, tau = del_time[i+1],
               span = 0.1, alpha = 0.05,
               n.grid = 1000, cut.off = 0.5)
    c(fit$main_res$AUC.empirical,as.numeric(fit$calibration_res[1]))
  })
  
  rowMeans(twROC.est,na.rm = TRUE)
})

cpTime<-c(cpTime,as.numeric(proc.time()-cp_time)[3])


data.frame(rbind(t(auc_bs_jm),t(auc_bs_jm1),t(auc_bs_ts),t(auc_bs_ts1),t(auc_bs_cx))) %>%
  set_names(c("iAUC","iBS")) %>%
  mutate("LmT"=rep(c(lm_time),5),
         "Method" = rep(c("JM","JM_S","TS","TS_S","CPH"),each=length(lm_time))) %>%
  write.csv(file = paste(res_dir,"auc_bs_res_",seed_n,".csv",sep = ))


data.frame(cpTime) %>%
  set_names("Time") %>%
  mutate("Type"=rep(c("MF","DP"),5),
         "Method" = rep(c("JM","JM_S","TS","TS_S","CPH"),each=2)) %>%
  write.table(file = paste(res_dir,"comp_time_",seed_n,".csv",sep = ))

