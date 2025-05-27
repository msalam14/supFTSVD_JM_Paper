# necessary packages
library(tidyverse)
library(dplyr)
library(rsvd)
library(matrixcalc)
library(survival)
library(MASS)
# directory for for codes
#code_dir<-"~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/"
code_dir<-"/hpc/home/ma521/hdld_surv_jm_adni/"
source(paste(code_dir,"high_dimensional_long_surv_jm.R",sep=""))
source(paste(code_dir,"supFTSVD.R",sep = ""))

# data directory
#data_dir<-"~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/adni_data_files_hdld_surv_jm/"
#res_dir<-"~/Academic/PostDocWorks/Research_PostDocs/High_Dimensional_Long_Surv_JM/adni_da_results/"

data_dir<-"/hpc/home/ma521/hdld_surv_jm_adni/adni_data_files_processed/"
res_dir<-"/cwork/ma521/adni_da_hdld_jm/"


# loading full survial and bile acid data
file_name<-"lipidBL_surv_full_17_797.csv"
# surv_long_data<-read.csv(paste(data_dir,file_name,sep="")) %>%
#   dplyr::select_if(is.numeric) %>%
#     group_by(RID) %>%
#       mutate(n_obs=n()) %>%
#         filter(n_obs>1) %>%
#           ungroup()
study_period<-c(0,120)
surv_long_data<-read.csv(paste(data_dir,file_name,sep="")) %>%
  dplyr::select_if(is.numeric) %>%
  filter(visit_time<=max(study_period)) %>%
  group_by(RID) %>%
  mutate(n_obs=n()) %>%
  filter(n_obs>1) %>%
  mutate(censor_status=ifelse(surv_time>max(visit_time),1,censor_status)) %>%
  mutate(surv_time=max(visit_time)) %>%
  ungroup()


#tst_id<-sort(sample(unique(surv_long_data$RID),100,replace = FALSE))

#surv_long_data<- surv_long_data %>%
#  filter(RID %in% tst_id)

# construction of meta data (first 16 columns are survival and clinical variables)
meta_data<-surv_long_data %>%
  dplyr::select(1:16) %>%
  distinct() 

# design matrix for supervised decomposition  
res_vec<-meta_data  %>%
  mutate(Interp=1) %>%
  dplyr::select(RID,Interp,AGE,PTGENDER,PTEDUCAT,APOE4) %>%
  distinct() %>%
  dplyr::select(-RID) %>%
  as.matrix()

# survival covariates
surv_pred<-meta_data  %>%
  dplyr::select(RID,AGE,PTGENDER,PTEDUCAT,APOE4) %>%
  distinct() %>%
  dplyr::select(-RID) %>%
  as.matrix()


# survival data 
surv_tdat<-meta_data %>%
  dplyr::select(RID,surv_time,censor_status) %>%
  distinct() 

# omics data preparation
omics_datF<-surv_long_data %>%
  dplyr::select(-c(1:16))

ls_omics_dat<-split.data.frame(omics_datF,surv_long_data$RID)
ls_obs_dat<-split(surv_long_data$visit_time,surv_long_data$RID)

n_sub<-length(ls_obs_dat)

omics_dat<-lapply(1:n_sub, function(i){
  xx<-t(ls_omics_dat[[i]])
  yy<-ifelse(xx==0,0.01,xx)
  rbind(ls_obs_dat[[i]],log(yy))
})


# cross-validation of rank selection
cp_time<-proc.time()
cv_res_lipid_blab<-cv_rank_supFTSVD_JM(datlist = omics_dat,
                              response=res_vec,
                              interval = study_period,
                              ranks = 1:10,
                              resolution=50,
                              CVPhi=TRUE,
                              K=5,
                              cvT=5,
                              smooth=round(exp(seq(-5,2,length.out=20)),3),
                              maxiter=200,
                              epsilon=1e-6,
                              KInd=NULL,
                              surv_time = surv_tdat$surv_time,
                              censor_status = surv_tdat$censor_status,
                              survX = surv_pred,
                              mc_ss = 1000,
                              constant_hazard=FALSE,
                              rsvd_seed=100,
                              conv_criteria = "cond_like",
                              scale=TRUE,
                              change_points=NULL,
                              rank_fold = 5,
                              rfolds_seed = 20,
                              rc_thresh = 0.05,
                              stratum = surv_tdat$censor_status,
                              par_computing=TRUE,
                              n_core=5)

save(cv_res_lipid_blab,file = paste(res_dir,"cv_res_lipid_blab_",max(study_period),".RData",sep=""))

load(file = paste(res_dir,"cv_res_lipid_blab_",max(study_period),".RData",sep=""))

# fitting of the joint model
jm_sup_res_lipid_blab<-supFTSVD_JM(datlist = omics_dat, 
                                  response=res_vec, 
                                  interval = study_period, 
                                  r = cv_res_lipid_blab$opt_r,
                                  resolution=50, 
                                  CVPhi=TRUE, 
                                  K=5, 
                                  cvT=5, 
                                  smooth=round(exp(seq(-5,2,length.out=20)),3),
                                  maxiter=200, 
                                  epsilon=1e-6,
                                  KInd=NULL,
                                  rsvd_seed=100,
                                  conv_criteria = "cond_like",
                                  surv_time = surv_tdat$surv_time,
                                  censor_status = surv_tdat$censor_status,
                                  survX = surv_pred,
                                  constant_hazard=FALSE)

save(jm_sup_res_lipid_blab,file = paste(res_dir,"jm_sup_res_lipid_blab_",max(study_period),".RData",sep=""))
proc.time()-cp_time

