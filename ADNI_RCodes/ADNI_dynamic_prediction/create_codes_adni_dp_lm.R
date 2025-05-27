# setwd('')
library(tidyverse)
n <- 120 ## number of repeats of simulation
raw.dir = "/hpc/home/ma521/hdld_surv_jm_adni/"
code.dir = "/cwork/ma521/adni_da_hdld_jm/RCodes/"

for (i in 1:n){
  fnm<-"adni_dp_lip_blab_lm.R"
  source_code <- paste0("source('", raw.dir, fnm,"')")
  source_code <- as.data.frame(source_code) %>% mutate(V1=source_code) %>% dplyr::select(V1)
  seed_code <- paste0("seed_n <- ", i+100)
  seed_code <- as.data.frame(seed_code) %>% mutate(V1=seed_code) %>% dplyr::select(V1)
  lmT<-paste0("lm_time<-",ifelse(i<=30,12,ifelse(i<=60,15,ifelse(i<=90,18,24))))
  lm_code<-as.data.frame(lmT) %>% mutate(V1=lmT) %>% dplyr::select(V1)
  code_file <- paste0(code.dir, i, '.R')
  R_code <- rbind(lmT,seed_code, source_code)
  write_delim(R_code, code_file, col_names = F, delim = '')
}
