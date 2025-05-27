# setwd('')
library(tidyverse)
n <- 300 ## number of repeats of simulation
raw.dir = "/hpc/home/ma521/hdld_surv_jm/"
code.dir = "/cwork/ma521/hdld_surv_jm/RCodes/"

indir <- code.dir 
for (i in 1:n){
  fnm<-ifelse(i%/%100==0,"iter_hdld_surv_jm_100.R",
              ifelse(i%/%100==1,"iter_hdld_surv_jm_200.R","iter_hdld_surv_jm_300.R"))
  if(i==100){
    fnm<-"iter_hdld_surv_jm_100.R"
  } 
  
  if(i==200){
    fnm<-"iter_hdld_surv_jm_200.R"
  } 
  
  if(i==300){
    fnm<-"iter_hdld_surv_jm_300.R"
  } 
  
   # if(i==400){
   #   fnm<-"iter_hdld_surv_jm_500.R"
   # }
  
  source_code <- paste0("source('", raw.dir, fnm,"')")
  source_code <- as.data.frame(source_code) %>% mutate(V1=source_code) %>% dplyr::select(V1)
  seed_code <- paste0("seed_n <- ", i)
  seed_code <- as.data.frame(seed_code) %>% mutate(V1=seed_code) %>% dplyr::select(V1)
  code_file <- paste0(code.dir, i, '.R')
  R_code <- rbind(seed_code, source_code)
  write_delim(R_code, code_file, col_names = F, delim = '')
}
