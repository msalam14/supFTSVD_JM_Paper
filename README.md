# supFTSVD-JM: a joint modeling approach


# Introduction

supFTSVD-JM is a joint modeling framework for high-dimensional
longitudinal and time-to-event data with right censoring. It integrates
a low-rank tensor decomposition approach, namely, supFTSVD, with Cox
proportional hazard model. A maximum likelihood estimation is developed
through a penalized Monte Carlo Expectation-Maximization (MCEM)
algorithm. The associated *R* package is available at . This vignette
illustrates implementation of supFTSVD-JM on a simulated data.

- Call necessary R packages

``` r
library(splines)
library(rstan)
library(parallel)
library(foreach)
library(tdROC)
library(survival)
library(supFTSVDJM)
library(supFTSVD)
```

- For installing the *supFTSVDJM* package, the following codes can be
  used

``` r
devtools::install_github("https://github.com/msalam14/supFTSVD_JM")
```

# Simulation study

- We carry out the simulation study using the Duke Cluster
  Computing (DCC) facilities. The 100 replicates of our Monte Carlo
  study are submitted as batch jobs.

- Codes used in the simulation study are available in the directory
  *SimulationStudy_RCodes*. Two settings investigated are named as
  *Censoring_25* and *Censoring_40*, representing censoring scenarios
  $25\%$ and $40\%$, respectively.

- There are three scripts in each censoring directory, namely,
  *iter_hdld_surv_jm_100.R*, *iter_hdld_surv_jm_200.R*, and
  *iter_hdld_surv_jm_300.R* for three sample sizes $100$, $200$, and
  $300$, respectively.

- For parallel running, we create three hundred scripts from these
  scripts, 100 copies for each of them. In each censoring directory, the
  script *create_codes_hdld_jm.R* accomplishes this task.

- For batch job submission, we used the file named as *sim_jm_hdld.sh*.

- Note that, files *high_dimensional_long_surv_jm.R* and *supFTSVD.R*
  contain source codes *R* packages *supFTSVDJM* and *supFTSVD*,
  respectively. One can either source these two files or call the
  corresponding *R* packages.

- For illustration, we obtain results for one Monte Carlo sample of size
  50 with feature dimension $100$.

``` r
seed_n<-25
source("SimulationStudy_RCodes/iter_hdld_surv_jm_50.R")
```

- Each run save three files, *hdld_tr_res_n_seed_n.txt*,
  *hdld_prd_res_n_seed_n.txt*, and *hdld_cmp_time_n_seed_n.txt*, where
  $n$ and $seed_n$ will be replaced by the used sample size and
  seed_n. For example in our demonstration the file names are
  *hdld_tr_res_50_25.txt*, *hdld_prd_res_50_25.txt*, and
  *hdld_cmp_time_50_25.txt*.

- For saving estimation results, we used *hdld_tr_res_50_25.txt*, which
  has following results

``` r
tr_res<-read.table("SimulationStudy_Results/hdld_tr_res_50_25.txt")
tr_res
```

                  JM           Two     Comps      Stat
    1   7.694925e-04  8.397033e-03       k=1  SingFunc
    2   1.826231e+00  1.129564e+00       k=2  SingFunc
    3   1.859482e+00  1.939306e+00       k=3  SingFunc
    4   5.600643e-03  5.168577e-02       k=1  FeatLoad
    5   4.972057e+00  3.456814e+00       k=2  FeatLoad
    6   5.163944e+00  5.280614e+00       k=3  FeatLoad
    7   0.000000e+00  1.417769e-02       all     difR2
    8   9.431084e+00  1.064758e+01    survLL       sLL
    9   2.446234e+04  2.184387e+04   jointLL   jointLL
    10  3.890144e+00  3.169484e+00 surv_coef surv_coef
    11  3.298456e+00  2.435853e+00 surv_coef surv_coef
    12 -1.853669e-01 -2.466834e-01 surv_coef surv_coef
    13  6.653620e-02  2.339719e-01 surv_coef surv_coef
    14 -1.285900e-01  1.699916e-01 surv_coef surv_coef

- On the other hand, for saving prediction results, we used
  *hdld_prd_res_50_25.txt*, which has following results

``` r
prd_res<-read.table("SimulationStudy_Results/hdld_prd_res_50_25.txt")
prd_res
```

       LM   prp_auc    tw_auc    tr_auc    prp_bs     tw_bs     tr_bs     prp_err
    1 0.2 0.8436173 0.8312768 0.8877225 0.1647860 0.1647806 0.1803827 0.115128479
    2 0.3 0.7641416 0.6928202 0.7903771 0.1869062 0.1776233 0.1917367 0.022728664
    3 0.4 0.7697522 0.7368826 0.7897728 0.1448116 0.1404785 0.1518835 0.009914312
    4 0.5 0.7708601 0.6622223 0.8269552 0.1385599 0.1380585 0.1428449 0.007090410
        sprp_err     prp_scrv     ts_scrv
    1 0.50207630 0.0069935049 0.008043215
    2 0.15625021 0.0016558160 0.006389419
    3 0.07158648 0.0012116711 0.007004464
    4 0.04670663 0.0008544647 0.006977670

- All original results obtained are available in the
  *SimulationStudy_Results*

# ADNI data analysis:

Due to unavailability of permission, we cannot share the data in this
repository. Thus, demonstration on the data is skipped here. However,
codes used and results obtained for ADNI data are available in the
directories *ADNI_RCodes* and *ADNI_Results*, respectively.
