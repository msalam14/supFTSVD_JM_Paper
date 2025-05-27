#' Supervised functional singular value decomposition for high-dimensional 
#' multivariate functional data.
#' 
#' This function performs low-rank decomposition of high-dimensional 
#' multivariate function following the methodology proposed in Alam (2024+)
#' 
#' @param datlist a list with n elements, each is a matrix of dimension p+1 
#' times m_i. The first row represents time points where measurements were 
#' obtained. Next p rows represent the function values observed at the time 
#' points
#' @param response a matrix of dimension n by q contains the design matrix for 
#' supervision of the subject loadings
#' @param interval represents the domain of the function. The default is NULL. 
#' @param r number of low-rank components involved in the assumed model
#' @param resolution grid size on which the singular functions will be estimated
#' @param CVPhi a logical scalar representing whether cross-validation (CV) will be 
#' performed for determining the smoothness parameters involved with singular 
#' functions
#' @param K number of folds for the CV step
#' @param cvT number of initial iterations at which CV will be performed
#' @param smooth Smoothing parameter for RKHS norm when CVPhi is FALSE. With 
#' CVPhi=TRUE, a vector of numeric values for smoothness parameters for grid
#' search
#' @param maxiter Maximum number of iteration. Default: 20.
#' @param epsilon Convergence threshold for objective function value
#' @param KInd is a vector indices to construct the folds for CV. The default is
#' NULL that indicates that folds will be constructed randomly inside
#' @importFrom stats approx coef cor lm rnorm
#' @return a list with following objects:
#' \itemize{
#'  \item A.hat : a matrix of dimension n by r containing the values of 
#'  conditional means for subjects
#'  \item B.hat : a matrix of dimension p by r representing the subject loadings
#'  \item Phi.hat : a matrix of dimension resolution by r representing the 
#'  estimated singular functions
#'  \item time : a vector of time points where Phi.hat values are estimated
#'  \item Lambda : a vector of r lambda values for weights (CP decomposition) of
#'   low-rank components
#'  \item r.square : component wise R^2 for tensor model by every rank-1 
#'  component constructed using conditional mean
#'  \item accum.r.square : cumulative R^2 values of tensor model by estimated 
#'  components constructed using conditional means
#'  \item r.square.Xb : component wise R^2 for tensor model by every rank-1 
#'  component constructed using the deterministic part of conditional mean
#'  \item accum.r.square.Xb : cumulative R^2 values of tensor model by estimated
#'   components constructed using deterministic parts of conditional means
#'   \item supR2 : total R^2 achieved by all estimated rank-1 components
#'   \item supR2.Xb : total R^2 achieved by all estimated rank-1 components 
#'   involving the deterministic parts of subject loadings
#'   \item Gamma : a matrix of q by r representing the estimated beta parameters
#'   \item Xb=Xb : deterministic parts of the conditional means
#'   \item Sigma2R : r estimated variance parameters involving with subject 
#'   loading models
#'   \item Sigma2 : estimated tensor noise variance
#'   \item cnvrgnt : change of objective function value at different iterations
#' }
#' @example examples/example_supFTSVD.R
#' @export
supFTSVD <- function(datlist, response, interval = NULL, r = 3,resolution=100, CVPhi=FALSE, K=5, cvT=5, smooth=1e-8,
                     maxiter=20, epsilon=1e-4,KInd=NULL,rsvd_seed=100,conv_criteria="par"){
  n = length(datlist)
  p = nrow(datlist[[1]])-1
  
  # response in matrix form
  if(is.null(response)){
    response<-matrix(rep(1,n),ncol=1)
  } else if(!is.matrix(response)){
    response<-matrix(response,ncol=1)
  }
  
  # storing original response
  ires<-response
  
  
  # Calculate range.
  timestamps.all = do.call(c,lapply(datlist, FUN=function(u){u[1,]}))
  
  timestamps.all = sort(unique(timestamps.all))
  if (is.null(interval)){
    interval = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  }
  
  # rescale the time to 0-1.
  input.time.range = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  for (i in 1:n){
    datlist[[i]][1,] = (datlist[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  }
  interval = (interval - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  
  ti = lapply(1:n,FUN=function(u){
    temp = 1 + round((resolution-1) * (datlist[[u]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] = 0
    temp
  })
  
  y0=do.call(c,lapply(1:n,FUN=function(u){
    as.vector(t(datlist[[u]][-1,which(ti[[u]]>0)]))
  }))
  
  y<-y0
  
  tipos<-lapply(ti,function(u){u>0})
  
  Lt<-lapply(datlist,FUN = function(u){u[1,]})
  ind_vec<-rep(1:n,sapply(Lt,length))
  
  tm <- unlist(Lt)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)
  
  # Initial values
  ## initialization of b
  data.unfold = do.call(cbind,lapply(datlist,FUN=function(u){u[-1,]}))
  
  set.seed(rsvd_seed)
  b.initials <- rsvd(data.unfold, k=r)$u
  b.hat = as.matrix(b.initials[,1:r])
  
  ## initialization of a
  if(r==1){
    a.hat<-as.matrix(sapply(1:n, function(i){
      colMeans(t(datlist[[i]][-1,])%*%b.hat)
    }))
  } else{
    a.hat<-t(sapply(1:n, function(i){
      colMeans(t(datlist[[i]][-1,])%*%b.hat)
    }))
  }
  
  
  ## initialize U
  Ahat<-a.hat
  
  # Initialize xi function
  # compress the data and apply function PCA.
  phi.hat<-sapply(1:r, function(k){
    Ly<-lapply(1:n,FUN = function(i){
      Ahat[i,k]*as.numeric(b.hat[,k]%*%datlist[[i]][2:(p+1),])
    })
    
    if(CVPhi){
      phiH = cv_freg_rkhs(Ly, Ahat[,k], ind_vec, Kmat, Kmat_output, smooth=smooth, kfold = K,KInd=KInd)[[2]]
    } else{
      phiH = freg_rkhs(Ly, Ahat[,k], ind_vec, Kmat, Kmat_output, smooth=smooth)
    }
    #phiH = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
    (phiH / sqrt(sum(phiH^2)))*sqrt(resolution)
  })
  
  t_arg<-seq(interval[1],interval[2],length.out = resolution)
  
  allT<-do.call(c,Lt)
  pXI<-apply(phi.hat,2,function(u){
    approx(t_arg,u,xout = allT,rule=2)$y
  })
  
  
  # initial values for error variance
  mi<-as.numeric(sapply(Lt,length))
  Mi<-cumsum(mi)
  
  prdXI<-lapply(1:n,function(w){
    if(w ==1){
      as.matrix(pXI[1:Mi[1],])
    } else{
      as.matrix(pXI[(Mi[w-1]+1):(Mi[w]),])
    }
  })
  
  cwY2<-lapply(1:r, function(k){
    do.call(cbind,lapply(1:n, function(i){
      outer(b.hat[,k],prdXI[[i]][,k])*Ahat[i,k]
    }))
  })
  
  prdY2<-Reduce(`+`,cwY2)
  
  sig<-mean(as.numeric((data.unfold-prdY2)^2))
  
  # Initial values for sigma^2_k
  # lambda adjusted a_i
  LAhat<-sapply(1:r,function(k){
    (as.numeric(coef(lm(as.numeric(data.unfold)~as.numeric(cwY2[[k]])-1))))*Ahat[,k]
  })
  
  Ahat<-LAhat
  
  
  ## Initialize for gamma
  Resp<-response #cbind(1,response)
  gammaP<-sapply(1:r,function(k){
    as.numeric((solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(Ahat[,k],ncol = 1)))
  })
  
  if(ncol(Resp)==1){
    gammaP<-t(as.matrix(gammaP))
  }
  
  
  # for all subject, rows correspond to subjects
  Xb<-Resp%*%gammaP
  
  
  ## initialize for sigma_r
  
  sigR<-colMeans((Ahat-Xb)^2)
  
  sigma_R<-diag(r)
  diag(sigma_R)<-1/sigR
  
  # Necessary matrices
  Hmat<-lapply(1:n, function(i){
    sapply(1:r, function(k){
      as.numeric(outer(b.hat[,k],prdXI[[i]][,k]))
    })
  })
  
  # for all subject, rows correspond to subjects
  HiB<-Resp%*%gammaP
  
  
  ## Subject wise super X matrix
  Hbeta<-lapply(1:n, function(i){
    as.numeric(Hmat[[i]]%*%matrix(HiB[i,],ncol=1))
  })
  
  # Subject wise vectorized data 
  Yvec<-lapply(1:n, function(i){
    as.numeric(datlist[[i]][-1,])
  })
  
  # Subject-specific conditional covariances
  sw_vcov<-lapply(Hmat,function(u){
    solve(crossprod(u)*(1/sig)+sigma_R)
  })
  
  
  if(r==1){
    varU<-as.matrix(sapply(1:n, function(i){
      diag(sw_vcov[[i]])
    }))
  } else{
    varU<-t(sapply(1:n, function(i){
      diag(sw_vcov[[i]])
    }))
  }
  
  # Prediction of conditional mean using initial values
  conU<-Ahat-Xb
  
  # another necessary matrix
  HUmat<-lapply(1:n, function(i){
    as.numeric(Hmat[[i]]%*%matrix(conU[i,],ncol=1))
  })
  
  
  
  ## objective function value using initial values
  obj_val<-Reduce(`+`,lapply(1:n, function(i){
    t1<-(-((p*mi[i])/2)*log(sig))-sum((1/2)*log(sigR))
    t2<-sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2/(sig))
    t3<-sum(diag((crossprod(Hmat[[i]])*(1/sig))%*%sw_vcov[[i]]))
    t4<-sum(diag(sigma_R%*%sw_vcov[[i]]))
    t5<-sum((conU[i,]^2/sigR))
    t1-(t2+t3+t4+t5)/2
  }))
  
  # Starting of the EM based estimation
  t<-1
  iter_dif<-NULL
  while(t<=(maxiter+cvT)){
    dif<-NULL
    # E step for EM algorithm
    # conditional mean for all subjects
    if(r==1){
      conU<-as.matrix(sapply(1:n,FUN = function(i){
        sw_vcov[[i]]%*%(t(Hmat[[i]]*(1/sig))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))
    } else{
      conU<-t(sapply(1:n,FUN = function(i){
        sw_vcov[[i]]%*%(t(Hmat[[i]]*(1/sig))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))
    }
    
    # M step
    # update of beta_k for all k
    eta_val<-NULL
    for(k in 1:r){
      #Xb<-Resp%*%gammaP
      Anew<-conU+Xb
      if(r>1){
        indx<-c(1:r)[-k]
      }
      
      # residual data
      if(r==1){
        rY<-lapply(1:n, function(i){
          datlist[[i]][-1,]
        })
      } else{
        rY<-lapply(1:n, function(i){
          datlist[[i]][-1,]-Reduce(`+`,lapply(indx, function(l){
            Anew[i,l]*outer(b.hat[,l],prdXI[[i]][,l])
          }))
        })
      }
      
      
      # design matrix for beta update
      Xs<-do.call(rbind,lapply(1:n, function(i){
        do.call(rbind,lapply(prdXI[[i]][,k], function(j){
          t(outer(Resp[i,],b.hat[,k])*j)
        }))
      }))
      
      ## Creating Z vector
      Zs<-matrix(do.call(c,lapply(1:n, function(i){
        as.numeric(rY[[i]]-(conU[i,k]*outer(b.hat[,k],prdXI[[i]][,k])))
      })))
      
      ngamma<-solve(t(Xs)%*%(Xs))%*%(t(Xs)%*%Zs)
      dif<-c(dif,sum((ngamma-gammaP[,k])^2)/sum(gammaP[,k]^2))
      gammaP[,k]<-ngamma
      
      # updated Xb
      Xb<-Resp%*%gammaP
      Anew<-conU+Xb
      
      # Scaling factor
      scD<-sqrt((Anew[,k]^2)+sapply(sw_vcov,function(u){u[k,k]}))
      
      if(r==1){
        srY<-lapply(1:n,function(i){
          (matrixcalc::hadamard.prod(outer(rep(Anew[i,k],p),rep(1,nrow(prdXI[[i]]))),rY[[i]]))*(1/scD[i])
        })
      } else{
        srY<-lapply(1:n,function(i){
          (matrixcalc::hadamard.prod(outer(rep(Anew[i,k],p),rep(1,nrow(prdXI[[i]]))),rY[[i]])-Reduce(`+`,lapply(indx,function(l){
            outer(b.hat[,l],prdXI[[i]][,l])*sw_vcov[[i]][k,l]
          })))*(1/scD[i])
        })
      }
      
      # update of xi_b k
      
      # scaled Residual data
      sY<-sapply(1:p,function(b){
        do.call(c,lapply(srY, function(u){u[b,]}))
      })
      
      
      nbhat<-NULL
      for(b in 1:p){
        bX<-do.call(c,lapply(1:n,function(i){scD[i]*prdXI[[i]][,k]}))
        
        nbhat<-c(nbhat,sum(sY[,b]*bX)/sum(bX^2))
      }
      n_bhat<-nbhat/sqrt(sum(nbhat^2))
      dif<-c(dif,sum((n_bhat-b.hat[,k])^2))
      b.hat[,k]<-n_bhat
      
      # Update of psi function
      Lty = lapply(1:n,function(i){
        list(scD[i]*as.numeric(b.hat[,k]%*%srY[[i]]))
        #(matrixcalc::hadamard.prod(scF[[i]],rY[[i]]))))
      })
      
      if(CVPhi & t<=cvT){
        cvfit<-cv_freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=smooth,kfold = K,KInd=KInd)
        Smv<-cvfit[[1]]
        phiH = cvfit[[2]]
        eta_val<-c(eta_val,Smv)
      } else{
        phiH = freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=ifelse(CVPhi,Smv,smooth))
        eta_val<-c(eta_val,ifelse(CVPhi,Smv,smooth))
      }
      phi.new = (phiH / sqrt(sum(phiH^2)))*sqrt(resolution)
      dif<-c(dif,sum((phi.hat[,k] - phi.new)^2))
      phi.hat[,k] <- phi.new
      
      kpXI<-approx(t_arg,phi.new,xout = allT,rule=2)$y
      
      for(i in 1:n){
        if(i ==1){
          prdXI[[i]][,k]<-kpXI[1:Mi[1]]
        } else{
          prdXI[[i]][,k]<-kpXI[(Mi[i-1]+1):(Mi[i])]
        }
      }
    }
    
    # update of sigma2_k
    sk2_new<-colMeans(conU^2+varU)
    
    # to avoid sigma2k=0
    if(any(sk2_new<1e-10)){
      sk2_new[sk2_new<1e-10]<-1e-10
    }
    
    
    dif<-c(dif,((sk2_new-sigR)^2)/sigR^2)
    
    
    sigR<-sk2_new
    
    diag(sigma_R)<-1/sigR
    
    
    # update of sigma2
    # Necessary matrices (updated)
    Hmat<-lapply(1:n, function(i){
      sapply(1:r, function(k){
        as.numeric(outer(b.hat[,k],prdXI[[i]][,k]))
      })
    })
    
    # for all subject, rows correspond to subjects
    HiB<-Resp%*%gammaP
    
    
    ## Subject wise super X matrix
    Hbeta<-lapply(1:n, function(i){
      as.numeric(Hmat[[i]]%*%matrix(HiB[i,],ncol=1))
    })
    
    HUmat<-lapply(1:n, function(i){
      as.numeric(Hmat[[i]]%*%matrix(conU[i,],ncol=1))
    })
    
    nsig2<-sum(sapply(1:n, function(i){
      sum(diag((t(Hmat[[i]])%*%Hmat[[i]])%*%sw_vcov[[i]]))+
        sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2)
    }))/sum(sapply(Yvec,length))
    
    dif<-c(dif,((nsig2-sig)^2)/sig^2)
    sig<-nsig2
    
    # objective function
    nobj_val<-Reduce(`+`,lapply(1:n, function(i){
      t1<-(-((p*mi[i])/2)*log(sig))-sum((1/2)*log(sigR))
      t2<-sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2/(sig))
      t3<-sum(diag((crossprod(Hmat[[i]])*(1/sig))%*%sw_vcov[[i]]))
      t4<-sum(diag(sigma_R%*%sw_vcov[[i]]))
      t5<-sum((conU[i,]^2/sigR))
      t1-(t2+t3+t4+t5)/2
    }))+sum(eta_val)
    
    
    par_dif<-max(dif)
    dif<-c(dif,nobj_val)
    iter_dif<-rbind(iter_dif,dif)
    
    if(CVPhi & t<=cvT){
      relC<-epsilon+0.05
    } else{
      relC<-(nobj_val-obj_val)/abs(obj_val)
    }
    
    if(conv_criteria=="par"){
      c_crit<-par_dif
    } else{
      c_crit<-relC
    }
    
    if(t>cvT & c_crit<=epsilon) 
      break
    
    obj_val<-nobj_val
    
    t <- t+1
    
    # Subject-specific conditional covariances
    sw_vcov<-lapply(Hmat,function(u){
      solve(crossprod(u)*(1/sig)+sigma_R)
    })
    
    if(r==1){
      varU<-as.matrix(sapply(1:n, function(i){
        diag(sw_vcov[[i]])
      }))
    } else{
      varU<-t(sapply(1:n, function(i){
        diag(sw_vcov[[i]])
      }))
    }
    
  }
  
  # calculate lambda
  tenY<-do.call(c,Yvec)
  
  # using Ahat
  cX<-do.call(rbind,lapply(1:n,function(i){
    sapply(1:r, function(k){
      as.numeric(Anew[i,k]*outer(b.hat[,k],prdXI[[i]][,k]))
    })
  }))
  
  cl.fit = lm(tenY~cX-1)
  cLambda = as.numeric(cl.fit$coefficients)
  cR2<-sapply(1:r, function(k){
    summary(lm(tenY~cX[,k]-1))$r.squared
  })
  
  cRsq<-sapply(1:r, function(k){
    summary(lm(tenY~cX[,1:k]-1))$r.squared
  })
  
  # using Xb
  X<-do.call(rbind,lapply(1:n,function(i){
    sapply(1:r, function(k){
      as.numeric(Xb[i,k]*outer(b.hat[,k],prdXI[[i]][,k]))
    })
  }))
  
  
  
  l.fit = lm(tenY~X-1)
  Lambda = as.numeric(l.fit$coefficients)
  mRsq=sapply(1:r, function(k){
    summary(lm(tenY~X[,k]-1))$r.squared
  })
  
  Rsq<-sapply(1:r, function(k){
    summary(lm(tenY~X[,1:k]-1))$r.squared
  })
  
  
  # elements of returning object
  PCname <- paste('Component', 1:r)
  colnames(Anew) = PCname
  colnames(Xb) = PCname
  colnames(b.hat) = PCname
  colnames(phi.hat) = PCname
  rownames(Anew) = names(datlist)
  rownames(Xb) = names(datlist)
  rownames(b.hat) = rownames(datlist[[1]])[-1]
  
  ## Time point where phi functions are estiamted
  
  time.return = seq(interval[1],interval[2],length.out = resolution)
  time.return = time.return * (input.time.range[2] - input.time.range[1]) + input.time.range[1]
  
  # return object
  results = list("A.hat" = Anew, "B.hat" = b.hat, 
                 "Phi.hat" = phi.hat, "time" = time.return,
                 "Lambda" = cLambda, "r.square" = cR2, "accum.r.square" = cRsq, 
                 "r.square.Xb" = mRsq, "accum.r.square.Xb" = Rsq,"supR2"=cRsq[r],"supR2.Xb"=Rsq[r],
                 "Gamma"=gammaP,"Xb"=Xb,"Sigma2R"=sigR,"Sigma2"=sig,cnvrgnt=iter_dif)
  return(results)
}

#' Functional tensor singular value decomposition (FTSVD)
#' 
#' This function models high-dimensional multivariate functional data using the
#' FTSVD approach (Han et. al. 2023). 
#' 
#' @param datlist a list with n elements, each is a matrix of dimension p+1 
#' times m_i. The first row represents time points where measurements were 
#' obtained. Next p rows represent the function values observed at the time 
#' points
#' @param interval represents the domain of the function. The default is NULL. 
#' @param r number of low-rank components involved in the assumed model
#' @param resolution grid size on which the singular functions will be estimated
#' @param CVPhi a logical scalar representing whether cross-validation (CV) will be 
#' performed for determining the smoothness parameters involved with singular 
#' functions
#' @param K number of folds for the CV step
#' @param cvT number of initial iterations at which CV will be performed
#' @param smooth Smoothing parameter for RKHS norm when CVPhi is FALSE. With 
#' CVPhi=TRUE, a vector of numeric values for smoothness parameters for grid
#' search
#' @param maxiter Maximum number of iteration. Default: 20.
#' @param epsilon Convergence threshold for objective function value
#' @param KInd is a vector indices to construct the folds for CV. The default is
#' NULL that indicates that folds will be constructed randomly inside
#' @return a list with following objects:
#' \itemize{
#'  \item A.hat : a matrix of dimension n by r containing the values of 
#'  conditional means for subjects
#'  \item B.hat : a matrix of dimension p by r representing the subject loadings
#'  \item Phi.hat : a matrix of dimension resolution by r representing the 
#'  estimated singular functions
#'  \item time : a vector of time points where Phi.hat values are estimated
#'  \item Lambda : a vector of r lambda values for weights (CP decomposition) of
#'   low-rank components
#'  \item r.square : component wise R^2 for tensor model by every rank-1 
#'  component constructed using conditional mean
#'  \item accum.r.square : cumulative R^2 values of tensor model by estimated 
#'  components constructed using conditional means
#' }
#' @example examples/example_FTSVD.R
#' @export
ftsvd<-function(datlist, interval = NULL, r = 3, resolution = 251, CVPhi=FALSE, 
                K=5, cvT=5, smooth=1e-8,maxiter=20, epsilon=1e-4,KInd=NULL){
  n = length(datlist)
  p = nrow(datlist[[1]])-1
  
  Lambda = rep(0, r)
  A = matrix(0, n, r)
  B = matrix(0, p, r)
  Phi = matrix(0, resolution, r)
  PCname <- paste('Component', 1:r)
  colnames(A) = PCname
  colnames(B) = PCname
  colnames(Phi) = PCname
  rownames(A) = names(datlist)
  rownames(B) = rownames(datlist[[1]])[-1]
  
  
  # Calculate range.
  timestamps.all = NULL
  for (i in 1:n){
    timestamps.all = c(timestamps.all, datlist[[i]][1,])
  }
  
  timestamps.all = sort(unique(timestamps.all))
  if (is.null(interval)){
    interval = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  }
  
  # rescale the time to 0-1.
  input.time.range = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  for (i in 1:n){
    datlist[[i]][1,] = (datlist[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  }
  interval = (interval - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  
  res = NULL
  Lambda = rep(0, r)
  X = NULL
  y0 <- NULL
  Rsq <- accumRsq <- rep(0, r)
  
  ti <- vector(mode='list', length=n)
  for (i in 1:n){
    temp = 1 + round((resolution-1) * (datlist[[i]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] = 0
    ti[[i]] <- temp
  }
  
  tipos <- vector(mode='list', length=n)
  for (i in 1:n){
    keep <- ti[[i]]>0
    tipos[[i]] <- keep
    y0 <- c(y0, as.vector(t(datlist[[i]][2:(p+1),keep])))
  }
  
  Lt = list()
  ind_vec <- NULL
  for (i in 1:n){
    Lt = c(Lt, list(datlist[[i]][1,]))
    ind_vec <- c(ind_vec, rep(i,length(Lt[[i]])))
  }
  
  tm <- unlist(Lt)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)
  
  for (s in 1:r){ 
    # calculate rank-1 component sequentially.
    # Step 1: initialization.
    print(sprintf("Calculate the %dth Component", s))
    
    # intialization of b
    data.unfold = NULL
    y <- NULL
    for (i in 1:n){
      data.unfold = cbind(data.unfold, datlist[[i]][2:(p+1),])
      y <- c(y, as.vector(t(datlist[[i]][2:(p+1),tipos[[i]]])))
    }
    b.initials <- svd(data.unfold, nu=r, nv=r)$u
    b.hat = b.initials[,1]
    # initialization of a
    a.hat <- rep(1,n)/sqrt(n)
    
    # iteratively update a,b,phi
    t <- 0
    dif <- 1
    while(t<=(maxiter+cvT) & dif>epsilon){
      # update phi:
      Ly = list()
      for (i in 1:n){
        Ly = c(Ly, list(a.hat[i]*as.numeric(b.hat%*%datlist[[i]][2:(p+1),])))
      }
      
      if(CVPhi & t<=cvT){
        cvfit<-cv_freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth,kfold = K,KInd=KInd)
        Smv<-cvfit[[1]]
        phi.hat = cvfit[[2]]
      } else{
        phi.hat = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=ifelse(CVPhi,Smv,smooth))
      }
      
      #phi.hat = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
      phi.hat = (phi.hat / sqrt(sum(phi.hat^2)))*sqrt(resolution)
      
      # update a:
      a.tilde <- rep(0,n)
      for (i in 1:n){
        t.temp <- tipos[[i]]
        a.tilde[i] <- b.hat %*% datlist[[i]][2:(p+1),t.temp] %*% phi.hat[ti[[i]][t.temp]] 
        a.tilde[i] <- a.tilde[i] / sum((phi.hat[ti[[i]][t.temp]])^2)
      }
      a.new <- a.tilde / sqrt(sum(a.tilde^2))
      if(t>cvT)
        dif <- sum((a.hat - a.new)^2)
      a.hat <- a.new
      
      # update b:
      temp.num <- matrix(0,p,n)
      temp.denom <- rep(0,n)
      for (i in 1:n){
        t.temp <- tipos[[i]]
        temp.num[,i] <- datlist[[i]][2:(p+1),t.temp] %*% phi.hat[ti[[i]][t.temp]]
        temp.denom[i] <-sum((phi.hat[ti[[i]][t.temp]])^2)
      }
      b.tilde <- as.numeric(temp.num%*%a.hat) / as.numeric(temp.denom%*%(a.hat^2))
      b.new <- b.tilde / sqrt(sum(b.tilde^2))
      if(t>cvT)
        dif <- max(dif, sum((b.hat - b.new)^2))
      b.hat <- b.new
      
      t <- t+1
    }
    
    # calculate lambda
    x = NULL
    for (i in 1:n){
      t.temp = ti[[i]]
      t.temp <- t.temp[t.temp>0]
      x <- c(x,as.vector(t(a.hat[i]*b.hat%o%phi.hat[t.temp])))
    }
    X = cbind(X, x)
    l.fit = lm(y~x-1)
    lambda = as.numeric(l.fit$coefficients)
    A[,s] = a.hat
    #P.A = P.A - a.hat %*% t(a.hat)
    B[,s] = b.hat
    Phi[,s] = t(phi.hat)
    Lambda[s] = lambda
    Rsq[s] <- summary(l.fit)$r.squared
    accumRsq[s] <- summary(lm(y0~X-1))$r.squared
    
    # update datlist
    for (i in 1:n){
      temp <- tipos[[i]]
      datlist[[i]][2:(p+1),which(temp)] = datlist[[i]][2:(p+1),which(temp)] - 
        Lambda[s] * A[i,s] * (B[,s] %*% t(Phi[ti[[i]][temp],s])) 
    }
    print(paste0("Convergence reached at dif=", dif, ', iter=', t-cvT))
  }
  l.fit = lm(y0~X-1)
  Lambda = as.numeric(l.fit$coefficients)
  
  # revise the sign of Lambda
  for (r in length(Lambda)){
    if (Lambda[r]<0){
      Lambda[r] = -Lambda[r]
      B[,r] = -B[,r]
    }
  }
  
  time.return = seq(interval[1],interval[2],length.out = resolution)
  time.return = time.return * (input.time.range[2] - input.time.range[1]) + input.time.range[1]
  results = list("A.hat" = A, "B.hat" = B, 
                 "Phi.hat" = Phi, "time" = time.return,
                 "Lambda" = Lambda, "r.square" = Rsq, "accum.r.square" = accumRsq)
  return(results)
}

#' Prediction using the fitted objects from supFTSVD
#' 
#' @param obj fitted object by supFTSVD
#' @param designM matrix of covariates for all subjects (rows represent subject)
#' @param new_dat is a list objects; each element is observed data for a subjects. 
#' The first row contains the time-points where measurements were taken. Thus, 
#' each element has (p+1) rows. The default is NULL; it is only required when 
#' subject specific prediction is the interest. 
#' @param TimeG Grid in the time domain where prediction will be made.
#' @return a list with following objects:
#' \itemize{
#'  \item meanTRJ : a list of matrices; every element is a matrix and 
#'  corresponds to a subject. Each matrix has p rows that represent p predicted 
#'  features. This is the prediction used by Xb as subject loading
#'  \item subTRJ : a list of matrices; every element is a matrix and 
#'  corresponds to a subject. Each matrix has p rows that represent p predicted 
#'  features. This is the prediction used by Xb as subject loading
#' }
#' @example examples/example_predict_supFTSVD.R
#' @export
predict.supFTSVD<-function(obj,designM,new_dat=NULL,TimeG=NULL){
  # number of components in the fitted model
  r<-ncol(obj$A.hat)
  # number of subjects for which prediction will be made
  ns<-nrow(designM)
  # Prediction of subject loading
  Xb<-designM%*%obj$Gamma
  # Predicted singular functions
  if(!is.null(new_dat)){
    newT=lapply(new_dat,function(u){
      u[1,]
    })
    
    newSF<-lapply(newT, function(nT){
      do.call(cbind,lapply(1:r, function(k){
        approx(x=obj$time,y=obj$Phi.hat[,k],xout = nT,rule = 2)$y
      }))
    })
    
    Hmat<-lapply(1:ns, function(i){
      as.matrix(sapply(1:r,function(k){
        as.numeric(outer(obj$B.hat[,k],newSF[[i]][,k]))
      }))
    })
    
    Hbeta<-lapply(1:ns, function(i){
      as.numeric(Hmat[[i]]%*%matrix(Xb[i,],ncol=1))
    })
    
    sw_vcov<-lapply(1:ns, function(i){
      solve(((t(Hmat[[i]])%*%Hmat[[i]])*(1/obj$Sigma2))+(diag(r)*(1/obj$Sigma2R)))
    })
    
    ## Vectorization of data
    Yvec<-lapply(1:ns, function(i){
      as.numeric(new_dat[[i]][-1,])
    })
    
    # conditional mean for all subjects
    if(r==1){
      conU<-as.matrix(sapply(1:ns,FUN = function(i){
        sw_vcov[[i]]%*%(t(Hmat[[i]]*(1/obj$Sigma2))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))
    } else{
      conU<-t(sapply(1:ns,FUN = function(i){
        sw_vcov[[i]]%*%(t(Hmat[[i]]*(1/obj$Sigma2))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))
    }
    Ahat<-Xb+conU
  } else{
    Ahat<-Xb
  }
  
  if(is.null(TimeG)){
    list("meanTRJ"=lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],obj$Phi.hat[,k])*Xb[i,k]*obj$Lambda[k]
      }))
    }),
    "subTRJ"=lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],obj$Phi.hat[,k])*Ahat[i,k]*obj$Lambda[k]
      }))
    }),
    "subEFF" = Ahat) 
  } else{
    predSF<-sapply(1:r, function(k){
      approx(x=obj$time,y=obj$Phi.hat[,k],xout = TimeG,rule = 2)$y
    })
    list("meanTRJ"=lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],predSF[,k])*Xb[i,k]*obj$Lambda[k]
      }))
    }),
    "subTRJ"=lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],predSF[,k])*Ahat[i,k]*obj$Lambda[k]
      }))
    }),
    "subEFF" = Ahat) 
  }
}


#' Prediction using the fitted objects from supFTSVD
#' 
#' @param obj fitted object by supFTSVD
#' @param new_dat is a list objects; each element is observed data for a 
#' subjects. The first row contains the time-points where measurements were 
#' taken. Thus, each element has (p+1) rows. 
#' @param TimeG Grid in the time domain where prediction will be made.
#' @returns a list with predicted trajectories. Elements correspond to the 
#' subjects
#' @example examples/example_predict_FTSVD.R
#' @export
predict.ftsvd<-function(obj,new_dat,TimeG=NULL){
  # number of components in the fitted model
  r<-ncol(obj$A.hat)
  # number of subjects for which prediction will be made
  ns<-length(new_dat)
  
  # Time points where measurement were taken
  newT=lapply(new_dat,function(u){
    u[1,]
  })
  
  
  # Predicted singular functions
  newSF<-lapply(newT, function(nT){
    sapply(1:r, function(k){
      approx(x=obj$time,y=obj$Phi.hat[,k],xout = nT,rule = 2)$y
    })
  })
  
  # matrix for estimating the subject loading
  Hmat<-lapply(1:ns, function(i){
    as.matrix(sapply(1:r,function(k){
      as.numeric(outer(obj$B.hat[,k],newSF[[i]][,k])*obj$Lambda[k])
    }))
  })
  
  ## Vectorization of data
  Yvec<-lapply(1:ns, function(i){
    as.numeric(new_dat[[i]][-1,])
  })
  
  ## Subject loading
  if(r==1){
    Ahat<-as.matrix(sapply(1:ns, function(i){
      as.numeric((solve(t(Hmat[[i]])%*%Hmat[[i]]))%*%(t(Hmat[[i]])%*%Yvec[[i]]))
    }))
  } else{
    Ahat<-t(sapply(1:ns, function(i){
      as.numeric((solve(t(Hmat[[i]])%*%Hmat[[i]]))%*%(t(Hmat[[i]])%*%Yvec[[i]]))
    }))
  }
  
  
  # Predicted data
  if(is.null(TimeG)){
    lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],obj$Phi.hat[,k])*Ahat[i,k]
      }))
    }) 
  } else{
    predSF<-sapply(1:r, function(k){
      approx(x=obj$time,y=obj$Phi.hat[,k],xout = TimeG,rule = 2)$y
    })
    
    lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],predSF[,k])*Ahat[i,k]
      }))
    })
  }
}


#' Format data table into input of ftsvd   
#'
#' @param taxon_table A table of read counts, with n rows for samples and p 
#' columns for taxa.
#' @param time_point The time stamp of each sample, relative to the start of 
#' the study. 
#' A length n vector.
#' @param subjectID The subject ID of each sample. A length n vector.
#' @param threshold A threshold for taxon filtering. 
#' Taxa with zero counts percentage >= threshold will be excluded.
#' @param pseudo_count A small number to add to all the counts before 
#' normalizing into proportions and log transformation.
#' @return Input for ftsvd. A list of matrices, each representing a subject data
#' @example examples/example_format_ftsvd.R
#' @export
format_ftsvd <- function(taxon_table, time_point, subjectID, threshold=0.95, 
                         feature_names=NULL, 
                         pseudo_count=0.5, transform="log_comp"){
  ntm <- which(table(subjectID)==1)
  if(length(ntm)>0)
    stop(paste('Please remove these subjects with only one time point:', 
               paste(names(ntm), collapse=', ')))
  if (length(subjectID)!=nrow(taxon_table)) 
    stop('length of subjectID does not match taxon_table!')
  if (length(time_point)!=nrow(taxon_table)) 
    stop('length of time_point does not match taxon_table!')
  # keep taxon that has non-zeros in >=1-threshold samples
  if (is.null(feature_names)){
    taxon_table <- taxon_table[,colMeans(taxon_table==0)<threshold]
  }else{
    taxon_table <- taxon_table[,feature_names]
  }
  if(transform=='log_comp'){
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- t(log(taxon_table/rowSums(taxon_table)))
  }else if(transform=='comp'){
    taxon_table <- taxon_table
    taxon_table <- t(taxon_table/rowSums(taxon_table))
  }else if(transform=='ast'){
    taxon_table <- taxon_table
    taxon_table <- t(asin(sqrt(taxon_table/rowSums(taxon_table))))
  }else if(transform=='clr'){
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- log(taxon_table/rowSums(taxon_table))
    taxon_table <- t(taxon_table-rowMeans(taxon_table))
  }else if(transform=='logit'){
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- t(taxon_table/rowSums(taxon_table))
    taxon_table <- log(taxon_table/(1-taxon_table))
  }else if(transform=='none'){
    taxon_table <- t(taxon_table)
  }else{
    print('Input transformation method is wrong! log_comp is applied instead')
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- t(log(taxon_table/rowSums(taxon_table)))
  }
  taxon_table <- rbind(time_point, taxon_table)
  rownames(taxon_table)[1] <- 'time_point'
  subID <- unique(subjectID)
  nsub <- length(subID)
  
  # construct list of data matrices, each element representing one subject
  datlist <- vector("list", length = nsub)
  names(datlist) <- subID
  
  # Each slice represents an individual (unequal sized matrix).
  for (i in 1:nsub){
    # print(i)
    datlist[[i]] <- taxon_table[, subjectID==subID[i]]
    datlist[[i]] <- datlist[[i]][,order(datlist[[i]][1,])]
    datlist[[i]] <- datlist[[i]][,!duplicated(datlist[[i]][1,])]
  }
  return(datlist)
}




#' Construct the Bernoulli kernel at a given set of knots and tiem points
#' 
#' @param x represents the observed time points
#' @param y set of reference points
#' @noRd
#' @export
bernoulli_kernel <- function(x, y){
  k1.x = x-0.5
  k1.y = y-0.5
  k2.x = 0.5*(k1.x^2-1/12)
  k2.y = 0.5*(k1.y^2-1/12)
  xy = abs(x %*% t(rep(1,length(y))) - rep(1,length(x)) %*% t(y))
  k4.xy = 1/24 * ((xy-0.5)^4 - 0.5*(xy-0.5)^2 + 7/240)
  kern.xy = k1.x %*% t(k1.y) + k2.x %*% t(k2.y) - k4.xy + 1
  return(kern.xy)
}

#' Perform a reproducing kernel hilbert space regression and being used by both 
#' ftsvd and supFTSVD function
#' 
#' @param Ly a list object contains the data
#' @param a.hat current estimate of subject loading
#' @param ind_vec a vector represents the ID
#' @param Kmat a matrix constructed using the bernoulli kernel
#' @param Kmat_output a matrix constructed using the bernoulli kernel and 
#' used for prediction
#' @param smooth a scalar represents the value of tuning parameter
#' @keywords internal
#' @noRd
#' @export
freg_rkhs <- function(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=1e-8){
  A <- Kmat
  for (i in 1:length(Ly)){
    #A[ind_vec==i,] <- A[ind_vec==i,]*a.hat[i]
    A[ind_vec==i,] <- A[ind_vec==i,]*a.hat[i]^2
  }
  #cvec <- Kmat%*%unlist(Ly)
  cvec <- unlist(Ly)
  
  #A.temp = crossprod(A) + smooth * Kmat
  #A.temp.eig <- eigen(A.temp, symmetric = TRUE)
  #A.d <- A.temp.eig$value
  #A.d[A.d<1e-10] <- 1e-10
  #beta <- ( (A.temp.eig$vector)%*%(t(A.temp.eig$vector)/A.d) ) %*% cvec
  
  A.temp <- A + smooth*diag(ncol(A))
  beta <- solve(A.temp)%*%cvec
  
  phi.est <- Kmat_output %*% beta
  return(phi.est)
}


#' Perform a reproducing kernel hilbert space regression that determines the 
#' smoothness parameter via g-fold cross-validation approach. It is used by both 
#' ftsvd and supFTSVD function when CVphi is set to be TRUE. 
#' 
#' @param Ly a list object contains the data
#' @param a.hat current estimate of subject loading
#' @param ind_vec a vector represents the ID
#' @param Kmat a matrix constructed using the bernoulli kernel
#' @param Kmat_output a matrix constructed using the bernoulli kernel and 
#' used for prediction
#' @param smooth a vector containing the grid of candidate values for the 
#' smoothness parameters
#' @param kfold number of folds to used
#' @param KInd a vector of index for making the folds; the function does the 
#' splitting automatically if set to NULL
#' @keywords internal
#' @noRd
#' @export
cv_freg_rkhs<-function(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth,kfold=5,
                       KInd=NULL){
  A <- Kmat
  for (i in 1:length(Ly)){
    A[ind_vec==i,] <- A[ind_vec==i,]*a.hat[i]^2
  }
  cvec <- unlist(Ly)
  if(is.null(KInd)){
    int<-length(cvec)%%kfold
    if(int==0){
      Kfold<-rep(1:kfold,each=length(cvec)/kfold)
    } else{
      Kfold<-c(rep(1:kfold,each=length(cvec)%/%kfold),1:(length(cvec)%%kfold))
    }
    
    KInd<-sample(Kfold,length(cvec),replace = FALSE)
  }
  
  KFres<-sapply(1:length(smooth),function(j){
    sapply(1:kfold,function(i){
      aind<-which(KInd!=i)
      A.temp <- A[aind,aind] + smooth[j]*diag(length(aind))
      beta <- solve(A.temp)%*%cvec[aind]
      cor(cvec[KInd==i],as.numeric((Kmat[KInd==i,aind]%*%beta)))^2
      #mean((cvec[KInd==i]-as.numeric((Kmat[KInd==i,aind]%*%beta)))^2)
    })
  })
  Sm<-smooth[which.max(colMeans(KFres))]
  print(Sm)
  #CV fit
  A.temp <- A + Sm*diag(ncol(A))
  beta <- solve(A.temp)%*%cvec
  
  phi.est <- Kmat_output %*% beta
  list(Sm,phi.est)
}


#' Function for generating longitudinal microbiome data
#' 
#' @param m_i number of measurements (repeated data) from subject i
#' @param M number of modalities
#' @param XMatrix design matrix for subject specific mean
#' @param Beta regression parameters for subject specific mean
#' @param Xi feature singular vector for all modalities ( a list of M elements)
#' @param PsiF is a list of M elements; each element is also a list with r 
#' elements which are functions
#' @param SubE_Var a list of M elements; each is a vector of r elements 
#' represent the variance parameters for subject-loading model
#' @param Data_Var a vector of M elements represent measurement error variances 
#' for different data modalities
#' @return a list with following iterms:
#' \itemize{
#'  \item data : a list of n elements; each correspond to data from a subject 
#'  and a matrix (p+1) by m_i matrix
#'  \item Xb : a n by r matrix of subject-specific mean for different low-rank 
#'  components
#'  \item zeta : a n by r matrix of subject loading for different low-rank 
#'  components 
#'  \item singF : singular functions used to generate data
#' }
#' @example examples/example_data_gen_supFTSVD.R
#' @export
data_gen_supFTSVD<-function(m_i,Xmatrix,Beta,Xi,PsiF,sing_val,SubE_Var,
                            Data_Var,Tgrid=seq(0,1,length.out=101)){
  
  n<-length(m_i)
  
  r<-ncol(Xi)
  
  EAval<-Xmatrix%*%as.matrix(Beta)
  
  miv_subE<-MASS::mvrnorm(n,mu = rep(0,r),Sigma=diag(SubE_Var))
  
  
  Zeta<-EAval+miv_subE
  
  obs_time_pos<-lapply(1:n,function(i){
    sort(sample(1:length(Tgrid),m_i[i],replace = FALSE))
  })
  
  
  obs_time<-lapply(1:n,function(i){
    Tgrid[obs_time_pos[[i]]] #sort(round(runif(n_i[i],TInt[1],TInt[2]),2))
  })
  
  
  singF<-lapply(1:n, function(i){
    sapply(1:r,function(k){
      PsiF[obs_time_pos[[i]],k] #PsiF[[m]][[k]](obs_time[[i]])
    })
  })
  
  
  dat<-lapply(1:n,function(i){
    rbind(obs_time[[i]],
          Reduce(`+`,lapply(1:r,function(k){
            (sing_val[k]*Zeta[i,k])*outer(Xi[,k],singF[[i]][,k])
          })) + rnorm(nrow(Xi)*length(obs_time[[i]]),mean = 0,sd = sqrt(Data_Var)))
  })
  
  list("data"=dat,
       "Xb" = EAval,
       "zeta"=Zeta,
       "singF"=singF)
}

#' Select a rank for supFTSVD by a modified version of cross-validation approach
#' in Lock and Li (2018)
#' 
#' @param datlist a list with n elements, each is a matrix of dimension p+1 
#' times m_i. The first row represents time points where measurements were 
#' obtained. Next p rows represent the function values observed at the time 
#' points
#' @param response a matrix of dimension n by q contains the design matrix for 
#' supervision of the subject loadings
#' @param interval represents the domain of the function. The default is NULL. 
#' @param ranks a vector of ranks on which grid search will be made 
#' @param resolution grid size on which the singular functions will be estimated
#' @param CVPhi a logical scalar representing whether cross-validation (CV) will be 
#' performed for determining the smoothness parameters involved with singular 
#' functions
#' @param K number of folds for the CV step
#' @param cvT number of initial iterations at which CV will be performed
#' @param smooth Smoothing parameter for RKHS norm when CVPhi is FALSE. With 
#' CVPhi=TRUE, a vector of numeric values for smoothness parameters for grid
#' search
#' @param maxiter Maximum number of iteration. Default: 20.
#' @param epsilon Convergence threshold for objective function value
#' @param KInd is a vector indices to construct the folds for CV. The default is
#' NULL that indicates that folds will be constructed randomly inside
#' @param rsvd_seed supFTSVD uses random svd for computational advantage while
#' setting inial values for feature loading vector. This seed number will be used
#' before setting the initial values. It confirms the reproducibility of your results.
#' @param conv_criteria the default is par which sets the threshold on the 
#' difference between updated and current estimate at M step to stop the iteration. 
#' It is computationally expensive; thus we recommend to use "cond_lik" instead of "par". 
#' The result is not sensitive. 
#' @param rank_fold number of folds to be used in the cross-validation
#' @param rfolds_seed seed number to be used in random splitting of the data
#' @param rc_thresh threshold on test-data likelihood for rank selection
#' @return a list with following objects:
#' \itemize{
#'  \item cv_res : detailed results for different folds
#'  \item opt_r : optimal rank
#' }
#' @example examples/example_cv_rank_supFTSVD.R
#' @export
cv_rank_supFTSVD<-function(datlist, response, interval = NULL, 
                           ranks = c(1:5),resolution=100, CVPhi=FALSE, K=5, cvT=5, 
                           smooth=1e-8,maxiter=20, epsilon=1e-4,KInd=NULL,
                           rsvd_seed=100,conv_criteria="par",
                           rank_fold=5,rfolds_seed=50,rc_thresh=0.05,stratum=NULL){
  n<-length(datlist)
  if(!is.null(stratum)){
    grp_frq<-as.numeric(table(stratum))
    lv_grp<-as.character(unique(stratum))
    set.seed(rfolds_seed)
    SJ_ind<-do.call(c,lapply(seq_len(length(lv_grp)),function(i){
      sample(1:grp_frq[i],grp_frq[i],replace=FALSE)
    }))
    rfolds_indx<-do.call(c,lapply(seq_len(length(lv_grp)),function(i){
      sj_ind<-sample(1:grp_frq[i],grp_frq[i],replace=FALSE)
      c(rep(1:rank_fold,each=grp_frq[i]%/%rank_fold),sample(1:rank_fold,grp_frq[i]%%rank_fold,replace = FALSE))
    }))
  } else{
    set.seed(rfolds_seed)
    SJ_ind<-sample(1:n,n,replace=FALSE)
    rfolds_indx<-c(rep(1:rank_fold,each=n%/%rank_fold),sample(1:rank_fold,n%%rank_fold,replace = FALSE))
  }
  rank_res<-NULL
  for(r in ranks){
    crank<-sapply(1:rank_fold, function(g){
      train_dat<-datlist[sort(SJ_ind[which(rfolds_indx!=g)])]
      test_dat<-datlist[sort(SJ_ind[which(rfolds_indx==g)])]
      fit_model<-supFTSVD(datlist=train_dat, response=response[sort(SJ_ind[which(rfolds_indx!=g)]),], interval = interval, 
                          r = r,resolution=resolution, CVPhi=CVPhi, K=K, cvT=cvT,
                          smooth=smooth,maxiter=maxiter, epsilon=epsilon,KInd=KInd,
                          rsvd_seed=rsvd_seed,conv_criteria=conv_criteria)
      predFIT<-predict.supFTSVD(obj = fit_model,
                                designM = response[sort(SJ_ind[which(rfolds_indx==g)]),],
                                new_dat = test_dat)
      
      test_Xb<-response[sort(SJ_ind[which(rfolds_indx==g)]),]%*%fit_model$Gamma
      test_conU<-predFIT$subEFF-test_Xb
      
      test_r_llik<-norm_llik_supFTSVD(datlist = test_dat,
                                      x_matrix = response[sort(SJ_ind[which(rfolds_indx==g)]),],
                                      comp_beta = fit_model$Gamma,
                                      subj_dev = test_conU,
                                      feat_loading = fit_model$B.hat,
                                      sing_func = fit_model$Phi.hat,
                                      sing_func_arg = fit_model$time,
                                      sing_val = fit_model$Lambda,
                                      sigR = fit_model$Sigma2R,
                                      sig = fit_model$Sigma2)
      
    })
    if(r==1){
      rc<-1
    } else{
      rc<-(mean(crank)-mean(orank))/abs(mean(orank))
    }
    if(rc<rc_thresh)
      break
    orank<-crank
    rank_res<-rbind(rank_res,crank)
  }
  list("cv_res"=rank_res,"opt_r"=ranks[nrow(rank_res)])
}


#' Compute the value of normal likelihood for observed high-dimensional multivariate functional data
#' 
#' @param datalist a list with n elements, each is a matrix of dimension p+1 
#' times m_i. The first row represents time points where measurements were 
#' obtained. Next p rows represent the function values observed at the time 
#' points
#' @param iterval represents the domain of the function. The default is NULL. 
#' @param x_matrix a matrix of dimension n by q contains the design matrix for 
#' supervision of the subject loadings
#' @param comp_beta parameter vector associated with x_matrix
#' @param subj_dev a matrix of n by r contains subject-specific-deviation
#' @param feat_loading a matrix of p by r contains feature loadings
#' @param sing_func a matrix represents the value of r singular functions on a 
#' smooth grid in the interval
#' @param sing_val a vector of singular values
#' @param sigR a vector of dimension r with variances of subject-loading model 
#' as elements
#' @param sig error variance
#' @param smooth_par smoothing parameters used to compute the singular functions;
#' it is an optional argument
#' @returns the value of the complete data likelihood for the given data and 
#' parameters based on normality assumption
#' @importFrom stats approx coef cor lm rnorm
#' @example examples/example_norm_join_llik_supFTSVD.R
#' @export
norm_llik_supFTSVD<-function(datlist,x_matrix,comp_beta,subj_dev,
                             feat_loading,sing_func,sing_func_arg,sing_val,sigR,sig,smooth_par=0){
  pj<-nrow(datlist[[1]])-1
  n = length(datlist)
  r<-ncol(subj_dev)
  obs_time<-lapply(datlist,function(u){u[1,]})
  m_i<-sapply(obs_time,length)
  
  
  # response in matrix form
  if(!is.matrix(x_matrix)){
    x_matrix<-matrix(x_matrix,ncol=1)
  }
  
  if(!is.matrix(comp_beta)){
    comp_beta<-matrix(comp_beta,ncol=1)
  }
  
  
  # Subject wise vectorized data 
  Yvec<-lapply(1:n, function(i){
    as.numeric(datlist[[i]][-1,])
  })
  
  Lt<-obs_time
  
  
  allT<-do.call(c,Lt)
  
  pXI<-apply(sing_func,2,function(u){
    approx(sing_func_arg,u,xout = allT,rule=2)$y
  })
  
  
  # initial values for error variance
  Mi<-cumsum(m_i)
  
  prdXI<-lapply(1:n,function(w){
    if(w ==1){
      as.matrix(pXI[1:Mi[1],])
    } else{
      as.matrix(pXI[(Mi[w-1]+1):(Mi[w]),])
    }
  })
  
  
  Ahat<-(x_matrix%*%comp_beta)+subj_dev
  Xb<-x_matrix%*%comp_beta
  
  # Necessary matrices (updated)
  Hmat<-lapply(1:n, function(i){
    sapply(1:r, function(k){
      as.numeric(outer(feat_loading[,k],prdXI[[i]][,k]))*sing_val[k]
    })
  })
  
  
  ## Subject wise super X matrix
  Hbeta<-lapply(1:n, function(i){
    as.numeric(Hmat[[i]]%*%matrix(Xb[i,],ncol=1))
  })
  
  # another necessary matrix
  HUmat<-lapply(1:n, function(i){
    as.numeric(Hmat[[i]]%*%matrix(subj_dev[i,],ncol=1))
  })
  
  ####
  sigma_R<-diag(ncol(subj_dev))
  diag(sigma_R)<-1/sigR
  
  # Subject-specific conditional covariances
  sw_vcov<-lapply(Hmat,function(u){
    solve(crossprod(u)*(1/sig)+sigma_R)
  })
  
  # objective function
  Reduce(`+`,lapply(1:n, function(i){
    t1<-((-pj*m_i[i]/2)*(log(sig)))-sum((1/2)*log(sigR))
    t2<-sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2/(sig))
    t3<-sum(diag((crossprod(Hmat[[i]])*(1/sig))%*%sw_vcov[[i]]))
    t4<-sum(diag(sigma_R%*%sw_vcov[[i]]))
    t5<-sum((subj_dev[i,]^2/sigR))
    t1-(t2+t3+t4+t5)/2
  }))+sum(smooth_par)
}