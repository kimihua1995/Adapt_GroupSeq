update.graph=function(S1=c(2,3),W0=c(0.5,0.5,0,0),G0=rbind(c(0,0,1,0),c(0,0,0,1),c(0,1,0,0),c(1,0,0,0)),S0=seq(1,length(W0),by=1)){
  #S1 new set of hypotheses #S1 must be a proper and non-empty subset of S0, S1 must be sorted increasingly
  #W0 initial weights
  #G0 initial matrix of length(W0)*length(W0)
  #S0 initial set of hypotheses from 1 to n
  SS=setdiff(S0,S1) 
  nss=length(SS)
  WT=W0;GT=G0
  if (nss>0){
    for (j in 1:nss){
      sj=SS[j]+1-j # this ensures we select the correct index after the graph is updated 
      WT=WT + WT[sj] * GT[sj,]  # update W 
      GT.new=GT
      for (l in 1:nrow(GT.new)){  # update G
        for (k in 1:ncol(GT.new)){
          GT.new[l,k] <- (GT[l,k] + GT[l,sj]*GT[sj,k])/(1 - GT[l,sj]*GT[sj,l])*(1 - (l == k))
        }
      }
      WT <- WT[-sj]
      GT <- as.matrix(GT.new[-sj,-sj])
    }
  }
  list(S1=S1,W1=WT,G1=GT)
}


graphical=function(X,W,G,zcorr,salpha,tIA){
  #X: vector of test-statistics
  #W: initial weights of length(p)
  #G: initial matrix of length(p)*length(p)
  #tIA: column of indicator for critical value (at which IA)
  nX=length(W)
  I=1:nX
  have_rej <- TRUE
  if (nX > 1){
    Xt <- X[,tIA]
  }else{
    Xt <- X[tIA]
  }
  
  while (have_rej){
    alpha_star <- rep(NA, length(I))
    for (i in 1:length(I)){
      if (W[i] > 0){
        if (length(I) > 1){
          temp <- crit.value.g(salpha[i,]*W[i], zcorr[[i]], sided=1)
        }else{
          temp <- crit.value.g(salpha*W[i], zcorr[[i]], sided=1)
        }
        if ((I[i] %in% 1:4) & tIA == 3){
          alpha_star[i] <- temp[tIA-1]
        }else{
          alpha_star[i] <- temp[tIA]
        }
      }else if (W[i] == 0){
        alpha_star[i] <- Inf
      }
    }
    rej_judge <- (Xt > alpha_star)
    have_rej <- sum(rej_judge) > 0
    if (have_rej){
      rej_pos <- which(rej_judge == TRUE)
      j <- rej_pos[1]  # pick the first rejected one
      W <- W + W[j] * G[j,]  # update W
      G_new <- G
      for (l in 1:nrow(G_new)){  # update G
        for (k in 1:ncol(G_new)){
          G_new[l,k] <- (G[l,k] + G[l,j]*G[j,k])/(1 - G[l,j]*G[j,l])*(1 - (l == k))
        }
      }
      # update I <- I/{j}
      I <- I[-j]
      W <- W[-j]
      if (length(I) > 0){
        Xt <- Xt[-j]
        G <- as.matrix(G_new[-j,-j])
        salpha <- salpha[-j,]
        zcorr <- zcorr[-j]
      }else {have_rej <- FALSE}
    }
  }
  rej_h <- rep(1,nX)
  rej_h[I] <- 0
  # return a vector of length of hypotheses number
  # 1: rejected, 0: not rejected
  if (length(I) > 0){
    return(list(rej = rej_h, W = W, G = G))
  }else{
    return(list(rej = rej_h, W = NULL, G = NULL))
  }
}


