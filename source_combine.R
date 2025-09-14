sim.design <- function(data.list, w, salpha_PFS, salpha_OS, B, nPFS, nOS, W, G, alpha){
# w: pre-specified vector of weight functions in test statistics (can be extended to be different for each hypothesis)
# salpha_PFS, salpha_OS: alpha-spending list for each hypothesis for PFS or OS
# B: number of simulated data
# nPFS, nOS: planned number of events for PFS or OS
# W: initial weight for the graph
# G: transition weight matrix of the graph
# alpha: one-sided significance level for type-1 error rate
    
  rej_multi1 <- rej_multi2 <- rej_multi3 <- array(0, dim = c(8, 3, B))
  adapt <- rep(1,B) # 1: all, 2: drop T1, 3: drop BM-, 4: stop all
  t_IAs_list <- matrix(NA, B, 4)
  
  
  for (b in 1:B){
    #if (!b%%1000) print(b)
    # 1) Generate the data
    dat <- data.list[[b]]
    
    # 2) Calculate the test statistics X_ij and its variance for each hypothesis and IA
    t_IAs <- c(get.IA(dat, c(1), c(0,2), c(ceiling(0.25*nPFS),ceiling(0.75*nPFS), nPFS), "PFS"),
               get.IA(dat, c(1), c(0,2), c(nOS), "OS"))
    t_IAs_list[b,] <- t_IAs
    X_PFS <- X_OS <- matrix(NA, 4, length(t_IAs)) # Test statistics include Stage-1
    var_PFS <- var_OS <- matrix(NA, 4, length(t_IAs) - 1) # Correlation matrix
    group.list <- list(c(1), c(1), c(0,1), c(0,1))
    trt.list <- list(c(0,2), c(0,1), c(0,2), c(0,1))
    ## T2 vs. T0 in BM+
    ## T1 vs. T0 in BM+
    ## T2 vs. T0 in overall
    ## T1 vs. T0 in overall
    for (i in 1:4){
      X_IAi <- X.IA(dat, t_IAs, group=group.list[[i]], trt = trt.list[[i]])
      X_PFS[i,] <- X_IAi$X_PFS
      X_OS[i,] <- X_IAi$X_OS
      var_PFS[i,] <- X_IAi$var_PFS
      var_OS[i,] <- X_IAi$var_OS
    }
    X_IAs <- rbind(X_PFS, X_OS)
    var_IAs <- rbind(var_PFS, var_OS)
    
    # 3) Adaptation in Stage-1
    ## This part can be modified based on specific study objective
    ## The return vector I2 contains the undropped hypotheses to Stage-2
    ## adapt: vector of adaption condition
    dat1 <- dat[dat$t.enroll <= t_IAs[1],]
    dat1_plus <- dat1[dat1$group==1,]
    ORR_20_plus <- mean(dat1_plus$y1[dat1_plus$trt==2])-mean(dat1_plus$y1[dat1_plus$trt==0])
    ORR_20 <- mean(dat1$y1[dat1$trt==2]) - mean(dat1$y1[dat1$trt==0])
    ORR_10_plus <- mean(dat1_plus$y1[dat1_plus$trt==1])-mean(dat1_plus$y1[dat1_plus$trt==0])
    if (ORR_20_plus < 0.15){
      adapt[b] <- 5; next
    }else if (ORR_20 < 0.1){
      if (ORR_10_plus < 0.1){
        adapt[b] <- 4; I2 <- c(1,5)
      }else{
        adapt[b] <- 3; I2 <- c(1,2,5,6)
      }
    }else if (ORR_10_plus < 0.1){
      adapt[b] <- 2; I2 <- c(1,3,5,7)
    }else{
      adapt[b] <- 1; I2 <- 1:8
    }
  
    
    
    
    # 4) Proposed method: calculate the group-sequential p-values in Stage-2
    ## i) Calculate the nominal p-values in Stage-1
    pI <- 1 - pnorm(X_IAs[,1])
    
    ## ii) Calculate the critical values g_ij for each hypothesis
    ### Calculate the correlation matrix for Z_ij
    zcorr <- list()
    for (i in 1:8){
      zcorr[[i]] <- create.m.corr(w, var_IAs[i,])
    }
    ### Calculate critical values
    g_crit <- matrix(NA, 8, length(t_IAs)-1)
    salpha <- rbind(salpha_PFS, salpha_OS)
    for (i in 1:8){
      if (W[i] > 0){
        g_crit[i,] <- crit.value.g(salpha[i,]*W[i], zcorr[[i]], sided=1)
      }else if (W[i] == 0){
        g_crit[i,] <- Inf
      }
    }
    
    
    ## iii) Calculate the critical values c_ij and prepare the sub-distribution function
    ### Calculate the critical values c_ij
    c_crit <- (g_crit - w[1]*matrix(X_IAs[,1], 8, length(t_IAs) - 1))/w[2]
    
    ### Determine the correlation matrix list for X_ij^[2]
    xcorr <- list()
    for (i in 1:8){
      xcorr[[i]] <- create.m.corr(c(0,1), var_IAs[i,])
    }
    
    ### Calculate the cumulative errors alpha^c using the sub-distribution
    cum_error <- matrix(NA, 8, length(t_IAs)-1)
    for (i in 1:8){
      if (W[i] > 0){
        cum_error[i,] <- cum.error(c_crit[i,], xcorr[[i]], 1)
      }else if (W[i] == 0){
        cum_error[i,] <- 0
      }
    }
    
    ## iv) Calculate the group-seq p-values for each hypothesis in Stage II
    xcorr_PFS <- list()
    for (k in 1:4){
      xcorr_PFS[[k]] <- xcorr[[k]][1:2,1:2]
    }
    gsp_PFS <- calgsp2(xm = X_IAs[1:4,2:3],
                       critm = c_crit[1:4,1:2],
                       alpham = cum_error[1:4,1:2],
                       matrix.list = xcorr_PFS[1:4],
                       sided = rep(1,4))$pm
    gsp_PFS <- cbind(gsp_PFS, gsp_PFS[,2])
    gsp_OS <- calgsp2(xm = X_IAs[5:8,-1],
                       critm = c_crit[5:8,],
                       alpham = cum_error[5:8,],
                       matrix.list = xcorr[5:8],
                       sided = rep(1,4))$pm
    gsp_all <- rbind(gsp_PFS, gsp_OS)
    
    
    
  
    ## v) Group-Sequential procedure using group-seq p-values
    ### Weighted Simes and Bonferroni tests are used for closed testing procedure
    I1 <- 1:8
    pset=powerset(I1)
    rej_h_B <- rej_h_S <- matrix(0, 8, length(t_IAs)-1)  # Matrix of rejection indicator for each hypothesis (row) at each IA in Stage-2 (column), based on weighted Bonferroi or Simes
    rej_h_B[I2,] <- rej_h_S[I2,] <- 1 # the dropped hypotheses always = 0, the initial value for undropped hypotheses = 1.
    for (h_index in I2){
      for (k in 1:length(pset)){
        if (is.element(h_index,pset[[k]])) {
          # For Stage-1
          W_kI <- update.graph(S1 = pset[[k]], W0 = W, G0 = G, S0 = I1)$W1
          pI_k <- pI[pset[[k]]]
          p1_adj_k_S=weighted.simes.BH(pv=pI_k,wv=W_kI,alpha=alpha)$adj.p
          p1_adj_k_B=weighted.bonferroni(pv=pI_k,wv=W_kI,alpha=alpha)$adj.p
          p1_adj_k_S=min(p1_adj_k_S,1)
          p1_adj_k_B=min(p1_adj_k_B,1)
          # For Stage-2
          W_kII <- update.graph(S1 = intersect(pset[[k]],I2), W0 = W, G0 = G, S0 = I1)$W1
          for (j in 2:length(t_IAs)){
            gsp_k <- gsp_all[intersect(pset[[k]],I2),j-1]
            p2_adj_k_S=weighted.simes.BH(pv=gsp_k,wv=W_kII,alpha=alpha)$adj.p
            p2_adj_k_B=weighted.bonferroni(pv=gsp_k,wv=W_kII,alpha=alpha)$adj.p
            p2_adj_k_S=min(p2_adj_k_S,1)
            p2_adj_k_B=min(p2_adj_k_B,1)
            # P-value combination
            Z_k_S <- w[1]*qnorm(1-p1_adj_k_S) + w[2]*qnorm(1-p2_adj_k_S) # Using Simes
            Z_k_B <- w[1]*qnorm(1-p1_adj_k_B) + w[2]*qnorm(1-p2_adj_k_B) # Using Bonferroni
            # Update the rejection matrix through closed testing
            rej_h_S[h_index,j-1] <- rej_h_S[h_index,j-1] * (Z_k_S >= qnorm(1-alpha))
            rej_h_B[h_index,j-1] <- rej_h_B[h_index,j-1] * (Z_k_B >= qnorm(1-alpha))
          }
        }
      }
    }

    
    
    # 5) Group-Sequential procedure using Sugitani et al. 2016
    ## i) Update the graph after adaption, using sub-graph of G (i.e., G_I2 and W_I2)
    ### This step controls the FWER, but is too conservative
    G_I2 <-  as.matrix(G[I2, I2])
    W_I2 <- W[I2]
    
    ## ii) Generate the weighted test statistics
    Z_G <- w[1]*X_IAs[,1] + w[2]*X_IAs[,-1]
    Z_G[1:4,3] <- Z_G[1:4,2] 
    Z_I2 <- Z_G[I2,]
    zcorr_I2 <- zcorr[I2]
    salpha_I2 <- salpha[I2,]
    
    rej_h_G <- matrix(0, 8, length(t_IAs)-1)  # Matrix of rejection indicator for each hypothesis (row) at each IA in Stage-2 (column)
    for (j in 1:(length(t_IAs)-1)){
      res_G <- graphical(Z_I2,W_I2,G_I2,zcorr_I2,salpha_I2,j)
      rej_j <- res_G$rej
      W_I2 <- res_G$W
      G_I2 <- res_G$G
      
      rej_h_G[I2,j:(length(t_IAs)-1)] <- rej_j
      I2 <- I2[rej_j == 0]
      
      if (length(I2) > 0){
        Z_I2 <- Z_G[I2,]
        salpha_I2 <- salpha[I2,]
        zcorr_I2 <- zcorr[I2]
      }else if (length(I2) == 0){ break }
    }
    
    
    
    # 6) Update in the rejection array for each approach
    rej_multi1[,,b] <- rej_h_S
    rej_multi2[,,b] <- rej_h_B
    rej_multi3[,,b] <- rej_h_G
    
  
    
  }
  
  return(list(rej1 = rej_multi1, 
              rej2 = rej_multi2,
              rej3 = rej_multi3,
              adapt = adapt,
              t_IAs = t_IAs_list))
}
