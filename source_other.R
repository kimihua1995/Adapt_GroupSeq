crosslist=function(b=list(a1=c(2,3),a2=c(2,4),a3=c(0,1))){
  bn=length(names(b))
  an=1
  for (i in 1:bn){an=an*length(b[[i]])}
  amatrix=matrix(0,nrow=an,ncol=bn)
  i=0;nsa=an;nsb=1
  while (i<bn){
    i=i+1
    nsa=nsa/length(b[[i]])
    if (i>=2)nsb=nsb*length(b[[i-1]])
    amatrix[,i]=rep(rep(b[[i]],each=nsb),times=nsa)
  }
  df=data.frame(amatrix)
  colnames(df)=names(b)
  
  return(df)
}




## create the correlation matrix
create.m.corr <- function(w,I){
  l <- length(I)
  smatrix <- diag(1, l)
  for (i in 1:(l-1)){
    for (j in (i+1):l){
      smatrix[i,j] <- smatrix[j,i] <- sqrt(I[i]/I[j])
    }
  }
  smatrix <- w[1]^2 + w[2]^2*smatrix
  return(smatrix)
}






## calculate the test statistics for each hypothesis
X.IA <- function(dat, t_IAs, group, trt){
  l <- length(t_IAs)
  X_PFS <- X_OS <- rep(NA, l)
  var_PFS <- var_OS <- rep(NA, l-1)
  dat0 <- dat[dat$group %in% group & dat$trt %in% trt, ]
  
  dat1 <- dat0[dat0$t.enroll <= t_IAs[1], ]
  dat1 <- dat.modify(dat1, t_IAs[1])
  LR1_PFS <- logrank.test(dat1$t2, dat1$status2, factor(dat1$trt))
  LR1_OS <- logrank.test(dat1$t3, dat1$status3, factor(dat1$trt))
  X_PFS[1] <- LR1_PFS$test$z
  X_OS[1] <- LR1_OS$test$z
  
  for (i in 2:l){
    dat1 <- dat.modify(dat1, t_IAs[i])
    LR1_PFS2 <- logrank.test(dat1$t2, dat1$status2, factor(dat1$trt))
    LR1_OS2 <- logrank.test(dat1$t3, dat1$status3, factor(dat1$trt))
    
    dat2 <- dat0[dat0$t.enroll <= t_IAs[i] & dat0$t.enroll > t_IAs[1], ]
    dat2 <- dat.modify(dat2, t_IAs[i])
    LR2_PFS <-  logrank.test(dat2$t2,dat2$status2,factor(dat2$trt))
    LR2_OS <-  logrank.test(dat2$t3,dat2$status3,factor(dat2$trt))
    
    # see section 4.2 for details
    var_PFS[i-1] <- LR1_PFS2$var - LR1_PFS$var + LR2_PFS$var
    X_PFS[i] <- (LR1_PFS2$test$z*sqrt(LR1_PFS2$var) - 
                 LR1_PFS$test$z*sqrt(LR1_PFS$var) +
                 LR2_PFS$test$z*sqrt(LR2_PFS$var))/sqrt(var_PFS[i-1])
    var_OS[i-1] <- LR1_OS2$var - LR1_OS$var + LR2_OS$var
    X_OS[i] <- (LR1_OS2$test$z*sqrt(LR1_OS2$var) - 
                LR1_OS$test$z*sqrt(LR1_OS$var) +
                LR2_OS$test$z*sqrt(LR2_OS$var))/sqrt(var_OS[i-1])
  }
  
  return(list(X_PFS = X_PFS, X_OS = X_OS,
              var_PFS = var_PFS, var_OS = var_OS))
}



# Create sets of all intersection hypotheses from I_1
powerset <- function(x) {
  sets <- lapply(1:(length(x)), function(i) combn(x, i, simplify = F))
  unlist(sets, recursive = F)
}




# generate the IA time given planned number of events
get.IA <- function(dat, group, trt, nevents, event){
  dats <- dat[dat$group %in% group & dat$trt %in% trt, ]
  if (event == "PFS") {
    y <- dats$y2
  }else if (event == "OS"){
    y <- dats$y3
  }
  t.time <- sort(dats$t.enroll + y)
  n <- NULL
  for (i in nevents[1]:nrow(dats)){
    dats$t.censor <- pmin(dats$t.loss, t.time[i] - dats$t.enroll)
    n <- c(n, sum(y <= dats$t.censor))
  }
  t_IAs <- NULL
  for (j in 1:length(nevents)){
    pos <- min(which(n >= nevents[j]))
    t_IAs <- c(t_IAs, t.time[nevents[1] + pos - 1])
  }
  return(t_IAs)
}

