# calculate the critical values of g_is
crit.value.g <- function(salpha, smatrix, sided, tol=1e-10, alpha.tol=tol/10){
  #salpha: cumulative alpha levels
  #smatrix: general correlation matrix
  #sided: 1: reject if large; -1: reject if small; 0: reject is absolute value is large
  ns=length(salpha)
  sc=rep(0,ns)
  
  if (sided==1) {
    sc[1]=qnorm(1-salpha[1])
  }else if(sided==0) {
    sc[1]=qnorm(1-salpha[1]/2)
  }else if (sided==-1) {
    sc[1]=qnorm(salpha[1])
  }
  
  if (ns >= 2){
    for (j in 2:ns){
      if (salpha[j]-salpha[j-1]>alpha.tol){
        xfunc=function(x,xcrit=rep(0,j-1)){
          if (sided==0){
            temp=OpenMx::omxMnor(covariance=smatrix[1:j,1:j], means=rep(0,j), lbound=c(-xcrit,x), ubound=c(xcrit,Inf))*2-(salpha[j]-salpha[j-1])
          }
          else if (sided==1){
            temp=OpenMx::omxMnor(covariance=smatrix[1:j,1:j], means=rep(0,j), lbound=c(rep(-Inf,j-1),x), ubound=c(xcrit,Inf))-(salpha[j]-salpha[j-1])
          }
          else if (sided==-1){
            temp=OpenMx::omxMnor(covariance=smatrix[1:j,1:j], means=rep(0,j), lbound=c(xcrit,-Inf), ubound=c(rep(Inf,j-1),x))-(salpha[j]-salpha[j-1])
          }
          temp
        }
        if (sided==0) tinterval=c(0,qnorm(1-alpha.tol/2))
        else if (sided==1) tinterval=c(0,qnorm(1-alpha.tol))
        else if (sided==-1) tinterval=c(-qnorm(1-alpha.tol),0)
        sc[j]<- uniroot(xfunc, interval=tinterval, tol = tol, xcrit=sc[1:(j-1)])$root
      }
      else if (salpha[j]-salpha[j-1]<=alpha.tol){
        sc[j]=(sided==0)*qnorm(1-alpha.tol/2)+(sided==1)*qnorm(1-alpha.tol)+(sided==-1)*qnorm(alpha.tol)
      }
    }
  }
  
  return(sc)
}



cum.error <- function(crit, smatrix, sided = 1){
  nx <- length(crit)
  if (sided==1){
    crit=-crit
  }
  calpha <- rep(0,nx)
  if (sided ==-1 | sided == 1){
    calpha[1] <- pnorm(crit[1])
    for (s in 2:nx){
      calpha[s] <- OpenMx::omxMnor(covariance=smatrix[1:s,1:s], means=rep(0,s), 
                        lbound=c(crit[1:(s-1)],-Inf), ubound=c(rep(Inf,s-1),crit[s])) +
        calpha[s-1]
    }
  }else if (sided == 0){
    calpha[1] <- 2*pnorm(-crit[1])
    for (s in 2:nx){
      calpha[s] <- OpenMx::omxMnor(covariance=smatrix[1:s,1:s], means=rep(0,s), 
                        lbound=c(-crit[1:(s-1)],crit[s]), ubound=c(crit[1:(s-1)],Inf))*2 +
        calpha[s-1]
    }
  }
  return(calpha)
}



# Calculate the group-sequential p-values
calgsp2=function(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
                 critm=matrix(rep(qnorm(c(0.02,0.03,0.05)),each=2),ncol=3,nrow=2),
                 alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
                 matrix.list=list(diag(3),diag(3)),sided=rep(-1,2)){
  #xm: matrix of test statistics, each row is for one hypothesis and is assumed to be multivariate normal
  #alpham: matrix of alpha-levels, each row is for one hypothesis at different time points, for each row, alpha levels must be non-decreasing
  #critm: matrix of critical values, each row is for one hypothesis at different time points
  #matrix.list: list of correlation matrix corresponding to each hypotheses
  #sided: 1: reject if large; -1: reject if small; 0: reject is absolute value is large
  #This program is faster as the critical values are provided
  n=nrow(xm)
  s=ncol(xm)
  pm=xm
  
  for (j in 1:n){
    crit=critm[j,]
    dcorr=matrix.list[[j]]
    pm[j,]=pmin(calgsp.fast(tx=xm[j,],
                            tcrit=critm[j,],
                            calpha=alpham[j,],
                            smatrix=matrix.list[[j]],
                            sided=sided[j])$pa,1)
  }
  list(pm=pm)
}



calgsp.fast=function(tx=qnorm(c(0.03,0.04,0.01)),
                     tcrit=qnorm(c(0.01,0.02,0.025)),
                     calpha=c(0.01,0.02,0.025),
                     smatrix=diag(3),sided=-1){
  #tx: test statistics, assumed to be multivariate normal
  #tcrit: critical values
  #calpha: cumulative alpha levels
  #smatrix: correlation matrix of the test statistics
  #sided: 1: reject if large; -1: reject if small; 0: reject is absolute value is large
  
  nx=length(tx)
  
  aa=pa=rep(0,nx)
  if (sided==-1){
    x=tx;crit=tcrit
  }else if (sided==1){
    x=-tx;crit=-tcrit
  }else if (sided==0){
      x=abs(tx);crit=tcrit
  }
  
  #calpha <- rep(0,nx)
  #if (sided==-1 | sided == 1){
  #  calpha[1] <- pnorm(crit[1])
  #  for (s in 2:nx){
  #    calpha[s] <- OpenMx::omxMnor(covariance=smatrix[1:s,1:s], means=rep(0,s), 
  #                      lbound=c(crit[1:(s-1)],-Inf), ubound=c(rep(Inf,s-1),crit[s])) +
  #      calpha[s-1]
  #  }
  #}else if (sided == 0){
  #  calpha[1] <- 2*pnorm(-crit[1])
  #  for (s in 2:nx){
  #    calpha[s] <- OpenMx::omxMnor(covariance=smatrix[1:s,1:s], means=rep(0,s), 
  #                      lbound=c(-crit[1:(s-1)],crit[s]), ubound=c(crit[1:(s-1)],Inf))*2 +
  #      calpha[s-1]
  #  }
  #}
  
  
  if (sided==-1|sided==1){
    aa[1]=pnorm(x[1])
    pa[1]=(x[1]<=crit[1])*aa[1]+(x[1]>crit[1])*1
    if (nx>=3){
      for (i in 2:(nx-1)){
        aa[i]=OpenMx::omxMnor(covariance=smatrix[1:i,1:i], means=rep(0,i), 
                              lbound=c(crit[1:(i-1)],-Inf), ubound=c(rep(Inf,i-1),x[i]))+
        calpha[i-1]
        #aa[i]=mvtnorm::pmvnorm(lower=c(crit[1:(i-1)],-Inf),upper=c(rep(Inf,i-1),x[i]),sigma=smatrix[1:i,1:i])+salpha[i-1]
        pa[i]=pa[i-1]-prod(x[1:(i-1)]>crit[1:(i-1)])*(x[i]<=crit[i])*(1-aa[i])
      }
    }
    aa[nx]=OpenMx::omxMnor(covariance=smatrix[1:nx,1:nx], means=rep(0,nx), 
                           lbound=c(crit[1:(nx-1)],-Inf), ubound=c(rep(Inf,nx-1),x[nx]))+
      calpha[nx-1]
    #aa[nx]=mvtnorm::pmvnorm(lower=c(crit[1:(nx-1)],-Inf),upper=c(rep(Inf,nx-1),x[nx]),sigma=smatrix[1:nx,1:nx])+salpha[nx-1]
    pa[nx]=pa[nx-1]-prod(x[1:(nx-1)]>crit[1:(nx-1)])*(1-aa[nx])
  }
  else if (sided==0){
    aa[1]=2*pnorm(-x[1])
    pa[1]=(x[1]>=crit[1])*aa[1]+(x[1]<crit[1])*1
    if (nx>=3){
      for (i in 2:(nx-1)){
        aa[i]=OpenMx::omxMnor(covariance=smatrix[1:i,1:i], means=rep(0,i), 
                              lbound=c(-crit[1:(i-1)],x[i]), ubound=c(crit[1:(i-1)],Inf))*2+
          calpha[i-1]
        #aa[i]=mvtnorm::pmvnorm(lower=c(-crit[1:(i-1)],x[i]),upper=c(crit[1:(i-1)],Inf),sigma=smatrix[1:i,1:i])*2+salpha[i-1]
        pa[i]=pa[i-1]-prod(x[1:(i-1)]<crit[1:(i-1)])*(x[i]>=crit[i])*(1-aa[i])
      }
    }
    aa[nx]=OpenMx::omxMnor(covariance=smatrix[1:nx,1:nx], means=rep(0,nx), 
                           lbound=c(-crit[1:(nx-1)],x[nx]), ubound=c(crit[1:(nx-1)],Inf))*2+
      calpha[nx-1]
    #aa[nx]=mvtnorm::pmvnorm(lower=c(-crit[1:(nx-1)],x[nx]),upper=c(crit[1:(nx-1)],Inf),sigma=smatrix[1:nx,1:nx])*2+salpha[nx-1]
    pa[nx]=pa[nx-1]-prod(x[1:(nx-1)]<crit[1:(nx-1)])*(1-aa[nx])
  }
  list(pa=pa)
}









weighted.bonferroni=function(pv=c(0.01,0.01,0.02),wv=c(0.5,0.4,0.1),alpha=0.025){
  if(!all(pv>=0)) stop('p-value must be non-negative')
  if(!all(wv>=0)) stop('weight must be non-negative')
  if(sum(wv)>1) {wv <- wv/sum(wv)}
  adj.p=min(pv/wv)
  rej=(adj.p<=alpha)*1
  list(adj.p=adj.p,rej=rej,pv=pv,wv=wv,alpha=alpha)    
}


weighted.simes.BH=function(pv=c(0.01,0.01,0.02),wv=c(0.5,0.4,0.1),alpha=0.025){
  #Benjamini & Hochberg (1997)
  if(!all(pv>=0)) stop('p-value must be non-negative')
  if(!all(wv>=0)) stop('weight must be non-negative')
  if(sum(wv)>1) {wv <- wv/sum(wv)}
  pdata=data.frame(p=pv,w=wv)
  pdata.sort=pdata[with(pdata, order(p, w)), ]
  pdata.sort$cw=cumsum(pdata.sort$w)
  pdata.sort$q=pdata.sort$p/pdata.sort$cw
  adj.p=min(pdata.sort$q)
  rej=(adj.p<=alpha)*1
  list(adj.p=adj.p,rej=rej,pv=pv,wv=wv,alpha=alpha,pdata.sort=pdata.sort)    
}

weighted.simes.HL=function(pv=c(0.01,0.01,0.02),wv=c(0.5,0.4,0.1),alpha=0.025){
  #Hochberg & Liberman (1993)
  if(!all(pv>=0)) stop('p-value must be non-negative')
  if(!all(wv>=0)) stop('weight must be non-negative')
  if(sum(wv)!=1) stop('sum of weights must be one')
  pdata=data.frame(p=pv,w=wv,q=pv/wv)
  pdata.sort=pdata[with(pdata, order(q, p, w)), ]
  pdata.sort$qj=pdata.sort$q/seq(1,length(pdata.sort$q),by=1)
  adj.p=min(pdata.sort$qj)
  rej=(adj.p<=alpha)*1
  list(adj.p=adj.p,rej=rej,pv=pv,wv=wv,alpha=alpha,pdata.sort=pdata.sort)    
}  
  
  
