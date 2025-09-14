# generate outcomes (Y1,Y2,Y3) per treatment/subgroup
## Y1: ORR, Y2: PFS, Y3: OS
y.sim <- function(N,p1,lambda2,lambda3,betas){
  # p1: response rate for ORR
  # lambda2: event rate in the exponential distribution for PFS
  # lambda3: event rate in the exponential distribution for OS
  u3 <- runif(N); x3 <- -log(u3)/lambda3
  u2 <- runif(N); x2 <- -log(u2)/(lambda2 - lambda3)
  y3 <- x3 # time for OS
  y2 <- pmin(x2, x3) # time for PFS
  prob <- 1/(1+exp(-p1-betas[1]*log(y2)-betas[2]*log(y3)))
  #y1 <- exp(p1 + betas[1]*log(y2) + betas[2]*log(y3) + rnorm(N,0,0.1))
  y1 <- rbinom(N, 1, prob)
  y <- data.frame(y1 = y1, y2 = y2, y3 = y3)
  return(y)
}



# generate data
# assume the number of total enrolled patients is equivalent (i.e. N) for each treatment & subgroup
# assume each treatment & subgroup has unique p1, lambda2, and lambda3
# for "group" matrix:
## T2, T1, T0
## group = 1: biomarker-positive, group = 0: biomarker-negative
dat.sim <- function(N, p1.m, lambda2.m, lambda3.m, betas, group,
                    t.enroll.end, drop.rate){
  # t.enroll.end: enrollment end time
  # drop.rate: dropout rate for loss of follow-up time
  J <- nrow(group)
  y.m <- NULL
  for (j in 1:J){
    yj <- y.sim(N, p1.m[j], lambda2.m[j], lambda3.m[j], betas)
    yj$trt <- group$trt[j]
    yj$group <- group$group[j]
    y.m <- rbind(y.m, yj)
  }
  y.m$t.enroll <- runif(J*N, 0, t.enroll.end)
  y.m$t.loss <- rexp(J*N, drop.rate)
  
  return(y.m)
}




# Data modification given IA time
# Generate censoring time and indicator
dat.modify <- function(dat, t.IA){
  dat$t.censor <- pmin(dat$t.loss, t.IA - dat$t.enroll)
  dat$t2 <- pmin(dat$y2, dat$t.censor)
  dat$status2 <- as.numeric(dat$y2 <= dat$t.censor)
  dat$t3 <- pmin(dat$y3, dat$t.censor)
  dat$status3 <- as.numeric(dat$y3 <= dat$t.censor)
  return(dat)
}





dat.prepare <- function(seed, B, N, q.m, lambda2.m, lambda3.m, betas, group,
                        t.enroll.end, drop.rate){
  set.seed(seed)
  
  p1.m <- rep(NA, length(q.m))
  for (i in 1:length(q.m)){
    p1.m[i] <- log(q.m[i]/(1-q.m[i])) - 
      betas[1]*log(log(2)/lambda2.m[i]) - 
      betas[2]*log(log(2)/lambda3.m[i])
  }
  
  dat.list <- list()
  for (i in 1:B){
    dat.list[[i]] <- dat.sim(N, p1.m, lambda2.m, lambda3.m, betas, group, t.enroll.end, drop.rate)
  }
  
  return(dat.list)
}


