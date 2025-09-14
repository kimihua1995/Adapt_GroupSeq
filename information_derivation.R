source("source_datsim.R")
source("source_other.R")

betas <- c(0.2,0.2)
q.m <- c(0.2,0.25,0.3,0.2,0.35,0.5)
lambda2.m <- log(2)/c(1.0,1.1,1.2,1.0,1.6,1.8)
lambda3.m <- log(2)/c(2.0,2.2,2.4,2.0,3.4,3.8)
group <- crosslist(list(trt = c(0,1,2), group = c(0,1)))
t.enroll.end <- 3
drop.rate <- 0.1
N <- 200
nPFS <- 166
nOS <- 160
B <- 1000


p1.m <- rep(NA, length(q.m))
for (i in 1:length(q.m)){
  p1.m[i] <- log(q.m[i]/(1-q.m[i])) - 
    betas[1]*log(log(2)/lambda2.m[i]) - 
    betas[2]*log(log(2)/lambda3.m[i])
}

# Derive the approximate interim timing by simulation
set.seed(212)
t_IAs.list <- NULL
for (b in 1:B){
  if (!b%%500) print(b)
  dat <- dat.sim(N, p1.m, lambda2.m, lambda3.m, betas, group, t.enroll.end, drop.rate)
  t_IAs <- c(get.IA(dat, c(1), c(0,2), c(ceiling(0.25*nPFS),ceiling(0.75*nPFS), nPFS), "PFS"),
             get.IA(dat, c(1), c(0,2), c(nOS), "OS"))
  t_IAs.list <- rbind(t_IAs.list, t_IAs)
}
colMeans(t_IAs.list) * 12  # 15, 28, 34, 47 months
t_IAs_m <- colMeans(t_IAs.list)


# Simulate a data with very large sample size, estimate the information fraction for each hypothesis
set.seed(212)
dat_L <- dat.sim(20000, p1.m, lambda2.m, lambda3.m, betas, group, t.enroll.end, drop.rate)
frac_PFS <- frac_OS <- matrix(NA, 4, 3)
group.list <- list(c(1), c(1), c(0,1), c(0,1))
trt.list <- list(c(0,2), c(0,1), c(0,2), c(0,1))
for (i in 1:4){
  dat_i <- dat_L[dat_L$trt %in% trt.list[[i]] & dat_L$group %in% group.list[[i]],]
  for (j in 2:length(t_IAs_m)){
    dat_it <- dat.modify(dat_i, t_IAs_m[j])
    frac_PFS[i,j-1] <- sum(dat_it$status2)
    frac_OS[i,j-1] <- sum(dat_it$status3)
  }
}
frac_PFS <- (frac_PFS/frac_PFS[,2])^3
frac_OS <- (frac_OS/frac_OS[,3])^3












