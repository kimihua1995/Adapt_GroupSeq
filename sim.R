library(MASS)
library(nph)
library(tidyverse)
library(gsDesign)
library(survival)
library(writexl)

source("source_combine.R")
source("source_datsim.R")
source("source_proposed_design.R")
source("source_graphical.R")
source("source_other.R")


#############################################################
# calculate the sample size required for a fixed design
nSurv(lambdaC = log(2)/1,
      hr = 1/1.8,
      eta = 0.1,
      T = 3 + 0,
      minfup = 0,
      ratio = 1, alpha = .025*2/3, beta = .05, sided = 1)
nSurv(lambdaC = log(2)/2.0, # does not affect the results
      hr = 2/3.8,  # by simulation
      eta = 0.1,
      T = 3 + 1,
      minfup = 1,
      ratio = 1, alpha = .025/3, beta = .05, sided = 1)
#############################################################


#############################################################
# Parameters for data generation
group <- crosslist(list(trt = c(0,1,2), group = c(0,1)))
## Group: 0 for BM-, 1 for BM+
## trt: 0 for T0, 1 for T1, 2 for T2
t.enroll.end <- 3
drop.rate <- 0.1
N <- 200  # sample size per treatment and biomarker subgroup
nPFS <- 166  # required # of PFS for H1
nOS <- 160   # required # of OS for H5
w <- c(sqrt(1/4), sqrt(3/4))
alpha <- 0.025
# The alpha spending is derived based on power spending function:
# alpha*I^3, where I is the information fraction
## For PFS, as the final analysis is conducted on IA3 instead of FA,
## the alpha should spend 100% at IA3., i.e., I = 1 at IA3.
salpha_PFS <- matrix(c(0.75^3, 1, 1.01, ## PFS: T2 vs. T0 in BM+
                       0.75^3, 1, 1.01, ## PFS: T1 vs. T0 in BM+ 
                       0.75^3, 1, 1.01, ## PFS: T2 vs. T0 in overall 
                       0.75^3, 1, 1.01), ## PFS: T1 vs. T0 in overall
                       nrow = 4, ncol = 3, byrow = TRUE) * alpha

## For OS, the information is derived based on simulation,
## see details in "information_derivation.R"
salpha_OS <- matrix(c(0.1, 0.25, 1, ## OS: T2 vs. T0 in BM+
                       0.1, 0.25, 1, ## OS: T1 vs. T0 in BM+ 
                       0.1, 0.25, 1, ## OS: T2 vs. T0 in overall 
                       0.1, 0.25, 1), ## OS: T1 vs. T0 in overall
                     nrow = 4, ncol = 3, byrow = TRUE) * alpha

## Weight in the graph
W <- c(2/3,0,0,0,1/3,0,0,0)
## Transition weight matrix in the graph
G <- matrix(c(0,1/3,2/3,0,0,0,0,0,
              0,0,1/2,1/2,0,0,0,0,
              0,1/2,0,1/2,0,0,0,0,
              0.5,0,0,0,0.5,0,0,0,
              0,0,0,0,0,1/3,2/3,0,
              0,0,0,0,0,0,1/2,1/2,
              0,0,0,0,0,1/2,0,1/2,
              1,0,0,0,0,0,0,0),
            nrow = 8, ncol = 8, byrow = T)


B <- 5000 # number of simulations
# Parameters for ORR
betas <- c(0.2,0.2)
q.m <- c(0.2,0.25,0.3,0.2,0.35,0.5)



# Scenario 1
## Parameters for PFS and OS distributions, same for other scenarios
lambda2.m <- log(2)/c(1.0,1.1,1.2,1.0,1.6,1.8)
lambda3.m <- log(2)/c(2.0,2.2,2.4,2.0,3.4,3.8)
dat.list1 <- dat.prepare(333, B, N, q.m, lambda2.m, lambda3.m, betas, group,
                        t.enroll.end, drop.rate)
res1 <- sim.design(dat.list1, w, salpha_PFS, salpha_OS,
                  B, nPFS, nOS, W, G, alpha)


# Scenario 2, null H1 & H5
lambda2.m <- log(2)/c(1.0,1.1,1.2,1.0,1.6,1.0)
lambda3.m <- log(2)/c(2.0,2.2,2.4,2.0,3.4,2.0)
dat.list2 <- dat.prepare(333, B, N, q.m, lambda2.m, lambda3.m, betas, group,
                         t.enroll.end, drop.rate)
res2 <- sim.design(dat.list2, w, salpha_PFS, salpha_OS,
                   B, nPFS, nOS, W, G, alpha)


# Scenario 3, null H2 & H6
lambda2.m <- log(2)/c(1.0,1.1,1.2,1.0,1.0,1.8)
lambda3.m <- log(2)/c(2.0,2.2,2.4,2.0,2.0,3.8)
dat.list3 <- dat.prepare(333, B, N, q.m, lambda2.m, lambda3.m, betas, group,
                         t.enroll.end, drop.rate)
res3 <- sim.design(dat.list3, w, salpha_PFS, salpha_OS,
                   B, nPFS, nOS, W, G, alpha)


# Scenario 4, null H1, H2, H5 & H6
lambda2.m <- log(2)/c(1.0,1.1,1.2,1.0,1.0,1.0)
lambda3.m <- log(2)/c(2.0,2.2,2.4,2.0,2.0,2.0)
dat.list4 <- dat.prepare(333, B, N, q.m, lambda2.m, lambda3.m, betas, group,
                         t.enroll.end, drop.rate)
res4 <- sim.design(dat.list4, w, salpha_PFS, salpha_OS,
                   B, nPFS, nOS, W, G, alpha)


# Scenario 5, null H1, H3, H5 & H7
lambda2.m <- log(2)/c(1.0,1.1,1.0,1.0,1.6,1.0)
lambda3.m <- log(2)/c(2.0,2.2,2.0,2.0,3.4,2.0)
dat.list5 <- dat.prepare(333, B, N, q.m, lambda2.m, lambda3.m, betas, group,
                         t.enroll.end, drop.rate)
res5 <- sim.design(dat.list5, w, salpha_PFS, salpha_OS,
                   B, nPFS, nOS, W, G, alpha)

# Scenario 6, null H2, H4, H6 & H8
lambda2.m <- log(2)/c(1.0,1.0,1.2,1.0,1.0,1.8)
lambda3.m <- log(2)/c(2.0,2.0,2.4,2.0,2.0,3.8)
dat.list6 <- dat.prepare(333, B, N, q.m, lambda2.m, lambda3.m, betas, group,
                         t.enroll.end, drop.rate)
res6 <- sim.design(dat.list6, w, salpha_PFS, salpha_OS,
                   B, nPFS, nOS, W, G, alpha)

# Scenario 7, null H1 to H8
lambda2.m <- log(2)/c(1.0,1.0,1.0,1.0,1.0,1.0)
lambda3.m <- log(2)/c(2.0,2.0,2.0,2.0,2.0,2.0)
dat.list7 <- dat.prepare(333, B, N, q.m, lambda2.m, lambda3.m, betas, group,
                         t.enroll.end, drop.rate)
res7 <- sim.design(dat.list7, w, salpha_PFS, salpha_OS,
                   B, nPFS, nOS, W, G, alpha)



custom_colSums <- function(x, na.rm = FALSE) {
  apply(x, 2, function(col) {
    if (na.rm && all(is.na(col))) {
      return(NA)
    }
    return(sum(col, na.rm = na.rm))
  })
}



rej.prob <- function(res,hset){
  res.S.FA <- res.S.FA2 <- res$rej1[,3,]
  res.B.FA <- res.B.FA2 <- res$rej2[,3,]
  res.G.FA <- res.G.FA2 <- res$rej3[,3,]
  adapt <- res$adapt
  res.S.FA2[c(2,4),adapt == 2] <- NA
  res.S.FA2[c(3,4),adapt == 3] <- NA
  res.S.FA2[,adapt == 4] <- NA
  res.B.FA2[c(2,4),adapt == 2] <- NA
  res.B.FA2[c(3,4),adapt == 3] <- NA
  res.B.FA2[,adapt == 4] <- NA
  res.G.FA2[c(2,4),adapt == 2] <- NA
  res.G.FA2[c(3,4),adapt == 3] <- NA
  res.G.FA2[,adapt == 4] <- NA
  
  if (is.null(hset)){
    prob1.S <- c(rowMeans(res.S.FA[1:8,]),0)
    prob1.B <- c(rowMeans(res.B.FA[1:8,]),0)
    prob1.G <- c(rowMeans(res.G.FA[1:8,]),0)
    prob2.S <- c(rowMeans(res.S.FA2[1:8,], na.rm = T),0)
    prob2.B <- c(rowMeans(res.B.FA2[1:8,], na.rm = T),0)
    prob2.G <- c(rowMeans(res.G.FA2[1:8,], na.rm = T),0)
  }else{
    prob1.S <- c(rowMeans(res.S.FA[1:8,]),
             mean(colSums(res.S.FA[hset,]) > 0))
    prob1.B <- c(rowMeans(res.B.FA[1:8,]),
               mean(colSums(res.B.FA[hset,]) > 0))
    prob1.G <- c(rowMeans(res.G.FA[1:8,]),
               mean(colSums(res.G.FA[hset,]) > 0))
    prob2.S <- c(rowMeans(res.S.FA2[1:8,], na.rm = T),
             mean(custom_colSums(res.S.FA2[hset,], na.rm = T) > 0, na.rm = T))
    prob2.B <- c(rowMeans(res.B.FA2[1:8,], na.rm = T),
               mean(custom_colSums(res.B.FA2[hset,], na.rm = T) > 0, na.rm = T))
    prob2.G <- c(rowMeans(res.G.FA2[1:8,], na.rm = T),
               mean(custom_colSums(res.G.FA2[hset,], na.rm = T) > 0, na.rm = T))
  }
  
  return(list(unadj.S = paste0(round(prob1.S*100,1),"%"), 
              adj.S = paste0(round(prob2.S*100,1),"%"),
              unadj.B = paste0(round(prob1.B*100,1),"%"), 
              adj.B = paste0(round(prob2.B*100,1),"%"),
              unadj.G = paste0(round(prob1.G*100,1),"%"), 
              adj.G = paste0(round(prob2.G*100,1),"%")))
}



write_xlsx(list(unadj_S = as.data.frame(cbind(rej.prob(res1,NULL)$unadj.S,
                               rej.prob(res2,c(1,5))$unadj.S,
                               rej.prob(res3,c(2,6))$unadj.S,
                               rej.prob(res4,c(1,2,5,6))$unadj.S,
                               rej.prob(res5,c(1,3,5,7))$unadj.S,
                               rej.prob(res6,c(2,4,6,8))$unadj.S,
                               rej.prob(res7,1:8)$unadj.S)),
                unadj_B = as.data.frame(cbind(rej.prob(res1,NULL)$unadj.B,
                                              rej.prob(res2,c(1,5))$unadj.B,
                                              rej.prob(res3,c(2,6))$unadj.B,
                                              rej.prob(res4,c(1,2,5,6))$unadj.B,
                                              rej.prob(res5,c(1,3,5,7))$unadj.B,
                                              rej.prob(res6,c(2,4,6,8))$unadj.B,
                                              rej.prob(res7,1:8)$unadj.B)),
                unadj_G = as.data.frame(cbind(rej.prob(res1,NULL)$unadj.G,
                                              rej.prob(res2,c(1,5))$unadj.G,
                                              rej.prob(res3,c(2,6))$unadj.G,
                                              rej.prob(res4,c(1,2,5,6))$unadj.G,
                                              rej.prob(res5,c(1,3,5,7))$unadj.G,
                                              rej.prob(res6,c(2,4,6,8))$unadj.G,
                                              rej.prob(res7,1:8)$unadj.G)),
                adj_S = as.data.frame(cbind(rej.prob(res1,NULL)$adj.S,
                                    rej.prob(res2,c(1,5))$adj.S,
                                    rej.prob(res3,c(2,6))$adj.S,
                                    rej.prob(res4,c(1,2,5,6))$adj.S,
                                    rej.prob(res5,c(1,3,5,7))$adj.S,
                                    rej.prob(res6,c(2,4,6,8))$adj.S,
                                    rej.prob(res7,1:8)$adj.S)),
                adj_B = as.data.frame(cbind(rej.prob(res1,NULL)$adj.B,
                                            rej.prob(res2,c(1,5))$adj.B,
                                            rej.prob(res3,c(2,6))$adj.B,
                                            rej.prob(res4,c(1,2,5,6))$adj.B,
                                            rej.prob(res5,c(1,3,5,7))$adj.B,
                                            rej.prob(res6,c(2,4,6,8))$adj.B,
                                            rej.prob(res7,1:8)$adj.B)),
                adj_G = as.data.frame(cbind(rej.prob(res1,NULL)$adj.G,
                                            rej.prob(res2,c(1,5))$adj.G,
                                            rej.prob(res3,c(2,6))$adj.G,
                                            rej.prob(res4,c(1,2,5,6))$adj.G,
                                            rej.prob(res5,c(1,3,5,7))$adj.G,
                                            rej.prob(res6,c(2,4,6,8))$adj.G,
                                            rej.prob(res7,1:8)$adj.G)),
                as.data.frame(rbind(paste0(table(res1$adapt)/B*100,"%"),
                                    paste0(table(res2$adapt)/B*100,"%"),
                                    paste0(table(res3$adapt)/B*100,"%"),
                                    paste0(table(res4$adapt)/B*100,"%"),
                                    paste0(table(res5$adapt)/B*100,"%"),
                                    paste0(table(res6$adapt)/B*100,"%"),
                                    paste0(table(res7$adapt)/B*100,"%")))
                ), "rejection probability.xlsx")