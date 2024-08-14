## SIR-Hawkes: Linking Epidemic Models and Hawkes Processes to Model Diffusions in Finite Populations
setwd("/Users/sarahmasri/Desktop/Research/MSc/SIR-Hawkes-Comparison/")


## 1: Preliminary (required packages)
library(parallel)
source('scripts/functions/functions-SIR-HawkesN-Rizoiu2018.R')
source('scripts/functions/functions-size-distribution-Rizoiu2018.R')





## 2:  Stochastic R simulation
## simulate 20 stochastic SIR relizations
params.S <- c(N = 1300, I.0 = 300, gamma = 0.2, beta = 1)
nsim <- 20
simdat <- replicate(
  n = nsim,
  generate.stochastic.sir(params = params.S, Tmax = 11, hide.output = T)    
)

# look at simulated data
as.data.frame(simdat[,1])[1:20,]





## 3: Fit stochastic SIR on simulated cascades
## fit stochastic SIR by 
## - choosing a starting point for all parameters, 
## - apply LBFGS algorithm for optimizing the likelihood function of SIR model

# initial fitting point for each execution
params.fit.start <- c(N = 0.1, I.0 = 0.1, gamma = 0.1, beta = 0.1)

.cl <- makeCluster(spec = min(nsim, detectCores()), type = 'FORK')
results <- parSapply(cl = .cl, X = 1:nsim, FUN = function(i) {
  mysim <- as.data.frame(simdat[, i])
  return(fit.stochastic.sir(mysim, params.fit.start))
})
stopCluster(.cl)

# reconstruct result data format
res <- as.data.frame(results[1,])
names(res) <- 1:nsim             
res <- as.data.frame(t(res))     
res$ll <- unlist(results[2,])
complete_res <- res

# let's see how well parameters were retreived
prnt <- rbind(params.S[c('N', 'I.0', 'gamma', 'beta')], 
              apply(X = complete_res[, c('N', 'I.0', 'gamma', 'beta')], MARGIN = 2, FUN = median),
              apply(X = complete_res[, c('N', 'I.0', 'gamma', 'beta')], MARGIN = 2, FUN = sd))
rownames(prnt) <- c('theoretical', 'median', 'sd')
print(prnt[, c('N', 'I.0', 'gamma', 'beta')], digits = 2)





## 4: Fit HawkesN on simulated cascades
## pull out the infective events for HawkesN model

# get the means at given time points, to be able to compare to deterministic
simhistory <- sapply(X = 1:nsim, FUN = function(i) {
  history.S <- SIR2HAWKES.stochastic.event.series(state = simdat[,i])  
})

# model HawkesN on the simulated data following the same steps as SIR

# start point 
params.init <- list(K = 1, c = 0.1, theta = 0.1)





lambda <- function(params, t) {
  (1 - Nt/N) * (mu + sum(triggering.kernel(t)))
}


triggering.kernel(t)









