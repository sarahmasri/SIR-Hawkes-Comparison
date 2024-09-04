## SIR-Hawkes: Linking Epidemic Models and Hawkes Processes to Model Diffusions in Finite Populations
setwd("/Users/sarahmasri/Desktop/Research/MSc/SIR-Hawkes-Comparison/")


## 1: Preliminary (required packages)
library(parallel)
source('scripts/functions/functions-SIR-HawkesN-Rizoiu2018.R')
source('scripts/functions/functions-size-distribution-Rizoiu2018.R')

N <- getN()



## 2:  Stochastic R simulation
## simulate 20 stochastic SIR relizations
params.S <- c(N = N , I.0 = 300, gamma = 0.2, beta = 1)
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
params.fit.start <- c(I.0 = 0.1, gamma = 0.1, beta = 0.1)

.cl <- makeCluster(spec = min(nsim, detectCores()), type = 'FORK')
results <- parSapply(cl = .cl, X = 1:nsim, FUN = function(i) {
  mysim <- as.data.frame(simdat[, i])
  return(fit.stochastic.sir(mysim, params.fit.start, N))
})
stopCluster(.cl)

# reconstruct result data format
res <- as.data.frame(results[1,])
names(res) <- 1:nsim             
res <- as.data.frame(t(res))     
res$ll <- unlist(results[2,])
complete_res <- res

# let's see how well parameters were retrieved
prnt <- rbind(params.S[c('I.0', 'gamma', 'beta')], 
              apply(X = complete_res[, c('I.0', 'gamma', 'beta')], MARGIN = 2, FUN = median),
              apply(X = complete_res[, c('I.0', 'gamma', 'beta')], MARGIN = 2, FUN = sd))
rownames(prnt) <- c('theoretical', 'median', 'sd')
print(prnt[, c('I.0', 'gamma', 'beta')], digits = 2)



















 
## 4: Fit HawkesN on simulated cascades
## pull out the infective events for HawkesN model

# get the means at given time points, to be able to compare to deterministic
simhistory <- sapply(X = 1:nsim, FUN = function(i) {
  history.S <- SIR2HAWKES.stochastic.event.series(state = simdat[,i])  
})

# model HawkesN on the simulated data following the same steps as SIR

# start point 
params.fit.start <- c(K = 0.1, c = 0.1, theta = 0.1)

# fit the event series with HawkesN
.cl <- makeCluster(spec = min(nsim, detectCores()), type = 'FORK')
results <- parSapply(cl = .cl, X = 1:nsim, FUN = function(i) {
  history.S <- as.data.frame(simhistory[,i])
  fitted.model <- fitSeries(history = history.S, params.fit.start, N)
})
stopCluster(.cl)


#res <- as.data.frame(results['par',])
res <- as.data.frame(matrix(nrow = length(results), ncol = length(results[[1]]$par)))


## TODO: I think there's a more efficient way to do this
for(i in 1:length(results)) {
  if(sum(is.na(results[[i]])) == 0) {
    res[i,] <- results[[i]]$par
  }
}

#names(res) <- 1:nsim
names(res) <- names(results[[1]]$par)
#res <- data.frame(t(res))
 
# these are the theoretical parameters
params.H <- c(K = 5, c = 0.001, theta = 0.2)
 
prnt <- rbind(params.H, 
              apply(X = res, MARGIN = 2, FUN = median, na.rm = T),
              apply(X = res, MARGIN = 2, FUN = sd, na.rm = T))
rownames(prnt) <- c('theoretical', 'median', 'sd')
print(prnt[, c('K', 'theta', 'c')], digits = 2)


res$gamma <- res$theta
res$beta <- res$K * res$theta
prnt <- rbind(params.S[c('gamma', 'beta')], 
              apply(X = res[, c('gamma', 'beta')], MARGIN = 2, FUN = mean, na.rm = T),
              apply(X = res[, c('gamma', 'beta')], MARGIN = 2, FUN = sd, na.rm = T))
rownames(prnt) <- c('theoretical', 'median', 'sd')
print(prnt, digits = 2)






## 5: Plotting size distribution

# theoretical parameters shown previously
params.S <- c(N = 1300, I.0 = 300, gamma = 0.2, beta = 1)
params.H <- c(K = 5, c = 0.001, theta = 0.2, N = 1300)

.transition <- construct.transition.matrix.and.states(params = params.S)

seen_perc <- 0.5
history <- as.data.frame(simhistory[,1])
seenEvents <- round(nrow(history) * seen_perc)

# compute size probs at event 1 and current event
size.est.at.zero <- get.size.distribution(params = params.H, .transition = .transition)
size.est.at.end <- get.size.distribution(params = params.H, .transition = .transition, history = history[seq(seenEvents),])

# plot our cascade
matplot(cbind(size.est.at.zero$final.state, size.est.at.end$final.state), 
        col = c('black', 'blue'), type = 'l', log = 'y', lty = c(1, 1), lwd = 3,
        xlab = 'Cascade final size', ylab = 'Probability',
        main = sprintf('Probability distribution of cascade final size\nfitted on %d seen events (N = %.2f)\n(seen %%: %.2f)', 
                       seenEvents, params.H['N'], seen_perc) )
abline(v = seenEvents, lty = 3, col = 'gray40')
abline(v = size.est.at.end$theo.mean, lty = 2, col = 'darkmagenta')
abline(v = nrow(history), lty = 1, col = 'red')
legend('bottomleft', 
       legend = c('Stoch. size distribution (apriori)', 
                  'Stoch. size distribution (aposteriori)',
                  sprintf('Observed events (%d)', seenEvents), 
                  sprintf('Deterministic size (%.2f)', size.est.at.zero$theo.mean), 
                  sprintf('Observed final size (%d)', nrow(history)) ), 
       lty = c(1, 1, 3, 2, 1), lwd = c(3, 3, 1, 1, 1), col = c('black', 'blue', 'gray40', 'darkmagenta', 'red'), bty = 'n')







## 6: Plotting Log Likelihood

gamma.vals = seq(0.05, 1, by=0.01)
beta.vals = seq(0.2, 5, by=0.05)
K.vals = seq(0.5, 10, by=0.1)
theta.vals = seq(0.01, 1, by=0.01)

N=getN()

plot.ll(model="SIR", param="gamma", values=gamma.vals, data=as.data.frame(simdat[, 1]), N=N)
plot.ll(model="SIR", param="beta", values=beta.vals, data=as.data.frame(simdat[, 1]), N=N)
plot.ll(model="HawkesN", param="K", values=K.vals, data=as.data.frame(simhistory[,1]), N=N)
plot.ll(model="HawkesN", param="theta", values=theta.vals, data=as.data.frame(simhistory[,1]), N=N)


gamma.vals.contour = seq(1/100, 1, by=1/100)
beta.vals.contour = seq(5/100, 5, by=5/100)
K.vals.contour = seq(10/100, 10, by=10/100)
theta.vals.contour = seq(1/100, 1, by=1/100)

contour.plot.ll(model="SIR", param.x="gamma", param.y="beta", 
                values.x=gamma.vals.contour, values.y=beta.vals.contour, data=as.data.frame(simdat[, 1]), N=N)
contour.plot.ll(model="HawkesN", param.x="K", param.y="theta", 
                values.x=K.vals.contour, values.y=theta.vals.contour, data=as.data.frame(simhistory[, 1]), N=N)

