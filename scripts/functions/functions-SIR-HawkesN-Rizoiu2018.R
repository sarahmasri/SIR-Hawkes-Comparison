# this is a library with the methods for the SIR and the HawkesN models
library(nloptr)
library(ggplot2)
library(plotly)

##################### GET HELPERS ######################

## helper function to get N
getN <- function() {1300}

get.theoretical.params <- function(model) {
  if(model == "SIR") {
    params = c(I.0 = 300, gamma = 0.2, beta = 1)
  }
  
  if(model == "HawkesN") {
    params = c(K = 5, c = 0.001, theta = 0.2)
  }
  
  return(params)
}




##################### SIMULATING A STOCHASTIC SIR MODEL ######################

#' Generates the stochastic Continuous Time Markov Chain (discreete states) SIR
#' model. In the stochastic version, time is not given in moments at which the
#' populations are to be computed, it evolves by itself. Tmax gives the end time
#' at which simulation halts.
generate.stochastic.sir <- function(params = c(N = 1300, I.0 = 300, gamma = 0.2, beta = 1),  Tmax = 10, hide.output = F) {
  params["b"] <- 0 ## b is the birth = death rate. In a fixed population it is fixed to zero, but I'm puting it here for later dev.
  
  ## start at time 0
  state <- data.frame(t = 0, I = params["I.0"], S = params["N"] - params["I.0"], C = params["I.0"])
  rownames(state) <- NULL
  j <- 1
  while (state$I[j]>0 & state$t[j]<Tmax) {
    ## first, compute next state rates
    ps1 <- params["beta"] * state$I[j] * state$S[j] / params["N"] ## new infection
    ps2 <- params["gamma"] * state$I[j] ## new recovery
    ps3 <- params["b"] * state$I[j] ## death in infected AND birth from infected
    ps4 <- params["b"] * ( params["N"] - state$I[j] - state$S[j]) ## birth from recovered
    stay <- 1 - ps1 - ps2 - ps3 - ps4 ## staying in the same state -> will determine time of next event
    
    ## get time of next event
    u1 <- runif(n = 1, min = 0, max = 1) # uniform random number
    a <-  ps1 + ps2 + ps3 + ps4
    newTime <- state$t[j] - log(u1) / a
    
    ## draw the next state. Each of the next states has a chance proportional to its ps1, 2, 3, 4.
    ## they need to amount to 1 (since there is no chance of staying in the current state -- accounted by the next time here above)
    u2 <- runif(n = 1, min = 0, max = 1)  # uniform random number
    
    if (u2 <= ps1/a ) { ## new infection -> S-1, I+1
      nextState <- data.frame(t = newTime, I = state$I[j]+1, S = state$S[j]-1, C = state$C[j]+1)
    } else if (u2 <= (ps1 + ps2)/a) {  ## new recovery -> S, I-1
      nextState <- data.frame(t = newTime, I = state$I[j]-1, S = state$S[j], C = state$C[j])
    } else if (u2 <- (ps1 + ps2 + ps3)/a) {  ## death in infected AND birth from infected -> S+1, I-1
      nextState <- data.frame(t = newTime, I = state$I[j]-1, S = state$S[j]+1, C = state$C[j])
    } else { ## birth from recovered -> S+1, I
      nextState <- data.frame(t = newTime, I = state$I[j], S = state$S[j]+1, C = state$C[j])
    }
    state <- rbind(state, nextState)
    j <- j + 1
    if (!hide.output) {
        cat(sprintf("\rCurrent simulation time: %.3f / %.3f (S=%d, I=%d, R=%d, C=%d).", state$t[j], Tmax, state$S[j], state$I[j], params["N"]-state$S[j]-state$I[j], state$C[j]))
    }
  }
  
  rownames(state) <- NULL
  state$R <- state$C - state$I
  names(state) <- c("time", "I", "S", "C", "R")
  state <- state[c("time", "S", "I", "R", "C")]
  
  if (!hide.output) {
      cat(sprintf("\n--> Simulation done!\n"))
  }
  return(state)
}

#' From a stochastic state data.frame (generated using generate.stochastic.sir),
#' this function extracts the event series. For compatibility with Hawkes, it
#' also creates dummy event magnitudes of 1.1.
SIR2HAWKES.stochastic.event.series <- function(state) {
  state <- as.data.frame(state)
  ## generate as many events with time 0 as there are I.0 (infections at time 0)
  history <- as.data.frame(matrix(data = rep(x = c(1.1, 0), times = state$I[1]), nrow = state$I[1], ncol = 2, byrow = T))
  colnames(history) <- c("magnitude", "time")
  
  ## compute where we had a new infection event
  state$new <- F
  state$new[-1] <- state$C[-1] > state$C[-nrow(state)]
  evnt <- data.frame(magnitude = 1.1, time = state$time[state$new])
  
  ## concatenate the two
  history <- rbind(history, evnt)
  return(history)
}

############################# FITTING A STOCHASTIC SIR MODEL ################################


# fit SIR using loglikelhood maximization
fit.stochastic.sir <- function(mysim, params.start, N, iterations = 5000) {
  # if (params.start["N"] <  max(mysim[, -1])) {
  #   params.start["N"] <-  max(mysim[, -1])
  # }
  
  params.start["I.0"] <- mysim$I[1]
  # res <- lbfgs(x0 = unlist(params.start),                            
  #              fn = stochastic.sir.complete.neg.log.likelihood,          
  #              lower =  c(I.0 = mysim$I[1], gamma = 0, beta = 0),
  #              upper = c(I.0 = mysim$I[1], gamma = Inf, beta = Inf),
  #              state = mysim,
  #              N = N)
  
  
  
  ## Error when useing L-BFGS-B:
  ## Error in optim(par = unlist(params.start), fn = stochastic.sir.complete.neg.log.likelihood, : non-finite finite-difference value [1]
  
  res <- optim(par = unlist(params.start),
               fn = stochastic.sir.complete.neg.log.likelihood,
              # method = "L-BFGS-B",
               method = "Nelder-Mead",
               hessian = FALSE,
              # lower = c(I.0 = mysim$I[1], gamma = 0, beta = 0),
              # upper = c(I.0 = mysim$I[1], gamma = Inf, beta = Inf),
               state = mysim,
               control = list(maxit = iterations, factr = "1e-8"),
                #list(trace = 3,
              #                maxit = 1000,
              #               ndeps = c(1e-4, 1e-4, 1e-4)),
               N = N)
  
  names(res$par) <- names(params.start)
  
  return(res)
}

#' Complete log-likelihood function for a stochastic SIR model, computed for a
#' set of params. It assumes that the state matrix is in the SIR format (with
#' "time", "I", "S", "C", "R" headers) and it contains all events in the
#' process: new infection, new recovery and (optionally) birth and death. Such a
#' matrix can be obtained using the "generate.stochastic.sir" function.
stochastic.sir.complete.neg.log.likelihood <- function(params, state, N) {
  state <- as.data.frame(state)
  names(params) <- c("I.0", "gamma", "beta")
  params <- unlist(params)
  params["b"] <- 0 ## b is the birth = death rate. In a fixed population it is fixed to zero, but I'm puting it here for later dev.
  
  ## we will sum here the log-likelihood
  ll <- 0
  
  contribs <- c(0, 0, 0, 0, 0); last_infection_time <- 0; last_recovery_time <- 0
  contrib_Lambda_inf <- 0; contrib_Lambda_rec <- 0
  ## go through events (rows) one by one. Note that the first row is the initial
  ## state (cannot compute any probs there)
  for (j in 2:nrow(state)) {
    ## first, compute next state probabilities/rates
    ps1 <- params["beta"] * state$I[j-1] * state$S[j-1] / N ## new infection 
    ps2 <- params["gamma"] * state$I[j-1] ## new recovery
    ps3 <- params["b"] * state$I[j-1] ## death in infected AND birth from infected
    ps4 <- params["b"] * ( N - state$I[j-1] - state$S[j-1]) ## birth from recovered
    a <-  ps1 + ps2 + ps3 + ps4 ## the inter-event times are exponentially distributed with param a. 
    
    ## the likelihood of having waited as much as we have
    ll <-  sum(ll, dexp(x = (state$time[j] - state$time[j-1]), rate = a, log = T), na.rm = T)
    
    new_event_prob_timing <- a * exp(-(state$time[j] - state$time[j-1])*a )
    contribs[3] <- sum(contribs[3], log(new_event_prob_timing), na.rm = T)
    
    ## compute contributions to Lambda by current component
    contrib_Lambda_inf <- sum(contrib_Lambda_inf, (state$time[j] - state$time[j-1]) * ps1)
    contrib_Lambda_rec <- sum(contrib_Lambda_rec, (state$time[j] - state$time[j-1]) * ps2)
    
    ## and now, how likely to have observed an event such as what we've seen
    if (state$I[j] == state$I[j-1]+1 && state$S[j] == state$S[j-1]-1) { ## new infection -- S-1, I+1
      ll <- sum(ll, log(ps1/a), na.rm = T)
      
      contribs[1] <- sum(contribs[1], log(ps1/a), na.rm = T)
      
      new_inf_prob_timing <- ps1 * exp(-contrib_Lambda_inf)
      contribs[4] <- sum(contribs[4], log(new_inf_prob_timing), na.rm = T)
      contrib_Lambda_inf <- 0
    } else if (state$I[j] == state$I[j-1]-1 && state$S[j] == state$S[j-1]) {  ## new recovery -- S, I-1
      ll <- sum(ll, log(ps2/a), na.rm = T)
      contribs[2] <- sum(contribs[2], log(ps2/a), na.rm = T)

      
      new_rec_prob_timing <- ps2 * exp(-contrib_Lambda_rec)
      contribs[5] <- sum(contribs[5], log(new_rec_prob_timing), na.rm = T)
      contrib_Lambda_rec <- 0
    } else if (state$I[j] == state$I[j-1]-1 && state$S[j] == state$S[j-1]+1) {## death in infected AND birth from infected -- S+1, I-1
      ll <- sum(ll, log(ps3/a), na.rm = T)
    } else if (state$I[j] == state$I[j-1] && state$S[j] == state$S[j-1]+1) { ## birth from recovered -- S+1, I
      ll <- sum(ll, log(ps4/a), na.rm = T)
    }
  }
  
  names(ll) <- "neg.ll"
  return(nll = -ll)
}

########################### FITTING HawkesN #############################





fitSeries <- function(history, params_init, N) {
  model <- NA
  tryCatch({
    # ## optimize with NLOPT
    model <- find.optimal.parameters(history = history, init_params = params_init, N = N)
  }, error = function(err) {
    print(paste("[fitSeries] Error in optim:  ", err))
    model <- list(par = c(K = NA, c = NA, theta = NA), value = NA, iter = 0, convergence = 0, message = 'Error')
  })
  
  # if (sum(is.na(model) | is.na(model$value)) > 0) {
  #   model <- list(par = c(K = NA, c = NA, theta = NA), value = NA, iter = 0, convergence = 0, message = 'Error')
  # }

  return(model)
}

find.optimal.parameters <- function(history, init_params, # = c(K = 0.1, c = 0.1, theta = 0.1), ## Note: uncomment?
                                    lowerBound = c(K = 0, c = .Machine$double.eps, theta = 0), 
                                    upperBound = c(K = Inf, c = 300, theta = Inf),
                                    iterations = 5000,
                                    N = getN(), 
                                    ...) {

  if ( sum(! names(init_params) %in% names(lowerBound)) > 0 ) {
    newParamName <- names(init_params)[! names(init_params) %in% names(lowerBound)]
    warning(sprintf("You gave me the param '%s', but no boundry for it. I'll assume entire real range.", newParamName ))
    lowerBound[newParamName] <- -Inf
    upperBound[newParamName] <- Inf
  }
  lowerBound <- lowerBound[names(init_params)]
  upperBound <- upperBound[names(init_params)]
  
  ## correct initial parameters out of bounds -- should not happen, but...
  init_params[init_params > upperBound] <- upperBound[init_params > upperBound]
  init_params[init_params < lowerBound] <- lowerBound[init_params < lowerBound]
  
  ## if my initial N (population size) is below what I observe, I know it is a
  ## bad starting point. Therefore, I will correct it to start at least where I
  ## see it.
  # if ("N" %in% names(init_params) && init_params["N"] < nrow(history))
  #   init_params["N"] <- nrow(history)
  
  #res <- lbfgs(x0 = unlist(init_params),
  #             fn = neg.log.likelihood,
  #             gr = closedGradient,
  #             lower = lowerBound,
  #             upper = upperBound,
  #             control = list(maxeval = iterations, xtol_rel = "1e-8"),
  #             history = history,  ...)
  
  
  
  res <- optim(par = unlist(init_params),
               fn = neg.log.likelihood,
               gr = closedGradient,
               method = "L-BFGS-B",
              # method = "Nelder-Mead",
              # hessian = FALSE,
               lower = lowerBound,
               upper = upperBound,
               control = list(maxit = iterations, factr = "1e-8"),
               history = history,  
               N = N) 
  
  
  
  
  
  print(res)
  
  names(res$par) <- names(init_params)
  res$value <- neg.log.likelihood(params = res$par, history = history, N) #, ... )
  
  return(res)
}

#' Next function calculates the negative
#' log-likelihood. This is given that the optim function minimizes a function. '
#' Minimizing the -1 * log-likelihood amounts to maximizing the log-likehood. '
#' @param .cl - ugly hack to pass the cluster parameter to the gradient function
#'  in the optimization procedure. Not necessary for this function. '
neg.log.likelihood <- function(params, history, N = getN()) { 
  names(params) <- c("K", "c", "theta")
  params <- as.list(unlist(params))
  
  ## if all good, continue
  bigT <- max(history$time)
  ## startEvent is the first event generated endogeneously (not in the initial pool of events)
  startEvent <- which.max(history$time > 0)
  endEvent <- nrow(history)
  
  ## get the formula of the log-likelihood * -1.
  ## that in the summation, we eliminated the first event, as there is no background rate in our lambda(t)
  lambdaInt <- integrateLambda(lower = 0, upper = bigT, history = history, params = params, N)
  lambdaSums <- sum(log(lambda(t = history$time[startEvent:endEvent], history = history, params = params)))
  retVal <- lambdaInt - lambdaSums 
  
  # if ("N" %in% names(params)) {
  #   ## apply ridge regularizer (with alpha = 10, desired N val is the observed
  #   ## size of the cascade) in the HawkesN case
  #   retVal <- retVal
  # }
  
  ## some optimizers don't like Inf, NA or NaN. Set it to the maximum real number
  if (!is.finite(retVal)) {
    warning(paste("Following params got an error value (", retVal, "): ", paste(params, collapse = ", "), sep = ""))
    retVal <- NA
  }
  
  return(retVal)
}

#' A function to calculate the value of the integral of lambda for 
#' Exponential Kernel with finite population
integrateLambda <- function(lower, upper, history, params, N = getN(), mmin = 1) { 
  names(params) <- c("K", "c", "theta")
  params <- as.list(unlist(params))
  endEvent <- nrow(history)
  
  ## closed form
  res <- sum(sapply(X = 1:nrow(history), FUN = function(i) {
    j <- i:endEvent
    ## make sure that the counter of available events never goes below zero,
    ## otherwise we get a negative integral
    avail_events <- N - i:(endEvent-1)
    avail_events[avail_events < 0] <- NA
    
    terms <- exp(-params$theta * (history$time[j] - history$time[i]) )
    
    val <- sum((avail_events / N) * (terms[-length(terms)] - terms[-1]))
    return(val)
  })) * params$K
  
  return(res)
}

#' Computes the conditional intensity of the non-homogenous Poisson process (
#' lambda(t) ) It is the equivalent of the CIF function in "simulation-model.R",
#' but updated for passing parameters as a list and calculating the intensity
#' for a vector on input (needed for numerical integration).
lambda <-  function(t, history, params = c(K = 0.024, c = 0.001, theta = 0.2), N = getN()) {
  params <- as.list(unlist(params))
  
  res <- sapply(X = t, FUN = function(ti) {
    subst <- history[history$time <= ti,]
    return(sum(kernelFct(event = subst, t = ti, params = params, N = N)))
  })
  
  res[res == 0] <- .Machine$double.eps
  
  return(res)
}

#' calculates the influence of a given event at the givent time the event is a 2
#' elements list (mi, ti). Inclusive parameter means that if an event is present
#' at time t, then it contributes to the conditional intensity. If inclusive ==
#' F, then events at time t are removed.
kernelFct <- function(event, t, params = c(K = 0.024, theta = 0.2), alpha = 2.016, mmin = 1, N = getN()) {
  params <- unlist(params)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]
  Nt <- min(sum(ti <= t), N)

  
  # f(p_j) part - virality of a diffusion. Constant for a given diffusion. Furthermore, discount for available events.
  fun_f <- params["K"] * (1 - Nt / N)
  
  fun_psi <- params["theta"] * (exp(-params["theta"] * (t - ti)))
  
  val = fun_f * fun_psi
  val[t<ti] = 0
  val[mi<mmin] = 0
  
  return(val)
}

#' A function to calculate the dervative in closed form for Finite Exponential Kernel.
closedGradient <- function(params, history, N = getN()){
  names(params) <- c("K", "c", "theta")
  
  ## the wrapping insures that we always return as many values as parameters
  ## compose first the return vector with zero everywhere
  ret <- rep(x = 0, times = length(params))
  names(ret) <- names(params)

  # variables collection
  params <- as.list(unlist(params))
  bigT <- max(history$time)
  n <- nrow(history)
  k <- params$K
  theta <- params$theta
  
  ## number of initial events 
  i_0 <- nrow(history[history$time <= 0,])
  
  # In all the following calcuaion, res is the part of derivative coming 
  # from bigLambda part and first is derivative coming from 
  # summation of log(lambda(ti))
  
  #calculating derivative wrt k in closed form
  res <- sum(sapply(X = 1:n, FUN = function(i) {
    ## compute all the exponentials e^{-\theta (t_j - t_i)}
    j <- i:n
    terms <- exp(-params$theta * (history$time[j] - history$time[i]) )
    
    val <- sum(((N - i:(n-1)) / N) * (terms[-length(terms)] - terms[-1]))
    return(val)
  }))
  derivK <- ((n-i_0) / k) - res
  
  #calculating derivative wrt theta in closed form
  res <- sum(sapply(X = 1:n, FUN = function(i) {
    ## compute all the exponentials e^{-\theta (t_j - t_i)}
    j <- i:n
    terms <- exp(-params$theta * (history$time[j] - history$time[i]) )
    times <- history$time[j] - history$time[i]
    
    val <- sum(((N - i:(n-1)) / N) * ((times[-1] * terms[-1]) - (times[-length(times)] * terms[-length(terms)])))
    return(val)
  })) * k
  
  first <- sum(sapply(X = (i_0+1):(n), FUN = function(i) {
    
    numerator <- sum(exp(-theta * (history$time[i] - history$time[1:i-1])) * 
                       - (history$time[i] - history$time[1:i-1]))
    
    denominator <- sum(exp(-theta * (history$time[i] - history$time[1:i-1])))
    return(numerator / denominator)
  }))
  derivTheta <- ((n-i_0) / theta) + first - res
  
  #calculating derivative wrt N in closed form
  res <- sum(sapply(X = 1:n, FUN = function(i) {
    ## compute all the exponentials e^{-\theta (t_j - t_i)}
    j <- i:n
    terms <- exp(-params$theta * (history$time[j] - history$time[i]) )
    
    val <- sum( ( ( i:(n-1) ) / N^2 ) * ( terms[-length(terms)] - terms[-1] ) )
    return(val)
  })) * k
  
  ## n! = gamma(n+1), therefore deriv of log(n!) is deriv lof log(gamma(n+1)) which is digamma(n+1)
  # first <- digamma(N-i_0+1) - ((n-i_0)/N)
  # first <- sum(sapply(X = (i_0+1):n, FUN = function(i) {return(1/(N-i+1))})) - ((n-i_0)/N)
  # derivN <- first - res
  
  
  
  ## added a -1 because we are minimizing log-likelihood
  retval <- c(K = derivK, theta = derivTheta)
  
  ## put the vals we got into the return array
  ret[names(retval)] <- -1 * retval
  return(ret)
}



########################### PLOTS #############################

#' plot likelihood against a set of parameter values
plot.ll <- function(model, param, values, data, N=getN()) {
  if(!(model %in% c("SIR", "HawkesN")) ) {
    return("Model choice not SIR or HawkesN")
  }
  
  if(!(param %in% c("gamma", "beta", "K", "theta"))) {
    return("Parameter choice not one of: gamma, beta, K, or theta")
  }
  
  n = length(values)
  ll = numeric(n)
  params = get.theoretical.params(model)
  
  if(model == "SIR") {
    for(i in 1:n) {
      params[param] = values[i]
      ll[i] = stochastic.sir.complete.neg.log.likelihood(params=params, state=data, N=N)
    }
  }
  
  if(model == "HawkesN") {
    for(i in 1:n) {
      params[param] = values[i]
      ll[i] =  neg.log.likelihood(params=params, history=data, N = N) 
    }
  }
  
  
  df = data.frame(values, ll)
  

  
  filename = paste0("/Users/sarahmasri/Desktop/Research/MSc/SIR-Hawkes-Comparison/figures/", 
                    model, 
                    "_neg-log-likelihood-plot_",
                    param,
                    ".png"
  )
  

  
  png(
    filename,
    width     = 5,
    height    = 5,
    units     = "in",
    res       = 1200,
    pointsize = 4
  )
  par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
  )
  
  plot <- ggplot(df, aes(x=values, y=ll)) + geom_point() + 
    geom_hline(yintercept = min(ll), linetype = "dashed") +
    #geom_vline(xintercept = get.theoretical.params(model)[param], color = "blue", linewidth = 1.5) +
    geom_vline(xintercept = values[which.min(ll)], color = "firebrick", linewidth = 1) +
    xlab(param) + 
    ylab("Negative Log Likelihood") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  
  print(plot)
 
  dev.off()
  
  
  return(ll)
}


## Create a contour plot of negative log likelihoods for two specified params
contour.plot.ll <- function(model, param.x, param.y, values.x, values.y, data, N=getN()) {
  if(!(model %in% c("SIR", "HawkesN")) ) {
    return("Model choice not SIR or HawkesN")
  }
  
  if(!(param.x %in% c("gamma", "beta", "K", "theta"))) {
    return("x-axis parameter choice not one of: gamma, beta, K, or theta")
  }
  
  if(!(param.y %in% c("gamma", "beta", "K", "theta"))) {
    return("y-axis parameter choice not one of: gamma, beta, K, or theta")
  }
  
  
  
  n.x = length(values.x)
  n.y = length(values.y)
  
  ll = matrix(numeric(n.x * n.y), nrow = n.y, ncol = n.x)  
  params = get.theoretical.params(model)
  
  
  if(model == "SIR") {
    for(i in 1:n.y) {
      params[param.y] = values.y[i]
      for(j in 1:n.x) {
        params[param.x] = values.x[j]
        ll[i,j] = stochastic.sir.complete.neg.log.likelihood(params=params, state=data, N=N)
      }
    }
  }
  
  if(model == "HawkesN") {
    for(i in 1:n.y) {
      params[param.y] = values.y[i]
      for(j in 1:n.x) {
        params[param.x] = values.x[j]
        ll[j,i] = neg.log.likelihood(params=params, history=data, N=N)
      }
    }
  }
  
  
  
  
  filename = paste0("/Users/sarahmasri/Desktop/Research/MSc/SIR-Hawkes-Comparison/figures/", 
                 model, 
                 "_contour-plot_",
                 param.x,
                 "_",
                 param.y,
                 ".png"
                 )
  
  
  png(
    filename,
    width     = 4,
    height    = 4,
    units     = "in",
    res       = 1200,
    pointsize = 4
  )
 
  
  filled.contour(
    x = values.x,
    y = values.y,
    z = ll,
    plot.axes={
      axis(1,cex.axis=1)
      axis(2,cex.axis=1)
    },
    plot.title={
      title(xlab=param.x,cex.lab=1.5)
      mtext(param.y,2,cex=1.5,line=3,las=0)
    },
    key.title = {par(cex.main=0.7);title(main="Negative\n Log\n-Likelihood")},
  )

  
  dev.off()
  
  
  

  
  
  return(ll)
}
