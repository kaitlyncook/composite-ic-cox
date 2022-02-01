#-----------------------------------------------------------------------------------------------------#
#         Poisson-Augmented Estimation of Stratified Cox Models for Interval-Censored Data:           #
#           Functions for Data Generation with Time-Varying Covariates/Covariate Effects              #
#-----------------------------------------------------------------------------------------------------#

# Notes:
#   All included data-generation functions adopt the structure of the Botswana Combination Prevention Project,
#   with a single binary intervention and a time-varying exogenous event indicator. These functions also simulate
#   observations assuming either independence within each cluster or membership within exchangeable sub-networks
#   within each cluster

library(copula)

#-----------------------------------------------------------------------------------------------------
# Data Generation Functions: Random Baseline Survival (Baseline Hazard) Function Generation
#-----------------------------------------------------------------------------------------------------

### Returns list of randomly generated baseline survival and corresponding hazard functions
  # maxT            numeric             upper bound of the support for the failure time distribution
  # knots           numeric             number of internal knots for spline representation of the baseline survival function
return_baseline <- function(maxT, knots, max.grid.width=0.001, min.grid.number=16000){
  
  # Selecting random anchor points for the construction of the baseline functions
  partition <- c(0, sort(runif(knots, 0, maxT)), maxT)
  surv <- c(1, sort(runif(knots), decreasing = TRUE), 0)
  
  # Fitting a monotonic cubic spline model to the anchor points
  survSpline <- stats::splinefun(x=partition, y=surv, method="hyman")
  
  # Piecewise-constant approximations to the baseline survival, density, and hazard funcitons
  timeGrid <- seq(0, maxT, by=min(max.grid.width, maxT/min.grid.number))
  survGrid <- survSpline(timeGrid)
  pdfGrid <- -survSpline(timeGrid, 1)
  hazGrid <- pdfGrid/survGrid
  
  # Quality assurance on the estimated components
  if (sum(survGrid > 1) > 0){
    index <- which(survGrid > 1)
    survGrid[index] <- 1
  }
  if (sum(survGrid < 0) > 0){
    index <- which(survGrid < 0)
    survGrid[index] <- 0
  }
  if (sum(is.na(hazGrid)) > 0){
    index  <- which(is.na(hazGrid))
    hazGrid[index] <- 0
  }
  if (sum(hazGrid < 0) > 0){
    index <- which(hazGrid < 0)
    hazGrid[index] <- 0
  }
  
  return(list(survival=data.frame(time=timeGrid, surv=survGrid), 
              density=data.frame(time=timeGrid, density=pdfGrid),
              hazard=data.frame(time=timeGrid, hazard=hazGrid)))
}

### Returns list of randomly generated baseline survival and corresponding hazard functions (under an exponential form for the baseline hazard)
  # maxT            numeric             upper bound of the support for the failure time distribution
  # lambda0         numeric             baseline hazard
  # alpha           numeric             random effects variance
  # max.grid.width  numeric             maximum size of the intervals in the time grid
  # min.grid.number numeric             minimum number of intervals in the time grid
return_exp_baseline <- function(maxT, lambda0, alpha, max.grid.width=0.001, min.grid.number=16000){
  
  timeGrid <- seq(0, maxT, by=min(max.grid.width, maxT/min.grid.number))
  b <- rnorm(1, 0, sqrt(alpha))
  hazGrid <- rep(lambda0*exp(b), length(timeGrid))
  survGrid <- c(1, exp(-cumsum(hazGrid[1:(length(hazGrid)-1)]*diff(timeGrid))))
  pdfGrid <- hazGrid*survGrid
  
  return(list(survival=data.frame(time=timeGrid, surv=survGrid), 
              density=data.frame(time=timeGrid, density=pdfGrid),
              hazard=data.frame(time=timeGrid, hazard=hazGrid)))
}


#-----------------------------------------------------------------------------------------------------
# Data Generation Functions: Time-Varying Binary Indicator of UTT
#-----------------------------------------------------------------------------------------------------

### Returns individual survival curve
  # baseline.i      list                list of baseline survival, density, and hazard for stratum i
  # tau.i           numeric             internal time of the UTT policy implementation  
  # beta            numeric             coefficient vector
  # x.i             numeric             indicator of intervention assignment
get_surv <- function(baseline.i, tau.i, beta, x.i){
  
  surv.i <- baseline.i$survival
  index <- sum(surv.i$time <= tau.i)
  
  # Calculating survival at times prior to the policy implementation
  s.prior <- surv.i$surv[1:index]^(exp(beta[1]*x.i))
  
  if (index == nrow(surv.i)){
    s.curve <- data.frame(time = surv.i$time, surv=s.prior)
  } else{
    # Calculating survival at times post the policy implementation
    s.tau.i <- surv.i$surv[index]
    s.post <- s.tau.i^(exp(beta[1]*x.i))*s.tau.i^(-exp(beta[1]*x.i + beta[2] + beta[3]*x.i))*surv.i$surv[(index+1):nrow(surv.i)]^(exp(beta[1]*x.i + beta[2] + beta[3]*x.i))
    
    s.curve <- data.frame(time=surv.i$time, surv=c(s.prior, s.post))
  }
  
  # Replacing missing/NaN survival times with zero
  if (sum(is.na(s.curve$surv)) > 0){
    na.index <- which(is.na(s.curve$surv))
    s.curve$surv[na.index] <- 0
  }
  
  return(s.curve)
}

### Returns failure time associated with an individual survival probability
  # s               numeric             individual survival probability
  # tau.i           numeric             internal time of the UTT policy implementation  
  # beta            numeric             coefficient vector
  # x.i             numeric             indicator of intervention assignment
  # baseline.i      list                list of baseline survival, density, and hazard for stratum i
get_inv_surv <- function(s, tau.i, beta, x.i, baseline.i){
  
  # Obtaining piecewise-approximated individual survival function
  s.curve.i <- get_surv(baseline.i, tau.i, beta, x.i)
  
  # Returning the failure time corresponding to s
  index <- which(s.curve.i$surv <= s)[1]
  t <- s.curve.i$time[index]
  return(t)
}

### Generates individual long-form data vector with time-varying UTT indicator
  # u.t             numeric             uniform random variate
  # tau.i           numeric             internal time of the UTT policy implementation  
  # beta            numeric             coefficient vector
  # x.i             numeric             indicator of intervention assignment
  # baseline.i      list                list of baseline survival, density, and hazard for stratum i  
  # visits.i        numeric             internal time of all study-mandated HIV testing dates
  # cens            numeric             loss-to-follow-up rate
  # stratum         numeric             ID number of the asociated stratum 
  # group           numeric             ID number of the associated community
  # id              numeric             ID number of the individual
gen_individual_tvc <- function(u.t, tau.i, beta, x.i, baseline.i, visits.i, cens, stratum, group, id){
  
  # Generating survival and censoring times
  s.time <- get_inv_surv(u.t, tau.i, beta, x.i, baseline.i)
  if (cens==0){
    c.time <- visits.i[length(visits.i)] + 1
  } else{
    u.c <- runif(1)
    c.time <- -log(u.c)/cens
  }
  
  # Deriving final interval censored observation and return (left, right)
  visits.i[visits.i > c.time] <- NA
  event <- as.numeric(visits.i > s.time)
  if (sum(is.na(event))==length(event)){
    out <- data.frame(stratum=stratum, group=group, treat=x.i, left=0, right=Inf)
    return(out)
  }
  event.status <- max(event, na.rm=TRUE)
  if (event.status==1){
    first.event.time <- min(which(event==1))
    l <- ifelse(first.event.time==1, 0, 
                ifelse(sum(!is.na(visits.i[1:(first.event.time-1)]))==0, 0,
                       max(visits.i[1:(first.event.time-1)], na.rm=TRUE)))
    r <- visits.i[first.event.time]
  } else{
    l <- max(visits.i, na.rm=TRUE)
    r <- Inf
  }
  
  # Returning the final long-form data for individual i
  out <- data.frame(stratum=rep(stratum, 2), group=rep(group, 2), id=rep(id, 2), treat=rep(x.i, 2), utt=c(0, 1), 
                    start=c(0, tau.i), stop=c(tau.i, Inf), left=rep(l, 2), right=rep(r, 2))
  return(out)
}

### Generates stratified dataset with a cluster- and time-constant binary variable and time-varying UTT indicator
  # M               numeric             number of communities 
  # ni              numeric             number (range) of individuals in each community 
  # beta            numeric             coefficient vector for (1) the time-constant variable, (2) the time-varying variable, and (3) their interaction
  # gen_visits      function            function for generating internal visit times
  # cens            numeric             loss-to-follow-up rate
  # rho             numeric             Kendall's tau for the dependence within community sub-networks
  # baseline        list                list containing the baseline distribution functions for each stratum
  # S               function            either (a) a function for generating cluster-varying stratum membership or (b) a vector indicating stratum membership for each cluster
gen_time_varying_cov_tx <- function(M, ni, beta, gen_visits, cens, rho=0, baseline, S){
  
  # Generating the design matrix if observations are independent within each cluster
  if (rho==0){
    
    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      n.i <- sample(ni[1]:ni[2], 1)
      u <- runif(n.i, 0, 1)
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }
    
    # Generating the design matrix if observations belong to an exchangeable subnetwork within each cluster
  } else {
    
    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      networks <- rpois(ni[2]/2, 2) + 2
      networks.cm <- cumsum(networks)
      n.i <- sample(ni[1]:ni[2], 1)
      sub.ni <- networks[networks.cm < n.i]
      copula.fn.list <- lapply(sub.ni, FUN=function(x) claytonCopula(-2*rho/(rho-1), dim = x))
      u <- do.call(c, lapply(copula.fn.list, FUN=function(x) as.numeric(rCopula(1, x))))
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }
    
  }
  
  # Determining stratum membership
  if (is.function(S)) {
    data.list.2 <- lapply(data.list, FUN=function(dat){
      n.i <- nrow(dat)
      s <- S(n.i)
      cbind(s, dat)
    })
  } else{
    data.list.2 <- mapply(FUN=function(s, dat){
      cbind(rep(s, nrow(dat)), dat)
    }, s=S, dat=data.list, SIMPLIFY=FALSE)
  }
  
  data <- do.call(rbind, data.list.2)
  n <- nrow(data)
  
  # Generating (internal) inspection times
  visits <- lapply(data[,4], FUN=gen_visits)
  
  # Generating (internal) time of UTT implementation
  tau <- runif(n, 0, max(do.call(c, visits)))
  
  data <- cbind(seq(1, n, by=1), data, tau)
  data.list.3 <- split(data, seq(nrow(data)))
  
  # Generating long-form dataset
  out <- do.call(rbind, mapply(FUN=function(y, v){
    stratum <- y[2]; group <- y[3]; id <- y[1]
    u.t <- y[4]; x.i <- y[5]; tau.i <- y[6]
    baseline.i <- baseline[[stratum]]
    gen_individual_tvc(u.t, tau.i, beta, x.i, baseline.i, v, cens, stratum, group, id)
  }, 
  y=data.list.3, v=visits, SIMPLIFY=FALSE))
  
  # Returning the final dataset
  return(out)
}

### Generates stratified dataset with a cluster- and time-constant binary variable and time-varying UTT indicator
  # M               numeric             number of communities 
  # ni              numeric             number (range) of individuals in each community 
  # beta            numeric             coefficient vector for (1) the time-constant variable, (2) the time-varying variable, and (3) their interaction
  # gen_visits      function            function for generating internal visit times
  # cens            numeric             loss-to-follow-up rate
  # rho             numeric             Kendall's tau for the dependence within community sub-networks
  # baseline        list                list containing the baseline distribution functions for each stratum
  # S               function            either (a) a function for generating cluster-varying stratum membership or (b) a vector indicating stratum membership for each cluster
gen_time_varying_cov <- function(M, ni, beta, gen_visits, cens, rho=0, baseline, S){
  
  # Generating the design matrix if observations are independent within each cluster
  if (rho==0){
    
    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      n.i <- sample(ni[1]:ni[2], 1)
      u <- runif(n.i, 0, 1)
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }
    
    # Generating the design matrix if observations belong to an exchangeable subnetwork within each cluster
  } else {
    
    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      networks <- rpois(ni[2]/2, 2) + 2
      networks.cm <- cumsum(networks)
      n.i <- sample(ni[1]:ni[2], 1)
      sub.ni <- networks[networks.cm < n.i]
      copula.fn.list <- lapply(sub.ni, FUN=function(x) claytonCopula(-2*rho/(rho-1), dim = x))
      u <- do.call(c, lapply(copula.fn.list, FUN=function(x) as.numeric(rCopula(1, x))))
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }
    
  }
  
  # Determining stratum membership
  if (is.function(S)) {
    data.list.2 <- lapply(data.list, FUN=function(dat){
      n.i <- nrow(dat)
      s <- S(n.i)
      cbind(s, dat)
    })
  } else{
    data.list.2 <- mapply(FUN=function(s, dat){
      cbind(rep(s, nrow(dat)), dat)
    }, s=S, dat=data.list, SIMPLIFY=FALSE)
  }
  
  data <- do.call(rbind, data.list.2)
  n <- nrow(data)
  
  # Generating (internal) inspection times
  visits <- gen_visits(n)
  
  # Generating (internal) time of UTT implementation
  tau <- runif(n, 0, max(visits[, ncol(visits)]))
  
  data <- cbind(seq(1, n, by=1), data, tau, visits)
  
  # Generating long-form dataset
  out <- do.call(rbind, apply(data, 1, FUN=function(y){
    stratum <- y[2]; group <- y[3]; id <- y[1]
    u.t <- y[4]; x.i <- y[5]; tau.i <- y[6]; visits.i <- y[7:length(y)]
    baseline.i <- baseline[[stratum]]
    gen_individual_tvc(u.t, tau.i, beta, x.i, baseline.i, visits.i, cens, stratum, group, id)
  }))
  
  # Returning the final dataset
  return(out)
}

#-----------------------------------------------------------------------------------------------------
# Data Generation Functions: Time-Dependent Effect of Combination Prevention
#-----------------------------------------------------------------------------------------------------

### List of possible forms for the time-dependent covariate effect 
  # u               numeric             time at which to evaluate the covariate effect
  # beta            numeric             coefficient vector for the (1) time constant and (2) time dependent components of the effect of x.i
  # tau             numeric             internal time of the discontinuity assuming a piecewise covariate effect 
tvbeta_fns <- list(
  log = list(
    name = "log",
    args = c("u", "beta"),
    beta.t = function(u, beta) beta[1] + beta[2]*log(u)
  ),
  linear = list(
    name = "linear",
    args = c("u", "beta"),
    beta.t = function(u, beta) beta[1] + beta[2]*u
  ),
  piecewise = list(
    name = "piecewise",
    args = c("u", "beta", "tau"),
    beta.t = function(u, beta, tau) beta[1] + ifelse(u > tau, 1, 0)*beta[2]
  )
)

### Returns individual survival curve
  # baseline.i      list                list of baseline survival, density, and hazard for stratum i
  # bfns            list                list containing the function for the time-dependent covariate effect, beta(t)
  # aux             list                arguments for the time-dependent covariate effect
  # x.i             numeric             indicator of intervention assignment
get_surv_tvb <- function(baseline.i, bfns, aux, x.i){
  
  h.i <- baseline.i$hazard
  
  # Obtaining the linear predictor at each timepoint
  aux$u <- h.i$time
  lp <- do.call(bfns$beta.t, aux)*x.i
  if (sum(is.na(lp)) > 0){
    na.index <- which(is.na(lp))
    if (length(na.index)==1){
      if (na.index != length(lp)){
        lp[na.index] <- lp[na.index + 1]
      } else if (na.index == length(lp)){
        lp[na.index] <- lp[na.index -1]
      }
    } else{
      for (i in na.index){
        if (i==length(lp)){
          lp[i] <- lp[tail(seq(1, i-1, by=1)[which(!is.na(lp[1:(i-1)]))], 1)]
        } else if (length(seq(i+1, length(lp), by=1)[which(!is.na(lp[(i+1):length(lp)]))]) > 0){
          lp[i] <- lp[head(seq(i+1, length(lp), by=1)[which(!is.na(lp[(i+1):length(lp)]))], 1)]
        } else{
          lp[i] <- lp[tail(seq(1, i-1, by=1)[which(!is.na(lp[1:(i-1)]))], 1)]
        }
      }
    }
  }
  if (sum(lp==Inf) > 0){
    inf.index <- which(lp==Inf)
    if ((length(inf.index)==1) & (inf.index != length(lp))){
      lp[inf.index] <- lp[inf.index + 1]
    } else if ((length(inf.index)==1) & (inf.index == length(lp))){
      lp[inf.index] <- lp[inf.index -1]
    } else {
      for (i in inf.index){
        lp[i] <- ifelse(i==length(lp), lp[tail(seq(1, i-1, by=1)[which(!(lp[1:(i-1)]==Inf))], 1)], 
                        lp[head(seq(i+1, length(lp), by=1)[which(!(lp[(i+1):length(lp)]==Inf))], 1)])
      }
    }
  }
  h.ij <-  h.i$hazard*exp(lp)
  
  # Approximating the cumulative hazard
  H.ij <- cumsum(h.ij[-length(h.ij)]*diff(h.i$time))
  
  # Approximating the corresponding survival
  S.ij <- exp(-H.ij)
  # Replacing missing/NA/infinite survival times with zero
  if (sum(is.na(S.ij)) > 0){
    na.index <- which(is.na(S.ij))
    S.ij[na.index] <- 0
  }
  if (sum(S.ij==Inf) > 0){
    inf.index <- which(S.ij==Inf)
    S.ij[inf.index] <- 0
  }
  
  s.curve <- data.frame(time=h.i$time[-length(h.i$time)], surv=S.ij)
  return(s.curve)
}

### Returns failure time associated with an individual survival probability
  # s               numeric             individual survival probability
  # baseline.i      list                list of baseline survival, density, and hazard for stratum i
  # bfns            list                list containing the function for the time-dependent covariate effect, beta(t)
  # aux             list                arguments for the time-dependent covariate effect
  # x.i             numeric             indicator of intervention assignment
get_inv_surv_tvb <- function(s, baseline.i, bfns, aux, x.i){
  
  # Obtaining piecewise-approximated individual survival function
  s.curve.i <- get_surv_tvb(baseline.i, bfns, aux, x.i)
  
  # Returning the failure time corresponding to s
  index <- which(s.curve.i$surv <= s)[1]
  t <- s.curve.i$time[index]
  return(t)
}


### Generates individual data vector with time-dependent combination prevention effect
  # u.t             numeric             uniform random variate
  # baseline.i      list                list of baseline survival, density, and hazard for stratum i
  # bfns            list                list containing the function for the time-dependent covariate effect, beta(t)
  # aux             list                arguments of the underlying time-to-event distribution
  # x.i             numeric             indicator of intervention assignment
  # visits.i        numeric             internal time of all study-mandated HIV testing dates
  # cens            numeric             loss-to-follow-up rate
  # stratum         numeric             ID number of the associated stratum
  # group           numeric             ID number of the associated community
gen_individual_tvbeta <- function(u.t, baseline.i, bfns, aux, x.i, visits.i, cens, stratum, group){
  
  # Generating survival and censoring times
  s.time <- get_inv_surv_tvb(u.t, baseline.i, bfns, aux, x.i)
  if (cens==0){
    c.time <- visits.i[length(visits.i)] + 1
  } else{
    u.c <- runif(1)
    c.time <- -log(u.c)/cens
  }
  
  # Deriving final interval censored observation and return (left, right)
  visits.i[visits.i > c.time] <- NA
  event <- as.numeric(visits.i > s.time)
  if (sum(is.na(event))==length(event)){
    out <- data.frame(stratum=stratum, group=group, treat=x.i, left=0, right=Inf)
    return(out)
  }
  event.status <- max(event, na.rm=TRUE)
  if (event.status==1){
    first.event.time <- min(which(event==1))
    l <- ifelse(first.event.time==1, 0, 
                ifelse(sum(!is.na(visits.i[1:(first.event.time-1)]))==0, 0,
                       max(visits.i[1:(first.event.time-1)], na.rm=TRUE)))
    r <- visits.i[first.event.time]
  } else{
    l <- max(visits.i, na.rm=TRUE)
    r <- Inf
  }
  
  # Returning the final observation for individual i
  out <- data.frame(stratum=stratum, group=group, treat=x.i, left=l, right=r)
  return(out)
}

### Generates stratified data with a time-dependent covariate effect for combination prevention
  # M               numeric             number of communities 
  # ni              numeric             number (range) of individuals in each community 
  # beta            list                list containing the function for the time-dependent covariate effect, beta(t)
  # gen_visits      function            function for generating internal visit times
  # cens            numeric             loss-to-follow-up rate
  # rho             numeric             Kendall's tau for the dependence within community sub-networks
  # args            list                optional beta, lambda, and tau arguments for the underlying time-to-event distribution
  # baseline        list                list containing the baseline distribution functions for each stratum
  # S               function            either (a) a function for generating cluster-varying stratum membership or (b) a vector indicating stratum membership for each cluster
gen_time_dependent_beta_tx <- function(M, ni, beta, beta.args, gen_visits, cens, rho=0, baseline, S){
  
  # Determining beta(t) parametrization
  if (missing(beta)){
    stop("Time-dependent covariate effect \"beta\" not specified")
  } else if (is.list(beta)){
    bfns <- beta
  } else {
    beta <- match.arg(tolower(beta), tolower(names(tvbeta_fns)))
    bfns <- tvbeta_fns[[beta]]
  }
  
  # Generating the design matrix if observations are independent within each stratum
  if (rho==0){
    
    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      n.i <- sample(ni[1]:ni[2], 1)
      u <- runif(n.i, 0, 1)
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }
    
    # Generating the design matrix if observations belong to an exchangeable subnetwork within each stratum
  } else {
    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      networks <- rpois(ni[2]/2, 2) + 2
      networks.cm <- cumsum(networks)
      n.i <- sample(ni[1]:ni[2], 1)
      sub.ni <- networks[networks.cm < n.i]
      copula.fn.list <- lapply(sub.ni, FUN=function(x) claytonCopula(-2*rho/(rho-1), dim = x))
      u <- do.call(c, lapply(copula.fn.list, FUN=function(x) as.numeric(rCopula(1, x))))
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }
    
  }
  
  # Determining stratum membership
  if (is.function(S)) {
    data.list.2 <- lapply(data.list, FUN=function(dat){
      n.i <- nrow(dat)
      s <- S(n.i)
      cbind(s, dat)
    })
  } else{
    data.list.2 <- mapply(FUN=function(s, dat){
      cbind(rep(s, nrow(dat)), dat)
    }, s=S, dat=data.list, SIMPLIFY=FALSE)
  }
  
  data <- do.call(rbind, data.list.2)
  n <- nrow(data)
  
  # Generating (internal) inspection times
  visits <- lapply(data[,4], FUN=gen_visits)
  data.list.3 <- split(data, seq(nrow(data)))
  
  # Generating dataset
  out <- do.call(rbind, mapply(FUN=function(y, v){
    stratum <- y[1]; group <- y[2]; u.t <- y[3]; x.i <- y[4]
    baseline.i <- baseline[[stratum]]
    gen_individual_tvbeta(u.t, baseline.i, bfns, beta.args, x.i, v, cens, stratum, group)
  }, 
  y=data.list.3, v=visits, SIMPLIFY=FALSE))
  
  # Returning the final dataset
  return(out)
}

### Generates stratified data with a time-dependent covariate effect for combination prevention
  # M               numeric             number of communities 
  # ni              numeric             number (range) of individuals in each community 
  # beta            list                list containing the function for the time-dependent covariate effect, beta(t)
  # gen_visits      function            function for generating internal visit times
  # cens            numeric             loss-to-follow-up rate
  # rho             numeric             Kendall's tau for the dependence within community sub-networks
  # args            list                optional beta, lambda, and tau arguments for the underlying time-to-event distribution
  # baseline        list                list containing the baseline distribution functions for each stratum
  # S               function            either (a) a function for generating cluster-varying stratum membership or (b) a vector indicating stratum membership for each cluster
gen_time_dependent_beta <- function(M, ni, beta, beta.args, gen_visits, cens, rho=0, baseline, S){
  
  # Determining beta(t) parametrization
  if (missing(beta)){
    stop("Time-dependent covariate effect \"beta\" not specified")
  } else if (is.list(beta)){
    bfns <- beta
  } else {
    beta <- match.arg(tolower(beta), tolower(names(tvbeta_fns)))
    bfns <- tvbeta_fns[[beta]]
  }
  
  # Generating the design matrix if observations are independent within each stratum
  if (rho==0){
    
    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      n.i <- sample(ni[1]:ni[2], 1)
      u <- runif(n.i, 0, 1)
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }
    
    # Generating the design matrix if observations belong to an exchangeable subnetwork within each stratum
  } else {
    
    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      networks <- rpois(ni[2]/2, 2) + 2
      networks.cm <- cumsum(networks)
      n.i <- sample(ni[1]:ni[2], 1)
      sub.ni <- networks[networks.cm < n.i]
      copula.fn.list <- lapply(sub.ni, FUN=function(x) claytonCopula(-2*rho/(rho-1), dim = x))
      u <- do.call(c, lapply(copula.fn.list, FUN=function(x) as.numeric(rCopula(1, x))))
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }
    
  }
  
  # Determining stratum membership
  if (is.function(S)) {
    data.list.2 <- lapply(data.list, FUN=function(dat){
      n.i <- nrow(dat)
      s <- S(n.i)
      cbind(s, dat)
    })
  } else{
    data.list.2 <- mapply(FUN=function(s, dat){
      cbind(rep(s, nrow(dat)), dat)
    }, s=S, dat=data.list, SIMPLIFY=FALSE)
  }
  
  data <- do.call(rbind, data.list.2)
  n <- nrow(data)
  
  # Generating (internal) inspection times
  visits <- gen_visits(n)
  data <- cbind(data, visits)
  
  # Generating dataset
  out <- do.call(rbind, apply(data, 1, FUN=function(y){
    stratum <- y[1]; group <- y[2]; u.t <- y[3]; x.i <- y[4]
    visits.i <- y[5:length(y)]
    baseline.i <- baseline[[stratum]]
    gen_individual_tvbeta(u.t, baseline.i, bfns, beta.args, x.i, visits.i, cens, stratum, group)
  }))
  
  # Returning the final dataset
  return(out)
}


#-----------------------------------------------------------------------------------------------------
# Data Generation Functions: Possible Monitoring Schemes
#-----------------------------------------------------------------------------------------------------

### Covariate-Dependent Monitoring Scheme
  # x               numeric             treatment assignment
gen_visits_tx <- function(x){
  if (x==1){
    temp <- cumsum(runif(15, 12, 24))
    temp <- c(temp, sapply(1:4, FUN=function(y)runif(1, 52*y-4, 52*y+4)))
    temp <- sort(temp)
    return(temp)
  } else if (x==0){
    return(sapply(1:4, FUN=function(y)runif(1, 52*y-4, 52*y+4)))
  }
}

### Yearly Monitoring Scheme
  # n               numeric             total sample size
gen_visits_bcpp <- function(n){
  visits <- do.call(rbind, lapply(seq_len(n), FUN=function(x)sapply(1:4, FUN=function(y)runif(1, 52*y-4, 52*y+4))))
  return(visits)
}

### Frequent Monitoring Scheme
  # n               numeric             total sample size
gen_visits_freq <- function(n){
  visits <- do.call(rbind, lapply(seq_len(n), FUN=function(x)cumsum(runif(20, 0, 16))))
  return(visits)
}
