#-----------------------------------------------------------------------------------------------------#
#         Poisson-Augmented Estimation of Stratified Cox Models for Interval-Censored Data:           #
#                 Functions for Point and Interval Estimation of Model Parameters                     #
#-----------------------------------------------------------------------------------------------------#

library(rootSolve)
library(abind)
library(ggplot2)


#-----------------------------------------------------------------------------------------------------
# Data Processing Functions
#-----------------------------------------------------------------------------------------------------

### Processes the input dataset
  # formula         formula             formula object for the stratified Cox model
  # data            data.frame          full set of observations 
  # strata          character           vector indicating the stratification variable(s)
  # id              character           individual subject identification variable (if there are time-varying covariates)
  # clus            character           vector indicating the clustering variable(s)
process_data <- function(formula, data, strata, id, clus){
  
  row.names(data) <- 1:nrow(data)
  
  # Creating list of datasets for each stratum
  if (!is.null(strata)){
    data.est <- split(data, data[, as.character(strata)])
  } else{
    data.est <- list(data)
  }
  
  # Creating processed data structures for each stratum
  processed.list <- lapply(data.est, FUN=function(d)process_data_s(formula, d, id))
  
  # Forming processed data object for point estimation
  processed <- list(l=lapply(processed.list, '[[', 'l.s'), 
                    u=lapply(processed.list, '[[', 'u.s'), 
                    u.star=lapply(processed.list, '[[', 'u.s.star'), 
                    tau=lapply(processed.list, '[[', 'tau.s'), 
                    risk=lapply(processed.list, '[[', 'risk.vr'),
                    x=lapply(processed.list, '[[', 'x.s'))
  
  if (is.null(clus) & is.null(strata)){
    
    # Unstratified model with model-based variance estimation
    mapping <- NULL
    
  } else if (is.null(clus)){
    
    # Creating list of datasets for each cluster
    data.clus <- list(data)
    
    # Creating the cluster-to-stratum mappings for each cluster
    mapping.list <- lapply(data.clus, FUN=function(d)cluster_to_stratum_mapping_i(d, data.est, strata, id))
    
    # Forming the final mapping object
    mapping <- list(stratum=lapply(mapping.list, '[[', 'stratum.id'),
                    indiv=lapply(mapping.list, '[[', 'indiv.id'))
    
  } else if (is.null(strata)){
    
    # Creating list of datasets for each cluster
    data.clus <- split(data, data[, as.character(clus)])
    
    # Creating the cluster-to-stratum mappings for each cluster
    mapping.list <- lapply(data.clus, FUN=function(d)cluster_to_stratum_mapping_i(d, data.est, strata, id))
    
    # Forming the final mapping object
    mapping <- list(stratum=lapply(mapping.list, '[[', 'stratum.id'),
                    indiv=lapply(mapping.list, '[[', 'indiv.id'))
    
  } else if (clus == strata){
    
    # Stratification and clustering variable(s) coincide
    mapping <- NULL
    
  } else {
    
    # Creating list of datasets for each cluster
    data.clus <- split(data, data[, as.character(clus)])
    
    # Creating the cluster-to-stratum mappings for each cluster
    mapping.list <- lapply(data.clus, FUN=function(d)cluster_to_stratum_mapping_i(d, data.est, strata, id))
    
    # Forming the final mapping object
    mapping <- list(stratum=lapply(mapping.list, '[[', 'stratum.id'),
                    indiv=lapply(mapping.list, '[[', 'indiv.id'))
    
  }
  
  return(list(pointEstimation=processed, varMapping=mapping))
}

### Processes the data for stratum s
  # formula         formula             formula object for the stratified Cox model
  # data            data.frame          observations in stratum s
  # id              character           individual subject identification variable (if there are time-varying covariates)
process_data_s <- function(formula, data, id){
  
  # Determining relevant column names and indices
  vars <- all.vars(formula); l <- vars[1]; u <- vars[2]
  
  # Identifying one unique row for each individual
  if (!is.null(id)){
    id.s <- data[, as.character(id)]
    row.id.s <- match(unique(id.s), id.s) 
  } else{
    row.id.s <- 1:nrow(data)
  }
  
  # Creating outcome, time, and riskset variables
  l.s <- data[row.id.s, l]; u.s <- data[row.id.s, u]
  u.s.star <- ifelse(u.s==Inf, l.s, u.s)
  tau.s <- sort(unique(c(l.s[l.s > 0], u.s[u.s < Inf])))
  risk.vr <- do.call(rbind, lapply(u.s.star, FUN=function(x)as.numeric(x >= tau.s)))
  
  # Creating covariate array
  mat <- model.matrix(formula[c(1, 3)], data=data)
  mat <- mat[, -which(colnames(mat)=="(Intercept)"), drop=FALSE]
  if (is.null(id)){
    x.s <- array(rep(mat, length(tau.s)), dim=c(nrow(mat), ncol(mat), length(tau.s)))
  } else{
    temp <- matrix(NA, nrow=length(unique(id.s))*ncol(mat), ncol=length(tau.s))
    temp.index <- seq(1, length(unique(id.s))*(ncol(mat)), by=length(unique(id.s)))
    for (v in unique(id.s)){
      index <- which(id.s==v)
      temp.time <- data[index, c("start", "stop")]
      temp.mat <- mat[index, ]
      temp[temp.index, ] <- do.call(cbind, lapply(tau.s, FUN=function(r){
        if (length(index)==1){
          return(temp.mat)
        } else{
          cov.index <- sum(temp.time[, "start"] < r)
          if (!is.matrix(temp.mat)){
            return(temp.mat[cov.index])
          } else{
            return(temp.mat[cov.index, ])
          }
        }
      }))
      temp.index <- temp.index + 1
    }
    x.s <- array(temp, dim=c(length(unique(id.s)), ncol(mat), length(tau.s)))
  }
  
  dimnames(x.s)[[2]] <- as.list(colnames(mat))
  
  # Returning the processed data structures for stratum s
  return(list(l.s=l.s, u.s=u.s, u.s.star=u.s.star, tau.s=tau.s, risk.vr=risk.vr, x.s=x.s))
  
}

### Determines the mapping between cluster membership and stratum membership
  # data            data.frame          observations in cluster i
  # stratum.list    list                list of data frames corresponding to each stratum
  # strata          character           vector indicating the stratification variable(s)
  # id              character           individual subject identification variable (if there are time-varying covariates)
cluster_to_stratum_mapping_i <- function(data, stratum.list, strata, id){
  
  # Identifying one unique row for each individual
  if (!is.null(id)){
    id.i <- data[, as.character(id)]
    row.id.i <- match(unique(id.i), id.i) 
  } else{
    row.id.i <- 1:nrow(data)
  }
  
  # Mapping each individual in cluster i to their corresponding stratum
  if (is.null(strata)){
    stratum.id <- rep(1, nrow(data[row.id.i, ]))
  } else {
    stratum.id <- apply(data[row.id.i, as.character(strata), drop=FALSE], 1, FUN=function(r){
      which(names(stratum.list) == paste0(as.character(r), collapse="."))
    })
  }
  
  # Mapping each individual in cluster i to their individual index in the corresponding stratum
  indiv.id <- mapply(function(row.name.ij, stratum.ij){
    strat <- stratum.list[[stratum.ij]]
    if (!is.null(id)){
      strat.id <- strat[, as.character(id)]
      strat.row.id <- match(unique(strat.id), strat.id)
    } else {
      strat.row.id <- 1:nrow(strat)
    }
    which(row.names(strat[strat.row.id, ]) == row.name.ij)
  }, row.name.ij = row.names(data[row.id.i, ]), stratum.ij = stratum.id)
  
  
  return(list(stratum.id=stratum.id, indiv.id=indiv.id))
  
}


#-----------------------------------------------------------------------------------------------------
# Model Estimation Wrapper Function
#-----------------------------------------------------------------------------------------------------

### Obtains composite maximum likelihood estimators for marginal Cox regression models with clustered interval-censored data
  # formula         formula             formula object for the marginal Cox model
  # data            data.frame/list     full set of observations; may be either an unprocessed data frame (in which case an accompanying formula object is also required) or a pre-processed list containing:
  #                                       $pointEstimation:
  #                                         l: list of length S containing l.s for each stratum
  #                                         u: list of length S containing u.s for each stratum
  #                                         u.star: list of length S containing u.star.s for each stratum
  #                                         tau: list of length S containing tau.s for each stratum
  #                                         x: list of length S containing x.s for each stratum
  #                                         risk: list of length S containing risk.vr for each stratum
  #                                       $varMapping: may be NULL, otherwise:
  #                                         stratum: list of length M containing the stratum id s for each observation in cluster i
  #                                         indiv: list of length M containing the stratum-specific subject id v for each observation in cluster i
  # strata          character           vector indicating the stratification variable(s)
  # clus            character           vector indicating the clustering variable(s)
  # id              character           individual subject identification variable (if there are time-varying covariates)
  # variance        logical             should the variance for the regression parameter estimators be computed and returned?
  # method          character           indicates whether the variance should be calculated using the "profile" composite likelihood or "bootstrap"
  # control         list                control parameters for the composite EM algorithm and variance estimation, namely:
  #                                       init: list of initial values for the coefficients and the baseline hazard(s);
  #                                       tol: stopping criterion (in terms of the l1 norm of the coefficient estimates) for the composite EM algorithm;
  #                                       maxit: maximum number of composite EM iterations;
  #                                       c: perturbation constant for the finite differences approximations;
  #                                       trace: should the full sequence of composite EM iterates be returned?
  # boot.control    list                control parameters for the bootstrap variance estimation procedure, namely:
  #                                       perturb: Boolean indicating whether a perturbation bootstrap (TRUE) or nonparametric bootstrap (FALSE) should be performed
  #                                       resample.clusters: Boolean indicating whether, if clus is non-null, clusters (TRUE) or individuals within clusters (FALSE) should be resampled
  #                                       n.boot: number of bootstrap resamples to use for variance calculation
  # boot.data       list                optional list of pre-processed bootstrap datasets (for use if data object is also processed and boot.control$perturb==FALSE)
composite_coxIC <- function(formula=NULL, data, strata=NULL, clus=NULL, id=NULL, variance=TRUE, method="profile", control=NULL, boot.control=NULL, boot.data=NULL){
  
  call <- match.call()
  
  args <- list(init=NULL, tol=1e-8, maxit.fit=1000, maxit.var=1000, c=1, trace=FALSE)
  if (!is.null(control)){
    if ("maxit" %in% names(control)){
      args$maxit.fit <- control$maxit
      args$maxit.var <- control$maxit
      args[setdiff(names(control), "maxit")] <- control[-which(names(control)=="maxit")]
    } else {
      args[names(control)] <- control
    }
  }
  
  boot.args <- list(perturb=FALSE, resample.clusters=FALSE, n.boot=100)
  if (!is.null(boot.control)){
    boot.args[names(boot.control)] <- boot.control
  }
  if (!is.data.frame(data) & method=="bootstrap" & is.null(boot.data) & boot.args$perturb==FALSE){
    stop("To conduct bootstrap inference for analyses with pre-processed data, provide a list of pre-processed bootstrapped datasets using the boot.data option.")
  }
  
  # Processing the input data
  if (is.data.frame(data)){
    temp <- process_data(formula, data, strata, id, clus)
    processed <- temp$pointEstimation
    var.mapping <- temp$varMapping
  } else{
    processed <- data$pointEstimation
    var.mapping <- data$varMapping
  }
  l <- processed$l; u <- processed$u; u.star <- processed$u.star; tau <- processed$tau
  risk <- processed$risk
  x <- processed$x
  
  # Initializing the point estimates
  if (is.null(args$init)){
    init <- em_init(dim(x[[1]])[2], tau)
  } else {
    init <- args$init
  }
  beta.current <- init$beta; lambda.current <- init$lambda
  
  # Initializing lists to store point estimates at each iteration
  if (args$trace){
    trace.beta <- list(beta.current)
    trace.lambda <- list(lambda.current)
  }
  
  # Iterating until convergence 
  fit.start <- Sys.time()
  l1.norm <- 1; n.iter <- 0
  while (l1.norm > args$tol & n.iter < args$maxit.fit){
    update <- point_iter(beta.current, lambda.current, l, u, u.star, tau, x, risk)
    l1.norm <- sum(abs(update$beta - beta.current))
    beta.current <- update$beta; lambda.current <- update$lambda
    if (args$trace){
      trace.beta[[length(trace.beta)+1]] <- beta.current
      trace.lambda[[length(trace.lambda)+1]] <- lambda.current
    }
    n.iter <- n.iter + 1
  }
  fit.end <- Sys.time()
  
  # Estimating the covariance matrix for beta
  var.start <- Sys.time()
  if (variance){
    if (method=="profile"){
      beta.cov <- variance_beta(beta.current, processed, var.mapping, args)
      colnames(beta.cov) <- dimnames(x[[1]])[[2]]
    } else if (method=="bootstrap"){
      # Creating the bootstrap distribution from user-provided processed datasets
      if (!is.null(boot.data)){
        boot.processed <- lapply(boot.data, FUN=function(b) b$pointEstimation)
        boot.weights <- lapply(boot.processed, FUN=function(b){
          lapply(b$l, FUN=function(b.l.i) rep(1, length(b.l.i)))
        })
      # Creating the bootstrap distribution by perturbing individuals 
      } else if (boot.args$perturb==TRUE & (is.null(var.mapping) & length(processed$l)==1)){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          boot.processed[[length(boot.processed) + 1]] <- processed
          boot.weights[[length(boot.weights) + 1]] <- list(rexp(sum(do.call(c, lapply(processed$x, FUN=function(x.s)dim(x.s)[1]))), 1))
        }
      } else if (boot.args$perturb==TRUE & length(var.mapping$stratum)==1){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          boot.processed[[length(boot.processed) + 1]] <- processed
          b.weights.vec <- rexp(sum(do.call(c, lapply(processed$x, FUN=function(x.s)dim(x.s)[1]))), 1)
          b.weights <- list(); n.s.vec <- c(0, do.call(c, lapply(processed$l, FUN=function(b.l.s)length(b.l.s))))
          for (s in 1:length(processed$l)){
            b.weights[[length(b.weights) + 1]] <- b.weights.vec[(n.s.vec[s] + 1):n.s.vec[s+1]]
          }
          boot.weights[[length(boot.weights) + 1]] <- b.weights
        }
      # Creating the bootstrap distribution by perturbing clusters
      } else if (boot.args$perturb==TRUE & (is.null(var.mapping) & length(processed$l) > 1)){
        warning("Use of the perturbation bootstrap is discouraged when the clustering and stratification factors coincide.")
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          boot.processed[[length(boot.processed) + 1]] <- processed
          boot.weights[[length(boot.weights) + 1]] <- lapply(processed$l, FUN=function(b.l.s){
            b.weight.s <- rexp(1, 1)
            rep(b.weight.s, length(b.l.s))
          })
        }
      } else if (boot.args$perturb==TRUE & (length(var.mapping$stratum) > 1)){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          boot.processed[[length(boot.processed) + 1]] <- processed
          b.weights <- lapply(processed$l, FUN=function(l.s) rep(NA, length(l.s)))
          b.weights.i <- rexp(length(var.mapping$stratum), 1)
          for (i in 1:length(var.mapping$stratum)){
            for (j in 1:length(var.mapping$stratum[[i]])){
              b.weights[[var.mapping$stratum[[i]][j]]][var.mapping$indiv[[i]][j]] <- b.weights.i[i]
            }
          }
          boot.weights[[length(boot.weights) + 1]] <- b.weights
        }
      # Creating the bootstrap distribution by resampling individuals
      } else if (boot.args$perturb==FALSE & is.null(clus)){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          b.index <- sample(nrow(data), nrow(data), replace=TRUE)
          b.data <- data[b.index, ]
          b.temp <- process_data(formula, b.data, strata, id, clus)
          boot.processed[[length(boot.processed) + 1]] <- b.temp$pointEstimation
          boot.weights[[length(boot.weights) + 1]] <- lapply(b.temp$pointEstimation$l, FUN=function(b.l.s) rep(1, length(b.l.s)))
        }
      # Creating the bootstrap distribution by resampling individuals within clusters: perturb is FALSE, there are clusters, resample.clusters=FALSE
      } else if (boot.args$perturb==FALSE & !is.null(clus) & boot.args$resample.clusters==FALSE){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          data.list <- split(data, data[, as.character(clus)])
          b.data.list <- lapply(data.list, FUN=function(b.d){
            b.index <- sample(nrow(b.d), nrow(b.d), replace=TRUE)
            b.d[b.index, ]
          })
          b.data <- do.call(rbind, b.data.list)
          b.temp <- process_data(formula, b.data, strata, id, clus)
          boot.processed[[length(boot.processed) + 1]] <- b.temp$pointEstimation
          boot.weights[[length(boot.weights) + 1]] <- lapply(b.temp$pointEstimation$l, FUN=function(b.l.s) rep(1, length(b.l.s)))
        }
      # Creating the bootstrap distribution by resampling clusters: perturb is FALSE, there are clusters, resample.clusters=TRUE
      } else if (boot.args$perturb==FALSE & !is.null(clus) & boot.args$resample.clusters==TRUE){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          data.list <- split(data, data[, as.character(clus)])
          b.index <- sample(length(data.list), length(data.list), replace=TRUE)
          b.data.list <- data.list[b.index]
          b.data <- do.call(rbind, b.data.list)
          b.temp <- process_data(formula, b.data, strata, id, clus)
          boot.processed[[length(boot.processed) + 1]] <- b.temp$pointEstimation
          boot.weights[[length(boot.weights) + 1]] <- lapply(b.temp$pointEstimation$l, FUN=function(b.l.s) rep(1, length(b.l.s)))
        }
      }
      beta.boot <- list()
      for (b in 1:length(boot.processed)){
        b.l <- boot.processed[[b]]$l; b.u <- boot.processed[[b]]$u; b.u.star <- boot.processed[[b]]$u.star; b.tau <- boot.processed[[b]]$tau
        b.risk <- boot.processed[[b]]$risk
        b.x <- boot.processed[[b]]$x
        b.init <- em_init(dim(b.x[[1]])[2], b.tau)
        b.beta.current <- b.init$beta; b.lambda.current <- b.init$lambda
        b.l1.norm <- 1; b.n.iter <- 0
        while (b.l1.norm > args$tol & b.n.iter < args$maxit.var){
          b.update <- point_iter(b.beta.current, b.lambda.current, b.l, b.u, b.u.star, b.tau, b.x, b.risk, weights=boot.weights[[b]])
          b.l1.norm <- sum(abs(b.update$beta - b.beta.current))
          b.beta.current <- b.update$beta; b.lambda.current <- b.update$lambda
          b.n.iter <- b.n.iter + 1
        }
        beta.boot[[length(beta.boot)+1]] <- b.beta.current
      }
      beta.cov <- var(do.call(rbind, beta.boot))
      colnames(beta.cov) <- dimnames(x[[1]])[[2]]
    }
  } else {
    beta.cov <- NULL
  }
  var.end <- Sys.time()
  
  # Formatting the final results
  lambda.mat <- lapply(seq_len(length(lambda.current)), FUN=function(x){
    temp <- cbind(tau[[x]], lambda.current[[x]])
    colnames(temp) <- c("tau.sr", "lambda.sr")
    return(temp)
  })
  names(lambda.mat) <- names(lambda.current)
  names(beta.current) <- dimnames(x[[1]])[[2]]
  trace.out <- NULL
  if (args$trace){
    trace.lambda.2 <- lapply(seq_len(length(trace.lambda[[1]])), FUN=function(x)lapply(trace.lambda, `[[`, x))
    names(trace.lambda.2) <- names(trace.lambda[[1]])
    trace.out <- list('beta'=trace.beta, 'lambda'=trace.lambda.2)
  }
  if (is.null(var.mapping) & length(processed$l)==1){
    var.type <- "model-based"
  } else if (is.null(var.mapping) & length(processed$l) > 1){
    var.type <- "robust"
  } else if (length(var.mapping$stratum)==1){
    var.type <- "model-based"
  } else {
    var.type <- "robust"
  }
  fit.time <- difftime(fit.end, fit.start, units="min")
  var.time <- ifelse(variance==TRUE, difftime(var.end, var.start, units="min"), NA)
  
  out <- list(call=call, coefficients=beta.current, var=beta.cov, var.type=var.type, var.method=method,
              baseline.hazard=lambda.mat, n.iter=n.iter, l1.norm=l1.norm, trace=trace.out,
              fit.times=list('model'=fit.time, 'variance'=var.time))
  class(out) <- 'compCoxIC'
  return(out)
}


#-----------------------------------------------------------------------------------------------------
# Initialization of the Composite EM Algorithm
#-----------------------------------------------------------------------------------------------------

### Initializes the coefficient and baseline hazard terms
  # p               integer             dimension of the coefficient vector
  # tau             list                list of length S containing tau.s for each stratum
em_init <- function(p, tau){
  
  # Setting the coefficient vector to zero and the baseline hazard function to 1/rho.s for each stratum
  beta.init <- rep(0, p)
  lambda.init <- lapply(tau, FUN=function(x)rep(1/length(x), length(x)))
  
  return(list(beta=beta.init, lambda=lambda.init))
}


#-----------------------------------------------------------------------------------------------------
# Single Iteration of the Composite EM Algorithm
#-----------------------------------------------------------------------------------------------------

### Performs a single iteration of the composite EM algorithm
  # beta            numeric             p-dimensional vector of the current fixed effect estimates 
  # lambda          list                list of length S containing the current estimate of lambda.s for each stratum
  # l               list                list of length S containing l.s for each stratum
  # u               list                list of length S containing u.s for each stratum
  # u.star          list                list of length S containing u.star.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # risk            list                list of length S containing risk.vr for each stratum
  # weights         list                list of length S containing the perturbation weights weight.s for each stratum
point_iter <- function(beta, lambda, l, u, u.star, tau, x, risk, weights=NULL){
  
  # Expectation step
  w <- calc_w(l, u, tau, x, lambda, beta)
  
  # Maximization step
  beta.iter <- update_beta(beta, u.star, tau, x, w, risk, weights)
  lambda.iter <- update_lambda(beta.iter, x, w, risk, weights)
  
  return(list(beta=beta.iter, lambda=lambda.iter))
}


#-----------------------------------------------------------------------------------------------------
# Expectation Step of the Composite EM Algorithm
#-----------------------------------------------------------------------------------------------------

### Calculates the conditional expectation of the Poisson augmentation variables
  # l               list                list of length S containing l.s for each stratum
  # u               list                list of length S containing u.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # lambda          list                list of length S containing the current estimate of lambda.s for each stratum
  # beta            numeric             p-dimensional vector of the current fixed effect estimates 
calc_w <- function(l, u, tau, x, lambda, beta){
  
  # Taking the projection of w.s for each stratum
  mapply(function(l.s, u.s, tau.s, x.s, lambda.s){
    w_s(l.s, u.s, tau.s, x.s, lambda.s, beta)
  }, l.s=l, u.s=u, tau.s=tau, x.s=x, lambda.s=lambda, SIMPLIFY=FALSE)
  
}

### Returns the matrix of projected w.s for stratum s
  # l.s             numeric             n.s-dimensional vector of left endpoints for stratum s
  # u.s             numeric             n.s-dimensional vector of right endpoints for stratum s
  # tau.s           numeric             rho.s-dimensional vector of sorted unique left and right endpoints for stratum s
  # x.s             array               covariate array so that x.s[v, , ] returns a p x rho.s matrix, x.s[ , p, ] an n.s x rho.s matrix, and x.s[ , , r] an n.s x p matrix
  # lambda.s        numeric             rho.s-dimensional vector of the current estimated baseline hazard function for stratum s
  # beta            numeric             p-dimensional vector of the current fixed effect estimates 
w_s <- function(l.s, u.s, tau.s, x.s, lambda.s, beta){
  
  # Taking the projection of w.sv for each individual
  do.call(rbind, mapply(FUN=function(l.sv, u.sv, x.pr) w_sv(l.sv, u.sv, tau.s, x.pr, lambda.s, beta), 
                        l.sv=l.s, u.sv=u.s, x.pr=lapply(seq_len(dim(x.s)[1]), FUN=function(x)matrix(x.s[x, , ], nrow=dim(x.s)[2])), SIMPLIFY=FALSE))
  
}

### Returns the vector w.sv for individual v in stratum s
  # l.sv            numeric             the left endpoint of the censoring interval for individual v in stratum s
  # u.sv            numeric             the right-endpoint of the censoring interval for individual v in stratum s
  # tau.s           numeric             rho.s-dimensional vector of sorted unique left and right endpoints for stratum s
  # x.pr            matrix              matrix (p x rho.s) containing the values of each covariate at each timepoint for individual v in stratum s
  # lambda.s        numeric             rho.s-dimensional vector of the current estimated baseline hazard function for stratum s
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
w_sv <- function(l.sv, u.sv, tau.s, x.pr, lambda.s, beta){
  
  # Initializing the vector of projections for all W_svr such that tau_sr <= L_sv
  w.sv <- rep(0, sum(tau.s <= l.sv))
  
  # Computing the conditional expectation for all W_svr such that tau_sr \in (L_sv, U_sv] for U_sv < infinity
  if (u.sv != Inf){
    r.ids <- (tau.s > l.sv & tau.s <= u.sv)
    lin.pred <- beta %*% x.pr[, r.ids]
    denom <- 1 - exp(-sum(lambda.s[r.ids] * exp(lin.pred)))
    w.sv <- c(w.sv, lambda.s[r.ids]*exp(lin.pred)/denom)
  }
  
  # Padding with additional zeroes for all W_svr such that tau_sr > U_sv* (to aid with computation and matrix algebra)
  w.sv <- c(w.sv, rep(0, length(tau.s) - length(w.sv)))
  
  return(w.sv)
}


#-----------------------------------------------------------------------------------------------------
# Maximization Step of the Composite EM Algorithm
#-----------------------------------------------------------------------------------------------------

### Solves the profile composite likelihood score equation for beta
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
  # u.star          list                list of length S containing u.star.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # w               list                list of length S containing w.s for each stratum
  # risk            list                list of length S containing risk.vr for each stratum
  # weights         list                list of length S containing the perturbation weights weight.s for each stratum
update_beta <- function(beta, u.star, tau, x, w, risk, weights=NULL){
  
  if (is.null(weights)){
    weights <- lapply(x, FUN=function(x.s) rep(1, dim(x.s)[1]))
  }
  
  # Finding and returning zeroes of the profile composite score equation for beta
  multiroot(beta_score, start=beta, u.star=u.star, tau=tau, x=x, w=w, risk=risk, weights=weights, maxiter=1)$root
  
}  

### Evaluates the profile composite likelihood score equation for beta
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
  # u.star          list                list of length S containing u.star.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # w               list                list of length S containing w.s for each stratum
  # risk            list                list of length S containing risk.vr for each stratum 
  # weights         list                list of length S containing the perturbation weights weight.s for each stratum
beta_score <- function(beta, u.star, tau, x, w, risk, weights){
  
  ### Evaluates the profile composite likelihood score equation in stratum s
    # beta            numeric             p-dimensional vector of the current fixed effect estimates
    # u.s.star        numeric             n.s-dimensional vector of u.sv.star for stratum s
    # tau.s           numeric             rho.s-dimensional vector of sorted unique left and right endpoints for stratum s
    # x.s             array               covariate array so that x.s[v, , ] returns a p x rho.s matrix, x.s[ , p, ] an n.s x rho.s matrix, and x.s[ , , r] an n.s x p matrix
    # w.s             matrix              matrix (n.s x rho.s) from the expectation step for stratum s
    # risk.vr         matrix              matrix (n.s x rho.s) indicating risk set membership for each individual at each timepoint in stratum s
    # weights.s       numeric             n.s-dimensional vector of perturbation weights for stratum s
  score_s <- function(beta, u.s.star, tau.s, x.s, w.s, risk.vr, weights.s){
    
    # Determining the weight matrix for each individual at each timepoint
    weight.vr <- risk.vr * w.s * weights.s
    
    # Calculating stratum-level score contribution for each covariate
    x.bar <- x_bar_s(beta, u.s.star, tau.s, x.s, risk.vr, weights.s)
    list.x.p <- lapply(seq_len(dim(x.s)[2]), function(x) x.s[, x, ])
    list.x.bar <- lapply(seq_len(ncol(x.bar)), function(x)x.bar[, x])
    score_out <- mapply(function(x.vr, x.bar){
      sum(rowSums(weight.vr * t(apply(x.vr, 1, FUN=function(x)x-x.bar))))
    }, x.vr=list.x.p, x.bar=list.x.bar)
    
    return(score_out)
  }
  
  # Calculating and aggregating the score contributions across strata
  rowSums(do.call(cbind, mapply(function(u.s.star, tau.s, x.s, w.s, risk.vr, weights.s) score_s(beta, u.s.star, tau.s, x.s, w.s, risk.vr, weights.s),
                                u.s.star=u.star, tau.s=tau, x.s=x, w.s=w, risk.vr=risk, weights.s=weights, SIMPLIFY=FALSE)))
}

### Computes the mean covariate vector at each timepoint in stratum s
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
  # u.s.star        numeric             n.s-dimensional vector of u.sv.star for stratum s
  # tau.s           numeric             rho.s-dimensional vector of sorted unique left and right endpoints for stratum s
  # x.s             array               covariate array so that x.s[v, , ] returns a p x rho.s matrix, x.s[ ,p, ] an n.s x rho.s matrix, and x.s[,,r] an n.s x p matrix
  # risk.vr         matrix              matrix (n.s x rho.s) indicating risk set membership for each individual at each timepoint in stratum s
  # weights.s       numeric             n.s-dimensional vector of perturbation weights for stratum s
x_bar_s <- function(beta, u.s.star, tau.s, x.s, risk.vr, weights.s){
  
  # Returns a matrix with rows corresponding to individuals and columns to timepoints
  constant.vr <- do.call(rbind, lapply(seq_len(dim(x.s)[1]), FUN=function(x) exp(beta %*% x.s[x, , ])))
  
  # Returns a matrix with rows corresponding to timepoints and columns to covariates
  out <- matrix(NA, nrow=length(tau.s), ncol=length(beta))
  for (p in 1:length(beta)){
    x.vr <- x.s[ , p, ]
    out[, p] <- colSums(risk.vr * x.vr * constant.vr * weights.s)/colSums(risk.vr * constant.vr * weights.s)
  }
  
  return(out)
}

### Updates all baseline hazard components
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
  # x               list                list of length S containing x.s for each stratum
  # w               list                list of length S containing w.s for each stratum
  # risk            list                list of length S containing risk.vr for each stratum  
  # weights         list                list of length S containing the perturbation weights weight.s for each stratum
update_lambda <- function(beta, x, w, risk, weights=NULL){
  
  if (is.null(weights)){
    weights <- lapply(x, FUN=function(x.s) rep(1, dim(x.s)[1]))
  }
  
  # Updating the estimated baseline hazard function for each stratum and aggregating as a list
  mapply(function(x.s, w.s, risk.vr, weights.s){
    update_lambda_s(beta, x.s, w.s, risk.vr, weights.s)
  }, x.s=x, w.s=w, risk.vr=risk, weights.s=weights, SIMPLIFY=FALSE)
  
}

### Updates the baseline hazard components for stratum s
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
  # x.s             array               covariate array so that x.s[v, , ] returns a p x rho.s matrix, x.s[ ,p, ] an n.s x rho.s matrix, and x.s[,,r] an n.s x p matrix
  # w.s             matrix              matrix (n.s x rho.s) from the expectation step for stratum s
  # risk.vr         matrix              matrix (n.s x rho.s) indicating risk set membership for each individual at each timepoint in stratum s
  # weights.s       numeric             n.s-dimensional vector of perturbation weights for stratum s
update_lambda_s <- function(beta, x.s, w.s, risk.vr, weights.s){
  
  # Calculating the denominator at each timepoint
  constant.vr <- do.call(rbind, lapply(seq_len(dim(x.s)[1]), FUN=function(x) exp(beta %*% x.s[x, , ])))
  denom <- colSums(risk.vr * constant.vr * weights.s)
  
  # Calculating the numerator at each timepoint
  num <- colSums(risk.vr * w.s * weights.s)
  
  return(num/denom)
  
}


#-----------------------------------------------------------------------------------------------------
# Composite Log-Likelihood Functions
#-----------------------------------------------------------------------------------------------------

### Evaluates the independence composite log-likelihood function
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
  # lambda          list                list of length S containing the current estimate of lambda.s for each stratum
  # l               list                list of length S containing l.s for each stratum
  # u               list                list of length S containing u.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # mapping         list                optional list mapping the (i, j) indexing of cluster assignment (for robust variance estimation) to the (s, v) indexing of stratum assignment (for point estimation)
composite <- function(beta, lambda, l, u, tau, x, mapping=NULL){
  
  # If no mapping is provided, the clustering and stratification variable(s) is(are) assumed to be the same
  if (is.null(mapping)){
    do.call(sum, mapply(FUN=function(lambda.i, l.i, u.i, tau.i, x.i) composite_i_1(beta, lambda.i, l.i, u.i, tau.i, x.i),
                        lambda.i=lambda, l.i=l, u.i=u, tau.i=tau, x.i=x, SIMPLIFY=FALSE))
    
  # Otherwise, the mapping between cluster and stratum memberships are used to construct the composite log-likelihood
  } else{
    do.call(sum, mapply(FUN=function(stratum.i, indiv.i) composite_i_2(beta, lambda, l, u, tau, x, stratum.i, indiv.i),
                        stratum.i=mapping$stratum, indiv.i=mapping$indiv, SIMPLIFY=FALSE))
  }
  
}

### Calculates the independence composite log-likelihood contribution of cluster i (clustering and stratification variable(s) is(are) the same)
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
  # lambda.i        numeric             rho.i-dimensional vector of the current estimated baseline hazard function for cluster (stratum) i
  # l.i             numeric             n.i-dimensional vector of left endpoints for cluster (stratum) i
  # u.i             numeric             n.i-dimensional vector of right endpoints for cluster (stratum) i
  # tau.i           numeric             rho.i-dimensional vector of sorted unique left and right endpoints for cluster (stratum) i
  # x.i             array               covariate array so that x.i[j, , ] returns a p x rho.i matrix, x.i[ , p, ] an n.i x rho.i matrix, and x.i[ , , r] an n.i x p matrix
composite_i_1 <- function(beta, lambda.i, l.i, u.i, tau.i, x.i){
  
  # Computing and aggregating the individual contributions to the cluster-specific composite likelihood
  do.call(sum, mapply(FUN=function(l.ij, u.ij, x.pr) composite_ij(l.ij, u.ij, tau.i, x.pr, lambda.i, beta),
                      l.ij=l.i, u.ij=u.i, x.pr=lapply(seq_len(dim(x.i)[1]), FUN=function(x)matrix(x.i[x, , ], nrow=dim(x.i)[2])),
                      SIMPLIFY=FALSE))
  
}

### Calculates the independence composite log-likelihood contribution of cluster i (clustering and stratification variable(s) differ)
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
  # lambda          list                list of length S containing the current estimate of lambda.s for each stratum
  # l               list                list of length S containing l.s for each stratum
  # u               list                list of length S containing u.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # stratum.i       numeric             n.i-dimensional vector mapping subject j in cluster i to the index s of their corresponding stratum
  # indiv.i         numeric             n.i-dimensional vector mapping subject j in cluster i to their corresponding subject index v in stratum s
composite_i_2 <- function(beta, lambda, l, u, tau, x, stratum.i, indiv.i){
  
  # Computing and aggregating the individual contributions to the cluster-specific composite likelihood
  do.call(sum, mapply(FUN=function(stratum.ij, indiv.ij){
    tau.s <- tau[[stratum.ij]]; lambda.s <- lambda[[stratum.ij]]
    l.ij <- l[[stratum.ij]][indiv.ij]; u.ij <- u[[stratum.ij]][indiv.ij]
    x.pr <- adrop(x[[stratum.ij]][indiv.ij, , , drop=FALSE], drop=1)
    composite_ij(l.ij, u.ij, tau.s, x.pr, lambda.s, beta) 
  }, stratum.ij=stratum.i, indiv.ij=indiv.i, SIMPLIFY=FALSE))
  
}

### Calculates the independence composite log-likelihood contribution of individual j in cluster i
  # l.ij            numeric             the left endpoint of the censoring interval for individual j in cluster i
  # u.ij            numeric             the right-endpoint of the censoring interval for individual j in cluster i
  # tau.s           numeric             rho.s-dimensional vector of sorted unique left and right endpoints for stratum s
  # x.pr            matrix              matrix (p x rho.s) containing the values of each covariate at each timepoint for individual j in cluster i
  # lambda.s        numeric             rho.s-dimensional vector of the current estimated baseline hazard function for stratum s
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
composite_ij <- function(l.ij, u.ij, tau.s, x.pr, lambda.s, beta){
  
  # Computing the first term of the individual log-likelihood
  r.ids.1 <- (tau.s <= l.ij)
  term.1 <- -sum(lambda.s[r.ids.1] * exp(beta %*% x.pr[, r.ids.1]))
  
  # Computing the second term of the individual log-likelihood
  if (u.ij < Inf){
    r.ids.2 <- (tau.s > l.ij & tau.s <= u.ij)
    term.2 <- log(1-exp(-sum(lambda.s[r.ids.2] * exp(beta %*% x.pr[, r.ids.2]))))
  } else{
    term.2 <- 0
  }
  
  # Individual composite log-likelihood contribution
  term.1 + term.2
}


#-----------------------------------------------------------------------------------------------------
# Profile Composite Log-Likelihood Function
#-----------------------------------------------------------------------------------------------------

### Returns the estimated baseline hazard function(s) that maximize(s) the profile composite log-likelihood
  # beta            numeric             p-dimensional vector of the current fixed effect estimates
  # l               list                list of length S containing l.s for each stratum
  # u               list                list of length S containing u.s for each stratum
  # u.star          list                list of length S containing u.star.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # mapping         list                optional list mapping the (i, j) indexing of cluster assignment (for robust variance estimation) to the (s, v) indexing of stratum assignment (for point estimation)
  # control         list                control parameters for the composite EM algorithm and variance estimation
argmax_profile <- function(beta, l, u, u.star, tau, x, risk, mapping=NULL, control){
  
  # Initial values for the baseline hazard and composite likelihood
  lambda <- lapply(tau, FUN=function(x)rep(1/length(x), length(x)))
  cl <- composite(beta, lambda, l, u, tau, x, mapping)
  
  # Obtaining profile maximum likelihood estimators for lambda given the specified value of beta
  l1.norm <- 1; n.iter <- 0
  while (l1.norm > control$tol & n.iter < control$maxit.var){
    w <- calc_w(l, u, tau, x, lambda, beta)
    lambda.iter <- update_lambda(beta, x, w, risk)
    cl.iter <- composite(beta, lambda.iter, l, u, tau, x, mapping)
    l1.norm <- abs(cl.iter - cl)
    lambda <- lambda.iter
    cl <- cl.iter
    n.iter <- n.iter + 1
  }
  
  return(lambda)
}


#-----------------------------------------------------------------------------------------------------
# Variance Estimation Functions: Profile Likelihood
#-----------------------------------------------------------------------------------------------------

### Constructs the kth canonical vector in R^p
  # k               integer             indicates which canonical basis vector to construct
  # p               integer             dimension of the canonical basis
make_basis <- function(k, p) replace(numeric(p), k, 1)

### Returns the estimated (either model or robust) profile composite likelihood covariance matrix for the fixed effects estimators
  # beta            numeric             p-dimensional vector of the final fixed effect estimates
  # processed       list                processed data object used for point estimation
  # mapping         list                mapping of the (i, j) indexing of cluster assignment (for robust variance estimation) to the (s, v) indexing of stratum assignment (for point estimation)
  # control         list                control parameters for the composite EM algorithm and variance estimation
variance_beta <- function(beta, processed, mapping, control){
  
  # Determining the number of parameters
  p <- length(beta)
  n <- do.call(sum, lapply(processed$x, FUN=function(z)dim(z)[1]))
  h.n <- control$c/sqrt(n)
  
  # Computing profile composite likelihood maximizers
  lambda.list <- list(CMLE=argmax_profile(beta, processed$l, processed$u, processed$u.star, processed$tau, processed$x, processed$risk, control=control))
  lambda.list$first.order <- vector("list", p)
  for (j in 1:p){
    lambda.list$first.order[[j]] <- argmax_profile(beta + h.n*make_basis(j, p), processed$l, processed$u, processed$u.star, processed$tau, processed$x, processed$risk, control=control)
    lambda.list$second.order[[j]] <- vector("list", j)
    for (k in 1:j){
      lambda.list$second.order[[j]][[k]] <- argmax_profile(beta + h.n*make_basis(j, p) + h.n*make_basis(k, p), processed$l, processed$u, processed$u.star, processed$tau, processed$x, processed$risk, control=control)
    }
  }
  
  if (is.null(mapping) & length(processed$l)==1){
    
    # Computing the Hessian matrix
    H.hat <- H(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)
    
    # Returning the model-based covariance matrix
    return(solve(-H.hat))
    
  } else if (is.null(mapping) & length(processed$l) > 1){
    
    # Computing the Hessian matrix
    H.hat <- H(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)
    
    # Computing the covariance matrix for the profile composite likelihood score
    J.hat <- J(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)
    
    # Computing the inverse robust sandwich (Godambe) information
    H.hat.inv <- solve(H.hat)
    G.hat.inv <- H.hat.inv %*% J.hat %*% H.hat.inv
    
    # Returning the robust sandwich covariance matrix
    return(G.hat.inv)
    
  } else if (length(mapping$stratum)==1){
    
    # Computing the Hessian matrix
    H.hat <- H(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)
    
    # Returning the model-based covariance matrix
    return(solve(-H.hat))
    
  } else {
    
    # Computing the Hessian matrix
    H.hat <- H(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)
    
    # Computing the covariance matrix for the profile composite likelihood score
    J.hat <- J(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)
    
    # Computing the inverse robust sandwich (Godambe) information
    H.hat.inv <- solve(H.hat)
    G.hat.inv <- H.hat.inv %*% J.hat %*% H.hat.inv
    
    # Returning the robust sandwich covariance matrix
    return(G.hat.inv)
    
  }
  
}

### Computes the full-data sensitivity matrix
  # beta            numeric             p-dimensional vector of the final fixed effect estimates
  # lambda.list     list                profile maximum likelihood estimators for lambda under various first- and second-order perturbations of beta
  # l               list                list of length S containing l.s for each stratum
  # u               list                list of length S containing u.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # mapping         list                mapping of the (i, j) indexing of cluster assignment (for robust variance estimation) to the (s, v) indexing of stratum assignment (for point estimation)
  # p               integer             dimension of the coefficient vector
  # h.n             integer             perturbation constant for the finite differences approximations
H <- function(beta, lambda.list, l, u, tau, x, mapping, p, h.n){
  
  # Initializing output matrix
  out <- matrix(NA, nrow=p, ncol=p)
  
  # Computing each element of the Hessian matrix
  #   Can improve this coding to avoid multiple computations of the same composite likelihood value
  for (j in 1:p){
    for (k in 1:j){
      temp <- (composite(beta, lambda.list$CMLE, l, u, tau, x, mapping) - composite(beta + h.n*make_basis(j, p), lambda.list$first.order[[j]], l, u, tau, x, mapping) - composite(beta + h.n*make_basis(k, p), lambda.list$first.order[[k]], l, u, tau, x, mapping) + composite(beta + h.n*make_basis(j, p) + h.n*make_basis(k, p), lambda.list$second.order[[j]][[k]], l, u, tau, x, mapping))/h.n^2
      out[j, k] <- out[k, j] <- temp
    }
  }
  
  return(out)
}

### Computes the full-data variability matrix
  # beta            numeric             p-dimensional vector of the final fixed effect estimates
  # lambda.list     list                profile maximum likelihood estimators for lambda under various first- and second-order perturbations of beta
  # l               list                list of length S containing l.s for each stratum
  # u               list                list of length S containing u.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # mapping         list                mapping of the (i, j) indexing of cluster assignment (for robust variance estimation) to the (s, v) indexing of stratum assignment (for point estimation)
  # p               integer             dimension of the coefficient vector
  # h.n             integer             perturbation constant for the finite differences approximations
J <- function(beta, lambda.list, l, u, tau, x, mapping, p, h.n){
  
  if (is.null(mapping)){
    
    # Restructuring lambda.list for use with mapply
    lambda <- vector("list", length(x))
    for (i in 1:length(x)){
      temp1 <- lambda.list$CMLE[[i]]
      temp2 <- lapply(seq_len(length(lambda.list$first.order)), FUN=function(z) lambda.list$first.order[[z]][[i]])
      lambda[[i]] <- list(CMLE=temp1, first.order=temp2)
    }
    
    # Summing across cluster-specific outer products (cluster-specific contributions to the full-data variability matrix)
    Reduce("+", mapply(FUN=function(lambda.list.i, l.i, u.i, tau.i, x.i){
      J_i_1(beta, lambda.list.i, l.i, u.i, tau.i, x.i, p, h.n)},
      lambda.list.i=lambda, l.i=l, u.i=u, tau.i=tau, x.i=x, SIMPLIFY=FALSE))
    
  } else {
    
    # Summing across cluster-specific outer products (cluster-specific contributions to the full-data variability matrix)
    Reduce("+", mapply(FUN=function(stratum.i, indiv.i){
      J_i_2(beta, lambda.list, l, u, tau, x, stratum.i, indiv.i, p, h.n)},
      stratum.i=mapping$stratum, indiv.i=mapping$indiv, SIMPLIFY=FALSE))
    
  }

}

### Computes the variability matrix contribution from cluster i (clustering and stratification variable(s) is(are) the same)
  # beta            numeric             p-dimensional vector of the final fixed effect estimates
  # lambda.list.i   list                profile maximum likelihood estimators for lambda under various first- and second-order perturbations of beta for cluster (stratum) i
  # l.i             numeric             n.i-dimensional vector of left endpoints for cluster (stratum) i
  # u.i             numeric             n.i-dimensional vector of right endpoints for cluster (stratum) i
  # tau.i           numeric             rho.i-dimensional vector of sorted unique left and right endpoints for cluster (stratum) i
  # x.i             array               covariate array so that x.i[j, , ] returns a p x rho.i matrix, x.i[ , p, ] an n.i x rho.i matrix, and x.i[ , , r] an n.i x p matrix
  # p               integer             dimension of the coefficient vector
  # h.n             integer             perturbation constant for the finite differences approximations
J_i_1 <- function(beta, lambda.list.i, l.i, u.i, tau.i, x.i, p, h.n){
  
  # Profile likelihood contribution under the composite maximum likelihood estimator
  pl.cmle <- composite_i_1(beta, lambda.list.i$CMLE, l.i, u.i, tau.i, x.i)
  
  # Computing cluster-specific gradient
  score_i <- matrix(NA, nrow=p, ncol=1)
  for (j in 1:p){
    score_i[j,] <- (composite_i_1(beta + h.n*make_basis(j, p), lambda.list.i$first.order[[j]], l.i, u.i, tau.i, x.i) - pl.cmle)/h.n
  }
  
  score_i %*% t(score_i)  
}

### Computes the variability matrix contribution from cluster i (clustering and stratification variable(s) differ)
  # beta            numeric             p-dimensional vector of the final fixed effect estimates
  # lambda.list     list                profile maximum likelihood estimators for lambda under various first- and second-order perturbations of beta
  # l               list                list of length S containing l.s for each stratum
  # u               list                list of length S containing u.s for each stratum
  # tau             list                list of length S containing tau.s for each stratum
  # x               list                list of length S containing x.s for each stratum
  # stratum.i       numeric             n.i-dimensional vector mapping subject j in cluster i to the index s of their corresponding stratum
  # indiv.i         numeric             n.i-dimensional vector mapping subject j in cluster i to their corresponding subject index v in stratum s
  # p               integer             dimension of the coefficient vector
  # h.n             integer             perturbation constant for the finite differences approximations
J_i_2 <- function(beta, lambda.list, l, u, tau, x, stratum.i, indiv.i, p, h.n){
  
  # Profile likelihood contribution under the composite maximum likelihood estimator
  pl.cmle <- composite_i_2(beta, lambda.list$CMLE, l, u, tau, x, stratum.i, indiv.i)
  
  # Computing cluster-specific gradient
  score.i <- matrix(NA, nrow=p, ncol=1)
  for (j in 1:p){
    score.i[j,] <- (composite_i_2(beta + h.n*make_basis(j, p), lambda.list$first.order[[j]], l, u, tau, x, stratum.i, indiv.i) - pl.cmle)/h.n
  }
  
  score.i %*% t(score.i)  
}
