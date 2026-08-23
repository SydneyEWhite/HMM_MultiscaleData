##### Libraries
library(moveHMM)
library(tidyr)
library(circular)
library(cluster)
library(gtools)
library(combinat)
library(future)
library(furrr)
library(dplyr)


set.seed(1)

################################################## Define Parameters ################################################## 
# 1: Travel
# 2: Territory / hunting (Marking trees and waterholes)
# 3: Rest

# fine scale: data collected 1 Hz --> once per second 
# coarse scale: data collected every 10 minutes

N <- 3 # number of states
coarse_gamma <- matrix(c(0.602, 0.119, 0.279, 0.364, 0.418, 0.218, 0.383, 0.081, 0.536), nrow = 3, byrow = TRUE)
# solve for stationary distribution 
delta <- solve(t(diag(3) - coarse_gamma + 1), rep(1, 3))

# coarse scale states
num_coarse_states <-  900 # 30 mins # 1800 # 1 hour #300 # 10 mins
S <- numeric(num_coarse_states) # state holder
S[1] <- sample(c(1,2,3), 1, prob = delta) # initial state based on initial distribution
for (i in 2:num_coarse_states) {
  S[i] <- sample(c(1,2,3), 1, prob = coarse_gamma[S[i-1],])
}
# fine scale states
fine_scale_points <- 50
s <- rep(S, each = fine_scale_points)
df <- data.frame(s)
### parameters
## step parameters: normal, hard
stepPar_mean <- c(10,0.1,0.1) # mean1, mean2, mean3,
stepPar_std <- c(5,0.1,0.1) # sd1, sd2, sd3; meters per second
## step parameters: easy
# stepPar_mean <- c(10,5,0.1) # mean1, mean2, mean3,
# stepPar_std <- c(5,3,0.1) # sd1, sd2, sd3; meters per second
## angle parameters
anglePar_mean <- c(0,0,0) # mean1, mean2, mean3,
anglePar_k <- c(4,0.8,0.2) # k1, k2, k3 -> k is like concentration around the mean
## acceleration parameters: easy, normal
dba_mean <- c(0.25, 0.45, 0.12)
dba_std   <- c(0.2, 0.3, 0.1)
## acceleration parameters: hard
# dba_mean <- c(0.25, 0.25, 0.12)
# dba_std   <- c(0.2, 0.2, 0.1)
# fill in parameters based on state
df$stepmean <- stepPar_mean[df$s]
df$stepstd <- stepPar_std[df$s]
df$anglemean <- anglePar_mean[df$s]
df$anglek <- anglePar_k[df$s]
df$dbamean <- dba_mean[df$s]
df$dbastd <- dba_std[df$s]

################################################## Probability and Likelihood Functions ################################################## 

unpack_pars <- function(theta.star, N, scenario, k = NULL) {
  ### create initial parameter values based on theta.star input
  # note we take the exponentials after inputting logs because it's easier to have unrestricted optimization. The exponential, then restricts the parameters to positive values
  # transition probability matrix
  Gamma <- diag(N) # creates a diagonal matrix of dim N with 1s on the diagonals 
  Gamma[!Gamma] <- exp(theta.star[1:((N-1) * N)]) # fills the off-diagonal elements with exponentiated gamma estimate values (ensure positivity)
  Gamma <- Gamma / rowSums(Gamma) # for each row divide the entries by the row sum to obtain row sums equal to one
  # stationary distribution
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  
  ###### not common but when solving for mean summary with only 10 minutes of time, there ended up being a negative delta value -5.7e-17 which caused an NaN in the log. This appears to just be noise rather than an actual issue so I added the floor for noise 
  if (any(delta < 0)) {
    if (any(delta < -1e-8)) {
      warning("non-negligible negative delta value. Delta = ", delta)
    }
    delta[delta < 0] <- 0
    delta <- delta / sum(delta)
  }
  
  # distribution parameters
  base <- (N - 1) * N
  out <- list(Gamma = Gamma, delta = delta)
  
  if (!scenario %in% c('base_acc')) {
    deltas <- exp(theta.star[base + 1:N]) 
    out$mu_step <- rev(cumsum(rev(deltas)))
    
    out$sig_step <- exp(theta.star[base + N + 1:N])
    out$mu_angle <- theta.star[base + 2*N + 1:N]        # not exponentiated, can be negative
    kappa <- exp(theta.star[base + 3*N + 1:N])   # must be positive
    out$kappa <- pmax(pmin(kappa, 1e4), 1e-6) # stops extreme k values
  }
  
  
  if (scenario %in% c('ideal', 'coarse')) {
    out$mu_acc <- exp(theta.star[base + 4*N + 1:N])
    out$sd_acc <- exp(theta.star[base + 5*N + 1:N])
  } else if (scenario %in% c('base_acc')) {
    deltas <- exp(theta.star[base + 1:N])
    out$mu_acc <- rev(cumsum(rev(deltas)))
    out$sd_acc <- exp(theta.star[base + N + 1:N])
  } else if (scenario == 'mean_summary') {
    out$mu_acc <- exp(theta.star[base + 4*N + 1:N])
    sd_acc <- exp(theta.star[base + 5*N + 1:N])
    out$sd_acc_mean <- sd_acc
  } else if (scenario == 'var_summary') {
    out$mu_acc_var <- exp(theta.star[base + 4*N + 1:N])
    out$sd_acc_var <- exp(theta.star[base + 5*N + 1:N])
  } else if (scenario == 'med_summary') {
    out$mu_acc <- exp(theta.star[base + 4*N + 1:N])
    sd_acc <- exp(theta.star[base + 5*N + 1:N])
    out$sd_acc_med <- sd_acc #* sqrt(pi/2) 
  } else if (scenario == 'stats_summary') {
    out$mu_acc <- theta.star[base + 4*N + 1:N] # unrestricted since normalized
    out$sd_acc <- exp(theta.star[base + 5*N + 1:N])
    out$sd_acc_mean <- out$sd_acc
    
    out$mu_acc_var <- theta.star[base + 6*N + 1:N] # unrestricted since normalized
    out$sd_acc_var <- exp(theta.star[base + 7*N + 1:N])
    
    out$mu_acc_med <- theta.star[base + 8*N + 1:N]
    out$sd_acc_med <- exp(theta.star[base + 9*N + 1:N])
  } else if (scenario == 'coarse_abs') {
    out$mu_acc <- exp(theta.star[base + 4*N + 1:N]) # exp because it must be positive
    out$sd_acc <- exp(theta.star[base + 5*N + 1:N])
  } else if (scenario %in% c('density_summary')) {
    K <- k  # number of density vector elements
    # create a mu (mean) and phi (concentration) which can be used later to recover alpha
    # alpha = (alpha_1, alpha_2, ...) is a vector of positive numbers that parameterizes the Dirichlet distribution
    # mu = alpha / sum(alpha)
    # phi = sum(alpha)
    # therefore alpha = phi * mu
    # since the rows have to sum to 1, there are only K-1 free parameters (the last can be calculated from the first 6)
    mu_raw <- matrix(theta.star[base + 4*N + 1:((K-1)*N)], nrow = N, ncol = K-1, byrow = TRUE)
    # additive log ratio transform
    # to keep the optimization unconstrained, we are actually estimating the log ratio values of mu
    # each of the 6 raw_mu are the log ratio of the mu_comp_i to mu_comp_7
    # since log(mu_comp_7/mu_comp_7) = log(1) = 0, we create a column of 0s
    mu_raw <- cbind(mu_raw, 0)  
    # now we exponentiate to get just the ratio of mu_comp_i to mu_comp_7
    mu_comp <- exp(mu_raw)
    # now normalize so each row sums to one
    out$mu_comp <- mu_comp / rowSums(mu_comp)
    
    # concentration parameter per state (controls spread around mu_j)
    out$phi <- exp(theta.star[base + 4*N + (K-1)*N + 1:N])
    out$alpha <- out$phi * out$mu_comp
  } else if  (scenario %in% c('quantile_summary')){
    K <- k
    quant_cols <- 3:(3+K-1)
    
    # model each quantile with a gamma distribution
    out$mu_q <- matrix(exp(theta.star[base + 4*N + 1:(K*N)]), nrow = N, ncol = K, byrow = TRUE)
    out$sd_q <- matrix(exp(theta.star[base + 4*N + K*N + 1:(K*N)]), nrow = N, ncol = K, byrow = TRUE)
  }
  out
}

### function to calculate probabilities of each obs being in each state
get_allprobs <- function(theta.star, x, N, scenario, k = NULL) {
  p <- unpack_pars(theta.star, N, scenario, k = k)
  Gamma <- p$Gamma
  delta <- p$delta

  if (!scenario %in% c('base_acc')) {
    mu_step <- p$mu_step
    sig_step <- p$sig_step
    mu_angle <- p$mu_angle
    kappa <- p$kappa
  }
  
  if (scenario %in% c('ideal', 'coarse')) {
    mu_acc <- p$mu_acc
    sd_acc <- p$sd_acc
  } else if (scenario %in% c('base_acc')) {
    mu_acc <- p$mu_acc
    sd_acc <- p$sd_acc
  } else if (scenario == 'mean_summary') {
    mu_acc <- p$mu_acc
    sd_acc <- p$sd_acc
    sd_acc_mean <- sd_acc
  } else if (scenario == 'var_summary') {
    mu_acc_var <- p$mu_acc_var
    sd_acc_var <- p$sd_acc_var
  } else if (scenario == 'med_summary') {
    mu_acc <- p$mu_acc
    sd_acc <- p$sd_acc
    sd_acc_med <- sd_acc 
  } else if (scenario == 'stats_summary') {
    mu_acc <- p$mu_acc
    sd_acc <- p$sd_acc
    sd_acc_mean <- sd_acc
    
    mu_acc_var <- p$mu_acc_var
    sd_acc_var <- p$sd_acc_var
    
    mu_acc_med <- p$mu_acc_med
    sd_acc_med <- p$sd_acc_med
  } else if (scenario == 'coarse_abs') {
    mu_acc <- p$mu_acc
    sd_acc <- p$sd_acc
  } else if (scenario %in% c('density_summary')) {
    K <- k  
    mu_comp <- p$mu_comp
    phi <- p$phi
  } else if  (scenario %in% c('quantile_summary')){
    K <- k
    mu_q <- p$mu_q
    sd_q <- p$sd_q
  }
  
  
  ### initialization
  n_obs <- nrow(x)
  log_allprobs <- matrix(0, n_obs, N) # T x N matrix of 0s (because log(1) = 0, and we are working in log space)
  
  ### iterate though the log probability of each observation (x) being in each state (j)
  for (j in 1:N) {
    
    if (!scenario %in% c('base_acc')) {
      # find indices of non-NA values
      ind_step  <- which(!is.na(x[, 1]))
      ind_angle <- which(!is.na(x[, 2]))
      # find the probability of state j given step lengths (ind_step)
      log_allprobs[ind_step, j]  <- log_allprobs[ind_step, j] + dgamma(x[ind_step, 1], shape = mu_step[j]^2 / sig_step[j]^2, scale = sig_step[j]^2 / mu_step[j], log = TRUE)   # log = TRUE is more stable in case values are close to 0
      # find the probability of state j given turning angle (ind_ang)
      log_allprobs[ind_angle, j] <- log_allprobs[ind_angle, j] + dvonmises(circular(x[ind_angle, 2], units = 'radians'), mu = circular(mu_angle[j], units = 'radians'), kappa = kappa[j], log = TRUE)
    }
    
    if (scenario %in% c('ideal', 'coarse')) {
      ind_acc <- which(!is.na(x[, 3]))
      # find the probability of state j given acceleration (ind_acc)
      log_allprobs[ind_acc, j] <- log_allprobs[ind_acc, j] + dgamma(x[ind_acc, 3], shape = mu_acc[j]^2 / sd_acc[j]^2, scale = sd_acc[j]^2 / mu_acc[j], log = TRUE)
    } else if (scenario %in% c('base_acc')) {
      ind_acc <- which(!is.na(x[, 1]))
      log_allprobs[ind_acc, j] <- log_allprobs[ind_acc, j] + dgamma(x[ind_acc, 1], shape = mu_acc[j]^2 / sd_acc[j]^2, scale = sd_acc[j]^2 / mu_acc[j], log = TRUE)
    } else if (scenario %in% c('mean_summary')) {
      ind_acc <- which(!is.na(x[, 3]))
      # find the probability of state j given acceleration mean (ind_acc)
      log_allprobs[ind_acc, j] <- log_allprobs[ind_acc, j] + dgamma(x[ind_acc, 3], shape = mu_acc[j]^2 / sd_acc_mean[j]^2, scale = sd_acc_mean[j]^2 / mu_acc[j], log = TRUE)
    } else if (scenario %in% c('med_summary')) {
      ind_acc <- which(!is.na(x[, 3]))
      log_allprobs[ind_acc, j] <- log_allprobs[ind_acc, j] + dgamma(x[ind_acc, 3], shape = mu_acc[j]^2 / sd_acc_med[j]^2, scale = sd_acc_med[j]^2 / mu_acc[j], log = TRUE)
    } else if (scenario == 'var_summary') {
      ind_acc_var <- which(!is.na(x[, 3]))
      log_allprobs[ind_acc_var, j] <- log_allprobs[ind_acc_var, j] + dgamma(x[ind_acc_var, 3], shape = mu_acc_var[j]^2 / sd_acc_var[j]^2, scale = sd_acc_var[j]^2 / mu_acc_var[j], log = TRUE)
    } else if (scenario %in% c('stats_summary')) {
      ind_acc <- which(!is.na(x[, 3]))
      ind_acc_var <- which(!is.na(x[, 4]))
      ind_acc_med <- which(!is.na(x[, 5]))
      
      # find the probability of state j given acceleration mean (ind_acc)
      log_allprobs[ind_acc, j] <- log_allprobs[ind_acc, j] + dnorm(x[ind_acc, 3], mean = mu_acc[j], sd = sd_acc[j], log = TRUE)

      # find the probability of state j given acceleration variance (ind_acc)
      log_allprobs[ind_acc_var, j] <- log_allprobs[ind_acc_var, j] + dnorm(x[ind_acc_var, 4], mean = mu_acc_var[j], sd = sd_acc_var[j], log = TRUE)
      
      # find the probability of state j given acceleration median (ind_acc)
      log_allprobs[ind_acc_med, j] <- log_allprobs[ind_acc_med, j] + dnorm(x[ind_acc_med, 5], mean = mu_acc_med[j], sd = sd_acc_med[j], log = TRUE)
    } else if (scenario %in% c('coarse_abs')) {
      ind_acc <- which(!is.na(x[, 3]))
      
      log_allprobs[ind_acc, j] <- log_allprobs[ind_acc, j] + dgamma(x[ind_acc, 3], shape = mu_acc[j]^2 / sd_acc[j]^2, scale = sd_acc[j]^2 / mu_acc[j], log = TRUE) # gamma distributed because absolute value
    } else if (scenario %in% c('density_summary')) {
      K <- k
      dens_cols <- 3:(3 + K - 1)  # x columns need to be step, angle, dens_1...dens_K
      ind_dens <- which(rowSums(is.na(x[, dens_cols])) == 0) # no column is na
      
      alpha_j <- phi[j] * mu_comp[j, ]  # calculate Dirichlet shape parameters for state j
      
      log_allprobs[ind_dens, j] <- log_allprobs[ind_dens, j] + log(gtools::ddirichlet(x[ind_dens, dens_cols, drop = FALSE], alpha_j)) # drop = False keeps it as a matrix if for some reason there would only be one column (which shouldn't happen)
      #browser()
    } else if  (scenario %in% c('quantile_summary')){
      K <- k
      dens_cols <- 3:(3 + K - 1)
      ind_dens <- which(rowSums(is.na(x[, dens_cols, drop = FALSE])) == 0)
      
      for (q in 1:K) {
        log_allprobs[ind_dens, j] <- log_allprobs[ind_dens, j] + dgamma(x[ind_dens, 2 + q], shape = mu_q[j, q]^2 / sd_q[j, q]^2, scale = sd_q[j, q]^2 / mu_q[j, q], log   = TRUE)
      }
    }
  }

  list(
    allprobs = exp(log_allprobs), # move from log probabilities to actual probabilities
    Gamma    = Gamma,
    delta    = delta
  ) 
}

### forward probability algorithm
mllk <- function(theta.star, x, N, scenario='base', k = NULL) {
  p <- get_allprobs(theta.star, x, N, scenario=scenario, k = k)
  n_obs <- nrow(x)
  
  # initiate forward probability
  alpha_1 <- p$delta %*% diag(p$allprobs[1, ]) # initial distribution * probability of seeing x given j
  l <- log(sum(alpha_1)) # log likelihood at t = 1
  phi <- alpha_1 / sum(alpha_1) # normalize
  # forward probability algorithm iteration
  for (t in 2:n_obs){
    alpha_t <- phi %*% p$Gamma %*% diag(p$allprobs[t, ])
    l <- l + log(sum(alpha_t)) # accumulated log likelihood (adding logs = multiplying probabilities)
    phi <- alpha_t / sum(alpha_t)
  }
  
  return(-l)
}

################################################## Parameter Initialization Functions ################################################## 
### generate clustering-based start
#possible_perm <- permutations(n = N, r = N)^

generate_cluster_start <- function(x, N, scenario, iter, k = NULL) {
  print("using cluster")
  if (scenario %in% c('ideal', 'coarse', 'mean_summary', 'med_summary', 'coarse_abs')) {
    cluster_data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), scale_acc = scale(x[,3]), step = x[,1], angle = x[,2], acc = x[,3]))
  } else if (scenario %in% c('base')){
    cluster_data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), step = x[,1], angle = x[,2]))
  } else if (scenario %in% c('var_summary')){
    cluster_data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), scale_acc_var = scale(x[,3]), step = x[,1], angle = x[,2], acc_var = x[,3]))
  } else if (scenario %in% c('stats_summary')){
    cluster_data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), scale_acc_mean = scale(x[,3]), scale_acc_var = scale(x[,4]), step = x[,1], angle = x[,2], acc_mean = x[,3], acc_var = x[,4], acc_med = x[,5]))
  } else if (scenario %in% c('density_summary', 'quantile_summary')) {
    K <- k
    dens_cols <- 3:(3 + K - 1)
    cluster_data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), step = x[,1], angle = x[,2], x[, dens_cols]))
    colnames(cluster_data)[5:(5+K-1)] <- paste0("dens_", 1:K)
  } else {
    print('which scenario?????')
  }
  
  set.seed(1)
  
  if (scenario %in% c('base', 'density_summary', 'quantile_summary')){
    km <- kmeans(cluster_data[, c('scale_step', 'scale_angle')], centers = 3, nstart = 100)
  } else if (scenario %in% c('var_summary')){
    km <- kmeans(cluster_data[, c('scale_step', 'scale_acc_var')], centers = 3, nstart = 100)
  } else if (scenario %in% c('stats_summary')){
    km <- kmeans(cluster_data[, c('scale_step', 'scale_acc_mean', 'scale_acc_var')], centers = 3, nstart = 100)
  } else {   
    km <- kmeans(cluster_data[, c('scale_step', 'scale_acc')], centers = 3, nstart = 100)
  }
  
  # order clusters by step length since states are travel, hunt, rest
  print(tapply(cluster_data$step, km$cluster, mean))
  cluster_order <- order(tapply(cluster_data$step, km$cluster, mean), decreasing = TRUE)
  
  #cluster_order <- possible_perm[iter,]
  cluster_labels <- match(km$cluster, cluster_order)
  
  # estimate starting values 
  init_step_mean <- tapply(cluster_data$step, cluster_labels, mean)
  init_step_sd <- tapply(cluster_data$step, cluster_labels, sd)
  init_angle_mean <- tapply(cluster_data$angle, cluster_labels, function(x) as.numeric(mean.circular(circular(as.numeric(x), type = 'angles', units = "radians"))))

  if (scenario %in% c('ideal', 'coarse', 'mean_summary', 'med_summary', 'coarse_abs')) {
    init_acc_mean <- pmax(tapply(cluster_data$acc, cluster_labels, mean), 1e-4) # ensuring no 0
    init_acc_sd <- tapply(cluster_data$acc, cluster_labels, sd)
    
    init_acc_sd <- pmax(init_acc_sd, 0.5 * sd(cluster_data$acc))   # Ensure sd values aren't too small
  } else if (scenario %in% c('var_summary')) {
    init_acc_var_mean <- tapply(cluster_data$acc_var, cluster_labels, mean)
    init_acc_var_sd <- tapply(cluster_data$acc_var, cluster_labels, sd)
    init_acc_var_sd <- pmax(init_acc_var_sd, 0.01)
  } else if (scenario %in% c('stats_summary')) {
    init_acc_mean <- pmax(tapply(cluster_data$acc_mean, cluster_labels, mean), 1e-4) # ensuring no 0
    init_acc_sd <- tapply(cluster_data$acc_mean, cluster_labels, sd)
    init_acc_sd <- pmax(init_acc_sd, 0.5 * sd(cluster_data$acc_mean))   # Ensure sd values aren't too small
    
    init_acc_var_mean <- tapply(cluster_data$acc_var, cluster_labels, mean)
    init_acc_var_sd <- tapply(cluster_data$acc_var, cluster_labels, sd)
    init_acc_var_sd <- pmax(init_acc_var_sd, 0.01)
    
    init_acc_med_mean <- tapply(cluster_data$acc_med, cluster_labels, mean)
    init_acc_med_sd <- tapply(cluster_data$acc_med, cluster_labels, sd)
    init_acc_med_sd <- pmax(init_acc_med_sd, 0.5 * sd(cluster_data$acc_med))
  } else if (scenario %in% c('density_summary')) {
    K <- k
    dens_mat <- as.matrix(cluster_data[, paste0("dens_", 1:K)])
    
    init_mu_comp <- t(sapply(1:N, function(j) {colMeans(dens_mat[cluster_labels == j, , drop = FALSE])}))
    init_mu_comp <- init_mu_comp / rowSums(init_mu_comp)  # ensure rows sum to 1
    
    # convert to additive-log-ratio scale relative to last component (component K as reference)
    init_mu_raw <- log(init_mu_comp[, 1:(K-1), drop = FALSE] / init_mu_comp[, K])
    
    # rough concentration estimate: inverse of average within-cluster variance -- note this is an approximation
    init_phi <- sapply(1:N, function(j) {
      # create matrix of just values in cluster j
      sub <- dens_mat[cluster_labels == j, , drop = FALSE]
      v <- mean(apply(sub, 2, var)) # mean of each columns variance 
      # phi is the inverse of variance (since concentration), and we add the small value to prevent division by 0
      max(1 / (v + 1e-6), 1)
    })
  } else if (scenario %in% c('quantile_summary')) {
    K <- k
    dens_mat <- as.matrix(cluster_data[, paste0("dens_", 1:K)])
    
    init_q_mean <- t(sapply(1:N, function(j) {colMeans(dens_mat[cluster_labels == j, , drop = FALSE])}))
    init_q_mean <- pmax(init_q_mean, 1e-6)
    
    init_q_sd <- t(sapply(1:N, function(j) {
      sub <- dens_mat[cluster_labels == j, , drop = FALSE]
      if (nrow(sub) < 2) return (rep(0.1, K))
      apply(sub, 2, sd) }))
    init_q_sd <- pmax(init_q_sd, 1e-6)
  }
  
  
  # Ensure sd values aren't too small
  init_step_sd <- pmax(init_step_sd, 0.5 * sd(cluster_data$step))
  
  if (scenario %in% c('ideal', 'coarse', 'mean_summary', 'med_summary')) {
    return(c(
      rep(-2, (N-1)*N),
      c(log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3]), log(init_step_mean[3])),
      #c(log(init_step_mean[1]), log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3])),
      log(init_step_sd),
      init_angle_mean,
      log(c(4, 0.8, 0.2)),
      log(init_acc_mean),
      log(init_acc_sd) 
    ))
  } else if (scenario %in% c('var_summary')) {
    return(c(
      rep(-2, (N-1)*N),
      c(log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3]), log(init_step_mean[3])),
      #c(log(init_step_mean[1]), log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3])),
      log(init_step_sd),
      init_angle_mean,
      log(c(4, 0.8, 0.2)),
      log(init_acc_var_mean),
      log(init_acc_var_sd) 
    ))
  } else if (scenario %in% c('base')){
    return(c(
      rep(-2, (N-1)*N),
      c(log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3]), log(init_step_mean[3])),
      #c(log(init_step_mean[1]), log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3])),
      log(init_step_sd),
      init_angle_mean,
      log(c(4, 0.8, 0.2))
    ))    
  } else if (scenario %in% c('stats_summary')) {
    return(c(
      rep(-2, (N-1)*N),
      c(log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3]), log(init_step_mean[3])),
      #c(log(init_step_mean[1]), log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3])),
      log(init_step_sd),
      init_angle_mean,
      log(c(4, 0.8, 0.2)),
      init_acc_mean,
      log(init_acc_sd), 
      init_acc_var_mean,
      log(init_acc_var_sd),
      init_acc_med_mean,
      log(init_acc_med_sd) 
    ))
  } else if (scenario %in% c('coarse_abs')) {
    return(c(
      rep(-2, (N-1)*N),
      c(log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3]), log(init_step_mean[3])),
      #c(log(init_step_mean[1]), log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3])),
      log(init_step_sd),
      init_angle_mean,
      log(c(4, 0.8, 0.2)),
      log(init_acc_mean),
      log(init_acc_sd) 
    ))
  } else if (scenario %in% c('density_summary')) {
    return(c(
      rep(-2, (N-1)*N),
      c(log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3]), log(init_step_mean[3])),
      log(init_step_sd),
      init_angle_mean,
      log(c(4, 0.8, 0.2)),
      as.vector(t(init_mu_raw)),  # flatten N x (K-1) matrix row-wise
      log(init_phi)
    ))
  } else if (scenario %in% c('quantile_summary')) {
    return(c(
      rep(-2, (N-1)*N),
      c(log(init_step_mean[1] - init_step_mean[2]), log(init_step_mean[2] - init_step_mean[3]), log(init_step_mean[3])),
      log(init_step_sd),
      init_angle_mean,
      log(c(4, 0.8, 0.2)),
      log(as.vector(t(init_q_mean))),
      log(as.vector(t(init_q_sd)))
    ))
  }
}

### generate random start
generate_random_start <- function(x, N, scenario, k = NULL) {
  print("using random number gen")
  if (scenario %in% c('ideal', 'coarse', 'mean_summary', 'med_summary', 'coarse_abs')) {
    data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), scale_acc = scale(x[,3]), step = x[,1], angle = x[,2], acc = x[,3]))
  } else if (scenario %in% c('base')){
    data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), step = x[,1], angle = x[,2]))
  } else if (scenario %in% c('base_acc')){
    data <- na.omit(data.frame(acc = x[,1]))
  } else if (scenario %in% c('var_summary')){
    data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), scale_acc_var = scale(x[,3]), step = x[,1], angle = x[,2], acc_var = x[,3]))
  } else if (scenario %in% c('stats_summary')){
    data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), scale_acc_mean = scale(x[,3]), scale_acc_var = scale(x[,4]), step = x[,1], angle = x[,2], acc_mean = x[,3], acc_var = x[,4], acc_med = x[,5]))
  } else if (scenario %in% c('density_summary', 'quantile_summary')) {
    K <- k
    dens_cols <- 3:(3 + K - 1)
    data <- na.omit(data.frame(scale_step = scale(x[,1]), scale_angle = scale(x[,2]), step = x[,1], angle = x[,2], x[, dens_cols]))
  } else {
    print('which scenario?????')
  }

  # use quantiles to inform random step ranges
  if (!scenario %in% c('base_acc')){
    step_range <- range(data$step)
    step_sd_range <- c(0.1, 2) * sd(data$step)
    
    sorted_means <- sort(runif(N, min = step_range[1], max= step_range[2]), decreasing = TRUE)
    ordered_step_params <- c(log(sorted_means[1] - sorted_means[2]), log(sorted_means[2] - sorted_means[3]), log(sorted_means[3])) #c(log(sorted_means[1]), log(sorted_means[1] - sorted_means[2]), log(sorted_means[2] - sorted_means[3]))
  }
  
  
  if (scenario %in% c('ideal', 'coarse', 'mean_summary', 'med_summary', 'coarse_abs')) {
    acc_range <- range(data$acc)
    acc_range[1] <- max(acc_range[1], 1e-4) # ensure no 0
    acc_sd_range <- c(0.1, 2) * sd(data$acc)
    
    return(c(
      rnorm((N-1)*N, mean = -2, sd = 1), # TPM: random around -2 (gives ~12% switching probability)
      ordered_step_params,
      log(runif(N, min = step_sd_range[1], max = step_sd_range[2])), # Step SDs: random but reasonable
      runif(N, -pi, pi), # Angle means: uniform over circle
      log(runif(N, 0.1, 5)), # Kappa: random but reasonable
      log(runif(N, min = acc_range[1], max = acc_range[2])),
      log(runif(N, min = acc_sd_range[1], max = acc_sd_range[2])) 
    ))
  } else if (scenario %in% c('var_summary')){
    acc_var_range <- range(data$acc_var)
    acc_var_sd_range <- c(0.1, 2) * sd(data$acc_var)
    
    return(c(
      rnorm((N-1)*N, mean = -2, sd = 1), # TPM: random around -2 (gives ~12% switching probability)
      ordered_step_params,
      log(runif(N, min = step_sd_range[1], max = step_sd_range[2])), # Step SDs: random but reasonable
      runif(N, -pi, pi), # Angle means: uniform over circle
      log(runif(N, 0.1, 5)), # Kappa: random but reasonable
      log(runif(N, min = acc_var_range[1], max = acc_var_range[2])), 
      log(runif(N, min = acc_var_sd_range[1], max = acc_var_sd_range[2])) 
    ))   
  } else if (scenario %in% c('stats_summary')) {
    acc_range <- range(data$acc_mean)
    acc_sd_range <- c(0.1, 2) * sd(data$acc_mean)
    
    acc_var_range <- range(data$acc_var)
    acc_var_sd_range <- c(0.1, 2) * sd(data$acc_var)
    
    acc_med_range <- range(data$acc_med)
    acc_med_sd_range <- c(0.1, 2) * sd(data$acc_med)
    
    return(c(
      rnorm((N-1)*N, mean = -2, sd = 1), # TPM: random around -2 (gives ~12% switching probability)
      ordered_step_params,
      log(runif(N, min = step_sd_range[1], max = step_sd_range[2])), # Step SDs: random but reasonable
      runif(N, -pi, pi), # Angle means: uniform over circle
      log(runif(N, 0.1, 5)), # Kappa: random but reasonable
      runif(N, min = acc_range[1], max = acc_range[2]),
      log(runif(N, min = acc_sd_range[1], max = acc_sd_range[2])), 
      runif(N, min = acc_var_range[1], max = acc_var_range[2]), 
      log(runif(N, min = acc_var_sd_range[1], max = acc_var_sd_range[2])) 

    ))
  } else if (scenario %in% c('base')) {
    return(c(
      rnorm((N-1)*N, mean = -2, sd = 1), # TPM: random around -2 (gives ~12% switching probability)
      ordered_step_params,
      log(runif(N, min = step_sd_range[1], max = step_sd_range[2])), # Step SDs: random but reasonable
      runif(N, -pi, pi), # Angle means: uniform over circle
      log(runif(N, 0.1, 5)) # Kappa: random but reasonable
    ))   
  } else if (scenario %in% c('base_acc')) {
    acc_range <- range(data$acc)
    acc_sd_range <- c(0.1, 2) * sd(data$acc)
    
    # similar to step length, draw ordered acceleration means as cumulative increments to guarantee state decoding matching
    raw_means <- sort(runif(N, min = acc_range[1], max = acc_range[2]), decreasing = TRUE)
    
    return(c(
      rnorm((N-1)*N, mean = -2, sd = 1), # TPM: random around -2 (gives ~12% switching probability)
      log(c(raw_means[1] - raw_means[2], raw_means[2] - raw_means[3], raw_means[3])),
      log(runif(N, min = acc_sd_range[1], max = acc_sd_range[2]))
    ))   
  } else if (scenario %in% c('coarse_abs')) {
    acc_range <- range(data$acc)
    acc_sd_range <- c(0.1, 2) * sd(data$acc)
    
    return(c(
      rnorm((N-1)*N, mean = -2, sd = 1), # TPM: random around -2 (gives ~12% switching probability)
      ordered_step_params,
      log(runif(N, min = step_sd_range[1], max = step_sd_range[2])), # Step SDs: random but reasonable
      runif(N, -pi, pi), # Angle means: uniform over circle
      log(runif(N, 0.1, 5)), # Kappa: random but reasonable
      log(runif(N, min = acc_range[1], max = acc_range[2])),
      log(runif(N, min = acc_sd_range[1], max = acc_sd_range[2])) 
    ))
  } else if (scenario %in% c('density_summary')) {
    K <- k
    return(c(
      rnorm((N-1)*N, mean = -2, sd = 1),
      ordered_step_params,
      log(runif(N, min = step_sd_range[1], max = step_sd_range[2])),
      runif(N, -pi, pi),
      log(runif(N, 0.1, 5)),
      rnorm((K-1)*N, mean = 0, sd = 1),  # random ALR-scale composition means
      log(runif(N, min = 1, max = 50))   # random concentration
    ))
  } else if (scenario %in% c('quantile_summary')) {
    K <- k
    
    q_range <- range(x[, 3:(2+K)], na.rm = TRUE)
    q_range[1] <- max(q_range[1], 1e-6)
    q_sd_range <- c(0.1, 2)*sd(x[,3], na.rm = TRUE)
    
    return(c(
      rnorm((N-1)*N, mean = -2, sd = 1),
      ordered_step_params,
      log(runif(N, min = step_sd_range[1], max = step_sd_range[2])),
      runif(N, -pi, pi),
      log(runif(N, 0.1, 5)),
      log(runif(K*N, min = q_range[1], max = q_range[2])),  
      log(runif(K*N, min = q_sd_range[1], max = q_sd_range[2]))   
    ))
  }
}
# note runif pulls N samples uniformly from the range min, max

################################################## Optimization Function (multiple initializations) ################################################## 

run_multiple_starts <- function(x, N, n_starts = 20, scenario, use_clustering = TRUE, k = NULL) {
  # store all results
  all_results <- vector("list", n_starts)
  best_mod <- NULL
  best_ll <- Inf
  
  for (i in 1:n_starts) {
    # generate starting values
    if (i == 1 && use_clustering) {   # first iteration: use clustering
      theta_start <- generate_cluster_start(x, N, scenario, iter = i, k = k)
    } else {                          # subsequent iterations: random perturbations
      theta_start <- generate_random_start(x, N, scenario, k = k) 
    }
    
    # run optimization
    mod <- tryCatch({
      nlm(mllk, theta_start, x = x, N = N, scenario = scenario, k = k, iterlim = 10000, stepmax = 10, print.level = 0) #### try a different optimizer -- optim 
    }, error = function(e) {
      print("error in one optimization round")
      list(minimum = Inf, estimate = theta_start, code = 999, iterations = 0, gradient = rep(NA, length(theta_start)))
    })
    
    # store results
    all_results[[i]] <- list(
      log_likelihood = -mod$minimum,
      estimate = mod$estimate,
      code = mod$code,
      iterations = mod$iterations,
      gradient_norm = sqrt(sum(mod$gradient^2)),
      starting_values = theta_start
    )
    
    # update best model
    converged <- mod$code %in% c(1, 2, 5)  # 1 = converged, 2 = iteration limit, 5 = step size limit
    if (converged && mod$minimum < best_ll) {
      best_ll <- mod$minimum
      best_mod <- mod
    }
  }
  
  cat("Best log-likelihood:", -best_ll, "\n")
  
  return(list(
    best = best_mod,
    all_results = all_results
  ))
}

################################################## State Decoding Function ################################################## 
#### State decoding using the viterbi algorithm
viterbi <- function(theta.star, x, N, scenario, k = NULL) {
  p <- get_allprobs(theta.star, x, N, scenario, k = k)
  n_obs <- nrow(x)
  
  state_prob_i <- matrix(-Inf, n_obs, N)  # stores max log-probabilities of being in state i given given all x up to time t (information about past x_t are accumulated)
  ps_i <- matrix(0, n_obs, N)  # stores most likely previous state for each possible current state (needed for backtracking)
  
  # initialization 
  state_prob_i[1, ] <- log(p$delta) + log(p$allprobs[1, ]) # non-matrix, non-log form: stationary probability of S1 = i * probability of x1 given S1 = i
  # solve for most likely previous state given each possible current state
  for (t in 2:n_obs) {
    for (j in 1:N) {
      probs <- state_prob_i[t-1, ] + log(p$Gamma[, j]) # prob of previous state * prob of transition to current state (j)
      
      if (length(probs) == 0 || all(is.na(probs))) {
        cat("Failure at t =", t, "j =", j, "\n")
        print(probs)
        browser()  # or stop() to inspect
      }
      
      ps_i[t, j] <- which.max(probs) # stores the most likely previous state
      state_prob_i[t, j] <- max(probs) + log(p$allprobs[t, j]) # max probability of being in state j at t
    }
  }
  
  # backtrack/decode to find most likely state sequence
  # initialization
  states <- integer(n_obs) # creates a list with the same length as the number of observations
  states[n_obs] <- which.max(state_prob_i[n_obs, ]) # most likely state at the final timestep 
  # work recursively
  for (t in (n_obs - 1):1) {
    states[t] <- ps_i[t + 1, states[t + 1]] # at each step, look up which previous state was most likely to lead to the current state
  }
  return(states)
}

################################################## Accuracy Calculation ################################################## 
accuracy <- function(decoded_states, true_states, N) {
  # create confusion matrix
  cm <- table(decoded = factor(decoded_states, levels = 1:N), true = factor(true_states,    levels = 1:N))
  # vector of all possible permutations of states
  perms <- permn(1:N)
  # find the permutation with the highest accuracy (aka highest sum of diagonal elements)
  acc <- sapply(perms, function(p) {sum(diag(cm[p, ]))})
  best <- which.max(acc)
  # create the best confusion matrix and calculate best accuracy
  best_cm <- cm[perms[[best]], ]
  accuracy <- max(acc) / sum(cm)
  return(list(
    best_cm = best_cm,
    best_acc = accuracy
  ))
}


#######################################################################################################################################
################################################## parallelization of multiple simulations ################################################## 
#######################################################################################################################################

run_one_iter <- function(i, df) {
  
  ## parallelization with help of Claude
  # one worker = one core: stop BLAS/OpenMP from spawning extra threads
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  set.seed(i)  
  
  param_i <- list()
  
  cat("Starting iteration:", i) 
  ################################################## Simulate Data ################################################## 
  # pick step lengths and turning angles
  df$step <- rgamma(n=nrow(df), shape = df$stepmean^2 / df$stepstd ^2, scale = df$stepstd ^2 / df$stepmean)
  df$angle <- mapply(function(mu, kappa) {rvonmises(1, mu = circular(mu, units = "radians"), kappa = kappa)}, df$anglemean, df$anglek)
  df$acc <- rgamma(n=nrow(df), shape = df$dbamean^2 / df$dbastd ^2, scale = df$dbastd ^2 / df$dbamean)
  
  # get coordinates
  n <- nrow(df)
  x <- numeric(n)
  y <- numeric(n)
  heading <- numeric(n)
  # initial conditions
  x[1] <- 0
  y[1] <- 0
  heading[1] <- 0
  # loop through to create coordinates
  for (t in 2:n) {
    heading[t] <- (heading[t - 1] + df$angle[t])%% (2*pi)
    x[t] <- x[t - 1] + df$step[t] * cos(heading[t])
    y[t] <- y[t - 1] + df$step[t] * sin(heading[t])
  }
  # add to df
  df$x <- x
  df$y <- y
  
  names(df)[names(df) == 's'] <- 'states'
  full_data <- df
  full_data$time <- as.integer(rownames(full_data)) # in seconds
  
  # sensor data -- remove gps except for at markers
  censored_data <- full_data
  censored_gps_data <- censored_data[which(censored_data$time %% fine_scale_points == 0), c('time', 'x', 'y')] 
  
  # compute step length and turning angle 
  censored_gps_data <- prepData(censored_gps_data, type = 'UTM', coordNames = c("x", "y"))
  
  # align coordinates (prepData shows the future step length instead of the past)
  censored_gps_data$step  <- c(NA, head(censored_gps_data$step,  -1)) 
  censored_gps_data$angle <- c(NA, head(censored_gps_data$angle, -1))
  
  # rejoin data
  censored_data <- merge(censored_gps_data, censored_data[,c('states', 'time', 'acc')], by = c('time'), all.y = TRUE) # right outer join
  
  # GPS data only at the coarse-scale
  baseline_data <- censored_data
  baseline_data <- baseline_data[which(baseline_data$time %% fine_scale_points == 0), c('time', 'x', 'y', 'step', 'angle', 'states')]
  
  empirical <- baseline_data %>%
    group_by(states) %>%
    summarise(
      mean_step = mean(step, na.rm = TRUE),
      mean_angle = mean(angle, na.rm = TRUE)
    )
  
  ################################################## Ideal Scenario ##################################################
  ##### Solve HMM with fine-scale GPS and acceleration data
  cat("\nIdeal Scenario for iteration ", i)
  # run optimization multiple times
  results <- run_multiple_starts(
    x = cbind(full_data$step, full_data$angle, full_data$acc), # maybe remove acceleration since it's just a derivation of the step and angle?
    N = 3,
    n_starts = 20, ### currently just set to 1 since it takes so long for the ideal scenario to run
    scenario = 'ideal',
    use_clustering = TRUE
  )

  # Use the best result
  mod <- results$best

  # Building the transition probability matrix from the estimated parameters
  thetastar <- mod$estimate[1:((N-1) * N)] # limit to the gamma estimates
  Gamma <- diag(N) # creates a diagonal matrix of dim N with 1s on the diagonals
  Gamma[!Gamma] <- exp(thetastar) # fills the off-diagonal elements with exponentiated thetastar values; taking the exponential allows us to obtain positive entries
  Gamma <- Gamma / rowSums(Gamma) # for each row divide the entries by the row sum to obtain row sums equal to one

  # decoding
  decoded_states <- viterbi(mod$estimate, x = cbind(full_data$step, full_data$angle, full_data$acc), N = N, scenario = 'ideal')

  a <- accuracy(decoded_states = decoded_states, true_states = full_data$states, N = N)
  acc_ideal_cm <- a$best_cm
  acc_ideal_cm
  acc_ideal <- a$best_acc
  acc_ideal

  print(acc_ideal_cm)

  param_i[["ideal"]] <- unpack_pars(mod$estimate, N, 'ideal')

  ################################################## Baseline Scenario ################################################## 
  ##### Solve HMM with coarse-scale GPS data
  cat("\nBaseline Scenarios GPS for iteration ", i)
  # run optimization multiple times
  results <- run_multiple_starts(
    x = cbind(baseline_data$step, baseline_data$angle),
    N = 3,
    scenario = 'base',
    n_starts = 20,
    use_clustering = TRUE
  )
  
  # use the best result
  mod <- results$best
  
  # decoding
  decoded_states <- viterbi(mod$estimate, x = cbind(baseline_data$step, baseline_data$angle), N = N, scenario = 'base')
  
  a <- accuracy(decoded_states = decoded_states, true_states =  baseline_data$states, N = N)
  acc_base_gps_cm <- a$best_cm
  acc_base_gps <- a$best_acc

  param_i[["base_gps"]] <- unpack_pars(mod$estimate, N, 'base')
  
  cat("\nBaseline Scenarios Acceleration for iteration ", i)
  #### baseline acceleration only
  # run optimization multiple times
  results <- run_multiple_starts(
    x = cbind(full_data$acc),
    N = 3,
    scenario = 'base_acc',
    n_starts = 20,
    use_clustering = FALSE
  )

  # use the best result
  mod <- results$best

  # decoding
  decoded_states <- viterbi(mod$estimate, x = cbind(full_data$acc), N = N, scenario = 'base_acc')

  a <- accuracy(decoded_states = decoded_states, true_states =  full_data$states, N = N)
  acc_base_acc_cm <- a$best_cm
  print(acc_base_acc_cm)
  acc_base_acc <- a$best_acc
  print(acc_base_acc)

  param_i[["base_acc"]] <- unpack_pars(mod$estimate, N, 'base_acc')

  ################################################## Coarse Scenario ##################################################
  ##### Solve HMM with coarse-scale GPS and acceleration data
  cat("\nCoarse Scenario for iteration ", i)
  # sensor data -- only keep data at 10 minute markers
  coarse_data <- full_data
  coarse_data <- coarse_data[which(coarse_data$time %% fine_scale_points == 0), c('time', 'x', 'y', 'states', 'acc')] # 10 minutes is 600 seconds
  
  # compute step length and turning angle 
  coarse_data <- prepData(coarse_data, type = 'UTM', coordNames = c("x", "y"))
  #plot(coarse_data$x, coarse_data$y, type = "l")
  # align coordinates (prepData shows the future step length instead of the past)
  coarse_data$step  <- c(NA, head(coarse_data$step,  -1)) 
  coarse_data$angle <- c(NA, head(coarse_data$angle, -1))
  
  # run optimization multiple times
  results <- run_multiple_starts(
    x = cbind(coarse_data$step, coarse_data$angle, coarse_data$acc),
    N = 3,
    scenario = 'coarse',
    n_starts = 20,
    use_clustering = TRUE
  )
  
  # use the best result
  mod <- results$best
  
  # decoding
  decoded_states <- viterbi(mod$estimate, x = cbind(coarse_data$step, coarse_data$angle, coarse_data$acc), N = N, scenario = 'coarse')
  
  a <- accuracy(decoded_states = decoded_states, true_states =  coarse_data$states, N = N)
  acc_coarse_cm <- a$best_cm
  acc_coarse <- a$best_acc

  param_i[["coarse"]] <- unpack_pars(mod$estimate, N, 'coarse')
  
  ################################################## Summary Scenario: Mean ##################################################
  ##### Solve HMM with coarse-scale GPS and a summary statistic of acceleration data
  cat("\nSummary Mean Scenario for iteration ", i)
  ##### mean as summary statistic for acceleration
  summary_data <- coarse_data[,c('step', 'angle', 'time') ]
  colnames(summary_data) <- c('coarse_step', 'coarse_angle', 'time')
  full_data1 <- merge(full_data, summary_data, by = c('time'), all.x = TRUE) 
  
  full_data1$block <- (full_data1$time - 1) %/% fine_scale_points # whole number group for every 10 minutes
  agg_mean <- data.frame(aggregate(acc ~ block, data = full_data1, FUN = mean))
  colnames(agg_mean) <- c("block", "acc_mean")
  agg_mean <- merge(full_data1, agg_mean, by = c('block'), all.x = TRUE) 
  agg_mean <- agg_mean[which(agg_mean$time %% fine_scale_points == 0), c('time', 'x', 'y', 'coarse_step', 'coarse_angle', 'states', 'acc_mean')] # 10 minutes is 600 seconds
  
  # run optimization multiple times
  results <- run_multiple_starts(
    x = cbind(agg_mean$coarse_step, agg_mean$coarse_angle, agg_mean$acc_mean),
    N = 3,
    scenario = 'mean_summary',
    n_starts = 20,
    use_clustering = TRUE
  )
  
  # use the best result
  mod <- results$best
  
  # decoding
  decoded_states <- viterbi(mod$estimate, x = cbind(agg_mean$coarse_step, agg_mean$coarse_angle, agg_mean$acc_mean), N = N, scenario = 'mean_summary')
  
  a <- accuracy(decoded_states = decoded_states, true_states =  agg_mean$states, N = N)
  acc_summary_mean_cm <- a$best_cm
  acc_summary_mean <- a$best_acc

  param_i[["mean_summary"]] <- unpack_pars(mod$estimate, N, 'mean_summary')
  
  ################################################## Summary Scenario: Variance ##################################################
  cat("\nSummary Variance Scenario for iteration ", i)
  agg_var <- data.frame(aggregate(acc ~ block, data = full_data1, FUN = var))
  colnames(agg_var) <- c("block", "acc_var")
  agg_var <- merge(full_data1, agg_var, by = c('block'), all.x = TRUE) 
  agg_var <- agg_var[which(agg_var$time %% fine_scale_points == 0), c('time', 'x', 'y', 'coarse_step', 'coarse_angle', 'states', 'acc_var')] # 10 minutes is 600 seconds
  n_agg<- fine_scale_points
  
  # run optimization multiple times
  results <- run_multiple_starts(
    x = cbind(agg_var$coarse_step, agg_var$coarse_angle, agg_var$acc_var),
    N = 3,
    scenario = 'var_summary',
    n_starts = 20,
    use_clustering = TRUE
  )
  
  # use the best result
  mod <- results$best
  
  # decoding
  decoded_states <- viterbi(mod$estimate, x = cbind(agg_var$coarse_step, agg_var$coarse_angle, agg_var$acc_var), N = N, scenario = 'var_summary')
  
  a <- accuracy(decoded_states = decoded_states, true_states =  agg_var$states, N = N)
  acc_summary_var_cm <- a$best_cm
  print(acc_summary_var_cm)
  acc_summary_var <- a$best_acc
  print(acc_summary_var)
  
  param_i[["var_summary"]] <- unpack_pars(mod$estimate, N, 'var_summary')
  
  ################################################## Summary Scenario: Median ##################################################
  cat("\nSummary Median Scenario for iteration ", i)
  agg_med <- data.frame(aggregate(acc ~ block, data = full_data1, FUN = median))
  colnames(agg_med) <- c("block", "acc_med")
  agg_med <- merge(full_data1, agg_med, by = c('block'), all.x = TRUE) 
  agg_med <- agg_med[which(agg_med$time %% fine_scale_points == 0), c('time', 'x', 'y', 'coarse_step', 'coarse_angle', 'states', 'acc_med')] # 10 minutes is 600 seconds
  n_agg<- fine_scale_points
  
  # run optimization multiple times
  results <- run_multiple_starts(
    x = cbind(agg_med$coarse_step, agg_med$coarse_angle, agg_med$acc_med),
    N = 3,
    scenario = 'med_summary',
    n_starts = 20,
    use_clustering = TRUE
  )
  
  # use the best result
  mod <- results$best
  
  # decoding
  decoded_states <- viterbi(mod$estimate, x = cbind(agg_med$coarse_step, agg_med$coarse_angle, agg_med$acc_med), N = N, scenario = 'med_summary')
  
  a <- accuracy(decoded_states = decoded_states, true_states =  agg_med$states, N = N)
  acc_summary_med_cm <- a$best_cm
  print(acc_summary_med_cm)
  acc_summary_med <- a$best_acc
  print(acc_summary_med)
  
  param_i[["med_summary"]] <- unpack_pars(mod$estimate, N, 'med_summary')
  
  ################################################## Summary Scenario: Summary Statistics ##################################################
  cat("\nSummars Stats Scenario for iteration ", i)
  agg_stats_temp <- merge(agg_mean, agg_var, by = c('time', 'x', 'y', 'coarse_step', 'coarse_angle', 'states'), all.x = TRUE) 
  agg_stats <- merge(agg_stats_temp, agg_med, by = c('time', 'x', 'y', 'coarse_step', 'coarse_angle', 'states'), all.x = TRUE) 
  
  agg_stats$acc_sd <- sqrt(agg_stats$acc_var)
  agg_stats <- agg_stats[which(agg_stats$time %% fine_scale_points == 0), c('time', 'x', 'y', 'coarse_step', 'coarse_angle', 'states', 'acc_mean', 'acc_var', 'acc_med')] # 10 minutes is 600 seconds
  n_agg<- fine_scale_points
  
  # try scaling, so neither dominates during likelihood
  mean_scale <- sd(agg_stats$acc_mean, na.rm = TRUE)
  var_scale <- sd(agg_stats$acc_var,  na.rm = TRUE)
  med_scale <- sd(agg_stats$acc_med,  na.rm = TRUE)
  
  agg_stats$acc_mean_input <- agg_stats$acc_mean / mean_scale
  agg_stats$acc_var_input <- agg_stats$acc_var / var_scale
  agg_stats$acc_med_input <- agg_stats$acc_med / med_scale
  
  # run optimization multiple times
  results <- run_multiple_starts(
    x = cbind(agg_stats$coarse_step, agg_stats$coarse_angle, agg_stats$acc_mean_input, agg_stats$acc_var_input, agg_stats$acc_med_input),
    N = 3,
    scenario = 'stats_summary',
    n_starts = 20,
    use_clustering = TRUE
  )
  
  # use the best result
  mod <- results$best
  
  # decoding
  decoded_states <- viterbi(mod$estimate, x = cbind(agg_stats$coarse_step, agg_stats$coarse_angle, agg_stats$acc_mean_input, agg_stats$acc_var_input, agg_stats$acc_med_input), N = N, scenario = 'stats_summary')

  a <- accuracy(decoded_states = decoded_states, true_states =  agg_stats$states, N = N)
  acc_summary_stats_cm <- a$best_cm
  acc_summary_stats <- a$best_acc

  param_i[["stats_summary"]] <- unpack_pars(mod$estimate, N, 'stats_summary')
  
  
  ################################################## Density Scenario ##################################################
  cat("\nDensity Scenario for iteration ", i)
  den_data <- censored_data
  den_data <- den_data[!is.na(den_data$acc), ]  # remove NAs from acc
  den_data$block <- (den_data$time - 1) %/% fine_scale_points # whole number group for every 600 datapoints, start at 1 not 0
  
  acc_min <- min(den_data$acc)
  acc_max <- max(den_data$acc)
  
  # define fixed evaluation points where we will evaluate the density estimate
  k <- 5 # number of grid points
  grid_pts <- quantile(den_data$acc, probs = seq(0.2, 0.8, length.out = k))
  
  # fit KDE per block and evaluate at fixed grid points
  density_df <- den_data %>%
    group_by(block) %>%
    summarise(
      density_vector = list({
        # Fit KDE on raw observations for this block
        # bw = "SJ" is more robust than the default "nrd0" for non-normal data
        kde <- density(acc,
                       bw     = "SJ",
                       kernel = "gaussian",
                       from   = acc_min,
                       to     = acc_max,
                       n      = 512)      # internal grid resolution, not output length
        
        # Evaluate at your fixed grid points by interpolation
        vals <- approx(kde$x, kde$y, xout = grid_pts)$y
        # prevent 0s since Dirichlet density cant handle 0
        vals <- vals + 1e-6
        # normalize
        vals / sum(vals)
      }),
      n = n(),
      .groups = "drop"
    )
  
  density_df
  
  # unnest density_vector into k separate columns:
  density_wide <- density_df %>%
    mutate(row_id = block) %>%
    tidyr::unnest_wider(density_vector, names_sep = "_") %>%
    rename_with(~ gsub("density_vector", "dens", .x), starts_with("density_vector"))
  
  density_wide$time <- density_wide$block * fine_scale_points + fine_scale_points
  
  agg_density <- merge(censored_data, density_wide, by = c("time"), all.x = TRUE)
  agg_density <- agg_density[which(agg_density$time %% fine_scale_points == 0), c('time', 'x', 'y', 'step', 'angle', 'states', paste0("dens_", 1:k))]
  
  results <- run_multiple_starts(
    x = cbind(agg_density$step, agg_density$angle, as.matrix(agg_density[, paste0("dens_", 1:k)])),
    N = 3,
    scenario = 'density_summary',
    n_starts = 20,
    use_clustering = TRUE, 
    k = k     
  )
  
  mod <- results$best
  
  decoded_states <- viterbi(mod$estimate, x = cbind(agg_density$step, agg_density$angle, as.matrix(agg_density[, paste0("dens_", 1:k)])), N = N, scenario = 'density_summary', k=k)

  a <- accuracy(decoded_states = decoded_states, true_states =  agg_density$states, N = N)
  acc_summary_density_cm <- a$best_cm
  acc_summary_density <- a$best_acc

  param_i[["density"]] <- unpack_pars(mod$estimate, N, 'density_summary', k = k)
  
  ################################################## Quantile Scenario ##################################################
  cat("\nQuantile Scenario for iteration ", i)
  k <- 7
  
  # find quantiles for each block
  quantile_df <- den_data %>%
    group_by(block) %>%
    summarise(
      density_vector = list({
        # Evaluate at quantiles
        vals <- quantile(acc, probs = seq(0.2, 0.8, length.out = k))
        # prevent 0s since Dirichlet density cant handle 0
        pmax(vals, 1e-6)
      }),
      n = n(),
      .groups = "drop"
    )
  
  # unnest density_vector into k separate columns:
  quantile_wide <- quantile_df %>%
    mutate(row_id = block) %>%
    tidyr::unnest_wider(density_vector, names_sep = "_") %>%
    rename_with(~ paste0("dens_", seq_along(.x)), starts_with("density_vector"))
  
  quantile_wide$time <- quantile_wide$block * fine_scale_points + fine_scale_points
  
  agg_quantile <- merge(censored_data, quantile_wide, by = c("time"), all.x = TRUE)
  agg_quantile <- agg_quantile[which(agg_quantile$time %% fine_scale_points == 0), c('time', 'x', 'y', 'step', 'angle', 'states', paste0("dens_", 1:k))]
  
  results <- run_multiple_starts(
    x = cbind(agg_quantile$step, agg_quantile$angle, as.matrix(agg_quantile[, paste0("dens_", 1:k)])),
    N = 3,
    scenario = 'quantile_summary',
    n_starts = 20,
    use_clustering = TRUE,
    k = k
  )
  
  mod <- results$best
  
  decoded_states <- viterbi(mod$estimate, x = cbind(agg_quantile$step, agg_quantile$angle, as.matrix(agg_quantile[, paste0("dens_", 1:k)])), N = N, scenario = 'quantile_summary', k = k)

  a <- accuracy(decoded_states = decoded_states, true_states =  agg_quantile$states, N = N)
  acc_summary_quantile_cm <- a$best_cm
  acc_summary_quantile <- a$best_acc

  param_i[["quantile"]] <- unpack_pars(mod$estimate, N, 'quantile_summary', k = k)
  
  
  ################################################## Accuracy Summary ################################################## 
  # Store results
  list(
    acc = c(
      #acc_ideal = acc_ideal,
      acc_base_gps = acc_base_gps,
      #acc_base_acc = acc_base_acc,
      acc_coarse = acc_coarse,
      acc_summary_mean = acc_summary_mean,
      acc_summary_var = acc_summary_var,
      acc_summary_med = acc_summary_med,
      acc_summary_stats = acc_summary_stats,
      acc_summary_density = acc_summary_density,
      acc_summary_quantile = acc_summary_quantile
    ),
    cm = list(
      #ideal = acc_ideal_cm,
      base_gps = acc_base_gps_cm,
      #base_acc = acc_base_acc_cm,
      coarse = acc_coarse_cm,
      summary_mean = acc_summary_mean_cm,
      summary_var = acc_summary_var_cm,
      summary_med = acc_summary_med_cm,
      summary_stats = acc_summary_stats_cm,
      summary_density = acc_summary_density_cm,
      summary_quantile = acc_summary_quantile_cm
    ), 
    params = param_i,
    empiricals = empirical
  )
  
}


##### for testing one iteration
# n_sims <- 1
# plan(sequential)
# one <- run_one_iter(1, df)
# names(one)                     # must be: acc, cm, params, empiricals
# str(one$acc)                   # named numeric length 8
# str(one$params$stats_summary)  # non-NULL Gamma, sd fields present
# str(one$empiricals)            # tibble with mean_step per state


### parallelization done with help from Claude 
n_sims <- 20

plan(multisession, workers = 4)

results_list <- future_map(
  1:n_sims, 
  run_one_iter, 
  df= df,
  .options = furrr_options(
    seed = TRUE,                                   # correct, reproducible parallel RNG
    packages = c("moveHMM", "tidyr", "dplyr", "circular",
                 "cluster", "gtools", "combinat", "RhpcBLASctl"),
    globals = TRUE                                 # ship df, N, fine_scale_points, and your functions to workers
  )
)

### reshape results
output <- as.data.frame(do.call(rbind, lapply(results_list, `[[`, "acc")))
cm_results <- lapply(results_list, `[[`, "cm")
param_results <- lapply(results_list, `[[`, "params")
empirical_results <- lapply(results_list, `[[`, "empiricals") ## note this isn't output in the files below


# acc <- colMeans(output)
# write.csv(output, "output_normal.csv")
# saveRDS(cm_results, "cm_results_normal.rds")
# saveRDS(param_results, "param_results_longer_normal.rds")

plan(sequential)
