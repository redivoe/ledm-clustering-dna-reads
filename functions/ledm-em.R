
deterministic_X <- function(Z, E, U){
  m <- length(Z)
  n <- nrow(E)
  
  X <- vector(mode = "list", length = n)
  l <- m + rowSums(E == 4) - rowSums(E == 2)
  for(i in 1:n){
    X[[i]] <- integer(length = l[i])
    h <- 1
    for(j in 1:m){
      if(E[i, j] == 1){
        X[[i]][h] <- Z[j]
        h <- h + 1
      }
      if(E[i, j] == 2){
        h <- h
      }
      if(E[i, j] == 3){
        X[[i]][h] <- U[i, j]
        h <- h + 1
      }
      if(E[i, j] == 4){
        X[[i]][h] <- Z[j]
        h <- h + 1
        X[[i]][h] <- U[i, j]
        h <- h + 1
      }
    }
  }
  return(X)
}

sample_led <- function(n,
                       p,
                       eta = rep(1/4, 4),
                       theta = c(0.7, 0.1, 0.1, 0.1),
                       epsilon = 1){
  
  Z <- sample.int(n = 4, size = p,  prob = eta, replace = TRUE)
  E <- array(sample.int(n = 4, size = n * p, replace = TRUE, prob = theta), dim = c(n, p))
  substitute_sets <- lapply(1:4, \(x) setdiff(1:4, x))
  
  U <- array(dim = c(n, p))
  for(i in 1:n){
    is_insert <- E[i, ] == 4
    U[i, is_insert] <- sample.int(n = 4, size = sum(is_insert), replace = TRUE, prob = eta)
    which_substitute <- which(E[i, ] == 3)
    if(length(which_substitute) > 0){
      for(j in which_substitute){
        U[i, j] <- sample(x = c(Z[j], substitute_sets[[Z[j]]]), 1, prob = c(1 - epsilon, rep(epsilon/3, 3)))
      }
    }
  }
  X <- deterministic_X(Z, E, U)
  
  return(list("X" = X,
              "Z" = Z,
              "E" = E,
              "U" = U))
}


led_em <- function(X,
                   p = 50,
                   B = 20,
                   max_em_iter = 10,
                   Z_init = NULL,
                   seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- length(X)
  theta <- c("K" = 0.85, "D" = 0.05, "S" = 0.05, "I" = 0.05)
  
  if (is.null(Z_init)) {
    Z <- sample.int(n = 4, size = p, replace = TRUE, prob = rep(1/4, 4))
  } else {
    Z <- Z_init
  }
  
  Z_move_trace <- integer(max_em_iter)
  loglik <- numeric(max_em_iter)
  
  for (iter in seq_len(max_em_iter)) {
    log_theta <- log(theta)
    log_Nc <- rep(-Inf, 4)
    
    for (i in seq_len(n)) {
      l <- length(X[[i]])
      
      log_prob_emissions <- get_log_prob_emissions_cpp(X[[i]], Z, B)
      log_alpha <- get_log_alpha_cpp(log_theta, log_prob_emissions, B)
      log_beta <- get_log_beta_cpp(log_theta, log_prob_emissions, B)
      log_gamma_num <- get_log_gamma_num_cpp(log_theta, log_alpha, log_beta, log_prob_emissions, B)
      
      loglik_i <- log_alpha[p + 1, l + 1]
      loglik[iter] <- loglik[iter] + loglik_i
      
      log_gamma <- log_gamma_num - loglik_i
      
      log_Nc[1] <- matrixStats::logSumExp(c(log_Nc[1], log_gamma[, , 1]))
      log_Nc[2] <- matrixStats::logSumExp(c(log_Nc[2], log_gamma[, , 2]))
      log_Nc[3] <- matrixStats::logSumExp(c(log_Nc[3], log_gamma[, , 3]))
      log_Nc[4] <- matrixStats::logSumExp(c(log_Nc[4], log_gamma[, , 4]))
    }
    
    theta <- exp(log_Nc - log(n * p))
    
    Z_old <- Z
    Z <- update_Z_cpp_fast(X, Z, log_theta, rep(1, n), B)
    Z_move_trace[iter] <- sum(Z != Z_old)
    
    cat("Iterazione", iter, "\n",
        "- loglik = ", loglik[iter], "\n",
        "- Z changed =", Z_move_trace[iter], "\n",
        "- theta = ", theta, "\n")
  }
  
  return(list(
    "theta" = theta,
    "Z" = Z,
    "loglik" = loglik
  ))
}



init_Z_consensus_raw <- function(X, p) {
  Z_init <- integer(p)
  
  for (j in seq_len(p)) {
    vals_j <- sapply(X, function(x) { if (length(x) >= j) x[j] else NA})
    vals_j <- vals_j[!is.na(vals_j)]
    
    if (length(vals_j) == 0) {
      Z_init[j] <- sample.int(n = 4, size = 1)
    } else {
      tab <- tabulate(vals_j, nbins = 4)
      Z_init[j] <- which.max(tab)
    }
  }
  return(Z_init)
}



##-----------
##  Mixture  
##-----------

sample_ledm <- function(nk,
                        p,
                        theta = c(0.8, 0.05, 0.1, 0.05),
                        eta = rep(1/4, 4),
                        epsilon = 1){
  
  K <- length(nk)
  Xk <- vector(mode = "list", length = K)
  Z <- matrix(NA, nrow = K, ncol = p)
  for(k in 1:K){
    out_k <- sample_led(n = nk[k], p = p, theta = theta, epsilon = epsilon)
    Xk[[k]] <- out_k$X
    Z[k, ] <- out_k$Z
  }
  X <- do.call(c, Xk)
  cl <- rep(1:K, times = nk)
  return(list("X" = X, "Z" = Z, "cl" = cl))
}

init_cl_ledm <- function(X, K, stringdist_method = "lv") {
  seqs <- sapply(X, paste0, collapse = "")
  dist_matrix <- stringdist::stringdistmatrix(a = seqs, method = stringdist_method, q = 3)
  hc <- hclust(d = dist_matrix, method = "complete")
  cl <- cutree(hc, k = K)
  return(cl)
}

init_Z_consensus_mixture <- function(X, p, cl) {
  K <- length(unique(cl))
  Z_init <- matrix(NA, nrow = K, ncol = p)
  for (k in seq_len(K)) {
    Xk <- X[cl == k]
    Z_init[k, ] <- init_Z_consensus_raw(Xk, p)
  }
  return(Z_init)
}


ledm_em <- function(X,
                    p = 50,
                    K,
                    B = 20,
                    Z_init_method = "consensus",
                    cl_init_method = "stringdist",
                    stringdist_method = "lv",
                    init = NULL,
                    seed = NULL,
                    max_em_iter = 100,
                    eps = 1e-3) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- length(X)
  
  # initialization
  if(is.null(init$theta)){
    theta <- c("K" = 0.85, "D" = 0.05, "S" = 0.05, "I" = 0.05)
  }else{
    theta <- init$theta
  }
  log_theta <- log(theta)
  
  stopifnot(cl_init_method %in% c("stringdist", "random"))
  if(is.null(init$cl)){
    if(K == 1){
      cl_init <- rep(1, n)
    }else{
      if(cl_init_method == "stringdist"){
        cl_init <- init_cl_ledm(X = X, K = K, stringdist_method = stringdist_method)
      }else if(cl_init_method == "random"){
        cl_init <- sample.int(n = K, size = n, replace = TRUE)
      }
    }
  }else{
    if(length(unique(init$cl)) != K){
      stop("The initial clustering does not have K groups.")
    }
    cl_init <- init$cl
  }
  
  pi_k <- tabulate(cl_init, K) / n
  
  stopifnot(Z_init_method %in% c("consensus", "random"))
  if(is.null(init$Z)){
    if(Z_init_method == "consensus"){
      Z_init <- init_Z_consensus_mixture(X, p, cl_init)
    }else if(Z_init_method == "random"){
      Z_init <- matrix(sample.int(n = 4, size = K * p, replace = TRUE, prob = rep(1/4, 4)), nrow = K, ncol = p)
    }
  }else{
    Z_init <- init$Z
  }
  Z <- Z_init
  
  Z_move_trace <- integer(max_em_iter)
  loglik <- rep(0, max_em_iter)
  tau <- ll <- matrix(NA, nrow = n, ncol = K)
  log_gamma_edit <- vector(mode = "list", length = 4)
  names(log_gamma_edit) <- names(theta)
  for(h in 1:4){
    log_gamma_edit[[h]] <- matrix(NA, nrow = n, ncol = K)
  }
  
  B_vec <- pmax(B, abs(sapply(X, length) - p))
  
  cat("EM iteration: ")
  for (iter in seq_len(max_em_iter)) {

    for(k in 1:K){
      for (i in seq_len(n)) {
        l <- length(X[[i]])
        log_prob_emissions <- get_log_prob_emissions_cpp(X[[i]], Z[k, ], B_vec[i])
        log_alpha <- get_log_alpha_cpp(log_theta, log_prob_emissions, B_vec[i])
        log_beta <- get_log_beta_cpp(log_theta, log_prob_emissions, B_vec[i])
        log_gamma_num <- get_log_gamma_num_cpp(log_theta, log_alpha, log_beta, log_prob_emissions, B_vec[i])

        ll[i, k] <- log_alpha[p + 1, l + 1]
        log_gamma <- log_gamma_num - ll[i, k]

        if(all(log_gamma_num == -Inf) & (ll[i, k] == -Inf)){
          for(h in 1:4){
            log_gamma_edit[[h]][i, k] <- -Inf
          }
        }else{
          for(h in 1:4){
            log_gamma_edit[[h]][i, k] <- matrixStats::logSumExp(log_gamma[, , h])
          }
        }
      }
    }

    log_fx <- matrixStats::rowLogSumExps(t(t(ll) + log(pi_k)))
    log_tau <-  t(t(ll) + log(pi_k)) - log_fx
    tau <- exp(t(t(ll) + log(pi_k)) - log_fx)
    
    log_Nc <- sapply(log_gamma_edit, \(x) matrixStats::logSumExp(x + log_tau))
    theta <- exp(log_Nc - log(n * p))
    log_theta <- log(theta)
    
    pi_k <- colMeans(tau)
    
    loglik[iter] <- sum(log_fx)
    
    Z_old <- Z
    for(k in seq_len(K)){
      Z[k, ] <- update_Z_cpp_fast(X, Z[k, ], log_theta, tau[, k], B_vec)
    }
    Z_move_trace[iter] <- sum(Z != Z_old)
    
    cat(iter, "/ ")
    # cat(loglik[iter], " ")
    if (iter > 10 && abs((loglik[iter] - loglik[iter - 1]) / loglik[iter - 1]) < eps) {
      break
    }
  }
  
  cl <- max.col(tau)
  
  n_params <- (K - 1) + K * p + 3
  ell_max <- mean(loglik[(iter - 2):iter])
  bic <- -2 * ell_max + n_params * log(n)
  aic <- -2 * ell_max + n_params * 2
  
  return(
    list(
      "theta" = theta,
      "Z" = Z,
      "cl" = cl,
      "loglik" = loglik[1:iter],
      "bic" = bic,
      "aic" = aic,
      "tau" = tau,
      "pi_k" = pi_k,
      "cl_init" = cl_init,
      "Z_init" = Z_init,
      "Z_move_trace" = Z_move_trace[1:iter]
    )
  )
}

ledm_em_sel <- function(X,
                        p = 50,
                        K_seq = 1:10,
                        selection_criterion = "bic",
                        B = 20,
                        Z_init_method = "consensus",
                        cl_init_method = "stringdist",
                        stringdist_method = "lv",
                        init = NULL,
                        seed = NULL,
                        max_em_iter = 100,
                        eps = 1e-3) {
  
  out_sel <- vector(mode = "list", length = length(K_seq))
  for(i in 1:length(K_seq)){
    cat("Estimating with K = ", K_seq[i], "\n", sep = "")
    out_sel[[i]] <- ledm_em(
      X = X,
      p = p,
      K = K_seq[i],
      B = B,
      Z_init_method = Z_init_method,
      cl_init_method = cl_init_method,
      stringdist_method = stringdist_method,
      init = init,
      seed = seed,
      max_em_iter = max_em_iter,
      eps = eps
    )
    cat("\n")
  }
  
  bic_seq <- sapply(out_sel, \(x) x$bic)
  aic_seq <- sapply(out_sel, \(x) x$aic)
  
  which_selected_bic <- which.min(bic_seq)
  which_selected_aic <- which.min(aic_seq)
  
  if(selection_criterion == "bic"){
    out <- out_sel[[which_selected_bic]]
  }else if(selection_criterion == "aic"){
    out <- out_sel[[which_selected_aic]]
  }
  
  out$K_selected_bic <- K_seq[which_selected_bic]
  out$K_selected_aic <- K_seq[which_selected_aic]
  
  out$bic_seq <- bic_seq
  out$aic_seq <- aic_seq
  
  return(out)
}
