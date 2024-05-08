# Tempering algorithm.

if (TRUE) {
  .libPaths("~/Rlib")
  setwd("~/Code/ising_model/scripts")
  source("algorithms.r")
  source("tools.r")
}

tempering <- function(A, beta_vector, init_matrix,
                      n_exchange, n_mc_iter, B = 0,
                      f1, f2 = f1, 
                      which_sampler = rep(1, length(beta_vector)),
                      n_states = 2,
                      is_potts = FALSE) {
  # Returns samples for graph A, across all temperatures
  # specified by the beta_vector. Uses tempering (or exchange)
  # Monte Carlo to improve sampling at high betas, i.e. cold temperatures.
  #
  # A: graph
  # beta_vector: betas for each replica
  # init_matrix: inits for each replica
  # n_exchange: number of exchanges.
  # n_mc_iter: number of sampling iterations between exchanges.
  # B: exterior magnetic field
  # f1: function for first sampler. Should accept A, beta, B,
  #     an init vector, n_mc_iter, and n_states as arguments.
  # f2: function for second sampler. Same as above.
  # which_sampler: which sampler to use for each replica. This
  #                should be a vector with the same length as
  #                beta_vector with entries 1 and 2.
  #
  n_replica <- length(beta_vector)
  n_iter <- n_exchange * n_mc_iter
  n_particles <- nrow(A)

  X <- array(NA, dim = c(n_replica, n_iter, n_particles))

  current_init <- init_matrix

  exchange_rate <- rep(0, n_replica - 1)
  
  for (i in 1:n_exchange) {
    first_index <- (i - 1) * n_mc_iter + 1
    last_index <- i * n_mc_iter
    for (b in 1:n_replica) {
      if (which_sampler[b] == 1) {
        X[b, first_index:last_index, ] <-
          f1(A, beta_vector[b], B, current_init[b, ],
             n_mc_iter, n_states)
      } else {
        X[b, first_index:last_index, ] <-
          f2(A, beta_vector[b], B, current_init[b, ],
             n_mc_iter, n_states)
      }
    }

    # exchange step
    # NOTE: pass 1 instead of beta to ising kernel to not double-count beta.
    for (b in 1:(n_replica - 1)) {
      if (b == 1) {
        if (!is_potts) {
          E_low = - ising_kernel(1, A, X[b, last_index, ], B, T)
        } else {
          E_low = - 0.5 * potts_log_kernel(1, A, X[b, last_index, ], 
                                           n_states = n_states)
        }
      } else { 
        E_low = E_high
      }

      if (!is_potts) {
        E_high = - ising_kernel(1, A, X[b + 1, last_index, ], B, T)
      } else {
        E_high = - 0.5 * potts_log_kernel(1, A, X[b + 1, last_index, ],
                                          n_states = n_states)
      }

      delta = (beta_vector[b + 1] - beta_vector[b]) * (E_high - E_low)

      if (runif(1) < exp(delta)) {
        saved_state = X[b + 1, last_index, ]
        X[b + 1, last_index, ] = X[b, last_index, ]
        X[b, last_index, ] = saved_state

        E_high = E_low  # update energy before next exchange

        exchange_rate[b] <- exchange_rate[b] + 1
      }

      current_init[b, ] = X[b, last_index, ]    }

    current_init[n_replica, ] = X[n_replica, last_index, ]
  }
  print(exchange_rate / n_exchange)

  X
}

###############################################################################
## Wrappers for samplers to be passed to tempering()

f_AG <- function(A, beta, B, init, n_iter, n_states) {
  ag_potts_simple(A = A, beta = beta, init = init, n_states = n_states,
                  n_iter = n_iter)
}

f_AG_lowrank <- function(A, beta, B, init, n_iter, n_states) {
  ag_potts_lowrank(A = A, beta = beta, init = init, n_states = n_states,
                   n_iter = n_iter, epsilon = 1e-15)
}

f_HB <- function(A, beta, B, init, n_iter, n_states) {
  hb_sampler(A = A, beta = beta, init = init,
             n_iter = n_iter, is_potts = TRUE, 
             n_states = n_states,
             sub_sample = 1)  
}

sub_adjust <- 10
f_HB_long <- function(A, beta, B, init, n_iter, n_states) {
  hb_sampler(A = A, beta = beta, init = init,
             n_iter = n_iter, is_potts = TRUE, 
             n_states = n_states,
             sub_sample = ceiling(graph_size / sub_adjust))
}

f_BHB <- function(A, beta, B, init, n_iter, n_states) {
  hb_sampler(A = A, beta = beta, init = init,
             n_iter = n_iter, is_potts = TRUE, 
             n_states = n_states,
             sub_sample = 1, black_box = TRUE)  
}

f_GWG <- function(A, beta, B, init, n_iter, n_states) {
  gwg_potts(A = A, beta = beta, init = init,
            n_states = n_states, n_iter = n_iter)
}

###############################################################################
## Unit test
unit_test <- FALSE
# WARNING: unit tests have not been updated with latest signature
# from tempering (with which_sampler argument)
if (unit_test) {
  n_part <- 64
  A <- adjacency_graph(n_part, type = "gaussian")
  beta_vector <- c(0.5, 1.8, 2.3, 2.6, 3.0, 3.25, 3.5, 3.70, 3.90, 4.0) 
    # c(seq(from = 0.5, to = 4, by = 0.35))
    # c(1.0, 1.0, 1.5, 2.0, 2.5)  # c(0.5, 0.75, 1.0, 1.25)
  B <- 0
  
  n_replica <- length(beta_vector)
  init_matrix <- matrix(sample(c(-1, 1), n_part * n_replica, replace = T),
                        nrow = n_replica, ncol = n_part)
  
  n_mc_iter <- 5e3
  n_exchange <- 40
  
  f1 <- sumit_sampler
  f2 <- function(A, beta, B, init, n_iter) { 
    hb_sampler(A, beta, B, init, n_iter, 
               sub_sample = ceiling(nrow(A) / 10))
  }
  threshold = 4.1

  sample <- tempering(A, beta_vector, init_matrix,
                      n_exchange, n_mc_iter, B = 0,
                      f1, f2, threshold)

  temp_estimates <- monte_carlo_estimate(sample[, -(1:1000), ])

  threshold = 0
  sample_hb <- tempering(A, beta_vector, init_matrix,
                         n_exchange, n_mc_iter, B, f1, f2, threshold)
  
  temp_hb_estimates <- monte_carlo_estimate(sample_hb[, -(1:1000), ])
  
  ## Benchmark against sumit's sampler without tempering
  n_iter <- n_mc_iter * n_exchange
  sample_const <- array(NA, c(n_replica, n_iter, n_part))
  for (i in 1:n_replica) {
    sample_const[i, , ] <- sumit_sampler(A, beta_vector[i], n_iter = n_iter,
                                       init = init_matrix[i, ])
  }

  const_estimates <- 
    monte_carlo_estimate(sample_const[, -(1:1000), ])

  temp_estimates
  temp_hb_estimates
  
  const_estimates
}
