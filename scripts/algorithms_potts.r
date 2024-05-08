
# Algorithms_potts
# Extend algorithms to work on pott sampler.

# library(extraDistr)  # CHECK -- do we need this?
source("tools.r")
source("algorithms.r")

# Simulate states from the Potts model, using auxiliary Gaussian variables.
#
# A: adjacency matrix
# beta: inverse-temperatures
# init: intial state for the chain
# n_states: number of states each particle can adopt.
# n_iter: number of sampling iteration.
# B (not implemented): external magnetic field.
# alpha_modif: modification along diagonal of A to make it semi-positive def.
# retain_alpha: whether to use alpha_modif or instead use min eigenvalue.
# pertub (not implemented): if TRUE, adds flip transition.
#
# return: matrix with samples.
ag_potts_simple <- function (A, beta, init, n_states = 2, n_iter = 1000, B = 0,
                             alpha_modif = 1e-7, retain_alpha = FALSE,
                             perturb = FALSE) {
  if (B != 0) print("WARNING: non-zero magnetic field is not supported.")
  if (perturb) print("WARNING: pertub = TRUE is not supported.")
  
  if (retain_alpha) {
    alpha = alpha_modif
  } else {
    alpha = abs(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
    alpha = alpha + alpha_modif
  }

  n_particles = nrow(A)
  C = 2 * beta * (A + alpha * diag(n_particles))
  L = t(chol(solve(C)))
  X = matrix(NA, nrow = n_iter, ncol = n_particles)
  x = init
  X[1, ] = x

  # Auxiliary Gaussians
  Q = matrix(NA, nrow = n_particles, ncol = n_states)

  for (i in 2:n_iter) {
    for (l in 1:n_states) {
      mu_q = rep(0, n_particles)
      mu_q[X[i - 1, ] == l] = 1
      Q[, l] = mvrnorm_cholesky(n = 1, mu = mu_q, L = L)
    }

    P = C %*% Q

    u = runif(n_particles, 0, 1)

    for (j in 1:n_particles) {
      p = rep(NA, n_states)
      for (l in 1:n_states) {
        p[l] = exp(P[j, l])
      }
      norm_const = sum(p)

      state = 1
      p_current = p[state] / norm_const
      while (u[j] > p_current) {
        state = state + 1
        p_current = p_current + p[state] / norm_const
      }
      
      X[i, j] = state
    }
  }
  X
}


###############################################################################

# Simulate states from the Potts model, using auxiliary Gaussian variables,
# and a low rank approximation. Takes in the modified adjacency matrix,
# C_tilde. See the wrapper below, which takes in A and beta.
# 
#
# C_tilde: modified adjacency matrix, which is low rank, semi+ def,
#          and scaled by the inverse-temperature beta.
# init: intial state for the chain
# n_states: number of states each particle can adopt.
# n_iter: number of sampling iteration.
# epsilon: value below which an eigenvalue is considered to be 0.
# eigenvalue_rounding: determines how different two eigen values need to
#                      be in order to be considered non-identical.
#                      Required to handle numerical precision of LAPACK. 
#
# return: matrix with samples.
ag_potts_lowrank__ <- function (C, init, n_states = 2, n_iter = 1000,
                                epsilon = 1e-12, eigenvalue_rounding = 1e12) {
  n_particles = nrow(C)
  spectral_decomposition = eigen(C, symmetric = TRUE)

  n_rank = length(
    spectral_decomposition$values[spectral_decomposition$values > epsilon])
  eigen_values_sqrt = sqrt(spectral_decomposition$values[1:n_rank])

  scaled_eigenvectors <- matrix(NA, nrow = n_particles, ncol = n_rank)
  for (i in 1:n_rank) {
    scaled_eigenvectors[, i] = spectral_decomposition$values[i] *
      spectral_decomposition$vectors[, i]
  }
  
  X = matrix(NA, nrow = n_iter, ncol = n_particles)
  x = init
  X[1, ] = x

  # Auxiliary Gaussians
  Q = matrix(NA, nrow = n_rank, ncol = n_states)

  for (i in 2:n_iter) {
    for (l in 1:n_states) {
      Y = rep(0, n_particles)
      Y[X[i - 1, ] == l] = 1
      mu_q = t(spectral_decomposition$vector[, 1:n_rank]) %*% Y

      # set.seed(1954)
      standard_normal = rnorm(n_rank)
      Q[, l] = standard_normal / eigen_values_sqrt + mu_q
    }

    P = scaled_eigenvectors %*% Q

    u = runif(n_particles, 0, 1)

    for (j in 1:n_particles) {
      p = exp(P[j, ])
      norm_const = sum(p)

      state = 1
      p_current = p[state] / norm_const
      while (u[j] > p_current) {
        state = state + 1
        p_current = p_current + p[state] / norm_const
      }

      X[i, j] = state
    }
  }
  X
}

# Wrapper around the above function, which takes A and beta.
ag_potts_lowrank <- function(A, beta, init, n_states = 2, n_iter = 1000,
                             alpha_modif = 1e-7, epsilon = 1e-12) {
  n_particles = nrow(A)
  alpha = abs(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
  alpha = alpha + alpha_modif
  C = 2 * beta * (A + alpha * diag(n_particles))

  ag_potts_lowrank__(C, init, n_states, n_iter, epsilon)
}

# Wrapper function for all AG algorithms
ag_sampler <- function(A, beta, init, n_states = 2, n_iter = 1000,
                       alpha_modif = 1e-7, lowrank = FALSE, epsilon = 1e-12,
                       B = 0, retain_alpha = FALSE, perturb = FALSE) {
  if (n_states == 2) {
    if (lowrank == TRUE) {
      X = ag_ising_lowrank(A, beta, init, n_iter, epsilon)
    } else {
      # Magnetic field is set to 0.
      X = sumit_sampler(A, beta, B, init, n_iter, alpha_modif, retain_alpha,
                        perturb)
    }
  } else {
    if (lowrank == TRUE) {
      X = ag_potts_lowrank(A, beta, init, n_states, n_iter, alpha_modif,
                           epsilon)
    } else {
      X = ag_potts_simple(A, beta, init, n_states, n_iter, B, alpha_modif,
                          retain_alpha, perturb)
    }
  }
  return(X)
}


###############################################################################
## Gibbs with gradient (Grathwohl et al 2021)
# NOTE: hand-coded categorical distribution is faster than R's built-in.

compute_d_tilde_half <- function(J, W_prime, index_row = NA, W = NA,
                                 d_tilde_half = NA, w_diff = NA,
                                 full_computation = TRUE) {
  # Follow Equation 4 from Grathwohl et al. but excludes the possibility
  # of proposing the same state again.
  # For efficiency, update previous d_tilde_half!
  # If full_computation is TRUE, then previous d_tilde_half is not used.
  # NOTE: based on unit tests, the speed gain from using stored d_tilde
  # is marginal at best. TRUE is recommended.
  #
  # J: Scaled adjacency matrix
  # W_prime: one-hoc representation of the Potts particles.
  # index_row: which particle is updated in the proposal.
  # W: one-hoc representation of the Potts particles at previous iter.
  # d_tilde_half: obtained using W at previous iteration.

  infinity = 10^10
  n = nrow(W_prime)
  q = ncol(W_prime)
  if (!full_computation) index_col_prime = which_one(W_prime[index_row, ])
  if (!full_computation) index_col = which_one(W[index_row, ])
  d_tilde_half_prime = array(NA, dim = c(n, q))
  d_tilde_half_corrected = array(NA, dim = c(n, q))
  
  for (i in 1:n) {
    if (full_computation) {
      JW_prime = J[i, ] %*% W_prime
      d_tilde_half_prime[i, ] = JW_prime - rep(JW_prime %*% W_prime[i, ], q)      
    } else {
      if (i != index_row) {
        d_tilde_half_prime[i, ] = d_tilde_half[i, ] +
           J[i, index_row] * (w_diff -
                              rep(W[i, index_col_prime] - W[i, index_col], q))
      } else {
        d_tilde_half_prime[i, ] = d_tilde_half[i, ] +
          rep(J[i, ] %*% (W_prime[, index_col_prime] - W[, index_col]), q)
      }
    }

    d_tilde_half_corrected[i, ] = d_tilde_half_prime[i, ]

    # exclude possibility of proposing the same state
    l_prime = which_one(W_prime[i, ])
    d_tilde_half_corrected[i, l_prime] = - infinity
  }

  # CHECK -- should there be a 0.5?
  return(list(d_tilde_half_prime = d_tilde_half_prime, 
              d_tilde_half_corrected = d_tilde_half_corrected))
}

gwg_potts <- function(A, beta, init, n_states = 2, n_iter = 1000) {
  n_particles = nrow(A)
  J = beta * A # J = 0.5 * beta * A -- CHECK

  W = array(0, dim = c(n_particles, n_states))
  for (l in 1:n_states) W[init == l, l] = 1

  X = array(NA, dim = c(n_iter, n_particles))

  # initialize variables which carry over if a proposal is rejected.
  d_list = compute_d_tilde_half(J, W, full_computation = TRUE)
  d_tilde_half = d_list$d_tilde_half_prime
  d_tilde_half_corrected = d_list$d_tilde_half_corrected

  P = softmax(d_tilde_half_corrected)
  f = 0
  for (i in 2:n_particles) {
    for (j in 1:(i - 1)) {
      f = f + J[i, j] * t(W[i, ]) %*% W[j, ]
    }
  }
  f = 2 * f
  x_current = init

  for (iter in 1:n_iter) {
    # draw index to update from categorical
    u = runif(1, 0, 1)
    index_row = 1; index_col = 1
    p_current = P[index_row, index_col]
    while (u > p_current) {
      if (index_row < n_particles) {
        index_row = index_row + 1
      } else {
        index_col = index_col + 1
        index_row = 1
      }
      p_current = p_current + P[index_row, index_col]
    }

    W_prime = W
    W_prime[index_row, ] = 0
    W_prime[index_row, index_col] = 1
    
    # Compute W' - W for flipped particle
    w_diff = rep(0, n_states)
    w_diff[index_col] = 1
    w_diff[which(W[index_row, ] == 1)] = -1

    d_list = compute_d_tilde_half(J, W_prime, index_row = index_row,
                                  W = W, d_tilde_half = d_tilde_half,
                                  w_diff = w_diff,
                                  full_computation = FALSE)
    d_tilde_half_prime = d_list$d_tilde_half_prime
    d_tilde_half_corrected = d_list$d_tilde_half_corrected

    exp_d_tilde_half = exp(d_tilde_half_corrected)
    P_prime = exp_d_tilde_half / sum(exp_d_tilde_half)

    # Reuse terms calculated for f to get f_proposal.
    f_proposal = f
    for (j in 1:n_particles) {
      f_proposal = f_proposal + 2 * J[index_row, j] *
                     t(w_diff) %*% W_prime[j, ]
    }

    index_col_past = which(W[index_row, ] == 1)
    prob_accept = exp(f_proposal - f) *
                     P_prime[index_row, index_col_past] / P[index_row, index_col]

    u = runif(1, 0, 1)
    if (u < prob_accept) {
      x_current = rep(NA, n_particles)
      for (n in 1:n_particles) x_current[n] = which(W_prime[n, ] == 1)
      W = W_prime
      P = P_prime
      f = f_proposal
      d_tilde_half = d_tilde_half_prime
    }

    X[iter, ] = x_current
  }
  
  return (X)
}

gwg_sampler <- function(A, beta, init, n_states, n_iter) {
  if (n_states == 2) {
    X = gwg_ising(A, beta, init, n_iter)
  } else {
    X = gwg_potts(A, beta, init, n_states, n_iter)
  }
  
  return(X)
}


###############################################################################
# compute_d_tilde <- function(By, n_states, n_particles, B, y) {
#   d_tilde = rep(NA, n_particles * n_states)
#   for (l in 1:n_states) {
#     index = (n_particles * (l - 1) + 1):(n_particles * l)
#     By[, l] = B %*% y[index]
#     d_tilde[index] = - 2 * y[index] * By[, l]
#   }
# 
#   d_tilde
# }
# 
# gwg_potts <- function(A, beta, init, n_states = 2, n_iter = 1000) {
#   n_particles = nrow(A)
#   
#   y = array(0, dim = c(n_particles, n_states))
#   for (l in 1:n_states) y[init == l, l] = 1
#   y = c(y)
#   
#   Y = matrix(NA, nrow = n_iter, ncol = n_particles * n_states)
#   Y[1, ] = y
#   B = 0.5 * beta * A
#   
#   for (i in 1:n_iter) {
#     # Use sparsity of augmented B to compute d_tilde
#     By = array(NA, dim = c(n_particles, n_states))
#     d_tilde = compute_d_tilde(By, n_states, n_particles, B, y)
# 
#     p = softmax(0.5 * d_tilde)
# 
#     # draw index to update from categorical
#     u = runif(1, 0, 1)
#     index = 1
#     p_current = p[index]
#     while (u > p_current) {
#       index = index + 1
#       p_current = p_current + p[index]
#     }
# 
#     y_proposal = y
#     y_proposal[index] = - y[index]
#     
#     By_proposal = array(NA, dim = c(n_particles, n_states))
#     d_tilde_prop = compute_d_tilde(By_proposal, n_states, n_particles, B, 
#                                    y_proposal)
#     p_proposal = softmax(0.5 * d_tilde_prop)
#     
#     # TODO: reuse the calculate unnormalized density from the previous proposal!!
#     log_f = 0
#     log_f_proposal = 0
#     for (l in 1:n_states) {
#       index = (n_particles * (l - 1) + 1):(n_particles * l)
#       log_f = log_f + y[index] %*% By[, l]
#       log_f_proposal = log_f_proposal + y_proposal[index] %*% By_proposal[, l]
#     }
# 
#     prob_accept = exp(lof_f_proposal - log_f) * p_proposal[index] / p_current[index]
#     
#     u = runif(1, 0, 1)
#     if (u < prob_accept) {
#       
#     }
#     
#   }
# }





## Test code

# x <- init_potts
# H = 0
# for (i in 1:n_part) {
#   for (j in 1:n_part) {
#     H = H + (x[i] == x[j]) * C[i, j]
#   }
# }


# 1) Reduce rank
# n_particles <- 16
# A <- adjacency_graph(n_particles, "complete")
# eigen_output <- eigen(A, symmetric = TRUE)
# epsilon <- 1e-12
# eigenvalue_rounding <- 1e12
# n_rank = length(round(
#   eigen_output$values[abs(eigen_output$values) > epsilon] *
#     eigenvalue_rounding) / eigenvalue_rounding)
# 
# # Nake sure we have eigenvalues and eigenvectors
# # index <- 4
# # t(A %*% eigen_output$vectors[, index])
# # t(eigen_output$values[index] * eigen_output$vectors[, index])
# 
# A_tilde = matrix(0, nrow = n_particles, ncol = n_particles)
# for (i in 1:n_rank) {
#   A_tilde = A_tilde +
#     eigen_output$values[i] * eigen_output$vectors[, i] %*% t(eigen_output$vectors[, i])
# }
# # QUESTION: why does A_tilde not recover A??
# 
# alpha <- abs(min(eigen(A_tilde, symmetric = TRUE, only.values = TRUE)$values))
# C_tilde <- A_tilde + alpha * diag(n_particles)
# init <- sample(1:n_states, n_particles, replace = T)
# n_states = 4

# B_tilde <- eigen_output$values + eigen_output$vectors[, 1] %*% t(eigen_output$vectors[, 1])
# eigen_output$values[3] * eigen_output$vectors[, 3]
# A %*% eigen_output$vectors[, 3]


  # # NOTE: The Auxiliary variable I introduce is an overcomplication.
  # # The ag_potts_simple for a simpler version.
  # # DEPRECATE.
  # ag_potts <- function (A, beta, init, n_states = 2, n_iter = 1000, B = 0,
  #                       alpha_modif = 1e-7, retain_alpha = FALSE,
  #                       perturb = FALSE) {
  #   print("WARNING: We recommend using ag_potts_simple.")
  #   if (retain_alpha) {
  #     alpha = alpha_modif
  #   } else {
  #     # CHECK ME!
  #     alpha = abs(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
  #     alpha = alpha + alpha_modif
  #   }
  #   
  #   n_particles = nrow(A)
  #   C = beta * (A + alpha * diag(n_particles))
  #   L = 2 * t(chol(solve(C)))
  #   X = matrix(NA, nrow = n_iter, ncol = n_particles)
  #   x = init
  #   X[1, ] = x
  #   
  #   # Auxiliary Gaussians
  #   Q = matrix(NA, nrow = n_particles, ncol = n_states)
  #   
  #   for (i in 2:n_iter) {
  #     for (l in 1:n_states) {
  #       mu_q = rep(-1, n_particles)
  #       mu_q[X[i - 1, ] == l] = 1
  #       Q[, l] = mvrnorm_cholesky(n = 1, mu = mu_q, L = L)
  #     }
  #     
  #     P = C %*% Q
  #     P_row_sum_term = 0.25 * rowSums(P)
  #     P_half = 0.5 * P
  #     
  #     u = runif(n_particles, 0, 1)
  #     
  #     for (j in 1:n_particles) {
  #       p = rep(NA, n_states)
  #       for (l in 1:n_states) {
  #         p[l] = exp(P_half[j, l] - P_row_sum_term[j])
  #       }
  #       norm_const = sum(p)
  #       
  #       state = 1
  #       p_current = p[state] / norm_const
  #       while (u[j] > p_current) {
  #         state = state + 1
  #         p_current = p_current + p[state] / norm_const
  #       }
  #       
  #       X[i, j] = state
  #     }
  #   }
  #   X
  # }  
