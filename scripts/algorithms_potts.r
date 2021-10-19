
# Algorithms_potts
# Extend algorithms to work on pott sampler.

# library(extraDistr)  # CHECK -- do we need this?
source("tools.r")

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
  # spectral_decomposition = eigen(2 * beta * A, symmetric = TRUE)
  # n_rank = length(
  #   spectral_decomposition$values[abs(spectral_decomposition$values) > epsilon])
  # 
  # A_tilde = matrix(0, nrow = n_particles, ncol = n_particles)
  # for (i in 1:n_rank) {
  #   A_tilde = A_tilde +
  #     spectral_decomposition$values[i] * spectral_decomposition$vectors[, i] %*% 
  #       t(spectral_decomposition$vectors[, i])
  # }
  # 
  # alpha <- abs(min(eigen(A_tilde, symmetric = TRUE, only.values = TRUE)$values))
  # C_tilde <- A_tilde + alpha * diag(n_particles)

  alpha = abs(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
  alpha = alpha + alpha_modif
  C = 2 * beta * (A + alpha * diag(n_particles))

  ag_potts_lowrank__(C, init, n_states, n_iter, epsilon)
}


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
