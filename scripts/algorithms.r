# Algorithms
# Source this file to load simulation algorithms for Ising model.
# Includes:
#  - Auxiliary Gaussian algorithms
#  - Heat Bath / Metropolis-Hastings
#  - Swendsen-Wang
#  - Wolff
#  - Exact (brute force, does not scale)
#
# Remark: this file uses functions in tools.r
# Source this file before running tests.

source("tools.r")

#####################################################################
## Sumit's algorithm (more formally, Auxiliary Gaussian algorithm
# for Ising model.)
# Note: this special case reduces to the algorithm by Martens & Sutskever (2010))

construct_C <- function(A, beta, alpha_modif = 1e-7, retain_alpha = FALSE) {
  # Make an adjacency matrix A positive definite but incrementing
  # the diagonal element. Returns energy of Ising model.
  #
  # A: adjacency matrix
  # beta: inverse temperature
  # alpha_modif: amount by which to increment diagonal element if
  #              retain_alpha = TRUE
  # retain_alpha: if TRUE, use alpha_modif. If FALSE, calculate the smallest
  #               absolute eigenvalue and use this to increment the diagonal
  #               elements of A. alpha_modif is also added (but can be set to
  #               a jitter).
  
  if (retain_alpha) { alpha = alpha_modif } else {
    alpha = abs(min(eigen(A, symmetric = TRUE,
                          only.values = TRUE)$values))
    alpha = alpha + alpha_modif 
  }
  
  n_particles = nrow(A)
  C = beta * (A + alpha * diag(n_particles))
  
  return(C)
}


sumit_sampler <- function (A, beta, B = 0, init, n_iter = 1000,
                           alpha_modif = 1e-7,
                           retain_alpha = FALSE,
                           perturb = TRUE) {
  # Also known as the Auxiliary Gaussian sampler.
  # Simulate states from the Ising model, using an auxiliary
  # Gaussian variable, q, for the Ising model.
  #
  # A: adjacency matrix
  # beta: inverse temperature
  # init: starting point for Markov chains
  # n_iter: length of Markov chain
  # alpha_modif: jitter or amount by which to increment diagonal element of A.
  # retain_alpha: if TRUE, only use alpha_modif, else use eigenvalue of A.
  # perturb: if TRUE, add a flipping step at each iteration.
  #
  # NOTE: this algorithm requires x in {-1, 1}.

  n_particles = nrow(A)
  C = construct_C(A, beta, alpha_modif, retain_alpha)
  L = t(chol(solve(C)))
  X = matrix(NA, nrow = n_iter, ncol = n_particles)
  x = init
  X[1, ] = x

  for (i in (1:(n_iter - 1))) {
    # update q
    q = mvrnorm_cholesky(n = 1, mu = x, L = L)

    # update x
    q_C = q %*% C
    u = runif(n_particles, 0, 1)
    for (j in 1:n_particles) {
      p = 1 / (1 + exp(-2 * (q_C[j] + B)))
      if (p > u[j]) { x[j] = 1 } else { x[j] = -1 }
    }

    # Random perturbator.
    if (perturb) {
      prob_flip = exp(- 2 * B * sum(x))
      u = runif(1, 0, 1)
      perturbator = - 2 * (u <= prob_flip) + 1
      x = perturbator * x
      X[i + 1, ] = x
    } else {
      X[i + 1, ] = x
    }
  }
  X
}

ag_ising_lowrank <- function (A, beta, init, n_iter, epsilon = 1e-12) {
  # Specialized AG algorithm using lowrank structure of A.
  #
  # A: adjacency matrix
  # beta: inverse temperature
  # init: starting point for Markov chains
  # n_iter: length of Markov chain
  # epsilon: threshold under which an eigenvalue is considered to be 0,
  #          when doing low-rank approximation.
  
  n_particles = nrow(A)
  alpha = abs(min(eigen(A, symmetric = TRUE,
                        only.values = TRUE)$values))
  C = beta * (A + alpha * diag(n_particles))
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
  
  q = rep(NA, n_rank)

  for (i in 2:n_iter) {
    mu_q = t(spectral_decomposition$vectors[, 1:n_rank]) %*% x
    q = rnorm(n_rank) / eigen_values_sqrt + mu_q

    P = scaled_eigenvectors %*% q

    u = runif(n_particles, 0, 1)
    for (j in 1:n_particles) {
      # p = 1 / (1 + exp(-2 * scaled_eigenvectors[j, ] %*% q))
      p = 1 / (1 + exp(-2 * P[j]))
      if (p > u[j]) { x[j] = 1 } else { x[j] = -1 }
    }
    X[i, ] = x
  }

  X
}

#####################################################################
## Metropolis-Hastings

mh_sampler <- function (A, beta, B = 0, 
                        init, n_iter, perturb = FALSE) {
  # A: adjacency matrix
  # beta: inverse temperature
  # B: magnetic field
  # init: starting point for Markov chains
  # n_iter: length of Markov chain
  # perturb: if TRUE, add a flipping step at each iteration.

  n_particles = nrow(A)
  X = matrix(NA, nrow = n_iter, ncol = n_particles)
  x = init
  X[1, ] = x

  kernel_x = ising_kernel(beta, A, x, B)

  for (i in 2:n_iter) {
    # randomly flip one of the spins
    x_new = x
    index <- sample(length(x), 1)
    x_new[index] = - x[index]

    kernel_new = ising_kernel(beta, A, x_new, B)
    a = kernel_new / kernel_x
    u = runif(1, 0, 1)
    if (u < a) {
      x = x_new
      kernel_x = kernel_new
    }

    if (perturb) {
      prob_flip = exp(- 2 * B * sum(x))
      u = runif(1, 0, 1)
      perturbator = - 2 * (u < prob_flip) + 1
      x = perturbator * x
      kernel_x = exp(log(kernel_x) - 2 * B * sum(x))
    } else {
      x = x
    }

    X[i, ] = x
  }
  X
}

#####################################################################
## Heat bath

hb_sampler <- function(A, beta, B_field = 0, init, n_iter,
                       sub_sample = nrow(A),
                       is_potts = FALSE, n_states = 2,
                       black_box = FALSE) {
  # A: adjacency matrix
  # beta: inverse temperature
  # B_field: magnetic field
  # init: starting point for Markov chains
  # n_iter: length of Markov chain
  # is_potts: if TRUE, the state space is over {0, 1,...}, and if
  #           FALSE, use {-1, 1}.
  # n_states: number of states q (2 for Ising, > 2 for Potts)
  # black_box: if TRUE, do not use any details of the structure in
  #            the energy to update the energy from one iteration to
  #            the other (i.e. the model is treated as a black box).
  #            This makes the algorithm slower. 
  
  n_particles = nrow(A)
  X = matrix(NA, nrow = n_iter, ncol = n_particles)
  x = init
  X[1, ] = x

  if (!is_potts) {
    B = 0.5 * beta * A
    if (black_box) kernel_x = ising_kernel(beta, A, x, B_field, log = T)
  } else {
    B = beta * A  # CHECK
    if (black_box) kernel_x = potts_log_kernel(beta = beta, A = A, x = x, B = B_field)
  }

  for (i in 2:n_iter) {
    # randomly sample which particles get updated
    update_particles = sample(1:n_particles, sub_sample, replace = FALSE)

    if (!is_potts) {  # Ising model with x in {-1, 1}
      for (j in update_particles) {
        x_new = x
        x_new[j] = - x_new[j]
        
        if (black_box) {
          kernel_new = ising_kernel(beta, A, x_new, B_field, log = T)
          a = 1 / (1 + exp(kernel_x - kernel_new))
        } else {
          # Trick to avoid computing the kernel at each iteration
          proposal_diff = x_new[j] - x[j]
          kernel_diff = 2 * sum(B[j, ] * x * proposal_diff)
          a = 1 / (1 + exp(- kernel_diff))
        }

        u = runif(1, 0, 1)
        if (u < a) {
          x = x_new
          if (black_box) kernel_x = kernel_new
        }
      }
    } else {  # Potts model with x in {1, 2, ..., n_states}
      for (j in update_particles) {
        x_new = x
        other_states <- (1:n_states)[-x[j]]
        x_new[j] = other_states[sample(n_states - 1, 1)]

        if (black_box) {
          kernel_new = potts_log_kernel(beta = beta, A = A, x = x_new,
                                        n_states = n_states)
          a = 1 / (1 + exp(kernel_x - kernel_new))
        } else {
          # Trick to avoid computing the kernel at each iteration
          kernel_diff = 2 * sum(B[j, ] * ((x_new[j] == x_new) - (x[j] == x)))
          a = 1 / (1 + exp(-kernel_diff))
        }

        u = runif(1, 0, 1)
        if (u < a) {
          x = x_new
          if (black_box) kernel_x = kernel_new
        }
      }
    }
    X[i, ] = x
  }
  X
}

#####################################################################
## Swendsen-Wang

sw_sampler <- function (A, beta, B = 0, init, n_iter, wolff = FALSE,
                        is_potts = FALSE, n_states = 2) {
  # Simulates samples from an ising model, using the
  # Swendsen-Wang algorithm. This algorithm introduces
  # an auxiliary random cluster.
  # The algorithm is written for a complete graph, meaning
  # every particle can be connected to any other particle,
  # when constructing the random cluster (i.e. the algorithm is
  # agnostic to whether we're sampling on a grid or not).
  # 
  # If Wolff = TRUE, use Wolff's algorithm instead of SW.
  # That is, for each iteration, construct only one cluster,
  # started at a random point.
  #
  # A: adjacency matrix
  # beta: inverse temperature
  # B: magnetic field
  # init: starting point for Markov chains
  # n_iter: length of Markov chain
  # wolff: if TRUE, only construct one auxliary cluster, per
  #        the Wolff algorithm.
  # is_potts: if TRUE, the state space {1, 2, ..., n_states}.
  #           if FALSE, the state space is {-1, 1}
  # n_states: number of states q.
  #
  # IDEA - use a lower triangular matrix instead of a full
  # matrix for the edges.
  n_particles = nrow(A)
  n_edges = n_particles * (n_particles - 1) / 2
  X = matrix(NA, nrow = n_iter, ncol = n_particles)
  x = init
  X[1, ] = x
  edges = matrix(NA, nrow = n_particles, ncol = n_particles)

  # Remark: p[i, j] = 0 if A[i, j] = 0.
  p = 1 - exp(-2 * beta * A)
  
  if (is_potts) p_states = 1 / n_states  # Assuming B = 0

  for (i in 1:(n_iter - 1)) {
    # Compute auxiliary random clusters.
    # j and k respectively index rows and columns in the
    # edge matrix.
    for (j in 2:n_particles) {
      for (k in 1:(j - 1)) {
        if (x[j] != x[k] | A[j, k] == 0) { edges[j, k] = 0
        } else {
          u = runif(1, 0, 1)
          if (u < p[j, k]) { edges[j, k] = 1 } else { edges[j, k] = 0 }
        }
      }
    }

    # Identify clusters and assign a spin.
    free_nodes = 1:n_particles
    j = 1
    if (wolff == T) wolff_node = sample(1:n_particles, 1)
    while (length(free_nodes) != 0) {
      if (wolff == T) j = wolff_node

      # Can only start a new cluster with a free node.
      cluster = c(j)
      cluster_grows = TRUE
      free_nodes = free_nodes[-match(j, free_nodes)]

      while (cluster_grows) {
        cluster_grows = FALSE
        potential_connection = free_nodes

        # Only expand the cluster if there are potential connections.
        if (length(potential_connection) != 0) {
          for (k in potential_connection) {
            if (edges[max(k, j), min(k, j)] == 1) {
              cluster = c(cluster, k)
              free_nodes = free_nodes[-match(k, free_nodes)]
            }
          }
        }
        if (length(cluster) > match(j, cluster)) {
          cluster_grows = T
          j = cluster[match(j, cluster) + 1]
        }
      }

      u = runif(1, 0, 1)
      
      if (!is_potts) {  # Ising model with x in {-1, 1}
        p_flip = 1 / (1 + exp(- 2 * length(cluster) * B))
        if (u < p_flip) { x[cluster] = 1 } else { x[cluster] = -1 }
        # print(paste0("cluster: ", cluster))  # use for tests
      } else {  # Potts models with x in {1, 2, ..., n_states}
        p_current = p_states
        state = 1
        while (u > p_current) {
          state = state + 1
          p_current = p_current + p_states
        }
        x[cluster] = state
      }

      if (length(free_nodes) != 0) j = min(free_nodes)
      if (wolff == T) break
    }

    X[i + 1, ] = x
  }
  X
}

wolff_sampler <- function (A, beta, B = 0, init, n_iter,
                           is_potts = FALSE, n_states = 2) {
  # Wolff algorithm. This is a variation of S-W that only
  # pertubs one cluster per iteration (and as a consequence,
  # it is easier to code).
  #
  # A: adjacency matrix
  # beta: inverse temperature
  # B: magnetic field
  # init: starting point for Markov chains
  # n_iter: length of Markov chain
  # is_potts: if TRUE, the state space {1, 2, ..., n_states}.
  #           if FALSE, the state space is {-1, 1}
  # n_states: number of states q.

  n_particles = nrow(A)
  X = matrix(NA, nrow = n_iter, ncol = n_particles)
  x = init
  X[1, ] = x

  # Remark: p[i, j] = 0 if A[i, j] = 0.
  p = 1 - exp(-2 * beta * A)
  
  if (is_potts) p_states = 1 / n_states  # Assuming B = 0

  for (i in 1:(n_iter - 1)) {
    source_node = sample(1:n_particles, 1)
    free_nodes = (1:n_particles)[-source_node]
    cluster = c(source_node)
    j = source_node
    cluster_grows = T

    while (cluster_grows) {
      cluster_grows = F
      for (k in free_nodes) {
        # compute whether there is an edge or not
        u = runif(1, 0, 1)
        edge = (x[j] == x[k]) & (u < p[j, k])
        if (edge) {
          # cluster_grows = T
          cluster = c(cluster, k)
          free_nodes = free_nodes[-match(k, free_nodes)]
        }
      }
      if (length(cluster) > match(j, cluster)) {
        cluster_grows = T
        j = cluster[match(j, cluster) + 1]
      }
      # if (cluster_grows) j = cluster[match(j, cluster) + 1]
    }
    u = runif(1, 0, 1)
    # p_flip = 1 / (1 + exp(- 2 * length(cluster) * B))
    if (!is_potts) {  # Ising model with x in {-1, 1}
      p_flip = exp(- 2 * B * sum(x[cluster]))
      if (u < p_flip) { x[cluster] = -x[cluster]}
    } else {  # Potts models with x in {1, 2, ..., n_states}
      p_current = p_states
      state = 1
      while (u > p_current) {
        state = state + 1
        p_current = p_current + p_states
      }
      x[cluster] = state
    }

    X[i + 1, ] = x
  }
  X
}

###############################################################################

gwg_ising <- function(A, beta, init, n_iter) {
  # Gibbs with gradient (Grathwohl et al 2021, https://arxiv.org/abs/2102.04509),
  # implemented with no exterior magnetic field.
  #
  # A: adjacency matrix
  # beta: inverse temperature
  # init: starting point for the Markov chains
  # n_iter: length of Markov chain
  #
  # NOTE: hand-coded categorical distribution is more efficient than R's
  #       implementation.
  
  n_particles = nrow(A)
  X = matrix(NA, nrow = n_iter, ncol = n_particles)
  x = init
  X[1, ] = x
  B = 0.5 * beta * A
  
  # initialize variables which carry over if a proposal is rejected
  Bx = B %*% x
  d_tilde_half = - x * Bx
  p = softmax(d_tilde_half)
  f = x %*% Bx

  for (i in 2:n_iter) {
    # draw index to update from categorical
    u = runif(1, 0, 1)
    index = 1
    p_current = p[index]
    while (u > p_current) {
      index = index + 1
      p_current = p_current + p[index]
    }

    x_proposal = x
    x_proposal[index] = - x[index]

    Bx_proposal = Bx + B[, index] * (x_proposal[index] - x[index])
    d_tilde_half_prop = - x_proposal * Bx_proposal
    exp_d_tilde_half = exp(d_tilde_half_prop)
    p_proposal = exp_d_tilde_half / sum(exp_d_tilde_half)
    # p_proposal = softmax(d_tilde_half_prop)

    f_proposal = x_proposal %*% Bx_proposal

    prob_accept = exp(f_proposal - f) * p_proposal[index] / p[index]

    u = runif(1, 0, 1)
    if (u < prob_accept) {
      x[index] = - x[index]
      p = p_proposal
      f = f_proposal
      Bx = Bx_proposal
    }

    X[i, ] = x
  }

  X
}

exact_sampler <- function (p, n_particles, n_iter) {
  # Brute force simulation, using the exact probability for each state.
  # The states are ordered using binary basis.
  #
  # p: probabilities for each state
  # n_particles: number of particles
  # n_iter: length of Markov chains.
  
  X = matrix(NA, nrow = n_iter, ncol = n_particles)
  for (i in 2:length(p)) p[i] = p[i - 1] + p[i]
  
  for (i in 1:n_iter) {
    u = runif(1, 0, 1)
    state = which(abs(p - u) == min(abs(p - u)))
    if (p[state] - u < 0) state = state + 1
    
    # remember we count states, starting at 0!
    X[i, ] = to_binary(state - 1, n_particles, ising_state = TRUE)
  }
  
  X
}

###############################################################################
## DRAFT CODE

###############################################################################
## Potts sampler
potts_sampler <- function(A, theta, init, n_sample) {
  # Uses the R package Potts sampler by (Okabayashi, Johnson and Geyer, 2011)
  # Assume we have a squared grid.
  # Convert from {-1, 1} to {2, 1} to respect Potts convention.
  # NOTE -- we may be waisting time doing a lot of copying back
  # and forth, so don't use this measure algorithm run time.
  #
  # theta is canonical parameter, equals 2 * beta.
  ncolor <- as.integer(2)
  n_particles <- nrow(A)
  init_color <- init
  init_color[init_color == -1] <- 2
  x_init <- packPotts(matrix(init_color,
                             nrow = sqrt(n_particles),
                             ncol = sqrt(n_particles)),
                      ncolor = 2)

  samples <- matrix(NA, n_sample, n_particles)

  # set canonical parameter. No B-field, so "color parameters"
  # need to be equal.
  parameter <- c(rep(0, 2), theta)

  for (i in 1:n_sample) {
    if (i %% 100 == 0) print(paste0("iteration: ", i, " / ", n_sample))
    if (i == 1) out <- x_init
    out <- potts(out, parameter, nbatch = 1,
                 boundary = "free")

    samples[i, ] <- as.vector(t(unpackPotts(out$final)))
    samples[i, samples[i, ] == 2] <- -1
  }
  samples
}
