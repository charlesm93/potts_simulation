# Convenient functions for evaluating ising model simulators.
library(knitr)

#####################################################################
## Functions to order Ising states.

to_binary <- function(n, n_digits, ising_state = FALSE) {
  # Convert base 10 number to base 2 number, with digits
  # stored as an n_digits-vector.
  #
  # n: the state number (start counting at 0)
  # n_digitis: the number of digits, hence of particles in the graph.
  # ising_state: if TRUE, returns a vector 1 and -1, instead of 1 and 0's.
  j = 0
  x = rep(ifelse(ising_state, -1, 0), n_digits)
  while (n > 0) {
    x[n_digits - j] = n %% 2
    if (ising_state && (x[n_digits - j] == 0)) x[n_digits - j] = ifelse(ising_state, -1, 0)
    n = as.integer(n/2)
    j = j + 1
  }
  x
}

to_decimal <- function(samples) {
  # convert samples to decimal numbers
  n_iter = nrow(samples)
  n_particles = ncol(samples)
  decimals <- rep(NA, n_iter)
  for (i in 1:n_iter) {
    decimals[i]  = 0
    for (j in 1:ncol(samples)) decimals[i] = decimals[i] +
        2^{n_particles - j} * samples[i, j] * (samples[i, j] == 1)
  }
  decimals
}

# returns the number of samples with a spin up (+1)
count_states <- function(samples, n_part) {
  0.5 * (rowSums(samples) + n_part)
}

# same as above, but for a single sample
count_states_single <- function(samples, n_part) {
  0.5 * (sum(samples) + n_part)
}

###############################################################################
## Functions to do an exact computation of the partition function.

# Returns the kernel for an ising model.
ising_kernel <- function (beta, A, x, B_field = 0, log = FALSE) {
  hamiltonian = 0.5 * beta * (t(x) %*% A %*% x) + B_field * sum(x)
  if (log) { hamiltonian } else { exp(hamiltonian) }
}

# use this in the tempering algorithm to not double-count the beta.
raw_hamiltonian <- function (A, x) {
  0.5 * t(x) %*% A %*% x
}

potts_log_kernel <- function (beta, A, x, B = 0, n_states = 2) {
  n_particle = nrow(A)
  exp_term = 0
  Y = matrix(NA, nrow = n_particle, ncol = n_states)
  for (l in 1:n_states) {
    Y[, l] = rep(0, n_particle)
    Y[x == l, l] = 1
    exp_term = exp_term + t(Y[, l]) %*% A %*% Y[, l]
  }

  # for (j in 2:n_particle) {
  #   for (i in 1:(j - 1)) {
  #     if (x[i] == x[j]) exp_term = exp_term + A[i, j]
  #   }
  # }

  beta * exp_term
}

partition <- function (beta, A, B = 0, prob = FALSE) {
  # brute-force computes the partition function for
  # an Ising model. Useful for simple configurations
  # and unit tests.
  # If prob = TRUE, return a list which contains the 
  # partition function and the probability mass for
  # each configuration.
  n_particles = nrow(A)
  x = rep(-1, n_particles)
  Z = 0
  if (prob) p = rep(NA, 2^n_particles)
  
  for (i in 0:((2^n_particles) - 1)) {
    # update state of the system
    x = to_binary(i, n_particles, ising_state = TRUE)
    local_kernel = ising_kernel(beta, A, x, B)
    Z = Z + local_kernel
    if (prob) p[i + 1] = local_kernel
  }
  
  Z = as.numeric(Z)
  
  if (!prob) { Z } else {
    list(Z = Z, p = p / Z)
  }
}

partition_complete <- function (beta, A, B = 0, prob = FALSE) {
  # a specialized function that computes the partition function
  # for a complete graph.
  # Remark: the probability are over the summary count space.
  n_particles = nrow(A)
  d = A[1, 2]

  # number of positive spins
  # n = 0:(floor(n_particles / 2))
  n = 0:n_particles
  length_n = length(n)
  # double_count = rep(2, length_n)
  # if (n_particles %% 2 == 0) {double_count[length_n] = 1}

  p = rep(NA, n_particles + 1)
  for (i in n) {
    j = i + 1
    p[j] = choose(n_particles, i) * exp(0.5 * beta * d *
                                        (4 * i^2 - 4 * i * n_particles +
                                        n_particles^2 - n_particles) +
                                        B * (2 * i - n_particles))
    # if (double_count[j] == 2) p[n_particles + 1 - (j - 1)] = p[j]
  }
  
  Z = sum(p)
  if (!prob) { Z } else {
    list(Z = Z, p = p / Z)
  }
}

count_probability <- function(p_analytical, n_particles) {
  # Returns probability over the summary count space, given
  # the exact probability for each state.
  # Useful in non-symmetric cases.
  
  n_states = length(p_analytical)
  # for (i in 0:n_states) {
  prob = rep(0, n_particles + 1)
  for (i in 0:(n_states - 1)) {
    index = sum(to_binary(i, n_particles, ising_state = FALSE)) + 1
    prob[index] = prob[index] + p_analytical[i + 1]
  }
  prob
}

#####################################################################
# Functions which return performance metrics, based on the
# produced samples.

summary <- function (samples, p_analytical, discard_burn_in = TRUE) {
  # Returns a few summaries statistics, including accuracy ratio
  # and variational distance.
  #
  # samples: simulated samples with x in {-1, 1}^n
  # p_analytical: a vector of probability
  # discard_burn_in: if TRUE, only use the second half of the samples.
  
  n_particles = ncol(samples)
  n_iter = nrow(samples)
  n_states = 2^n_particles
  
  # convert results to decimal numbers
  decimals = to_decimal(samples)
  if (discard_burn_in) decimals = decimals[(n_iter / 2):n_iter]
  
  p_empirical = rep(NA, n_states)
  for (i in 0:(n_states - 1)) p_empirical[i + 1] = sum(decimals == i)
  p_empirical = p_empirical / length(decimals)
  
  # p_empirical = table(decimals) / length(decimals)
  
  accuracy_ratio = as.vector(abs(p_empirical / p_analytical - 1))
  max_ratio = max(accuracy_ratio)
  total_distance = 0.5 * sum(abs(p_empirical - p_analytical))
  variational_distance = max(abs(p_empirical - p_analytical))
  
  list(decimals = decimals,
       accuracy_ratio = accuracy_ratio,
       max_ratio = max_ratio,
       total_distance = total_distance,
       variational_distance = variational_distance)
}

gen_frequency_iter <- function (decimals, n_particles) {
  # Returns a 2-array with the frequency of each state
  # after the i^th iteration.
  # WARNING: this function takes some time to run!
  # For loops in R are very inefficient.
  # FIX ME -- find a slicker way of writing this.

  n_iter = length(decimals)
  n_states = 2^n_particles

  frequency_iter = matrix(0, nrow = n_iter, ncol = 2^n_particles)
  
  for (j in 1:n_iter) {
    if (j %% 100 == 0) print(paste0("iteration: ", j, " / ", n_iter))
    unique_elements = sort(unique(decimals[1:j])) + 1  # adjust index
    current_frequency = table(decimals[1:j]) / j
    
    if (length(current_frequency) == n_states) {
      frequency_iter[j, ] = current_frequency
    } else {
      for (i in 1:length(unique_elements)) {
        frequency_iter[j, unique_elements[i]] = current_frequency[i]
      }
    }
  }
  
  frequency_iter
}

ttv_total <- function (samples, p_analytical, discard_burn_in = TRUE) {
  # Return the total variational distance (ttv). 
  n_particles = ncol(samples)
  n_iter = nrow(samples)
  n_states = 2^n_particles
  
  # convert results to decimal numbers
  decimals = to_decimal(samples)
  if (discard_burn_in) decimals = decimals[(n_iter / 2):n_iter]
  
  p_empirical = rep(NA, n_states)
  for (i in 0:(n_states - 1)) p_empirical[i + 1] = sum(decimals == i)
  p_empirical = p_empirical / length(decimals)
  
  0.5 * sum(abs(p_empirical - p_analytical))
}

ttv_count <- function (samples, p_analytical, discard_burn_in = TRUE) {
  # Return the variational distance over the count space
  # (i.e. the number of up spins)

  n_particles = ncol(samples)
  n_iter = nrow(samples)
  n_states = n_particles + 1

  # convert results to spin up counts.
  if (discard_burn_in) samples = samples[(n_iter / 2):n_iter, ]
  count_sample = count_states(samples, n_particles)

  p_empirical = rep(NA, n_states)
  for (i in 0:(n_states - 1)) p_empirical[i + 1] = sum(count_sample == i)
  p_empirical = p_empirical / length(count_sample)
  
  0.5 * sum(abs(p_empirical - p_analytical))
}

ttv_comp <- function(samples, p_analytical, type = "total", 
                     discard_burn_in = T) {
  if (type == 'total') 
    return(ttv_total(samples, p_analytical, discard_burn_in))
  if (type == 'count') 
    return(ttv_count(samples, p_analytical, discard_burn_in))
}

ESS <- function(gen_quant, batch_size, burn_in, n_algorithms,
                sd_gen_quant, gq_index = 1) {
  # Do multiple estimates of the generated quantity, by splitting
  # the sample into multiple batches. This allows us to estimate
  # the Monte Carlo error.
  n_sample <- dim(gen_quant)[1]
  batch <- seq(from = burn_in + 1, to = n_sample, by = batch_size)
  batch_means <- matrix(NA, nrow = length(batch), ncol = n_algorithms)
  for (b in 1:length(batch)) {
    for (i in 1:n_algorithms) {
      batch_means[b, i] =
        mean(gen_quant[batch[b]:(batch[b] + batch_size - 1), i, gq_index])
    }
  }
  mc_sd <- rep(NA, n_algorithms)
  ratio <- rep(NA, n_algorithms)
  for (i in 1:n_algorithms) {
    mc_sd[i] <- sd(batch_means[, i])
    ratio[i] <- (sd_gen_quant / mc_sd[i])^2 / 1e3
  }

  list(mc_sd = mc_sd, ratio = ratio)
}

# Compute means for sufficient statistics for an array
# of samples (whihch a n x m x p object)
monte_carlo_estimate <- function(sample_all) {
  n_replica <- dim(sample_all)[1]
  n_iter <- dim(sample_all)[2]

  M_square <- rep(NA, n_replica)
  for (i in 1:n_replica) M_square[i] <- mean(rowSums(sample_all[i, , ])^2)
  
  mean_log_kernel <- rep(NA, n_replica)
  log_kernel <- array(NA, c(n_replica, n_iter))
  for (i in 1:n_replica) {
    for (j in 1:n_iter) {
      log_kernel[i, j] <- ising_kernel(beta_vector[i], A, 
                                       sample_all[i, j, ],
                                       log = T)
    }
    mean_log_kernel[i] <- mean(log_kernel[i, ])
  }

  list(M_square = M_square, mean_log_kernel)
}

###############################################################################
## Function to compute a multivariate normal distribution, given
# a mean vector and cholesky decomposition of the covariance matrix.
mvrnorm_cholesky <- function (n = 1, mu, L) {
  # Samples from a multivariate normal distribution, given a mean vector
  # a (lower triangular) Cholesky decomposition of the covariance matrix.
  dimension <- length(mu)
  X = matrix(NA, nrow = n, ncol = dimension)
  for (i in 1:n) {
    z = rnorm(dimension)
    X[i, ] = L %*% z + mu 
  }
  X
}

###############################################################################
## Function to generate the adjancency matrix for special grids

grid_graph <- function (m, anti_corr = FALSE, ignore_degrees = FALSE) {
  # The nodes row-ordered.
  # For example:
  #  1 2 3
  #  4 5 6
  #  7 8 9
  # The nodes are connected using the "+m" and
  # "non-edge neighbor" rules.
  # 
  # m: the number of column in the graph, and in the square grid setting
  # the square root of the total number of particles.
  # anti_corr: if true, then the connection is -1 instead of 1 with
  # probability 0.5.

  # FIX ME: condition sometimes fail -- numerical error?
  # if (m %% 1 != 0) stop("The number of particles must be a square!")

  n_particle = m^2
  A = matrix(0, nrow = n_particle, ncol = n_particle)
  for (i in 1:n_particle) {
    # edge neihbor rule
    if (i %% m != 0) {
      if (anti_corr) { connection = 1 - 2 * rbinom(1, 1, 0.5)
      } else { connection = 1}
      A[i + 1, i] = connection
      A[i, i + 1] = connection
    }

    # "+m" rule
    if (i + m <= n_particle) {
      if (anti_corr) { connection = 1 - 2 * rbinom(1, 1, 0.5)
      } else { connection = 1}
      A[i + m, i] = connection
      A[i, i + m] = connection
    }
  }

  degree_of_connect = mean(colSums(abs(A)))
  if (ignore_degrees) { d_divide = 1 } else 
    { d_divide = degree_of_connect}
  list(A = A / d_divide, d = degree_of_connect)
}

cube_graph <- function (m, anti_corr = FALSE) {
  # The nodes are row-ordered, then column-ordered.
  # The nodes are connected by first creating the block-diagonal
  # adjacency matrix with grid connection, and then connecting
  # each to the node "in front" of it.
  #
  # m is the lattice number, hence the total number of particles
  # is m x m x m.
  
  # FIX ME: condition sometimes fail -- numerical error?
  # if (m %% 1 != 0) stop("The number of particles must be a cube!")

  n_particle = m^3
  A = matrix(0, nrow = n_particle, ncol = n_particle)
  for (i in 1:m) {
    grid = grid_graph(m, anti_corr, ignore_degrees = TRUE)
    cube_index = ((i - 1) * m^2 + 1):(i * m^2)
    A[cube_index, cube_index] = grid$A
  }
  
  for (i in 1:(n_particle - m^2)) {
    j = i + m^2
    if (anti_corr) { connection = 1 - 2 * rbinom(1, 1, 0.5)
    } else { connection = 1}
    A[i, j] = connection
    A[j, i] = connection
  }

  degree_of_connect = mean(colSums(abs(A)))
  list(A = A / degree_of_connect, d = degree_of_connect)  
}

cube_graph(m = round(216^{1/3}), anti_corr = TRUE)

# Test the method
# m <- 3
# alpha <- 2
# A <- grid_graph(m)
# C <- A + alpha * diag(m^2)
# solve(C)


adjacency_graph <- function(n_part, type, anti_corr = FALSE,
                            do_round = TRUE, m_hopfield = 1) {
# n_part: number of particles
# type: complete, grid, gaussian
# anti_corr: for complete and grid, make some of the connections
#            negative.
# do_round: round the cube or square of n_part to make sure we
#           we get an integer. Useful for numerical stability.
  if (type == "complete") {
    A <- matrix(1, nrow = n_part, ncol = n_part)
    diag(A) <- 0
    A <- A / (n_part - 1)  # CHECK -- normalize?
  } else if (type == "grid") {
    lattice <- sqrt(n_part)
    if (do_round) lattice <- round(lattice)
    A <- grid_graph(lattice, anti_corr = anti_corr)$A
  } else if (type == "cube") {
    lattice <- n_part^{1/3}
    if (do_round) lattice <- round(lattice)
    A <- cube_graph(lattice, anti_corr = anti_corr)$A
  } else if (type == "gaussian") {
    A <- matrix(NA, nrow = n_part, ncol = n_part)
    for (i in 1:n_part) {
      for (j in 1:i) {
        A[i, j] = rnorm(1, 0, 1)
        if (i != j) A[j, i] = A[i, j]
      }
    }
    diag(A) <- 0
    A <- A / sqrt(n_part)
  } else if (type == "hopfield") {
    eta <- matrix(sample(c(-1, 1), n_part * m_hopfield, replace = TRUE),
                 nrow = n_part, ncol = m_hopfield)
    # CHECK: how do we normalize this?
    # A <- eta %*% t(eta) / m_hopfield
    A <- eta %*% t(eta) / max(n_part, m_hopfield)
    diag(A) <- 0
  } else {
    print("Could not recognize the type of graph.")
  }
  A
}


###############################################################################
## Softmax function
# from https://rpubs.com/FJRubio/softmax

softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

softmax_brute <- function(par) {
  
}

###############################################################################
## custom which function

which_one <- function (x) {
  i = 1
  while (x[i] != 1) i = i + 1
  i
}







