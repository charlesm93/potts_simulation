
rm(list = ls())
gc()
set.seed(1954)

# Adjust to your settings.
.libPaths("~/Rlib/")
setwd("~/Code/ising_model/scripts/")

library(ggplot2)
library(posterior)
source("tools.r")
source("algorithms.r")
source("algorithms_potts.r")

## Tuning parameters for experiments
# For Potts model with q = 4.
n_states <- 4
# beta_range <- c(0.25, 0.45, 0.5493062, 0.65, 0.85)  # grid graph
# beta_range <- c(1.35, 1.55, 1.647918, 1.75, 1.95)  # complete graph (q = 4)
# beta_range <- c(0.5, 0.8, 1, 1.2, 2)  # gaussian graph
# beta_range <- c(0.25, 0.44, 0.63)
# beta_range <- c(0.5, 1, 1.5) # complete graph for q = 2
beta_range <- rep(1, 14)  # hopfield model

n_beta <- length(beta_range)
# Options for graph types: "complete", "grid", "gaussian", "hopfield"
# m_hopfield parameter is only relevant for Hopfield model
m_hopfield <- c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 256)
  # c(1) # c(150, 200) # c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
graph_types <- rep("hopfield", n_beta)
anti_corr <- FALSE
graph_sizes <- rep(256, n_beta)
anti_corr <- rep(FALSE, n_beta)

n_iter <- 1e4
n_chains <- 4  # need multiple chains for diagnostics
total_sample <- n_iter * n_chains

# For Gaussian graphs, Wolff and SW cannot be used.
# algo_names <- c("AG_potts", "HB", "Wolff", "SW", "AG_potts_lowrank")  # complete
# algo_names <- c("AG_potts", "HB")  # Gaussian
algo_names <- c("AG_potts", "HB", "AG_potts_lowrank")
# algo_names <- c("AG_potts", "HB", "Wolff", "SW", "AG_potts_lowrank")
n_algorithms <- length(algo_names)

# Number of generated quantities, based on samples
# (for now only do one)
n_gen_quant <- 2
n_gen_sample <- n_iter / 2

burn_in <- 1e3

length_graph_size <- length(graph_sizes)
if (FALSE) {
  beta_range[1] <- 1.60
  length_graph_size <- 1  # only run one loop.
}

for (k in 1:length_graph_size) {
  graph <- graph_types[k]
  n_part <- graph_sizes[k]
  beta <- beta_range[k]
  print(paste0("graph size: ", n_part, ", type: ", graph_types[k],
               ", beta: ", beta_range[k], ", anti_corr: ", anti_corr[k]))
  
  # Build adjacency graph
  A <- adjacency_graph(n_part, type = graph, anti_corr = anti_corr[k],
                       m_hopfield = m_hopfield[k])
  
  init <- matrix(NA, nrow = n_chains, ncol = n_part)
  for (c in 1:n_chains) init[c, ] <- sample(1:n_states, n_part, replace = T)
  
  samples <- array(NA, c(n_iter, n_part, n_algorithms, n_chains))
  time <- rep(NA, n_algorithms)
  
  print("Doing AG sampling")
  time[1] <- system.time(
    for (c in 1:n_chains) {
      samples[, , 1, c] <- ag_potts_simple(A = A, beta = beta, init = init[c, ],
                                           n_iter = n_iter, n_states = n_states)
    })[3]

  print("Doing HB sampling")
  time[2] <- system.time(
    for (c in 1:n_chains) {
      samples[, , 2, c] <- hb_sampler(A = A, beta = beta, init = init[c, ],
                                      n_iter = n_iter, is_potts = TRUE, 
                                      n_states = n_states,
                                      sub_sample = ceiling(n_part / 10))
    })[3]

  if (graph != "gaussian" & graph != "hopfield") {
    print("Doing Wolff sampling")
    time[3] <- system.time(
      for (c in 1:n_chains) {
        samples[, , 3, c] <- wolff_sampler(A = A, beta = beta, init = init[c, ],
                                           n_iter = n_iter, is_potts = TRUE,
                                           n_states = n_states)
      })[3]

    print("Doing SW sampling")
    time[4] <- system.time(
      for (c in 1:n_chains) {
        samples[, , 4, c] <- sw_sampler(A = A, beta = beta, init = init[c, ],
                                        n_iter = n_iter, is_potts = TRUE,
                                        n_states = n_states)
      })[3]
  }

  if (graph == "complete" || graph == "hopfield") {
    print("Doing AG lowrank sampling")
    if (graph == "complete") index = 5
    if (graph == "hopfield") index = 3
    
    # REMARK: default epsilon = 1e-12 doesn't remove approximatively 0 eigenvalues.
    time[index] <- system.time(
      for (c in 1:n_chains) {
        samples[, , index, c] <- ag_potts_lowrank(A = A, beta = beta, init = init[c, ],
                                              n_states = n_states, n_iter = n_iter,
                                              alpha_modif = 0,
                                              epsilon = 1e-5)
      })[3]
  }

  all_samples <- array(NA, dim = c(total_sample, n_part, n_algorithms))
  for (c in 1:n_chains) {
    iteration_index <- ((c - 1) * n_iter + 1):(c * n_iter)
    all_samples[iteration_index, , ] <- samples[, , , c]
  }
  
  gen_quant <- array(NA, c(n_gen_sample, n_algorithms, n_chains, n_gen_quant))
  index_start <- n_iter - n_gen_sample + 1
  print("Computing generated quantities.")
  for (j in 1:n_algorithms) {
    for (c in 1:n_chains) {
      for (i in 1:n_gen_sample) {
        gen_quant[i, j, c, 1] <- potts_log_kernel(beta, A, 
                                  samples[index_start - 1 + i, , j, c])
        gen_quant[i, j, c, 2] <- gen_quant[i, j, c, 1]
      }
    }
  }

  # Compute the effective samples size and results
  mean <- array(NA, dim = c(n_gen_quant, n_algorithms))
  mcse <- array(NA, dim = c(n_gen_quant, n_algorithms))
  rhat <- array(NA, dim = c(n_gen_quant, n_algorithms))
  ess <- array(NA, dim = c(n_gen_quant, n_algorithms))
  eff <- array(NA, dim = c(n_gen_quant, n_algorithms))

  for (j in 1:n_algorithms) {
    draws_formatted <- as_draws_df(gen_quant[, j, , ])
    names(draws_formatted)[1] <- c("log_kernel")
    summary <- summarize_draws(draws_formatted,
                               measure = c("mean", "mcse_mean", "ess_bulk", "rhat"))
    names(summary) <- c("variable", "mean", "mcse_mean", "ess_bulk", "rhat")
    mean[, j] <- summary$mean
    mcse[, j] <- summary$mcse_mean
    rhat[, j] <- summary$rhat
    ess[, j] <- summary$ess_bulk   # / total_sample * 1e3
    eff[, j] <- ess[, j] / time[j]   # (time[j] * 1e3 / total_sample)
    # Note: returns ESS / time, and corrects for the fact ess is the number
    # of effectively independent samples per 1000 iterations.
  }

  eff_data <- data.frame(cbind(mean, mcse, rhat, ess, eff))
  names(eff_data) <- c(paste0("mean_", algo_names), paste0("mcse_", algo_names),
                       paste0("rhat_", algo_names),
                       paste0("ess_", algo_names), paste0("eff_", algo_names))

  output_name <- paste0("deliv/measurements/potts/", n_part, graph,
                        "_", beta, "_", n_states, "_", m_hopfield[k], "_efficiency")
  write.csv(eff_data, paste0(output_name, ".csv"))
}
