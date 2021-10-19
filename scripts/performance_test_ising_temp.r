
# Script to test the performance of tempering algorithms.

rm(list = ls())
gc()
set.seed(1955)  # change seed to check if behavior is consistent

# Adjust to your settings.
.libPaths("~/Rlib/")
setwd("~/Code/ising_model/scripts/")

library(ggplot2)
library(tidyverse)
library(posterior)
source("tools.r")
source("algorithms.r")
source("algorithm_tempering.r")

## Tuning parameters for the problem
beta_vector <- c(0.25, 0.44, 0.63)
graph_type <- "grid"
# options: "complete", "grid", "cube", and "gaussian"
anti_corr <- FALSE  # allow bonds to be negative for grid and cubic graphs
B <- 0  # magnetic field
graph_size <- 16
A <- adjacency_graph(graph_size, type = graph_type, anti_corr = anti_corr)

## Tuning parameters for the algorithm
n_mc_iter <- 1e4  # 1e3
n_exchange <- 1  # 40  # pick an exchange such that r_hat goes to 1.0
n_burnin <- 1e3
total_sample <- n_mc_iter * n_exchange
n_chains <- 4
algo_names <- c("temp_AG", "AG", "Wolff")
# algo_names <- c("temp_AG", "temp_HB", "temp_hybrid", "HB")

# Settings for hybrid method (not used much)
f1 <- function(A, beta, B, init, n_iter) {
  sumit_sampler(A, beta, B, init, n_iter, perturb = TRUE)
}
sub_adjust <- 10  # 10
f2 <- function(A, beta, B, init, n_iter) {
  hb_sampler(A, beta, B, init, n_iter, 
             sub_sample = ceiling(nrow(A)/ sub_adjust))
}
n_temp <- length(beta_vector)
n_algorithms <- length(algo_names)
which_sampler <- array(c(rep(1, n_temp),
                         rep(2, n_temp),
                         c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
                         rep(NA, n_temp)), dim = c(n_temp, n_algorithms))

## Tuning parameters for generated quantities
n_gen_quant <- 3
n_gen_sample <- total_sample - n_burnin

##########################################################################
## Run experiment
print(paste0("graph_size: ", graph_size, ", B-field: ", B,
             ", type: ", graph_type, ", beta: varying",
             ", anti_corr: ", anti_corr))
print("beta_vector: ")
print(beta_vector)

init <- array(NA, dim = c(n_chains, n_temp, graph_size))
for (c in 1:n_chains) {
  for (b in 1:n_temp) init[c, b, ] <- sample(c(-1, 1), graph_size, replace = T)
}

## Draw samples
samples <- 
  array(NA, c(n_algorithms, n_chains, n_temp, total_sample, graph_size))
time <- rep(NA, n_algorithms)

print("Doing temp_AG sampling")
time[1] <- system.time (
  for (c in 1:n_chains) {
    samples[1, c, , , ] <- tempering(A, beta_vector, init[c, , ], n_exchange,
                                   n_mc_iter, B, f1, f2, which_sampler[, 1])
  }
)[3]

# print("Doing temp_HB sampling")
# time[2] <- system.time (
#   for (c in 1:n_chains) {
#     samples[2, c, , , ] <- tempering(A, beta_vector, init[c, , ], n_exchange,
#                                    n_mc_iter, B, f1, f2, which_sampler[, 2])
#   }
# )[3]
# 
# print("Doing temp_hybrid sampling")
# time[3] <- system.time (
#   for (c in 1:n_chains) {
#     samples[3, c, , , ] <- tempering(A, beta_vector, init[c, , ], n_exchange,
#                                    n_mc_iter, B, f1, f2, 
#                                    which_sampler[, 3])
#   }
# )[3]

print("Doing AG sampling")
time[2] <- system.time (
  for (c in 1:n_chains) {
    for (b in 1:n_temp) {
      samples[2, c, b, , ] <- 
        sumit_sampler(A, beta_vector[b], B, init[c, b, ], total_sample)
    }
  }
)[3]

print("Doing Wolff sampling")
time[3] <- system.time (
  for (c in 1:n_chains) {
    for (b in 1:n_temp) {
      samples[3, c, b, , ] <- 
        wolff_sampler(A, beta_vector[b], B, init[c, b, ], total_sample)
    }
  }
)[3]


# Compute generated quantities
# NOTE: this part takes quite a long time...
print("Computing generated quantities")
n_gen_quant <- 3
gen_quant <-
  array(NA, c(n_gen_quant, n_algorithms, n_chains, n_temp,
              total_sample - n_burnin))
for (j in 1:n_algorithms) {
  print(paste0("Chain: ", 0, " / ", n_chains))
  for (c in 1:n_chains) {
    for (b in 1:n_temp) {
      gen_quant[1, j, c, b, ] <- rowSums(samples[j, c, b, -(1:n_burnin), ])
      gen_quant[2, j, c, b, ] <- gen_quant[1, j, c, b, ]^2
      
      for (i in (n_burnin + 1):total_sample) {
        gen_quant[3, j, c, b, i - n_burnin] <-
          ising_kernel(beta_vector[b], A, samples[j, c, b, i, ],
                       log = T)
      }
    }
    print(paste0("Chain: ", c, " / ", n_chains))
  }
}

# Compute ESS, eff, and other diagnostics using the package Posterior.
dim_summary <- c(n_gen_quant, n_algorithms, n_temp)
sample_mean <- array(NA, dim = dim_summary)
sample_mcse <- array(NA, dim = dim_summary)
sample_r_hat <- array(NA, dim = dim_summary)
sample_ess <- array(NA, dim = dim_summary)
sample_eff <- array(NA, dim = dim_summary)

for (j in 1:n_algorithms) {
  for (b in 1:n_temp) {
    perm_gen_quant <- aperm(gen_quant[, j, , b, ], perm = c(3, 2, 1))
    draws_formatted <- as_draws_df(perm_gen_quant)
    names(draws_formatted)[1:3] <- c("spin", "spin_squared", "corr")
    
    summary_sample <- summarize_draws(draws_formatted,
      measure = c("mean", "mcse_mean", "rhat", "ess_bulk"))
    names(summary_sample) <- 
      c("variable", "mean", "mcse_mean", "rhat", "ess_bulk")
    
    sample_mean[, j, b] <- summary_sample$mean
    sample_mcse[, j, b] <- summary_sample$mcse_mean
    sample_r_hat[, j, b] <- summary_sample$rhat
    sample_ess[, j, b] <- summary_sample$ess_bulk
    sample_eff[, j, b] <- sample_ess[, j, b] / time[j]
  }
}

# create an output file for each replica
for (b in 1:n_temp) {
  saved_summary <- data.frame(cbind(sample_mean[, , b],
                                    sample_mcse[, , b],
                                    sample_r_hat[, , b],
                                    sample_ess[, , b],
                                    sample_eff[, , b]))
  names(saved_summary) <-
    c(paste0("mean_", algo_names),
      paste0("mcse_", algo_names),
      paste0("rhat_", algo_names),
      paste0("ess_", algo_names),
      paste0("eff_", algo_names))

  print(paste0("beta = ", beta_vector[b]))
  print(paste0(saved_summary$mean_temp_AG, " ", saved_summary$mean_temp_HB))
  
  saved_summary
  # Can see from the summary that the heat bath algorithm is strongly biased.

  output_name <- paste0("deliv/measurements_temp/", graph_size, graph_type, 
                        "_", B, "_", beta_vector[b], "_summary")
  write.csv(saved_summary, paste0(output_name, ".csv"))
}


# Examine the efficiency over betas.
gen_quant_index <- 2  # 1, 2, 3

# Don't include heat bath because the estimator is biased.
plot_eff <- as.vector(sample_eff[gen_quant_index, 1:3, ]) # %*%
                        # diag(1 /sample_eff[gen_quant_index, 1, ]))
plot_algo <- rep(algo_names, n_temp)
plot_beta <- rep(beta_vector, each = n_algorithms)

plot <- ggplot(data = data.frame(eff = plot_eff,
                                 algo = plot_algo,
                                 beta = plot_beta),
               aes(x = beta, y = eff, color = algo)) +
  theme_bw() + geom_point() + geom_line() +
  scale_y_continuous(trans = 'log10')
plot

# Examine the distributions of the generated quantities
# gen_data <- as.vector(gen_quant[1, , , , ])
# 
# algorithm_data <- rep(algo_names, n_chains * n_temp * n_gen_sample)
# beta_data <- rep(rep(beta_vector, each = n_algorithms * n_chains),
#                  n_gen_sample)

dens_data <- data.frame(m_squared = as.vector(gen_quant[2, , , , ]),
                    corr = as.vector(gen_quant[3, , , , ]),
                    algo = rep(algo_names, n_chains * n_temp * n_gen_sample),
                    chain = as.factor(rep(rep(1:n_chains, each = n_algorithms),
                                n_temp * n_gen_sample)),
                    beta = rep(rep(beta_vector, 
                                   each = n_algorithms * n_chains),
                               n_gen_sample)
                    )

## Compare posterior samples generated by both algorithms
dens_plot <- 
  ggplot(data = dens_data[dens_data$algo %in% c("temp_AG", "AG", "Wolff"), ],
        aes(x = m_squared, fill = algo)) + theme_bw() +
  geom_histogram(color = "black",alpha = 0.5,
                 bins = 30, position = "identity") +
  facet_wrap(~beta, scale = "free")
dens_plot

dens_plot2 <-
  ggplot(data = dens_data[dens_data$algo %in% c("temp_AG", "AG", "Wolff"), ],
         aes(x = corr, fill = algo)) + theme_bw() +
  geom_histogram(color = "black",alpha = 0.5,
                 bins = 30, position = "identity") +
  facet_wrap(~beta, scale = "free")
dens_plot2


## Compare samples between chains
dens_plot3 <- 
  ggplot(data = dens_data[dens_data$algo == "temp_AG", ],
         aes(x = corr, fill = chain)) + theme_bw() +
  geom_histogram(color = "black",alpha = 0.5,
                 bins = 30, position = "identity") +
  facet_wrap(~beta, scale = "free")
dens_plot3

