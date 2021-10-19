
# Test Ising model simulation algorithms under various
# configurations. This file can be run on the cluster
# using the run_R_script.sh file.
#
# This is version 2. This code makes it easier to
# add or remove algorithms being tested.
# Version 1 has some experimental code that is still
# of interest. Some comments:
#
# ttv_type: type of total variation being computed.
#           Over large complete graphs, it is best to use
#           variation over count space to insure scalability.
# graph_time: measures and plot ttv against (average estimated)
#             time. This is limited by the  time taken by the 
#             shortest algorithm to run.
# iterations: set of points where we calculate the ttv. If we
#             produce large samples, it is best to take big
#             steps.

rm(list = ls())
gc()
set.seed(1954)

# Adjust to your settings.
.libPaths("~/Rlib/")
setwd("~/Code/ising_model/scripts/")
# setwd("/rigel/home/ccm2172/ising_model/scripts")  # for cluster
# .libPaths("/rigel/home/ccm2172/Rlib")

library(ggplot2)
library(tidyverse)
library(posterior)
source("tools.r")
source("algorithms.r")
source("algorithms_potts.r")


## Tuning parameters for experiments.

# NOTE: for grid graphs, beta gets adjusted after normalization
# beta_range <- c(log(1 + sqrt(2)) / 2, log(1 + sqrt(2)) / 2)
# beta_range <- c(0.5, 1 / sqrt(2), 0.9) # c(3, 5, 6)  # Grid graph
beta_range <- c(0.5, 1, 2)  # c(0.5, 1, 1.5) # Complete graph
n_beta <- length(beta_range)

graph_types <- rep("complete", n_beta)  
# options: "grid", "complete", "cube", and "gaussian"
anti_corr <- FALSE  # allows for negative bonds.

B_range <- rep(0, n_beta)  # magnetic field
graph_sizes <- rep(256, n_beta)  # 16^2

n_sample <- 1e4 # 1e5
n_chains <- 4  # need multiple chains for diagnostics
total_sample <- n_sample * n_chains

graph_ttv <- FALSE
ttv_type <- "total"  # options: "total" or "count"
iterations <- c(10, seq(from = 500, to = n_sample * n_chains, by = 500))
n_strat = length(iterations)

graph_time <- FALSE
time_strat <- 40

# algo_names <- c("Sumit", "Sumit_p", "MH", "MH_p", "SW", "Wolff") #, "Exact")
algo_names <- c("AG_lowrank")
# algo_names <- c("Sumit", "HB", "Wolff", "AG_lowrank")
# algo_names <- c("Sumit", "Sumit_p", "MH", "MH_p", "HB")
# algo_names <- c("Sumit", "HB", "AG_potts")
n_algorithms <- length(algo_names)
do_exact_sample <- FALSE

# Number of generated quantities, based on samples
# Do total spin and number of identical neighbors.
n_gen_quant <- 3
n_gen_sample <- n_sample / 2
n_plot_gen_sample <- 2e2
index_gen_plot <- (n_gen_sample - n_plot_gen_sample + 1):n_gen_sample

batch_size <- 1e3
burn_in <- 1e3

for (k in 1:length(graph_sizes)) {
  graph <- graph_types[k]
  n_part <- graph_sizes[k]
  B <- B_range[k]
  beta <- beta_range[k]
  print(paste0("graph size: ", n_part, ", B-field: ", B, ", type: ", graph,
               ", beta: ", beta, ", anti_corr: ", anti_corr))

  # Build adjacency graph
  A <- adjacency_graph(n_part, type = graph, anti_corr = anti_corr)
  
  # if (graph == "complete") {
  #   A <- matrix(1, nrow = n_part, ncol = n_part)
  #   diag(A) <- 0
  #   A <- A / (n_part - 1)  # CHECK -- normalize?
  # } else {  # graph == "grid"
  #   A_prop <- grid_graph(sqrt(n_part))
  #   A <- A_prop$A
  #   beta <- beta *  A_prop$d
  # }

  init <- matrix(NA, nrow = n_chains, ncol = n_part)
  for (c in 1:n_chains) init[c, ] <- sample(c(-1, 1), n_part, replace = T)

  ## Draw samples
  samples <- array(NA, c(n_sample, n_part, n_algorithms, n_chains))
  time <- rep(NA, n_algorithms)

  print("Doing Sumit sampling")
  print(paste0("Chain: ", 0, " / ", n_chains))
  time[1] <- system.time(
    for (c in 1:n_chains) {
      samples[, , 1, c] <- sumit_sampler(A, beta, B, init[c, ],
                                      n_iter = n_sample,
                                      perturb = F)
      print(paste0("Chain: ", c, " / ", n_chains))
      })[3]
  
  # print("Doing Sumit sampling with perturbation")
  # time[2] <- system.time(
  #   for (c in 1:n_chains) {
  #     samples[, , 2, c] <- sumit_sampler(A, beta, B, init[c, ],
  #                                     n_iter = n_sample,
  #                                     perturb = T)
  #   })[3]

  # print("Doing MH sampling.")
  # time[3] <- system.time(
  #   for (c in 1:n_chains) {
  #   samples[, , 3, c] <- mh_sampler(A, beta, B, init[c, ],
  #                                n_iter = n_sample,
  #                                perturb = F)
  #   })[3]

  # print("Doing MH sampling with perturbation.")
  # time[4] <- system.time(
  #   for (c in 1:n_chains) {
  #     samples[, , 4, c] <- mh_sampler(A, beta, B, init[c, ], 
  #                                     n_sample, perturb = T)
  #     })[3]

  print("Doing heatbath sampling.")
  print(paste0("Chain: ", 0, " / ", n_chains))
  time[2] <- system.time(
    for (c in 1:n_chains) {
      samples[, , 2, c] <- hb_sampler(A, beta, B, init[c, ],
                                      n_sample,
                                      sub_sample = ceiling(n_part / 10))
      print(paste0("Chain: ", c, " / ", n_chains))
    })[3]

  # print("Doing potts sampling")
  # init_potts <- matrix(NA, nrow = n_chains, ncol = n_part)
  # for (c in 1:n_chains) init_potts[c, ] <- init[c, ]
  # init_potts[init_potts == -1] = 2
  # 
  # time[3] <- system.time(
  #   for (c in 1:n_chains) {
  #     samples[, , 3, c] <- ag_potts(A = A, beta = 2 * beta, n_states = 2,
  #                                   init = init_potts[c, ], n_iter = n_sample
  #                                   )
  #  })[3]
  # 
  #  samples[samples == 2] <- -1
  
  # print("Doing SW sampling.")
  # time[5] <- system.time(
  #   for (c in 1:n_chains) {
  #     samples[, , 5, c] <- sw_sampler(A, beta, B, init[c, ], n_sample)
  #   })[3]

  print("Doing Wolff sampling.")
  print(paste0("Chain: ", 0, " / ", n_chains))
  time[3] <- system.time(
    for (c in 1:n_chains) {
      samples[, , 3, c] <- wolff_sampler(A, beta, B, init[c, ], n_sample)
      print(paste0("Chain: ", c, " / ", n_chains))
    })[3]
  
  print("Doing AG lowrank sampling.")
  print(paste0("Chain: ", 0, " / ", n_chains))
  time[4] <- system.time(
    for (c in 1:n_chains) {
      samples[, , 4, c] <- ag_ising_lowrank(A = A, beta = beta,
                                            init = init[c, ],
                                            n_iter = n_sample,
                                            epsilon = 1e-12)
      print(paste0("Chain: ", c, " / ", n_chains))
    })[3]

  # print("Doing AG lowrank sampling.")
  # print(paste0("Chain: ", 0, " / ", n_chains))
  # time[1] <- system.time(
  #   for (c in 1:n_chains) {
  #     samples[, , 1, c] <- ag_ising_lowrank(A = A, beta = beta,
  #                                           init = init[c, ],
  #                                           n_iter = n_sample,
  #                                           epsilon = 1e-12)
  #     print(paste0("Chain: ", c, " / ", n_chains))
  #   })[3]
  
  if (do_exact_sample) {
    exact_solution <- partition(beta, A, B, prob = T)
    prob_analytical <- exact_solution$p
    time[7] <- system.time(
      for (c in 1:n_chains) {
        samples[, , 7, c] <- exact_sampler(prob_analytical,
                                           n_part, n_sample)
      })[3]
  }

  # TO DO -- find a way to not duplicate the data set. Are there
  # such things as C++ pointers in R?
  all_samples <- array(NA, dim = c(total_sample, n_part, n_algorithms))
  for (c in 1:n_chains) {
    iteration_index <- ((c - 1) * n_sample + 1):(c * n_sample)
    all_samples[iteration_index, , ] <- samples[, , , c]
  }

  need_exact <- ((!do_exact_sample & (ttv_type == "total" | graph == "grid")) &
                   graph_ttv)
  if (need_exact) prob_analytical <- partition(beta, A, prob = T)$p

  ## NOTE:
  # Variational distance over the count space is computed 
  # efficiently for the complete graph. 
  # For the grid graph, the exact solution is first computed
  # and the count sapce deduced from it (not faster than doing ttv)
  if (graph_ttv) {
    if (ttv_type == "count") {
      if (graph == "complete") prob_analytical <-
          partition_complete(beta, A, B, prob = T)$p

      if (graph == "grid") prob_analytical <-
          count_probability(prob_analytical, n_particles)
    }
  
    tv <- array(NA, c(n_strat, n_algorithms))
    print("Computing variational distance.")
    for (i in 1:n_strat) {
      if (i %% 5 == 0) print(paste0(i, " / ", n_strat))
      for (j in 1:n_algorithms) {
        tv[i, j] = ttv_comp(all_samples[1:iterations[i], , j],
                            prob_analytical, ttv_type,
                            discard_burn_in = T)
      }
    }
  }
  
  gen_quant <- array(NA, c(n_gen_sample, n_algorithms, n_chains, n_gen_quant))
  index_start <- n_sample - n_gen_sample + 1
  print("Computing generated quantities.")
  for (j in 1:n_algorithms) {
    for (c in 1:n_chains) {
      gen_quant[, j, c, 1] <- rowSums(samples[index_start:n_sample, , j, c])
      gen_quant[, j, c, 2] <- rowSums(samples[index_start:n_sample, , j, c])^2
      for (i in 1:n_gen_sample) {
        gen_quant[i, j, c, 3] <-
          ising_kernel(beta, A, samples[index_start - 1 + i, , j, c], log = T)
      }
    }
  }
  
  ## Summary plot
  # 1. total variation vs iterations
  # 2. total variation vs (varying) time
  # 3. total variation vs fixed time.
  if (graph_ttv) {
    tv_data <- data.frame(iterations, tv)
    tv_time_tot <- array(NA, c(n_strat, n_algorithms))
    for (j in 1:n_algorithms) {
      tv_time_tot[, j] <- iterations * time[j] / n_sample
    }
    tv_data <- cbind(tv_data, tv_time_tot)
    tv_names <- c(paste0("tv_", algo_names), paste0("time_", algo_names))
    names(tv_data) <- c("iterations", tv_names)
  
    tv_data <- tv_data %>%
      gather(key = key, value = time,
             -c(iterations, paste0("tv_", algo_names)))
  
    tv_data <- tv_data %>%
      gather(key = key, value = tv,
             -c(iterations, time, key))
  
    plot <- ggplot(tv_data) + aes(iterations, tv, color = key) +
      geom_line() + scale_y_continuous(trans = 'log10') + theme_bw() +
      ylab("total variation distance")
    plot
  }
  
  if (graph_time) {
    print("Computing total variation distance accross strata.")
    tv_time <- array(NA, c(time_strat, n_algorithms))
    time_min <- min(time)
    time_iter <- seq(from = 0, to = time_min, by = time_min / time_strat)
    n_produced <- time_min / time * (n_sample * n_chains)
    
    for (j in 1:n_algorithms) {
      step <- round(n_produced[j] / time_strat)
      strat <- c(seq(from = step, to = round(n_produced[j]), by = step),
                 round(n_produced[j]))
      if (length(strat) < time_strat) strat <- c(strat, round(n_produced[j]))
      for (i in 1:time_strat) {
        tv_time[i, j] = ttv_comp(all_samples[1:strat[i], , j],
                                 prob_analytical, type = ttv_type,
                                 discard_burn_in = T)
      }
    }
  }
  # TO DO: do full plot against time, instead of iterations.
  # Can't do this using the above method.
  
  ## Do tv vs time plots
  if (graph_time) {
    tv_time_data <- data.frame(time_iter[-1], tv_time)
    names(tv_time_data) <- c("time", paste0("tv_", algo_names))
    tv_time_data <- tv_time_data %>%
      gather(key = key, value = tv, -time)
    plot2 <- ggplot(tv_time_data) + aes(time, tv, color = key) +
      geom_line() + scale_y_continuous(trans = 'log10') + theme_bw()
    plot2
  }

  ## plot markov chains for chain 1.
  gen_data <-
    data.frame(gen_quant[index_gen_plot, , 1, 1],
               gen_quant[index_gen_plot, , 2, 1],
               gen_quant[index_gen_plot, , 3, 1],
               1:n_plot_gen_sample)
  names(gen_data) <- c(paste0("count_", algo_names),
                       paste0("sqr_count_", algo_names),
                       paste0("corr_", algo_names), "iteration")
  gen_data <- gen_data %>%
    gather(key = key, value = count,
           -c(iteration, paste0("sqr_count_", algo_names),
                         paste0("corr_", algo_names))) %>%
    gather(key = key2, value = sqr_count,
           -c(iteration, count, key,
              paste0("corr_", algo_names))) %>%
    gather(key = key3, value = corr,
           -c(iteration, count, sqr_count, key, key2))
  
  # TO DO -- find slicker way to do this
  # gen_count_data$corr <- as.numeric(gen_count_data$corr)

  plot3 <- ggplot(gen_data) +
    aes(iteration, count, color = key) +
    geom_line(size = 0.2) + facet_wrap(~key) + theme_bw()
  plot3
  
  # gen_count_data <- gen_count_data %>%
  #   gather(key = key2, value = corr,
  #          -c(iteration, key, count))
  # gen_corr_data <- data.frame(gen_count_data,
  #                             rep(1:n_gen_sample, n_algorithms))
  # names(gen_corr_data) <- c(paste0("corr_", algo_names), "iteration")
  # gen_corr_data <- gen_corr_data %>%
  #   gather(key = key, value = corr, -iteration)
  
  plot4 <- ggplot(gen_data) +
    aes(iteration, sqr_count, color = key2) +
    geom_line(size = 0.2) + facet_wrap(~key2) + theme_bw()
  plot4
  
  plot5 <- ggplot(gen_data) +
    aes(iteration, corr, color = key3) +
    geom_line(size = 0.2) + facet_wrap(~key3) + theme_bw()
  plot5
  
  # Compute the effective samples size and results
  mean <- array(NA, dim = c(n_gen_quant, n_algorithms))
  mcse <- array(NA, dim = c(n_gen_quant, n_algorithms))
  ess <- array(NA, dim = c(n_gen_quant, n_algorithms))
  eff <- array(NA, dim = c(n_gen_quant, n_algorithms))
  for (j in 1:n_algorithms) {
    draws_formatted <- as_draws_df(gen_quant[, j, , ])
    names(draws_formatted)[1:3] <- c("spin", "spin_squared", "corr")
    
    summary <- summarize_draws(draws_formatted)
    summary_mcse <- summarize_draws(draws_formatted,
                                    measure = c("mcse_mean"))

    # summary <- summarize_draws(draws_formatted,
    #    measure = c("mean", "mcse_mean", "ess_bulk", "rhat"))
    # names(summary) <- c("variable", "mean", "mcse_mean", "ess_bulk", "rhat")

    mean[, j] <- summary$mean
    mcse[, j] <- summary_mcse$measure
    ess[, j] <- summary$ess_bulk   # / total_sample * 1e3
    eff[, j] <- ess[, j] / time[j]   # (time[j] * 1e3 / total_sample)
    # Note: returns ESS / time, and corrects for the fact ess is the number
    # of effectively independent samples per 1000 iterations.
  }

  ## Save results and figures.
  if (graph_ttv) {
    output_name <- paste0("deliv/measurements/", n_part, graph,
                          "_", ttv_type, "_tv.data")
    write.csv(tv_data, paste0(output_name, ".csv"))

    output_name <- paste0("deliv/measurements/", n_part, graph,
                          "_", ttv_type, "_time_tv.data")
    write.csv(tv_time_data, paste0(output_name, ".csv"))
  }

  # gen_corr_data$count <- gen_count_data$count
  output_name <- paste0("deliv/measurements/", n_part, graph, "_", B,
                        "_", beta, "_mcmc.data")
  write.csv(gen_data, paste0(output_name, ".csv"))
  
  output_name <- paste0("deliv/measurements/", n_part, graph, "_", B,
                         "_", beta,"_efficiency")
  
  eff_data <- data.frame(cbind(mean, mcse, ess, eff))
  names(eff_data) <- c(paste0("mean_", algo_names), paste0("mcse_", algo_names),
                       paste0("ess_", algo_names), paste0("eff_", algo_names))

  write.csv(eff_data, paste0(output_name, ".csv"))

  pdf(file = paste0("deliv/figures/", n_part, graph, "_", B,
                    "_performance%03d.pdf"),
      width = 6, height = 6, onefile = F)
  if (graph_ttv) print(plot)
  if (graph_time) print(plot2)
  print(plot3)
  print(plot4)
  print(plot5)
  dev.off()
}
