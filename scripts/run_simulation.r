
rm(list = ls())
gc()
seed <- 1954
set.seed(seed)

# Adjust to your settings.
.libPaths("~/Rlib/")
setwd("~/Code/ising_model/scripts/")

# library(ggplot2)
# library(tidyverse)
# library(posterior)
source("tools.r")
# source("algorithms.r")
source("algorithms_potts.r")
source("algorithms_hmc.r")
source("algorithm_tempering.r")

###############################################################################
## Utility functions

generate_quantities <- function(samples, beta, A, n_chains, q) {
  # Transform samples into quantities of interest:
  #  (i) sum of spins
  #  (ii) sum of squared spins
  #  (iii) Hamiltonian
  #
  # samples: input samples with shape (n_sample, n_particle, algorith, n_chain)
  # algo_index: which third index to use in samples
  # beta: inverse temperature
  # A: adjacency graph.
  n_gen_quant <- 3
  n_sample <- dim(samples)[1]
  gen_quant <- array(NA, dim = c(n_sample, n_chains, n_gen_quant))
  for (c in 1:n_chains) {
    gen_quant[, c, 1] <- rowSums(samples[, , c])
    gen_quant[, c, 2] <- rowSums(samples[, , c]^2)
    for (i in 1:n_sample) {
      if (q == 2) {
        gen_quant[i, c, 3] <- ising_kernel(beta, A, samples[i, , c],
                                           log = TRUE)
      } else {
        gen_quant[i, c, 3] <- potts_log_kernel(beta, A, samples[i, , c])
      }
    }
  }
  return(gen_quant)
}

save_output <- function(samples, time, algo_index, beta, A, n_chains,
                        output_name_root, algo_names,
                        perturbation_step=FALSE, q=2, m_hopfield=0) {
  # Compute generated quantities and save output to a file.
  gen_quant <- generate_quantities(samples, beta, A, n_chains, q)
  
  sample_output <- list(samples = samples,
                        time = time,
                        gen_quant = gen_quant)
  
  output_name <- paste0(output_name_root, algo_names[algo_index], "_perturb_",
                        perturbation_step)
  if (q != 2) {
    output_name <- paste0(output_name, "_q=", q)
  }
  
  if (m_hopfield != 0) {
    output_name <- paste0(output_name, "_m=", m_hopfield)
  }

  output_name <- paste0(output_name, ".RData")

  deliv_dir <- file.path(getwd(), "deliv_v3")
  saveRDS(sample_output, file = paste0(deliv_dir, "/", output_name))
}


###############################################################################
## Function to run experiments

run_simulations <- function(beta_range, graph_types, graph_sizes, n_samples_algo,
                            algo_names, burnin_fraction = 0.5, n_chains = 4, 
                            anti_corr = FALSE, perturbation_step = FALSE, q = 2,
                            m_hopfield_range = NULL,
                            beta_vector = seq(from = 0.5, to = 3, by = 0.25),
                            n_exchange = 40, n_mc_iter = 1e3) {
  # Runs and saves sampling algorithms for various Ising and Potts models.
  #
  # beta_range: vector of inverse temperatures. For grid graphs, beta gets
  #             adjusted after normalization.
  # graph_types: adjacency graph for model. Options: "grid", "complete", "cube"
  #              and "gaussian".
  # graph_sizes: vector of sizes for graph. If grid, needs to be a square.
  # n_samples_algo: number of iterations per chain for each algorithm.
  # algo_names: vector containing the name of the sampling algorithms. Options
  #             of algorithms to include:
  #             - "AG": Gibbs with an auxiliary Gaussian (AG)
  #             - "AG_lowrank": Gibbs with AG, using a low-rank approximation.
  #             - "AG_potts": Gibbs with AG for potts models. Here, it is used
  #                           for the case where there are two "colors", although
  #                           unlike AG, it is not specialized for this case.
  #             - "HB": Heat bath with single site update.
  #             - "HB_long": Heat bath with sequential update (of a subset)
  #             - "BHB": black-box heat bath, with non-specialized calculations
  #                      of the energy.
  #             - "SW": Swendsen-Wang algorithm.
  #             - "Wolff": Wolff algorithm.
  #             - "GWG": Gibbs with gradient-based proposal.
  #             - "dHMC": discrete HMC using AG.
  # burnin_fraction: fraction of n_sample discarded for burnin phase.
  # n_chains: number of chains. Needs to be greater than 2 for diagnostics.
  # anti_corr: if TRUE, allows negative edges in the graph.
  # perturbation_step: if TRUE, adds the flip step to all algorithms.
  # q: number of colors (if 2, reduces to Ising model)
  # m_hopfield_range: vector of ranks of adjacency matrix when graph_type == "hopfield".
  # beta_vec: for tempering algorithms, the temperatures at which the different replicas
  #           are run. NOTE: tempering algorithms do not use the beta_range argument.
  # n_exchange: for tempering algorithms, number of attempted exchanges.
  # n_mc_iter: number of MCMC iterations between exchanges for each tempered sampler.

  B_range = rep(0, length(graph_sizes))
  if (q >= 3) {
    is_potts = TRUE
  } else {
    is_potts = FALSE
  }

  for (k in 1:length(graph_sizes)) {
    graph <- graph_types[k]
    n_part <- graph_sizes[k]
    B <- B_range[k]
    beta <- beta_range[k]
    print(paste0("graph size: ", n_part, ", B-field: ", B, ", type: ", graph,
                 ", beta: ", beta, ", anti_corr: ", anti_corr))

    if (graph == "hopfield") {
      m_hopfield <- m_hopfield_range[k]
      A <- adjacency_graph(n_part, type = graph, anti_corr = anti_corr,
                           m_hopfield = m_hopfield_range[k])
    } else {
      m_hopfield = 0
      A <- adjacency_graph(n_part, type = graph, anti_corr = anti_corr) 
    }
    
    init <- matrix(NA, nrow = n_chains, ncol = n_part)
    
    if (q == 2) {
      for (c in 1:n_chains) init[c, ] <- sample(c(-1, 1), n_part, replace = T)
    } else {
      for (c in 1:n_chains) init[c, ] <- sample(1:q, n_part, replace = T)
    }

    output_name_root <- paste0(n_part, graph, "_", beta, "_")

    if ("AG" %in% algo_names) {
      set.seed(seed)
      algo_index = which(algo_names == "AG")
      samples <- array(NA, c(n_samples_algo[algo_index], n_part, n_chains))
      print("Doing AG sampling")
      print(paste0("Chain: ", 0, " / ", n_chains))
      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , c] <- ag_sampler(A, beta, init[c, ], q,
                                       n_samples_algo[algo_index],
                                       perturb = perturbation_step)
          # samples[, , c] <- sumit_sampler(A, beta, B, init[c, ],
          #                                 n_iter = n_samples_algo[algo_index],
          #                                 perturb = perturbation_step)
          print(paste0("Chain: ", c, " / ", n_chains))
          })[3]

      save_output(samples, time, algo_index, beta, A, n_chains, output_name_root,
                  algo_names, perturbation_step, q, m_hopfield)
    }

    if ("AG_lowrank" %in% algo_names) {
      if (graph == "hopfield") {
        epsilon = 1e-6
      } else {
        epsilon = 1e-6
      }
      
      set.seed(seed)
      algo_index = which(algo_names == "AG_lowrank")
      samples <- array(NA, c(n_samples_algo[algo_index], n_part, n_chains))
      print("Doing AG lowrank sampling.")
      print(paste0("Chain: ", 0, " / ", n_chains))
      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , c] <- ag_sampler(A, beta, init[c, ], q,
                                       n_samples_algo[algo_index],
                                       lowrank = TRUE,
                                       perturb = perturbation_step,
                                       epsilon = epsilon)
          # samples[, , c] <- ag_ising_lowrank(A = A, beta = beta,
          #                                    init = init[c, ],
          #                                    n_iter = n_samples_algo[algo_index],
          #                                    epsilon = 1e-12)
          print(paste0("Chain: ", c, " / ", n_chains))
        })[3]

      save_output(samples, time, algo_index, beta, A, n_chains, output_name_root,
                  algo_names, perturbation_step, q, m_hopfield)
    }
    
    if ("HB" %in% algo_names) {
      set.seed(seed)
      algo_index <- which(algo_names == "HB")
      print("Doing heatbath sampling.")
      print(paste0("Chain: ", 0, " / ", n_chains))
      samples <- array(NA, c(n_samples_algo[algo_index], n_part, n_chains))
      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , c] <- hb_sampler(A, beta, B, init[c, ],
                                       n_samples_algo[algo_index],
                                       sub_sample = 1, is_potts = is_potts,
                                       n_states = q)
          print(paste0("Chain: ", c, " / ", n_chains))
        })[3]

      save_output(samples, time, algo_index, beta, A, n_chains, output_name_root,
                  algo_names, q=q, m_hopfield=m_hopfield)
    }
    
    if ("HB_long" %in% algo_names) {
      set.seed(seed)
      algo_index <- which(algo_names == "HB_long")
      print("Doing heatbath long sampling.")
      print(paste0("Chain: ", 0, " / ", n_chains))
      samples <- array(NA, c(n_samples_algo[algo_index], n_part, n_chains))
      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , c] <- hb_sampler(A, beta, B, init[c, ],
                                       n_samples_algo[algo_index],
                                       sub_sample = ceiling(n_part / 10),
                                       is_potts = is_potts,
                                       n_states = q)
          print(paste0("Chain: ", c, " / ", n_chains))
        })[3]
      
      save_output(samples, time, algo_index, beta, A, n_chains, output_name_root,
                  algo_names, q=q, m_hopfield=m_hopfield)
    }
  
    if ("BHB" %in% algo_names) {
      set.seed(seed)
      algo_index <- which(algo_names == "BHB")
      print("Doing blackbox heatbath sampling.")
      print(paste0("Chain: ", 0, " / ", n_chains))
      samples <- array(NA, c(n_samples_algo[algo_index], n_part, n_chains))
      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , c] <- hb_sampler(A, beta, B, init[c, ],
                                       n_samples_algo[algo_index],
                                       sub_sample = 1,
                                       is_potts = is_potts,
                                       n_states = q,
                                       black_box = TRUE)
          print(paste0("Chain: ", c, " / ", n_chains))
        })[3]
      save_output(samples, time, algo_index, beta, A, n_chains, output_name_root,
                  algo_names, q=q, m_hopfield=m_hopfield)
    }
  
    # if ("SW" %in% algo_names) {
    #   set.seed(seed)
    #   algo_index <- which(algo_names == "SW")
    #   print("Doing SW sampling.")
    #   time <- c(time, system.time(
    #     for (c in 1:n_chains) {
    #       print(paste0("Chain: ", 0, " / ", n_chains))
    #       samples[, , algo_index, c] <- sw_sampler(A, beta, B, init[c, ], n_sample)
    #     })[3])
    # }
  
    if ("Wolff" %in% algo_names) {
      set.seed(seed)
      algo_index <- which(algo_names == "Wolff")
      print("Doing Wolff sampling.")
      print(paste0("Chain: ", 0, " / ", n_chains))
      samples <- array(NA, c(n_samples_algo[algo_index], n_part, n_chains))
      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , c] <- wolff_sampler(A, beta, B, init[c, ],
                                          n_samples_algo[algo_index],
                                          is_potts, q)
          print(paste0("Chain: ", c, " / ", n_chains))
        })[3]
      
      save_output(samples, time, algo_index, beta, A, n_chains, output_name_root,
                  algo_names, q=q, m_hopfield=m_hopfield)
    }
  
    if ("GWG" %in% algo_names) {
      set.seed(seed)
      algo_index <- which(algo_names == "GWG")
      print("Doing Gibbs with Gradient sampling")
      print(paste0("Chain: ", 0, " / ", n_chains))
      samples <- array(NA, c(n_samples_algo[algo_index], n_part, n_chains))
      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , c] <- gwg_sampler(A = A, beta = beta,
                                        init = init[c, ],
                                        n_states = q,
                                        n_iter = n_samples_algo[algo_index])
          print(paste0("Chain: ", c, " / ", n_chains))
        })[3]
  
      save_output(samples, time, algo_index, beta, A, n_chains, output_name_root,
                  algo_names, q=q, m_hopfield=m_hopfield)
    }
    
    if ("dHMC" %in% algo_names) {
      # WARNING: dHMC is only implemented for the Ising model!
      set.seed(seed)
      algo_index <- which(algo_names == "dHMC")
      print("Doing discrete HMC sampling")
      print(paste0("Chain: ", 0, " / ", n_chains))
      samples <- array(NA, c(n_samples_algo[algo_index], n_part, n_chains))
      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , c] <- discrete_hmc(A = A, beta = beta,
                                                    init = init[c, ], 
                                                    n_iter = n_samples_algo[algo_index])
          print(paste0("Chain: ", c, " / ", n_chains))
        })[3]
  
      save_output(samples, time, algo_index, beta, A, n_chains, output_name_root,
                  algo_names, q=q, m_hopfield=m_hopfield)
    }
    
    
    # Initialization for tempering algorithms
    n_beta_temp <- length(beta_vector)
    init_temp <- array(NA, c(n_beta_temp, n_chains, n_part))
    
    for (i in 1:n_beta_temp) {
      if (q == 2) {
        for (c in 1:n_chains) init_temp[, c, ] <- sample(c(-1, 1), n_part, replace = T)
      } else {
        for (c in 1:n_chains) init_temp[, c, ] <- sample(1:q, n_part, replace = T)
      }
    }

    if ("temp_AG" %in% algo_names) {
      set.seed(seed)
      algo_index <- which(algo_names == "temp_AG")
      print("Doing temp AG sampling")
      print(paste0("Chain: 0 / ", n_chains))
      samples <- array(NA, c(n_beta_temp, n_mc_iter[algo_index] * n_exchange,
                             n_part, n_chains))

      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , , c] <- tempering(A, beta_vector, init_temp[, c, ],
                                        n_exchange,  n_mc_iter[algo_index], B,
                                        f_AG, is_potts = is_potts,
                                        n_states = q)
        }
      )

      for (i in 1:length(beta_vector)) {
        output_name_root <- paste0(n_part, graph, "_", beta_vector[i], "_")
        save_output(samples[i, , , ], time, algo_index, beta_vector[i], A,
                    n_chains, output_name_root, algo_names, q=q)
      }
    }
    
    
    if ("temp_HB_long" %in% algo_names) {
      set.seed(seed)
      algo_index <- which(algo_names == "temp_HB_long")
      print("Doing temp HB_long sampling")
      print(paste0("Chain: 0 / ", n_chains))
      samples <- array(NA, c(n_beta_temp, n_mc_iter[algo_index] * n_exchange,
                             n_part, n_chains))
      
      sub_adjust <- 10
      f_HB_long <- function(A, beta, B, init, n_iter, n_states) {
        hb_sampler(A = A, beta = beta, init = init,
                   n_iter = n_iter, is_potts = TRUE, 
                   n_states = n_states,
                   sub_sample = ceiling(n_part / sub_adjust))
      }

      time <- system.time(
        for (c in 1:n_chains) {
          samples[, , , c] <- tempering(A, beta_vector, init_temp[, c, ],
                                        n_exchange,  n_mc_iter[algo_index], B,
                                        f_HB_long, is_potts = is_potts,
                                        n_states = q)
        }
      )
      
      for (i in 1:length(beta_vector)) {
        output_name_root <- paste0(n_part, graph, "_", beta_vector[i], "_")
        save_output(samples[i, , , ], time, algo_index, beta_vector[i], A,
                    n_chains, output_name_root, algo_names, q=q)
      }      
    }
    
  }
}

###############################################################################
## Run experiments

# Ising complete
run_simulations(beta_range = c(0.5, 1, 1.5),
                graph_types = rep("complete", 3),
                graph_sizes = rep(24^2, 3),
                n_samples_algo = c(5e4, 5e4, 1e5, 5e4, 5e4, 5e3, 1e5, 1e3),
                algo_names = c("AG", "AG_lowrank", "HB", "HB_long", "BHB",
                               "Wolff", "GWG", "dHMC"))

# Ising grid
run_simulations(beta_range = c(0.5, 0.65, 1 / sqrt(2), 0.75, 0.8),
                graph_types = rep("grid", 5),
                graph_sizes = rep(24^2, 5),
                n_samples_algo = c(5e4, 1e5, 5e4, 5e4, 5e3, 1e5, 1e3),
                algo_names = c("AG", "HB", "HB_long", "BHB",
                               "Wolff", "GWG", "dHMC"))

# Potts grid
run_simulations(beta_range = c(0.6),  # c(0.5, 0.55, 0.6) 
                graph_types = rep("grid", 1),
                graph_sizes = rep(24^2, 2),
                n_samples_algo = c(5e4, 1e5, 1e5, 5e4, 5e4, 1e4),
                algo_names = c("AG", "HB", "HB_long", "BHB", "Wolff", "GWG"),
                q = 4)

# Potts complete
run_simulations(beta_range = c(1.7),  # c(1.2, 1.64),
                graph_types = rep("complete", 1),
                graph_sizes = rep(24^2, 1),
                n_samples_algo = c(5e4, 5e4, 1e5, 1e5, 5e4, 5e4, 1e4, 1e3),
                algo_names = c("AG", "AG_lowrank", "HB", "HB_long", "BHB", "Wolff", "GWG",
                               "dHMC"),
                q = 4)



# Hopfield
m_hopfield_range <- c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200)
run_simulations(beta_range = rep(1, length(m_hopfield_range)),
                graph_types = rep("hopfield", length(m_hopfield_range)),
                graph_sizes = rep(16^2, length(m_hopfield_range)),
                n_samples_algo = c(5e4, 5e4, 5e4), # c(5e4, 5e4, 5e4),
                algo_names = c("AG", "AG_lowrank", "HB_long"),
                q = 4,
                anti_corr = TRUE,
                m_hopfield_range = m_hopfield_range)


# Potts spin glass
beta_vector <- seq(from = 0.5, to = 3, by = 0.25)
run_simulations(beta_range = c(3),
                graph_types = c("gaussian"),
                graph_sizes = c(128),
                n_samples_algo = c(0, 0),
                algo_names = c("temp_AG"), #, "temp_HB_long"),
                q = 4,
                beta_vector = beta_vector,
                n_exchange = 40,
                n_mc_iter = c(2e3, 2e3))


# Potts spin glass: additional runs for non-tempered algorithm
beta_range <- seq(from = 0.5, to = 3, by = 0.25)
run_simulations(beta_range = beta_range,
                graph_types = rep("gaussian", length(beta_range)),
                graph_sizes = rep(128, length(beta_range)),
                n_samples_algo = c(8e4, 0, 0),
                algo_names = c("AG"),
                q = 4)


###############################################################################
## Test and patchups...
run_simulations(beta_range = c(0.55),
                graph_types = rep("grid", 1),
                graph_sizes = rep(24^2, 1),
                n_samples_algo = c(1e5, 1e4),
                algo_names = c("HB", "GWG"),
                q = 4)


run_simulations(beta_range = c(0.25),
                graph_types = rep("grid", 1),
                graph_sizes = rep(24^2, 1),
                n_samples_algo = rep(1e2, 7),
                algo_names = c("AG", "HB", "HB_long", "BHB",
                               "Wolff", "GWG"),
                q = 3)



run_simulations(beta_range = c(0.5, 0.55, 0.6),
                graph_types = rep("grid", 3),
                graph_sizes = rep(24^2, 3),
                n_samples_algo = c(1e3),
                algo_names = c("dHMC"),
                q = 4)


run_simulations(beta_range = c(1.2, 1.64, 1.7),  # c(1.2, 1.64),
                graph_types = rep("complete", 3),
                graph_sizes = rep(24^2, 3),
                n_samples_algo = c(1e4),
                algo_names = c("GWG"),
                q = 4)


run_simulations(beta_range = c(1.2, 1.64, 1.7),  # c(1.2, 1.64),
                graph_types = rep("complete", 3),
                graph_sizes = rep(24^2, 3),
                n_samples_algo = c(5e3),
                algo_names = c("Wolff"),
                q = 4)

run_simulations(beta_range = c(1.64, 1.7),  # c(1.2, 1.64),
                graph_types = rep("complete", 2),
                graph_sizes = rep(24^2, 2),
                n_samples_algo = c(5e4),
                algo_names = c("AG_lowrank"),
                q = 4)

beta_range = c(1.2)  # c(1.2, 1.64),
graph_types = rep("complete", 1)
graph_sizes = rep(24^2, 1)
n_samples_algo = c(1e4)
algo_names = c("AG_lowrank")
q = 4

anti_corr = FALSE
perturbation_step = FALSE
m_hopfield_range = NULL
beta_vector = seq(from = 0.5, to = 3, by = 0.25)
n_exchange = 40
n_mc_iter = 1e3


###############################################################################
## Draft code

# run_simulations <- function(beta_range, graph_types, graph_sizes, n_samples_algo,
#                             algo_names, burnin_fraction = 0.5, n_chains = 4, 
#                             anti_corr = FALSE, perturbation_step = FALSE, q = 2)
# 
# 
# # algo_names <- c("AG",  "HB", "HB_long", "BHB", "Wolff", "GWG", "dHMC")
# # n_sample <- 5e4
# # n_samples_algo <- rep(n_sample, length(algo_names))
# # n_samples_algo[which(algo_names == "GWG")] <- 10e4
# # n_samples_algo[which(algo_names == "HB")] <- 10e4
# # n_samples_algo[which(algo_names == "Wolff")] <- 5e3
# # n_samples_algo[ which(algo_names == "dHMC")] <- 1e3
# 
#   
# 
# # beta_range <- c(log(1 + sqrt(2)) / 2, log(1 + sqrt(2)) / 2)
# # beta_range <- c(0.5, 0.65, 0.69, 1 / sqrt(2), 0.72, 0.75, 0.9)  # Grid graph
# # beta_range <- c(0.5, 0.9, 1, 1.1, 2) # Complete graph
# 
# # beta_range <- c(0.5, 1, 1.5)  # Complete graph
# # beta_range <- c(0.5, 0.65, 1/sqrt(2), 0.8, 1)
# beta_range <- c(0.75)
# n_beta <- length(beta_range)
# B_range <- rep(0, n_beta) 
# graph_types <- rep("grid", n_beta)
# anti_corr <- FALSE
# graph_sizes <- rep(24^2, n_beta)  # 16^2
# n_sample <- 100 # 5e4
# n_chains <- 4
# total_sample <- n_sample * n_chains
# burnin_fraction <- 0.5
# 
# 
# algo_names <- c("AG",  "HB", "HB_long", "BHB", "Wolff", "GWG", "dHMC")
# n_algorithms <- length(algo_names)
# 
# perturbation_step <- FALSE
# q = 4
# 
# # Some algorithms have a large ESS/iteration, but each iteration takes
# # a long time, so it makes sense to run less iterations.
# n_samples_algo <- rep(n_sample, n_algorithms)
# n_samples_algo[which(algo_names == "GWG")] <- 10e4
# n_samples_algo[which(algo_names == "HB")] <- 10e4
# n_samples_algo[which(algo_names == "Wolff")] <- 5e3
# n_samples_algo[ which(algo_names == "dHMC")] <- 1e3
# 
# n_gen_quant <- 3
# 
# deliv_dir <- file.path(getwd(), "deliv_v3")



