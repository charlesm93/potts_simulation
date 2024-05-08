
rm(list = ls())
gc()
set.seed(1954)

# Adjust to your settings.
.libPaths("~/Rlib/")
setwd("~/Code/ising_model/scripts/")

library(ggplot2)
library(tidyverse)
library(posterior)
library(latex2exp)



deliv_dir <- file.path(getwd(), "deliv_v3")

# options: "Ising grid", "Potts grid", "Ising complete", "Potts complete",
#          "Spin Glass"
experiment <- "Spin Glass"

plot_temp <- FALSE
if (experiment == "Ising grid") {
  cbPalette <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#999999")
  n_part <- 24^2; graph <- "grid"; q <- 2;
  beta_range <- c(0.65, 1/sqrt(2), 0.75)
  algo_names <- c("AG", "HB", "HB_long", "BHB", "Wolff", "GWG", "dHMC")
  n_sample_vec <- c(5e4, 10e4, 5e4, 5e4, 5e3, 10e4, 1e3)
  
  algo_to_emphasize <- c("AG")
}

if (experiment == "Potts grid") {
  cbPalette <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#999999")
  n_part <- 24^2; graph <- "grid"; q <- 4;
  beta_range <- c(0.5, 0.55, 0.6)
  algo_names <- c("AG", "HB", "HB_long", "BHB", "Wolff", "GWG", "dHMC")
  n_sample_vec = c(5e4, 1e5, 1e5, 5e4, 5e4, 1e4, 1e3)
  
  algo_to_emphasize <- c("AG", "AG_lowrank")
}

if (experiment == "Ising complete") {
  cbPalette <-c("#D55E00", "#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                "#999999")
  n_part <- 24^2; graph <- "complete"; q <- 2;
  beta_range <- c(0.5, 1, 1.5)
  algo_names <- c("AG", "AG_lowrank", "HB", "HB_long", "BHB",
                  "Wolff", "GWG", "dHMC")
  n_sample_vec <- c(5e4, 5e4, 1e5, 5e4, 5e4, 5e3, 1e5, 1e3)
  
  algo_to_emphasize <- c("AG", "AG_lowrank")  
}

if (experiment == "Potts complete") {
  cbPalette <-c("#D55E00", "#CC79A7", "#E69F00", "#009E73", "#F0E442", "#0072B2",
                "#999999")
  n_part <- 24^2; graph <- "complete"; q <- 4;
  beta_range <- c(1.2, 1.64, 1.7)
  algo_names <- c("AG", "AG_lowrank", "HB", "HB_long", "BHB", "Wolff", "GWG")
  n_sample_vec <- c(5e4, 5e4, 1e5, 1e5, 5e4, 1e3, 5e3)

  algo_to_emphasize <- c("AG", "AG_lowrank", "Wolff")  
}


if (experiment == "Hopfield") {
  cbPalette <-c("#D55E00", "#CC79A7", "#0072B2")
  n_part <- 16^2; graph <- "hopfield"; q <- 4;
  beta_range <- c(1)
  algo_names <- c("AG", "AG_lowrank", "HB_long")
  n_sample_vec <- c(5e4, 5e4, 5e4)
  
  algo_to_emphasize <- c("AG", "AG_lowrank")  
}

if (experiment == "Spin Glass") {
  cbPalette <- c("#D55E00", "#F0E442", "#0072B2")
  n_part <- 128; graph <- "gaussian"; q <- 4;
  beta_range <- seq(from = 0.5, to = 3, by = 0.25)
  algo_names <- c("AG", "temp_AG", "temp_HB_long")
  n_sample_vec <- c(8e4, 8e4, 8e4)
  
  algo_to_emphasize <- c("AG", "temp_AG")
  
  plot_temp <- TRUE
}

n_beta <- length(beta_range)
m_hopfield_range <- c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200)
n_algorithms <- length(algo_names)
perturbation <- rep(FALSE, n_algorithms)
n_sample <- max(n_sample_vec)
n_chains <- 4

# n_part <- 24^2 # options: 128, 16^2, 24^2
# graph <- "grid" # options: "grid", "complete", "hopfield", "gaussian
# q <- 4

# beta_range <- c(0.5, 1, 1.5)
# beta_range <- c(0.5, 0.65, 1/sqrt(2), 0.75, 0.8, 1)  # Ising grid
# beta_range <- c(0.55)  # Potts grid
# beta_range <- c(0.71)
# beta_range <- seq(from = 0.5, to = 3, by = 0.25)




# algo_names <- c("AG", "AG_lowrank", "HB", "HB_long", "BHB", "Wolff", "GWG", "dHMC")  # complete
# algo_names <- c("AG", "HB", "HB_long", "BHB", "Wolff", "GWG", "dHMC")  # Ising grid
# algo_names <- c("AG", "HB", "HB_long", "BHB", "Wolff", "GWG")  # potts grid
# algo_names <- c("AG", "AG_lowrank", "HB_long")
# algo_names <- c("AG", "temp_AG", "temp_HB_long")

# n_algorithms <- length(algo_names)
# perturbation <- rep(FALSE, n_algorithms)

# n_sample_vec <- c(5e4, 5e4, 10e4, 5e4, 5e4, 5e3, 10e4, 1e3)  # complete
# n_sample_vec <- c(5e4, 10e4, 5e4, 5e4, 5e3, 10e4, 1e3) # Ising grid
# n_sample_vec <- c(5e4, 1e5, 1e5, 5e4, 5e4, 1e4)
# n_sample_vec <- c(5e4, 5e4, 5e4, 5e4, 5e4, 5e4)
# n_sample_vec <- c(5e4, 1e5, 1e5, 5e4, 5e4, 1e4)
# n_sample_vec <- c(8e4, 8e4, 8e4)

###############################################################################

if (graph == "hopfield") {
  n_cases <- length(m_hopfield_range)
} else {
  n_cases <- n_beta
}

energy <- array(NA, c(n_cases, n_algorithms, n_sample, n_chains))
time <- array(NA, c(n_cases, n_algorithms))
time_per_iter <- array(NA, c(n_cases, n_algorithms))

for (j in 1:n_cases) {
  if (graph == "hopfield") {
    beta_char <- beta_range[1]
  } else {
    beta_char <- beta_range[j]
  }

  output_name_root <- paste0(n_part, graph, "_", beta_char, "_")

  for (i in 1:length(algo_names)) {
    file_name <- paste0(output_name_root, algo_names[i], "_perturb_",
                        perturbation[i])
    if (q >= 3) file_name <- paste0(file_name, "_q=", q)
    if (graph == "hopfield") file_name <- paste0(file_name, "_m=", m_hopfield_range[j])
    file_name <- paste0(file_name, ".RData")

    saved_list <- readRDS(file = paste0(deliv_dir, "/", file_name))
    energy[j, i, 1:n_sample_vec[i],] <-
      saved_list$gen_quant[1:n_sample_vec[i], 1:n_chains, 3]

    if (algo_names[i] %in% c("temp_AG", "temp_HB_long")) {
      time[j, i] <- saved_list$time[3]
    } else {
      time[j, i] <- saved_list$time
    }
    
    time_per_iter[j, i] <- time[j, i] / n_sample_vec[i]
  }
}

if (plot_temp) {
  algo_index <- which(algo_names == "AG")
  time_total <- sum(time[, algo_index])
  for (j in 1:n_cases) {
    time[j,  algo_index] <- time_total
    time_per_iter[j, algo_index] <- time[j, algo_index] / n_sample_vec[algo_index]
  }
}


checkpoints <- c(4,  seq(from = 10, to = 400, by = 10), 
                 seq(from = 500, to = n_sample, by = 500))
n_checkpoints <- length(checkpoints)

checkpoints_time <- array(NA, dim = c(n_cases, n_checkpoints, n_algorithms))
for (j in 1:n_cases) {
  for (a in 1:n_algorithms) {
    checkpoints_time[j, , a] <- checkpoints * time_per_iter[j, a]
  }
}

rhat_vec <- array(NA, c(n_cases, n_checkpoints, n_algorithms))
for (j in 1:n_cases) { 
  for (a in 1:n_algorithms) {
    for (i in 1:n_checkpoints) {
      if (checkpoints[i] <= n_sample_vec[a]) {
        n_burnin <- floor(checkpoints[i] / 2)
        rhat_vec[j, i, a] <- rhat_nested(energy[j, a, n_burnin:checkpoints[i], ],
                                         superchain_ids = 1:n_chains)
        
        last_index = i
      } else {
        rhat_vec[j, i, a] = rhat_vec[j, last_index, a]
        checkpoints_time[j, i, a] = checkpoints_time[j, last_index, a]
      }
    }
  }
}

plot_data <-
  data.frame(time = c(checkpoints_time),
             iteration = rep(rep(checkpoints, each = n_cases), n_algorithms),
             rhat = c(rhat_vec) - 1,
             algorithm = rep(algo_names, each = n_cases * n_checkpoints))


if (graph == "hopfield") {
  plot_data$beta <- paste0("b=", rep(beta_range[1], n_checkpoints * n_cases * n_algorithms))
  plot_data$rank <- paste0("m=", rep(m_hopfield_range, n_algorithms * n_checkpoints))
} else {
  plot_data$beta <- paste0("b=", rep(trunc(beta_range * 10^3) / 10^3,
                                     n_algorithms * n_checkpoints))
  plot_data$rank <- paste0("m=", rep(m_hopfield_range[1], n_checkpoints * n_cases * n_algorithms))
}

plot_data$alpha <- rep(0.75, n_algorithms * n_checkpoints)

# algo_to_emphasize <- c("AG", "AG_lowrank")  # complete
# algo_to_emphasize <- c("AG")  # grid ising, grid potts
# algo_to_emphasize <- c("temp_AG")  # tempering

plot_data[plot_data$algorithm %in% algo_to_emphasize, ]$alpha <- 1

scaleFUN <- function(x) sprintf("%.3f", x)

# Plot by inverse temperature beta
p <- ggplot(data = plot_data, aes(x = time, y = rhat, color = algorithm)) +
  geom_line(aes(alpha=alpha), linewidth=1) + 
  scale_y_continuous(trans='log10') +
  facet_wrap(~beta) +
  theme_bw() + geom_hline(yintercept=0.01, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.1, linetype="dashed", color = "grey") +
  scale_color_manual(values= cbPalette) +
  ylab(TeX("$\\widehat{R}$ - 1")) + xlab("Wall time (s)") +
  scale_x_continuous(trans='log2', labels=scaleFUN) +
  scale_alpha(range = c(0.5, 1)) +  guides(alpha = "none") +
  theme(text = element_text(size = 15))
p

# Plot by rank of hopfield graph
if (graph == "hopfield") {
  p <- ggplot(data = plot_data, aes(x = time, y = rhat, color = algorithm)) +
    geom_line(aes(alpha=alpha), linewidth=1) + 
    scale_y_continuous(trans='log10') +
    facet_wrap(~rank) +
    theme_bw() + geom_hline(yintercept=0.01, linetype="dashed", color = "black") +
    geom_hline(yintercept=0.05, linetype="dashed", color = "grey") +
    scale_color_manual(values= cbPalette) +
    ylab(TeX("$\\widehat{R}$ - 1")) + xlab("Wall time (s)") +
    scale_x_continuous(trans='log2', labels=scaleFUN) +
    scale_alpha(range = c(0.5, 1)) +  guides(alpha = "none")
  p
}



energy_mean <- array(NA, c(n_cases, n_algorithms))
ess <- array(NA, c(n_cases, n_algorithms))
ess_s <- array(NA, c(n_cases, n_algorithms))
rhat_final <- array(NA, c(n_cases, n_algorithms))
for (j in 1:n_cases) {
  for (i in 1:n_algorithms) {
    energy_mean[j, i] <- mean(energy[j, i, 1:n_sample_vec[i], ]) 
    ess[j, i] <- ess_bulk(energy[j, i, (n_sample_vec[i] / 2):n_sample_vec[i], ])
    ess_s[j, i] <- ess[j, i] / time[j, i]
    rhat_final[j, i] <- rhat_vec[j, n_checkpoints, i]
  }
}

plot_data <- data.frame(ess_s = c(ess_s),
                        algorithm = rep(algo_names, each = n_cases),
                        rhat_final = c(rhat_final))

if (graph == "hopfield") {
  plot_data$beta <- rep(beta_range[1], n_cases * n_algorithms)
  plot_data$rank <- rep(m_hopfield_range, n_algorithms)
} else {
  plot_data$beta <- rep(beta_range, n_algorithms)
  plot_data$rank <- rep(m_hopfield_range[1], n_cases * n_algorithms)
}

## Tuning parameters for what constitutes an acceptable convergence
if (plot_temp) {
  rhat_threshold <- 1.15
  min_ess <- 1e-4
} else {
  rhat_threshold <- 1.15
  min_ess <- 1e-2
}

plot_data$eff_if_conv <- plot_data$ess_s * (plot_data$rhat_final <= rhat_threshold) + min_ess
plot_data$beta <- as.factor(plot_data$beta)
if (graph == "grid") {
  levels(plot_data$beta)[which(beta_range == 1 / sqrt(2))] <- "0.71"
}


p_bar_ising <- ggplot(plot_data, aes(x = beta, y = eff_if_conv, fill = algorithm)) +
  geom_col(color = 'black', position = position_dodge()) + theme_bw() +
  scale_fill_manual(values= cbPalette) +
  scale_y_continuous(trans='sqrt') +
  xlab(" ") + xlab(TeX("$\\beta$")) +
  theme(text = element_text(size = 13)) + ggtitle(paste0("q = ", q)) +
  ylab(TeX("$ESS / s"))
p_bar_ising

plot_data$beta <- rep(beta_range, n_algorithms)
p_line <- ggplot(plot_data[plot_data$rhat_final < rhat_threshold, ], 
                 aes(x = beta, y = ess_s, color = algorithm)) +
  geom_line(linewidth=1) + geom_point(size = 2) + theme_bw() +
  scale_color_manual(values= cbPalette) +
  scale_y_continuous(trans='log2') +
  xlab(" ") + xlab(TeX("$\\beta$")) +
  theme(text = element_text(size = 13)) + ggtitle(paste0("q = ", q)) +
  ylab(TeX("$ESS / s"))

p_line

if (graph == "hopfield") {
  p_bar_ising <- ggplot(plot_data, aes(x = rank, y = eff_if_conv, fill = algorithm)) +
    geom_col(color = 'black', position = position_dodge()) + theme_bw() +
    scale_fill_manual(values= cbPalette) +
    scale_y_continuous(trans='sqrt') +
    xlab(" ") + xlab(TeX("$\\beta$")) +
    theme(text = element_text(size = 13)) + ggtitle(paste0("q = ", q)) +
    ylab(TeX("$ESS / s"))

  p_bar_ising
  
  p_line <- ggplot(plot_data, aes(x = rank, y = eff_if_conv, color = algorithm)) +
    geom_line(linewidth=2) + geom_point(size = 3) + theme_bw() +
    scale_color_manual(values= cbPalette) +
    xlab(" ") + xlab("rank of A") +
    theme(text = element_text(size = 13)) + ggtitle(paste0("q = ", q)) +
    ylab("ESS / s") # + scale_y_continuous(trans = 'log2')
  p_line
}

if (TRUE) {
  print(algo_names)
  print(energy_mean)
  print(ess)
  print(ess_s)
}
