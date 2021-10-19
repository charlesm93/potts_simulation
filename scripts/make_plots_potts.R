
# generate plots for the manuscript

rm(list = ls())
gc()

setwd("~/Code/ising_model/scripts")
.libPaths("~/Rlib/")

library(ggplot2)
library(latex2exp)
library(gridExtra)
library(scales)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 <- c("#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

################################################################
## complete "Curry Weiss" model
deliv_dir <- file.path("deliv", "numerical_experiment", "potts_complete")

n_algorithm <- 4
n_states <- 4
# beta_range <- c(1.35, 1.55, 1.6, 1.609, 1.647918, 1.75, 1.95)
beta_range <- c(1.35, 1.55, 1.647918, 1.75, 1.95)
beta_range <- beta_range[1:3]  # After 1.65, results become inacurate for
                               # random walk methods.

eff_data <- c()
for (i in 1:length(beta_range)) {
  eff_file <- read.csv(file = file.path(deliv_dir,
                                paste0("256complete_", beta_range[i],
                                       "_", n_states, "_efficiency.csv")))
  eff_file <- eff_file[1, ]

  eff_data <- c(eff_data, eff_file$eff_AG_potts,
                eff_file$eff_HB,
                eff_file$eff_Wolff,  # , WOLFF
                # eff_file$eff_SW)
                eff_file$eff_AG_potts_lowrank)
}

method <- rep(c("AG", "Heat Bath", "Wolff", "AG (low)"), length(beta_range))
method <- factor(method, level = c("AG", "Wolff", "Heat Bath", "AG (low)"))

# beta <- rep(c(2, 1, 0.5), each = 6)
beta <- rep(beta_range, each = n_algorithm)
parm <- rep("K(x)", length(eff_data))
parm_labels <- rep(c(TeX("$K(x)$")))
parm <- factor(parm, label = parm_labels)

p_potts <- ggplot(data.frame(eff = eff_data, beta = beta,
                       method = method, parm = parm),
            aes(x = beta, y = eff, linetype = method,
                color = method)) +
  # geom_point() + 
  geom_line(size = 1.5) + theme_bw() +
  # facet_wrap(~parm, labeller = "label_parsed") + 
  scale_y_continuous(trans='log2') +
  xlab(TeX("$\\beta$")) +
  ylab("") +
  theme(text = element_text(size = 13)) +
  scale_color_manual(values= cbPalette) + ggtitle("q = 4")
p_potts

# Ising model
beta_range <- c(0.5, 1, 2)
eff_data <- c()
deliv_dir <- file.path("deliv", "numerical_experiment", "ising_complete")
for (i in 1:length(beta_range)) {
  eff_file <- read.csv(file = file.path(deliv_dir,
                                        paste0("256complete_0_", beta_range[i],
                                               "_efficiency.csv")))
  eff_file <- eff_file[3, ]
  
  # NOTE: Had to rerun the lowrank algorithm after further optimizing
  # and stored the results in a seperate file.
  eff_file_lowrank <- read.csv(file = file.path(deliv_dir,
                                 paste0("256complete_0_", beta_range[i],
                                 "_efficiency_lowrank.csv")))
  
  eff_file_lowrank <- eff_file_lowrank[3, ]
  
  eff_data <- c(eff_data, eff_file$eff_Sumit,
                eff_file$eff_HB,
                eff_file$eff_Wolff,
                eff_file_lowrank$eff_AG_lowrank)
}

method <- rep(c("AG", "Heat Bath", "Wolff", "AG (low)"), length(beta_range))
method <- factor(method, level = c("AG", "Wolff", "Heat Bath", "AG (low)"))

beta <- rep(beta_range, each = n_algorithm)
p_ising <- ggplot(data.frame(eff = eff_data, beta = beta,
                             method = method),
                  aes(x = beta, y = eff, linetype = method,
                      color = method)) +
  geom_line(size = 1.5) + theme_bw() +
  scale_y_continuous(trans='log2') +
  xlab(TeX("$\\beta$")) +
  ylab(TeX("$m_{eff} / s")) +
  theme(text = element_text(size = 13)) +
  scale_color_manual(values= cbPalette) + ggtitle("q = 2") +
  theme(legend.position = "none")
p_ising

# p_ising <- p_potts + theme(legend.position = "none") +
#   ylab(TeX("$m_{eff} / s"))

grid.arrange(p_ising, p_potts, widths = c(0.8, 1),nrow = 1)  

################################################################
## Lattice (grid) model for q = 2 (Ising) and Potts (q = 4)
deliv_dir <- file.path("deliv", "numerical_experiment", "potts_grid")

## (1) create plot for Potts experiment
n_states <- 4  # Options: 3, 4, 5
n_algorithm <- 3  # (AG, HB, Wolff)
beta_range <- c(0.25, 0.45, 0.5493062, 0.65, 0.85)  # for n_states = 4
# beta_range <- c(0.25, 0.44, 0.63, 1, 1.41421356237309, 1.8)  # n_states = 5
# beta_range <- c(1, 1.41421356237309, 1.8)

eff_data <- c()
for (i in 1:length(beta_range)) {
  eff_file <- read.csv(file = file.path(deliv_dir,
                                        paste0("256grid_", beta_range[i],
                                               "_", n_states, "_efficiency.csv")))
  eff_file <- eff_file[1, ]
  
  eff_data <- c(eff_data, eff_file$eff_AG_potts,
                eff_file$eff_HB,
                eff_file$eff_Wolff)  # , WOLFF
                # eff_file$eff_SW)
}

method <- rep(c("AG", "Heat Bath", "Wolff"), length(beta_range))
method <- factor(method, level = c("AG", "Wolff", "Heat Bath"))

# beta <- rep(c(2, 1, 0.5), each = 6)
beta <- rep(beta_range, each = n_algorithm)
parm <- rep("K(x)", length(eff_data))
parm_labels <- rep(c(TeX("$K(x)$")))
parm <- factor(parm, label = parm_labels)

p_potts <- ggplot(data.frame(eff = eff_data, beta = beta,
                       method = method, parm = parm),
            aes(x = beta, y = eff, linetype = method,
                color = method)) +
  # geom_point() + 
  geom_line(size = 1.5) + theme_bw() +
  # facet_wrap(~parm, labeller = "label_parsed") + 
  scale_y_continuous(trans='log2') +
  ylab("") +
  xlab(TeX("$\\beta$")) +
  # ylab(TeX("$m_{eff} / s")) +
  theme(text = element_text(size = 13)) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  ggtitle("q = 4")
p_potts


## (2) Create plots for Ising model
# TODO: use a smarter way to store the results, as was done for
# Potts models.
deliv_dir <- file.path("deliv", "numerical_experiment")
filename <- "256grid_0_efficiency.csv"
beta_cold <- 
  read.csv(file = file.path(deliv_dir, "grid_beta_cold",
                            filename))
beta_crit <- 
  read.csv(file = file.path(deliv_dir, "grid_beta_crit",
                            filename))
beta_warm <- 
  read.csv(file = file.path(deliv_dir, "grid_beta_warm",
                            filename))

# Only examine calculations of the Hamiltonian
eff_data <- c(beta_cold$eff_Sumit_p[3],
              beta_cold$eff_Wolff[3],
              beta_cold$eff_MH[3],
              beta_crit$eff_Sumit_p[3],
              beta_crit$eff_Wolff[3],
              beta_crit$eff_MH[3],
              beta_warm$eff_Sumit_p[3],
              beta_warm$eff_Wolff[3],
              beta_warm$eff_MH[3])

method <- rep(c("AG", "Wolff", "Metropolis"), 3)
method <- factor(method, level = c("AG", "Wolff", "Metropolis"))

# beta <- rep(c(0.5, 1 / sqrt(2), 0.9), each = 3)
beta <- rep(c(0.9, 1 / sqrt(2), 0.5), each = 3)
# parm <- rep(c("M^2", "K"), 9)
# parm_labels <- rep(c(TeX("$M(X)^2$"), TeX("$K(X)$")))
# parm <- factor(parm, label = parm_labels)

p_ising <- ggplot(data.frame(eff = eff_data, beta = beta,
                       method = method),
            aes(x = beta, y = eff, linetype = method,
                color = method)) +
  # geom_point() + 
  geom_line(size = 1.5) + theme_bw() +
  scale_y_continuous(trans='log2') +
  xlab(TeX("$\\beta$")) +
  ylab(TeX("$m_{eff} / s")) +
  theme(text = element_text(size = 13)) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme(legend.position = "none") + ggtitle("q = 2")
p_ising

grid.arrange(p_ising, p_potts, widths = c(0.8, 1),nrow = 1)


#####################################################################
## Hopfield model
deliv_dir <- file.path("deliv", "numerical_experiment", "potts_hopfield")
n_states <- 4
n_algorithm <- 3
m_range <- c(1, seq(from = 10, to = 100, by = 10), 150, 200, 256)
# m_range <- seq(from = 10, to = 100, by = 10)
beta_range <- rep(1, length(m_range))

eff_data <- c()
for (i in 1:length(m_range)) {
  eff_file <- read.csv(file = file.path(deliv_dir,
                                  paste0("256hopfield_", beta_range[i], "_", n_states,
                                         "_", m_range[i], "_efficiency.csv")))
  eff_file <- eff_file[1, ]
  
  eff_data <- c(eff_data, eff_file$eff_AG_potts,
                eff_file$eff_HB, eff_file$eff_AG_potts_lowrank)
}

method <- rep(c("AG", "Heat Bath", "AG (low)"), length(beta_range))
method <- factor(method, level = c("AG", "Heat Bath", "AG (low)"))

m_hopfield <- rep(m_range, each = n_algorithm)

p_hopfield <- ggplot(data.frame(eff = eff_data, method = method,
                                m_hopfield = m_hopfield),
                     aes(x = m_hopfield, y = eff, linetype = method,
                         color = method)) +
  # geom_point() + 
  geom_line(size = 1.5) + theme_bw() +
  theme(text = element_text(size = 13)) +
  # scale_y_continuous(trans = "log2") +
  xlab(TeX("Rank of $A_n$")) + ylab(TeX("$m_{eff} / s")) +
  scale_color_manual(values= cbPalette2) +
  theme(legend.position = c(0.8, 0.75)) +
  theme(legend.title = element_blank())

p_hopfield

ess_data <- c()
for (i in 1:length(m_range)) {
  ess_file <- read.csv(file = file.path(deliv_dir,
                paste0("256hopfield_", beta_range[i], "_", n_states,
                "_", m_range[i], "_efficiency.csv")))
  ess_file <- ess_file[1, ]
  
  ess_data <- c(ess_data, ess_file$ess_AG_potts,
                ess_file$eff_HB, ess_file$ess_AG_potts_lowrank)
}

p_hopfield <- ggplot(data.frame(ess = ess_data, method = method,
                                m_hopfield = m_hopfield),
                     aes(x = m_hopfield, y = ess, linetype = method,
                         color = method)) +
  # geom_point() + 
  geom_line(size = 1.5) + theme_bw() +
  theme(text = element_text(size = 13)) +
  # scale_y_continuous(trans = "log2") +
  xlab(TeX("Rank of $A_n$")) + ylab(TeX("$m_{eff}")) +
  scale_color_manual(values= cbPalette2) +
  theme(legend.position = c(0.8, 0.75)) +
  theme(legend.title = element_blank())

p_hopfield


#####################################################################
## SK Model (cold, using tempering)
n_algorithm <- 3  # 3
beta_range <- seq(from = 0.5, to = 3, by = 0.25)
deliv_dir_temp <- file.path("deliv", "numerical_experiment", "potts",
                       "gaussian_cold", "temp")

deliv_dir_ag <- file.path("deliv", "numerical_experiment", "potts",
                          "gaussian_cold", "regAG")

eff_data <- c()
for (i in 1:length(beta_range)) {
  eff_file_temp <- read.csv(file = file.path(deliv_dir_temp,
    paste0("128gaussian_0_", beta_range[i], "_summary_lr.csv")))

  eff_file_ag <- read.csv(file = file.path(deliv_dir_ag,
     paste0("128gaussian_0_", beta_range[i], "_summary.csv")))
  
  # eff_file_lr <- read.csv(file = file.path(deliv_dir,
  #   paste0("128gaussian_0_", beta_range[i], "_summary_lr.csv")))

  eff_data <- c(eff_data, eff_file_temp$eff_temp_AG[1],
                eff_file_temp$eff_temp_HB[1],
                eff_file_ag$eff_AG[1])
}

beta <- rep(beta_range, each = n_algorithm)
# parm <- rep(c("K(x)"), length(beta_range))
# parm_labels <- c(TeX("$K(x)$"))
# parm <- factor(parm, label = parm_labels)
# method <- rep(c("AG (temp)", "Heat Bath (temp)", "AG",
#                 "AG (low, temp)", "AG (low)"), length(beta_range))
# method <- factor(method, level = c("AG", "AG (temp)", "Heat Bath (temp)",
#                                    "AG (low, temp)", "AG (low)"))
method <- rep(c("AG (temp)", "Heat Bath (temp)", "AG"))
method <- factor(method, level = c("AG", "AG (temp)", "Heat Bath (temp)"))

p <- ggplot(data.frame(eff = eff_data, beta = beta,
                       method = method),
            aes(x = beta, y = eff, linetype = method,
                color = method)) +
  geom_line(size = 1.5) + theme_bw() +
  # facet_wrap(~parm, labeller = "label_parsed", ) + 
  scale_y_continuous(trans='log2') +  #, labels = trans_format("log10", math_format(10^.x))) +
  xlab(TeX("$\\beta$")) +
  ylab(TeX("$m_{eff} / s")) +
  theme(text = element_text(size = 13)) +
  scale_color_manual(values=c("#999999", "#D55E00", "#CC79A7", "#F0E442", "#0072B2")) +
  theme(legend.position = c(0.75, 0.75)) +
  theme(legend.title = element_blank()) +
  geom_point(aes(x = 3, y = eff_data[33]), shape = 4, size = 4, colour = "red") +
  geom_point(aes(x = 2.75, y = eff_data[30]), shape = 4, size = 4, colour = "red") +
  geom_point(aes(x = 2.5, y = eff_data[27]), shape = 4, size = 4, colour = "red")
p
