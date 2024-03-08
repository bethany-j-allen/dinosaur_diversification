# Bethany Allen (updated 19th Feb 2024)
# Adapted from code by Maria Volkova
# Code to plot results of BDSKY analysis

library(tidyverse)
library(coda)
library(ggthemes)

all_trees_birth <- data.frame(); all_trees_death <- data.frame()
all_trees_div <- data.frame(); all_trees_turn <- data.frame()
all_trees_samp <- data.frame(); all_trees_ori <- data.frame()

trees <- c("Benson1", "Benson2", "Lloyd1", "Lloyd2")

for (i in 1:4){
  #Create path string
  log_path <- paste0("BDSKY_analysis/results/BDSKY_8.", trees[i], "clean.log")
  
  #Read in BDSKY log and trim 10% burn-in
  data <- read.table(log_path, header = T) %>% slice_tail(prop = 0.9)
  
  #Calculate diversification and turnover
  birth_rates <- select(data, starts_with("birthRate"))
  death_rates <- select(data, starts_with("deathRate"))
  
  div_rates <- birth_rates - death_rates
  colnames(div_rates) <- paste0("divRate.",
                                      seq(1:ncol(div_rates)))
  
  TO_rates <- birth_rates + death_rates
  colnames(TO_rates) <- paste0("TORate.",
                                      seq(1:ncol(TO_rates)))
  
  #Create function to calculate proportion of posterior above 0
  prop.positive <- function(x) length(which(x > 0))/length(x)
  
  #Calculate median and 95% HPD values for births
  birth_mcmc <- as.mcmc(birth_rates)
  birth_data <- as.data.frame(HPDinterval(birth_mcmc))
  birth_data$median <- apply(birth_rates, 2, median)
  
  #Add interval names
  birth_data$interval <- c("Early-Mid Triassic", "Late Triassic",
                           "Early Jurassic", "Middle Jurassic",
                           "Late Jurassic", "Berriasian- Barremian",
                           "Aptian- Turonian", "Coniacian- Maastrichtian")
  #Add tree number
  birth_data$tree <- trees[i]
  
  #Bind to overall results
  all_trees_birth <- rbind(all_trees_birth, birth_data)
  
  #Calculate median and 95% HPD values for deaths
  death_mcmc <- as.mcmc(death_rates)
  death_data <- as.data.frame(HPDinterval(death_mcmc))
  death_data$median <- apply(death_rates, 2, median)
  
  #Add interval names
  death_data$interval <- c("Early-Mid Triassic", "Late Triassic",
                           "Early Jurassic", "Middle Jurassic",
                           "Late Jurassic", "Berriasian- Barremian",
                           "Aptian- Turonian", "Coniacian- Maastrichtian")
  #Add tree number
  death_data$tree <- trees[i]
  
  #Bind to overall results
  all_trees_death <- rbind(all_trees_death, death_data)
  
  #Calculate median and 95% HPD values for diversification
  div_mcmc <- as.mcmc(div_rates)
  div_data <- as.data.frame(HPDinterval(div_mcmc))
  div_data$median <- apply(div_rates, 2, median)
  div_data$prop_positive <- apply(div_rates, 2, prop.positive)
  
  #Add interval names
  div_data$interval <- c("Early-Mid Triassic", "Late Triassic",
                         "Early Jurassic", "Middle Jurassic",
                         "Late Jurassic", "Berriasian- Barremian",
                         "Aptian- Turonian", "Coniacian- Maastrichtian")
  #Add tree number
  div_data$tree <- trees[i]
  
  #Bind to overall results
  all_trees_div <- rbind(all_trees_div, div_data)
  
  #Calculate median and 95% HPD values for turnover
  turn_mcmc <- as.mcmc(TO_rates)
  turn_data <- as.data.frame(HPDinterval(turn_mcmc))
  turn_data$median <- apply(TO_rates, 2, median)
  
  #Add interval names
  turn_data$interval <- c("Early-Mid Triassic", "Late Triassic",
                          "Early Jurassic", "Middle Jurassic",
                          "Late Jurassic", "Berriasian- Barremian",
                          "Aptian- Turonian", "Coniacian- Maastrichtian")
  #Add tree number
  turn_data$tree <- trees[i]
  
  #Bind to overall results
  all_trees_turn <- rbind(all_trees_turn, turn_data)
  
  #Calculate median and 95% HPD values for sampling
  samp_log <- select(data, starts_with("samplingRate"))
  samp_mcmc <- as.mcmc(samp_log)
  samp_data <- as.data.frame(HPDinterval(samp_mcmc))
  samp_data$median <- apply(samp_log, 2, median)
  
  #Add interval names
  samp_data$interval <- c("Early-Mid Triassic", "Late Triassic",
                          "Early Jurassic", "Middle Jurassic",
                          "Late Jurassic", "Berriasian- Barremian",
                          "Aptian- Turonian", "Coniacian- Maastrichtian")
  #Add tree number
  samp_data$tree <- trees[i]

  #Bind to overall results
  all_trees_samp <- rbind(all_trees_samp, samp_data)
  
  #Extract origin data
  origin_mcmc <- as.mcmc(pull(data, "originFBD"))
  origin_data <- as.data.frame(HPDinterval(origin_mcmc))
  origin_data$median <- median(data$originFBD)
  origin_data <- origin_data + 66

  #Add tree number
  origin_data[4] <- trees[i]
  
  #Bind to overall results
  all_trees_ori <- rbind(all_trees_ori, origin_data)
}

all_trees_birth$interval <- factor(all_trees_birth$interval,
                                 levels = c("Early-Mid Triassic", "Late Triassic",
                                            "Early Jurassic", "Middle Jurassic",
                                            "Late Jurassic", "Berriasian- Barremian",
                                            "Aptian- Turonian", "Coniacian- Maastrichtian"))

all_trees_birth$tree <- factor(all_trees_birth$tree,
                         levels = c("Lloyd1", "Benson1", "Benson2",
                                    "Lloyd2"))

all_trees_death$interval <- factor(all_trees_death$interval,
                                 levels = c("Early-Mid Triassic", "Late Triassic",
                                            "Early Jurassic", "Middle Jurassic",
                                            "Late Jurassic", "Berriasian- Barremian",
                                            "Aptian- Turonian", "Coniacian- Maastrichtian"))

all_trees_death$tree <- factor(all_trees_death$tree,
                         levels = c("Lloyd1", "Benson1", "Benson2",
                                    "Lloyd2"))

all_trees_div$interval <- factor(all_trees_div$interval,
                                 levels = c("Early-Mid Triassic", "Late Triassic",
                                            "Early Jurassic", "Middle Jurassic",
                                            "Late Jurassic", "Berriasian- Barremian",
                                            "Aptian- Turonian", "Coniacian- Maastrichtian"))

all_trees_div$tree <- factor(all_trees_div$tree,
                         levels = c("Lloyd1", "Benson1", "Benson2",
                                    "Lloyd2"))

all_trees_turn$interval <- factor(all_trees_turn$interval,
                                  levels = c("Early-Mid Triassic", "Late Triassic",
                                             "Early Jurassic", "Middle Jurassic",
                                             "Late Jurassic", "Berriasian- Barremian",
                                             "Aptian- Turonian", "Coniacian- Maastrichtian"))

all_trees_turn$tree <- factor(all_trees_turn$tree,
                         levels = c("Lloyd1", "Benson1", "Benson2",
                                    "Lloyd2"))

all_trees_samp$interval <- factor(all_trees_samp$interval,
                                  levels = c("Early-Mid Triassic", "Late Triassic",
                                             "Early Jurassic", "Middle Jurassic",
                                             "Late Jurassic", "Berriasian- Barremian",
                                             "Aptian- Turonian", "Coniacian- Maastrichtian"))

all_trees_samp$tree <- factor(all_trees_samp$tree,
                         levels = c("Lloyd1", "Benson1", "Benson2",
                                    "Lloyd2"))

## Error plots
ggplot(data = all_trees_birth, aes(x = interval, y = median, ymin = lower,
                                 ymax = upper, group = tree, col = tree)) +
  annotate(geom = "rect", xmin = c(1.5, 3.5, 5.5, 7.5),
           xmax = c(2.5, 4.5, 6.5, 8.5),
           ymin = -Inf, ymax = Inf,
           colour = "grey", alpha = 0.1, linewidth = 0)  +
  geom_point(size = 1.5, position = position_dodge(1)) +
  geom_errorbar(linewidth = 1, width = 0.8, position = position_dodge(1)) +
  scale_colour_colorblind() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Speciation rate") +
  theme_classic(base_size = 17)
ggsave("figures/BDSKY_all_birth.pdf", scale = 3.5, dpi = 600)

ggplot(data = all_trees_death, aes(x = interval, y = median, ymin = lower,
                                   ymax = upper, group = tree, col = tree)) +
  annotate(geom = "rect", xmin = c(1.5, 3.5, 5.5, 7.5),
           xmax = c(2.5, 4.5, 6.5, 8.5),
           ymin = -Inf, ymax = Inf,
           colour = "grey", alpha = 0.1, linewidth = 0)  +
  geom_point(size = 1.5, position = position_dodge(1)) +
  geom_errorbar(linewidth = 1, width = 0.8, position = position_dodge(1)) +
  scale_colour_colorblind() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Extinction rate") +
  theme_classic(base_size = 17)
ggsave("figures/BDSKY_all_death.pdf", scale = 3.5, dpi = 600)

ggplot(data = all_trees_div, aes(x = interval, y = median, ymin = lower,
                             ymax = upper, group = tree, col = tree)) +
  annotate(geom = "rect", xmin = c(1.5, 3.5, 5.5, 7.5),
           xmax = c(2.5, 4.5, 6.5, 8.5),
           ymin = -Inf, ymax = Inf,
           colour = "grey", alpha = 0.1, linewidth = 0)  +
  geom_point(size = 1.5, position = position_dodge(1)) +
  geom_errorbar(linewidth = 1, width = 0.8, position = position_dodge(1)) +
  scale_colour_colorblind() +
  geom_hline(aes(yintercept = 0), colour = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  #scale_y_continuous(limits = c(-0.1, 0.3)) +
  labs(x = "Interval", y = "Diversification rate") +
  theme_classic(base_size = 17)
ggsave("figures/BDSKY_all_div.pdf", scale = 3.5, dpi = 600)

ggplot(data = all_trees_turn, aes(x = interval, y = median, ymin = lower,
                             ymax = upper, group = tree, col = tree)) +
  annotate(geom = "rect", xmin = c(1.5, 3.5, 5.5, 7.5),
           xmax = c(2.5, 4.5, 6.5, 8.5),
           ymin = -Inf, ymax = Inf,
           colour = "grey", alpha = 0.1, linewidth = 0)  +
  geom_point(size = 1.5, position = position_dodge(1)) +
  geom_errorbar(linewidth = 1, width = 0.8, position = position_dodge(1)) +
  scale_colour_colorblind() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Turnover rate") +
  theme_classic(base_size = 17)
ggsave("figures/BDSKY_all_turn.pdf", scale = 3.5, dpi = 600)

ggplot(data = all_trees_samp, aes(x = interval, y = median, ymin = lower,
                                  ymax = upper, group = tree, col = tree)) +
  annotate(geom = "rect", xmin = c(1.5, 3.5, 5.5, 7.5),
           xmax = c(2.5, 4.5, 6.5, 8.5),
           ymin = -Inf, ymax = Inf,
           colour = "grey", alpha = 0.1, linewidth = 0)  +
  geom_point(size = 1.5, position = position_dodge(1)) +
  geom_errorbar(linewidth = 1, width = 0.8, position = position_dodge(1)) +
  scale_colour_colorblind() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  #scale_y_continuous(limits = c(0, 0.015)) +
  labs(x = "Interval", y = "Sampling rate") +
  theme_classic(base_size = 17)
ggsave("figures/BDSKY_all_samp.pdf", scale = 3.5, dpi = 600)

all_trees_ori$median <- as.numeric(all_trees_ori$median)
all_trees_ori$lower <- as.numeric(all_trees_ori$lower)
all_trees_ori$upper <- as.numeric(all_trees_ori$upper)
all_trees_ori$V4 <- factor(all_trees_ori$V4,
                             levels = c("Lloyd1", "Benson1", "Benson2",
                                        "Lloyd2"))

ggplot(data = all_trees_ori, aes(x = V4, y = median, ymin = lower,
                                  ymax = upper, group = V4, col = V4)) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf,
           ymin = 251.9, ymax = Inf,
           colour = "grey", alpha = 0.1, linewidth = 0)  +
  geom_point(size = 3) +
  geom_errorbar(linewidth = 1, width = 0.8) +
  scale_colour_colorblind() +
  labs(x = "Interval", y = "Origin") +
  theme_classic(base_size = 14)
ggsave("figures/BDSKY_all_ori.pdf", dpi = 600)

## Results tables
birth_table <- cbind(tree = as.character(all_trees_birth$tree),
                      interval = as.character(all_trees_birth$interval),
                      lower_HPD = signif(all_trees_birth$lower, 3),
                      median = signif(all_trees_birth$median, 3),
                      upper_HPD = signif(all_trees_birth$upper, 3))
write.csv(birth_table, "tables/FBD_speciation.csv",
          row.names = FALSE)

death_table <- cbind(tree = as.character(all_trees_death$tree),
                     interval = as.character(all_trees_death$interval),
                     lower_HPD = signif(all_trees_death$lower, 3),
                     median = signif(all_trees_death$median, 3),
                     upper_HPD = signif(all_trees_death$upper, 3))
write.csv(death_table, "tables/FBD_extinction.csv",
          row.names = FALSE)

div_table <- cbind(tree = as.character(all_trees_div$tree),
                     interval = as.character(all_trees_div$interval),
                     lower_HPD = signif(all_trees_div$lower, 3),
                     median = signif(all_trees_div$median, 3),
                     upper_HPD = signif(all_trees_div$upper, 3),
                     prop_positive = signif(all_trees_div$prop_positive, 3))
write.csv(div_table, "tables/FBD_diversification.csv",
          row.names = FALSE)

samp_table <- cbind(tree = as.character(all_trees_samp$tree),
                   interval = as.character(all_trees_samp$interval),
                   lower_HPD = signif(all_trees_samp$lower, 3),
                   median = signif(all_trees_samp$median, 3),
                   upper_HPD = signif(all_trees_samp$upper, 3))
write.csv(samp_table, "tables/FBD_sampling.csv",
          row.names = FALSE)

origin_table <- cbind(tree = as.character(all_trees_ori$V4),
                        lower_HPD = signif(all_trees_ori$lower, 5),
                        median = signif(all_trees_ori$median, 5),
                        upper_HPD = signif(all_trees_ori$upper, 5))
write.csv(origin_table, "tables/FBD_origin.csv",
          row.names = FALSE)
