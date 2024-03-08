# Bethany Allen (updated 19th Feb 2024)
# Adapted from code by Maria Volkova
# Code to plot results of coalescent analysis

library(tidyverse)
library(coda)
library(ggthemes)

all_trees <- data.frame(); all_trees_pop <- data.frame()

trees <- c("Benson1", "Benson2", "Lloyd1", "Lloyd2")

for (i in 1:4){
  #Create path string
  log_path <- paste0("Coalescent_analysis/results/PiecewiseCoalescent_8.", trees[i], "clean.log")
  
  #Read in coalescent log and trim 10% burn-in
  log_table <- read.table(log_path, header = T) %>% slice_tail(prop = 0.9)
  log_mcmc <- as.mcmc(log_table)
  
  #Calculate median and 95% HPD values
  summary_data <- as.data.frame(HPDinterval(log_mcmc))
  summary_data$median <- apply(log_table, 2, median)
  
  #Create function to calculate proportion of posterior above 0
  prop.positive <- function(x) length(which(x > 0))/length(x)
  summary_data$prop_positive <- apply(log_table, 2, prop.positive)
  
  diversification_data <- summary_data[6:13,]

  #Add interval names
  diversification_data$interval <- c("Coniacian- Maastrichtian", "Aptian- Turonian",
                                     "Berriasian- Barremian", "Late Jurassic",
                                     "Middle Jurassic", "Early Jurassic",
                                     "Late Triassic", "Early-Mid Triassic")
  
  #Add tree number
  diversification_data$tree <- trees[i]

  #Bind to overall results
  all_trees <- rbind(all_trees, diversification_data)
  
  #Extract origin data
  pop_data <- summary_data[5,]
  
  #Add tree number
  pop_data[4] <- trees[i]
  
  #Bind to overall results
  all_trees_pop <- rbind(all_trees_pop, pop_data)
}

all_trees$interval <- factor(all_trees$interval,
                             levels = c("Early-Mid Triassic", "Late Triassic",
                                        "Early Jurassic", "Middle Jurassic",
                                        "Late Jurassic", "Berriasian- Barremian",
                                        "Aptian- Turonian", "Coniacian- Maastrichtian"))

all_trees$tree <- factor(all_trees$tree,
                             levels = c("Lloyd1", "Benson1", "Benson2",
                                        "Lloyd2"))

## Error plot
ggplot(data = all_trees, aes(x = interval, y = median, ymin = lower,
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
  labs(x = "Interval", y = "Diversification rate") +
  theme_classic(base_size = 17)
ggsave("figures/Coalescent_all.pdf", scale = 3.5, dpi = 600)

## Results tables
divers_table <- cbind(tree = as.character(all_trees$tree),
                      interval = as.character(all_trees$interval),
                      lower_HPD = signif(all_trees$lower, 3),
                      median = signif(all_trees$median, 3),
                      upper_HPD = signif(all_trees$upper, 3),
                      prop_positive = signif(all_trees$prop_positive, 3))
write.csv(divers_table, "tables/Coalescent_diversification.csv",
          row.names = FALSE)

startpop_table <- cbind(tree = all_trees_pop$V4,
                      lower_HPD = signif(all_trees_pop$lower, 5),
                      median = signif(all_trees_pop$median, 5),
                      upper_HPD = signif(all_trees_pop$upper, 5))
write.csv(startpop_table, "tables/Coalescent_popsize.csv",
          row.names = FALSE)
