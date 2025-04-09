# This file generates the figures that visualize the outputs from "inference on minimizer for l2svm.R"
library(PropCIs)
library(gridExtra)
library(ggplot2)

sgd_sub_results <- read.csv("sgd_sub_rand_unbalanced_a0.55-0.75_eta0.4_n_iter40000_m10000_K50_b600_lambda0.csv")
sgd_plug_results <- read.csv("sgd_plug_in_unbalanced_a0.55-0.75_eta0.4_n_iter40000_m10000_K50_b600_lambda0.csv")
momentum_sub_results <- read.csv("momentum_sub_rand_unbalanced_a0.55-0.75_eta0.4_momentum_coef1_n_iter40000_m10000_K50_b600_lambda0.csv")
momentum_plug_results <- read.csv("momentum_plug_in_unbalanced_a0.55-0.75_eta0.4_momentum_coef1_n_iter40000_m10000_K50_b600_lambda0.csv")


# Initialize coverage results for 
grid_a = seq(0.55, 0.75, 0.05)
sim = 500
m = 10000  
K = 50
b = 600
n_iter = m+K*b
eta = 0.4
momentum_coef = 1


calculate_metrics <- function(conf, model_weight, sim) {
  # Split the confidence intervals
  even <- seq(2, ncol(conf), 2)
  odd <- seq(1, ncol(conf) - 1, 2)
  right <- conf[, even]
  left <- conf[, odd]
  
  # Calculate coverage ratio
  accept <- (left < model_weight) & (model_weight < right)
  cov0 <- apply(accept, 2, sum)
  coverage <- cov0 / sim
  
  # Calculate interval lengths
  lengths <- right - left
  mean_lengths <- colMeans(lengths)
  sd_lengths <- apply(lengths, 2, sd)
  # Function to calculate Clopper-Pearson interval
  ci <- function(x) {
    exactci(x, sim, 0.95)$conf.int[1:2]
  }
  
  # Calculate Clopper-Pearson intervals
  int <- rep(0,  ncol(conf))
  for (j in 1:(ncol(conf)/2)) {
    int[(2*j-1):(2*j)] <- round(ci(cov0[j]), digits = 3)
  }
  
  # Return a list of results
  list(coverage = coverage, mean_lengths = mean_lengths, sd_lengths = sd_lengths, clopper_pearson_intervals = int)
}



# from simulation with n = 500000, and GD with 1000 step, and unbalanced class 0.2, 0.8
model_deterministic<-NULL
model_deterministic$w<-c(0.603725412,  0.599636710, -0.001170211, -0.603290833, -0.601349460)


# Calculate metrics for both SGD, Momentum, and Plug-in
sgd_metrics <- calculate_metrics(sgd_sub_results, model_deterministic$w[1], sim)
momentum_metrics <- calculate_metrics(momentum_sub_results, model_deterministic$w[1], sim)
sgd_plug_metrics <- calculate_metrics(sgd_plug_results, model_deterministic$w[1], sim)
momentum_plug_metrics <- calculate_metrics(momentum_plug_results, model_deterministic$w[1], sim)

# Create data frame for all methods
plot_data <- data.frame(
  a = rep(grid_a, 4),
  method = factor(rep(c("SGD Sub-randomization", "Momentum Sub-randomization", "SGD Plug-in", "Momentum Plug-in"), each = length(grid_a))),
  coverage = c(sgd_metrics$coverage, momentum_metrics$coverage, sgd_plug_metrics$coverage, momentum_plug_metrics$coverage),
  lower_coverage = c(
    sgd_metrics$clopper_pearson_intervals[seq(1, length(sgd_metrics$clopper_pearson_intervals), 2)],
    momentum_metrics$clopper_pearson_intervals[seq(1, length(momentum_metrics$clopper_pearson_intervals), 2)],
    sgd_plug_metrics$clopper_pearson_intervals[seq(1, length(sgd_plug_metrics$clopper_pearson_intervals), 2)],
    momentum_plug_metrics$clopper_pearson_intervals[seq(1, length(momentum_plug_metrics$clopper_pearson_intervals), 2)]
  ),
  upper_coverage = c(
    sgd_metrics$clopper_pearson_intervals[seq(2, length(sgd_metrics$clopper_pearson_intervals), 2)],
    momentum_metrics$clopper_pearson_intervals[seq(2, length(momentum_metrics$clopper_pearson_intervals), 2)],
    sgd_plug_metrics$clopper_pearson_intervals[seq(2, length(sgd_plug_metrics$clopper_pearson_intervals), 2)],
    momentum_plug_metrics$clopper_pearson_intervals[seq(2, length(momentum_plug_metrics$clopper_pearson_intervals), 2)]
  ),
  length = c(sgd_metrics$mean_lengths, momentum_metrics$mean_lengths, sgd_plug_metrics$mean_lengths, momentum_plug_metrics$mean_lengths),
  length_lower = c(
    sgd_metrics$mean_lengths - sgd_metrics$sd_lengths,
    momentum_metrics$mean_lengths - momentum_metrics$sd_lengths,
    sgd_plug_metrics$mean_lengths - sgd_plug_metrics$sd_lengths,
    momentum_plug_metrics$mean_lengths - momentum_plug_metrics$sd_lengths
  ),
  length_upper = c(
    sgd_metrics$mean_lengths + sgd_metrics$sd_lengths,
    momentum_metrics$mean_lengths + momentum_metrics$sd_lengths,
    sgd_plug_metrics$mean_lengths + sgd_plug_metrics$sd_lengths,
    momentum_plug_metrics$mean_lengths + momentum_plug_metrics$sd_lengths
  )
)

# Create the coverage plot with all methods
coverage_plot <- ggplot(plot_data, aes(x = a, y = coverage, color = method, linetype = method, shape = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_coverage, ymax = upper_coverage), width = 0.02, size = 1) +
  geom_hline(yintercept = 0.9, linetype = "longdash", color = "black", size = 1.2) +
  labs(
    x = "a",
    y = "Coverage",
  ) +
  theme_minimal(base_size = 18) +
  ylim(0.5, 1) +
  scale_color_manual(values = c("SGD Sub-randomization" = "blue", "Momentum Sub-randomization" = "darkgreen", "SGD Plug-in" = "orange", "Momentum Plug-in" = "purple")) +
  scale_linetype_manual(values = c("SGD Sub-randomization" = "twodash", "Momentum Sub-randomization" = "longdash", 
                                   "SGD Plug-in" = "dashed", "Momentum Plug-in" = "solid")) +
  scale_shape_manual(values = c("SGD Sub-randomization" = 3, "Momentum Sub-randomization" = 15, 
                                "SGD Plug-in" = 17, "Momentum Plug-in" = 16)) +
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text = element_text(size = 18))

coverage_plot
# Create the length plot with all methods
length_plot <- ggplot(plot_data, aes(x = a, y = length, color = method, linetype = method, shape = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = length_lower, ymax = length_upper), width = 0.02, size = 1) +
  labs(
    x = "a",
    y = "Confidence Interval Length",
  ) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("SGD Sub-randomization" = "blue", "Momentum Sub-randomization" = "darkgreen", "SGD Plug-in" = "orange", "Momentum Plug-in" = "purple")) +
  scale_linetype_manual(values = c("SGD Sub-randomization" = "twodash", "Momentum Sub-randomization" = "longdash", 
                                   "SGD Plug-in" = "dashed", "Momentum Plug-in" = "solid")) +
  scale_shape_manual(values = c("SGD Sub-randomization" = 3, "Momentum Sub-randomization" = 15, 
                                "SGD Plug-in" = 17, "Momentum Plug-in" = 16)) +
  coord_cartesian(ylim = c(0, max(plot_data$length) * 1.1)) +
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text = element_text(size = 18))

length_plot

# Function to extract legend from a ggplot
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


plot_with_legend <- coverage_plot + 
  theme(legend.position = "right",legend.title = element_blank(),
        legend.text = element_text(size = 14))


legend <- get_legend(plot_with_legend)
combined_plot <- grid.arrange(
  coverage_plot + theme(legend.position = "none"),
  length_plot + theme(legend.position = "none"),
  legend,
  ncol = 3,  # Set to 2 columns for side-by-side arrangement
  widths = c(1, 1, 0.5)  # Equal widths for both plots
)

final_file_name <- paste0("final_figure",min(grid_a), "-", max(grid_a), 
                          "_eta", eta, "_momentum_coef", momentum_coef, "_n_iter", n_iter, "_m", m, "_K", K, "_b", b, "_lambda0_", "sim", sim, ".pdf")
ggsave(final_file_name, 
       plot = combined_plot, 
       width = 16,  # Increased width to accommodate side-by-side layout
       height = 6)

