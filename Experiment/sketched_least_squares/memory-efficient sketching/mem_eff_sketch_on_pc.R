# Required Libraries
library(phangorn) # for fhm function
library(data.table)
library(ggplot2)
library(gridExtra)


# SRHT Function
# Implements the Subsampled Randomized Hadamard Transform.
SRHT <- function(select, D, a) {
  Da <- D*a
  Sa <- fhm(Da)[select]  
  return(Sa)
}


# Large Matrix Generation
set.seed(123)
n <- 2^12  # number of rows
p <- 100 # number of columns
X <- matrix(rnorm(n * p), n, p)
grid_b = 200
write.csv(X, sprintf("data_matrix_n%d_p%d.csv", n, p), row.names = FALSE)  # CSV format


# Performance Measurement Function
measure_loading_performance <- function(method, file_path, n, p, b, K) {
  gc()
  start_time_total <- Sys.time()
  loading_time <- 0
  
  select_list <- lapply(1:K, function(k) which(rbinom(n, 1, b/n) != 0))
  D_list <- lapply(1:K, function(k) sample(c(1,-1), n, replace=TRUE))
  sketches_b <- vector("list", K)
  for(k in 1:K) {
    sketches_b[[k]] <- matrix(0, length(select_list[[k]]), p)
  }
  
  if (method == "full") {
    for (k in 1:K) {
      start_load <- Sys.time()
      X <- as.matrix(fread(file_path))
      end_load <- Sys.time()
      loading_time <- loading_time + as.numeric(difftime(end_load, start_load, units = "secs"))
      for(j in 1:p) {
        sketches_b[[k]][,j] <- SRHT(select_list[[k]], D_list[[k]], X[,j])
      }
      
      rm(X)
      gc()
    }
  } else if (method == "block") {
    for(k in 1:K) {
      sketches_b[[k]] <- matrix(0, length(select_list[[k]]), p)
    }
    
    num_blocks <- 10
    block_size <- p/num_blocks
    
    for(i in 1:num_blocks) {
      start_col <- ((i-1) * block_size + 1)
      end_col <- i * block_size
      cols <- start_col:end_col
      
      start_load <- Sys.time()
      block_data <- fread(file_path, 
                          select = cols,
                          header = TRUE,
                          data.table = FALSE)
      end_load <- Sys.time()
      loading_time <- loading_time + as.numeric(difftime(end_load, start_load, units = "secs"))
      #print(loading_time)
      for(j in 1:ncol(block_data)) {
        col_idx <- start_col + j - 1
        for(k in 1:K) {
          sketches_b[[k]][, col_idx] <- SRHT(select_list[[k]], 
                                             D_list[[k]], 
                                             block_data[,j])
        }
      }
      
      rm(block_data)
      gc()
    }
  }
  
  end_time_total <- Sys.time() 
  return(list(
    total_time = difftime(end_time_total, start_time_total, units = "secs"),
    loading_time = loading_time
  ))
}

# Performance Testing
n_reps <- 5
K_range <- seq(20, 100, by = 20)

loading_results <- data.frame(
  method = character(),
  K = numeric(),
  mean_total_time = numeric(),
  sd_total_time = numeric(),
  mean_loading_time = numeric(),
  sd_loading_time = numeric()
)

for (method in c("full", "block")) {
  for (K in K_range) {
    total_times <- numeric(n_reps)
    loading_times <- numeric(n_reps)
    
    cat(sprintf("\nTesting %s loading method with K=%d...\n", method, K))
    
    for (i in 1:n_reps) {
      cat(sprintf("Replication %d/%d\n", i, n_reps))
      perf <- measure_loading_performance(method, sprintf("data_matrix_n%d_p%d.csv", n, p), n, p, grid_b[1], K)
      total_times[i] <- as.numeric(perf$total_time)
      loading_times[i] <- as.numeric(perf$loading_time)
    }
    
    loading_results <- rbind(loading_results, data.frame(
      method = method,
      K = K,
      mean_total_time = mean(total_times),
      sd_total_time = sd(total_times),
      mean_loading_time = mean(loading_times),
      sd_loading_time = sd(loading_times)
    ))
  }
}

# Results Visualization
time_results_long <- rbind(
  data.frame(
    method = loading_results$method,
    K = loading_results$K,
    time = loading_results$mean_loading_time,
    sd_time = loading_results$sd_loading_time,
    type = "Loading Time"
  ),
  data.frame(
    method = loading_results$method,
    K = loading_results$K,
    time = loading_results$mean_total_time - loading_results$mean_loading_time,
    sd_time = sqrt(loading_results$sd_total_time^2 + loading_results$sd_loading_time^2),
    type = "Arithmetic Time"
  ),
  data.frame(
    method = loading_results$method,
    K = loading_results$K,
    time = loading_results$mean_total_time,
    sd_time = loading_results$sd_total_time,
    type = "Total Time"
  )
)


write.csv(loading_results, 
          sprintf("loading_results_n%d_p%d_b%d.csv", n, p, grid_b[1]), 
          row.names = FALSE)

# Detailed time results for different types
write.csv(time_results_long,
          sprintf("time_results_detailed_n%d_p%d_b%d.csv", n, p, grid_b[1]),
          row.names = FALSE)


loading_results <- read.csv(sprintf("loading_results_n%d_p%d_b%d.csv", n, p, grid_b[1]))
time_results_long <- read.csv(sprintf("time_results_detailed_n%d_p%d_b%d.csv", n, p, grid_b[1]))


# Loading Time visualization
p1 <- ggplot(subset(time_results_long, type == "Loading Time"), 
             aes(x = K, y = time, color = method, group = method, 
                 linetype = method, shape = method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = time - sd_time, 
                    ymax = time + sd_time),
                width = 2, linewidth = 1.5) +
  labs(title = "Loading Time by Method and K",
       x = "Number of Sketches (K)",
       y = "Time (seconds)") +
  scale_shape_manual(values = c(16, 17)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold")
  )

# Arithmetic Time visualization
p2 <- ggplot(subset(time_results_long, type == "Arithmetic Time"), 
             aes(x = K, y = time, color = method, group = method, 
                 linetype = method, shape = method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = time - sd_time, 
                    ymax = time + sd_time),
                width = 2, linewidth = 1.5) +
  labs(title = "Arithmetic Time by Method and K",
       x = "Number of Sketches (K)",
       y = "Time (seconds)") +
  scale_shape_manual(values = c(16, 17)) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 18)
  )

# Total Time visualization
p3 <- ggplot(subset(time_results_long, type == "Total Time"), 
             aes(x = K, y = time, color = method, group = method, 
                 linetype = method, shape = method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = time - sd_time, 
                    ymax = time + sd_time),
                width = 2, linewidth = 1.5) +
  labs(title = "Total Time by Method and K",
       x = "Number of Sketches (K)",
       y = "Time (seconds)") +
  scale_shape_manual(values = c(16, 17)) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 18)
  )

# Arrange all plots in a 2x2 grid
grid.arrange(p1,  p2, p3, ncol = 3)

# Print detailed results
print("\nDetailed Loading Performance Results (Mean ± SD):")
print(loading_results)


ggsave("loading_time.pdf", p1, width = 10,  height = 8,  dpi = 300)

ggsave("arithmetic_time.pdf", p2, width = 8, height = 6,  dpi = 300)

ggsave("total_time.pdf", p3, width = 8, height = 6, dpi = 300)




p_combined <- ggplot(subset(time_results_long, type != "Arithmetic Time"), 
                     aes(x = K, y = time, color = method, group = method, 
                         linetype = method, shape = method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = time - sd_time, 
                    ymax = time + sd_time),
                width = 2, linewidth = 1.5) +
  facet_wrap(~factor(type, levels = c("Loading Time", "Total Time")), 
             scales = "free_y", ncol = 2) +
  labs(x = "Number of Sketches (K)",
       y = "Time (seconds)") +
  scale_shape_manual(values = c(16, 17)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.key.size = unit(2, "cm"),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

ggsave(sprintf("combined_time_analysis_n%d_p%d_b%d.pdf", 
               n, p, grid_b[1]),
       p_combined,
       width = 12,
       height = 6,
       dpi = 300)


