# Experiment for inference in least squares support vector machine via vanilla SGD and stochastic heavy ball
library(MASS)

# Define the mixture distribution
mean_1 <- c(1, 1, 1, 0, 0)    
mean_2 <- c(0, 0, 1, 1, 1) 
mixture_weights <- c(0.2, 0.8) 
Sigma <- diag(rep(0.5, 5))    
lambda0 <- 0


# Function to generate a single sample from the mixture
generate_single_sample <- function() {
  component <- sample(1:2, 1, prob = mixture_weights)
  
  if (component == 1) {
    x <- mvrnorm(1, mu = mean_1, Sigma = Sigma)
    y <- 1
  } else {
    x <- mvrnorm(1, mu = mean_2, Sigma = Sigma)
    y <- -1
  }
  
  return(list(x = x, y = y))
}


# Obtain the minimizer of population loss via Gradient Descent
n_samples <- 500000
X <- matrix(0,n_samples,5)
y <- numeric(n_samples)
for(i in 1:n_samples){
  new_sample<-generate_single_sample()
  X[i,]<-new_sample$x
  y[i]<-new_sample$y
}
  
gd_l2svm <- function(X, y, learning_rate, n_iter = 1000, lambda = 0) {
    n_samples <- nrow(X)
    n_features <- ncol(X)
    w<-rep(0,5)
    w0 <- 0
    
    # Store solution path
    w_path <- matrix(0, nrow = n_iter, ncol = n_features)
    w0_path <- numeric(n_iter)
    loss_path <- numeric(n_iter)
    
    for(iter in 1:n_iter) {
      # Compute hinge loss gradient for the entire dataset
      margins <- y * (X %*% w + w0)
      indicator <- as.numeric(ifelse(margins < 1, 1, 0))
      # Compute gradient
      gradient_w <- -colSums((indicator  * pmax(0,1-margins)*y) * X) / n_samples + lambda * w
      gradient_w0 <- -sum(indicator *pmax(0,1-margins)*y) / n_samples
      
      # Update weights and bias
      w <- w - learning_rate * gradient_w
      w0 <- w0 - learning_rate * gradient_w0
      
      # Store current solution
      w_path[iter,] <- w
      w0_path[iter] <- w0
      
      hinge_losses <- cbind(0,1-margins)
      loss_path[iter] <- sum((apply(hinge_losses,1,max))^2)/(2*n_samples)+  (lambda/2) * sum(w^2)
    }
    
    list(w = w, w0 = w0, w_path = w_path, w0_path = w0_path, loss_path = loss_path)
}

lambda0 <- 0
model_deterministic <- gd_l2svm(X, y, learning_rate = 0.1, n_iter = 1000, lambda = lambda0)
model_deterministic$w

# Solution path for sgd/momentum-sgd
sgd_l2svm_online <- function(learning_rates, n_iter = 1000, lambda = 0, 
                             w_start = NULL, w0_start = NULL, burnin = 0, use_momentum = FALSE, momentum = NULL) {
  
  if (is.null(w_start)) {
    w_start <- rep(0, 5)
    w0_start <- 0   
    if (burnin > 0) {
      for (iter in 1:burnin) {
        new_sample <- generate_single_sample()
        x_i <- new_sample$x
        y_i <- new_sample$y
        
        margin <- y_i * (sum(w * x_i) + w0)
        
        indicator <- as.numeric(ifelse(margin < 1, 1, 0))
        # Compute gradient
        gradient_w <- -indicator  * pmax(0,1-margin)*y_i * x_i + lambda * w
        gradient_w0 <- -indicator *pmax(0,1-margin)*y_i
        
        # Standard SGD update
        w_start  <- w_start - current_lr * gradient_w
        w0_start <- w0_start - current_lr * gradient_w0   
      }
    }
  }
  
  w <- w_start
  w0 <- w0_start
  v_w <- rep(0, 5)
  v_w0 <- 0
  
  # Store solution path
  w_path <- matrix(0, nrow = n_iter, ncol = 5)
  w0_path <- numeric(n_iter)
  loss_path <- numeric(n_iter)
  accum_hinge_loss <- 0
  
  for (iter in 1:n_iter) {
    # Get current learning rate from schedule
    current_lr <- learning_rates[iter]
    
    # Randomly select one sample
    new_sample <- generate_single_sample()
    x_i <- new_sample$x
    y_i <- new_sample$y
    
    # Compute hinge loss
    margin <- y_i * (sum(w * x_i) + w0)
    indicator <- as.numeric(ifelse(margin < 1, 1, 0))
    # Compute gradient
    gradient_w <- -indicator  * pmax(0,1-margin)*y_i * x_i + lambda * w
    gradient_w0 <- -indicator *pmax(0,1-margin)*y_i
    if (use_momentum) {
      # Update velocities for momentum  
      v_w <- momentum[iter] * v_w + current_lr * gradient_w
      v_w0 <- momentum[iter] * v_w0 + current_lr * gradient_w0    
      
      # Update parameters
      w <- w - current_lr * v_w
      w0 <- w0 - current_lr * v_w0
    } else {
      # Standard SGD update
      w <- w - current_lr * gradient_w
      w0 <- w0 - current_lr * gradient_w0
    }
    
    
    accum_hinge_loss <- accum_hinge_loss + (max(0, 1 - margin))^2
    loss_path[iter] <- accum_hinge_loss / (2 * iter) + (lambda / 2) * sum(w^2)
    
    w_path[iter, ] <- w
    w0_path[iter] <- w0
  }
  
  list(w = w, w0 = w0, w_path = w_path, w0_path = w0_path, loss_path = loss_path)
}

# CI for sub-randomamization and multi-run plug-in 
l2svm_online_inference_ci<-function(n_iter,alpha, m, K, b, lambda0, eta, a, use_momentum=FALSE, momentum_coef = 1){
  learning_rates <- eta * (1:n_iter)^(-a) 
  momentum <- pmax(1 - momentum_coef * learning_rates, 0) 
  subske <- numeric(K)
  if(use_momentum){
    res <- sgd_l2svm_online(learning_rates, n_iter = m, lambda = lambda0, w_start = NULL, w0_start = NULL, burnin = 0, use_momentum = TRUE, momentum = momentum)
    center <- res$w[1]
    for (i in 1:K){
      res <- sgd_l2svm_online(learning_rates, n_iter = b, lambda = lambda0, w_start = NULL, w0_start = NULL, burnin = 0, use_momentum = TRUE, momentum = momentum)
      subske[i] <- res$w[1]
    }
  }else{
    res <- sgd_l2svm_online(learning_rates, n_iter = m, lambda = lambda0, w_start = NULL, w0_start = NULL, burnin = 0, use_momentum = FALSE)
    center <- res$w[1]
    for (i in 1:K){
      res <- sgd_l2svm_online(learning_rates, n_iter = b, lambda = lambda0, w_start = NULL, w0_start = NULL, burnin = 0, use_momentum = FALSE)
      subske[i] <- res$w[1]
    }
  }
  taum <- m^(a/2)  
  taub <- b^(a/2)
  
  lb <- quantile(taub*(subske-center), alpha/2)/(taum-taub)
  rb <- quantile(taub*(subske-center), (1-alpha/2))/(taum-taub)
  sub_conf <- c(center-rb, center-lb)
  
  
  plug_v <- var(subske)
  plug_rb <- qnorm(1-alpha/2, sd=sqrt(plug_v)*taub/taum)
  plug_lb <- qnorm(alpha/2, sd=sqrt(plug_v)*taub/taum)
  plug_conf<- c(center-plug_rb, center-plug_lb)
  
  
  return(list(sub_conf = sub_conf, plug_conf = plug_conf))
}



sim = 2
grid_a = seq(0.55,0.75,0.05)
sgd_sub_conf_list <- list()
sgd_plug_conf_list <- list()
momentum_sub_conf_list <- list()
momentum_plug_conf_list <- list()
for (a_idx in 1:length(grid_a)) {
  sgd_sub_conf_list[[a_idx]] <- matrix(0, sim, 2)
  sgd_plug_conf_list[[a_idx]] <- matrix(0, sim, 2)
  momentum_sub_conf_list[[a_idx]] <- matrix(0, sim, 2)
  momentum_plug_conf_list[[a_idx]] <- matrix(0, sim, 2)
}



m = 10000
K = 50
b = 600
n_iter = m + K*b
momentum_coef = 1
eta = 0.4

sgd_sub_conf_list <- list()
sgd_plug_conf_list <- list()
momentum_sub_conf_list <- list()
momentum_plug_conf_list <- list()
for (a_idx in 1:length(grid_a)) {
  sgd_sub_conf_list[[a_idx]] <- matrix(0, sim, 2)
  sgd_plug_conf_list[[a_idx]] <- matrix(0, sim, 2)
  momentum_sub_conf_list[[a_idx]] <- matrix(0, sim, 2)
  momentum_plug_conf_list[[a_idx]] <- matrix(0, sim, 2)
}
for(a_idx in 1:length(grid_a)){
  for(i in 1:sim){
    result <- l2svm_online_inference_ci(n_iter = n_iter, alpha = 0.1, m = m, K = K, b = b, lambda0 = 0, eta = eta, a = grid_a[a_idx], use_momentum=FALSE, momentum_coef = momentum_coef)
    sgd_sub_conf_list[[a_idx]][i,] <- result$sub_conf
    sgd_plug_conf_list[[a_idx]][i,] <- result$plug_conf
    
    result <- l2svm_online_inference_ci(n_iter = n_iter, alpha = 0.1, m = m, K = K, b = b, lambda0 = 0, eta = eta, a = grid_a[a_idx], use_momentum=TRUE, momentum_coef = momentum_coef)
    momentum_sub_conf_list[[a_idx]][i,] <- result$sub_conf
    momentum_plug_conf_list[[a_idx]][i,] <- result$plug_conf
  }
}

# Combine all results into a single matrix
final_sgd_sub_results <- do.call(cbind, sgd_sub_conf_list)
final_sgd_plug_results <- do.call(cbind, sgd_plug_conf_list)
final_momentum_sub_results <- do.call(cbind, momentum_sub_conf_list)
final_momentum_plug_results <- do.call(cbind, momentum_plug_conf_list)

# Define the file name for the combined results
final_file_name_1 <- paste0("sgd_sub_rand_unbalanced_a",min(grid_a), "-", max(grid_a),"_eta", eta, "_n_iter", n_iter, "_m", m, "_K", K, "_b", b, "_lambda0_", index, ".csv")
final_file_name_2 <- paste0("sgd_plug_in_unbalanced_a",min(grid_a), "-", max(grid_a),"_eta", eta, "_n_iter", n_iter, "_m", m, "_K", K, "_b", b, "_lambda0_", index, ".csv")
final_file_name_3 <- paste0("momentum_sub_rand_unbalanced_a",min(grid_a), "-", max(grid_a),"_eta", eta, "_momentum_coef", momentum_coef, "_n_iter", n_iter, "_m", m, "_K", K, "_b", b, "_lambda0_", index, ".csv")
final_file_name_4 <- paste0("momentum_plug_in_unbalanced_a",min(grid_a), "-", max(grid_a),"_eta", eta, "_momentum_coef", momentum_coef, "_n_iter", n_iter, "_m", m, "_K", K, "_b", b, "_lambda0_", index, ".csv")

# Save the combined results to a CSV file
write.csv(final_sgd_sub_results, file = final_file_name_1, row.names = FALSE)
write.csv(final_sgd_plug_results, file = final_file_name_2, row.names = FALSE)
write.csv(final_momentum_sub_results, file = final_file_name_3, row.names = FALSE)
write.csv(final_momentum_plug_results, file = final_file_name_4, row.names = FALSE)














