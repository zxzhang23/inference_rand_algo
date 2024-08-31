args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  index <- as.numeric(args[1])  #folder number
  set.seed(index)
} else {
  stop()
}

sigmoid <- function(z) {
  return(1 / (1 + exp(-z)))
}

###  Stochastic Gradient Descent, including vanilla SGD, heavy ball SGD, adam

logistic_regression_sgd_online <- function(learning_rate, num_iterations, batch_size, method = "sgd", average = FALSE, trajectory=FALSE,heavy_ball = FALSE, momentum, adam = FALSE, beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8,start=NULL,burnin=0) {
  weights_history <- matrix(0, nrow = num_iterations + 1, ncol = d)
  avg_weights_history <- matrix(0, nrow = num_iterations + 1, ncol = d)
  
  if(is.null(start)){
    start <- rep(0, d)
    if (burnin > 0) {
      for (i in 1:burnin) {
        X_batch <- t(rnorm(d))
        y_batch <- rbinom(1,1,sigmoid(sum(X_batch*theta)))
        z <- X_batch %*% start
        p <- sigmoid(z)
        grad <- t(X_batch) %*% (p - y_batch) / nrow(X_batch)
        start <- start - learning_rate[i] * grad
      }
    }
  }
  
  weights_history <- matrix(0, nrow = num_iterations + 1, ncol = d)
  weights_history[1, ] <- start
  
  w <- start
  avg_w <- start
  
  # Initialize velocity for heavy ball iteration
  velocity <- rep(0, d)
  
  for (i in 1:num_iterations) {
    X_batch <- t(rnorm(d))
    y_batch <- rbinom(1,1,sigmoid(sum(X_batch*theta)))
    # Compute the gradient for the batch
    z <- X_batch %*% w
    p <- sigmoid(z)
    grad <- t(X_batch) %*% (p - y_batch) 
    
    
    if (method == "sgd") {
      # Update the weights using standard SGD
      w <- w - learning_rate[i] * grad
    } else if (method == "heavy_ball") {
      velocity <- momentum[i] * velocity + grad
      w <- w - learning_rate[i] * velocity
    } else if (method == "adam") {
      # Update variables for Adam
      m <- beta1 * m + (1 - beta1) * grad
      v <- beta2 * v + (1 - beta2) * (grad^2)
      
      # Compute bias-corrected first and second moment estimates
      m_hat <- m / (1 - beta1^i)
      v_hat <- v / (1 - beta2^i)
      
      # Update the weights using Adam
      w <- w - learning_rate[i] * m_hat / (sqrt(v_hat) + epsilon)
    } else {
      stop("Invalid optimization method. Choose from 'sgd', 'heavy_ball', or 'adam'.")
    }
    
    if (average) {
      avg_w <- (avg_w * (i - 1) + w) / i
    }
    
    # Store the weights and average weights at the current iteration
    weights_history[i + 1, ] <- w
    if (average) {
      avg_weights_history[i + 1, ] <- avg_w
    }
  }
  
  if (average & trajectory) {
    return(list(weights_history = weights_history, avg_weights_history = avg_weights_history))
  } else if(trajectory & !average){
    return(list(weights_history = weights_history))
  }else if(!trajectory & average){list(avg_w=avg_w)
  }else{list(last_w=w)}
}


d <- 5
theta<-seq(0, d,d/4)/(d)
num_iterations <- 100000
batch_size <- 1 

alpha=0.1
grid_b<-seq(2000,10000,4000)
grid_m<-rep(num_iterations,length(grid_b))
grid_K=rep(50,length(grid_b))
sim=500

##parameters for learning rate eta/t^a
eta<-0.3
grid_a=c(0.5,0.9)
conf_list <- list()
for (a_idx in 1:length(grid_a)) {
  conf_list[[a_idx]] <- matrix(0, sim, 6*length(grid_m))
}


pivo_conf=sub_conf=boot_conf=plug_conf=multi_run_conf=multi_run_conf1=matrix(0,sim,2*length(grid_m))



for (a in grid_a) {
  learning_rate <- eta*(1:num_iterations)^(-a) 
  momentum<-1-3*sqrt(learning_rate)  #momentum can be asjusted here
  
  for (i in 1:sim) {
    for (j in 1:length(grid_m)) {
      m <- grid_m[j]
      results <- logistic_regression_sgd_online(learning_rate, m, batch_size,method='heavy_ball',momentum=momentum)
      center <- results$last_w[5]
      
      subske <- 0
      b <- grid_b[j]
      K <- grid_K[j]
      
      find_start <- logistic_regression_sgd_online(learning_rate, 1000, batch_size, method='heavy_ball', momentum = momentum, start=NULL,burnin=0)
      startb <- find_start$last_w
      
      for (l in 1:K) {
        resultb <- logistic_regression_sgd_online(learning_rate, b, batch_size,method='heavy_ball',momentum=momentum,start=startb,burnin = 0)
        subske[l] <- resultb$last_w[5]
      }
      
      taum <- m^(a/4)
      taub <- b^(a/4)
      lb <- quantile(taub*(subske-center), alpha/2)/(taum-taub)
      rb <- quantile(taub*(subske-center), (1-alpha/2))/(taum-taub)
      sub_conf[i, (2*j-1):(2*j)] <- c(center-rb, center-lb)
      
      plug_v <- var(subske)
      plug_rb <- qnorm(1-alpha/2, sd=sqrt(plug_v)*taub/taum)
      plug_lb <- qnorm(alpha/2, sd=sqrt(plug_v)*taub/taum)
      plug_conf[i, (2*j-1):(2*j)] <- c(center-plug_rb, center-plug_lb)
      
      multi_v <- var(subske)
      multi_rb <- qnorm(1-alpha/2, sd=sqrt(multi_v/K))
      multi_lb <- qnorm(alpha/2, sd=sqrt(multi_v/K))
      multi_center <- mean(subske)
      multi_run_conf[i, (2*j-1):(2*j)] <- c(multi_center-multi_rb, multi_center-multi_lb)
    }
  }
  
  conf <- cbind(sub_conf, plug_conf, multi_run_conf)
  conf_list[[which(grid_a == a)]] <- conf
}

write.csv(conf_list[[1]],paste("0logistic_sgd_hb_conf_eta0.3_a0.5_n1e6p5m1e5b2e3to10e3K50rep500_",index,".csv",sep=""))
write.csv(conf_list[[2]],paste("0logistic_sgd_hb_conf_eta0.3_a0.9_n1e6p5m1e5b2e3to10e3K50rep500_",index,".csv",sep=""))


