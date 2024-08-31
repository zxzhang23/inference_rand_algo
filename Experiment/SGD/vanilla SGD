
##############################################################
### SGD with Polyak-averaging for logistic regression ########
##############################################################
logistic_regression_asgd_online <- function(learning_rate, num_iterations, batch_size, average = TRUE, heavy_ball = FALSE, momentum = 0.9,burnin= round(nrow(X) / 1000),start=NULL) {
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
    
    if (heavy_ball) {
      # Update velocity and weights using heavy ball iteration
      velocity <- momentum * velocity + grad
      w <- w - learning_rate[i]* velocity
    } else {
      # Update the weights using standard SGD
      w <- w - learning_rate[i] * grad
    }
    
    # Update the average weights if averaging is enabled
    if (average) {
      avg_w <- (avg_w * (i - 1) + w) / i
    }
  }
  return(list(avg_w = avg_w))
}

###CI for sub-randomamization, multi-run plug-in and aggregation
conf_subrand_asgd_online<-function(learning_rate, m, b, K, batch_size, alpha=0.1,average = TRUE, heavy_ball = FALSE, momentum = 0.9,burnin= 0){
  results <- logistic_regression_asgd_online(learning_rate, m, batch_size, average = TRUE, heavy_ball = FALSE, momentum = 0.9,burnin=0)
  center<-results$avg_w[5]
  
  find_start <- logistic_regression_asgd_online(learning_rate, 1000, batch_size, average = TRUE, heavy_ball = FALSE, momentum = 0.9,burnin=0)
  startb <- find_start$avg_w
  
  subske=0
  for(l in 1:K){
    resultb <- logistic_regression_asgd_online(learning_rate, b,batch_size, average = TRUE, heavy_ball = FALSE, momentum = 0.9,burnin=0,start=startb)
    subske[l]<-resultb$avg_w[5]
  }
  
  taum<-sqrt(m)
  taub<-sqrt(b)
  lb<-unname(quantile(taub*(subske-center),alpha/2)/(taum-taub))
  rb<-unname(quantile(taub*(subske-center),(1-alpha/2))/(taum-taub))
  sub_conf<-rbind(center-rb,center-lb)
  
  
  plug_v<-var(subske)
  plug_rb=qnorm(1-alpha/2,sd=sqrt(plug_v)*taub/taum)
  plug_lb=qnorm(alpha/2,sd=sqrt(plug_v)*taub/taum)
  plug_conf=rbind(center-plug_rb,center-plug_lb)
  
  multi_v<-var(subske)
  multi_rb=qnorm(1-alpha/2,sd=sqrt(multi_v/K))
  multi_lb=qnorm(alpha/2,sd=sqrt(multi_v/K))
  multi_center<-mean(subske)
  multi_run_conf=rbind(multi_center-multi_rb,multi_center-multi_lb)
  
  out <- list(subrand_conf = sub_conf, plug_conf=plug_conf,multi_run_conf=multi_run_conf)
  return(out)
}

