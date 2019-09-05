# table
g = seq(from=0,to=20,length.out = 20) * 12
p = abs(seq(from=0.1, to=0.8, length.out = 20) + rnorm(mean=0,sd=0.1,20))
plot(g,p, col="blue", xlab = "months", ylab="freq")


pp = (p^2*(1+s) + p*(1-p)(1+h*s)) / W


# loss function
lm.loss <- function(par) {
  
  h.par = par[1]
  s.par = par[2]
  W.par = par[3]
  err.sigma = par[4]
  
  if(err.sigma < 0) {deviance <- 10000000}
  if(err.sigma > 0) {
    likelihoods <- dnorm(p, mean = g * h.par + s.par + W.par, sd = err.sigma)
    log.likelihoods <- log(likelihoods)
    deviance <- -2 * sum(log.likelihoods)
  }
  return(deviance)
}

lm.loss(c(1,5,20,0.1))

# optimise
parameter.fits <- optim(par = c(1,5,20,0.1),
                        fn = lm.loss, hessian = T
)
parameter.fits$par



stop("Ara")




### input data

N = 1000

# variable
data.x <- rnorm(n = N, mean = 10, sd = 2)



# true parameters
a.true <- 3
b.true <- 8
true.y <- data.x * a.true + b.true

# response
err.sd.true <- 1  # Set noise sd
noise <- rnorm(n = N, mean = 0, sd = 2)  # Generate noise
data.y <- true.y + noise  # Add noise to true (latent) responses

# loss function
lm.loss <- function(par) {
  
  a.par <- par[1]  # The current slope
  b.par <- par[2]  # The current intercept
  err.sigma <- par[3]  # The current error standard deviation
  
  # If the error standard deviation is invalid (i.e.; negative), then we need to return a very high deviance
  # This will tell the optimization procedure to stay away from invalid (either statistically or psychologically)
  #   parameter values.
  
  if(err.sigma < 0) {deviance <- 10000000}
  
  # If the error standard deviation is valid (i.e.; > 0), then calculate the deviance...
  
  if(err.sigma > 0) {
    
    # Calculate the likelihood of each data point.
    # Here is where you specify the model and how you calculate likelihoods.
    
    likelihoods <- dnorm(data.y, mean = data.x * a.par + b.par, sd = err.sigma)
    
    # Now let's convert the vector of likelihoods to a summary deviance score (-2 times sub log-lik)
    
    # Calculate log-likelihoods of each data point
    log.likelihoods <- log(likelihoods)
    
    # Calculate the deviance (-2 times sum of log-likelihoods)
    deviance <- -2 * sum(log.likelihoods)
    
  }
  
  return(deviance)
  
}


# tries how deviant is the function
dev.temp <- lm.loss(c(1,5,20))
dev.temp

dev.temp <- lm.loss(c(3,8,2))
dev.temp

# optimise
parameter.fits <- optim(par = c(0,0,1.9),
                        fn = lm.loss, hessian = T
)
parameter.fits$par

hessian <- parameter.fits$hessian
hessian.inv <- solve(hessian)
parameter.se <- sqrt(diag(hessian.inv))
parameter.se

