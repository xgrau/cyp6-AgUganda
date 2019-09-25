# libraries
library(stats4)


# input
# table
x = 2001:2020
y = abs(seq(from=0.1, to=0.9, length.out = 20) + rnorm(mean=0,sd=0.1,20))
plot(x,y, ylim=c(0,1), col="blue", xlab = "year", ylab="freq")


# Lynd et al MBE 2010
pp = (p^2*(1+s) + p*(1-p)(1+hs)) / W

# function to estimate
LL = function(beta0, beta1, mu, sigma) {
  # Find residuals
  R = y - x * beta1 - beta0
  # Calculate the likelihood for the residuals (with mu and sigma as parameters)
  R = suppressWarnings(dnorm(R, mu, sigma))
  # Sum the log likelihoods for all of the data points
  -sum(log(R))
}

fit <- mle(LL, start = list(beta0 = 2, beta1 = 1, mu = 0.5, sigma=0.1))




stop("ara")

# fit linear model
x <- runif(N)
y <- 5 * x^3 + x + 3 + rnorm(N)
y <- 5 * x + 3 + rnorm(N)
fit <- lm(y ~ x)
summary(fit)
plot(x, y)

# function to estimate
LL <- function(beta0, beta1, mu, sigma) {
  # Find residuals
  R = y - x * beta1 - beta0
  # Calculate the likelihood for the residuals (with mu and sigma as parameters)
  R = suppressWarnings(dnorm(R, mu, sigma))
  # Sum the log likelihoods for all of the data points
  -sum(log(R))
}


fit <- mle(LL, start = list(beta0 = 4, beta1 = 2, mu = 0, sigma=1))
summary(fit)
