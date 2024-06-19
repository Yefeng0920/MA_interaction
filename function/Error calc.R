###########
# Function to calculate power (two-tail) for meta-analysis
power.ma_Shinichi <- function(mu, SE, alpha = 0.05) {
  2-pnorm(qnorm(1-alpha/2)-abs(mu)/SE)-pnorm(qnorm(1-alpha/2)+abs(mu)/SE)
} # or power.ma_Shinichi1 <- function(mu,SE){1 - pnorm(qnorm(1-0.05/2)-abs(mu)/SE) + pnorm(-qnorm(1-0.05/2)-abs(mu)/SE)}


# Function for power analysis for empirical data point
power.individual_Shinichi <- function(mu, se, alpha = 0.05) {
  2-pnorm(qnorm(1-alpha/2)-abs(mu)/se)-pnorm(qnorm(1-alpha/2)+abs(mu)/se)} # two-tailed power


# Function for Type S error for empirical data point
error_S <- function(mu, se, alpha = 0.05){
  #z <- qnorm(1 - alpha/2) # Z-score or quantile
  p.u <- 1 - pnorm(qnorm(1 - alpha/2) - abs(mu)/se) # upper-tail probability
  p.l <- pnorm(-qnorm(1 - alpha/2) - abs(mu)/se) # lower-tail probability
  power <- p.u + p.l # upper + lower
  errorS <- p.l/power # percentage of the opposite direction
  return(errorS)
} 

# Function for Type M error for empirical data point
error_M <- function(mu, se, alpha = 0.05, N = 10000) {
  est.random <- rnorm(n=N, mean = mu, sd = se)
  # est.random <- mu + se*rnorm(n=N, mean=0, sd=1)
  sig.index <- abs(est.random) > se*qnorm(1 - alpha/2)
  overestimate <- mean(abs(est.random)[sig.index])/abs(mu) # ratio is regardnesss of sign, so we need absolute value
  absolute_error <- overestimate*abs(mu) - abs(mu)
  relative_error <- absolute_error/(overestimate*abs(mu))
  return(abs(overestimate) %>% round(3))
}