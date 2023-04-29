############# Code for Optional Stopping in a Bayesian Framework for Misspecified Models ##################
# Set seed
set.seed(123)

# Libraries
library(tidyverse)
############# Example 1: Sample size N = 100 ########################
# H0
tau <- 1
n <- 100
H0_reps <- replicate(50000, {
xbar <-  rnorm(1, mean = 0, sd = sqrt(1/n))
 num <- exp(((n^2)*(xbar)^2)/(2*(n+tau^2)))
 den <- sqrt(n+tau^2)
 BF <- num/den
})

tau <- 1
n <- 100
H1_reps <- replicate(50000, {
  delta <- rnorm(1, 0, sd = sqrt(tau^2))
  xbar <- rnorm(1, mean = delta, sd = sqrt(1/n))
  num <- exp(((n^2)*(xbar)^2)/(2*(n+tau^2)))
  den <- sqrt(n+tau^2)
  BF <- num/den
})



# Combine results in dataframe
dat1 <- as.data.frame(cbind(H0_reps, H1_reps))


# Histogram plot
png(filename = "~/Desktop/example1N.png", width = 20, height = 15, units = "cm", res = 600)
ggplot(data = dat1) + 
  geom_histogram(aes(x = log(H1_reps), y = ..count..), colour = "black", fill = "grey", bins = 20) +  
  geom_histogram(aes(x = log(H0_reps), y = -..count..), colour = "black", fill = "red", bins = 20) + theme_bw() + 
   labs(x = "Bayes Factor (Log Scale)", y = "Counts") + scale_x_continuous(limits = c(-0.1, 3))
dev.off()
###################### Example 1: optional stopping ###############################
# H0
tau <- 1
H0_reps <- replicate(50000, {
n <- seq(1, 100, 1)
BF <- c()
nlast <- 0
for (i in seq_along(n)){
  xbar <-  rnorm(1, mean = 0, sd = sqrt(1/n[i]))
  num <- exp(((n[i]^2)*(xbar)^2)/(2*(n[i]+tau^2)))
  den <- sqrt(n[i]+tau^2)
  BF[i] <- num/den
  if (BF[i] < 1/10 || BF[i] > 10) {
    break
  }
}
BF <- BF[i]
})



# H1
tau <- 1
H1_reps <- replicate(50000, {
  n <- seq(1, 100, 1)
  BF <- c()
  nlast <- 0
  for (i in seq_along(n)){
    delta <- rnorm(1, 0, sd = sqrt(tau^2))
    xbar <- rnorm(1, mean = delta, sd = sqrt(1/n[i]))
    num <- exp(((n[i]^2)*(xbar)^2)/(2*(n[i]+tau^2)))
    den <- sqrt(n[i]+tau^2)
    BF[i] <- num/den
    if (BF[i] < 1/10 || BF[i] > 10) {
      break
    }
  }
  BF <- BF[i]
})


# Combine results in dataframe
dat2 <- as.data.frame(cbind(H0_reps, H1_reps))


# Histogram plot
png(filename = "~/Desktop/example_stopping.png", width = 20, height = 15, units = "cm", res = 600)
ggplot(data = dat2) + 
  geom_histogram(aes(x = log(H1_reps), y = ..count..), colour = "black", fill = "grey", bins = 60) +  
  geom_histogram(aes(x = log(H0_reps), y = -..count..), colour = "black", fill = "red", bins = 60) + theme_bw() + 
   labs(x = "Bayes Factor (Log Scale)", y = "Counts") 
dev.off()
######################## Plot for Bayes factor change under null hypothesis simulation ############################
# H0
n <- seq(1, 2000, 1)
BF <- c()
nlast <- 0
tau <- 1
for (i in seq_along(n)){
xbar <-  rnorm(1, mean = 0, sd = sqrt(1/n[i]))
num <- exp(((n[i]^2)*(xbar)^2)/(2*(n[i]+tau^2)))
den <- sqrt(n[i]+tau^2)
BF[i] <- num/den
}

# Convert to dataframe 
dat3 <- as.data.frame(cbind(n , BF))

png(filename = "~/Desktop/bayes_over_n.png", width = 20, height = 15, units = "cm", res = 600)
ggplot(dat3, aes(n, log(BF))) + geom_line(colour = "blue") + 
  theme_bw() +  labs(x = "Sample Size (n)", y = "Bayes Factor") + scale_y_continuous(limits = c(0, 1.5)) + 
 geom_hline(yintercept = 1)
dev.off()

######################## Example 2: Misspecified prior in model #########################
# Generate data from null
n <- seq(1, 2000, 1)
BF <- c()
for (i in seq_along(n)) {
# Directionality hypotheses
delta1 <- 0.5
delta2 <- -0.5

# Generate data
x <- rnorm(n[i], mean = 0, sd = 1)

# Calculate sample mean
xbar <- mean(x)

BF[i] <- exp(-(xbar - delta1)^2/(2/n[i]))/exp(-(xbar + delta1)^2/(2/n[i]))
}

# Convert to dataframe 
dat4 <- as.data.frame(cbind(n , BF))

png(filename = "~/Desktop/model_misspecification.png", width = 20, height = 15, units = "cm", res = 600)
ggplot(dat4, aes(n, log(BF))) + geom_line(colour = "blue") + 
  theme_bw() +  labs(x = "Sample Size (n)", y = "Log Bayes Factor")  + 
  geom_hline(yintercept = 1)
dev.off()

##################################### Plot posterior odds with optional stopping ########################
# Without optional stopping (N = 50)
dat_N <- replicate(50000, {
  n <- 50
    # Generate data
    x <- rnorm(n, mean = 0, sd = 1)
    
    # Calculate sample mean
    xbar <- mean(x)
    
    BF <- exp(-(xbar - delta1)^2/(2/n))/exp(-(xbar + delta1)^2/(2/n))
})


# Histogram for optional stopping 
dat6 <- as.data.frame(cbind(x = dat_N))

png(filename = "~/Desktop/hist_N.png", width = 20, height = 15, units = "cm", res = 600)
ggplot(data = dat6) + 
  geom_histogram(aes(x = log(dat_N), y = ..count..), colour = "black", fill = "grey", bins = 60)+
 theme_bw() + labs(x = "Log Posterior Odds", y = "Count")
dev.off()


# Histogram for optional stopping 
dat5 <- as.data.frame(cbind(x = dat_opt))

png(filename = "~/Desktop/hist_opt.png", width = 20, height = 15, units = "cm", res = 600)
ggplot(data = dat5) + 
  geom_histogram(aes(x = log(dat_opt), y = ..count..), colour = "black", fill = "grey", bins = 60) + 
  scale_x_continuous(limits = c(0, 10)) + theme_bw() + labs(x = "Log Posterior Odds", y = "Count")
dev.off()

# With optional stopping
dat_opt <- replicate(50000, {
  n <- seq(1, 100, 1)
  BF <- c()
  nlast <- 0
  for (i in seq_along(n)){
    # Generate data
    x <- rnorm(n[i], mean = 0, sd = 1)
    
    # Calculate sample mean
    xbar <- mean(x)
    
    BF[i] <- exp(-(xbar - delta1)^2/(2/n[i]))/exp(-(xbar + delta1)^2/(2/n[i]))
    if (BF[i] < 1/10 || BF[i] > 10) {
      break
    }
  }
  BF <- BF[i]
})

# Histogram for optional stopping 
dat5 <- as.data.frame(cbind(dat_opt, dat_N))

png(filename = "~/Desktop/hist_opt.png", width = 20, height = 15, units = "cm", res = 600)
ggplot(data = dat5) + 
  geom_histogram(aes(x = log(dat_opt), y = ..count..), colour = "black", fill = "red", bins = 60) + 
  geom_histogram(aes(x = log(dat_N), y = -..count..), colour = "black", fill = "grey", bins = 60) + 
 theme_bw() + labs(x = "Log Bayes Factor", y = "Count")
dev.off()