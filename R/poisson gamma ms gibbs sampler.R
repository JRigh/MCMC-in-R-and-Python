#-----------------------------------------
# Multi-stage Poisson-Gamma Gibbs sampling
# in R
# Robert & Casella, 10.17
#-----------------------------------------

library(tidyverse)
library(xtable)

setwd("C:/Users/julia/OneDrive/Desktop/github/20. Gibbs sampling")

# write the data and export in LaTeX table 

dataset = data.frame('Pump' = 1:10,
                     'Failures' = c(5,1,5,14,3,19,1,1,4,22),
                     'Time' = c(94.32, 15.72, 62.88, 125.76, 5.24, 31.44, 1.05, 1.05, 2.10, 10.48))
print(xtable(dataset, type = "latex"), file = "PGMSGS_tables.tex")

# 1. multi-stage gibbs sampler
twostageGibbsSamplerPG <- function(nsim, beta.start, alpha, gamma, delta, y, t, burnin = 0,
                                   thin = 1) {
  beta.draws <- c()
  lambda.draws <- matrix(NA, nrow = nsim, ncol = length(y))
  beta.cur <- beta.start
  
  lambda.update <- function(alpha, beta, y, t) {
    rgamma(length(y), y + alpha, t + beta)
  }
  
  beta.update <- function(alpha, gamma, delta, lambda, y) {
    rgamma(1, length(y) * alpha + gamma, delta + sum(lambda))
  }
  for (i in 1:nsim) {
    lambda.cur <- lambda.update(alpha = alpha, beta = beta.cur, y = y, t = t)
    beta.cur <- beta.update(alpha = alpha, gamma = gamma, delta = delta,
                            lambda = lambda.cur, y = y)
    if (i > burnin & (i - burnin)%%thin == 0) {
      lambda.draws[(i - burnin), ] <- lambda.cur
      beta.draws[(i - burnin)] <- beta.cur
    }
  }
  return(list(lambda.draws = lambda.draws, beta.draws = beta.draws))
}


# posterior distributions
set.seed(2023)
posterior <- twostageGibbsSamplerPG(nsim = 10000, beta.start = 1, alpha = 1.8, gamma = 0.01,
                                    delta = 1, y = dataset[,2], t = dataset[,3])

# 2. posterior quantities
mean_lambda = mean(posterior$lambda.draws)
mean_beta = mean(posterior$beta.draws)
sd_lambda = sd(posterior$lambda.draws)
sd_beta = sd(posterior$beta.draws)
df1 <- matrix(c(mean_lambda, mean_beta, sd_lambda, sd_beta), byrow = TRUE, ncol =2)
colnames(df1) = c('lambda', 'beta')
rownames(df1) = c('post mean', 'post sd')

quantiles_lambda = quantile(posterior$lambda.draws, probs = c(0.025, 0.5, 0.975))
quantiles_beta = quantile(posterior$beta.draws, probs = c(0.025, 0.5, 0.975))
df2 <- matrix(c(quantiles_lambda, quantiles_beta), byrow = FALSE, ncol =2)
colnames(df2) = c('lambda', 'beta')
rownames(df2) = c('2.5%', '50%', '97.5%')

post <- rbind(df1, df2)
setwd("C:/Users/julia/OneDrive/Desktop/github/20. Gibbs sampling")
print(xtable(post, type = "latex", digits=4), file = "tables.tex")

# 3. plot
p1 <- ggplot(data.frame('x' = matrix(posterior$lambda.draws, ncol = 1)), aes(x = x)) +
  geom_histogram( color = 'black', fill = 'darkred', binwidth = 0.1, aes(y = ..density..)) +
  geom_density(color = 'black', lwd = 1.4) +
  labs(title = 'Histogram of posterior distribution of lambda parameter',
       y="count", x="lambda") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.subtitle=element_text(size=10, face="italic", color="darkred"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey90"))

p2 <- ggplot(data.frame('x' = matrix(posterior$beta.draws, ncol = 1)), aes(x = x)) +
  geom_histogram( color = 'black', fill = 'darkred', binwidth = 0.1, aes(y = ..density..)) +
  geom_density(color = 'black', lwd = 1.4) +
  labs(title = 'Histogram of posterior distribution of beta hyperparameter',
       y="count", x="beta") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.subtitle=element_text(size=10, face="italic", color="darkred"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey90"))

grid.arrange(p1, p2, nrow = 1)

#----
# end
#----
