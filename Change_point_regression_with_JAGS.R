### July 2021 - Kilian Eichenseer
### Change point regression analysis - For further description see keichenseer.com/...
###
#
### Generate the data
set.seed(10) # change the seed for a different sequence of random numbers
n <- 60 # number of total data points
n_shift <- 35 # the data point at which we introduce a change
x <- rnorm(n,0,1) # generate x
y <- rnorm(n,0,0.5) + 0.5 * x # generate y without a change
y[n_shift:n] <- rnorm(length(n_shift:n),0,1) + 1 * x[n_shift:n] + 0.75 # introduce change
#
#
### Plot the data
phase_1 <- 1:(n_shift-1)
phase_2 <- n_shift:n
phase_col <- rep(rgb(0,0.3,1,0.75), n)
phase_col[phase_2] <- rgb(0.9,0.4,0,0.75)
#
layout(matrix(c(1,1,1,1,1,2,2,2), nrow = 1, ncol = 8, byrow = TRUE))
par(mar = c(4,4,1,1), las = 1, mgp = c(2.25,0.75,0), cex = 1.35)
#
plot(x, type = "o", pch = 19, ylim = c(min(y),max(y)), cex = 0.6, xlab = "time", ylab = "")
abline(v = 34.5, lty = 3, lwd = 1.5)
points(y, type = "o", col = "red", pch = 19, cex = 0.6)
legend("topleft", legend = c("x","y"), col = c("black","red"), pch = 19, lwd = 1.5, pt.cex = 0.6, bty = "n")
#
plot(x,y, type = "n")
abline(h=0, v=0, lty = 3)
#points(c(min(x[phase_1]),max(x[phase_1])), c(min(x[phase_1]), max(x[phase_1])) *  coefficients(lm(y[phase_1] ~ x[phase_1]))[2] + coefficients(lm(y[phase_1] ~ x[phase_1]))[1], type = "l", col = rgb(0,0.3,1) , lwd = 2)
#points(c(min(x[phase_2]),max(x[phase_2])), c(min(x[phase_2]), max(x[phase_2])) *  coefficients(lm(y[phase_2] ~ x[phase_2]))[2] + coefficients(lm(y[phase_2] ~ x[phase_2]))[1], type = "l", col = rgb(0.9,0.4,0), lwd = 2)
points(x,y, bg = phase_col, pch = 21, cex = 0.9)
legend("topleft", legend = c(expression("t"[1]~"- t"[34]), expression("t"[35]~"- t"[60])), pt.bg  = c(rgb(0,0.3,1), rgb(0.9,0.4,0)), pch = 21,  pt.cex = 1, bty = "n")
#
#
### Create the Model for JAGS
model_CPR <- function(){
  
  ### Likelihood or data model part
  for(i in 1:n){
    
    y[i] ~ dnorm(mu[i], tau[i]) 
    
    mu[i] <- alpha_1 + 
      alpha_2 * step(i - n_change) +
      (beta_1 + beta_2 * step(i - n_change))*x[i]
    
    tau[i] <- exp(log_tau[i])
    
    log_tau[i] <- log_tau_1 + log_tau_2 * 
      step(i - n_change)
  } 
  
  ### Priors
  alpha_1 ~ dnorm(0, 1.0E-4)
  alpha_2 ~ dnorm(0, 1.0E-4)
  
  beta_1 ~ dnorm(0, 1.0E-4)
  beta_2 ~ dnorm(0, 1.0E-4)
  
  log_tau_1 ~ dnorm(0, 1.0E-4)
  log_tau_2 ~ dnorm(0, 1.0E-4)
  
  K ~ dcat(p)
  n_change <- possible_change_points[K]
  
}
#
#
### Prepare data for passing it to JAGS
# minimum number of the data points before and after the change
min_segment_length <- 5 
#
# assign indices to the potential change points we allow
possible_change_points <- (1:n)[(min_segment_length+1):(n+1-min_segment_length)] 
#
# number of possible change points
M <- length(possible_change_points)  
#
# probabilities for the discrete uniform prior on the possible change points, 
# i.e. all possible change points have the same prior probability
p <- rep(1 / M, length = M) 
# 
# save the data to a list for jags
data_CPR <- list("x", "y", "n", "possible_change_points", "p")
#
#
### Run the Regression with JAGS
require(R2jags) 
#
CPR  <- jags(data = data_CPR, 
             parameters.to.save = c("alpha_1", "alpha_2", 
                                    "beta_1","beta_2",
                                    "log_tau_1","log_tau_2",
                                    "n_change"), 
             n.iter = 2000, 
             n.chains = 3,
             model.file = model_CPR)
#
#
### Inspect the results
# First, a traceplot
library(ggmcmc)
CPR.ggs <- ggs(as.mcmc(CPR)) # convert to ggs object
ggs_traceplot(CPR.ggs, family = "n_change") 
#
# Posterior probabilities for change point:
ggplot(data = CPR.ggs %>% filter(Parameter == "n_change"),
       aes(x=value, y = 3*(..count..)/sum(..count..), fill = as.factor(Chain))) + 
  geom_vline(xintercept = 35,lty = 2) + geom_bar(position = "identity", alpha = 0.5) +
  ylab("posterior probability") + xlab("n_change") + labs(fill='Chain')
#
# "In which interval does the change point fall with 90 % probability?"
quantile(CPR$BUGSoutput$sims.list$n_change, probs = c(0.05, 0.95))
#
# The probability that the change point falls in the interval 34 to 38:
round(length(which(CPR$BUGSoutput$sims.list$n_change %in% 34:38))/(CPR$BUGSoutput$n.sims),2)
#
#
### Plot regression parameters:
CPRm <- CPR$BUGSoutput$mean
CPRs <- CPR$BUGSoutput$sims.list
#
par(mar = c(2.75,4,0,1), las = 1, mfrow = c(1,1))
plot(0,0,type = "n", xlim = c(0.5,6.5), xaxs = "i", ylim= c(-0.2,2.18), xaxt = "n", xlab = NA, ylab = "mean and 90 % CI")
axis(1,at = 1:6, label = expression(alpha[1],alpha[1]+alpha[2],
                                    beta[1],beta[1]+beta[2],
                                    sigma[1]^2,sigma[1]^2+sigma[2]^2))
points(1:6,c(CPRm$alpha_1,CPRm$alpha_1+CPRm$alpha_2, 
             CPRm$beta_1,CPRm$beta_1+CPRm$beta_2,
             1/(exp(CPRm$log_tau_1)),1/(exp(CPRm$log_tau_1+CPRm$log_tau_2))), 
       pch = 19, cex = 1.4)
points(c(1,1), c(quantile(CPRs$alpha_1,probs = c(0.05,0.95))), type = "l", lwd = 2)
points(c(2,2), c(quantile(CPRs$alpha_1+CPRs$alpha_2,probs = c(0.05,0.95))), type = "l", lwd = 2)
points(c(3,3), c(quantile(CPRs$beta_1,probs = c(0.05,0.95))), type = "l", lwd = 2)
points(c(4,4), c(quantile(CPRs$beta_1+CPRs$beta_2,probs = c(0.05,0.95))), type = "l", lwd = 2)
points(c(5,5), c(quantile(1/(exp(CPRs$log_tau_1)),probs = c(0.05,0.95))), type = "l", lwd = 2)
points(c(6,6), c(quantile(1/(exp(CPRs$log_tau_1+CPRs$log_tau_2)),probs = c(0.05,0.95))), type = "l", lwd = 2)
#
#
### Visualise the change point regression:
change_point <- as.numeric(names(sort(table(CPRs$n_change),decreasing = T)))[1] # mode as the change point
phase_1 <- 1:(change_point-1)
phase_2 <- change_point:n
phase_col <- rep(rgb(0,0.3,1,0.75), n)
phase_col[phase_2] <- rgb(0.9,0.4,0,0.75)
#
par(mar = c(4,4,0.3,0.3), las = 1, mgp = c(2.25,0.75,0), cex = 1.35)
#
reg1_seq <- seq(min(x[phase_1]),max(x[phase_1]),length.out = 100)
reg2_seq <- seq(min(x[phase_2]),max(x[phase_2]),length.out = 100)
#
reg1 <- CPRm$alpha_1 + CPRm$beta_1*reg1_seq
reg2 <- (CPRm$alpha_1+CPRm$alpha_2) + (CPRm$beta_1+ CPRm$beta_2)*reg2_seq
#
### Calculate confidence intervals
reg1_025 <- sapply(1:100, function(x) quantile(CPRs$alpha_1 + CPRs$beta_1*reg1_seq[x], probs = 0.025))
reg1_975 <- sapply(1:100, function(x) quantile(CPRs$alpha_1 + CPRs$beta_1*reg1_seq[x], probs = 0.975))
reg2_025 <- sapply(1:100, function(x) quantile(CPRs$alpha_1 + CPRs$alpha_2 + (CPRs$beta_1+CPRs$beta_2)*reg2_seq[x], probs = 0.025))
reg2_975 <- sapply(1:100, function(x) quantile(CPRs$alpha_1 + CPRs$alpha_2 + (CPRs$beta_1+CPRs$beta_2)*reg2_seq[x], probs = 0.975))
#
plot(x,y, type = "n")
abline(h=0, v=0, lty = 3)
#
error_polygon <- function(x,en,ep,color) { # A function to facilitate drawing credible intervals around the regression line
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}
#
error_polygon(reg2_seq,reg2_025,reg2_975,rgb(0.9,0.4,0,0.22))
error_polygon(reg1_seq,reg1_025,reg1_975,rgb(0,0.3,1,0.22))
#
points(reg1_seq, reg1, type = "l", col = rgb(0,0.3,1) , lwd = 2)
points(reg2_seq, reg2, type = "l", col = rgb(0.9,0.4,0), lwd = 2)
points(x,y, bg = phase_col, pch = 21, cex = 0.9)
legend("topleft", legend = c(expression("before "*italic(n["change"])), expression("after "*italic(n["change"]))), pt.bg  = c(rgb(0,0.3,1), rgb(0.9,0.4,0)), pch = 21,  pt.cex = 1, bty = "n", cex = 0.85)
