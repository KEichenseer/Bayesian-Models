######################################################
######################################################
### A Metropolis algorithm in R - Part 2: 
### Adaptive proposals
###
### by Kilian Eichenseer
### March 2022
###
######################################################
######################################################

######################################################
######################################################
### FIGURE 1 - illustration of Metropolis proposals

set.seed(2)
par(mfrow=c(2,1),mar = c(0.5,4,0,0), las = 1)
xd <- seq(-12.2,12.2,0.01)
yd <- dnorm(xd,-1.5,1)+dnorm(xd,1.5,1)
plot(xd,yd, xlab = "", ylab = "density",type = "n", yaxs = "i", xaxs = "i", ylim = c(0,0.43), xaxt = "n", yaxt = "n")
axis(2,seq(0,0.5,0.1), c("0",as.character(seq(0.1,0.4,0.1)),NA))
text(-11.5,0.39,"target distribution", adj = c(0,.5), cex = 1.1)

error_polygon <- function(x,en,ep,color) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}

col3 <- rgb(0.8,0.8,0,.75)
col2 <- rgb(0.8,0.5,0,.75)
col1 <- rgb(1,0,0,.75)


error_polygon(xd,rep(0,length(xd)),yd,rgb(0,0.35,0.7,0.33))
p1 <- -0.75
sd3 <- p1+rnorm(10,0,0.2)
points(sd3,dnorm(sd3,-1.5,1)+dnorm(sd3,1.5,1), pch = 4, cex = 1.2, lwd = 2, col = col1, xpd = T)

sd2 <- p1+rnorm(10,0,2.5)
points(sd2,dnorm(sd2,-1.5,1)+dnorm(sd2,1.5,1), pch = 4, cex = 1.2, lwd = 2, col = col2, xpd = T)

sd1 <- p1+rnorm(10,0,8)
points(sd1,dnorm(sd1,-1.5,1)+dnorm(sd1,1.5,1), pch = 4, cex = 1.2, lwd = 2, col = col3, xpd = T)

abline(v= -0.75, lty = 3, lwd = 2)
points(-0.75,dnorm(-0.75,-1.5,1)+dnorm(-0.75,1.5,1), pch = 4, cex = 1.75, lwd = 3, col = rgb(0,0,0,0.67), xpd = T)
sigmas <- c(8,2.5,0.2)
legend("topright", as.expression(c("current value",sapply(1:3, function(x) bquote(italic(sigma["proposal"]~"="~.(sigmas[x])))))), 
       col = c(rgb(0,0,0,0.67),col1,col2,col3), border = NA, bty = "n", pch = 4, 
       pt.cex = c(1.7,1.2,1.2,1.2), pt.lwd = c(3,2,2,2), cex = 0.95)


par(mar = c(2,4,0.5,0))

plot(xd,yd, xlab = "", ylab = "density (not to scale)",type = "n", yaxs = "i", xaxs = "i", ylim = c(0,0.43), yaxt = "n")
text(-11.5,0.39,"proposal distributions", adj = c(0,.5), cex = 1.1)
axis(2,seq(0,0.5,0.1), c("0",as.character(seq(0.1,0.4,0.1)),NA))

col1 <- rgb(0.8,0.8,0,0.35)
col2 <- rgb(0.8,0.5,0,0.4)
col3 <- rgb(1,0,0,0.4)

error_polygon(xd,rep(0,length(xd)),dnorm(xd,p1,8)*4,col1)
error_polygon(xd,rep(0,length(xd)),dnorm(xd,p1,2.5)*1.8,col2)
error_polygon(xd,rep(0,length(xd)),dnorm(xd,p1,0.2)*0.2,col3)


abline(v= -0.75, lty = 3, lwd = 2)
legend("topright", as.expression(sapply(1:3, function(x) bquote(italic(sigma["proposal"]~"="~.(sigmas[x]))))), 
       fill = c(col3,col2,col1), border = NA, bty = "n", cex = 0.95)
par(mfrow=c(1,1))
                                 
######################################################
######################################################
### MCMC algorithm
                                 
### FUNCTIONS

gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  return(A + max(c(K-A,0))/((1+(exp(Q*(x-M))))) + rnorm(length(x),0,sdy))
}

loglik <- function(x, y,  coeff, sdy) {
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  pred = A + max(c(K-A,0))/((1+(exp(Q*(x-M)))))
  return(sum(dnorm(y, mean = pred, sd = sdy, log = TRUE)))
}

logprior <- function(coeff) {
  return(sum(c(
    dunif(coeff[1], -4, 40, log = TRUE),
    dunif(coeff[2], -4, 40, log = TRUE),
    dnorm(coeff[3], 45, 10, log = TRUE),
    dlnorm(coeff[4], -2, 1, log = TRUE))))
}

logposterior <- function(x, y, coeff, sdy){
  return (loglik(x, y, coeff, sdy) + logprior(coeff))
}

MH_propose <- function(coeff, proposal_sd){
  rnorm(4,mean = coeff, sd= proposal_sd)
}

weighted_var <- function(x, weights, sum_weights) {
  sum(weights*((x-sum(weights*x)/sum_weights)^2))/(sum_weights)
}

# Main MCMCM function
run_MCMC <- function(x, y, coeff_inits, sdy_init, nIter, proposal_sd_init = rep(5,4), 
                     nAdapt = 5000, adaptation_decay = 500){
  ### Initialisation
  coefficients = array(dim = c(nIter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,nIter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy
  A_sdy = 3 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 0.1 # parameter for the prior on the inverse gamma distribution of sdy
  n <- length(y)
  shape_sdy <- A_sdy+n/2 # shape parameter for the inverse gamma
  sd_it <- 1 # iteration index for the proposal standard deviation
  proposal_sd <- array(NA_real_,dim = c(nAdapt,4)) # array to store proposal SDs
  proposal_sd[1:3,] <- proposal_sd_init # proposal SDs before adaptation
  # pre-define exp. decaying weights for weighted variance
  allWeights <- exp((-(nAdapt-2)):0/adaptation_decay) 
  accept <- rep(NA,nIter) # vector to store the acceptance or rejection of proposals
  ### The MCMC loop
  for (i in 2:nIter){

   ## 1. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(
      1,shape_sdy,B_sdy+0.5*sum((y-gradient(x,coefficients[i-1,],0))^2)))

   ## 2. Metropolis-Hastings step to estimate the regression coefficients
    proposal = MH_propose(coefficients[i-1,],proposal_sd[sd_it,]) # new proposed values
    if(any(proposal[4] <= 0)) HR = 0 else {# Q and nu need to be >0
      # Hasting's ratio of the proposal
      HR = exp(logposterior(x = x, y = y, coeff = proposal, sdy = sdy[i]) -
                 logposterior(x = x, y = y, coeff = coefficients[i-1,], sdy = sdy[i]))}

    #if(gradient(65, proposal,0) >10) HR = 0
    # accept proposal with probability = min(HR,1)
    if (runif(1) < HR){
      accept[i] <- 1
      coefficients[i,] = proposal
      # if proposal is rejected, keep the values from the previous iteration
    }else{
      accept[i] <- 0
      coefficients[i,] = coefficients[i-1,]
    }
    # Adaptation of proposal SD
    if(i < nAdapt){ # stop adaptation after nAdapt iterations
    if(i>=3) {
    weights = allWeights[(nAdapt-i+2):nAdapt-1] # select weights for current iteration
    sum_weights = sum(weights) 
    weighted_var_coeff <- apply(coefficients[2:i,], 2, # calculate weighted variance
          function(f) weighted_var(
          f, weights = weights, sum_weights = sum_weights))


    for(v in 1:4) {if(weighted_var_coeff[v]==0)   { # 
              proposal_sd[i+1,v] <- sqrt(proposal_sd[i,v]^2/5)
      } else  proposal_sd[i+1,v] <- 2.4/sqrt(4) * sqrt(weighted_var_coeff[v]) + 0.0001
    }                   
                           
    }
    sd_it <- i+1
    }
  } # end of the MCMC loop

  ###  Function output
  output = list(coeff = data.frame(A = coefficients[,1],
                           K = coefficients[,2],
                           M = coefficients[,3],
                           Q = coefficients[,4],
                           sdy = sdy),
                proposal_sd = data.frame(proposal_sd),
                accept = accept)
  return(output)
}
                                
# for plotting the 95 % CI shading
error_polygon <- function(x,en,ep,color) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}


### Taking samples
set.seed(10)
sample_lat <- runif(10,0,90)
sample_data <- data.frame(
  x = sample_lat, 
  y = gradient(x = sample_lat, coeff = c(-2.0, 28, 41, 0.1), sd = 2))

### Analysis
nIter <- 100000
print(system.time({m <- run_MCMC(x = sample_data$x, y = sample_data$y,
                                 coeff_inits = c(0,30,45,0.2), sdy_init = 4, 
                                 nIter = nIter, nAdapt = 5000, adaptation_decay =500,
                                 proposal_sd_init = c(5,5,5,5))}))


######################################################
######################################################
### PLOTS

### Traceplot and density
                                
library(ggplot2)
library(dplyr)
library(reshape)
library(cowplot)
mm <- m$coeff %>% dplyr::mutate(iteration = 1:nrow(.))
mm <- melt(mm,id.vars = "iteration")

facet.labs = c("A", "K", "M", "Q", expression(sigma["y"]))
names(facet.labs) <- c("A", "K", "M", "Q", "sdy")

p1 <- ggplot(mm %>% filter(iteration %% 10 == 0)) + geom_line(aes(x = iteration, y = value),
                             colour = rgb(0,0.35,0.7,1))+
      facet_grid(variable ~ ., switch = "y",scales = "free_y",
      labeller = labeller(variable = facet.labs)) + 
      theme(legend.position = "none")+
      theme_bw(base_size = 14)+
      theme(strip.text.y.left = element_text(angle = 0))+
      scale_y_continuous(position = "right")+
      ylab(NULL)+
      ggtitle(expression("traceplot (showing every 10th"~"iteration)"))+ 
      theme(plot.title = element_text(size = 13, face = "bold",hjust = 0.5),
            axis.title=element_text(size=13), axis.text.y=element_blank(),)


p2 <- ggplot(mm %>% filter(iteration %% 10 == 0 & iteration > 5000)) +
  #geom_vline(aes(xintercept = trueval), colour = "black", size = 1,linetype = "solid")+

 stat_density(aes(x = value, y =..scaled..), 
   fill = rgb(0,0.35,0.7,0.55), size = 1.5)+
      facet_grid(variable ~ ., scales = "free_y")+ 
     theme(legend.position = "none")+
     theme_bw(base_size = 14)+
     theme(strip.background = element_blank(),
   strip.text.y = element_blank())+ 
 scale_y_continuous(position = "left",breaks = c(0,0.5,1))+
  scale_x_continuous(position = "top")+
coord_flip()+
  ylab("scaled density")+
      xlab(NULL)+
      ggtitle(expression("density"))+ 
      theme(plot.title = element_text(size = 13, face = "bold",hjust = 0.5),
            axis.title=element_text(size=13))

plot_grid(p1, p2, ncol = 2, labels = c("", ""),
         rel_widths = c(7,2))

### Traceplot of proposal standard deviations
colnames(m$proposal_sd) <- c("A","K", "M", "Q")
mm <- m$proposal_sd %>% dplyr::mutate(iteration = 1:nrow(.))
mm <- melt(mm,id.vars = "iteration")
facet.labs = c("A", "K", "M", "Q")
names(facet.labs) <- c("A", "K", "M", "Q")

ggplot(mm) + geom_line(aes(x = iteration, y = value),
                             colour = rgb(0,0.35,0.7,1))+
      facet_grid(variable ~ ., switch = "y",scales = "free_y",
      labeller = labeller(variable = facet.labs)) + 
      theme(legend.position = "none")+
      theme_bw(base_size = 14)+
      theme(strip.text.y.left = element_text(angle = 0))+
      scale_y_continuous(position = "right")+
      ylab(NULL)+
      ggtitle(expression("traceplot of "*sigma["proposal"]))+ 
      theme(plot.title = element_text(size = 13, face = "bold",hjust = 0.5),
            axis.title=element_text(size=13))
                                
######################################################
######################################################
### effective sample size
round( mcmcse::multiESS(m$coeff[5001:nIter,]),0)
### increase in effective sample size (in percent)
round((mcmcse::multiESS(m$coeff[5001:nIter,])/2453 - 1)*100,0)
                                

