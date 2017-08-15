library(msm)
library(ggplot2)

simParamExp <- function(n=10^6){
  
  # randomly initialise
  theta1 <- runif(n,27,30)
  theta2 <- rtnorm(n,0.4,0.15,0,Inf)
  
  param <- data.frame(theta1,theta2)
  
}

simDataExp <- function(param,num_tp = 5,noise = 3,high_m = 4,high_sd = 0.5,
                       low_m = 1,low_sd = 2 ){
  
  n <- dim(param)[1]
  
  # double gaussian for start time
  start_time <- mat.or.vec(n,1)
  
  high_low <- sample(0:1,n,replace=TRUE)
  
  start_time[high_low==1] <- rnorm(length(which(high_low==1)),low_m,low_sd)
  start_time[high_low==0] <- rnorm(length(which(high_low==0)),high_m,high_sd)
  
  data <- mat.or.vec(n,num_tp)
  
  for (i in 1:n){
    
    data[i,] <- param[i,1] - exp((start_time[i] + (0:(num_tp-1))) * param[i,2]) + runif(num_tp,-noise,noise)
    
  }
  
  ## plot noisy data
  #plot(0:4,data[1,],type='l',ylim=c(0,1))
  
  #for (i in 2:1000){
  
  #lines(0:4,data[i,],col=i)
  
  #}
  
  ret1 <- vector('list',2)
  
  ret1[[1]] <- data.frame(data)
  ret1[[2]] <- start_time
  
  return(ret1)
  
}

# for exponential decline, add unit, ceiling and death effects
mmseStyle <- function(sim_data_exp){
  
  sim_mmse <- round(sim_data_exp)
  
  sim_mmse[sim_mmse>30] <- 30
  
  sim_mmse[sim_mmse<3] <- NA
  
  sim_mmse
  
}

# Monte Carlo point-wise average for generating curve, MMSE style data
monteCarloExp <- function(n=10^6,grid=30,modelLength = 30){
  # CAN BE COMPUTATIONALLY EXPENSIVE - I.E. MAY HAVE LONG RUNTIME
  #
  # may require shift to align with generating curve
  
  param <- simParamExp(n)
  
  full_data <- mat.or.vec(n,grid+1)
  
  for (i in 1:n){
    
    full_data[i,] <-  param[i,1] - exp(((0:grid)*(modelLength/grid))*param[i,2])
    
    print(i) # uncomment to see progress, but will slow progress
    
  }
  
  MCmeans_df <- data.frame(t=(0:grid)*(modelLength/grid),y=colMeans(full_data))
  
  MCmeans_df
  
}

# estimates midpoints and rates
slopes <- function(sim_data,mm=0){
  
  n <- dim(sim_data)[1]
  num_tp <- dim(sim_data)[2]
  
  rate <- mat.or.vec(n,1)
  mid <- mat.or.vec(n,1)
  
  # get rates and midpoints from linear models, or linear mixed models
  
  if (mm){
    
    # not updated for missing data
    
    simPanel <- matToPanel(sim_data)
    
    lmm_fit <- lmer(abeta ~ day + (day | person), simPanel)
    
    rate <- fixef(lmm_fit)[[2]] + ranef(lmm_fit)$person[,2]
    intercept <- fixef(lmm_fit)[[1]] + ranef(lmm_fit)$person[,1]
    
    for (i in 1:n){
      
      x <- 0:(num_tp-1)
      
      min_tp <- min(x)
      
      max_tp <- max(x)
      
      midpoint <- (max_tp-min_tp)/2
      
      mid[i] <- rate[i]*midpoint + intercept[i]
      
    }
    
    
    
  } else {
    
    for (i in 1:n){
      
      #print(i)
      
      y <- as.numeric(sim_data[i,])
      
      if (length(which(!is.na(y))) > 1) {
        
        #print(y)
        
        x <- 0:(num_tp-1)
        
        fit <- lm(y~x)
        
        rate[i] <- as.numeric(fit$coefficients[2])
        
        min_tp <- min(x[!is.na(y)])
        
        max_tp <- max(x[!is.na(y)])
        
        midpoint <- data.frame(x=min_tp+(max_tp-min_tp)/2)
        
        mid[i] <- predict(fit,midpoint)
        
      } else {
        
        mid[i] <- NA
        rate[i] <- NA
        
      }
      
    }
    
  }
  
  data.frame(mid,rate)
  
}


# set random seed to make reproducible
set.seed(1)


# Generate data
#

# generate parameters for simulation
mmse_param <- simParamExp(1000)

# generate untransformed data
sim_tmp <- simDataExp(mmse_param)

sim_data <- sim_tmp[[1]]
sim_starts <- sim_tmp[[2]]

# transform to make more MMSE-like
mmse_data <- mmseStyle(sim_data)

# Calculate monte carlo average
#

mmse_MCmeans <- monteCarloExp()


# Plots to illustrate simulation
#

pdf("Sim_shift.pdf",height=5,width=5)
plot(mmse_MCmeans[,'t'],mmse_MCmeans[,'y'],ylim=c(0,30),type='l',xlim=c(0,10),xlab="Years of cognitive decline",ylab="MMSE",lwd=2)

for (i in 300:350){

  lines((0:4) + sim_starts[i],mmse_data[i,],col=i) 
  
  
}

legend('bottomleft',"Simulation curve",lty=1,lwd=2)
dev.off()

pdf("Sim.pdf",height=5,width=5)
plot(0:4,mmse_data[300,],ylim=c(0,30),type='l',xlab="Years since first MMSE",ylab="MMSE")

for (i in 301:350){

  lines(0:4,mmse_data[i,],col=i) 
  
  
}

dev.off()

pdf("Sim_both.pdf",height=5,width=7)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))

plot(mmse_MCmeans[,'t'],mmse_MCmeans[,'y'],ylim=c(0,30),type='l',xlim=c(0,10),xlab="Years of cognitive decline",ylab="MMSE",lwd=2)

for (i in 300:350){
  
  lines((0:4) + sim_starts[i],mmse_data[i,],col=i) 
  
  
}

legend('bottomleft',"Average",lty=1,lwd=2)

plot(0:4,mmse_data[300,],ylim=c(0,30),type='l',xlab="Years since first MMSE",ylab="MMSE")

for (i in 301:350){
  
  lines(0:4,mmse_data[i,],col=i) 
  
}

dev.off()

mmse_slopes <- slopes(mmse_data)


Decline_parameter <- mmse_param[,2]
Estimated_slopes <- mmse_slopes[,2] 


First_MMSE <- mmse_data[,1]
df <- data.frame(Decline_parameter,Estimated_slopes,First_MMSE)
ggplot(data=df, aes(x=Decline_parameter,y=Estimated_slopes,col=First_MMSE)) + geom_point()
ggsave("Rate_sim_col.pdf",width = 5,height=5)


r2_table <- data.frame(multipleR2=NA,adjustedR2=NA,multipleR2_first=NA,adjustedR2_first=NA)

tmp <- lm(Decline_parameter~Estimated_slopes-1)

r2_table['Slope',1:2] <- c(summary(tmp)[[8]],summary(tmp)[[9]])

tmp <- lm(Decline_parameter~Estimated_slopes-1+First_MMSE)

r2_table['Slope',3:4] <- c(summary(tmp)[[8]],summary(tmp)[[9]])

