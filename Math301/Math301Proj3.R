##Math proj 3
###########packages##############
library(deSolve)
library(minpack.lm)
library(tidyverse)
library(readr)

##############Data##############
data <- read_csv("Desktop/infection_data.csv")
Data2 <- read_csv("Desktop/Data2.csv")
data <- infection_data
plot(x=data$Day,y=data$Infections)

########Functions#############
SIRmodel = function(t,y,parms){
  with(as.list(parms),{
    return(list(c(-r*y[1]*y[2], r*y[1]*y[2]-a*y[2], a*y[2])))
  })
}

ex1ODE = function(t,y,parms){
  list(c(y[1]-y[2], 1-exp(y[1])))
}

residFn=function(parValues){
  # initial populations
  yinit=c(762, 1,0)
  # time points to use when solving the ODE
  timespan=seq(0,max(data$Day),0.05)
  # parameters from the parameter estimation routine (only a for this case)
  r=parValues[1]
  # fixed values for b,m,n
  afix = 0.4761905
  parms=c(r,a = afix)
  # solve ODE for a given set of parameters
  output=ode(y=yinit,times=timespan,func=SIRmodel,parms=parms)
  
  ## Filter data that contains time points where data is available
  outdf=data.frame(output[,1],output[,3]) #put time, wolf pop into a data frame
  outdf=outdf[outdf$output...1 %in% data$Day,] #only include the time points from the dataset
  
  residFn = outdf$output...3 - data$Infections
  
  # return predicted vs experimental residual
  return(residFn)
}

########Solving#########
## s(0)=763
## 1/a = 2.1 --> a = 0.476190476190476
## r = 0.918787406326296 (from fitval below)
##fit data
parValues <- c(r = .001)
fitval = nls.lm(par = parValues, fn = residFn)

r <- fitval$par[1]
nums <- c(r = fitval$par[1], a = 0.47619)

sys <- ode(y=c(762,1,0), func = SIRmodel, times = seq(0,13,.1), parms = nums)
plot(sys)

#########Applying#########
## Fix model parameters ##
# Reaction rates
rr=0.00221
ra=0.4761904

# Initial populations 
S0 = 762
I0 = 1
R0 = 0

# Time to stop a simulation
maxt = 17
# Number of simulations
nsims = 100


## Start Simulations ##
# Store the populations at the end of each simulation
popEnd = matrix(0,nsims,3)

# Initialize a plot to show AB molecules over time
plot(NULL, xlim=c(0,maxt), ylim=c(0,800), ylab="Number of Infected People", xlab="Time")
par(ask=TRUE)

# Run through nsims simulations (FINISH FOR LOOP CONDITION)
for (i in 1:nsims ){
  t=0  #time variable
  #Initialize population counters (FILL IN BELOW)
  S=S0
  I=I0 
  R=R0
  #Track populations over time for a single simulation
  popMat = c(t,S,I,R)
  
  # Run through a single simulation (lasts until maxt)
  # FILL IN WHILE LOOP CONDITION
  while(t < maxt && I > 0){
    # FILL IN ALL INFORMATION BELOW
    # Calculate sum of current rates 
    rateSum = rr*S*I+ra*I
    # two uniformly distributed random numbers
    u = runif(2)
    # increase time by exponentially distributed RV tau
    tau = -log(u[1])/rateSum
    t = t + tau
    # determine which event occurs
    if (u[2] < (rr*S*I)/rateSum ){ #forward reaction occurs
      S = S-1
      I = I+1
    }else{ #backward reaction occurs
      R = R+1
      I = I-1
    }
    # add new row to matrix with this population change:
    popMat = rbind(popMat, c(t,S,I,R))
  }
  # record ending populations
  popEnd[i,] = popMat[nrow(popMat),2:4]
  # plot the number of AB molecules over time
  lines(popMat[,1],popMat[,3],col=sample(rainbow(10)))
}

par(ask=FALSE) #Stop adding to the current plot

write.csv(popMat, file="Desktop/popmat.csv")
## Plot histograms ##
# Run the following after you have saved the plot above
# Plot a histogram of the A,B, AB populations at the end of the simulations
###Plots
#

# INCLUDE A LINE TO ADD A HISTOGRAM FOR THE AB POPULATION IN PURPLE
legend( x= 15, y=25,c("S","I","R"),col=c("red", "blue","purple"),lwd=5)


####extensions
########Functions#############
SIRmodel = function(t,y,parms){
  with(as.list(parms),{
    return(list(c(-r*y[1]*y[2], r*y[1]*y[2]-a*y[2], a*y[2])))
  })
}

ex1ODE = function(t,y,parms){
  list(c(y[1]-y[2], 1-exp(y[1])))
}

residFn=function(parValues){
  # initial populations
  yinit=c(8500000, 20,0)
  # time points to use when solving the ODE
  timespan=seq(0,max(Data2$Day),0.05)
  # parameters from the parameter estimation routine (only a for this case)
  r=parValues[1]
  # fixed values for b,m,n
  afix = 0.03571429
  parms=c(r,a = afix)
  # solve ODE for a given set of parameters
  output=ode(y=yinit,times=timespan,func=SIRmodel,parms=parms)
  
  ## Filter data that contains time points where data is available
  outdf=data.frame(output[,1],output[,3]) #put time, wolf pop into a data frame
  outdf=outdf[outdf$output...1 %in% Data2$Day,] #only include the time points from the dataset
  
  residFn = outdf$output...3 - Data2$I
  
  # return predicted vs experimental residual
  return(residFn)
}

########Solving#########
## s(0)=763
## 1/a = 2.1 --> a = 0.476190476190476
## r = 0.918787406326296 (from fitval below)
##fit data
parValues <- c(r = .65/6000000)
fitval2 = nls.lm(par = parValues, fn = residFn)

r <- fitval2$par[1]
nums2 <- c(r = fitval2$par[1], a = .2)
sys2 <- ode(y=c(8500000,20,0), func = SIRmodel, times = seq(0,180,.01), parms = nums2)
plot(sys2)

## Fix model parameters ##
# Reaction rates
rr=fitval2$par[1]
ra=0.2

# Initial populations 
S0 = 8500000
I0 = 20
R0 = 0

# Time to stop a simulation
maxt = 150
# Number of simulations
nsims = 3


## Start Simulations ##
# Store the populations at the end of each simulation
popEnd = matrix(0,nsims,3)

# Initialize a plot to show AB molecules over time
plot(NULL, xlim=c(0,150), ylim=c(0,8750000), ylab="Number of Infected People", xlab="Time")
par(ask=TRUE)

# Run through nsims simulations (FINISH FOR LOOP CONDITION)
for (i in 1:nsims ){
  t=0  #time variable
  #Initialize population counters (FILL IN BELOW)
  S=S0
  I=I0 
  R=R0
  #Track populations over time for a single simulation
  popMat = c(t,S,I,R)
  
  # Run through a single simulation (lasts until maxt)
  # FILL IN WHILE LOOP CONDITION
  while(t < maxt && I > 0){
    # FILL IN ALL INFORMATION BELOW
    # Calculate sum of current rates 
    rateSum = rr*S*I+ra*I
    # two uniformly distributed random numbers
    u = runif(2)
    # increase time by exponentially distributed RV tau
    tau = -log(u[1])/rateSum
    t = t + tau
    # determine which event occurs
    if (u[2] < (rr*S*I)/rateSum ){ #forward reaction occurs
      S = S-1
      I = I+1
    }else{ #backward reaction occurs
      R = R+1
      I = I-1
    }
    # add new row to matrix with this population change:
    popMat = rbind(popMat, c(t,S,I,R))
  }
  # record ending populations
  popEnd[i,] = popMat[nrow(popMat),2:4]
  # plot the number of AB molecules over time
  lines(popMat[,1],popMat[,3],col=sample(rainbow(10)))
}

par(ask=FALSE)

