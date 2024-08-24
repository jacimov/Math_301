####Lab2
library(deSolve)
library(minpack.lm)
library(tidyverse)
#########################Exercise1############################
ex1ODE = function(t,y,parms){
  list(c(y[1]-y[2], 1-exp(y[1])))
}

ex1sys = ode(y=c(1.5, 3), func = ex1ODE, times = seq(0,3,0.1), parms=NA)
plot(ex1sys)

##graph
time = ex1sys[,1] # The first column in ex1Sys is t
x = ex1sys[,2]    # The second column holds x
y = ex1sys[,3]    # The last column holds y
plot(time,x,type = "l", lty = 1,ylim = range(x,y), xlab = "Time",ylab = "x,y",col = "red")   # Plot x versus time
par(new=TRUE)   # Tell R you want to plot on top of an existing plot
plot(time,y,type = "l",lty = 1,ylim = range(x,y),
     xlab = "",ylab = "",col = "blue")   # Plot y versus time
legend("right",legend=c("x", "y"),col=c("red","blue"),lty=1)   # Add a legend


####################Exercise2######################
Compmodel = function(t,y,parms){
  with(as.list(parms),{
    return(list(c(a*y[1]-b*y[1]*y[2], m*y[2]-n*y[1]*y[2])))
  })
}

parValues = c(a=0.2, b=0.1, m=0.3, n=0.1)
ex2sys = ode(y=c(1.5, 3), func = ex1ODE, times = seq(0,12,0.1), parms=NA)
time = ex2sys[,1] # The first column in ex1Sys is t
x = ex2sys[,2]    # The second column holds x
y = ex2sys[,3]    # The last column holds y
plot(time,x,type = "l", lty = 1,ylim = range(x,y), xlab = "Time",ylab = "x,y",col = "red")   # Plot x versus time
par(new=TRUE)   # Tell R you want to plot on top of an existing plot
plot(time,y,type = "l",lty = 1,ylim = range(x,y),
     xlab = "",ylab = "",col = "blue")   # Plot y versus time
legend("right",legend=c("x", "y"),col=c("red","blue"),lty=1)   # Add a legend

#######################Exercise3###################
tData = 0:12
wolfData = c(2, 2, 2.5, 3, 3, 3.5, 4, 4.5, 5, 6.5, 9, 13, 20)

residFn=function(parValues){
  # initial populations
  yinit=c(2, 3)
  # time points to use when solving the ODE
  timespan=seq(0,max(tData),.05)
  # parameters from the parameter estimation routine (only a for this case)
  a=parValues[1]
  # fixed values for b,m,n
  bfix = 0.1
  mfix = 0.3
  nfix = 0.1
  parms=c(a,b=bfix,m=mfix,n=nfix)
  # solve ODE for a given set of parameters
  output=ode(y=yinit,times=timespan,func=Compmodel,parms=parms)
  
  ## Filter data that contains time points where data is available
  outdf=data.frame(output[,1],output[,2]) #put time, wolf pop into a data frame
  outdf=outdf[outdf$output...1 %in% tData,] #only include the time points from the dataset
  
  residFn = outdf$output...2 - wolfData
  
  # return predicted vs experimental residual
  return(residFn)
}

parValues = c(a=0.2)
plot(tData,wolfData)
fitval <- nls.lm(par = parValues, lower = NULL, upper = NULL, fn = residFn, jac = NULL)
a <- fitval$par
nums <- c(a = 0.45,b = .1, m = .3, n = .1)

ex3sys <- ode(y=c(2, 3), func = Compmodel, times = seq(0,12,1), parms = nums)
plot(ex3sys)