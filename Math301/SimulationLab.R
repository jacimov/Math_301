#####SimulationLab
######start########
a=4 ##when a is 5, prints the second if
if( a<5 ){
  print("a is less than 5")
}
if( a<6 && a>5 ){
  print("a is a number between 5 and 6.")
}

##########Exercise1###########

if( a<6 && a>5 ){
  print("a is a number between 5 and 6.")
} else {
  print("a is not between 5 and 6.")
}
x=3.5
if( x<1 ){
  print("y=1")
} else if (x<2){
  print("y=2")
} else if (x<4){
  print("y=3")
} else{
  print("y=4")
}
## x=3.5,y=3 ; x=5,y=4

########WhileLoops#########
a=2
while( a<100000 ){
  a=a^(2)
}
exercise2 <- sqrt(a) ## 65536

##########ForLoops############
a=2
for( i in 1:4 ){
  a=a^2
}
a=2
for( i in 1:4 ){
  a=a + i
}

exercise3 <- for( i in 1:8){
  i = i^2
  print(i)
}

##############GeneratingRandomNumbers###############
x=runif(1)
x=runif(10)
x=runif(1,3,7)

##Exercise4
r = runif(1)
if (  r    < 1/6                ) {
  print("y = 1")
}    else if (1/6<r && r<(2/6)                ) {
  print("y = 2")
}    else if ( 2/6<r && r  < 3/6                 ) {
  print("y = 3")
}    else if (  3/6 < r && r  < 4/6                  ) {
  print("y = 4")
}    else if ( 4/6 < r && r < 5/6                   ) {
  print("y = 5")
}    else  {
  print("y = 6")
}   

##Exercise5
N=100 #number of trials

m = 10 # number of flips per trial
keeptrack = matrix(0,1,N)

for (i in 1:N){
  counter = 0 # initialize counter
  for (j in 1:m){
    x=runif(1)
    if (x <= 0.5){counter = counter+1} #increment counter if heads (x<= 0.5)
  }
  keeptrack[i] = counter #store current trial info
}
mean(keeptrack) #find mean number of heads over all trials
hist(keeptrack) #plot a histogram of the results

##Modiefied
N=100 #number of trials

m = 10 # number of flips per trial
keeptrack = matrix(0,1,N)

for (i in 1:N){
  counter = 0 # initialize counter
  for (j in 1:m){
    x=runif(1)
    if (x <= 0.35){counter = counter+1} #increment counter if heads (x<= 0.5)
  }
  keeptrack[i] = counter #store current trial info
}
mean(keeptrack) #find mean number of heads over all trials
hist(keeptrack) #plot a histogram of the results

#########################GillespieCode###################

# Gillespie Algorithm Example 

# This is a code outline for Exercise 6 from the Simulation Lab.
# Fill in the missing pieces of the code to simulate 
# the chemical reactions.

## Fix model parameters ##
# Reaction rates
rf=0.4
rb=0.25

# Initial populations 
N_A0 = 25
N_B0 = 35
N_AB0 = 5

# Time to stop a simulation
maxt = 1
# Number of simulations
nsims = 100


## Start Simulations ##
# Store the populations at the end of each simulation
popEnd = matrix(0,nsims,3)

# Initialize a plot to show AB molecules over time
plot(NULL, xlim=c(0,maxt), ylim=c(0,40), ylab="Number of AB molecules", xlab="Time")
par(ask=TRUE)

# Run through nsims simulations (FINISH FOR LOOP CONDITION)
for (i in 1:nsims ){
  t=0  #time variable
  #Initialize population counters (FILL IN BELOW)
  A=N_A0
    B=N_B0 
    AB= N_AB0
    #Track populations over time for a single simulation
    popMat = c(t,A,B,AB)
  
  # Run through a single simulation (lasts until maxt)
  # FILL IN WHILE LOOP CONDITION
  while(t < maxt ){
    # FILL IN ALL INFORMATION BELOW
    # Calculate sum of current rates 
    rateSum = rf*A*B+rb*AB
      # two uniformly distributed random numbers
      u = runif(2)
      # increase time by exponentially distributed RV tau
      tau = -log(u[1])/rateSum
    t = t + tau
      # determine which event occurs
      if (u[2] < (rf*A*B)/rateSum ){ #forward reaction occurs
        A = A-1
          B = B-1
          AB = AB+1
      }else{ #backward reaction occurs
        A = A+1
          B = B+1
          AB = AB-1
      }
    # add new row to matrix with this population change:
    popMat = rbind(popMat, c(t,A,B,AB))
  }
  # record ending populations
  popEnd[i,] = popMat[nrow(popMat),2:4]
  # plot the number of AB molecules over time
  lines(popMat[,1],popMat[,4],col=sample(rainbow(10)))
}

par(ask=FALSE) #Stop adding to the current plot

## Plot histograms ##
# Run the following after you have saved the plot above
# Plot a histogram of the A,B, AB populations at the end of the simulations
hist(popEnd[,1],main = "Final Populations",col="red",xlim=c(0,30),xlab="Population")
hist(popEnd[,2],add=TRUE,col="blue")
hist(popEnd[,3],add=TRUE, col="purple")
# INCLUDE A LINE TO ADD A HISTOGRAM FOR THE AB POPULATION IN PURPLE
legend( x= 15, y=25,c("A","B","AB"),col=c("red", "blue","purple"),lwd=5)

