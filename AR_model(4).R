#計量経済学１

#OFFICE HOURS : 火曜日4-5, 3-312

rm(list = ls(all=T)) # This code clears all.

y <- vector(length=1000) #to save simulated data 벡터 공간 생성

for (i in 1:length(y)){
  y[i] <- rnorm(1, mean=0, sd=1) 
}
dev.new()
plot(y, type="l", xlab="period", col="blue")
legend("bottomleft", inset= 0.05,  "white noise process", lwd=2, lty=1, col="blue")

#------------AR(1) model: simulation exercise-----------------
# See how the parameter phi1 affects the dynamics of y. 

set.seed(123)
a <- rnorm(length(y), mean=0, sd=1) 

y[1] <- 0 #최초의 y값은 y1, y1 = 0으로 가정하는게 best guess, 
          #u을 0으로 가정했기떄문

phi1 <- 0.9
for (i in 2:length(y)){
  y[i] <- phi1*y[i-1] + a[i]
} #y1 = phi1*y0 + a[1]인데 y0은 data가 존재하지 않기때문에 
  #y1을 구할 수 없음

dev.new()
plot(y, type="l", xlab="period", col="skyblue")
legend("bottomleft", inset=0.05, c(paste("Phi1=", phi1)), lwd=2, lty=1, col="skyblue")

#Autocorrelation Function
dev.new()
acf(y, lag=12, main="ACF") #corr가 지속적으로 있으면 AR model
                          #corr가 중간에 끈기면 MA model

#Partial Autocorrelation Function
dev.new()
pacf(y, lag=12, main="PACF") #스파이크가 1개면 AR(1) model
                              #2개면 AR(2) model
#파랑색 밴드 밖이면 0이 아님


#-------Given various values of phi1, generate y---------------------
dev.new()
par(mfcol=c(2,2)) # 2 by 2
#options: mfcol and mfrow

#y[1] <- 0 # initial value 
for (phi1 in c(0.0, 0.5, 0.9, 1)){
  
  for (i in 2:length(y)){
    y[i] <- phi1*y[i-1] + a[i]
  }
  
  plot(y, type="l", xlab="period", col="blue")
  legend("bottomleft", inset=.05, c(paste("Phi1=", phi1)), lwd=2, lty=1, col="blue")
}
par(mfcol=c(1,1))

#------------AR(1) model with a constant term: simulation exercise-----------------
# See how the parameter phi0 affects the dynamics of y. 

phi0 <- 0.3
phi1 <- 1 #phi1=0.5 

#y[1] <- 0 # initial value
for (i in 2:length(y)){
  y[i] <- phi0 + phi1*y[i-1] + a[i]
}

dev.new()
plot(y, type="l", xlab="period", col="blue")
legend("bottomright", inset=0.05,  "non-stationary process", lwd=2, lty=1, col="blue")


#------Simulation of a random walk process: phi=1 ----------------------------------
#y[1] <- 0 # initial value 
for (j in 1:1000){     #number of simulations
  
  a <- rnorm(length(y), mean=0, sd=1) # Increase sd to see what happens.  
  for (i in 2:length(y)){
    y[i] <- 0.1+y[i-1] + a[i] #random walk with a drift
  }
  
  if (j==1) {
    dev.new()
    plot(y, type="l", col="blue", xlab="period", ylim=c(-100, 200), main="stock price")
  } else {
    lines(y, type="l", col="blue")
  }
  
}


#------------AR(2) model--------------------------------------------

dev.new()
par(mfcol=c(2,2)) # 2 by 2 

y[1]=0 # initial value
y[2]=0
for (phi in c(0.0, 0.5, 0.9, 1)){
  
  if (phi==0) {
    #white noise
    phi1=0  #phi=phi1+phi2
    phi2=0  
  }else if (phi==0.5) {
    phi1=0.25
    phi2=0.25
  }else if (phi==0.9) {
    phi1=0.45
    phi2=0.45
  } else {
    phi1=0.5
    phi2=0.5
    #두개 합이 1에 가까울수록 random work
  }
  
  for (i in 3:length(y)){
    y[i] <- phi1*y[i-1]+phi2*y[i-2] + a[i]
  }
  
  plot(y, type="l", xlab="period", col="blue")
  legend("bottomleft", inset=.05, c(paste("phi1+phi2=", phi)), lwd=2, lty=1, col="blue")
}
par(mfcol=c(1,1))




#---------------MA(1)-------------------------
y <- vector(length=50)
set.seed(123)
a <- rnorm(length(y), mean=0, sd=1) 

dev.new()
par(mfcol=c(1,2)) # 2 by 2 

y[1]=0 # first value 
for (theta in c(0.8, -0.8)){
  
  for (i in 2:length(y)){
    y[i] <- a[i]+theta*a[i-1] # case mean = 0 
  }
  
  plot(y, type="l", xlab="period", col="blue")
  legend("topleft", inset=.05, c(paste("theta=", theta)), lwd=2, lty=1, col="blue")
}
par(mfcol=c(1,1))

#y[t] = a[t] + 0.8a[t-1] vs y[t] = a[t]  - 0.8a[t-1]
#in the case corr(y[t], y[t-1]) is positive and negative 
# theta access to zero, which become white noise.


 

