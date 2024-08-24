rm(list = ls(all=T))
graphics.off()

getwd()
setwd("F:\\teaching\\Time_Series_Analysis_2023\\R_code_examples_2023")

load("F:\\teaching\\Time_Series_Analysis_2023\\R_code_examples_2023\\sigma2garch.RData")

#============================================================================
#
#   Estimation of an GARCH Model by Maximum Likelihood
#
#============================================================================
library(quantmod) # Load the package
#install.packages("tseries")ㅌ
#library(tseries) # for garch function 
#install.packages("rugarch")
#library(rugarch) # for ugarchspec


###################################################################################
#------------------SP500-----------------------------------------------------------
###################################################################################

SP500<-getSymbols("^GSPC",from="1990-01-03", to="2024-04-01", auto.assign = FALSE)
#MSFT<-getSymbols("MSFT",from="1990-01-05", to="2024-04-01", auto.assign = FALSE)

#compute daily returns
SP500.return=dailyReturn(SP500, type='arithmetic')
#SP500.return=monthlyReturn(SP500, type='arithmetic')

#Eliminate date from the data
date <-index(SP500.return)
SP500.return <-coredata(SP500.return)
SP500.return <-SP500.return-mean(SP500.return)


#----------------------------------------------------------------------------
# Likelihood function for a TGARCH(1,1) model
#----------------------------------------------------------------------------
likeli_tgarch <- function( p,y ) {
  alpha <- abs(p)
  n <- length( y )
  e <- y  #y=SP500.return
  sigma2 <- sd(y)^2*rep(1, n) #sigma2=sigma^2
  
  for (i in 2:n) {
    
    if (e[i-1]>=0) {
      sigma2[i] <- alpha[1] + alpha[2]*e[i-1]^2 + alpha[3]*sigma2[i-1] 
    } else {
      sigma2[i] <- alpha[1] + (alpha[2]+alpha[4])*e[i-1]^2 + alpha[3]*sigma2[i-1]
    }
    
  }
  
  likeli <-0
  for (j in 2:n) {
    likeli <- likeli - 0.5*log( 2*pi ) - 0.5*log( sigma2[j] ) - 0.5*(e[j]^2/sigma2[j]) 
  }
  minus_like <- -likeli
  return(minus_like)
}


#--------------------------------------------------------------
# Estimating the TGARCH(1) model
p <- c(0.50, 0.90, 0.90, 0.2)
output <- optim(p, likeli_tgarch, y=SP500.return, method="Nelder-Mead", hessian=TRUE) 
#Algorithms: ??Nelder-Mead??, ??BFGS??, ??CG??, ??L-BFGS-B??, ??SANN??, ??Brent??
theta <- output$par
theta 

likelihood <- -output$val
likelihood

information <- output$hessian #Fisher Information          
var_cov <- solve(information) #Variance-Covariance Matrix
var_cov

se <-sqrt(diag(var_cov))  #The function diag takes the diagnal elements to compute se. 
se

t <-output$par/se
t

# Compute conditional variance at optimal parameters
y<-SP500.return
n <-length(y)

e <- y 
sigma2 <- sd( e )^2*rep(1,n)


for (i in 2:n) {
  
  if (e[i-1]>=0) {
    sigma2[i] <- theta[1] + theta[2]*e[i-1]^2 + theta[3]*sigma2[i-1] 
  } else {
    sigma2[i] <- theta[1] + (theta[2]+theta[4])*e[i-1]^2 + theta[3]*sigma2[i-1]
  }
}
head(sigma2)
#----------------Figure for the estimated conditional variance-----------------------------
dev.new()  
par(mar=c(5.0,5.0,1.0,1.0)) #sets the bottom, left, top and right margins respectively 
                            #of the plot region in number of lines of text. 
plot(date, sigma2, type="l", xlab="date", ylab=expression(sigma[t]^2), col="blue")


#----------------comparison of GARCH and TGARCH--------------------------------------------
load("F:\\teaching\\Time_Series_Analysis_2023\\R_code_examples_2023\\sigma2garch.RData")

dev.new()  
par(mar=c(5.0,5.0,1.0,1.0))  
plot(date, sigma2, type="l", xlab="date", ylab=expression(sigma[t]^2), col="blue")
lines(date, sigma2garch, type="l", xlab="date", ylab=expression(sigma[t]^2), col="red")
text(date[n-2000], 0.0035, labels="TGARCH", col="blue")
text(date[n-2000], 0.0033, labels="GARCH", col="red")
  

