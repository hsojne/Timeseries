rm(list = ls(all=T))
graphics.off()

getwd()

#============================================================================
#
#   Estimation of an ARCH(1) Model by Maximum Likelihood
#
#============================================================================
library(quantmod) # Load the package
#library(tseries) # for garch function 

###################################################################################
#------------------SP500-----------------------------------------------------------
###################################################################################
SP500<-getSymbols("^GSPC",from="1990-01-05", to="2024-04-01", auto.assign = FALSE)
tail(SP500)
#^GSPC
#compute monthly returns
SP500.return=dailyReturn(SP500, type='arithmetic') #산술평균

mean(annualReturn(SP500, type='arithmetic')) #average annual return

#Eliminate date from the data
date <-index(SP500.return)
SP500.return=coredata(SP500.return)
SP500.return <-SP500.return-mean(SP500.return) #demeaned returns, it is 'rt' 

#----------------------------------------------------------------------------
# Likelihood function for a ARCH(1) model
#----------------------------------------------------------------------------
likeli_arch <- function( p,y ) {
  #y = data
  alpha <- abs(p)
  n <- length( y )
  e <- y  #y=SP500.return
  sigma2 <- sd(y)^2*rep(1, n)  #sigma2 is set to the variance of y=SP500.return
  
  for (i in 2:n) {
    sigma2[i] <- alpha[1] + alpha[2]*e[i-1]^2 
    #sigma2[i] <- exp(alpha[1]) + alpha[2]*e[i-1]^2 #if alpha[1]<0
  }
  
  likeli <-0
  for (j in 2:n) {
    likeli <- likeli + (- 0.5*log( 2*pi ) - 0.5*log( sigma2[j] ) - 0.5*(e[j]^2/sigma2[j]) ) 
  }
  minus_like <- -likeli
  return(minus_like)
}


#--------------------------------------------------------------
# Estimating the ARCH(1) model
p <- c(1, 1)
output <- optim(p, likeli_arch, y=SP500.return, method="Nelder-Mead", hessian=TRUE) 
                           #Algorithms: ??Nelder-Mead??, ??BFGS??, ??CG??, ??L-BFGS-B??, ??SANN??, ??Brent??
theta <- output$par
theta 

likelihood <- -output$val
likelihood

information <- output$hessian #Fisher Information: I=-(H)         
var_cov <- solve(information) #Variance-Covariance Matrix
var_cov

se <-sqrt(diag(var_cov))  #The function diag takes the diagnal elements to compute se. 
se

t <-output$par/se
t
#print the estimates of the model parameters, standard errors, and t-values
summary.output<-data.frame(output$par, se, t)
#summary.output<-round(summary.output, digits = 4)
colnames(summary.output)<-c("parameter", "standard error", "t-value")
rownames(summary.output)<-c("alpha0", "alpha1")
summary.output

#-------------------------------------------------------------------
# Use the built-in function garch to estimate the ARCH(1) model 
#install.packages("tseries")
library(tseries)
arch1=garch(SP500.return, c(0,1)) #arch1
arch1$coef #The same results are obtained from the garch function. 
#-------------------------------------------------------------------

# Compute the conditional variance using the estimated parameters
y<-SP500.return
sd(y)

n <-length(y)
sigma2 <- sd( y )^2*rep(1,n)
e <- y 

for (i in 2:n) {
  sigma2[i] <- theta[1] + theta[2]*e[i-1]^2 
}

#----------------graph of the conditional variance-----------------------------
dev.new()  
par(mfrow=c(2,1), mar=c(5.0,5.0,1.0,1.0)) #set the bottom, left, top and right margins respectively 
                                          #of the plot region in number of lines of text. 
plot(date, sigma2,type="l", xlab="date", ylab=expression(sigma[t]^2), col="blue")
plot(date, SP500.return,type="l", xlab="date", ylab="return", col="blue")

