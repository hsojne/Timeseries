rm(list = ls(all=T))
graphics.off()

getwd()
setwd("F:\\teaching\\Time_Series_Analysis_2023\\R_code_examples_2023")

#============================================================================
#   Estimation of an GARCH Model by Maximum Likelihood
#============================================================================
#install.packages("quantmod")
library(quantmod) # Load the package
#install.packages("tseries")
library(tseries) # for garch function 




###################################################################################
#------------------SP500-----------------------------------------------------------
###################################################################################

SP500<-getSymbols("^GSPC",from="1990-01-03", to="2024-04-01", auto.assign = FALSE)
#SP500<-getSymbols("^GSPC",from="1990-01-01", to="2023-09-30",periodicity = 'monthly', auto.assign = FALSE)

#compute daily/monthly returns
SP500.return=dailyReturn(SP500, type='arithmetic')
#SP500.return=monthlyReturn(SP500, type='arithmetic')

#Eliminate date from the data
date <-index(SP500.return)
SP500.return <-coredata(SP500.return)
SP500.return <-SP500.return-mean(SP500.return)


#----------------------------------------------------------------------------
# Likelihood function for the GARCH(1,1) model
#----------------------------------------------------------------------------
likeli_garch <- function( p,y ) {
  alpha <- abs(p)
  n <- length( y )
  e <- y  #y=SP500.return
  sigma2garch <- sd(y)^2*rep(1, n) #sigma2garch=sigma^2
  
  for (i in 2:n) {
    sigma2garch[i] <- alpha[1] + alpha[2]*e[i-1]^2 + alpha[3]*sigma2garch[i-1] # Notice that sigma2garch[1] is set to sd(y)^2. 
  }
  
  likeli <-0
  for (j in 2:n) {
    likeli <- likeli - 0.5*log( 2*pi ) - 0.5*log( sigma2garch[j] ) - 0.5*(e[j]^2/sigma2garch[j]) 
  }
  minus_like <- -likeli
  return(minus_like)
}


#--------------------------------------------------------------
# Estimating the GARCH(1) model
p <- c(0.30, 0.30, 0.40)
output <- optim(p, likeli_garch, y=SP500.return, method="Nelder-Mead", hessian=TRUE) 
                            #Algorithms: ??Nelder-Mead??, ??BFGS??, ??CG??, ??L-BFGS-B??, ??SANN??, ??Brent??
theta <- output$par
theta 

likelihood <- output$val
likelihood <- -likelihood
likelihood 

information <- output$hessian #Fisher Information          
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
rownames(summary.output)<-c("alpha0", "alpha1", "beta1")
summary.output 


#-------------------------------------------------------------------
# Use the built-in function garch to estimate the GARCH(1) model 
garch1=garch(SP500.return, c(1,1)) #,trace=F
out<-summary(garch1)
out$coef

#-------------------------------------------------------------------


# Compute the conditional variance using the estimated parameters
y<-SP500.return
sd(y)

n <-length(y)
sigma2garch <- sd( y )^2*rep(1,n)
e <- y 

for (i in 2:n) {
  sigma2garch[i] <- theta[1] + theta[2]*e[i-1]^2 + theta[3]*sigma2garch[i-1]
}

#----------------graph for the conditional variance-----------------------------
dev.new()  
par(mfrow=c(2,1), mar=c(5.0,5.0,1.0,1.0)) #sets the bottom, left, top and right margins respectively 
                            #of the plot region in number of lines of text. 
plot(date, SP500.return,type="l", xlab="date", ylab="Stock return", col="blue")
plot(date, sigma2garch, type="l", xlab="date", ylab=expression(sigma[t]^2), col="blue")


save(sigma2garch, file="sigma2garch.RData")

