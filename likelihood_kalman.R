#log-likelihood function 
likelihood <- function(b,y) {
  
  # system matrices: Lambda, Phi, R, and Q
  
  Lambda  <- matrix(c(1, 1, 0, 0), nrow = 1)
  
  Phi  <- matrix(c(1,      0,     0,   1,
                   0,   b[1],  b[2],   0,
                   0,      1,     0,   0,
                   0,      0,     0,   1),    4,4, byrow=T)  
  R <- 0 
  
  Q <- matrix(c(b[3]^2,      0,       0,       0,
                0,      b[4]^2,       0,       0,
                0,           0,       0,       0,
                0,           0,       0, b[5]^2),  4,4, byrow=T)
  
  likeli <- -sum( kalman.filter(y,Phi,Lambda,R,Q) ) #minus log-likelihood
  return (likeli)
}

#-------------------------------------------------------------------
# Kalman filter algorithm
#-------------------------------------------------------------------
kalman.filter <- function(y,Phi,Lambda,R,Q) {
  
  t <- nrow(y)
  n <- ncol(y)           #number of variables
  k <- nrow(Q)           #number of state variables
  lnl <- rep(0, t)       #for log-likelihood values
  
  #Recursions of the Kalman Filter
  
  #Initial values
  st <- matrix(rep(0, k), nrow = k) #k=4,   #s{1|0}
  st[1,1]<-y[1]
  
  pt <- diag(k)
  pt[1,1]<-sd(diff(y))^2 #guess             #P{1|0}
  pt[2,2]<-sd(diff(y))^2
  pt[4,4]<-0.10*sd(diff(y))^2
  
  #Prediction of yt, its error and variance
  mt <- Lambda %*% st                       #y{1|0}  
  vt <- Lambda %*% pt %*% t(Lambda) + R     #V{1|0}    
  ut <- y[1,] - mt                          #prediction error
  
  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5* t(ut) %*% solve(vt) %*% (ut)
  
  #Kalman gain and updating of s and p
  Kal <- pt %*% t(Lambda) %*% solve(vt)     #Kalman gain
  s0 <- st + Kal %*% ut                     #s{1|1} 
  p0 <- pt - Kal %*% Lambda %*% pt          #P{1|1}
  
  
  # Loop to compute s{t|t}, s{t|t-1}, y{t|t}, y{t|t-1}, P{t|t}, and P{t|t-1} 
  for (i in 2:t) {
    
    # Prediction of st and pt
    st <- Phi %*% s0                           #s{t|t-1}
    pt <- Phi %*% p0 %*% t(Phi) + Q            #P{t|t-1}
    
    # Prediction of yt, its error and variance 
    mt <- Lambda %*% st                        #y{t|t-1}
    vt <- Lambda %*% pt %*% t(Lambda) + R      #V{t|t-1}
    ut <- y[i,] - mt                           #u{t|t-1}
    
    # Log-likelihood function
    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5* t(ut) %*% solve(vt) %*% (ut)
    
    # Updating equation        	
    Kal <- pt %*% t(Lambda) %*% solve(vt)      #Kalman gain
    s0 <- st + Kal %*% ut                      #s{t|t}
    p0 <- pt - Kal %*% Lambda %*% pt           #P{t|t}
  }    
  return(lnl)
}


#------------------------------------------------------------------
#       Smoothing
#------------------------------------------------------------------
smoothing <- function( y, para ) {
  
  phi1   <- para[1]
  phi2   <- para[2]
  sigtau <- para[3]
  sigc   <- para[4]
  sigmu  <- para[5]
  
  Lambda  <- matrix(c(1, 1, 0, 0), nrow = 1)
  
  Phi  <- matrix(c(1,      0,     0,   1,
                   0,   phi1,  phi2,   0,
                   0,      1,     0,   0,
                   0,      0,     0,   1),    4,4, byrow=T)  
  
  R <- 0 
  
  Q <- matrix(c(sigtau^2,      0,       0,        0,
                0,        sigc^2,       0,        0, 
                0,             0,       0,        0,
                0,             0,       0,   sigmu^2), 4,4, byrow=T)
  
  # Arrays to save s and p 
  t <- nrow(y)
  n <- ncol(y)
  
  lnl     <- rep(0, t)
  k       <- nrow(Q)
  
  s10     <- array(0, c(t,k))    #To save s{t|t-1}
  s11     <- array(0, c(t,k))    #To save s{t|t}  
  p10     <- array(0, c(k,k,t))  #To save P{t|t-1}
  p11     <- array(0, c(k,k,t))  #To save P{t|t}
  ss      <- array(0, c(t,k))    #To save s{t|T}
  
  # Initial values
  st <- matrix(rep(0, k),nrow=k) #s{1|0}
  st[1,1]<-y[1]
  pt <- diag(k)                  #P{1|0}
  pt[1,1]<-sd(diff(y))^2
  pt[2,2]<-sd(diff(y))^2
  pt[4,4]<-0.10*sd(diff(y))^2
  
  s10[1,] <- t(st)               #Save s{1|0}   #s10[1,] <- st  is also okay
  p10[,,1] <- pt                 #Save P{1|0}
  
  mt <- Lambda %*% st                    #y{1|0}
  vt <- Lambda %*% pt %*% t(Lambda) + R  #V{1|0}
  ut <- y[1,] - mt                       #u{1|0}
  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5* t(ut) %*% solve(vt) %*% (ut)
  
  Kal <- pt %*% t(Lambda) %*% solve(vt)  #Kalman gain
  s0 <- st + Kal %*% ut                  #s{1|1}
  p0 <- pt - Kal %*% Lambda %*% pt       #P{1|1}
  s11[1,] <- t(s0)                       #Save s{1|1}
  p11[,,1] <- p0                         #Save P{1|1}
  
  
  # Main loop over observations  
  for (i in 2:t) {
    # Prediction equations
    st <- Phi %*% s0                      #s{t|t-1}
    pt <- Phi %*% p0 %*% t(Phi) + Q       #P{t|t-1}
    s10[i,] <-  t(st)                     #Save s{t|t-1}
    p10[,,i] <- pt                        #Save P{t|t-1}
    
    # Observation equations
    mt <- Lambda %*% st                   #y{t|t-1}
    vt <- Lambda %*% pt %*% t(Lambda) + R #V{t|t-1}
    ut <- y[i,] - mt                      #u{t|t-1}
    
    # Log-likelihood function
    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*t(ut) %*% solve(vt) %*% (ut)
    
    # Updating equations         
    Kal <- pt %*% t(Lambda) %*% solve(vt) #Kalman gain
    s0 <- st + Kal %*% ut                 #s{t|t}
    p0 <- pt - Kal %*% Lambda %*% pt      #P{t|t}
    s11[i,] <- t(s0)                      #Save s{t|t}
    p11[,,i] <- p0                        #Save P{t|t}
  }
  
  # Smoothing    
  ss <- s11
  for (j in 1:(t-1)) {    
    Jt      <- p11[,,t-j] %*% t(Phi) %*% solve(p10[,,t-j])
    ss[t-j,] <- s11[(t-j),] + Jt %*% (ss[(t-j+1),] - s10[(t-j+1),])
  }  
  likeli   <- -sum(lnl) 
  return(list(ss=ss, lf=likeli))
}
