#Function to run the Kalman filter on simple C model
#CTS 28 July 2014
#This uses the same code as kalmanStep, but instead of returning NLL it returns
#estimates of the state at each time step

#Arguments are
#  p - Parameters to be estimated
#  Y - Observations of state variables. Rows are time points, columns are (DIC, DOC)
#  H - Estimate of the observation error covariance matrix
#  otherPar - Additional fixed parameters, not to be estimated, needed for calculations
#  drivers - A data.frame with the driver data needed for the ut matrix.

kalmanState<-function(p,Y,H,drivers){
  require(abind)
  with(drivers,{
    #Unpack estimated parameters
    r <- p[1]
    iota <- p[2]
    k <- p[3]
    Q <- matrix(c(exp(p[4]),0,0,exp(p[5])),nrow=2,byrow=TRUE)	#process error covariance matrix [2x2]
    
    #Initialize P (covariance matrix for estimate of state) from Q
    P <- Q
    #Initialize POut to hold P estimate at each time point
    POut <- P
      
    #Number of observations
    nObs <- nrow(Y)
    
    #Number of state variables
    nVar <- 2
    
    #Initialize a (estimate of state)
    a <- matrix(as.numeric(Y[1,]),nrow=2,byrow=T)
    #Initialize aOut to hold a estimate at each time point
    aOut <- a
    
    #Define matrix C, parameters of covariates [2x5]
    C=matrix(c(k,1,1,-iota,0,0,0,0,0,0,1,1),nrow=2,byrow=TRUE)
    
    #Iterate through time
    for(t in 2:nObs){
     
      #Define transition matrix B [2x2]
      B <- matrix(c(1-qOut[t]/v[t]-k/zMix[t],r,0,1-r-qOut[t]/v[t]),nrow=2,byrow=TRUE)
      
      #Define covariates matrix ut [6x1]
      ut <- matrix(c(DICeq[t]*v[t]/zMix[t],qIn[t]*DICIn[t],precip[t]*DICPrecip[t],0.2*light[t]*v[t],qIn[t]*DOCIn[t],precip[t]*DOCPrecip[t]),ncol=1)
      
      #New observation
      y=matrix(as.numeric(Y[t,]),nrow=2,byrow=TRUE)
      
      #Predictions (a is estimate of state, P is covariance matrix of estimation error)
      a = (B %*% a) + (C%*%ut)
      P = (B %*% P %*% t(B))+Q
      
      #Update
      if (all(is.na(y))==FALSE) {  #If no obs at this time step, skip updating and don't accumulate to likelihood. Otherwise...
        
        #Identify state variables for which obs are available at this time step, adjust y etc for missing obs
        #See Harvey Eq. 3.4.75 f.f. for this approach (weights matrix W etc)
        avail <- which(is.na(y)==FALSE)
        nT <- length(avail)
        W <- matrix(0,nrow=nT,ncol=nVar,byrow=TRUE) 
        W[cbind((1:length(avail)),avail)] <- 1
        y[which(is.na(y))] <- 0 #To make matrix math in next line work
        yT <- W %*% y
        
        #Update
        V = yT-a[avail]  # residuals
        F = P+H
        Finv_P = solve(F[avail,avail],P[avail,avail])
        Finv_V = solve(F[avail,avail],V)
        a[avail,] = a[avail,] + (P[avail,avail] %*% Finv_V)
        P[avail,avail] = P[avail,avail] - (P[avail,avail] %*% Finv_P)
      }
      
      #Save the estimates of a and from this time step
      aOut <- cbind(aOut,a)
      POut <- abind(POut,P,along=3)
      
    } # End iteration
    
    #Return the state estimates and their covariance matrices
    return(list(a=aOut,P=POut))
    
  }) #end 'with'
} #end function
