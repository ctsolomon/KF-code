#Function to run the Kalman filter on simple C model
#CTS 25 July 2014, based on earlier version by SEJ, all built from SRC code

#Arguments are
#  p - Parameters to be estimated
#  Y - Observations of state variables. Rows are time points, columns are (DIC, DOC)
#  H - Estimate of the observation error covariance matrix
#  otherPar - Additional fixed parameters, not to be estimated, needed for calculations
#  drivers - A data.frame with the driver data needed for the ut matrix.

kalmanStep<-function(p,Y,H,drivers){
  with(drivers,{
    #Unpack estimated parameters
    r <- p[1]
    iota <- p[2]
    k <- p[3]
    Q <- matrix(c(exp(p[4]),0,0,exp(p[5])),nrow=2,byrow=TRUE)	#process error covariance matrix [2x2]
    
    #Initialize P from Q
    P <- Q
  
    #Set up terms of NLL expression
    VFV <- 0
    logdetF <- 0
    
    #Number of observations
    nObs <- nrow(Y)
    
    #Number of state variables
    nVar <- 2
    
    #Initialize a (estimate of state)
    a <- matrix(as.numeric(Y[1,]),nrow=nVar,byrow=T)
       
    #Define matrix C, parameters of covariates [2x5]
    C=matrix(c(k,1,1,-iota,0,0,0,0,0,0,1,1),nrow=2,byrow=TRUE)
    
    #Iterate through time
    for(t in 2:nObs){
      
      #Define transition matrix B [2x2]
      B <- matrix(c(1-qOut[t]/v[t]-k/zMix[t],r,0,1-r-qOut[t]/v[t]),nrow=nVar,byrow=TRUE)
      
      #Define covariates matrix ut [6x1]
      ut <- matrix(c(DICeq[t]*v[t]/zMix[t],qIn[t]*DICIn[t],precip[t]*DICPrecip[t],0.2*light[t]*v[t],qIn[t]*DOCIn[t],precip[t]*DOCPrecip[t]),ncol=1)
      
      #New observation
      y=matrix(as.numeric(Y[t,]),nrow=nVar,byrow=TRUE)
             
      #Predictions
      a = (B %*% a) + (C%*%ut)
      P = (B %*% P %*% t(B))+Q
      
      #Update
      if (all(is.na(y))==FALSE) {  #If no obs at this time step, skip updating and don't accumulate to likelihood. Otherwise...
        
        #Identify state variables for which obs are available at this time step, adjust y etc for missing obs
        #See Harvey Eq. 3.4.75 f.f. for this approach (weights matrix W etc)
        avail <- which(is.na(y)==FALSE)
        nT <- length(avail)
        W <- matrix(0,nrow=nT,ncol=nVar,byrow=TRUE) 
        W[cbind((1:length(avail)),avail)] <- 1  ###>>>>double-check that this works as intended
        y[which(is.na(y))] <- 0 #To make matrix math in next line work
        yT <- W %*% y
        
        #Update
        V = yT-a[avail]  # residuals
        F = P+H #P[avail,avail]+H[avail,avail]  ###>>>>>>>>>Need to deal with dimensions here
        Finv_P = solve(F,P)
        Finv_V = solve(F,V)
        a = a + (P %*% Finv_V)
        P = P - (P %*% Finv_P)
        
        #Accumulate likelihood terms
        logdetF = c( logdetF, 0.5*log(det(F)) )
        VFV = c( VFV, 0.5*t(V)%*%Finv_V )
      }

    } # End iteration
    
    # Calculate negative log likelihood
    NLL = (0.5*nObs*log(2*pi)) + sum(logdetF,na.rm=TRUE) + sum(VFV,na.rm=TRUE)
    
    return(NLL)
  }) #end 'with'
} #end function
