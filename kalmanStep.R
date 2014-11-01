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
    Q <- matrix(c(exp(p[4]),0,0,exp(p[5])),nrow=2,byrow=TRUE)  #process error covariance matrix [2x2]
    
    #Initialize P from Q
    P <- Q
    
    #Set up terms of NLL expression
    VFV <- 0
    logdetF <- 0
    
    #Number of observations
    nObs <- nrow(Y)
    
    #Initialize a (estimate of state)
    a <- matrix(as.numeric(Y[1,]),nrow=2,byrow=T)
    
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
      
      #Predictions
      a = (B %*% a) + (C%*%ut)
      P = (B %*% P %*% t(B))+Q
      
      #Update (as long as there are no NA in the y vector at this time point)
      if (any(is.na(y))==FALSE) {
        V = y-a  # residuals
        F = P+H
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
