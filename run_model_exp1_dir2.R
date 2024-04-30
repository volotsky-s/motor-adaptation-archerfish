
#===============================================================================

runMCMC = function( datFrm , yName="y" , x1Name="x1" , x2Name="x2" ,
                    numSavedSteps=50000 , thinSteps=2 , saveName=NULL ,
                    runjagsMethod="parallel" , 
                    nChains=3 ) { 
  
  #------------------------------------------------------------------------------
  
  # Check that required packages are installed:
  want = c("parallel","rjags","runjags","compute.es")
  have = want %in% rownames(installed.packages())
  if ( any(!have) ) { install.packages( want[!have] ) }
  
  # Load rjags. Assumes JAGS is already installed.
  try( library(rjags) )
  # Load runjags. Assumes JAGS is already installed.
  try( library(runjags) )
  try( runjags.options( inits.warning=FALSE , rng.warning=FALSE ) )
  
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  y = as.numeric(datFrm[,yName])
  x1 = as.numeric(as.factor(datFrm[,x1Name]))
  x1levels = levels(as.factor(datFrm[,x1Name]))
  x2 = as.numeric(as.factor(datFrm[,x2Name]))
  x2levels = levels(as.factor(datFrm[,x2Name]))
  Ntotal = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  
  
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  yMean = mean(y , na.rm=TRUE)
  ySD = sd(y , na.rm=TRUE)
  
  
  yEpochMean = c(mean(y[x1==1]), mean(y[x1==2]), mean(y[x1==3]), mean(y[x1==4]), 
                 mean(y[x1==5]), mean(y[x1==6]))
  
  
  ySubjectMean = c(mean(y[x2==1]), mean(y[x2==2]),  mean(y[x2==3]), mean(y[x2==4]))
  
  
  yEpochSubjectMean = matrix(data=NA,nrow=Nx1Lvl,ncol=Nx2Lvl)
  for (e in 1:Nx1Lvl) {for (s in 1:Nx2Lvl) { yEpochSubjectMean[e,s] = mean(y[x1==e & x2==s]) } }
  
  
  # For hyper-prior on deflections:
  gammaRate =  ( ySD/2 + sqrt( (ySD/2)^2 + 4 * (2*ySD)^2 ) ) / ( 2 * (2*ySD)^2 )
  gammaShape = 1 + ySD/2 * gammaRate  
  
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    x1 = x1 ,
    x2 = x2 ,
    
    Ntotal = Ntotal ,
    Nx1Lvl = Nx1Lvl ,
    Nx2Lvl = Nx2Lvl ,
    
    # data properties for scaling the prior:
    yMean = yMean ,
    ySD = ySD ,
    yEpochMean = yEpochMean,
    ySubjectMean = ySubjectMean,
    yEpochSubjectMean = yEpochSubjectMean,
    gammaRate = gammaRate,
    gammaShape = gammaShape
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dnorm( mu[i] , 1/(ySigma^2) )
      mu[i] <- a0 + aE[x1[i]] + aS[x2[i]] + aES[x1[i],x2[i]]
    }

    a0 ~ dnorm( yMean , 1/(ySD*5)^2 ) 
    ySigma ~ dgamma(gammaShape,gammaRate) 
    #
    for ( j1 in 1:Nx1Lvl ) { aE[j1] ~ dnorm( yEpochMean[j1] , 1/aE_SD^2 ) }
    aE_SD ~ dgamma(gammaShape,gammaRate) 
    #
    for ( j2 in 1:Nx2Lvl ) { aS[j2] ~ dnorm( ySubjectMean[j2] , 1/aS_SD^2 ) }
    aS_SD ~ dgamma(gammaShape,gammaRate) 
    #
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) { 
      aES[j1,j2] ~ dnorm( yEpochSubjectMean[j1,j2] , 1/aES_SD^2 )
    } } 
    aES_SD ~ dgamma(gammaShape,gammaRate) 
    
    # Convert a0,a1[],a2[],a1a2[,] to sum-to-zero b0,b1[],b2[],b1b2[,] :
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      m[j1,j2] <- a0 + aE[j1] + aS[j2] + aES[j1,j2] # cell means 
    } }
    b0 <- mean( m[1:Nx1Lvl,1:Nx2Lvl] )
    for ( j1 in 1:Nx1Lvl ) { bE[j1] <- mean( m[j1,1:Nx2Lvl] ) - b0 }
    for ( j2 in 1:Nx2Lvl ) { bS[j2] <- mean( m[1:Nx1Lvl,j2] ) - b0 }
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      bES[j1,j2] <- m[j1,j2] - ( b0 + bE[j1] + bS[j2] )  
    } }
    
  }
  " # close quote for modelstring
  writeLines(modelstring,con="TEMPmodel.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS it automatically...
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  require(rjags)
  require(runjags)
  parameters = c( "b0" ,  "bE" ,  "bS" ,"bES" , "m" , "ySigma" , 
                  "ySD" , "aE_SD" , "aS_SD", "aES_SD" )
  adaptSteps = 2000 
  burnInSteps = 2000 
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
}

#===============================================================================
#===============================================================================

openGraph = function( width=7 , height=7 , mag=1.0 , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    tryInfo = try( X11( width=width*mag , height=height*mag , type="cairo" , 
                        ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      X11( width=width*mag , height=height*mag , type="cairo" , ... )
    }
  } else { # Windows OS
    tryInfo = try( windows( width=width*mag , height=height*mag , ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      windows( width=width*mag , height=height*mag , ... )    
    }
  }
}

#------------------------------------------------------------------------------
# Function(s) for plotting properties of mcmc coda objects.

DbdaAcfPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  for ( cIdx in 1:nChain ) {
    acfInfo = acf(codaObject[,c(parName)][[cIdx]],plot=FALSE) 
    xMat = cbind(xMat,acfInfo$lag)
    yMat = cbind(yMat,acfInfo$acf)
  }
  matplot( xMat , yMat , type="o" , pch=20 , col=plColors , ylim=c(0,1) ,
           main="" , xlab="Lag" , ylab="Autocorrelation" )
  abline(h=0,lty="dashed")
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  text( x=max(xMat) , y=max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
        labels=paste("ESS =",round(EffChnLngth,1)) )
}

DbdaDensPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject) # or nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  for ( cIdx in 1:nChain ) {
    densInfo = density(codaObject[,c(parName)][[cIdx]]) 
    xMat = cbind(xMat,densInfo$x)
    yMat = cbind(yMat,densInfo$y)
    hdiLims = cbind(hdiLims,HDIofMCMC(codaObject[,c(parName)][[cIdx]]))
  }
  matplot( xMat , yMat , type="l" , col=plColors , 
           main="" , xlab="Param. Value" , ylab="Density" )
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims) , 0 , "95% HDI" , adj=c(0.5,-0.2) )
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  MCSE = sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth) 
  text( max(xMat) , max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
        paste("MCSE =\n",signif(MCSE,3)) )
}

diagMCMC = function( codaObject , parName=varnames(codaObject)[1] ,
                     saveName=NULL , saveType="jpg" ) {
  DBDAplColors = c("skyblue","black","royalblue","steelblue")
  openGraph(height=5,width=7)
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot and gelman.plot are from CODA package:
  require(coda)
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors ) 
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
                       col=DBDAplColors )
  )  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning: coda::gelman.plot fails for ",parName))
    }
  }
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  
}

#------------------------------------------------------------------------------
# Functions for computing limits of HDI's:

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

