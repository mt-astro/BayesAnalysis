# Utility programs for use with the book,
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   This script has been lightly edited by Neil Frazer <---------------- NB

# A list of Kruschke's functions w/ approximate starting line numbers
# openGraph(), 67. 
# saveGraph(), 87.
# HDIofMCMC(), 96. Indep.
# HDIofICDF(), 136. Indep.
# HDIofGrid(), 163. Indep.
# DbdaAcfPlot(), 175. coda, 
# DbdaDensPlot(), 199. coda, HDIofMCMC,
# diagMCMC(), 230. coda, DbdaAcfPlot, DbdaDensPlot
# diagStanFit(), 264. coda, rstan, dbdaAcfPlot, dbdaDensPlot
# normalize(), 305. Indep.
# summarizePost(), 309. coda, HDIofMCMC,
# plotPost(), 375. coda, HDIofMCMC, 
# betaABfromMeanKappa(), 564.
# betaABfromModeKappa(), 572.
# gammaShRaFromMeanSD(), 590.
# gammaShRaFromModeSD(), 598.
# genYwithOut(), 618.

# Check that required packages are installed:
want  <-  c("parallel","rjags","runjags","compute.es")
have  <-  want %in% rownames(installed.packages())
if ( any(!have) ) 
  install.packages(want[!have], repos="http://cran.us.r-project.org" )

# Load rjags. Assumes JAGS is already installed.
try( library(rjags) )
# Load runjags. Assumes JAGS is already installed.
try( library(runjags) )
try( runjags.options( inits.warning=FALSE , rng.warning=FALSE ) )

# set default number of chains and parallelness for MCMC:
library(parallel) # for detectCores().
nCores <- detectCores() 
if ( !is.finite(nCores) ) { nCores <- 1 } 
if ( nCores > 4 ) { 
  nChainsDefault <- 4  # because JAGS has only 4 rng's.
  runjagsMethodDefault = "parallel"
}
if ( nCores == 4 ) { 
  nChainsDefault = 3  # save 1 core for other processes.
  runjagsMethodDefault = "parallel"
}
if ( nCores < 4 ) { 
  nChainsDefault = 3 
  runjagsMethodDefault = "rjags" # NOT parallel
}

# Functions for opening and saving graphics that operate the same for 
# Windows and Macintosh and Linux operating systems. At least, that's the hope!
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openGraph <- function( width=7 , height=7 , mag=1.0 , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    tryInfo <- 
      try( X11( width=width*mag, height=height*mag, type="cairo", ... ) )
    if ( class(tryInfo) == "try-error" ) {
      lineInput <- readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      X11( width=width*mag , height=height*mag , type="cairo" , ... )
    }
  } else { # Windows OS
    tryInfo <- 
      try( windows( width=width*mag , height=height*mag , ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput <- readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      windows( width=width*mag , height=height*mag , ... )    
    }
  }
} # end openGraph()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveGraph <- function( file="saveGraphOutput" , type="pdf" , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    if ( any( type == c("png","jpeg","jpg","tiff","bmp")) ) {
      sptype <- type    # for savePlot()
      if ( type == "jpg" ) { sptype <- "jpeg" } 
      savePlot( file=paste0(file,".",type) , type=sptype , ... )     
    }
    if ( type == "pdf" ) {
      dev.copy2pdf(file=paste0(file,".",type) , ... )
    }
    if ( type == "eps" ) {
      dev.copy2eps(file=paste0(file,".",type) , ... )
    }
  } else { # Windows OS
    file <- paste0(file,".",type) 
    savePlot( file=file , type=type , ... )
  }
} # end saveGraph()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions for computing limits of HDI's:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HDIofMCMC <- function( sampleVec , credMass=0.95 ) {
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
  sortedPts <- sort( sampleVec )
  ciIdxInc <- # How many sample elements are in the CI?
    ceiling( credMass * length( sortedPts ) )
  nCIs <-     # number of possible CIs
    length( sortedPts ) - ciIdxInc
  ciWidth <- rep( 0 , nCIs ) # create storage
  for ( i in 1:nCIs ) {
    ciWidth[ i ] <- # last value in CI - first value in CI
      sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  } # end for
  HDImin <- sortedPts[ which.min( ciWidth ) ]
  HDImax <- sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim <- c( HDImin , HDImax )
  return( HDIlim )
} # end HDIofMCMC()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HDIofICDF <- function( ICDFname , credMass=0.95 , tol=1e-8 , ... ) {
  # Arguments:
  #   ICDFname is R's name for the inverse cumulative density function
  #     of the distribution.
  #   credMass is the desired mass of the HDI region.
  #   tol is passed to R's optimize function.
  # Return value:
  #   Highest density iterval (HDI) limits in a vector.
  # Example of use: For determining HDI of a beta(30,12) distribution, type
  #   HDIofICDF( qbeta , shape1 = 30 , shape2 = 12 )
  #   Notice that the parameters of the ICDFname must be explicitly named;
  #   e.g., HDIofICDF( qbeta , 30 , 12 ) does not work.
  # Adapted and corrected from Greg Snow's TeachingDemos package.
  incredMass <- 1.0 - credMass
  intervalWidth <- function( lowTailPr , ICDFname , credMass , ... ) {
    ICDFname( credMass + lowTailPr , ... ) - ICDFname( lowTailPr , ... )
  }
  optInfo <- # list(minimum = <argmin>, objective = <minimum>)
    optimize( intervalWidth , # lowTailPr that minimizes intervalWidth
              c( 0 , incredMass ) , 
              ICDFname=ICDFname ,
              credMass=credMass , tol=tol , ... )
  HDIlowTailPr <- optInfo$minimum
  return( c( ICDFname( HDIlowTailPr , ... ) ,
             ICDFname( credMass + HDIlowTailPr , ... ) ) )
} # end HDIofICDF()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HDIofGrid <- function( probMassVec , credMass=0.95 ) {
  # Arguments:
  #   probMassVec is a vector of probability masses at each grid point.
  #   credMass is the desired mass of the HDI region.
  # Return value:
  #   A list with components:
  #   indices is a vector of indices that are in the HDI
  #   mass is the total mass of the included indices
  #   height is the smallest component probability mass in the HDI
  # Example of use: For determining HDI of a beta(30,12) distribution
  #   approximated on a grid:
  #   > probDensityVec = dbeta( seq(0,1,length=201) , 30 , 12 )
  #   > probMassVec = probDensityVec / sum( probDensityVec )
  #   > HDIinfo = HDIofGrid( probMassVec )
  #   > show( HDIinfo )
  sortedProbMass <- sort( probMassVec , decreasing=TRUE )
  HDIheightIdx   <- min( which( cumsum( sortedProbMass ) >= credMass ) )
  HDIheight      <- sortedProbMass[ HDIheightIdx ]
  HDImass        <- sum( probMassVec[ probMassVec >= HDIheight ] )
  return( list( indices = which( probMassVec >= HDIheight ) ,
                mass = HDImass , height = HDIheight ) )
} # end HDIofGrid()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function(s) for plotting properties of mcmc coda objects.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DbdaAcfPlot <- function( codaObject , 
                         parName=varnames(codaObject)[1] , 
                         plColors=NULL ) {
  # Calls coda::effectiveSize()
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject) # or coda::nchains()
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
} # end DbdaAcfPlot()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DbdaDensPlot <- function( codaObject , 
                          parName=varnames(codaObject)[1] , 
                          plColors=NULL 
                          ) {
  # Calls coda::effectiveSize(), HDIofMCMC()
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain <- length(codaObject) # or coda::nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat <- NULL
  yMat <- NULL
  hdiLims <- NULL
  for ( cIdx in 1:nChain ) {
    densInfo <- density( codaObject[,c(parName)][[cIdx]] ) 
    xMat <- cbind(xMat,densInfo$x)
    yMat <- cbind(yMat,densInfo$y)
    hdiLims <- cbind( hdiLims, HDIofMCMC(codaObject[,c(parName)][[cIdx]]) )
  }
  matplot( xMat , yMat , type="l" , col=plColors , 
           main="" , xlab="Param. Value" , ylab="Density" )
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims) , 0 , "95% HDI" , adj=c(0.5,-0.2) )
  EffChnLngth <- effectiveSize(codaObject[,c(parName)])
  MCSE <- sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth) 
  text( max(xMat) , max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
        paste("MCSE =\n",signif(MCSE,3)) )
} # end DbdaDensPlot()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diagMCMC <- 
  function( codaObject , parName=varnames(codaObject)[1] ) {
  # Calls coda::traceplot()&gelman.plot(), DbdaAcfPlot(), DbdaDensPlot()
  
  DBDAplColors <- c("skyblue","black","royalblue","steelblue")
  par( mar=0.5+c(3,4,1,0) , 
       oma=0.1+c(0,0,2,0) , 
       mgp=c(2.25,0.7,0) , 
       cex.lab=1.5 
       )
  
  layout(matrix(1:4,nrow=2))
  
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors ) 
  tryVal <-  try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
                       col=DBDAplColors ) )  
  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() # skip this plot
    print(paste0("Warning: coda::gelman.plot fails for ", parName))
  } else if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
    plot.new() # skip this plot
    print(paste0("Warning: coda::gelman.plot fails for ", parName))
  }

  DbdaAcfPlot(codaObject, parName, plColors=DBDAplColors)
  DbdaDensPlot(codaObject, parName, plColors=DBDAplColors)
  mtext( text=parName, outer=TRUE, adj=c(0.5,0.5), cex=2.0 )

} # end diagMCMC()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diagStanFit <- function( stanFit , parName ,
                        saveName=NULL , saveType="jpg" ) {
  # Calls coda::mcmc.list&gelman.plot, rstan::traceplot, dbdaAcfPlot, dbdaDensPlot
  codaFit <- 
    coda::mcmc.list( 
      lapply( 1:ncol(stanFit), function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
  
  DBDAplColors <- c("skyblue","black","royalblue","steelblue")
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))

  require(rstan)
  rstan::traceplot(stanFit,pars=parName,nrow=1,ncol=1)#,main="",ylab="Param. Value",col=DBDAplColors) 

  require(coda)
  tryVal <- try(
    coda::gelman.plot( codaObject[,c(parName)] , 
                       main="" , 
                       auto.layout=FALSE , 
                       col=DBDAplColors ) )
  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  } else if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  }

  DbdaAcfPlot(codaFit,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaFit,parName,plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"Diag",parName), type=saveType)
  }
} # end diagStanFit()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions for summarizing and plotting distribution of a large sample; 
# typically applied to MCMC posterior.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalize <- function( v ){ return( v / sum(v) ) }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
require(coda) # loaded by rjags, but redundancy doesn't hurt
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summarizePost <- function( # summarize a large sample
  paramSampleVec, compVal=NULL, ROPE=NULL, credMass=0.95 ) {
  ## Author: J.K.Kruschke, in DBDA2E-utilities.R
  # Needs coda, HDIofMCMC
  meanParam   <- mean(    paramSampleVec )
  medianParam <- median(  paramSampleVec )
  
  dres        <- density( paramSampleVec )
  ## stats::density() returns an object of class "density"
  ## which is a list with components x (the n coordinates of
  ## the points where the density was estimated), y (the n
  ## estimated density values, all ≥ 0), bw (the bandwidth
  ## used), n (the sample size after elimination of missing
  ## values), call (the call that produced the result), 
  ## data.name (the name of the x argument), and has.na
  ## (a logical argument included for compatibility; it is
  ## always FALSE).
  modeParam   <- dres$x[which.max(dres$y)]
  
  mcmcEffSz <- round( effectiveSize( paramSampleVec ) , 1 )
  ## coda::effectiveSize() takes an argument x of class mcmc,
  ## or mcmc.list, and returns a vector giving the effective
  ## sample size for each column of x. If x is of class
  ## mcmc.list, the effective sizes are summed across chains.
  ## To get the size for each chain individually use
  ## lapply(x,effectiveSize)
  names(mcmcEffSz) <- NULL # remove any names
  
  ## Get HDI limits
  hdiLim <- HDIofMCMC( paramSampleVec , credMass=credMass )
  
  ## Get percent > compVal
  if ( !is.null(compVal) ) {
    pcgtCompVal <- ( 100 * sum( paramSampleVec > compVal ) 
                    / length( paramSampleVec ) )
  } else {
    compVal <- NA
    pcgtCompVal <- NA
  } # end if
  
  ## Get pcltRope, pcgtRope, pcinRope
  if ( !is.null(ROPE) ) {
    pcltRope <- ( 100 * sum( paramSampleVec < ROPE[1] ) 
                 / length( paramSampleVec ) )
    pcgtRope <- ( 100 * sum( paramSampleVec > ROPE[2] ) 
                 / length( paramSampleVec ) )
    pcinRope <- 100-(pcltRope+pcgtRope)
  } else { 
    ROPE <- c(NA,NA)
    pcltRope <- NA 
    pcgtRope <- NA 
    pcinRope <- NA 
  } # end if
  
  ## Return a vector with named components
  return( c( 
      Mean=meanParam      , Median=medianParam , Mode=modeParam , 
      ESS=mcmcEffSz       ,
      HDImass=credMass    , HDIlow=hdiLim[1] , HDIhigh=hdiLim[2] , 
      CompVal=compVal     , PcntGtCompVal=pcgtCompVal , 
      ROPElow=ROPE[1]     , ROPEhigh=ROPE[2] ,
      PcntLtROPE=pcltRope , PcntInROPE=pcinRope , PcntGtROPE=pcgtRope ) 
      )

} # end summarizePost()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotPost <- function( # from J.K.Kruschke's DBDA2E-utilities.R
  paramSampleVec , 
  cenTend = c("mode","median","mean")[2] , 
  compVal = NULL, 
  ROPE = NULL, 
  credMass = 0.95, 
  HDItextPlace = 0.7, 
  xlab = "Param. Val." , 
  xlim = range( c( compVal , ROPE , paramSampleVec ) ) , 
  yaxt = "n" , 
  ylab = "" , 
  main = "" , 
  cex = 1.4 , 
  cex.lab = 1.5 ,
  col = "skyblue" , 
  border = "white" , 
  showCurve = FALSE , 
  breaks = NULL , 
  ... ) {
  
  # requires coda package and JKK's HDIofMCMC() <--------------------- NB
  
  # convert coda object to matrix:
  if ( class(paramSampleVec) == "mcmc.list" ) {
    paramSampleVec <- as.matrix(paramSampleVec)
  } # Class mcmc.list is described in the coda manual, page 23.
    # An object of class mcmc is a vector or else a matrix with with one
    #  column for each variable. An object of class mcmc.list is a list
    #  containing objects of class mcmc all with the same length.
    #  An object of class mcmc has a dimnames attribute that is a list of
    #  two character vectors. The first character vector consists of the
    #  row numbers in character form; the second character vector contains
    #  the names of the parameters. If xx is an object of class mcmc then
    #  as.matrix(xx) discards the character vector of rownames and keeps
    #  the character vector of colnames. If yy is an object of class 
    #  mcmc.list then as.matrix(yy) 'stacks' the matrices associated with
    #  the mcmc objects in the list to make a single matrix with the same
    #  number of columns, and the same column names as one would get from 
    #  the command as.matrix(yy[1]), or the command as.matrix(yy[[1]]) which
    #  gives the same thing. Naturally the stacking process results in a 
    #  matrix with many more rows: if there are nc chains in yy with nr rows
    #  in each chain then the number of rows in as.matrix(yy) will be 
    #  (nc)(nr).
  
  summaryColNames <- c("ESS","mean","median","mode",
                      "hdiMass","hdiLow","hdiHigh",
                      "compVal","pGtCompVal",
                      "ROPElow","ROPEhigh","pLtROPE","pInROPE","pGtROPE")
  
  postSummary <- # create an empty one-row matrix with column names
    matrix( NA , nrow=1 , ncol=length(summaryColNames) , 
            dimnames=list( c( xlab ) , summaryColNames ) 
            )
  # Now put numbers into that single row of the matrix
  postSummary[,   "ESS"] <- effectiveSize(paramSampleVec)
  ## coda::effectiveSize() takes an argument x of class mcmc,
  ## or mcmc.list, and returns a vector giving the effective
  ## sample size for each column of x. If x is of class
  ## mcmc.list, the effective sizes are summed across chains.
  ## To get the size for each chain individually use
  ## lapply(x,effectiveSize)
  
  postSummary[,  "mean"] <- mean(   paramSampleVec)
  
  postSummary[,"median"] <- median( paramSampleVec)
  
  mcmcDensity            <- density(paramSampleVec)
  ## stats::density() returns an object of class "density"
  ## which is a list with components: x (the n coordinates of
  ## the points where the density was estimated), y (the n
  ## estimated density values, all ≥ 0), bw (the bandwidth
  ## used), n (the sample size after elimination of missing
  ## values), call (the call that produced the result), 
  ## data.name (the name of the x argument), and has.na
  ## (a logical argument included for compatibility; it is
  ## always FALSE).
  
  postSummary[,  "mode"] <- mcmcDensity$x[which.max(mcmcDensity$y)]
  
  postSummary[,"hdiMass"] <- credMass
  
  HDI <- HDIofMCMC( paramSampleVec , credMass ) # <--------------------- NB

  postSummary[,"hdiLow" ] <- HDI[1]
  
  postSummary[,"hdiHigh"] <- HDI[2]
  
  # Plot histogram.
  cvCol   <- "darkgreen"
  ropeCol <- "darkred"
  
  if ( is.null(breaks) ) { # make breaks for histogram
    if ( max(paramSampleVec) > min(paramSampleVec) ) {
      breaks <- c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                       by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec) )
    } else {
      breaks <- c(min(paramSampleVec)-1.0E-6,max(paramSampleVec)+1.0E-6)
      border <- "skyblue"
    }
  }
  
  if ( !showCurve ) { # plot a histogram
    par(xpd=NA)
    histinfo <- hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks , ... )
  } else { # plot a density curve
    par(xpd=NA)
    histinfo  <- hist( paramSampleVec , plot=F )
    densCurve <- density( paramSampleVec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  
  cenTendHt  <- 0.90*max(histinfo$density)
  cvHt       <- 0.70*max(histinfo$density)
  ROPEtextHt <- 0.55*max(histinfo$density)
  
  # Display central tendency:
  mn  <- mean(paramSampleVec)
  med <- median(paramSampleVec)
  mcmcDensity <- density(paramSampleVec)
  mo <- mcmcDensity$x[which.max(mcmcDensity$y)]
  
  if ( cenTend=="mode" ){ 
    text( mo , cenTendHt ,
          bquote(mode==.(signif(mo,3))) , adj=c(.5,0) , cex=cex )
  } else if ( cenTend=="median" ){ 
    text( med , cenTendHt ,
          bquote(median==.(signif(med,3))) , adj=c(.5,0) , cex=cex , col=cvCol )
  } else if ( cenTend=="mean" ){ 
    text( mn , cenTendHt ,
          bquote(mean==.(signif(mn,3))) , adj=c(.5,0) , cex=cex )
  }
  
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    pGtCompVal <- sum( paramSampleVec > compVal ) / length( paramSampleVec ) 
    pLtCompVal <- 1 - pGtCompVal
    lines( c(compVal,compVal) , c(0.96*cvHt,0) , 
           lty="dotted" , lwd=2 , col=cvCol )
    text( compVal , cvHt ,
          bquote( .(round(100*pLtCompVal,1)) * "% < " *
                   .(signif(compVal,3)) * " < " * 
                   .(round(100*pGtCompVal,1)) * "%" ) ,
          adj=c(pLtCompVal,0) , cex=0.8*cex , col=cvCol )
    postSummary[,   "compVal"] <- compVal
    postSummary[,"pGtCompVal"] <- pGtCompVal
  }
  
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    pInROPE <- ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                / length( paramSampleVec ) )
    pGtROPE <- ( sum( paramSampleVec >= ROPE[2] ) / length( paramSampleVec ) )
    pLtROPE <- ( sum( paramSampleVec <= ROPE[1] ) / length( paramSampleVec ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol)
    text( mean(ROPE) , ROPEtextHt ,
          bquote( .(round(100*pLtROPE,1)) * "% < " * .(ROPE[1]) * " < " * 
                  .(round(100*pInROPE,1)) * "% < " * .(ROPE[2]) * " < " * 
                  .(round(100*pGtROPE,1)) * "%" 
                  ) ,
          adj=c(pLtROPE+.5*pInROPE,0) , cex=1 , col=ropeCol 
          )
    postSummary[, "ROPElow"] <- ROPE[1] 
    postSummary[,"ROPEhigh"] <- ROPE[2] 
    postSummary[, "pLtROPE"] <- pLtROPE
    postSummary[, "pInROPE"] <- pInROPE
    postSummary[, "pGtROPE"] <- pGtROPE
  }
  
  # Display the HDI.
  lines( HDI , c(0,0) , lwd=4 , lend=1 )
  text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
        adj=c(.5,-1.7) , cex=cex )
  text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
        adj=c(HDItextPlace,-0.5) , cex=cex )
  text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
        adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  par(xpd=FALSE)
  
  # Done
  return( postSummary )
} # end plotPost()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Shape parameters from central tendency and scale:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaABfromMeanKappa <- function( mean , kappa ) {
  if ( mean <=0 | mean >= 1) stop("must have 0 < mean < 1")
  if ( kappa <=0 ) stop("kappa must be > 0")
  a <- mean * kappa
  b <- ( 1.0 - mean ) * kappa
  return( list( a=a , b=b ) )
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaABfromModeKappa <- function( mode , kappa ) {
  if ( mode <=0 | mode >= 1) stop("must have 0 < mode < 1")
  if ( kappa <=2 ) stop("kappa must be > 2 for mode parameterization")
  a <- mode * ( kappa - 2 ) + 1
  b <- ( 1.0 - mode ) * ( kappa - 2 ) + 1
  return( list( a=a , b=b ) )
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaABfromMeanSD <- function( mean , sd ) {
  if ( mean <= 0 | mean >= 1) stop("must have 0 < mean < 1")
  if ( sd <= 0              ) stop("sd must be > 0")
  kappa <- mean*(1-mean)/sd^2 - 1
  if ( kappa <= 0 ) stop("invalid combination of mean and sd")
  a <- mean * kappa
  b <- ( 1.0 - mean ) * kappa
  return( list( a=a , b=b ) )
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaShRaFromMeanSD <- function( mean , sd ) {
  if ( mean <= 0 ) stop("mean must be > 0")
  if (  sd  <= 0 ) stop(" sd  must be > 0")
  shape <- mean^2/sd^2
  rate  <- mean/sd^2
  return( list( shape=shape , rate=rate ) )
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaShRaFromModeSD <- function( mode , sd ) {
  if ( mode <= 0 ) stop("mode must be > 0")
  if ( sd   <= 0 ) stop(" sd  must be > 0")
  rate   <- ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape  <- 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Make some data files for examples...
createDataFiles <- FALSE # <--------------------------- NB
if ( createDataFiles ) {
  source("HtWtDataGenerator.R")
  N <- 300
  m <- HtWtDataGenerator( N , rndsd=47405 )
  write.csv( file=paste0("HtWtData",N,".csv") , row.names=FALSE , m )
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Function for generating normal data with normal outliers:
  genYwithOut <- function( N , pcntOut=15 , sdOut=3.0 ) {
    inl <- rnorm( N-ceiling(pcntOut/100*N) )
    out <- rnorm(   ceiling(pcntOut/100*N) )
    inl <- (inl-mean(inl))/sd(inl)
    out <- (out-mean(out))/sd(out) * sdOut
    return( c(inl,out) )
  }
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Two-group IQ scores with outliers 
  set.seed(47405)
  y1 <- round(pmax(50,genYwithOut(63,20,3.5)*17.5+106))
  y2 <- round(pmax(50,genYwithOut(57,20,3.5)*10+100))
  write.csv( file="TwoGroupIQ.csv" , row.names=FALSE ,
             data.frame( Score=c(y1,y2) , 
                         Group=c(rep("Smart Drug",length(y1)),
                                 rep("Placebo",length(y2))) ) )
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # One-group log-normal
  set.seed(47405)
  z <- rnorm(123)
  logY <- (z-mean(z))/sd(z) * 0.5 + 5.5 # logY has mean 5.5 and sd 0.5
  y <- round( exp(logY) , 2 )
  write.csv( file="OneGroupLogNormal.csv" , row.names=FALSE ,
             cbind(y) )
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # One-group gamma
  desiredMode <- 250
  desiredSD <- 100
  desiredRate <- (desiredMode+sqrt(desiredMode^2+4*desiredSD^2))/(2*desiredSD^2)
  desiredShape <- 1+desiredMode*desiredRate
  set.seed(47405)
  y <- round( rgamma( 153 , shape=desiredShape , rate=desiredRate ) , 2 )
  write.csv( file="OneGroupGamma.csv" , row.names=FALSE , cbind(y) )
} # end if (createDataFiles)
