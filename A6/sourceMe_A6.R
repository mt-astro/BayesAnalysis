

HDIofICDF <- function( ICDFname , credMass=0.95 , tol=1e-8 , ... ) {
  # Arguments:
  #   ICDFname will be an inverse cumulative density function, e.g., qbeta.
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
#######################################################################
#######################################################################



BernBeta <- function( priorBetaAB, Data, plotType=c("Points","Bars")[2],
                      showCentTend=c("Mean","Mode","None")[3],
                      showHDI=c(TRUE,FALSE)[2], HDImass=0.95,
                      showpD =c(TRUE,FALSE)[2], ROPE=NULL) {
  # priorBetaAB is two-element vector of beta a,b shape parameters
  # Data is vector of 0's and 1's.
  # Check for input errors:
  if ( any( priorBetaAB < 0 ) ) {
    stop("priorBetaAB values must be nonnegative")
  }
  if ( length( priorBetaAB ) != 2 ) {
    stop("priorBetaAB must be a vector of two values")
  }
  if ( !all( Data == 1 | Data == 0 ) ) {
    stop("Data values must be 0 or 1")
  }
  
  # For notational convenience, rename beta shape parameters
  a <- priorBetaAB[1]
  b <- priorBetaAB[2]
  
  # Create summary values of Data:
  z <- sum( Data )      # number of 1's in Data
  N <- length( Data )   # number of trials (flips)
  
  Theta <-              # points for plotting
    seq(0.001,0.999, by=0.001) 
  pTheta <-             # prior pdf of theta
    dbeta( Theta, a, b)  
  pThetaGivenData <-    # posterior pdf of theta
    dbeta( Theta, a+z, b+N-z)  
  pDataGivenTheta <-    # likelihood (sampling distribution)
    Theta^z * (1-Theta)^(N-z)
  
  # Compute the evidence for optional display:
  # pData <- beta(z+a,N-z+b) / beta(a,b)  # may underflow for large a,b
  pData <-          # lbeta() is log( B(a,b) )
    # B(a,b) is Int_0^1 of theta^(a-1)(1-theta)^(b-1).
    exp( lbeta(z+a,N-z+b) - lbeta(a,b) )
  
  ## Plot results.
  layout( matrix( 1:3, nrow=3, ncol=1, byrow=FALSE ) ) # 3x1 panels
  par( mar=c(3,3,1,0),         # number of lines, each margin
       mgp=c(2,0.7,0),         # margin lines for title and axis
       mai=c(0.5,0.5,0.3,0.1)  # margin sizes in inches
  ) 
  cexAxis <- 1.33   # character expansion factor for axis annotation
  cexLab  <- 1.75   # character expansion factor for axis labels
  # convert plotType to notation used by plot():
  if ( plotType=="Points" ) { plotType <- "p" }
  if ( plotType=="Bars"   ) { plotType <- "h" }
  dotsize <- 5 # size of plotted dots
  barsize <- 5 # width of bar lines    
  yLim <- # y limits for plotting prior and posterior:
    c(0, 1.1*max( c(pTheta, pThetaGivenData) ))
  
  plot( # Plot the prior.
    Theta, pTheta, type=plotType, 
    pch=".", cex=dotsize, lwd=barsize,
    xlim=c(0,1), ylim=yLim, cex.axis=cexAxis,
    xlab=bquote(theta), 
    ylab=bquote(dbeta(theta*"|"*.(a),.(b))), 
    cex.lab=cexLab,
    main="Prior (beta)", cex.main=1.5, col="skyblue"
  )
  
  #if ( showCentTend != "None" ) {
  
  if ( showCentTend == "Mean" ) {
    meanTheta <- a/(a+b) 
    if ( meanTheta > .5 ) {
      textx <- 0; 
      textadj <- c(0,1) # left-justify,  top-justify
    } else { # meanTheta ≤ .5
      textx <- 1; 
      textadj <- c(1,1) # right-justify, top-justify
    } # end if ( meanTheta > .5 )
    text( textx, yLim[2],
          bquote( "mean=" * .(signif(meanTheta,3)) ),
          cex=2.0 , adj=textadj )
  } # end if ( showCentTend == "Mean" )
  
  if ( showCentTend == "Mode" ) {
    if ( a+b-2 > 0 ) {
      modeTheta <- (a-1)/(a+b-2)
      if ( modeTheta > .5 ) {
        textx <- 0 ; textadj <- c(0,1)
      } else { # modeTheta ≥ .5
        textx <- 1 ; textadj <- c(1,1)
      } # end if ( modeTheta > .5 )
      text( textx, yLim[2],
            bquote( "mode=" * .(signif(modeTheta,3)) ),
            cex=2.0 , adj=textadj )
    } else { # a+b-2 ≤ 0
      warning("Mode does not exist if a + b ≤ 2.")
    } # end if ( a+b-2 > 0 )
  } # end if ( showCentTend == "Mode" )
  
  # } # end if ( showCentTend != "None" )
  
  # Mark the HDI. HDI points are not thinned in the plot.
  if ( showHDI ) {
    if ( a+b-2 > 0 ) { # a mode exists
      HDIinfo <- HDIofICDF( qbeta, shape1=a, shape2=b, credMass=HDImass )
      HDIheight <- mean( dbeta(HDIinfo, shape1=a, shape2=b) )
      lines( HDIinfo, rep(HDIheight, length(HDIinfo)) )
      text( mean(HDIinfo), HDIheight,
            bquote( .(100*signif(HDImass,3)) * "% HDI" ) ,
            adj=c(0.5,-1.5) , cex=1.5 )
      
      # Mark the left and right ends of the waterline. 
      for ( i in 1:2 ) {
        lines( c(HDIinfo[i], HDIinfo[i]), 
               c(-0.5,HDIheight), type="l" , lty=2 , lwd=1.5 )
        text( HDIinfo[i], HDIheight, bquote(.(round(HDIinfo[i],3))),
              adj=c(0.5,-0.15), cex=1.2 )
      }
    }
  }
  # Mark the ROPE
  if ( !is.null(ROPE) ) {
    #pInRope = ( pbeta( ROPE[2] , shape1=a+z , shape2=b+N-z ) 
    #            - pbeta( ROPE[1] , shape1=a+z , shape2=b+N-z ) )
    pInRope <- (  pbeta( ROPE[2], shape1=a, shape2=b ) 
                  - pbeta( ROPE[1], shape1=a, shape2=b ) )
    ropeTextHt <- 0.7*yLim[2]
    ropeCol <- "darkred"
    lines( c(ROPE[1],ROPE[1]), c(-0.5,ropeTextHt), type="l", lty=2, 
           lwd=1.5, col=ropeCol )
    lines( c(ROPE[2],ROPE[2]), c(-0.5,ropeTextHt), type="l", lty=2, 
           lwd=1.5, col=ropeCol )
    text( mean(ROPE), ropeTextHt,
          paste0(ROPE[1], "<",round(pInRope,4)*100, "%<",ROPE[2]),
          adj=c(0.5,-0.15), cex=1.2, col=ropeCol )    
  }
  
  # Plot the likelihood: p(Data|Theta)
  plot( Theta, pDataGivenTheta, type=plotType, 
        pch=".", cex=dotsize, lwd=barsize,
        xlim=c(0,1), ylim=c(0,1.1*max(pDataGivenTheta)), 
        cex.axis=cexAxis,
        xlab=bquote(theta), ylab=bquote( "p(D|" * theta * ")" ), 
        cex.lab=cexLab,
        main="Likelihood (Bernoulli)", cex.main=1.5, col="skyblue" 
  )
  if ( z > .5*N ) { textx <- 0 ; textadj <- c(0,1) }
  else { textx <- 1 ; textadj <- c(1,1) }
  text( textx, 1.0*max(pDataGivenTheta), cex=2.0,
        bquote( "Data: z=" * .(z) * ",N=" * .(N) ), adj=textadj )
  
  
  if ( showCentTend != "None" ) {
    
    if ( showCentTend == "Mean" ) {
      meanTheta <- z/N
      if ( meanTheta > .5 ) {
        textx <- 0 ; textadj <- c(0,1)
      } else { # meanTheta ≤ .5
        textx <- 1 ; textadj <- c(1,1)
      } # end if ( meanTheta > .5 )
      text( textx, 0.7*max(pDataGivenTheta),
            bquote( "z/N = " * .(signif(meanTheta,3)) ),
            cex=2.0, adj=textadj )
    } # end if ( showCentTend == "Mean" )
    
    if ( showCentTend == "Mode" ) {
      modeTheta <- Theta[ which.max( pDataGivenTheta ) ]
      if ( modeTheta > .5 ) {
        textx <- 0 ; textadj <- c(0,1)
      } else { # modeTheta ≤ .5
        textx <- 1 ; textadj <- c(1,1)
      } # end if ( modeTheta > .5 )
      text( textx, 0.7*max(pDataGivenTheta),
            bquote( "max at " * .(signif(modeTheta,3)) ),
            cex=2.0 , adj=textadj )
    } # end if ( showCentTend == "Mode" )
    
  } # end if ( showCentTend != "None" )
  
  ## Plot the posterior.
  plot( Theta , pThetaGivenData , type=plotType , 
        pch="." , cex=dotsize , lwd=barsize ,
        xlim=c(0,1) , ylim=yLim , cex.axis=cexAxis ,
        xlab=bquote(theta) , 
        ylab=bquote(dbeta(theta*"|"*.(a+z),.(b+N-z))),
        cex.lab=cexLab,
        main="Posterior (beta)" , cex.main=1.5 , col="skyblue" 
  )
  
  if ( showCentTend != "None" ) {
    
    if ( showCentTend == "Mean" ) {
      meanTheta <- (a+z)/(a+b+N) 
      if ( meanTheta > .5 ) {
        textx = 0 ; textadj <- c(0,1)
      } else { # meanTheta ≤ .5
        textx <- 1 ; textadj <- c(1,1)
      } # end if ( meanTheta > .5 )
      text( textx , yLim[2] ,
            bquote( "mean=" * .(signif(meanTheta,3)) ) ,
            cex=2.0 , adj=textadj )
    } # end if ( showCentTend == "Mean" )
    
    if ( showCentTend == "Mode" ) {
      if ( a+z > 1 & b+N-z > 1 ) {
        modeTheta <- (a+z-1)/(a+b+N-2)
        
        if ( modeTheta > .5 ) {
          textx <- 0 ; textadj <- c(0,1)
        } else { # modeTheta ≤ 0.5
          textx <- 1 ; textadj <- c(1,1)
        } # end if ( modeTheta > .5 )
        
        text( textx , yLim[2] ,
              bquote( "mode=" * .(signif(modeTheta,3)) ) ,
              cex=2.0 , adj=textadj )
      } # end if ( a+z > 1 & b+N-z > 1 )
    } # if ( showCentTend == "Mode" )
    
  } # end if ( showCentTend != "None" )
  
  if ( showpD ) { # Plot marginal likelihood pData:
    meanTheta <- (a+z)/(a+b+N)
    if ( meanTheta > .5 ) {
      textx <- 0 ; textadj <- c(0,1)
    } else { # meanTheta ≤ .5
      textx <- 1 ; textadj <- c(1,1)
    } # end if ( meanTheta > .5 )
    text( textx, 0.55*yLim[2], 
          bquote( "p(D)=" * .(signif(pData,3)) ),
          cex=2.0, adj=textadj )
  } # end if ( showpD )
  
  
  if ( showHDI ) { # Mark the highest density interval.
    if ( a+b+N-2 > 0 ) {
      HDIinfo <- HDIofICDF( qbeta , shape1=a+z , shape2=b+N-z , credMass=HDImass)
      HDIheight <- mean( dbeta(HDIinfo,shape1=a+z,shape2=b+N-z) )
      lines( HDIinfo, rep(HDIheight,length(HDIinfo)) )
      text( mean(HDIinfo), HDIheight,
            bquote( .(100*signif(HDImass,3)) * "% HDI" ),
            adj=c(0.5,-1.5) , cex=1.5 )
      # Mark the left and right ends of the waterline. 
      for ( i in 1:2 ) {
        lines( c(HDIinfo[i], HDIinfo[i]) , c(-0.5,HDIheight), 
               type="l", lty=2, lwd=1.5 )
        text( HDIinfo[i], HDIheight, bquote(.(round(HDIinfo[i],3))),
              adj=c(0.5,-0.15), cex=1.2 )
      } # end for ( i in 1:2 )
    } # end if ( a+b+N-2 > 0 )
  } # if ( showHDI )
  
  if ( !is.null(ROPE) ) { # Mark the ROPE
    pInRope <- (  pbeta( ROPE[2], shape1=a+z, shape2=b+N-z ) 
                  - pbeta( ROPE[1], shape1=a+z, shape2=b+N-z ) )
    ropeTextHt <- 0.7*yLim[2]
    ropeCol <- "darkred"
    lines( c(ROPE[1],ROPE[1]) , c(-0.5,ropeTextHt), 
           type="l", lty=2, lwd=1.5 , col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(-0.5,ropeTextHt), 
           type="l", lty=2, lwd=1.5, col=ropeCol )
    text( mean(ROPE) , ropeTextHt ,
          paste0(ROPE[1],"<",round(pInRope,4)*100,"%<",ROPE[2]) ,
          adj=c(0.5,-0.15) , cex=1.2 , col=ropeCol )    
  } # end if ( !is.null(ROPE) )
  
  return( c(a+z,b+N-z) )
} # end BernBeta()