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