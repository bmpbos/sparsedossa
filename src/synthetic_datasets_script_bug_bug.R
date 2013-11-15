#!/usr/bin/env Rscript
#######################################################################################
# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the MaAsLin (Multivariate Associations Using Linear Models), 
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Timothy Tickle, ttickle@hsph.harvard.edu).
#######################################################################################


c_strBugBugAssociations = "BugToBugAssociations"
c_strMinimumSamples = "Minimum Spiked-in Samples:"
c_strCorrDomainBugs = "Number of bugs each correlated bug is correlated with:"
c_strCorrDomainBugsIdx ="Indices of the bugs each correlated bug is correlated with:"
c_strCorrRangeBugsIdx = "Indices of bugs correlated with others:"
c_strMaxCorrDomainBugs = "Maximum number of bugs with which one bug is correlated:"
c_strNumberOfAssociations = "Number of associations (bugs correlated with others):"
c_strNoiseScaling  = "Scaling parameter for variance of noise:"

funcAssociation = func_linear_association
lAssociationParams = list(dIntercept = 0, vdSlope = 1)

func_get_corr_indices = function(
### A function to obtain the indices of the range and domain bugs
iNumAssociations,
### The number of associations to generate
iMaxDomainNumber,
### The maximum number of features with which one feature can be correlated
iNumFeatures
### The total number of features
){

  # We can choose any index to begin with
  possibleIndices = seq(1,iNumFeatures)

  viIdxCorrRangeBugs = sample(possibleIndices,iNumAssociations)
  possibleIndices = setdiff(possibleIndices,viIdxCorrRangeBugs)

  # Initialize the holder for the domain indices
  liIdxCorrDomainBugs = vector("list",iNumAssociations)
  if(length(possibleIndices) < iNumAssociations){

    iNumAssociations.trunc <- floor(iNumFeatures/2)
    warning( paste( "The number of features",iNumFeatures,
                    "is insufficient to make",iNumAssociations,
                    "associations.  Making the maximum possible associations:",iNumAssociations.trunc ) )

    # Generate the maximum possible associations
    possibleIndices     <- seq(1,iNumFeatures)
    viIdxCorrRangeBugs  <- sample(possibleIndices,iNumAssociations.trunc)
    possibleIndices     <- setdiff(possibleIndices, viIdxCorrRangeBugs)
    liIdxCorrDomainBugs <- as.list( sample(possibleIndices,iNumAssociations.trunc) )

  } else {

    # The current association number and number of associations remaining
    iAssociation   <- 1 
    iNumAssociationsRemaining <- iNumAssociations

    while( length(possibleIndices) > iNumAssociationsRemaining && iAssociation <= iNumAssociations){

      iNumDomainBugs <- sample( seq( 1,
                                     min( iMaxDomainNumber,
                                          length(possibleIndices) ) ),
                                1 )
      liIdxCorrDomainBugs[[iAssociation]] <- sample(possibleIndices,iNumDomainBugs)
      possibleIndices                     <- setdiff(possibleIndices,liIdxCorrDomainBugs[[iAssociation]])
      iAssociation                        <- iAssociation + 1
      iNumAssociationsRemaining           <- iNumAssociationsRemaining - 1
      print(c(iAssociation,iNumAssociationsRemaining))

    }

    if( iAssociation <= iNumAssociations ){

      warning( paste( "Only generated", iNumAssociations - iNumAssociationsRemaining, 
                      "associations with random numbers of domain bugs.  The remaining",iNumAssociationsRemaining,
                      "associations will have 1 domain bug" ) )

      if( length(possibleIndices) < iNumAssociationsRemaining ){

        iAssociation    <- iAssociation - 1
        possibleIndices <- c(possibleIndices,liIdxCorrDomainBugs[[iAssociation]])
      }

      for( i in iAssociation:iNumAssociations ){
        liIdxCorrDomainBugs[[i]] <- sample(possibleIndices,1)
        possibleIndices          <- setdiff(possibleIndices,liIdxCorrDomainBugs[[i]])
      }
    }
  }

  return( list( viIdxCorrRangeBugs  = viIdxCorrRangeBugs,
                liIdxCorrDomainBugs = liIdxCorrDomainBugs ) )
}

func_linear_association = function(
### A function to spike in positive linear bug-bug associations
dfeatures.x,
### The data for the independent features (as a matrix; rows are features)
dIntercept,
### The intercept parameter for the relationship
vdSlope,
### The slope parameter(s) for the relationship, as a vector
### Elements are recycled if there are fewer than 1 per feature
dVarScale
### The scaling parameter for the variance
){
   if(is.matrix(dfeatures.x)) {
      n.features <- nrow(dfeatures.x)
      if( length(vdSlope) < n.features ) vdSlope <- rep( vdSlope,ceiling( n.features/length( vdSlope ) ) )
      vdfeature.y = dIntercept + apply(vdSlope[1:n.features]*dfeatures.x,2,sum) + dVarScale * sum( apply( dfeatures.x,1,var ) )
      vdSlopeUsed = vdSlope[1:n.features]
   }
   vdfeature.y = dIntercept + vdSlope[1] * dfeatures.x + dVarScale
   vdSlopeUsed = vdSlope[1]
   return( list(vdfeature.y = vdfeature.y,
                vdSlopeUsed = vdSlopeUsed) )
}