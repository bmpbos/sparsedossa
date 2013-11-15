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

### Not formally tested                                                                                                                                                                                          
func_generate_bug_bug_spikes = function(
### Add spiked bug-bug correlations                                                                                                                                                                              
mtrxData,
### The matrix into which to spike the associations                                                                                                                                                              
dVarScale,
### The scaling parameter for the VARIANCE (the noise added will have variance dVarScale*var(bug))                                                                                                               
funcAssociation,
### The function to be used to generate the association; takes two bugs and the variance scaling parameter                                                                                                       
lAssociationParams,
### The possible extra parameters for the association function                                                                                                                                                   
viIdxCorrRangeBugs,
### row indices of the range bugs as a vector                                                                                                                                                                    
liIdxCorrDomainBugs,
### The row indices of the domain bugs, arranged as a list                                                                                                                                                       
iDataset,
### The dataset number this matrix corresponds to                                                                                                                                                                
vdPercentZero = NA,
### Vector of percent zero parameters for features if not supplied, one will be generated by rlnorm                                                                                                              
fZeroInflate = TRUE,
### Controls if zero inflation is used.                                                                                                                                                                          
fVerbose = FALSE
### If true, plotting and logging occur                                                                                                                                                                          
){
  print("start func_generate_bug_bug_spikes")
  # Define some parameters                                                                                                                                                                                       
  iNumAssociations           <- length( liIdxCorrDomainBugs )
  iNumAssociationsToGenerate <- 0
  iNumFeatures               <- nrow( mtrxData )
  viNumberCorrDomainBugs     <- unlist( lapply( liIdxCorrDomainBugs, "length" ) )

  lAssociationParams$dVarScale <- dVarScale

  # Some generic error reporting                                                                                                                                                                                 
  strWarning = "No correlation was added to the bug-bug correlation matrix."

  if( length(viIdxCorrRangeBugs) < iNumAssociations ){
    print(paste( "Number of range relationships ",
                 length(viIdxCorrRangeBugs),
                 " is less than the number of domain relationships ",
                 iNumAssociations,
                 ". ",
                 strWarning,
                 sep="" ))
    iNumAssociations = 0
  }
  if( length(viIdxCorrRangeBugs) > iNumAssociations ){
    print(paste( "Number of range relationships ",
                 length(viIdxCorrRangeBugs),
                 " is greater than the number of domain relationships ",
                 iNumAssociations,
                 ". ", strWarning,
                 sep="" ))
    iNumAssociationsToGenerate <- length(viIdxCorrRangeBugs) - iNumAssociations
  }
  if( min(viNumberCorrDomainBugs) <=0 ){
    print(paste( "Some number of domain bugs", min(viNumberCorrDomainBugs), "is not possible.", strWarning ))
    iNumAssociations = 0
  }
  if( dVarScale < 0 ){
    paste( "The variance scale parameter,",dVarScale,", is < 0. ", strWarning, sep="" )
    iNumAssociations = 0
  }
  if( iNumAssociations > floor(iNumFeatures/2) ){
     print(paste( iNumAssociations,"associations cannot be formed with only", int_number_features,
                  "features: The number of associations must be less than half the number of features" ) )
    iNumAssociations = 0
  }

  if(iNumAssociations > 0){
    ## Converting vectors and lists to appropriately delimited strings                                                                                                                                           
    strNumberCorrDomainBugs = paste(viNumberCorrDomainBugs[1])
    if(length(viNumberCorrDomainBugs)>1){

      for(k in 2:length(viNumberCorrDomainBugs)){
        strNumberCorrDomainBugs = paste(strNumberCorrDomainBugs,viNumberCorrDomainBugs[k],sep='; ')
      }

    }

    strIdxCorrRangeBugs = paste(viIdxCorrRangeBugs[1])
    if(length(viIdxCorrRangeBugs)>1){
      for(k in 2:length(viIdxCorrRangeBugs)){
        strIdxCorrRangeBugs = paste(strIdxCorrRangeBugs,viIdxCorrRangeBugs[k],sep='; ')
      }
    }

    strIdxCorrDomainBugs = toString(liIdxCorrDomainBugs[[1]])
    if(length(liIdxCorrDomainBugs) > 1){
      for(k in 2:length(liIdxCorrDomainBugs)){
        strIdxCorrDomainBugs = paste(strIdxCorrDomainBugs,toString(liIdxCorrDomainBugs[[k]]),sep='; ')
      }
    }
    ## End converting to strings                                                                                                                                                                                 

    ## Generate each association                                                                                                                                                                                 
    mtrx_final <- mtrxData
    for(i in seq(1,iNumAssociations)){

      lAssociationParams$dfeatures.x   <- mtrxData[ liIdxCorrDomainBugs[[ i ]], ]
      vdfeature.y                      <- do.call( what = funcAssociation,
                                                   args = lAssociationParams )
      mtrx_final[viIdxCorrRangeBugs[i],] <- vdfeature.y$vdfeature.y
    }
    ## End generating associations                                                                                                                                                                               

  } else {
    viNumCorrDomainBugs = NA
    viCorrRangeBugsIdx  = NA
    liCorrDomainBugsIdx = NA

    strNumCorrDomainBugs = "NA"
    strCorrRangeBugsIdx  = "NA"
    strCorrDomainBugsIdx = "NA"
  }

  # This will hold the associations that you create and be placed in the truth file that is records association spike-ins for later assessment                                                                   
  # It starts with the name of the microbiome you are creating                                                                                                                                                   
  # Parameters of interest and then your feature associations                                                                                                                                                    
  mtrxParameters = matrix(data=NA, nrow=11, ncol=1)

  mtrxParameters[1,1]  = paste(c_strSyntheticMicrobiome,   c_strBugBugAssociations, "_d_", iDataset, sep='')
  mtrxParameters[2,1]  = paste(c_strNumberOfFeatures,      int_number_features)
  mtrxParameters[3,1]  = paste(c_strNumberOfSamples,       int_number_samples)
  mtrxParameters[4,1]  = paste(c_strNumberCounts,          iMinNumberCounts ) )
  mtrxParameters[5,1]  = paste(c_strMinimumSamples,        iMinNumberSamples ) )
  mtrxParameters[6,1]  = paste(c_strNoiseScaling,          dVarScale ) )
  mtrxParameters[7,1]  = paste(c_strNumberOfAssociations,  iNumAssociations ) )
  mtrxParameters[8,1]  = paste(c_strMaxCorrDomainBugs,     iMaxNumberCorrDomainBugs ) )
  mtrxParameters[9,1]  = paste(c_strCorrDomainBugs,        strNumCorrDomainBugs ) )
  mtrxParameters[10,1] = paste(c_strCorrRangeBugsIdx,      strCorrRangeBugsIdx ) )
  mtrxParameters[11,1] = paste(c_strCorrDomainBugsIdx,     strCorrDomainBugsIdx ) )

  print("stop func_generate_bug_bug_spikes_matrix")

  return(list( mtrxAssnParameters = mtrxParameters
               mat_bugs           = mtrxFinal))
}