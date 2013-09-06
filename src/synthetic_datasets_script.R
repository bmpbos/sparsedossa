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


# This is the core script for generating synthetic datasets. This will call each function in turn to generate the metadata, and the bugs, and combine them into the appropriate pcl files.

library(optparse)
library(getopt)
library(MASS)

# Plot constants
c_strDefaultMarkerColor = "black"
c_strOutlierColor = "red"
c_strFeatureOutlier = "cyan"

# Control flow constants
iLoopingControlIncrement = 1000
### Max number of looping in looping functions
c_iCountTypesOfMetadata = 4
### There are 4 types of metadata supported in metadata generation
c_fIgnoreZerosInOutliers = TRUE
### Will not allow zeros to be used as the min value in swapping unless they are needed to fulfill the number of
### swaps (if there are a whole bunch of zeros, some zeros may be needed or no swapping can be performed).
c_fDummyFactorData = TRUE
### For Factor data only, selects a level and sets all of that level to 1 and others to 0
c_dRunifMin = .01
c_dRunifMax = .99
### Binary and quarternary metadata are draw from uniform distributions, these are the bounds
c_dMinBinary = .3
c_dMaxBinary = .7
### Min and max probabilities for binary metadata levels
c_iTimesSDIsOutlier = 4
### The number of deviations until a value is truncated as an outlier
c_dALittleMoreThanZero = 0.00001

# Labeling Constants
c_strLevel = "Level"
c_strFeature = "Feature"
c_strMetadata = "Metadata"
c_strRandom = "Lognormal"
c_strOutlier = "Outlier"
c_strSpike = "spike"
c_strBugBugAssocations = "BugToBugAssociations"
c_strMinimumSamples = "Minimum Spiked-in Samples:"
c_strSyntheticMicrobiome = "SyntheticMicrobiome"
c_strNumberOfFeatures = 'Number of features:'
c_strNumberOfSamples = 'Number of samples:'
c_strPercentSpikes = 'Percent spikes:'
c_strMultiplier = 'Multiplier:'
c_strMultiplierParameter = 'Multivariate Parameter:'
c_strTotalSampleBugOccurrence = "Total Reads per Sample:"
c_strNumberCounts = "Minimum Number of Counts:"
c_strNumberSamples = "in Minimum Number of Samples:"
c_strPercentOutliers = "Max Percent Outliers in a Sample:"
c_strPercentSampleOutliers = "Percent Samples with Outliers:"
c_strOutlierParameter = "Outlier Swap:"
c_strSampleParameter = "Sample:"
c_strContinuous = "Continuous"
c_strFactor = "Factor Levels"
c_strMetadataDetails = "Metadata: Details"

### For Bug-bug spikin matrix (in alphabetical order)
c_strCorrDomainBugs = "Number of bugs each correlated bug is correlated with:"
c_strCorrDomainBugsIdx ="Indices of the bugs each correlated bug is correlated with:"
c_strCorrRangeBugsIdx = "Indices of bugs correlated with others:"
c_strMaxCorrDomainBugs = "Maximum number of bugs with which one bug is correlated:"
c_strNumberOfAssociations = "Number of associations (bugs correlated with others):"
c_strNoiseScaling  = "Scaling parameter for variance of noise:"

# Temporary control flags
c_dfFreezeSDFeatures = FALSE
c_dfFreezeSDGrandMu = FALSE
c_fPrintLognormalMatrix = TRUE

# Define the relatioship between different feature properties
# These variables are associated with settings for
# Calculating the SD and percent Zero based on the mu or expectation.
# The constants for these variables have been estimated from
# a real (IBD) data set. These constants are used unless a
# calibration file is given which estimates the values in the
# same manner as the constants were estimated.
c_dSDBeta = 0.242678502955986
### The estimate for the relationship between SD and exp
c_dBetaZero = -0.09194 # -0.119637999937095 # .86160
### The extimate for the relationship between exp and zero percent
c_dBetaGrandSD = 0.01699164
### The estimate for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)

################### Methods


# Plotting function. Not in regression suite
funcBoxPlotOutliers = function(
### Generic plotter for a series of boxplots.
mtrxData,
### Matrix of values
strTitle,
### Title for the plot
lviHighlight = list(),
### Optional. Points to highlight, list of vectors of indices. Each entry is to be highlighted, each list vector is a sample. Default no highlighting.
strColor = "red",
### Optional. Highlighting color. Default "red".
fWithZeros = TRUE,
### Optional. Plot with zeros. Default TRUE.
fBySamples = TRUE
### Optional. Indicates if the matrix columns are samples. By default (TRUE) samples are expected as columns. One can give FALSE if they want the matrix transposed.
){
  strXLabel = "Samples"
  if( !fBySamples ){
    strXLabel = "Features"
    mtrxData = t(mtrxData)
  }
  if( fWithZeros ){
    boxplot( mtrxData, main = paste( strTitle, " - With Zeros - By ", strXLabel, sep = "" ), xlab = strXLabel, ylab = "Distribution of Counts" )
    lvdData = list()
    viPosition = c()
    for( iCol in 1:ncol( mtrxData ) ){
      lvdData[[length( lvdData )+1]] = mtrxData[, iCol]
      viPosition = c( viPosition, length( lvdData ) )
      lstrColoring = rep( c_strDefaultMarkerColor, length( mtrxData[,iCol] ) )

      if( !is.null( lviHighlight[iCol][[1]] ) ){
        lstrColoring[lviHighlight[[iCol]]] = strColor
      }

      # Highlight for the different levels (highlight and not)
      for(strLevel in levels( as.factor( lstrColoring ) ) ){
        stripchart( mtrxData[, iCol][which( lstrColoring == strLevel )], at = iCol, add = TRUE, vertical = TRUE, method = "jitter", jitter = 0.3, col = strLevel )
      }
    }
 } else {
    liNoZeros = list()
    for( iColumn in 1:ncol( mtrxData ) ){
      viNotZeros = which( mtrxData[, iColumn] != 0 )
      liNoZeros[[length( liNoZeros ) + 1]] = mtrxData[, iColumn][viNotZeros]

      fTrimHighlights = !is.null( lviHighlight[iColumn][[1]] )
      # Trim down colors to the same indices as non zeros
      if( fTrimHighlights ){
        lviHighlight[[iColumn]] = intersect( lviHighlight[[iColumn]], viNotZeros )
        # Shift the indices down to account for removing zeros
        viShiftedIndices = c()
        for( iIndex in lviHighlight[[iColumn]] ){
          viShiftedIndices = c( viShiftedIndices, 1 + length( which( viNotZeros < iIndex ) ) )
        }
        if( is.null( viShiftedIndices ) ){
          lviHighlight[[iColumn]][1] = NULL
        } else {
          lviHighlight[[iColumn]] = viShiftedIndices
        }
      }
    }
    boxplot( liNoZeros, main = paste( strTitle, " - Without Zeros - By ", strXLabel, sep = ""), xlab = strXLabel, ylab = "Distribution of Counts" )
    for( iColumn in 1:ncol( mtrxData ) ){
      liNoZeros[[length( liNoZeros ) + 1]] = mtrxData[, iColumn][which( mtrxData[, iColumn] != 0 )]
      lstrColoring = rep( c_strDefaultMarkerColor, length( mtrxData[, iColumn] ) )
      if( !is.null( lviHighlight[iColumn][[1]] ) ){
        lstrColoring[lviHighlight[[iColumn]]] = strColor
      }

      # Highlight for the different levels (highlight and not)
      for( strLevel in levels( as.factor( lstrColoring ) ) ){
        stripchart( liNoZeros[[iColumn]][which( lstrColoring == strLevel )], at = iColumn, add = TRUE, vertical = TRUE, method = "jitter", jitter = 0.3, col = strLevel )
      }
    }
  }
}


# Reviewed but not tested.
funcCalibrateRLNormToMicrobiome = function(
### Given a TSV file
### Parameters for distribution generatio are given
### To be used in matrix generation
### Estimated parameters include
### SD excluding zeros
### Mus excluding zeros
### Percent zeros
### The beta for estimating the SD given the mu
### The beta for estimating the percent zero given the mu
### All Mus, SD, grand MU, and grand SD are ready for rlnorm (have been measured from a logged (rlnorm)
sCalibrationFile,
### File to be read in and used to calibrate constansts and relationships in the underlying data.
fVerbose = FALSE
### Flag to turn on logging and pdf creation
){
  print("start funcCalibrateRLNormToMicrobiome")
  # Read in file
  print("Reading file.")
  dfData = read.table(sCalibrationFile)
  row.names(dfData) = dfData[[1]]
  dfData = dfData[-1,-1]

  # Get read depth of the samples (From a tsv file, samples = rows)
  ldReadDepths = sapply(1:nrow(dfData), function(x) sum(as.numeric(as.matrix(dfData)[x,])))

  # Get the vector of Exps, Mus and SD ignoring 0s
  # Logged will be used to make features, they can be directly used in the rlnorm function
  # Not logged will be used to later be logged and estimate a grand mean and grand SD for the initial distribution,
  # Given that every point in the larger vector is the mu of the feature vectors.
  # Which could then be used by rlnorm to generate a each feature if needed.
  # This is not needed here but the calculation of this value is useful. It is used as an initial value for synthetic creation
  # of the initial vector of feature mus.

  # Mean of the nonzero data
  vdExp = c()
  # Mean of the logged nonzero data
  vdMu = c()
  # Standard deviation of the logged nonzero data
  vdLogSD = c()
  # Percent zeros in data
  vdPercentZero = c()

  # Calculate parameters for each feature (column)
  for(iIndex in 1:ncol(dfData))
  {
    # Get the percent zero before removing zeros for other measurements
    vdCur = as.numeric(as.vector(as.matrix(dfData[iIndex])))
    vdPercentZero = c(vdPercentZero, length(which(vdCur==0))/length(vdCur))

    # Remove zeros from data
    vdExp = c(vdExp,mean(vdCur))
    vdCur = vdCur[which(vdCur!=0)]
#    vdExp = c(vdExp,mean(vdCur))
    
    #### Note
    #### rlnorm needs a mean and sd from a logged rlnorm distribution which would match this
    #### without further manipulation. The "mean" in the formula is actually not the expectation
    #### The expectation is e^mu+.5*sd^2 so this is always more than mu.

    # Log nonzero data
    vdLogCur = log(vdCur)
    vdLogSD = c(vdLogSD, sd(vdLogCur))
#    vdMu = c(vdMu, mean(vdLogCur))
    vdMu = c( vdMu, funcGetMu( vdExp[iIndex], exp(sd( vdLogCur ) ) ) )
  }

  # Estimate the distribution parameters from the expectation vector
  # Includes the relationship between the grand mu and grand sd
  # The grand mu, grand expectation (logged) and the grand sd
  lParams = funcGenerateExpVectorParameters(vdExp,TRUE)

  ##### Get relationship between logSD and log(Exp)
  # Log to make a linear relationship
  # This is both values logged as SD is based on vdLogCur
  vdExpLog = log(vdExp)

  lmod = lm(as.formula(paste("vdLogSD","vdExpLog",sep="~")))
  dBetaSD = coef(lmod)["vdExpLog"]

  #### Percent Zero and Exp
  ### When using NLS, first supply with a guess (took from a lm)
  ### Then take the LM guess and do the NLS
  dBetaZero = NA
  if(sum(vdPercentZero)>0)
  {
    # Make a first guess
    lmod = lm(vdPercentZero ~ vdExpLog)
    dBetaZero = coef(lmod)["vdExpLog"]

    print("lmod")
    print(lmod)

    # Using nonlinear fit and giving the linear fit as the first choice.
#    modNLS = nls(vdPercentZero ~ vdExpLog*dBetaZZ+dBetaInter,start=list(dBetaZZ=dBetaZero))
#    print("modNLS")
#    print(modNLS)
#    print(coef(summary(modNLS)))
#    dBetaZero = coef(summary(modNLS))["dBetaZZ.vdExpLog","Estimate"]
  }

  if(fVerbose)
  {
    # Initial description of calibration data set
    barplot(ldReadDepths, main=paste("Read depth of Template Microbiome, Average",mean(ldReadDepths)))
    abline(mean(ldReadDepths),0,col="red")
    hist(as.numeric(as.vector(as.matrix(dfData))), main="Histogram of Template Microbiome Counts")
    hist(log(as.numeric(as.vector(as.matrix(dfData)))), main="Histogram of Logged Template Microbiome Counts")
    plot(as.numeric(as.vector(as.matrix(dfData))), main="Distribution of Template Microbiome")

    # Get relationship between SD and MU
    plot(vdExp, vdLogSD, main="Estimating SD with Exp: baseline")
    points(x=vdExp,y=funcEstimateFeatureSD(vdExpLog,dBetaSD),col="orange")
    plot(vdExpLog, vdLogSD, main="Estimating SD with Log Exp: predicted points")
    points(x=vdExpLog,y=funcEstimateFeatureSD(vdExpLog,dBetaSD),col="green")

    # Percent zero #!#
    plot(vdExp, vdPercentZero, main="Estimating Percent Zero: Original baseline relationship 1")
    points(x=vdExp,y=funcEstimatePercentZero(vdExpLog,dBetaZero),col="violet")
    plot(vdExpLog, vdPercentZero, main="Estimating Zero with Exp: NLS predicted points 1")
    points(x=vdExpLog,y=funcEstimatePercentZero(vdExpLog,dBetaZero),col="violet")
  }

  print("stop funcCalibrateRLNormToMicrobiome")

  return(list(exp=vdExp, mu=vdMu, sd=exp(vdLogSD), dSDBeta=dBetaSD, percentZero=vdPercentZero,
              dZeroBeta=dBetaZero, dGrandBeta=lParams$GrandLogSDBeta,
              dAverageReadDepth=mean(ldReadDepths), iFeatureCount=ncol(dfData)))

### When returning the grand Mu remember that you are returning the Mu that gives the expectation for the Mus
### given the rlnorm function so this is different than the mus measured in the logged distribution (rlnorm)
### exp: Not Logged expectation of distribution ( mean(x) )
### mu: The exponentiated logMu
### sd: The exponentiated logSD
### dSDBeta: Relationship between logSD and log(mean(x))
### percentZero: Percent of values which are zero
### dZeroBeta: Relationship between percent zero and log(mean(x))
### dGrandBeta: Relationship between the logSD and Exp
### dAverageReadDepth: The average read depth per samples
### iFeatureCount: The number of features measured
}


# 5 Tests 9/4/2013
funcCycleVectorIndices = function(
### Modified modulus function that returns iElementsToPick  number of indices to use to cycle through a
### vector of data a certain amount of times. Allows a vector to acts like 
### a list of length iVectorLength with ends linked together.
iElementsToPick,
### The number of indices or elements needed
iVectorLength
### The total number of the vector being cycled over. Will often be shorter than the iElementsToPick.
){
  viReturn = c()
  for( iElements in 1:iElementsToPick ){
    viReturn = c( viReturn, ( iElements %% iVectorLength ) )
  }
  viReturn[viReturn == 0] = iVectorLength
  return( viReturn )
}


# 5 Tests 9/4/2013
funcEstimateFeatureSD=function(
### Estimate the SD given the Log Exp and parameters modeling the relationship between the Log Exp and the SD
vdExpLog,
### The measured mean of feature values
dBetaSD
### The beta for the relationship between the mu and the SD
){
  return(vdExpLog * dBetaSD)
  ### Returns the LogSD for a feature
}


# 5 Tests 9/4/2013
funcEstimateGrandSD=function(
### Estimate the grand SD given the grand Mu and parameters modeling the relationship between the grand mu and grand sd
dExp,
dBetaGrandSD
){
  return(dExp * dBetaGrandSD)
  ### Returns the LogSD for the expecation vector
}


# 5 Tests 9/4/2013
funcEstimatePercentZero = function(
### Estimate the percent zero given the logged expectation and parameters modeling the relationship between the log Exp and the percent zero
vdExpLog,
### The measured mean of feature values
dBetaZero
){
  return((vdExpLog * dBetaZero)+.70)
}


# 2 Tests 9/4/2013
funcGenerateExpVectorParameters = function(
### Get the point estimate for the relationship between the mu and sd
vdExpectations,
### These are the mus of each of the features from the calibration file.
### These are untransformed (not logged) and the expectation
### This will not be a sparse/zero-inflated distribution given that it is the expectation of
### each feature.
fVerbose = FALSE
### If pdfs should be made
){
  print("Start funcGenerateExpVectorParameters")

  #Point estimates of the distributions
  #The exp you should get if you use these SD and Mu directly in the rlnorm
  # Not logged
  dExp = mean(vdExpectations)
  #This gives the logSD value
  dLogSD = sd(log(vdExpectations))

  # Get the Mu that would generate the Exp given the SD
  # Already logged (logMu)
  dLogMu = log(funcGetMu(dExp, exp(dLogSD)))

  # The relationship between the Exp and the LogSD
  dGrandBeta = dLogSD / dExp

  return(list(GrandLogMu=dLogMu, GrandExp=dExp, GrandLogSD=dLogSD, dReadDepth=sum(vdExpectations), GrandLogSDBeta=dGrandBeta))
### GrandMu: The logMus of the distributions
### GrandExp: The expectation of the distributions mean(x)
### GrandSD: The logSDs of the distributions
### dReadDepth: The total average read depth of the samples
### GrandSDBeta: The relationship between the logSD and the Exp
}

#!# TODO update exp business
func_generate_bug_bug_spiking_matrix = function(
### Stub for bug bug interactions as a matrix. This just gives a zero inflated lognormal distributrion.
### You'll need to add correlations
int_number_features,
### Number of features
int_number_samples,
### Number of samples
iMinNumberCounts,
### Minimum number of counts for a feature to be considered in a sample for the QC filtering
iMinNumberSamples,
### Created bugs must have a minimum number of samples (iMinNumberSamples ) that have a minimum number of counts (iMinNumberSamples)
iReadDepth,
### Simulated read depth
dVarScale,
### The scaling parameter for the VARIANCE (the noise added will have variance dVarScale*var(bug))
iNumAssociations,
### The number of correlation "structures" to introduce (ie. the number of bugs which will be correlated with others)
iMaxNumberCorrDomainBugs,
### The maximum number of bugs each of the iNumAssociations correlated bugs can be correlated with; the minimum is 1
vdExp = NA,
### Vector of expectation for the features (should add up to read depth)
vdMu = NA,
### Vector oflog Mu parameters for the original mu distribution (means of features) if not supplied, one will be generated by rlnorm
vdSD = NA,
dPercentZero = NA,
fZeroInflate = TRUE,
### Controls if zero inflation is used.
dBetaSD = c_dBetaSD,
### The beta for the relationship between the mu and sd in which generates the expectation vector
dBetaGrandSD = c_dBetaGrandSD,
### The beta for the relationship between the Exp and the SD
dBetaZero = c_dBetaZero,
### The beta for the relationship between the Exp and the percent zero
fVerbose = FALSE
### If true, plotting and logging occur
){
  print("start func_generate_bug_bug_spiking_matrix")
  # Some generic error reporting
  strWarning = "No correlation was added to the bug-bug correlation matrix."
  if( iNumAssociations < 0 ){
    print(paste( "The number of associations given, ", iNumAssociations, ", is < 0. ", strWarning,sep="" ))
  }
  if( iMaxNumberCorrDomainBugs <=0 ){
    print(paste( "Number of domain bugs", iMaxNumberCorrDomainBugs, "is not possible.", strWarning ))
    iNumAssociations = 0
  }
  if( dVarScale < 0 ){
    paste( "The variance scale parameter,",dVarScale,", is < 0. ", strWarning, sep="" )
    iNumAssociations = 0
  }
  if( iNumAssociations > floor(int_number_features/2) ){
     print(paste( iNumAssociations,"associations cannot be formed with only", int_number_features,
                  "features: The number of associations must be less than half the number of features" ) )
    iNumAssociations = 0
  }

  ### Generating the indices for all bugs concerned
  if( iNumAssociations > 0 ){

    viAvailableIndices  = seq(1,int_number_features)                                              # Can use any index for domain or range at initialization
    liCorrDomainBugsIdx = list()                                                                  # List of the indices of the bugs in the domain
    viCorrRangeBugsIdx  = sample(viAvailableIndices,iNumAssociations,replace=FALSE)               # The indices of the bugs in the range
    viAvailableIndices  = viAvailableIndices[-which(viAvailableIndices %in% viCorrRangeBugsIdx)]  # Remove indices already selected from the available pool

    viNumCorrDomainBugs = sample( seq(1,iMaxNumberCorrDomainBugs),
                                  iNumAssociations, replace=TRUE )  # Vector of the number of bugs with which each of the correlated bugs will be correlated
      
    # Get the domain indices for each association
    for(k in 1:iNumAssociations){
      if( length(viAvailableIndices) >= viNumCorrDomainBugs[k] ){                     # Check that the number of remaining indices is sufficient
         liCorrDomainBugsIdx[[k]] = sample( viAvailableIndices,
                                            viNumCorrDomainBugs[k],
                                            replace=FALSE )
         viAvailableIndices       = viAvailableIndices[-which(viAvailableIndices %in% liCorrDomainBugsIdx[[k]])]
      } else {                                                                        # Otherwise stop and print an error message
        stop( paste( "Insufficient Features: there remain only",length(viAvailableIndices),
                     "feature(s) left to generate domain bugs, but the number of domain bugs needed is",viNumCorrDomainBugs[k] ) )
      }
    }

    ## Converting vectors and lists to appropriately delimited strings
    strNumCorrDomainBugs = paste(viNumCorrDomainBugs[1])
    if(length(viNumCorrDomainBugs)>1){
      for(k in 2:length(viNumCorrDomainBugs)) strNumCorrDomainBugs = paste(strNumCorrDomainBugs,viNumCorrDomainBugs[k],sep='; ')
    }

    strCorrRangeBugsIdx = paste(viCorrRangeBugsIdx[1])
    if(length(viCorrRangeBugsIdx)>1){
      for(k in 2:length(viCorrRangeBugsIdx)) strCorrRangeBugsIdx = paste(strCorrRangeBugsIdx,viCorrRangeBugsIdx[k],sep='; ')
    }

    strCorrDomainBugsIdx = toString(liCorrDomainBugsIdx[[1]])
    if(length(liCorrDomainBugsIdx) > 1){
      for(k in 2:length(liCorrDomainBugsIdx)) strCorrDomainBugsIdx = paste(strCorrDomainBugsIdx,toString(liCorrDomainBugsIdx[[k]]),sep='; ')
    }
    ## End converting to strings

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
  vStrParameters = c(paste(c_strSyntheticMicrobiome, c_strBugBugAssocations, sep='_' ) )
  vStrParameters = c(vStrParameters, paste(c_strNumberOfFeatures,      int_number_features ) )
  vStrParameters = c(vStrParameters, paste(c_strNumberOfSamples,       int_number_samples ) )
  vStrParameters = c(vStrParameters, paste(c_strNumberCounts,          iMinNumberCounts ) )
  vStrParameters = c(vStrParameters, paste(c_strMinimumSamples,        iMinNumberSamples ) )
  vStrParameters = c(vStrParameters, paste(c_strNoiseScaling,          dVarScale ) )
  vStrParameters = c(vStrParameters, paste(c_strNumberOfAssociations,  iNumAssociations ) )
  vStrParameters = c(vStrParameters, paste(c_strMaxCorrDomainBugs,     iMaxNumberCorrDomainBugs ) )
  vStrParameters = c(vStrParameters, paste(c_strCorrDomainBugs,        strNumCorrDomainBugs ) )
  vStrParameters = c(vStrParameters, paste(c_strCorrRangeBugsIdx,      strCorrRangeBugsIdx ) )
  vStrParameters = c(vStrParameters, paste(c_strCorrDomainBugsIdx,     strCorrDomainBugsIdx ) )


  # Here is a method that will get you a mu_vector with default pushed through the script from the sfle call
  mu_vector = funcGenerateFeatureParameters( int_number_features = int_number_features, 
                                       int_number_samples = int_number_samples, 
                                       iMinNumberSamples = iMinNumberSamples, 
                                       iReadDepth = iReadDepth,
                                       vdExp = vdExp,
                                       vdMu = vdMu, 
                                       vdSD = vdSD, 
                                       vdPercentZero = dPercentZero, 
                                       dBetaSD = dBetaSD, 
                                       dBetaZero = dBetaZero,
                                       dBetaGrandSD = dBetaGrandSD,
                                       fVerbose = fVerbose)

  # This gets a lognormal distribution, you can do this with a predefined mu vector if needed.
  # After this step you should add the bug spike-ins
  # Note this is count data, normalization will happen automatically for you after you return the matrix
  mtrxBugs = func_generate_random_lognormal_matrix( int_number_features = int_number_features, 
                                                    int_number_samples = int_number_samples, 
                                                    iMinNumberCounts = iMinNumberCounts, 
                                                    iMinNumberSamples = iMinNumberSamples, 
                                                    iReadDepth = iReadDepth,
                                                    vdExp = mu_vector[["exp"]],
                                                    vdMu = mu_vector[["mu"]],
                                                    xPercentZero = mu_vector[["PercentZero"]],
                                                    vdSD = mu_vector[["sd"]],
                                                    fZeroInflate = fZeroInflate, 
                                                    dBetaSD = dBetaSD, 
                                                    dBetaZero = dBetaZero,
                                                    dBetaGrandSD = dBetaGrandSD,
                                                    fVerbose = fVerbose )[["mat_bugs"]]

  # Generating the correlation structure (very simple at the moment)
  if( iNumAssociations > 0 ){
    for(i in seq(1,iNumAssociations)){

      vdVarCorrDomainBugs = apply(mtrxBugs,1,var)[liCorrDomainBugsIdx[[i]]]  # Get the variance of the domain bugs
      if(length(liCorrDomainBugsIdx[[i]])==1){
        mtrxBugs[viCorrRangeBugsIdx[i],] = mtrxBugs[liCorrDomainBugsIdx[[i]],]+
                                           rnorm( int_number_samples,mean=0,sd = sqrt( dVarScale*sum( vdVarCorrDomainBugs ) ) )

      } else {
        # TODO Emma check, I was messin with your code here.
        vdCurDomainBug = mtrxBugs[liCorrDomainBugsIdx[[i]],]
        if(is.null(dim(vdCurDomainBug))){
          vdCurDomainBug = sum(vdCurDomainBug)
        } else {
          vdCurDomainBug = apply(mtrxBugs[liCorrDomainBugsIdx[[i]],],2,sum)
        }
        mtrxBugs[viCorrRangeBugsIdx[i],] = vdCurDomainBug + rnorm( int_number_samples,mean=0,sd = sqrt( dVarScale*sum( vdVarCorrDomainBugs ) ) )
      }
    }
  }
  print("stop func_generate_bug_bug_spiking_matrix")
  # This is the minimum return.
  return(list(mtrxBugs=mtrxBugs, vStrParameters=vStrParameters))
}


func_generate_lefse_matrices = function(
lefse_file,
metadata_parameters,
int_number_features,
### Number of features
int_number_samples,
### Number of samples
mat_metadata,
mat_bugs,
dataset_type
){
  # only have discrete metadata to pick from
  # identify which matrices are discrete from the truth file
  # Get the indices of factor metadata
  lsIndicesOfMicrobiomes = which(sapply(as.vector(metadata_parameters), function(x) grepl('Factor',x))==TRUE)

  l_lefse_matrix = list()
  lefse_matrix_1 = matrix(data=NA, nrow=(int_number_features+2), ncol=(int_number_samples+1))
  for (meta_index in lsIndicesOfMicrobiomes)
  {
    lefse_matrix_1[1,1] = 'SampleID'
    lefse_matrix_1[1,2:(int_number_samples + 1)] = paste('Sample',1:int_number_samples,sep='')
    lefse_matrix_1[2,1] = 'Metadata'
    lefse_matrix_1[2,2:(int_number_samples+1)] = mat_metadata[(meta_index-1),]
    lefse_matrix_1[3:(int_number_features+2),1] = paste("Bug", 1:int_number_features, sep='')
    lefse_matrix_1[3:(int_number_features+2),2:(int_number_samples+1)] = mat_bugs
    lefse_filenames = file.path(lefse_file, paste('LefSe_',dataset_type, '_metadata_',(meta_index-1), '.1-pcl', sep=''))
    write.table(lefse_matrix_1,file=lefse_filenames, quote=FALSE, row.names=FALSE, col.names=FALSE,sep='\t')
    l_lefse_matrix[[length(l_lefse_matrix)+1]] = lefse_matrix_1
  }
  return(l_lefse_matrix)
}

# 32 Tests 9/4/2013
# Reviewed but not tested
funcGenerateFeatureParameters = function(
### Generate the feature mu paramater vector and associated SD parameter, expectation, and percent zeros if needed
# If a mu parameter vector is given of the size of the int_number_samples then pass through
# Otherwise sample to that size with replacement.
int_number_features,
### Number of features
int_number_samples,
### Number of samples
iMinNumberSamples,
### Minimum number of samples not zero
iReadDepth,
### Simulated read depth
vdExp = NA,
### Vector of expectation for the features (should add up to read depth)
vdMu = NA,
### Vector of Mu parameters for the original expectation distribution (means of features) if not supplied, one will be generated by rlnorm
vdSD = NA,
### The vector of SD parameters matching the vdMus.
vdPercentZero = NA,
### The vector of percent zeros matching the vdMus
dBetaSD = c_dBetaSD,
### If vdSD is not given, SD will be generated by a relationship with Exp parameters using the beta given here
dBetaZero = c_dBetaZero,
### If vdPercentZero is not given, vdPercentZero will be generated by a relationship with Exp using this relationship between Mu and percent zero with an beta
dBetaGrandSD = c_dBetaGrandSD,
### The beta for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
fVerbose = FALSE
### Controls the plotting of graphic pdf logging (Default FALSE, TRUE indicates logging occurs)
){
  print("start funcGenerateFeatureParameters")

  if(fVerbose & !is.na(vdExp))
  {
    hist(vdExp,main="funcGenerateFeatureParameters: Original data expectation")
  }

  if(is.na(vdExp))
  {
    print("funcGenerateFeatureParameters: Generate vdExp Vector.")

    # Draw a vector of expectations
    # This will be the template distribution all bugs within a sample will be based on.
    # This allows to have structure within the data set of bugs more or less prevalent with a level of consistency
    # The mean of the distribution of bugs is derived from the max number of bug counts for a sample, making the
    # total number of bugs per sample the same within a small random margin
    lsParams = NULL
    if(c_dfFreezeSDGrandMu)
    {
      print("The Grand SD is frozen")
      lsParams = list(dBestExp=iReadDepth/int_number_features, dBestSD=1, dBestMu=iReadDepth/int_number_features, dBestDepth=iReadDepth, dTargetDepth=iReadDepth)
    } else {
      lsParams = funcGetParamsForReadDepth(dReadDepthPerFeature=iReadDepth/int_number_features, dBeta=dBetaGrandSD)
    }
    print("lsParams")
    print(lsParams)

    # Remember this is a vector of expectation not of log mu parameters.
    lExpVectorReturn = funcTruncatedRLNorm(iNumberMeasurements=int_number_features, dLogMean=lsParams$dLogMu, dLogSD=lsParams$dLogSD, iThreshold=exp((c_iTimesSDIsOutlier*lsParams$dLogSD)+lsParams$dLogMu))
    vdExp = lExpVectorReturn$Feature

    print("lExpVectorReturn")
    print(lExpVectorReturn)

    ### Update vdExp to read depth
    dCurReadDepthDifference = iReadDepth-sum(vdExp)

    ### Update  the distribution to the read depth requested.
    ### Depending on how many features are requested this is more or less needed
    ### This is not needed at the limit with many features.
    if(dCurReadDepthDifference > 0)
    {
      viSamplingIndices = sample(1:length(vdExp),abs(dCurReadDepthDifference),replace=TRUE)
      for(iIndex in viSamplingIndices)
      {
        vdExp[iIndex] = vdExp[iIndex]+1
      }
    } else if(dCurReadDepthDifference < 0) {
      for(iIndex in 1:round(abs(dCurReadDepthDifference)))
      {
        viGreaterThanOne = which(vdExp > 2)
        if(length( viGreaterThanOne ) > 0)
        {
          iGreaterThanOne = viGreaterThanOne
          if(length(iGreaterThanOne)>1){iGreaterThanOne = sample(iGreaterThanOne,1)}
          vdExp[iGreaterThanOne] = vdExp[iGreaterThanOne]-1
        }
      }
    }

    print("vdExp")
    print(vdExp)

    if(fVerbose)
    {
      barplot(vdExp, main=paste("Per feature count",mean(vdExp)), xlab="Feature", ylab="Expectation")
      abline(mean(vdExp),0,col="violet")
    }
    # Mus, SD, and percent zero depend on expectations. Since new expectations have been created these are reset to NA so they will be regenerated.
    vdMu = NA
    vdSD = NA
    vdPercentZero = NA
  }

  # Make sure there are no zeros in the expectation vector
  vdExp[ which( vdExp == 0 ) ] = c_dALittleMoreThanZero
  print("vdExp vdExp")
  print(vdExp)

  # TODO update
  if(is.na(vdSD))
  {
    # Generate vector of SD based on mu since it is not known
    print("funcGenerateFeatureParameters: Generate vdSD Vector.")

    vdSD = exp(funcEstimateFeatureSD(log(vdExp),dBetaSD))
    if(length(vdExp)==1){vdSD = 0}

    # Floor to close to 1 or exp(0)
    viLessThan1 = which(vdSD<1)
    if(length(viLessThan1)>0)
    {
      print(paste("funcGenerateFeatureParameters: Changing low SDs to 1. # occurences = ",length(viLessThan1)))
      vdSD[which(vdSD<1)] = 1
    }

    if(fVerbose)
    {
      vdVisExp = vdExp
      vdVisExp[which(vdVisExp==0)] = 0.00001
      vdVisSD = vdSD
      vdVisSD[which(vdVisSD==0)] = 0.00001
      plot(vdVisExp, vdVisSD, main="Generated Relationship of Exp and SD", col="orange")
    }
  }

  if(is.na(vdMu))
  {
    print("funcGenerateFeatureParameters: Generate vdMu Vector.")
    # We know the vdExp for each sample
    # We know the SD for each sample
    vdMu = sapply(1:length(vdExp), function(x) funcGetMu(vdExp[x],vdSD[x]))
  }

  if(is.na(vdPercentZero))
  {
    print("dBetaZero")
    print(dBetaZero)
    # Generate vector of percent zero based on exp since it is not known    print("funcGenerateFeatureParameters: Generate vdPercentZero Vector.")
    vdPercentZero = funcEstimatePercentZero(log(vdExp),dBetaZero)
    print("vdPercentZero")
    print(vdPercentZero)
    viLessThanZero = which(vdPercentZero < 0)
    if(length(viLessThanZero>0))
    {
      print(paste("funcGenerateFeatureParameters: Changing negative Percent Zeros to 0. # occurences = ",length(viLessThanZero)))
      vdPercentZero[vdPercentZero < 0] = 0
    }
    if(fVerbose){plot(log(vdExp), vdPercentZero, main="Generated Relationship of PercentZero and Log exp", col="purple")}
  }

  if(!length(vdMu)==int_number_features)
  {
    # This is the scenario that the calibration file is used and the number of the samples needed are not equal
    # This correct number of are selected with replacement.
    print("funcGenerateFeatureParameters: Reseting count of samples.")
    viWhich = sample(1:length(vdExp), size = int_number_features, replace = TRUE)
    vdExp = vdExp[viWhich]
    vdMu = vdMu[viWhich]
    vdSD = vdSD[viWhich]
    vdPercentZero = vdPercentZero[viWhich]   
  }

  print("stop funcGenerateFeatureParameters")

  # QC and contraints for percent zero
  # Make sure the percent zero passes the max
  # If there are not enough nonzeros, there is no signal to use.
  # This should be th max. Given there is a certain number of samples that have to have signal
  # The percentage of zeros must allow for those samples not to be zero and so restricts
  # the max the percent zero can be.
  dMaxPercent = 1-(iMinNumberSamples/int_number_samples)
  vdPercentZero[which(vdPercentZero>dMaxPercent)] = dMaxPercent
  return(list(exp=vdExp, mu=vdMu,sd=vdSD,PercentZero=vdPercentZero))
  ### exp Vector of expectations (the expectation vector)
  ### mu Vector of mu (not logMu) associated with the vdExp (by index)
  ### sd Vector of sd (not logSD) associated with the vdExp (by index)
  ### PercentZero Vector of percent zeros (0-1) associated with the vdExp (by index)
}


#!# TODO look over
func_generate_metadata = function(
int_base_metadata_number,
int_number_samples,
### Number of samples
dMinLevelPercent
### The minimum percent of samples a level can have
){
  # Preallocate matrix
  mat_metadata = matrix(data=NA,nrow=(int_base_metadata_number*c_iCountTypesOfMetadata),ncol=int_number_samples)

  # Used to report on metadata
  mtrxParameters = c(c_strMetadataDetails)

  # Continous metadata, generated means
  # Generating and padding the list of potential mean values
  li_mean_value_list = list(runif(1, c_dRunifMin, c_dRunifMax),1,100)
  if(length(li_mean_value_list) < (int_base_metadata_number*2))
  {
    for (k in (length(li_mean_value_list)+1):(int_base_metadata_number*2))
    {
      #TODO update the min and max
      li_mean_value_list = c(li_mean_value_list,runif(1, c_dRunifMin, c_dRunifMax))
    }
  } else {
    li_mean_value_list = li_mean_value_list[1:(int_base_metadata_number*2)]
  }
	
  # generating the continuous metadata
  # Should be random normal (is not bugs)
  i = 1
  for (mean_value in li_mean_value_list)
  {
    mat_metadata[i,] = rnorm(int_number_samples,mean=mean_value,sd=mean_value/5)
    mtrxParameters = c(mtrxParameters, paste(c_strMetadata, i, " ", c_strContinuous, sep =""))
    i = i+1
  }
	
  # generating the binary variables
  # Set up what the distributions for the different binary metadata will be
  li_probability_values_list = list(c(.5,.5))
  if(length(li_probability_values_list)+1 < int_base_metadata_number)
  {
    for (k in (length(li_probability_values_list)+1):int_base_metadata_number)
    {
      probability = runif(1, c_dMinBinary, c_dMaxBinary)
      li_probability_values_list = c(li_probability_values_list,list(c(probability, 1-probability)))
    }
  } else {
    li_probability_values_list = li_probability_values_list[1: int_base_metadata_number]
  }

  # Check to make sure minimal number of occurences is possible given the number of levels
  # for example, if the min percentage of occurences was originally 60%, it would not be possible to have 60% of samples
  # in both binary grouping, the max in this case would be 50%, the max for quarternery would be 25%
  # In this case the min value will be set to the max number of samples in a level in an even distribution * the percentage
  # So an original 60% samples for a binary case would be ceiling(60%*(#samples/2)*#samples) and a quarternery case 
  # would be ceiling(60%*(#samples/4)*#samples)
  iMinLevelBinaryCount = ceiling(int_number_samples*dMinLevelPercent)
  iMinLevelQuarterneryCount = iMinLevelBinaryCount
  if((dMinLevelPercent)>.25)
  {
    iMinLevelQuarterneryCount = ceiling(dMinLevelPercent*(int_number_samples/4)*int_number_samples)
    if((dMinLevelPercent)>.5)
    {
      iMinLevelBinaryCount = ceiling(dMinLevelPercent*(int_number_samples/2)*int_number_samples)
    }
  }

  # Draw from the binary values in binary_names given the previously generated binary distributions
  binary_names = c(1,2)
  for(probability in li_probability_values_list)
  {
    vsCurMetadata = sample(binary_names, size=int_number_samples, prob=probability, replace = TRUE)
    fMetadataFailed = !funcIsFactorMetadataValid(vsCurMetadata, iMinLevelBinaryCount)
    iBinaryMetadataLoop = 1
    while(fMetadataFailed)
    {
      vsCurMetadata = sample(binary_names, size=int_number_samples, prob=probability, replace = TRUE)
      fMetadataFailed = !funcIsFactorMetadataValid(vsCurMetadata, iMinLevelBinaryCount)
      iBinaryMetadataLoop = iBinaryMetadataLoop + 1
      if(iBinaryMetadataLoop > iLoopingControlIncrement)
      {
        fMetadataFailed = FALSE
        print(paste("Suboptional metadata was created, did not pass quality control, is too imbalanced. Minimum level is preferred to be ", iMinLevelBinaryCount,"."))
      }
    }
    mat_metadata[i,] = vsCurMetadata
    mtrxParameters = c(mtrxParameters, paste(c_strMetadata, i, " ", c_strFactor, " ", paste(levels(as.factor(mat_metadata[i,])), collapse = " "), sep=""))
    i = i+1
  }

  # generating the quarternary metadata with simple distributions and at least one uniform
  # names of the feature choices
  li_list_of_distributions = list(uniform=c(.25,.25,.25,.25), c(.2,.2,.3,.3), c(.3,.3,.2,.2), c(.2,.2,.2,.4))
  if(length(li_list_of_distributions) < int_base_metadata_number)
  {
    for (k in (length(li_list_of_distributions)+1):int_base_metadata_number)
    {
      li_list_of_distributions = c(li_list_of_distributions, list(c(.25,.25,.25,.25)))
    }
  } else {
    li_list_of_distributions = li_list_of_distributions[1:int_base_metadata_number]
  }
		
  # handpicked n-nomial distributions (metadata values)
  quarternary_names = c(1,2,3,4)

  for (distribution in li_list_of_distributions)
  {
    vsCurMetadata = sample(quarternary_names, size=int_number_samples, prob=distribution, replace = TRUE)
    fMetadataFailed = !funcIsFactorMetadataValid(vsCurMetadata, iMinLevelQuarterneryCount)
    iBinaryMetadataLoop = 1
    while(fMetadataFailed)
    {
      vsCurMetadata = sample(quarternary_names, size=int_number_samples, prob=distribution, replace = TRUE)
      fMetadataFailed = !funcIsFactorMetadataValid(vsCurMetadata, iMinLevelQuarterneryCount)
      iBinaryMetadataLoop = iBinaryMetadataLoop + 1
      if(iBinaryMetadataLoop > iLoopingControlIncrement)
      {
        fMetadataFailed = FALSE
        #print(paste("Suboptional metadata was created, did not pass quality control, is too imbalanced. Minimum level is preferred to be ", iMinLevelQuarterneryCount,"."))
        #print(vcCurMetadata)
      }
    }
    mat_metadata[i,] = vsCurMetadata

    mtrxParameters = c(mtrxParameters, paste(c_strMetadata, i, " ", c_strFactor, " ", paste(levels(as.factor(mat_metadata[i,])), collapse = " "), sep=""))
    i = i+1
  }
  return(list(mat_metadata=mat_metadata, mtrxParameters=mtrxParameters))
}

#!# TODO check for exp update
func_generate_random_lognormal_matrix = function(
int_number_features,
### Number of features
int_number_samples,
### Number of samples,
iMinNumberCounts,
### Minimum number of counts for a feature to be considered in a sample for the QC filtering
iMinNumberSamples,
### Created bugs must have a minimum number of samples (iMinNumberSamples) that have a minimum number of counts (iMinNumberSamples)
iReadDepth,
### Simulated read depth for sample creation
vdExp = NA,
### The vector of expectations for each feature. If not provided one will be generated and vdMu, vdSD, and xPercentZero will be reset
vdMu = NA,
### Vector of log Mu parameters for the original mu distribution (means of features) if not supplied, one will be generated
xPercentZero = c_dBetaZero,
### Either a vector of percent zero for each mu or a single parameter modeling the relationship between percent zeros and mus
vdSD = NA,
fZeroInflate = TRUE,
### Turns off Zero inflation if FALSE (default TRUE, zero inflation turned on)
dBetaSD = c_dBetaSD,
### The beta in the relationship between the sd and the feature exp
dBetaZero = c_dZeroBeta,
### The beta in the relationship between the percent zero and the exp
dBetaGrandSD = c_dBetaGrandSD,
### The beta in the relationship between mu and sd in the vector of feature expectations
fVerbose = FALSE
### Turns on logging (typically generates pdfs)
){
  print(paste("func_generate_random_lognormal_matrix::", "int_number_samples", int_number_samples, "int_number_features", int_number_features, "iReadDepth", iReadDepth, "dBetaSD", dBetaSD, "fZeroInflate", fZeroInflate))

  # Preallocating for speed
  mat_bugs = matrix(data=NA,nrow=int_number_features,ncol=int_number_samples)
  # Get the initial mu vector for generating features.
  lsInitialDistribution = funcGenerateFeatureParameters(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberSamples=iMinNumberSamples, iReadDepth=iReadDepth, vdExp=vdExp, vdMu=vdMu, vdSD=vdSD, vdPercentZero=xPercentZero, dBetaSD=dBetaSD, dBetaZero=dBetaZero, dBetaGrandSD = dBetaGrandSD, fVerbose=fVerbose)

  print("lsInitialDistribution")
  print(lsInitialDistribution)
  print("SUM exp")
  print(sum(lsInitialDistribution$exp))

#  # Need to set the distributions to the right magnitude, previously in the log form
#  lsInitialDistribution$mu = exp(1)^lsInitialDistribution$mu
#  lsInitialDistribution$sd = exp(1)^lsInitialDistribution$sd

  # Update the Mu, SD and Percent zero bugs and report on distributions
  vdMu = lsInitialDistribution[["mu"]]
  vdSD = lsInitialDistribution[["sd"]]
  xPercentZero = lsInitialDistribution[["PercentZero"]]
  vdExp = lsInitialDistribution[["exp"]]

  # TODO remove
  if(c_dfFreezeSDFeatures)
  { 
    print("Feature SDs are frozen.")
    if(length(vdSD)>0)
    {
      vdMu = rep(iReadDepth/int_number_features, length(vdSD))
      vdSD = rep(1,length(vdSD))
    }
  }

  # Number of samples needed to have signal as a constraint
  iNumberSamples = min(int_number_samples, iMinNumberSamples)

  # Record the current read depths of the columns (samples)
  vdLeftOver = rep(0, int_number_samples)
  print("*** 3 fZeroInflate")
  print(fZeroInflate)
  # Make features and assign feature samples to samples giving higher counts to lower read depth samples.
  for(iReset in 1:int_number_features)
  {
    # Create new feature
    lFeatureDetails = funcMakeFeature(dMu=vdMu[iReset], dSD=vdSD[iReset], dPercentZero=xPercentZero[iReset], iNumberSamples=int_number_samples, iMinNumberCounts=iMinNumberCounts, iMinNumberSamples=iMinNumberSamples, dTruncateThreshold=(c_iTimesSDIsOutlier*vdSD[iReset])+vdMu[iReset], fZeroInflate=fZeroInflate, fVerbose=fVerbose )

    # Store the left over counts per sample
    vdLeftOver = vdLeftOver + lFeatureDetails$LeftOver

    # Update the matrix with the new feature
    mat_bugs[iReset,] = lFeatureDetails$Feature
  }

  # Remove any fully zero sample
  lZeroCorrectionResults = funcZeroCorrectMatrix(mtrxData=mat_bugs, vdFeatureMus=vdExp, vdLeftOver=vdLeftOver)
  mat_bugs = lZeroCorrectionResults[["Data"]]
  vdLeftOver = lZeroCorrectionResults[["LeftOver"]]
  
  # Shuffle back in removed signal.
  mat_bugs = funcShuffleMatrix(mtrxData=mat_bugs, vdFeatureMus=vdMu, vdShuffleIn=vdLeftOver)

  # Floor to Counts
  lRoundedResult = funcRoundMatrix(mtrxData=mat_bugs, vdLeftOver=vdLeftOver)
  mat_bugs = lRoundedResult$Data
  vdLeftOver = lRoundedResult$LeftOver

  if(c_dfFreezeSDGrandMu ||c_dfFreezeSDFeatures||c_fPrintLognormalMatrix)
  {
#    print("bug counts")
#    print(mat_bugs)
    print("Read Depth")
    print(colSums(mat_bugs))
    print("Average read depth")
    print(mean(colSums(mat_bugs)))
    print("Feature mean")
    print(funcGetRowMetric(mat_bugs,mean))
    print("vdExp")
    print(vdExp)
    print("Sum exp")
    print(sum(vdExp))
    print("Sum exp, should be read depth")
    print(sum(vdExp))
    print("Feature mean summary")
    print(summary(funcGetRowMetric(mat_bugs,mean)))
  }

  if(fVerbose)
  {
    ## Plot
    plot(colSums(mat_bugs),main=paste("Read depth mean=",mean(colSums)), xlab="Samples")
    abline(mean(colSums(mat_bugs)),0)
  }

  # Truth table for log normal data
  mtrxParameters = matrix(data=NA, nrow=6, ncol=1)
  mtrxParameters[1,1] = paste(c_strSyntheticMicrobiome, c_strRandom, sep='')
  mtrxParameters[2,1] = paste(c_strNumberOfFeatures, int_number_features)
  mtrxParameters[3,1] = paste(c_strNumberOfSamples, int_number_samples)
  mtrxParameters[4,1] = paste(c_strTotalSampleBugOccurrence, iReadDepth)
  mtrxParameters[5,1] = paste(c_strNumberCounts, iMinNumberCounts)
  mtrxParameters[6,1] = paste(c_strNumberSamples, iMinNumberSamples)

  print("stop func_generate_random_lognormal_matrix")
  return(list(mat_bugs=mat_bugs, mtrxParameters=mtrxParameters))
  ### Returns a row major matrix of log-normal data.
}


func_generate_random_lognormal_with_outliers = function(
### Generates a random log normal distribution of data as a null matrix
### The option of using a matrix passed in as a parameter as the null matrix is provided (mtrxBugs)
### A percent of samples are given outliers based on the percent parameter
int_number_features,
### Number of features
int_number_samples,
### Number of samples
iMinNumberCounts,
### Minimum number of counts for a feature to be considered in a sample for the QC filtering
iMinNumberSamples,
### Created bugs must have a minimum number of samples (iMinNumberSamples) that have a minimum number of counts (iMinNumberSamples)
dMaxPercentOutliers,
### The maximum percent of outliers to create in each sample (0 =< dPercent =< 1.0)
dPercentSamples,
### Percent of samples given outliers (0 =< dPercent =< 1.0)
mtrxBugs,
### Precalculated null matrix
fVerbose = FALSE
){
  print("start func_generate_random_lognormal_with_outliers")
  # Plot Null matrix by sample and feature before outliers are generated
  if(fVerbose)
  {
    funcBoxPlotOutliers(mtrxData=mtrxBugs, strTitle="Random Log Normal - No Outliers")
    funcBoxPlotOutliers(mtrxData=mtrxBugs, strTitle="Random Log Normal - No Outliers", fWithZeros=FALSE)
    funcBoxPlotOutliers(mtrxData=mtrxBugs, strTitle="Random Log Normal - No Outliers", fBySamples=FALSE)
    funcBoxPlotOutliers(mtrxData=mtrxBugs, strTitle="Random Log Normal - No Outliers", fWithZeros=FALSE, fBySamples=FALSE)
  }

  # convert to dataframe
  mat_bugs_dataframe = as.data.frame(mtrxBugs)

  # Determine the number of samples
  iNumberSamples = round(int_number_samples*dPercentSamples)
  # Randomly select samples
  liIndices = sample(1:int_number_samples,iNumberSamples,replace=FALSE)

  # Detemine the max number of outliers
  iMaxNumberOfOutliers = round(int_number_features*dMaxPercentOutliers)

  lviSwapped = list()
  lviFeatureSwaps = list()
  if(iMaxNumberOfOutliers>0)
  {
    #for each column
    for(iSample in liIndices)
    {
      # Select the number of outliers for the sample
      iNumberOutliers = sample(1:iMaxNumberOfOutliers,1,replace=FALSE)

      # sort the column
      dfSorted = mat_bugs_dataframe[order(mat_bugs_dataframe[,iSample]),]

      viAlreadySelected = c()
      iBufferForZeros = 0
      if(c_fIgnoreZerosInOutliers)
      {
        iBufferForZeros = length(which(dfSorted[[iSample]]==0))
      }
      if((iNumberOutliers+iBufferForZeros)>length(dfSorted[[iSample]]))
      {
        iBufferForZeros = length(dfSorted[[iSample]])-iNumberOutliers
        print(paste("Allowing",iBufferForZeros,"zeros to be selected in outlier swapping so that all",iNumberOutliers,"swaps can be performed in sample",iSample))
      }

      for (iOutlier in 1:iNumberOutliers)
      {
        # Get the indices of all Min and max values ties
        iMaxValue = dfSorted[(int_number_features-(iOutlier-1)), iSample]
        iMinValue = dfSorted[iOutlier+iBufferForZeros, iSample]
        viMaxValueIndices = which(mat_bugs_dataframe[[iSample]]==iMaxValue)
        viMinValueIndices = which(mat_bugs_dataframe[[iSample]]==iMinValue)
        iRowMax = setdiff(viMaxValueIndices,viAlreadySelected)
        iRowMin = setdiff(viMinValueIndices,viAlreadySelected)

        # Sample from ties
        if(length(iRowMax) > 1)
        {
          iRowMax = sample(c(setdiff(viMaxValueIndices,viAlreadySelected)), size = 1)
        }
        if(length(iRowMin) > 1)
        {
          iRowMin = sample(c(setdiff(viMinValueIndices,viAlreadySelected)), size = 1)
        }

        mtrxBugs[c(iRowMin, iRowMax), iSample] = mtrxBugs[c(iRowMax, iRowMin), iSample]
        viAlreadySelected = c(viAlreadySelected, iRowMax, iRowMin)
      }
      # Record which features were swapped with lists as features, basically a transpose
      for(iSwapIndex in viAlreadySelected)
      {
        lviFeatureSwaps[[iSwapIndex]] = c(lviFeatureSwaps[iSwapIndex][[1]],iSample)
      }
      lviSwapped[[iSample]] = viAlreadySelected
    }
  }

  if(fVerbose)
  {
    # Plot box plot with outliers
    funcBoxPlotOutliers(mtrxData=mtrxBugs, strTitle = "Random Log Normal - With Outliers", lviHighlight=lviSwapped)
    funcBoxPlotOutliers(mtrxData=mtrxBugs, strTitle = "Random Log Normal - With Outliers", fWithZeros=FALSE, lviHighlight=lviSwapped)
    funcBoxPlotOutliers(mtrxData=mtrxBugs, strTitle = "Random Log Normal - With Outliers", strColor = c_strFeatureOutlier, fBySamples=FALSE, lviHighlight=lviFeatureSwaps)
    funcBoxPlotOutliers(mtrxData=mtrxBugs, strTitle = "Random Log Normal - With Outliers", strColor = c_strFeatureOutlier, fWithZeros=FALSE, fBySamples=FALSE, lviHighlight=lviFeatureSwaps)
  }

  # Count swaps
  iSwapCount = c()
  if(length(lviSwapped)>0)
  {
    for(iIndex in 1:length(lviSwapped))
    {
      iSwapCount=c(iSwapCount, lviSwapped[[iIndex]])
    }
  }
  iSwapCount = length(iSwapCount)

  # Truth table for log normal data
  mtrxParameters = c(paste(c_strSyntheticMicrobiome, c_strOutlier, sep=''))
  mtrxParameters = c(mtrxParameters, paste(c_strNumberOfFeatures, int_number_features))
  mtrxParameters = c(mtrxParameters, paste(c_strNumberOfSamples, int_number_samples))
  mtrxParameters = c(mtrxParameters, paste(c_strPercentOutliers, dMaxPercentOutliers))
  mtrxParameters = c(mtrxParameters, paste(c_strPercentSampleOutliers, dPercentSamples))
  if(length(lviSwapped)>0)
  {
    for(iIndex in 1:length(lviSwapped))
    {
      for(iItemSwapped in lviSwapped[[iIndex]])
      {
        mtrxParameters = c(mtrxParameters, paste(c_strOutlierParameter,paste(c_strFeature,c_strOutlier, iItemSwapped, sep="_"),c_strSampleParameter,iIndex))
      }
    }
  }

  print("Stop func_generate_random_lognormal_with_outliers")
  # And return
  return(list(mat_bugs=mtrxBugs, mtrxParameters=mtrxParameters))
}


func_generate_random_lognormal_with_multivariate_spikes = function(
### Spike in associations with 1 or more metadata
int_number_features,
### Number of features
int_number_samples,
### Number of samples
iMinNumberCounts,
### Minimum number of counts for a feature to be considered in a sample for the QC filtering
iMinNumberSamples,
### Created bugs must have a minimum number of samples (iMinNumberSamples ) that have a minimum number of counts (iMinNumberSamples)
percent_spikes,
### The percent of features (bugs) to spike
multiplier,
### Used to multiple the metadata before adding to a feature to strengthen the signal of the metadata
metadata_matrix,
### Matrix of metadata (by row) to spike in
multivariate_parameter,
### Number of metadata to spike in
dMinLevelCountPercent,
### Minimum number of samples allowed to be part of the spiked in relationship
mtrxBugs,
### Random log normal matrix
fZeroInflated,
### True indicates it is a zero inflated model.
lviFrozeMetadataIndices = NULL,
### If given, the method must select these specific metadata indicies.
### This allow selection of features to be the same when evaluating the multiplier
liFrozeDataIndicies = NULL,
### If given, the method must select these specific data indicies.
### This allows the selection of features to be the same when evaluating the multiplier
lsFrozeLevels = NULL,
### If given, the method must select these specific level indicies.
### This allows the slection of features to be the same when evaluating the multiplier
fVerbose = FALSE
){
  print("start func_generate_random_lognormal_with_multivariate_spikes")

  # Tracks the bug of interest
  iIndexSpikedFeature = NA

  # Initialize froze levels if need be
  if(is.null(lsFrozeLevels))
  {
    lsFrozeLevels = list()
  }

  # looping control
  iLoopingControl = min(iLoopingControlIncrement, int_number_samples*int_number_features)

  # Min number of samples in the relationship
  iMinSamples = floor(dMinLevelCountPercent*int_number_samples)

  # Creating the truth file which contains true positive spike-ins
  strMult = paste(multiplier)
  strMult = sub(".","_", strMult, fixed=TRUE)
  m_parameter_rec = c(paste(paste(c_strSyntheticMicrobiome, c_strSpike,sep=""), "n", multivariate_parameter,"m", strMult, sep='_'))
  m_parameter_rec = c(m_parameter_rec, paste(c_strNumberOfFeatures, int_number_features))
  m_parameter_rec = c(m_parameter_rec, paste(c_strNumberOfSamples, int_number_samples))
  m_parameter_rec = c(m_parameter_rec, paste(c_strPercentSpikes, percent_spikes))
  m_parameter_rec = c(m_parameter_rec, paste(c_strMultiplier, multiplier))
  m_parameter_rec = c(m_parameter_rec, paste(c_strMultiplierParameter, multivariate_parameter))
  m_parameter_rec = c(m_parameter_rec, paste(c_strMinimumSamples, iMinSamples))

  # Features to be spiked to select from
  viRemainingFeatures = 1:int_number_features

  # Holds the metadata selected for spikin
  lviMetadata = list()
  # holds the data selected for spikin
  liData = list()

  # The number of bugs to spike in
  iSpikeinCount = floor(int_number_features*percent_spikes)

  if(iSpikeinCount > 0)
  {
    # For each spiked in bug
    for(iSpikedBug in 1:iSpikeinCount)
    {
      vsCurFrozeLevels = c()
      iSpikeInLoopControl = 1
      # Holds the currently selected metadata indicies
      viSelectedMetadata = NULL
      # Currently selected metadata names
      vstrSpikedMetadata = NA
      # Bug to spike with
      vdCurData = NA
      # Indicator if the spikin passed quality control
      fSpikeInFailed = TRUE
      # Current best failed run
      # List including 
      ## Metadata = matrix of metadata that has been prepped for the spikin (row major),
      ## MetadataNames = vector of strings for the name sof each row of metadata,
      ## SpikinBug = vector of doubles (bug measurements),
      ## Count = vector of non-zeros elements overlapping both the bugs and each metadata,
      ## if not using metadata at a level but using all levels, the overlap is given for each level.
      lxCurrentBestRun = list(Metadata = c(), MetadataNames = c(), SpikinBug = c(), BugIndex = NA, Count = -1)

      # Metadata indices that are selected for spikin
      viSelectedMetadata = c()

      # Find valid spike-in scenario or best choice
      while(fSpikeInFailed)
      {
        # Get the bug to attempt the association
        # If previously associations have been made
        # liFrozenDataIndices makes the same associations happen here
        # This is so if multiple multipliers are given
        # The different matrices show the differences given increased size of effect
        # Not difference driven by selecting diffrent bugs
        if(!is.null(liFrozeDataIndicies) & length(lviFrozeMetadataIndices)>0)
        {
          iIndexSpikedFeature = liFrozeDataIndicies[[iSpikedBug]]
          iSpikeInLoopControl = iSpikeinCount + 1
        } else {
          iIndexSpikedFeature = sample(viRemainingFeatures,1)
        }

        # Select which of the metadatum we will be using to scale
        if(!is.null(lviFrozeMetadataIndices) & length(lviFrozeMetadataIndices)>0)
        {
          viSelectedMetadata = lviFrozeMetadataIndices[[iSpikedBug]]
        } else {
          viSelectedMetadata = sample(1:nrow(metadata_matrix), multivariate_parameter, replace=TRUE)
        }

        # Get matrix of metadata and feature
        vdCurMetadata = metadata_matrix[viSelectedMetadata,]
        vdCurData = mtrxBugs[iIndexSpikedFeature,]

        if(fVerbose)
        {
          plot(vdCurData,main=paste("Spikin: Original bug ", paste(viSelectedMetadata,sep=","),sep=""))
        }

        # Attempt to spike in a new bug
        lsMetadataInfo = funcPrepareMetadata(viSelectedMetadata=viSelectedMetadata, vdCurMetadata=vdCurMetadata, vsFrozeLevels=unlist(lsFrozeLevels[iSpikedBug]))
        vdCurMetadata = lsMetadataInfo[["metadata"]]
        vstrSpikedMetadata = lsMetadataInfo[["names"]]
        vsCurFrozeLevels = lsMetadataInfo[["vsLevels"]]

        # Spike in new bug
        vdSpikedBug = funcSpikeNewBug(vdCurMetadata=vdCurMetadata, vdCurData=vdCurData, multiplier=multiplier, fZeroInflated=fZeroInflated)

        # Check to see if a failure occured
        lxQCInfo = funcQCSpikin(vdCurMetadata, vdSpikedBug, iMinSamples)

        # lxQCInfo has the slots PASS = Boolean, CommonCounts = vector of integers
        fSpikeInFailed = !lxQCInfo[["PASS"]]

        # Check to see if the looping must end and no success is given
        # Also if the spikein failed, make sure to update the best failed scenario
        iSpikeInLoopControl = iSpikeInLoopControl + 1
        if(fSpikeInFailed)
        {
          # Keep the best failed scenario so far.
          if(!is.null(lxQCInfo[["CommonCounts"]]) & !is.null(lxQCInfo[["SpikinBug"]]) & sum(lxQCInfo[["CommonCounts"]]>iMinSamples) > lxCurrentBestRun[["Count"]])
          {
            lxCurrentBestRun = list(Metadata = vdCurMetadata, MetadataNames = vstrSpikedMetadata, SpikinBug = vdSpikedBug, BugIndex = iIndexSpikedFeature,  Count = sum(lxQCInfo[["CommonCounts"]]>iMinSamples), MetadataIndices = viSelectedMetadata, Levels = vsCurFrozeLevels)
          }
        
          # If we have ran out of iteration, use the best failed scenario and indicate this.
          if(iSpikeInLoopControl > iLoopingControl)
          {
            # Reset the current spikin variables to the best scenario
            vdCurMetadata = lxCurrentBestRun[["Metadata"]]
            vstrSpikedMetadata = lxCurrentBestRun[["MetadataNames"]]
            vdSpikedBug = lxCurrentBestRun[["SpikinBug"]]
            iIndexSpikedFeature = lxCurrentBestRun[["BugIndex"]]
            viSelectedMetadata = lxCurrentBestRun[["MetadataIndices"]]
            if(length(lsFrozeLevels) < iSpikedBug)
            {
              lsFrozeLevels[iSpikedBug] = lxCurrentBestRun[["Levels"]]
            }
            print(paste("While spiking in a relationship between metadata and bug was not able to meet the minimal percentage of spiked-in samples for the relationship. Min count = ", iMinSamples))
            print(paste("The following spike-in does not pass quality control: Bug=", iIndexSpikedFeature," Metadata=",paste(vstrSpikedMetadata,collapse=","),"Multiplier=", multiplier,"Count=", multivariate_parameter))
            # Break while to use the best case scenario
            break
          }
        } else {
          lxCurrentBestRun = list(Metadata = vdCurMetadata, MetadataNames = vstrSpikedMetadata, SpikinBug = vdSpikedBug, BugIndex = iIndexSpikedFeature,  Count = sum(lxQCInfo[["CommonCounts"]]>iMinSamples), MetadataIndices = viSelectedMetadata, Levels = vsCurFrozeLevels)
        }
      }

      # If the spike-in was successful, update
      if(!is.na(lxCurrentBestRun[["BugIndex"]]))
      {
        # Update metadata and data indices for the features used for spikins
        lviMetadata[[ length(lviMetadata)+1 ]] = viSelectedMetadata
        liData[[ length(liData)+1 ]] = iIndexSpikedFeature

        # Update spike-in
        mtrxBugs[iIndexSpikedFeature,] = vdSpikedBug
        if(length(lsFrozeLevels) < iSpikedBug)
        {
          lsFrozeLevels[[iSpikedBug]] = vsCurFrozeLevels
        }

        # Update the truth table with which feature is spiked with which metadata
        m_parameter_rec = c(m_parameter_rec, paste(c_strFeature,c_strSpike,"n", multivariate_parameter,"m", strMult,iIndexSpikedFeature,sep='_'))
        m_parameter_rec = c(m_parameter_rec, vstrSpikedMetadata)

        # Remove bug from pool of potential bugs to spike.
        viRemainingFeatures = setdiff(viRemainingFeatures, iIndexSpikedFeature)

        # Start making plots
        if(fVerbose)
        {
          plot(0, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
          text(1,1,paste("Start Association Bug", iIndexSpikedFeature,"with metadata",paste(viSelectedMetadata,sep=",")))
        }
        i = 1
        for(iSelectedMetadata in viSelectedMetadata)
        {
          # Plot spike in relationship by
          # Non-zero bug entries and metadata dummy level if there is one.
          vstrDummyLevel = rep(c_strOutlierColor,length(vdCurMetadata[i,]))
          vstrDummyLevel[which(vdCurData==0)] = c_strDefaultMarkerColor
          # Color original relationship by bug non zero entries
          vstrRelationship = rep(c_strDefaultMarkerColor,length(vdCurData))
          vstrRelationship[which(vdCurData!=0)] = c_strOutlierColor

          if(!sum(sapply(vdCurMetadata[i,],funcNumericIsInt))==0)
          {
            vstrDummyLevel[which(vdCurMetadata[i,]==1)] = c_strDefaultMarkerColor
            vstrRelationship[which(vdCurMetadata[i,]==1)] = c_strDefaultMarkerColor
          }
          if(fVerbose)
          {
            plot(vdCurData, col=vstrDummyLevel, main=paste("Spikin: Original bug ", iIndexSpikedFeature,sep=""))
            plot(as.vector(vdSpikedBug), col=vstrDummyLevel, main=paste("Spikin: New bug ", iIndexSpikedFeature,sep=""))
            plot(vdCurMetadata[i,], col=vstrDummyLevel, main="Spikin: Dummied metadata")
            plot(metadata_matrix[iSelectedMetadata,], col=vstrDummyLevel, main=paste("Spikin: Original metadata n=", iSelectedMetadata, " of ", multivariate_parameter,sep=""))
            plot(metadata_matrix[iSelectedMetadata,], vdCurData, col=vstrRelationship, main = paste("Spikin: Original Relationship", iSelectedMetadata," in ", paste(viSelectedMetadata,sep=","),sep=""))
            plot(metadata_matrix[iSelectedMetadata,], vdSpikedBug, col=vstrRelationship, main = paste("Spikin: Spikin Relationship", iSelectedMetadata," in ", paste(viSelectedMetadata,sep=","),sep=""))
          }
          i = i + 1
        }
        if(fVerbose)
        {
          plot(0, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
          text(1,1,paste("Stop Association Bug", iIndexSpikedFeature,"with metadata",paste(viSelectedMetadata,sep=",")))
        }
      }
    }
  }
  # Floor to counts (spike-ins will be real numbers)
  mtrxBugs = floor(mtrxBugs)

  print("stop func_generate_random_lognormal_with_multivariate_spikes")

  # Return normalized matrix and spike-in list
  return(list(mat_bugs=mtrxBugs, m_parameter_rec=m_parameter_rec, MetadataIndices=lviMetadata, DataIndices=liData, Levels=lsFrozeLevels))
}


# 3 Tests 9/4/2013
funcGetExp = function(
### Given the Mu and SD, get the expectation
### Will need to 
# Expects the values not logged
# Log is base exp(1)
# Returns a value that target the mean(data) not the logged data.
dMu,
dSD
){
  return( exp( log( dMu, exp( 1 ) ) + 0.5 * ( log( dSD, exp( 1 ) )^2 ) ) )
}


# 3 Tests 9/4/2013
funcGetMu = function(
### Get the mu given SD and an expectation
# Expects the Exp to be the mean(data) and SD to be the sd(log(data))
# so sd(log(rlnorm(10000,3,2))) to calculate the value
# output from these functions can be directly used in the rlnorm (is the logged mu)
# Log is base exp(1)
# To use these values in rlnorm you would log them
dEx,
dSD
){
  return( exp( -1 * ( ( ( log( dSD, exp( 1 ) )^2 ) /2 ) - log( dEx, exp( 1 ) ) ) ) )
}


# 2 Tests 9/4/2013
funcGetParamsForReadDepth = function(
### From the relationship between the grand mu and grand sd contained in dBetadGradSD
### A mu and sd is calculated that should satisfy the given read depth
dReadDepthPerFeature,
### The read depth of interest in the original scale
dBeta
### The beta describing the relationship
){
  # Get the associated SD from the EXP
  dLogSD = funcEstimateGrandSD( dReadDepthPerFeature, dBeta )

  # Get the associated mu
  dMu = funcGetMu( dReadDepthPerFeature, exp( dLogSD ) )

  return( list( dLogMu = log( dMu ), dLogSD = dLogSD ) )
  ### dLogMu: The logMu which, with the logSD, would give the read depth (exp)
  ### dLogSD: The logMu which, with the logSD, would give the read depth (exp)
}


# 3 Tests 9/4/2013
funcGetRowMetric = function(
### Get mean for vector or matrix
lxValues,
### Values of metadata which may be discrete (in which after the _ is a number) or a numeric
### Method to perform on values
funcMethod
){
  if( is.null( dim( lxValues )[1] ) )
  {
    return( mean( lxValues, na.rm=TRUE ) )
  }
  return( apply( lxValues, 1, funcMethod, na.rm=TRUE ) )
}


# 3 Tests 9/4/2013
funcGetSD = function(
### Get standard deviation for vector or matrix
lxValues
### matrix or vector
){
  if( is.null( dim( lxValues )[ 1 ] ) )
  {
    return( sd( lxValues, na.rm = TRUE ) )
  }
  return( apply( lxValues, 1, sd, na.rm = TRUE ) )
}


# 4 Tests 9/4/2013
funcIsFactorMetadataValid = function(
### Check to make sure a level of entrophy is in a discontinuous metadata so it has a minimal level to associate
vxMetadata,
### A metadata that can be factor data
iMin
### Minimum number of instances of a level
){
  vxFactorMetadata = as.factor( vxMetadata )
  lstrLevels = levels( vxFactorMetadata )
  for( strLevel in lstrLevels )
  {
    if( length( which( vxFactorMetadata == strLevel ) ) < iMin )
    {
      return( FALSE )
    }
  }
  return( TRUE )
}


funcMakeFeature = function(
### Create a feature given parameters
### Uses a zero inflated model if necessary
### Enforces a min number samples with signal if needed.
dMu,
### Mu of the rnorm distribution associate with the rlnorm draw that will occur, this will be logged
dSD,
### SD of the rnorm distribution associate with the rlnorm draw that will occur, this will be logged
dPercentZero,
### Percent of zeros to add to the feature
iNumberSamples,
### Number of measurements for the feature
iMinNumberCounts,
### Minimum number of counts to be considered signal
iMinNumberSamples,
### Minimum number of samples needed to have signal. If this is not fulfilled, signal will be generated and added.
dTruncateThreshold = NA,
### Threshold to truncate the underlying distribution
fZeroInflate = TRUE,
### If the feature should be zero inflated
fVerbose = FALSE
### If pdf logging of feaure should occur
){
#  print("funcMakeFeature")
  # If not zero inflated
  if(!fZeroInflate){dPercentZero = 0}

  # Generate feature
  lFeatureData = func_zero_inflate(log(dMu), dPercentZero, iNumberSamples, log(dSD), dTruncateThreshold)
  dFeature = lFeatureData$Feature

  # TODO I would rather have this smoother
  # Min is placed here because if you ask for 2 samples but this is evaluating to find atleast 3 samples at a certain level you will inf loop
  # Look to see how many measurements are above the min.
  # If there are not enough counts more the the min number of counts
  # then make sure there are by drawing more samples (not zero inflated)
  # And changing any measurement less than the min to the min
  # Then randomly assigning this to a location that is zero
  vdSignal = which(dFeature >= iMinNumberCounts)

  iNeededExtraValues = iMinNumberSamples - length(vdSignal)

  if(iNeededExtraValues > 0)
  {
    viNonZeroLocations = which(dFeature<iMinNumberCounts)
    if(length(viNonZeroLocations) < iNeededExtraValues)
    {
      viNonZeroLocations = which(dFeature!=0)
      if(length(viNonZeroLocations)<iNeededExtraValues)
      {
        viNonZeroLocations = 1:length(vdSignal)
      }
    }
    # TODO need to assign using sample not feature
    # Randomly select locations to add additional signal to
    # Would like this to be non-zero value that are less than the min but if they do not exist may just be any value.
    vdNewSignalLocations = sample(x=viNonZeroLocations, iNeededExtraValues, replace=FALSE)

    # Replace the feature with these nonzero measurements
    # Adding in the random element to smooth out the min value when added.
    dFeature[vdNewSignalLocations] = dFeature [vdNewSignalLocations] + sapply(1:length(vdNewSignalLocations), function(x) iMinNumberCounts * 1 + runif(1))
  }

  # Extra useful measurements, the true and expected means
  dMean = mean(dFeature)
  dExpCal = funcGetExp(dMu,dSD)
  dExp = mean(dMu,dSD)

  if(fVerbose)
  {
    plot( dFeature, main = paste("funcMakeFeature","Mean",round(dMean,2),"SD",round(sd(log(dFeature)),2),"Exp",round(dExp,2)))
  }
#  print("End funcMakeFeature dFeature and leftover")
  return(list(Feature=dFeature, Mean=dMean, Exp=dExp, ExpCal=dExpCal, LeftOver=lFeatureData$LeftOver))
}


# 4 Tests 9/4/2013
funcNormalizeMicrobiome = function(
### Normalize the data assuming that the 
### Samples are columns and the features (bugs) are rows.
### Zero columns are preserved as zero.
mtrxMicrobiome
### Matrix or dataframe microbiome
){
  mtrxRelAb = mtrxMicrobiome / unlist( lapply( colSums( mtrxMicrobiome ), function( x ) rep( x, nrow( mtrxMicrobiome ) ) ) )
  viZero = which( apply( mtrxMicrobiome, 2, sum ) == 0 )
  mtrxRelAb[, viZero] = rep( 0, nrow( mtrxRelAb ) )
  return( mtrxRelAb )
}


# 4 Tests 9/4/2013
# Will fail if all are integer if using to detect metadata
# Returns a TRUE (1) for each int value so if
# the return vector's sum is 0, the data is not numeric
funcNumericIsInt = function(
### Tests to see if a numeric is an interger
dValue
### The numeric value to check
){
  return( paste( as.integer( dValue ), sep = "" ) == paste( dValue, sep = "" ) )
}


# Plotting function. Not in regression suite
funcPlotSpikeins = function(
### Plots normalized spikeins
vsTruth,
mtrxFinal,
strFileDir,
strFileName
)
{
  # Print the spiked-in markers after normalization to make sure they are still there
  pdf(file.path(strFileDir,strFileName), useDingbats=FALSE)
  strSpikeinIndicator = paste(c_strSyntheticMicrobiome, c_strSpike,sep="")
  fPlotSpikeIn = FALSE
  strFeature = NULL
  strMetadata = NULL
  strMetadataNoLevel = NULL
  for(sLine in vsTruth)
  {
    if(!is.null(strFeature) && !is.null(strMetadata))
    {
      iFeatureIndex = which(mtrxFinal[,1] == strFeature)
      vsFeature = mtrxFinal[iFeatureIndex,]
      vdFeature = as.numeric(vsFeature[2:length(vsFeature)])

      if(strMetadata == strMetadataNoLevel)
      {
        iMetadataIndex = which(mtrxFinal[,1] == strMetadata)
        vsMetadata = mtrxFinal[iMetadataIndex,]
        vdMetadata = as.numeric(vsMetadata[2:length(vsMetadata)])
        plot(vdMetadata, vdFeature, xlab=strMetadata, ylab=strFeature)
        
      } else {

        iMetadataIndex = which(mtrxFinal[,1] == strMetadataNoLevel)
        vsMetadata = mtrxFinal[iMetadataIndex,]
        vdMetadata = as.factor(vsMetadata[2:length(vsMetadata)])
        plot(vdMetadata, vdFeature, xlab=strMetadata, ylab=strFeature)
      }
    }
    fInMicrobiome = grepl(c_strSyntheticMicrobiome,sLine,fixed=TRUE)
    fInSpikeinMicrobiome = grepl(strSpikeinIndicator,sLine,fixed=TRUE)
    # In the spike in microbiome
    if(fInMicrobiome && fInSpikeinMicrobiome){fPlotSpikeIn = TRUE}
    # Out of the spike in microbiome
    if(fInMicrobiome && !fInSpikeinMicrobiome)
    {
      fPlotSpikeIn = FALSE
      strFeature = NULL
      strMetadata = NULL
      strMetadataNoLevel = NULL
    }
    if(fPlotSpikeIn)
    {
      if(grepl(c_strFeature,sLine,fixed=TRUE))
      {
        strFeature = sLine
        strMetadata = NULL
        strMetadataNoLevel = NULL
      }
      else if(grepl(c_strMetadata,sLine,fixed=TRUE))
      {
        strMetadata = sLine
        strMetadataNoLevel = unlist(strsplit(strMetadata,"_Level_"))[1]
      }
      else
      {
         strFeature = NULL
         strMetadata = NULL
         strMetadataNoLevel = NULL
      }
    }
  }
  dev.off()  
}


funcPrepareMetadata = function(
### If the data is continuous just return it and the metadata name
### If the data is discontinuous and dummying return one randomly selected level
### and give it a name that reflects the level selected.
viSelectedMetadata,
### Indices of metadata to prepare
vdCurMetadata,
### Matrix (row major) of original metadata but only metadata to be spiked
fDummyData = TRUE,
vsFrozeLevels = NULL
### Holds the name and level of metadata that is discontinuous
### Allows one to force the selection of a specific level
){
  # Get the names of the metadata associated with this bug
  # If dummy is occuring this will be written over
  vstrSpikedMetadata = paste(c_strMetadata, viSelectedMetadata, sep='')

  if(is.null(vsFrozeLevels))
  {
    vsFrozeLevels = c()
  }
  vsCurLevels = c()

  # Dummy the data if need be
  if(fDummyData)
  {
    if(is.null(nrow(vdCurMetadata)))
    {
      vdCurMetadata = matrix(vdCurMetadata,nrow=1)
    }

    # Reset for levels
    vstrSpikedMetadata = c()
    for(iRowIndex in 1:nrow(vdCurMetadata))
    {
      # Get metadata
      vxMetadata = as.vector(vdCurMetadata[iRowIndex,])
      # If the factor data is not numeric
      if(!sum(sapply(vxMetadata,funcNumericIsInt))==0)
      {
        # Get levels of the metadata
        vsLevels = levels(as.factor(vxMetadata))

        # Add factor metadata with levels
        # Select from any level except for the first level
        # This is because the first level is always the reference level
        # in comparisons and will never be significant
        # (the other value compared to it will be).
        strLevel = sample(vsLevels[2:length(vsLevels)],size=1)
        if(length(vsFrozeLevels)>0)
        {
          strLevel = vsFrozeLevels[iRowIndex]
        } else {
          vsCurLevels = c(vsCurLevels,strLevel)
        }
        # Make the selected level 2 and those levels not selected 1
        vxMetadata[which(vxMetadata != strLevel)]=0
        vxMetadata[which(vxMetadata == strLevel)]=2
        vxMetadata[which(vxMetadata == 0)]=1
        vdCurMetadata[iRowIndex,] = vxMetadata
        vstrSpikedMetadata = c(vstrSpikedMetadata,paste(c_strMetadata,viSelectedMetadata[iRowIndex],"_",c_strLevel,"_",strLevel,sep=""))

      } else {
        if(length(vsFrozeLevels)==0)
        {
          vsCurLevels = c(vsCurLevels,NA)
        }
        # Add continuous metadata
        vstrSpikedMetadata = c(vstrSpikedMetadata,paste(c_strMetadata,viSelectedMetadata[iRowIndex],sep=""))
      }
    }
  }
  # Return which level was dummied
  if(!length(vsFrozeLevels)>0)
  {
    vsFrozeLevels = vsCurLevels
  }
  return(list(names=vstrSpikedMetadata, metadata=vdCurMetadata, vsLevels=vsFrozeLevels))
}

funcQCSpikin = function(
### Check to make sure the minimal number of spiked-in samples are met
### Asks the question is there a minimal number of overlap for metadata features with data features
### Works for dummied metadata or not
### In level based metadata, checks if there is a minimal # of nonzero bug entries with the level of interest
### In factor data not reduced to a level, checks if each level has at least  a min of non zero entries >= min # non-zeros / # metadata levels
vdCurMetadata,
### The metadata that was spiked-in, could be a vector or matrices (row major) depending on if the spikin is multivariate or not.
vdSpikedBug,
### The bug feature created with a spiked in relationship
iMinSpikedSamples,
### Minimal number of samples to be spiked in to pass QC
fDummyFactorData = TRUE
){
  if(is.null(vdSpikedBug)){return(list(PASS=FALSE, CommonCounts=c()))}

  # Will eventually hold which samples are non zero
  # Will eventually be a list of vectors
  # Multiple vectors will arise if the data is discontinuous
  # In this case each vector will be a level
  vNonZeroMetadata = c()

  # Get which bug samples are 0
  viNonZeroSpikedBugs = which(vdSpikedBug!=0)

  # Get number of row (the matrix is row major)
  dNumberMetadata = nrow(vdCurMetadata)
  if(is.null(dNumberMetadata))
  {
    dNumberMetadata = 1
    vdCurMetadata = matrix(vdCurMetadata,nrow=1)
  }

  # This is to return a pass or fail and the number of nonzero items so
  # even if one fails, the best choice can be made in a series of these calls
  fPass = TRUE
  viCommon = c()
  for(iMetadata in 1:dNumberMetadata)
  {
    vCurrentMetadata = vdCurMetadata[iMetadata,]

    # Test to see if the data is continuous or not
    if(sum(sapply(vCurrentMetadata,funcNumericIsInt))==0)
    {
      iCommon = intersect(viNonZeroSpikedBugs, which(vCurrentMetadata!=0))
      if( length(iCommon) < iMinSpikedSamples )
      {
        fPass = FALSE
      }
      viCommon = c(viCommon, length(iCommon))

    } else {
      # Check factor data
      # Dummied means there is going to be two values, 1s and 2s 2s being the level of interest.
      # So you are interested in those things not == 1 which is everything but that the level of interest.
      if(fDummyFactorData)
      {
        iCommon = intersect(viNonZeroSpikedBugs, which(vCurrentMetadata!=1))
        if( length(iCommon) < iMinSpikedSamples )
        {
          fPass = FALSE
        }
        viCommon = c(viCommon, length(iCommon) )

      } else {
        # Not dummied
        # Check for each level
        # Store each level and then set diff like normal.
        # At the very end each level will be checked if Non dummying is used.
        vstrLevels = levels(as.factor(vCurrentMetadata))
        iMinSamplesPerLevel = ceiling( iMinSpikedSamples/length(vstrLevels) )

        for( strLevel in vstrLevels )
        {
          iCommon = intersect(viNonZeroSpikedBugs, which(vCurrentMetadata==strLevel))
          if( length(iCommon) < iMinSamplesPerLevel)
          {
            fPass = FALSE
          }
          viCommon = c(viCommon, length(iCommon) )
        }
      }
    }
  }
  # All data / levels passes in this spiked relationship
  return(list(PASS=fPass, CommonCounts=viCommon))
}


funcRoundMatrix = function(
### Round a matrix
mtrxData,
### Matrix of data values to round
vdLeftOver
### Left over signal to use
){
  for( iSample in 1:ncol( mtrxData ) )
  {
    vdCur = mtrxData[,iSample]
    viRound = intersect( which( vdCur < 1 ), which( vdCur > 0 ) )
    iCompensate = length( viRound )
    if( iCompensate > 0)
    {
      mtrxData[viRound, iSample] = 1

      if( vdLeftOver[iSample] > 0 )
      {
        if( iCompensate > vdLeftOver[iSample] )
        {
          iCompensate = iCompensate - vdLeftOver[iSample]
          vdLeftOver[iSample] = 0
        } else {
          vdLeftOver[iSample] = vdLeftOver[iSample] - iCompensate
          iCompensate = 0
        }
      }
      if( iCompensate > 0 )
      {
        for( iRemove in 1:iCompensate )
        {
          viReplace = which( mtrxData[,iSample] > 1 )
          if( length(viReplace) > 0 )
          {
            iReplace = viReplace
            if( length(iReplace) > 1)
            {
              iReplace = sample(iReplace,1)
            }
            mtrxData[iReplace,iSample] = mtrxData[iReplace,iSample] - 1
          }
        }
      }
    }
  }
  return( list( Data = mtrxData, LeftOver = vdLeftOver ) )
}


funcShuffleMatrix = function(
### Shuffle in left over singal into the matrix
mtrxData,
### Matrix of values (rows are features, columns are samples)
vdFeatureMus,
### Vector of mus for the features
vdShuffleIn
### Vector of values to shuffle into each sample (positionally, first vdShuffleIn entry into first sample).
){
  print("Start funcShuffleMatrix")
  # Add back in left over counts
  vdFeatureMeans = funcGetRowMetric(mtrxData,mean)
#  vdFeatureMusPercent = vdFeatureMus/sum(vdFeatureMus)
  # Make left over an int
  viLeftOver = floor(vdShuffleIn)
  # Number of samples
  iNumberSamples = ncol(mtrxData)

  for(iShuffle in 1:iNumberSamples)
  {
    # Shuffle back in signal but not in zero locations
    if(viLeftOver[iShuffle] > 0)
    {
      # Find out which are zero locations
      viSampleZeros = which(mtrxData[,iShuffle]==0)

      for(iCountIndex in 1:viLeftOver[iShuffle])
      {
        # For every count
        # Select the feature location most differing from but not going over it's expected mean
        # If all go over just randomly select the rest
        # Get difference of expected and actual mean
        vdDifference = vdFeatureMus - vdFeatureMeans
      
        # If there are zero locations make sure they are not selected
        if(length(viSampleZeros)>0)
        {
          vdDifference[viSampleZeros] = min(vdDifference)-1
        }
        dMax = max(vdDifference)

        # Add a count to the feature location
        # If dMax is less than zero then all features have their means so, if any more counts are needed, added randomly
        # Otherwise do not add randomly, add by difference and
        iShuffleTo = 0
        if(dMax < 0)
        {
          # Do not select zero locations
          iShuffleTo = setdiff(1:length(vdDifference), viSampleZeros)
          if(length(iShuffleTo)>1)
          {
            iShuffleTo = sample(x=iShuffleTo, size=1)
          }
        } else {
          iShuffleTo = which(vdDifference == dMax)
          if(length(iShuffleTo)>1){iShuffleTo = sample(x=iShuffleTo, size=1)}
        }

        # Update the data and the means
        mtrxData[iShuffleTo,iShuffle] = mtrxData[iShuffleTo,iShuffle] + 1
        vdFeatureMeans[iShuffleTo] = vdFeatureMeans[iShuffleTo] + 1/iNumberSamples
        viLeftOver[iShuffle] = viLeftOver[iShuffle] - 1
      }
    }
    
    # Update the feature means
    vdFeatureMeans = funcGetRowMetric(mtrxData,mean)
  }
  print("Stop funcShuffleMatrix")
  return(mtrxData)
}


funcSpikeNewBug = function(
### Combine a bug data and metadata with a certain multiplier to create a known positive association
vdCurMetadata,
### Metadata to spike in
vdCurData,
### Data (bug to use in the spike-in)
multiplier,
### Multiplier to increase the signal of the metadata in the association.
### Larger values make metadata more of the signal in the spike-in.
### 1 makes an equal parts metadata and data relationship.
fZeroInflated = True
){
  # Get average and SD of each metadata
  metadata_average = funcGetRowMetric(vdCurMetadata,mean)
  metadata_sigma = funcGetSD(vdCurMetadata)

  # Get the average and SD of the feature (ignoring zeros if zero inflated)
  data_average = mean(vdCurData[which(vdCurData>0)])
  data_sigma = sd(vdCurData[which(vdCurData>0)])
  if(is.na(data_sigma)){return(NULL)}

  if(!fZeroInflated)
  {
    data_average = mean(vdCurData)
    data_sigma = sd(vdCurData)
  }

  # Feature as occurence (0,1)
  # This is used so that the features with zero keep those measurements at zero
  # Only used in a zero inflated environment so everything is set to 1 otherwise which cancels it's affect
  liOccurenceFilter = vdCurData
  liOccurenceFilter[liOccurenceFilter>0]=1
  if(!fZeroInflated)
  {
    liOccurenceFilter[1:length(liOccurenceFilter)] = 1
  }

  # Spike in the feature with the metadata (multivariate_parameter == 1 means 1 metadata spike-in for the feature)
  # Make the scaled metadata
  scaled_metadata = (vdCurMetadata - metadata_average)/metadata_sigma
  if(data_sigma!=0){ scaled_metadata = scaled_metadata * data_sigma }
  scaled_metadata = scaled_metadata + data_average
  if(!(is.null(nrow(vdCurMetadata))||(nrow(vdCurMetadata)==1)))
  {
    scaled_metadata = colSums(scaled_metadata)
  }
  # Spike in the metadata that is now scaled.
  # make the metadata sparse based on the feature ("and" function filter)
  vdSpikedBug = vdCurData + (multiplier*scaled_metadata*liOccurenceFilter)
  # Scale the spike in down so the bug count does not increase as a (whole) feature just because adding is happening
  vdSpikedBug = vdSpikedBug /((nrow(vdCurMetadata)*multiplier)+1)

  vdSpikedBug[which(vdSpikedBug<0)]=0
  return(vdSpikedBug)
}

# 4 Tests 9/4/2013
funcTruncatedRLNorm = function(
### Return draws from a random log normal distribution with a truncated tail so outliers are not introduced.
iNumberMeasurements,
### The number of measurements for this distribution
dLogMean,
### The mean of the logged distribution
dLogSD,
### The SD of the logged distribution
iThreshold = NA
### If FALSE, the outlier is given a max value and then the difference between the max value and the outlier values
### are randomly shuffled into the feature after all values have been defined.
### If TRUE, the outlier is replaced with the max value, and the difference is ignored.
){
  # Get a feature measurement one at a time checking for outliers
  # and redrawing outliers
  vdFeature = rlnorm( iNumberMeasurements, dLogMean, dLogSD )
  vdDifference = rep(0,iNumberMeasurements)

  # If the outliers are there
  if(!is.na(iThreshold))
  {
    viOutliers = which(vdFeature>iThreshold)

    # If there are outliers
    # Change the outliers to the max
    # and optionally shuffle back in the differences between the outliers and the max
    if(length(viOutliers))
    {
      vdDifference[viOutliers] = vdFeature[viOutliers]-iThreshold
      vdFeature[viOutliers] = iThreshold
    }
  }
  #Truncate negatives to zero
  vdFeature[vdFeature<0] = 0

  return(list(Feature=vdFeature, LeftOver=vdDifference))
### vdFeature: The signal (optionally truncated, log normal distribution)
### vdDifference: Left over signal that was reomved from the original vdFeature
}


funcZeroCorrectMatrix = function(
### Make sure there are no samples which sum to zero
### If there are take a measurement from a sample in the sample feature
### And adjust the left over values for read depth adjustments.
mtrxData,
### Matrix of values (rows are features, columns are samples).
vdFeatureMus,
### Vector of mus for the features
vdLeftOver
### Counts to be suffled into each sample (column) by position.
### First element in the vector for the first element in the matrix
### Here, this is used to store or take values from the samples when correcting
### for zeros.
){
  print("Start funcZeroCorrectMatrix")

  # Which samples have no signal in them
  vdZeroReadDepth = which(colSums(mtrxData)==0)
  # Find which samples have some signal in them (how many features have signal in them)
  viNotZeroOccurences = apply(mtrxData,2,function(x) length(which(x>0)))
  if(!length(viNotZeroOccurences)){viNotZeroOccurences = 0}
  # Find the average non-zero occurences in the samples, this is how many features will be given signal in a sample.
  iLowestMeasurementOccurence = floor(mean(viNotZeroOccurences))

  # If there are samples with read depth equal to zeros
  if(length(vdZeroReadDepth))
  {
    # If there are not enough features to move around to fill in zero samples
    # Add in features with signal based on the expected means
    # and remove the added signal from the feature in another sample (from the left over signal)
    if(iLowestMeasurementOccurence < 1)
    {
      # For each of the samples with zero sum
      # Give a random 1 to a random feature in reverse order of feature mean
      viFeaturesOrder = order(vdFeatureMus,decreasing=FALSE)

      # Incase there are less features than samples shift the indices so that the features
      # are cycled through the correct amount of times (which == the # samples)
      vdFeatureIndices = funcCycleVectorIndices(iElementsToPick=length(vdZeroReadDepth), iVectorLength=length(viFeaturesOrder))

      # For each zero sample
      for(iSampleIndex in 1:length(vdZeroReadDepth))
      {
        # Add a count to the zero sample and remove it from the left over signal to later shuffle in.
        mtrxData[viFeaturesOrder[vdFeatureIndices[iSampleIndex]], vdZeroReadDepth[iSampleIndex]] = 1
        vdLeftOver[vdZeroReadDepth[iSampleIndex]] = vdLeftOver[vdZeroReadDepth[iSampleIndex]] - 1
        # Get the sample with the max amount of signal
        # If there are ties randomly sample one from the group
        vdRds = colSums(mtrxData)
        viMaxRds = which(vdRds == max(vdRds))
        if(length(viMaxRds)>1){viMaxRds = sample(viMaxRds,1)}

        # Update the matrix and the left over vector
        mtrxData[viFeaturesOrder[vdFeatureIndices[iSampleIndex]],viMaxRds] = mtrxData[viFeaturesOrder[vdFeatureIndices[iSampleIndex]],viMaxRds] - 1
        vdLeftOver[viMaxRds] = vdLeftOver[viMaxRds] + 1
      }
    } else {
      # Start reordering features until there are samples with atleast one measure in each sample
      # For each sample that has 0 read depth
      for(iIndex in vdZeroReadDepth)
      {
        # Give counts to iLowestMeasurmentOccurence Features
        # Attempting to make these samples like an average sample
        for(iMeasurmentOccurenceIndex in 1:iLowestMeasurementOccurence)
        {
          # Find the sample with the most measurements
          viMeasurements = apply(mtrxData,2,function(x) length(which(x>0)))
          dMaxSample = which(viMeasurements == max(viMeasurements))
          if(length(dMaxSample)>1){ dMaxSample = sample(dMaxSample,1) }

          # Randomly select one of the signal features of the sample
          iTake = sample(which(mtrxData[,dMaxSample]>0),1)
          iTakeValue = min(mtrxData[iTake,dMaxSample], vdLeftOver[iIndex])

          # Take the value from the more measurement occuring sample and shuffle it to the 
          # zero measurement occuring sample in the same feature
          # Update the left over measurements accordingly
          vdLeftOver[dMaxSample] = vdLeftOver[dMaxSample]+ iTakeValue
          vdLeftOver[iIndex] = vdLeftOver[iIndex] - iTakeValue
          mtrxData[iTake,iIndex] = mtrxData[iTake,iIndex] + iTakeValue
          mtrxData[iTake,dMaxSample] = mtrxData[iTake,dMaxSample] - iTakeValue
        }
      }
    }
  }
  print("Stop funcZeroCorrectMatrix")
  return(list(Data=mtrxData,LeftOver=vdLeftOver))
}


# 9 Tests 9/4/2013
func_zero_inflate = function(
### Create a zero inflated log normal distribution with a specified mean and percentage of zeros.
### If you want to get the original values of mu and sd used in rlnorm, use mean(log(func_zero_inflate()))
### and sd(log(func_zero_inflate()))
dMean,
### Mean of the distribution (logged)
dPercentZeroInflated,
### Percentage of return which is zero
int_number_samples,
### The number of samples to create
dSD,
### The sd of the distribution (logged)
iThreshold = NA
### The threshold for outliers
){
  # Get feature given distribution parameters
  lFeatureData = funcTruncatedRLNorm( int_number_samples, dMean, dSD, iThreshold = iThreshold )
  vdFeature = lFeatureData$Feature

  # Zero inlate and store removed signal to the left over storage vector
  viZeroLocations = sample( 1:int_number_samples, floor( int_number_samples * dPercentZeroInflated ), replace = FALSE )
  vdLeftOver = lFeatureData$LeftOver
  if( length( viZeroLocations ) )
  {
    vdLeftOver[ viZeroLocations ] = vdLeftOver[ viZeroLocations ] + vdFeature[ viZeroLocations ]
    vdFeature[ viZeroLocations ] = 0
  }

  # Return zero-inflated truncated lognormal feature with left over signal to later be shuffled back in.
  return( list( Feature = vdFeature, LeftOver = vdLeftOver ) )
}




option_list = list(
make_option(c("-a","--variance_scale"), type="double", default = .01, help="Tuning parameter for noise in bug-bug associations"),
make_option(c("-b","--bug_to_spike"), type="integer", default=5, help="Number of bugs to correlate with others"),
make_option(c("-c","--calibrate"), type="character", default=NA, help="Calibration file for generating the random log normal data. TSV file (column = feature)"),
make_option(c("-d","--int_multiplier_delta"), type="double", default=1, help="Multiplier delta, a positive real number is expected."),
make_option(c("-e","--read_depth"), type="integer", default=8030, help="Simulated read depth for counts. A positive integer value is expected."),
make_option(c("-f","--number_features"), type="integer", default=300, help="The number of features per sample to create. A positive integer value is expected."),
make_option(c("-g","--int_max_multiplier_range"), type="double", default=2, help="Maximum value for the multiplier range. A positive real number is expected."),
make_option(c("-i","--collinear_increment"), type="integer", default=1, help="This focuses on the number of metadata spiked in to a feature. This is the delta used to create a range of values to indicate the number of spiked in metadata. This helps when creating multiple microbiomes with different levels of collinearity (or different numbers of metadata to spike in). A positive integer is expected."),
make_option(c("-j","--lefse_file"), type="character", default=NULL, help="Folder containing lefSe inputs."),
make_option(c("-k","--percent_spiked"), type="double", default=.03, help="The percent of features spiked-in. A real number between 0 and 1 is expected."),
make_option(c("-l","--minLevelPercent"), type="double", default=.1, help="Minimum percent of measurements out of the total a level can have in a discontinuous metadata (Rounded up to the nearest count). A real number between 0 and 1 is expected."),
make_option(c("-m","--max_domain_bugs"),type="integer", default=2, help="Maximum number of bugs with which each correlated bug can be associated with. A positive integer greater than 0 is expected."),
make_option(c("-n","--number_samples"), type="integer", default=50, help="The number of samples to generate. A positive integer greater than 0 is expected."),
make_option(c("-o","--max_percent_outliers"), type="double", default=.05, help="The maximum percent of outliers to spike into a sample. A real number between 0 and 1 is expected."),
make_option(c("-p","--number_metadata"), type="integer", default=5, dest='number_metadata',help="Indicates how many metadata are created. number_metadata*2 = number continuous metadata, number_metadata = number binary metadata, number_metadata = number quaternary metadata. A positive integer greater than 0 is expected."),
make_option(c("-q","--int_min_multiplier_range"), type="double", default=1, help="Minimum value for the multiplier range. A positive real number is expected."),
make_option(c("-r","--collinear_range"), type="integer", default=1, help="The maximum value for the range of number determining how many metadata to use in a feature spike-in. A positive integer greater than 0 is expected."),
make_option(c("-s","--seed"), type="integer", default=NA, help="A seed to freeze the random generation of counts/relative abundance. If left as default (NA), generation is random. If seeded, data generation will be random within a run but identical if ran again under the same settings. An integer is expected."),
make_option(c("-t","--percent_outlier_spikins"), type="double", default=.05, help="The percent of samples to spike in outliers. A real number between 0 to 1 is expected."),
make_option(c("-u","--minOccurence"), type="integer", default=4, help="Minimum counts a bug can have for the ocurrence quality control filter used when creating bugs. ( Filtering minimum number of counts in a minimum number of samples). A positive integer is expected."),
make_option(c("-v","--verbose"), action="store_false", default = TRUE, help="If True logging and plotting is made by the underlying methodology. This is a flag, it is either included or not included in the commandline, no value needed."),
make_option(c("-w","--minSample"), type="integer", default=4, help="Minimum samples a bug can be in for the ocurrence quality control filter used when creating bugs. ( Filtering minimum number of counts in a minimum number of samples). A positive integer is expected."),
make_option(c("-z","--noZeroInflate"), action="store_true", default = FALSE, help="If given, zero inflation is not used when generating a feature. This is a flag, it is either included or not included in the commandline, no value needed.")
)


main = function(
pArgs
){
lxArgs = parse_args(pArgs,positional_arguments = TRUE)

options = lxArgs[['options']]
seed = options[['seed']]
lefse_file = options[['lefse_file']]
int_base_metadata_number = options[['number_metadata']]
if(int_base_metadata_number<1) stop("Please provide the base number for metadata generation as 1 or greater.")
int_number_features = options[['number_features']]
if(int_number_features<1) stop("Please provide a number of features of atleast 1")
int_number_samples = options[['number_samples']]
if(int_number_samples<1) stop("Please provide a number of samples of atleast 1")
dPercentOutliers = options[['max_percent_outliers']]
if( (dPercentOutliers>1) | (dPercentOutliers<0) ) stop("Please provide a percent outliers in the range of 0 to 1")
dPercentOutlierSpikins = options[['percent_outlier_spikins']]
if( (dPercentOutlierSpikins>1) | (dPercentOutlierSpikins<0) ) stop("Please provide a percent spikins in the range of 0 to 1")
iReadDepth = options[['read_depth']]
if(iReadDepth < max(int_number_features, int_number_samples)) stop("Please provide a read depth of atleast equal to feature size or sample size (which ever is larger)")
int_max_multiplier_range = options[['int_max_multiplier_range']]
int_min_multiplier_range = options[['int_min_multiplier_range']]
int_multiplier_delta = abs(options[['int_multiplier_delta']])
if(int_max_multiplier_range<int_min_multiplier_range)
{
  iTemp = int_min_multiplier_range
  int_min_multiplier_range = int_max_multiplier_range
  int_max_multiplier_range = iTemp
}
if(int_min_multiplier_range<=0)stop("Please provide make sure the smallest value in the multiplier range is greater than 0.")
collinear_range = options[['collinear_range']]
if(collinear_range < 1) stop("Please provide a collinear increment greater than or equal to 1")
collinear_increment = options[['collinear_increment']]
if(collinear_increment < 0) stop("Please provide a collinear increment greater than 0")
iNumAssociations = options[['bug_to_spike']]
if(iNumAssociations<0) stop("Please provide a number of associations (bug-bug correlation) greater than or equal to 0")
dVarScale =options[['variance_scale']]
iMaxNumberCorrDomainBugs = options[['max_domain_bugs']]
dPercentMultSpiked = options[['percent_spiked']]
if( (dPercentMultSpiked>1) | (dPercentMultSpiked<0) ) stop("Please provide a percent multivariate spike in the range of 0 to 1")
strCalibrationFile = options[['calibrate']]
dMinLevelCountPercent = options[['minLevelPercent']]
if( (dMinLevelCountPercent>1) | (dMinLevelCountPercent<0) ) stop("Please provide a min level percent in the range of 0 to 1")
dMinOccurenceCount = options[['minOccurence']]
if(dMinOccurenceCount<0) stop("Please provide a min occurence greater than or equal to 0")
dMinOccurenceSample = options[['minSample']]
if(dMinOccurenceSample<0) stop("Please provide a min sample greater than or equal to 0")
if(dMinOccurenceSample>int_number_samples)
{
  dMinOccurenceSample = int_number_samples
  print(paste("The min sample (for QC) was larger than the actual sample size, reset the min sample to the sample size, minSample is now equal to number_samples which is ",int_number_samples))
}
fVerbose = options[['verbose']]
fZeroInflate = !options[['noZeroInflate']]
print("%%%%% fZeroInflate")
print(fZeroInflate)

# locational arguments
file_names = lxArgs[['args']]
strNormalizedFileName = file_names[1]
strCountFileName = file_names[2]
parameter_filename = file_names[3]
strBugBugInteractionFile = file_names[4]
if(is.na(strNormalizedFileName))
{
  strNormalizedFileName = "SyntheticMicrobiome.pcl"
}
if(is.na(strCountFileName))
{
  strCountFileName = "SyntheticMicrobiome-Counts.pcl"
}
if(is.na(parameter_filename))
{
  parameter_filename = "SyntheticMicrobiomeParameterFile.txt"
}

# seed the random number generator
if(!is.na(seed))
{
  set.seed(seed)
}

# List of associations and bugs
vParametersAssociations = c()
list_of_bugs = list()

# generate the metadata
lsMetadataInfo = func_generate_metadata(int_base_metadata_number,int_number_samples,dMinLevelCountPercent)
mat_metadata =  lsMetadataInfo[["mat_metadata"]]
metadata_parameters = lsMetadataInfo[["mtrxParameters"]]
vParametersAssociations = c(vParametersAssociations,lsMetadataInfo[["mtrxParameters"]])

# generate plain random lognormal bugs
pdf(file.path(dirname(strCountFileName),"FuncGenerateRLNorm.pdf"), useDingbats=FALSE)
# Get the fitted values for calibrating rlnorm
vdExp = NA
vdMu = NA
vdSD = NA
vdPercentZero = NA
dSDBeta = c_dSDBeta
dBetaZero = c_dBetaZero
dGrandBeta = c_dBetaGrandSD

print("Parameters BEFORE Calibration File")
print(paste("Length exp",NA,"Length vdMu", NA, "length vdSD", NA, "length vdPercentZero", NA, "dSDBeta", dSDBeta,"dBetaZero", dBetaZero, "dGrandBeta", dGrandBeta, "Read depth", iReadDepth))

if(!is.na(strCalibrationFile) & (strCalibrationFile!="NA"))
{
  # Get the fit for the data
  print("Calibrating...")
  lsFit = funcCalibrateRLNormToMicrobiome(strCalibrationFile,fVerbose)
  print("lsFit")
  print(lsFit)
  vdExp = lsFit[["exp"]]
  vdMu = lsFit[["mu"]]
  vdSD = lsFit[["sd"]]
  vdPercentZero = lsFit[["percentZero"]]
  dSDBeta  = lsFit[["dSDBeta"]]
  dBetaZero = lsFit[["dZeroBeta"]]
  dGrandBeta = lsFit[["dGrandBeta"]]
  iReadDepth = lsFit[["dAverageReadDepth"]]
#  int_number_features = lsFit[["iFeatureCount"]]
#  dMaxCount = lsFit[["dMaxCount"]]
}
print("Parameters AFTER Calibration File (if no calibration file is used, defaults are shown")
print(paste("Length exp",length(vdExp),"Length vdMu", length(vdMu), "length vdSD", length(vdSD), "length vdPercentZero", length(vdPercentZero), "dSDBeta", dSDBeta,"dBetaZero", dBetaZero, "dGrandBeta", dGrandBeta, "Read depth", iReadDepth))

mat_random_lognormal_bugs = func_generate_random_lognormal_matrix(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberCounts=dMinOccurenceCount, iMinNumberSamples=dMinOccurenceSample, iReadDepth=iReadDepth, vdExp=vdExp, vdMu=vdMu, xPercentZero=vdPercentZero, vdSD=vdSD, fZeroInflate=fZeroInflate, dBetaSD=dSDBeta,  dBetaZero=dBetaZero, dBetaGrandSD=dGrandBeta, fVerbose=fVerbose)

lefse_lognormal = NULL
if(!is.null(lefse_file))
{
  lefse_file = dirname(lefse_file)
  lefse_lognormal = func_generate_lefse_matrices(lefse_file, metadata_parameters, int_number_features, int_number_samples,mat_metadata, mat_random_lognormal_bugs[['mat_bugs']], 'lognormal')
}

vParametersAssociations = c(vParametersAssociations,mat_random_lognormal_bugs[["mtrxParameters"]])
list_of_bugs[[length(list_of_bugs) + 1]] = mat_random_lognormal_bugs[["mat_bugs"]]
hist(as.vector(mat_random_lognormal_bugs[["mat_bugs"]]), main="Final: Log normal matrix")
dev.off()

# generate random lognormal with outliers
pdf(file.path(dirname(strCountFileName),"LognormalWithOutliers.pdf"), useDingbats=FALSE)
mat_random_lognormal_outliers_bugs = func_generate_random_lognormal_with_outliers(int_number_features=int_number_features, int_number_samples=int_number_samples,dMaxPercentOutliers=dPercentOutliers,dPercentSamples=dPercentOutlierSpikins,mtrxBugs=mat_random_lognormal_bugs[["mat_bugs"]],fVerbose=fVerbose)
lefse_outliers = NULL
if(!is.null(lefse_file))
{
  lefse_outliers = func_generate_lefse_matrices(lefse_file, metadata_parameters,int_number_features,int_number_samples,mat_metadata, mat_random_lognormal_outliers_bugs[['mat_bugs']], 'outliers')
}
vParametersAssociations = c(vParametersAssociations,mat_random_lognormal_outliers_bugs[["mtrxParameters"]])
list_of_bugs[[length(list_of_bugs) + 1]] = mat_random_lognormal_outliers_bugs[["mat_bugs"]]
hist(as.vector(mat_random_lognormal_outliers_bugs[["mat_bugs"]]), main="Final: Log normal matrix with outliers")
dev.off()

# Holds key words for the feature names of the microbiomes
lsMicrobiomeKeys = c(c_strRandom,c_strOutlier)

# There are 4 groups of metadata (2 continuous, binary, and quarternery)
number_metadata = c_iCountTypesOfMetadata*int_base_metadata_number

# Generate random lognormal with varying amounts of spikes
for (how_many_multivariates in seq(1,collinear_range,collinear_increment))
{
  lviMetadata = NULL
  liData = NULL
  lsLevels = NULL

  for (mult in seq(int_min_multiplier_range, int_max_multiplier_range, int_multiplier_delta))
  {
    pdf(file.path(dirname(strCountFileName),paste("SpikeIn_n_", how_many_multivariates,"_m_", mult,".pdf",sep="")), useDingbats=FALSE)
    mat_random_lognormal_multivariate_spikes = func_generate_random_lognormal_with_multivariate_spikes(int_number_features=int_number_features, int_number_samples=int_number_samples,  percent_spikes=dPercentMultSpiked, multiplier=mult, metadata_matrix=mat_metadata, multivariate_parameter=how_many_multivariates, dMinLevelCountPercent=dMinLevelCountPercent, mtrxBugs=mat_random_lognormal_bugs[["mat_bugs"]],fZeroInflated=fZeroInflate, lviFrozeMetadataIndices=lviMetadata, liFrozeDataIndicies=liData, lsFrozeLevels=lsLevels, fVerbose=fVerbose)
    mat_random_lognormal_multivariate_spikes_bugs = mat_random_lognormal_multivariate_spikes[["mat_bugs"]]
    lviMetadata = mat_random_lognormal_multivariate_spikes[["MetadataIndices"]]
    liData = mat_random_lognormal_multivariate_spikes[["DataIndices"]]
    lsLevels = mat_random_lognormal_multivariate_spikes[["Levels"]]

    if(!is.null(lefse_file))
    {
      lefse_spike = func_generate_lefse_matrices(lefse_file, metadata_parameters,int_number_features,int_number_samples,mat_metadata, mat_random_lognormal_multivariate_spikes_bugs, paste('multivariate_n_', how_many_multivariates, '_m_', mult,sep=""))
    }
    # generate known associations for random lognormal with spikes
    vParametersAssociations = c(vParametersAssociations,mat_random_lognormal_multivariate_spikes[["m_parameter_rec"]])
    list_of_bugs[[length(list_of_bugs)+1]] = mat_random_lognormal_multivariate_spikes_bugs
    lsMicrobiomeKeys[[length(lsMicrobiomeKeys)+1]] = paste(c_strSpike,"n",how_many_multivariates,"m",sub(".","_",paste(mult),fixed=TRUE), sep="_")

    hist(as.vector(mat_random_lognormal_multivariate_spikes_bugs), main=paste("Final: Spiked matrix n_", how_many_multivariates,"_m_", mult,sep=""))
    dev.off()
  }
}

# Add bug associated bug microbiome
lsBugBugInfo = func_generate_bug_bug_spiking_matrix(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberSamples=dMinOccurenceSample, iReadDepth=iReadDepth, iMinNumberCounts=dMinOccurenceCount, vdExp=vdExp, vdMu=vdMu, vdSD=vdSD, dPercentZero=vdPercentZero, fZeroInflate=fZeroInflate, dBetaSD=dSDBeta, dBetaZero=dBetaZero, fVerbose=fVerbose, dVarScale=dVarScale,iNumAssociations=iNumAssociations,iMaxNumberCorrDomainBugs=iMaxNumberCorrDomainBugs)
list_of_bugs[[length(list_of_bugs)+1]] = lsBugBugInfo[["mtrxBugs"]]
vParametersAssociations = c(vParametersAssociations,lsBugBugInfo[["vStrParameters"]])
lsMicrobiomeKeys[[length(lsMicrobiomeKeys)+1]] = c_strBugBugAssocations

# preallocate final pcl  matrix
final_matrix = matrix(data=NA,nrow=(number_metadata+int_number_features*length(list_of_bugs))+1, ncol=(int_number_samples+1))
final_matrix[1,1] = '#SampleID'
final_matrix[1,2:(int_number_samples + 1)] = paste('Sample',1:int_number_samples,sep='')
final_matrix[2:(number_metadata+1),1] = paste(c_strMetadata,1:number_metadata,sep='')
vdDim = dim(mat_metadata)
mat_metadata[(floor(vdDim[1]/2)+1):vdDim[1],] = paste("Group_",mat_metadata[(floor(vdDim[1]/2)+1):vdDim[1],],sep="")
final_matrix[2:(number_metadata+1),2:(int_number_samples+1)] = mat_metadata

# Make a matrix for counts (use the other for normalized)
mtrxFinalCounts = final_matrix

start = 2 + number_metadata
end = (2+number_metadata) + (int_number_features-1)
iFirstData = start

for(iMatIndex in 1:length(list_of_bugs))
{
  final_matrix[start:end,1] = paste(paste(c_strFeature,lsMicrobiomeKeys[iMatIndex],sep="_"),1:int_number_features,sep='_')
  # Normalize each column by thier sum and add to output
  final_matrix[start:end,2:(int_number_samples+1)] = funcNormalizeMicrobiome(list_of_bugs[[iMatIndex]])
  mtrxFinalCounts[start:end,1] = paste(paste(c_strFeature,lsMicrobiomeKeys[iMatIndex],sep="_"),1:int_number_features,sep='_')
  mtrxFinalCounts[start:end,2:(int_number_samples+1)] = list_of_bugs[[iMatIndex]]
  start = start + int_number_features
  end = end + int_number_features
}

# Plot spike-ins before and after normalization
funcPlotSpikeins(vsTruth = vParametersAssociations, mtrxFinal = mtrxFinalCounts, strFileDir = dirname(strCountFileName), strFileName = "SpikeinsBeforeNormalization.pdf")
funcPlotSpikeins(vsTruth = vParametersAssociations, mtrxFinal = final_matrix, strFileDir = dirname(strCountFileName), strFileName = "SpikeinsAfterNormalization.pdf")

# Write the table as normalized counts
write.table(final_matrix, file=strNormalizedFileName, quote=FALSE,row.names=FALSE, col.names=FALSE,sep = '\t')

# Write the tables as counts
mtrxCounts = final_matrix
write.table(mtrxFinalCounts, file=strCountFileName, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')

# Write truth tables
write.table(as.matrix(vParametersAssociations), file=parameter_filename,quote=FALSE, row.names=FALSE, col.names=FALSE,sep='\t')
}

# This is the equivalent of __name__ == "__main__" in Python.
# That is, if it's true we're being called as a command line script;
# if it's false, we're being sourced or otherwise included, such as for
# library or inlinedocs.
if( identical( environment( ), globalenv( ) ) &&
	!length( grep( "^source\\(", sys.calls( ) ) ) ) {
	main(OptionParser(usage="synthetic_datasets_script.R [options] NormalizedFile(Optional) CountFile(Optional) TrueFile(Optional)",option_list=option_list) ) }
