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

iLoopingControlIncrement = 1000
c_iCountTypesOfMetadata = 4
# Simulates an occurence filter
c_strDefaultMarkerColor = "black"
c_strOutlierColor = "red"
c_strFeatureOutlier = "cyan"
c_fIgnoreZerosInOutliers = TRUE
### Will not allow zeros to be used as the min value in swapping unless they are needed to fulfill the number of
### swaps (if there are a whole bunch of zeros, some zeros may be needed or no swapping can be performed).
c_fDummyFactorData = TRUE
### For Factor data only, selects a level and sets all of that level to 1 and others to 0
c_fRarify = FALSE
### Rarifies the samples to the min sample count
### Binary and quarternary metadata are draw from uniform distributions, these are the bounds
c_dRunifMin = .4
c_dRunifMax = .5

# Todo need to be in a constants file
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
c_strFirstQuartile = "1st Qu."
c_strThirdQuartile = "3rd Qu."

# These variables are associated with settings for
# Calculating the SD and percent Zero based on the mu.
# The constants for these variables have been estimated from
# a real (IBD) data set. These constants are used unless a
# calibration file is given which estimates the values in the
# same manner as the constants were estimated.
# The beta for fitting
strBeta = "Beta"
# The beta estimate for the relationship between SD and mu
c_dSDBeta = 0.634072110412161
# The intecept estimate for the relationship between SD and mu
c_dSDIntercept = 0.777602075408279
# Beta for the relationship between mu and zero percent
c_dBetaZero = -0.0127906473815179
# Zero intercept
c_dZeroIntercept = -0.0518608591994982
# The grand mu of the mus of all the bug distributions
c_dGrandMu = 1.0624445524326
# The SD of the mus of all the bug distributions
c_dGrandSD = 2.53216541658158
### The beta for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
c_dBetaGrandSD = 0.3269432 # 0.281797820896239
### The beta for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
c_dInterceptGrandSD = 1.711547 # 1.82333891862979
### The max count for an entry
c_dMaxCount = 6135

# Tested 1 - Generates RLNorm data set
funcCreateRLNormSyntheticCalibrationFile = function(
### File to write the data to
strWriteFile,
### Indicates to create a sparse, zero-inflated data set
fZeroInflate,
### Indicates if pdfs are plotted for checking
fVerbose = FALSE
){
  c_iNumberSamples = 100
  c_iNumberFeatures = 1000

  # Create default matrix
  mtrxBugs <- matrix(data=NA,nrow=c_iNumberSamples,ncol=c_iNumberFeatures)

  # Expectation vector
  vdExpVector = rlnorm(c_iNumberFeatures,c_dGrandMu,c_dGrandSD)

  # Generate each feature based on the expectation vector
  # Here each vector entry is the expectation of a feature.
  # We do not have the mu or the sd that generates this expectation
  # But we do have the relationship between the mu and the sd
  # funcGetParamsForExpectation explores mus and sds that are related to
  # that mu to see what expectation they will generate.
  vdDifferences = c()
  for(iFeatureIndex in 1:c_iNumberFeatures)
  {
    iZeroCount = 0
    viZeroIndices = c()
    vdFeature = rep(NA,c_iNumberSamples)

    lExpResults = funcGetParamsForExpectation( dMuInit=log(1), dReadDepthTarget=vdExpVector[iFeatureIndex], dSampleCount=c_iNumberSamples, dFeatureBetaSD=c_dSDBeta, dFeatureInterceptSD=c_dSDIntercept )
    dFeatureMu = lExpResults[["dBestMu"]]
    dFeatureSD = lExpResults[["dBestSD"]]

    # Generating a tsv file, so features are columns
    vdFeature[which(is.na(vdFeature))] = rlnorm(c_iNumberSamples-iZeroCount, log(dFeatureMu), log(dFeatureSD))
    vdDifferences = c(vdDifferences,(mean(vdFeature[which(vdFeature != 0)])-vdExpVector[iFeatureIndex])/vdExpVector[iFeatureIndex])

    # Shuffle counts to make zeros
    if(fZeroInflate)
    {
      dPercentZero = funcEstimatePercentZero(dFeatureMu,c_dBetaZero,c_dZeroIntercept)
      iZeroCount = floor(dPercentZero*length(vdFeature))
      viZeroIndices = sample(1:length(vdFeature), iZeroCount, replace=FALSE)
      vdFeature = funcAdjSeqForZeros(vdFeature,viZeroIndices)
    }

    # Store feature
    mtrxBugs[,iFeatureIndex] = vdFeature
  }

  # Remove NA columns
  # I don't really like that I have to remove NAs, I would rather not generate them
  # This happends when atttempting to use rlnorm with a mu or sd < 1
  # These are so small they are not interesting but nonetheless I do not like this.
#  mtrxBugs = mtrxBugs[,-1*which(is.nan(mtrxBugs[1,]))]

  # Log the pdf for the data set
  if(fVerbose)
  {
    barplot(height=apply(mtrxBugs,1,sum), main=paste("LNORM Data Set, Exp Read Depth",sum(vdExpVector)))
    abline(mean(apply(mtrxBugs,1,sum)),0,col="orange")
    abline(sum(vdExpVector),0,col="grey")
  }

  # Add sample names
  row.names(mtrxBugs)=paste("Sample",1:c_iNumberSamples,sep="_")

  # Write to file
  write.table(mtrxBugs,sep="\t",file=strWriteFile, col.names=NA, row.names=TRUE)
}

funcStandardize = function(adValues)
### Returns a numeric vector as normal deviates
{
  return((adValues/sd(adValues))-mean(adValues))
}

funcIdentifyNonZeroOutliers = function(
### Returns the indices of outliers in a vector of values
vdValues
){
  #!# Should this be logged?
  tbleSum = summary(vdValues[which(vdValues!=0)])
  dThirdQuartile = tbleSum[c_strThirdQuartile]
  dFirstQuartile = tbleSum[c_strFirstQuartile]
  dMax = (dThirdQuartile-dFirstQuartile)+dThirdQuartile
  return(which(vdValues>dMax))
}

# Tested
funcNormalizeMicrobiome = function(
### Normalize the data assuming that the 
### Samples are columns and the features (bugs) are rows.
### Zero columns are preserved as zero.
mtrxMicrobiome
### Matrix or dataframe microbiome
){
  mtrxRelAb = mtrxMicrobiome/unlist(lapply(colSums(mtrxMicrobiome),function(x) rep(x,nrow(mtrxMicrobiome))))
  viZero = which(apply(mtrxMicrobiome,2,sum)==0)
  mtrxRelAb[,viZero]=rep(0,nrow(mtrxRelAb))
  return(mtrxRelAb)
}

# Tested
funcGetExp <- function(
### Given the Mu and SD, get the expectation
### Will need to 
# Expects the values not logged
# Log is base exp(1)
# Returns a value that is consistent and can be directly placed into
# the funcGetMu and funcGetExpSD functions. To use to describe the actual lrnorm, exp(return_value)
# Otherwise this is in norm terms
dMu,
dSD
){
  return(exp(log(dMu,exp(1))+.5*(log(dSD,exp(1))^2)))
}

# Tested
funcGetMu <- function(
### Get the mu given SD and an expectation
# Expects the values to be derived from a logged rlnorm
# so sd(log(rlnorm(10000,3,2))) to calculate the value
# output from these functions can be directly used in the rlnorm
# Log is base exp(1)
# To use these values in rlnorm you would log them
dEx,
dSD
){
  return(exp( -1*(( (log( dSD, exp(1) )^2) /2 )-log( dEx, exp(1) )) ))
}

#funcGetExpSD <- function(
#### Get the SD given the Mu and an expectation
## To use these values in rlnorm you would log them
## so mean(log(rlnorm(10000,3,2))) to calculate the value
## output from these functions can b directly used in the rlnorm
#dEx,
#dMu
#){
#  return( exp(sqrt( 2*( log(dEx,exp(1))-log( dMu,exp(1) ) ) )))
#}

#funcGetLogMu <- function(
#dMu,
#dSD
#){
#  return( log( dMu/(sqrt (dSD^2/ (dMu^2) +1) ) ) )
#}

#funcGetLogSD <- function(
#dMu,
#dSD
#){
#  return(sqrt(log( ( dSD^2/ (dMu^2) )+1 )))
#}

funcEstimateFeatureSD<-function(
### Estimate the SD given the Mu and parameters modeling the relationship between the mu and the SD
dMu,
### The measured mean of feature values
dBetaSD,
### The beta for the relationship between the mu and the SD
dSDIntercept
### The intercept for the relationship between the mu and the SD
){
#  return(max(c(log(dMu)*dBetaSD+dSDIntercept,funcEstimateMinimumFeatureSD(dBetaSD,dSDIntercept))))
  return(log(dMu)*dBetaSD+dSDIntercept)
}

funcEstimateMinimumFeatureSD<-function(
### This is the opposite of the funcEstimateSD function
dBetaSD,
### The beta for the relationship between the mu and the SD
dSDIntercept
### The intercept for the relationship between the mu and the SD
){
  return(exp(1)^(((1+c_dMinimumSD)-dSDIntercept)/dBetaSD))
}

funcEstimateGrandSD<-function(
### Estimate the grand SD given the grand Mu and parameters modeling the relationship between the grand mu and grand sd
dMu,
dBetaGrandSD,
dInterceptGrandSD
){
  return(dMu * dBetaGrandSD + dInterceptGrandSD)
}

funcEstimatePercentZero <- function(
### Estimate the percent zero given the mu and parameters modeling the relationship between the mu and the percent zero
dMu,
dBetaZero,
dInterceptZero
){
  return(exp(1)^(exp(1)^dMu*dBetaZero+dInterceptZero))
}

# Tested
funcAdjSeqForZeros <- function(
### Adjust sequence by removing counts from the given location to make them zeros
### and shuffling out the counts to other sequence locations
vdFeature,
### The measurements
viZeros
### The locations where zeros should occur
){
  viSignalLocations = setdiff(1:length(vdFeature),viZeros)
  for(iZeroLocation in viZeros)
  {
    viAddCount = sample(viSignalLocations,floor(vdFeature[iZeroLocation]),replace=TRUE)
    for(iLocation in viAddCount){ vdFeature[iLocation] = vdFeature[iLocation]+1 }
    vdFeature[iZeroLocation]=0
  }
  return(vdFeature)
}

func_generate_lefse_matrices <- function(
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

  l_lefse_matrix <- list()
  lefse_matrix_1 <- matrix(data=NA, nrow=(int_number_features+2), ncol=(int_number_samples+1))
  for (meta_index in lsIndicesOfMicrobiomes)
  {
    lefse_matrix_1[1,1] = 'SampleID'
    lefse_matrix_1[1,2:(int_number_samples + 1)] <- paste('Sample',1:int_number_samples,sep='')
    lefse_matrix_1[2,1] <- 'Metadata'
    lefse_matrix_1[2,2:(int_number_samples+1)] <- mat_metadata[(meta_index-1),]
    lefse_matrix_1[3:(int_number_features+2),1] <- paste("Bug", 1:int_number_features, sep='')
    lefse_matrix_1[3:(int_number_features+2),2:(int_number_samples+1)] <- mat_bugs
    lefse_filenames <- file.path(lefse_file, paste('LefSe_',dataset_type, '_metadata_',(meta_index-1), '.1-pcl', sep=''))
    write.table(lefse_matrix_1,file=lefse_filenames, quote=FALSE, row.names=FALSE, col.names=FALSE,sep='\t')
    l_lefse_matrix[[length(l_lefse_matrix)+1]] <- lefse_matrix_1
  }
  return(l_lefse_matrix)
}

# Tested
funcGetSD <- function(
### Get standard deviation for vector or matrix
lxValues
### matrix or vector
){
  if(is.null(dim(lxValues)[1]))
  {
    return(sd(lxValues,na.rm=TRUE))
  }
  return( apply(lxValues, 1, sd, na.rm=TRUE))
}

# Tested
# Will fail if all are integer if using to detect metadata
# Returns a TRUE (1) for each int value so if
# the return vector's sum is 0, the data is not numeric
funcNumericIsInt = function(
### Tests to see if a numeric is an interger
dValue
### The numeric value to check
){
  return(paste(as.integer(dValue),sep="")==paste(dValue,sep=""))
}

# Tested
funcGetRowMean <- function(
### Get mean for vector or matrix
lxValues
### Values of metadata which may be discrete (in which after the _ is a number) or a numeric
){
  if(is.null(dim(lxValues)[1]))
  {
    return(mean(lxValues,na.rm=TRUE))
  }
  return( apply(lxValues, 1, mean, na.rm=TRUE))
}

# Reviewed but not tested
funcGenerateExpVectorParameters <- function(
### Need several pieces of data to generate the expectation vector
### Need the SD and the Mu that will generate the vector using rlnorm
### Need to make sure that this reasonably models the original distribution
### Need to make sure the vectors of mus generated from the measured parameters
### create a distribution that maintains the read depth (sum of the original expectations or the generated values)
vdExpectations,
### These are the mus of each of the features from the calibration file.
### These are untransformed (not logged) and the expectation
### This will not be a sparse/zero-inflated distribution given that it is the expectation of
### each feature.
fVerbose = FALSE
### If pdfs should be made
){
  print("Start funcGenerateExpVectorParameters")
  # Bootstrap the Mu of the expectations and the Mu of the features
  # Do this by randomly selecting a number of features to select
  # Then randomly selecting that number of features from the calibration data set.
  # Then model with lm (because it is linear from the figure).
  # Remember when defining the relationship between the grand mu and the grand sd
  # directly taking the mu will just give you the relationship with the expectation and the sd,
  # not the mu. We want the relationship between the mus and the SD
  vdBootGrandMu = c()
  vdBootGrandSD = c()
  vdBootGrandExp = c()
  vdPredictedGrandExp = c()
  for(iIndex in 1:1000)
  {
    # The minimum sample for the feature length here is 5, this is because a sd of 1 is NA and
    # so the calculation would not be helpful. As well, approximating mean and sd with a very small
    # subsample will be highly volatile. I will use the other more stable measurements of more
    # samples to estimate this smaller more volatile group. 
    vdBootFeatures = sample(vdExpectations, sample(5:length(vdExpectations), 1), replace=TRUE)

    # Since all the features are not being selected, the case can arise that one selects very small
    # features without the other more prevalent features. This can cause a very small mean which can cause
    # problems with the log and rlnorm calls. These are also not important groupings of features so they are ignored and
    # the relationship in the low end can be inferred by the others.
    if(log(mean(vdBootFeatures))>0)
    {
      vdBootGrandExp = c(vdBootGrandExp,log(mean(vdBootFeatures)))
      vdBootGrandSD = c(vdBootGrandSD,sd(log(vdBootFeatures)))
      vdBootGrandMu = c(vdBootGrandMu,funcGetMu(vdBootGrandExp[length(vdBootGrandExp)],vdBootGrandSD[length(vdBootGrandSD)]))
      vdPredictedGrandExp = c(vdPredictedGrandExp, log(mean(rlnorm(length(vdBootFeatures),vdBootGrandMu[length(vdBootGrandMu)],vdBootGrandSD[length(vdBootGrandSD)]))))
      if(is.nan(vdPredictedGrandExp[iIndex])){print(paste("Feature c(",paste(vdBootFeatures, collapse=","),")"))}
    }
  }

  # Get the relationship between the Mu of Mus and the SD of Mus
  lmGrand = lm(vdBootGrandSD ~ vdBootGrandMu)
  dGrandBetaGuess = as.numeric(coef(lmGrand)["vdBootGrandMu"])
  dGrandInterceptGuess = as.numeric(coef(lmGrand)["(Intercept)"])

  # Parameters to generate the expection vector
  dMu = mean(log(vdExpectations))
  dSD = sd(log(vdExpectations))

  # Document function
  if(fVerbose)
  {
    # Distribution of Mus
    # Hypothetical guess at generating Mus
    vdGuess = rlnorm(length(vdExpectations), dMu, dSD)
    hist(vdExpectations,main=paste("Untransformed Mu Vector from Calibration File", "Sum", sum(vdExpectations)))
    plot(sort(log(vdGuess)),sort(log(vdExpectations)), xlab="Logged Predicted", ylab="Logged Measured", main=paste("Real Mus vs Estimated Mus", "Est Sum", sum(vdGuess)))

    # Plot the relationship between Mu and SD
    # Plot the expectation predicted from this relationship and the actual relationship
    plot(vdBootGrandMu, vdBootGrandSD, main="Boosted Grand Mu v Grand sd")
    points(x=vdBootGrandMu, y=funcEstimateGrandSD(vdBootGrandMu,dGrandBetaGuess,dGrandInterceptGuess), col="tan")
    plot(vdPredictedGrandExp,vdBootGrandExp, xlab="Log Estimated", ylab="Log Measured", main="Check difference between measured and estimated expectations")
  }

  # Return the SD and Mu which would generate this vector with rlnorm
  # Return the relationship between these
  # Return the original read depth
  # Returned values are already logged and can be used directly in rlnorm
  return(list(GrandMu=dMu, GrandSD=dSD, GrandExp=log(mean(vdExpectations)), dReadDepth=sum(vdExpectations), iFeatureCount=length(vdExpectations), GrandSDBeta=dGrandBetaGuess, GrandSDIntercept=dGrandInterceptGuess))
}

# Reviewed but not tested.
funcCalibrateRLNormToMicrobiome <- function(
### Given a TSV file
### Parameters for distribution generatio are given
### To be used in matrix generation
### Estimated parameters include
### SD excluding zeros
### Mus excluding zeros
### Percent zeros
### The beta for estimating the SD given the mu
### The beta for estimating the percent zero given the mu
### The intercept for estimating the SD given the mu
### The intercept for esitmating the percent zero given the mu
### All Mus, SD, grand MU, and grand SD are ready for rlnorm (have been measured from a logged (rlnorm)
sCalibrationFile,
### File to be read in and used to calibrate constansts and relationships in the underlying data.
fVerbose = FALSE
### Flag to turn on logging and pdf creation
){
  # Read in file
  print("Reading file.")
  dfData = read.table(sCalibrationFile)
  row.names(dfData) = dfData[[1]]
  dfData = dfData[-1,-1]

  # Get read depth of the samples (From a tsv file, samples = rows)
  ldReadDepths = sapply(1:nrow(dfData), function(x) sum(as.numeric(as.matrix(dfData)[x,])))

  # Get max count
  dMaxCount = max(as.numeric(as.vector(as.matrix(dfData))))

  # Get the vector of Mus and SD ignoring 0s
  # Mus are also given from the untransformed distribution
  # Logged will be used to make features, they can be directly used in the rlnorm function
  # Not logged will be used to later be logged and estimate a grand mean and grand SD for the initial distribution,
  # Given that every point in the larger vector is the mu of the feature vectors.
  # Which could then be used by rlnorm to generate a each feature if needed.
  # This is not needed here but the calculation of this value is useful. It is used as an initial value for synthetic creation
  # of the initial vector of feature mus.
  vdSD = c()
  vdMu = c()
  vdPercentZero = c()
  vdUntransformedExp = c()

  for(iIndex in 1:ncol(dfData))
  {
    # Get the percent zero before removing zeros for other measurements
    vdCur = as.numeric(as.vector(as.matrix(dfData[iIndex])))
    vdUntransformedExp = c(vdUntransformedExp,mean(vdCur))
    vdPercentZero = c(vdPercentZero, length(which(vdCur==0))/length(vdCur))
#    iLocations = which(vdCur!=0)
    vdCur = vdCur[which(vdCur!=0)]
    
    #### Note
    #### rlnorm needs a mean and sd from a logged rlnorm distribution which would match this
    #### without further manipulation. The "mean" in the formula is actually not the expectation
    #### The expectation is e^mu+.5*sd^2 so this is always more than mu.
    vdSD = c(vdSD, sd(log(vdCur)))
    vdMu = c(vdMu, mean(log(vdCur)))
#    vdCur[iLocations] = log(vdCur[iLocations])
#    vdSD = c(vdSD, sd(vdCur))
#    vdMu = c(vdMu, mean(vdCur))
  }

#  print("calibrated mu")
#  print(vdMu)
#  print("vdSD")
#  print(vdSD)
#  print("vdPercentZero")
#  print(vdPercentZero)

  # Estimate the distribution parameters from the expectation vector
  # Includes the relationship between the grand mu and grand sd
  # The grand mu, grand expectation (logged) and the grand sd
  lParams = funcGenerateExpVectorParameters(vdUntransformedExp)

  ##### SD and Mu
  # Log to make a linear relationship
  vdMuLog = log(vdMu)

  lmod = lm(as.formula(paste("vdSD","vdMuLog",sep="~")))
  dBetaSD = coef(lmod)["vdMuLog"]
  dInterceptSD = coef(lmod)["(Intercept)"]

 # modNLS = nls(vdSD ~ exp(1)^(vdMu*dBetaSDCoef),start=list(dBetaSDCoef=dBetaSD))
 # dBetaSD = coef(summary(modNLS))["dBetaSDCoef","Estimate"]
 # plot(vdSD, vdMu, main="Estimating Zero with Mu: NLS predicted points")
 # points(x=exp(1)^(vdMu*dBetaSD),y=vdMu,col="red")

  #### Percent Zero and Mu
  ### When using NLS, first supply witha guess (took from a transformed lm)
  ### Then take the LM guess and do the NLS
  dBetaZero = NA
  dInterceptZero = NA
  if(sum(vdPercentZero)>0)
  {
    lmod = lm(asin(sqrt(vdPercentZero)) ~ vdMu)
    dBetaZero = coef(lmod)["vdMu"]
    dInterceptZero = coef(lmod)["(Intercept)"]

    # Using nonlinear fit and giving the linear fit as the first choice.
    modNLS = nls(vdPercentZero ~ exp(1)^(exp(1)^vdMu*dBetaZZ+dBetaInter),start=list(dBetaZZ=dBetaZero,dBetaInter=dInterceptZero))
    dBetaZero = coef(summary(modNLS))["dBetaZZ","Estimate"]
    dInterceptZero = coef(summary(modNLS))["dBetaInter","Estimate"]
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
    plot(vdMu, vdSD, main="Estimating SD with Mu: baseline")
    points(x=vdMu,y=sapply(vdMu, function(x) funcEstimateFeatureSD(x,dBetaSD,dInterceptSD)),col="orange")
    plot(vdMuLog, vdSD, main="Estimating SD with Mu: Linear relationship, predicted points")
    points(x=vdMuLog,y=sapply(vdMu, function(x) funcEstimateFeatureSD(x,dBetaSD,dInterceptSD)),col="green")

    # Percent zero #!#
    plot(vdMu, vdPercentZero, main="Estimating Percent Zero: Original baseline relationship 1")
    plot(exp(1)^vdMu, vdPercentZero, main="Estimating Percent Zero: Original baseline relationship 2")
    plot(exp(1)^vdMu, vdPercentZero, main="Estimating Zero with Mu: NLS predicted points 1")
    points(x=exp(1)^vdMu,y=funcEstimatePercentZero(vdMu,dBetaZero,dInterceptZero),col="violet")
  }

  # Return
  # When returning the grand Mu remember that you are returning the Mu that gives the expectation for the Mus
  # given the rlnorm function so this is different than the mus measured in the logged distribution (rlnorm)
  return(list(mu=vdMu, sd=vdSD, dSDBeta=dBetaSD, dSDIntercept=dInterceptSD, percentZero=vdPercentZero,
              dZeroBeta=dBetaZero, dZeroIntercept=dInterceptZero, dGrandMu=lParams$GrandMu,
              dGrandSD=lParams$GrandSD, dGrandBeta=lParams$GrandSDBeta, dGrandIntercept=lParams$GrandSDIntercept,
              dAverageReadDepth=mean(ldReadDepths), dMaxCount = dMaxCount))
}

func_generate_bug_bug_spiking_matrix <- function(
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
vdMu = NA,
### Vector of Mu for the original mu distribution (means of features) if not supplied, one will be generated by rlnorm
vdSD = NA,
dPercentZero = NA,
fZeroInflate = TRUE,
### Controls if zero inflation is used.
dBetaSD = c_dBetaSD,
### The beta for the relationship between the Mu and the SD
dSDIntercept = c_dInterceptSD,
### The intercept for the relationship between the Mu and the SD
dBetaZero = c_dBetaZero,
### The beta for the relationship between the Mu and the percent zero
dInterceptZero = c_dZeroIntercept,
fVerbose = FALSE
### If true, plotting and logging occur
){
  # The parameters which will eventually be passed in; the current numbers represent arbitrary first tries
  dVarScale = .01  # The scaling parameter for the VARIANCE (the noise added will have variance dVarScale*var(bug))
  iNumberAssociations = 2 # The number of correlation "structures" to introduce
  # The names for the next two parameters come from viewing the correlation structure as a function, with the domain
  # consisting of independent bugs and the range being one bug dependent on the bugs in the domain
  iMaxNumberCorrRangeBugs  = 1   # The maximum number of correlated bugs to make; the minimum is 0
  iMaxNumberCorrDomainBugs = 2   # The maximum number of bugs each of the N correlated bugs can be correlated with; the minimum is 1
  iNumCorrRangeBugs = sample(seq(1,iMaxNumberCorrRangeBugs),1)  # The true number of correlated bugs which will be made
  viNumCorrDomainBugs = sample(seq(1,iMaxNumberCorrDomainBugs),iNumCorrRangeBugs) # The number of bugs each of the correlated bugs will be correlated with
  boolIndepCorr = TRUE      # A boolean indicating whether the correlated bugs should be independent of each other

  # Generating the indices for all bugs concerned
  if(iMaxNumberCorrRangeBugs>0){
    viAvailableIndices = seq(1,int_number_features)
    viCorrRangeBugsIdx = sample(viAvailableIndices,iNumCorrRangeBugs,replace=FALSE)
    viAvailableIndices = viAvailableIndices[-viCorrRangeBugsIdx]                       # Remove indices already selected from the available pool
    liCorrDomainBugsIdx = list()
    for(k in 1:iNumCorrRangeBugs){
      liCorrDomainBugsIdx[[k]] = sample(viAvailableIndices,viNumCorrDomainBugs[k],replace=FALSE)
      if(boolIndepCorr) viAvailableIndices[-liCorrDomainBugsIdx[[k]]]
    }

    # Converting vectors and lists to appropriately delimited strings
    strNumCorrDomainBugs = paste(viNumCorrDomainBugs[1])
    if(length(viNumCorrDomainBugs)>1){
      for(k in 2:length(viNumCorrDomainBugs)) strNumCorrDomainBugs <- paste(strNumCorrDomainBugs,viNumCorrDomainBugs[k],sep='; ')
    }

    strCorrRangeBugsIdx = paste(viCorrRangeBugsIdx[1])
    if(length(viCorrRangeBugsIdx)>1){
      for(k in 2:length(viCorrRangeBugsIdx)) strCorrRangeBugsIdx <- paste(strCorrRangeBugsIdx,viCorrRangeBugsIdx[k],sep='; ')
    }

    strCorrDomainBugsIdx = toString(liCorrDomainBugsIdx[[1]])
    if(length(liCorrDomainBugsIdx) > 1){
      for(k in 2:length(liCorrDomainBugsIdx[[1]])) strCorrDomainBugsIdx <- paste(strCorrDomainBugsIdx,toString(liCorrDomainBugsIdx[[k]]),sep='; ')
    }
  } 
  # Global parameters that should probably go outside
  c_strNoiseScaling = "Scaling parameter for variance of noise:"
  c_strMaxCorrRangeBugs = "Maximum number of bugs correlated with others:"
  c_strMaxCorrDomainBugs = "Maximum number of bugs with which one bug is correlated:"
  c_strCorrRangeBugs = "Number of bugs correlated with others:"
  c_strCorrDomainBugs = "Number of bugs each correlated bug is correlated with:"
  c_strCorrRangeBugsIdx = "Indices of bugs correlated with others:"
  c_strCorrDomainBugsIdx = "Indices of the bugs each correlated bug is correlated with:"
  c_strIndepCorr = "Whether the correlated bugs are independent of each other:"

  # This will hold the associations that you create and be placed in the truth file that is records association spik-ins for later acessment
  # It starts with the name of the microbiome you are creating
  # Parameters of interest and then your feature associations
  vStrParameters = c(paste(c_strSyntheticMicrobiome, c_strBugBugAssocations, sep='_'))
  vStrParameters = c(vStrParameters, paste(c_strNumberOfFeatures, int_number_features))
  vStrParameters = c(vStrParameters, paste(c_strNumberOfSamples, int_number_samples))
  vStrParameters = c(vStrParameters, paste(c_strNumberCounts, iMinNumberCounts))
  vStrParameters = c(vStrParameters, paste(c_strMinimumSamples, iMinNumberSamples))
  vStrParameters = c(vStrParameters, paste(c_strNoiseScaling, dVarScale))
  vStrParameters = c(vStrParameters, paste(c_strMaxCorrRangeBugs, iMaxNumberCorrRangeBugs))
  vStrParameters = c(vStrParameters, paste(c_strMaxCorrDomainBugs, iMaxNumberCorrDomainBugs))
  vStrParameters = c(vStrParameters, paste(c_strCorrRangeBugs,iNumCorrRangeBugs))
  vStrParameters = c(vStrParameters, paste(c_strCorrDomainBugs, strNumCorrDomainBugs))
  vStrParameters = c(vStrParameters, paste(c_strCorrRangeBugsIdx, strCorrRangeBugsIdx))
  vStrParameters = c(vStrParameters, paste(c_strCorrDomainBugsIdx, strCorrDomainBugsIdx))
  vStrParameters = c(vStrParameters, paste(c_strIndepCorr, toString(boolIndepCorr)))


  # Here is a method that will get you a mu_vector with default pushed through the script from the sfle call
  mu_vector = func_generate_mu_vector(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberSamples=iMinNumberSamples, iReadDepth=iReadDepth, vdMu=vdMu, vdSD=vdSD, vdPercentZero=dPercentZero, dBetaSD=dBetaSD, dSDIntercept=dSDIntercept, dBetaZero=dBetaZero)

  # This gets a lognormal distribution, you can do this with a predefined mu vector if needed.
  # After this step you should add the bug spike-ins
  # Note this is count data, normalization will happen automatically for you after you return the matrix
  mtrxBugs = func_generate_random_lognormal_matrix(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberCounts=iMinNumberCounts, iMinNumberSamples=iMinNumberSamples, iReadDepth=iReadDepth, vdMu=mu_vector[["mu"]], vdSD=mu_vector[['sd']], xPercentZero=dPercentZero, fZeroInflate=fZeroInflate, dBetaSD=dSDBeta, dSDIntercept=dSDIntercept, dBetaZero=dBetaZero)[["mat_bugs"]]

  # Generating the correlation strucutre (very simple at the moment)
  for(i in seq(1,iNumCorrRangeBugs)){
    vdVarCorrDomainBugs = apply(mtrxBugs,1,var)[liCorrDomainBugsIdx[[i]]]
    if(length(liCorrDomainBugsIdx[[i]])==1){
      mtrxBugs[viCorrRangeBugsIdx[i],] = mtrxBugs[liCorrDomainBugsIdx[[i]],]+
                                         rnorm(int_number_samples,mean=0,sd = sqrt(dVarScale*sum(vdVarCorrDomainBugs)))
    } else {
      mtrxBugs[viCorrRangeBugsIdx[i],] = apply(mtrxBugs[liCorrDomainBugsIdx[[i]],],2,sum)+
                                         rnorm(int_number_samples,mean=0,sd = sqrt(dVarScale*sum(vdVarCorrDomainBugs)))
    }
  }

  # This is the minimum return.
  return(list(mtrxBugs=mtrxBugs, vStrParameters=vStrParameters))
}

# Tested 3
func_zero_inflate <- function(
### Create a zero inflated log normal distribution with a specified mean and percentage of zeros.
### If you want to get the original values of mu and sd used in rlnorm, use mean(log(func_zero_inflate()))
### and sd(log(func_zero_inflate()))
dMean,
### Mean of the distribution (logged)
dPercentZeroInflated,
### Percentage of return which is zero
int_number_samples,
### The number of samples to create
dSD
### The sd of the distribution (logged)
){
  vdData = c()
  for( i in 1:int_number_samples)
  {
    fZeroInflated = sample(c(TRUE,FALSE), size=1, prob=c(dPercentZeroInflated, 1.0-dPercentZeroInflated), replace=TRUE)
    if(fZeroInflated)
    {
      vdData = c(vdData,0)
    } else {
      ### MUs and SD coming in should be from a logged rlnorm distribution
      vdData = c(vdData,rlnorm(1, dMean, dSD))
    }
  }

  #Truncate negatives to zero
  vdData[vdData<0]=0

#!#TODO should this be removed  # Remove extreme outliers
#  aiIndices = funcIdentifyNonZeroOutliers(vdData)
#  while(length(aiIndices)>0)
#  {
#    ### Mus and SD coming in should be from a logged rlnorm distribution
#    vdData[aiIndices] = rlnorm(length(aiIndices),dMean,dSD)
#    print("vdData")
#    print(vdData)
#    aiIndices = funcIdentifyNonZeroOutliers(vdData)
#  }
  return(vdData)
}

# TODO
nfunc_zero_inflate <- function(
### Create a zero inflated log normal distribution with a specified mean and percentage of zeros.
### If you want to get the original values of mu and sd used in rlnorm, use mean(log(func_zero_inflate()))
### and sd(log(func_zero_inflate()))
dMean,
### Mean of the distribution (logged)
dPercentZeroInflated,
### Percentage of return which is zero
int_number_samples,
### The number of samples to create
dSD
### The sd of the distribution (logged)
){
  vdData  = funcAdjSeqForZeros( rlnorm(int_number_samples, dMean, dSD), sample(1:int_number_samples, floor(int_number_samples*dPercentZeroInflated), replace=FALSE) ) 

  #Truncate negatives to zero
  vdData[vdData<0]=0

  return(vdData)
}

# Tested 2
#!# Replace with a numeric function
funcGetParamsForExpectation <- function(
### Searches for the best Mu and SD combination to match the Expectation
### for the data set.
dMuInit,
### Initial starting guess for the Mu
dReadDepthTarget,
### Target read depth of interest
dSampleCount,
### Describes the relationship between the Mu and SD (SD will be calculated). This is the coefficient of the relationship.
dFeatureBetaSD,
### Describes the relationship between the Mu and SD (SD will be calculated). This is the intercept of the reationship.
dFeatureInterceptSD,
### Number of features which will be generated
dError = .1,
### Minimal cost in the difference between the predicted ReadDepth and the target  read depth
iIterations = 10000,
### Maximal number of iterations before returning an answer even if the cost is less than the dError
dIncrement = .1
### The increment amount for searching for the optimal Mu
){
  # Initialize best to the first guess
  # Can not have a negative mu or SD
  dBestMu = max(dMuInit,0+dIncrement)
  dBestSD = funcEstimateFeatureSD(dBestMu,dFeatureBetaSD,dFeatureInterceptSD)
#  dBestSD = max(funcEstimateFeatureSD(dBestMu,dFeatureBetaSD,dFeatureInterceptSD),0+dIncrement)
  dBestDepth = funcGetExp(dBestMu,dBestSD)
  dTarget = dReadDepthTarget

  # Will search forwards and backwards for a mu and sd pair that meet the criteria for the read depth
  dMuForwards = dMuInit
  dMuBackwards = dMuInit

  # Will be flagged True if the interations run out without satistfying the dError criteria
  fNoConverge = FALSE

  # Iterate modifying the Mu less/more by the increment to find the optimal solution.
  for(iIndex in 1:iIterations)
  {
    dMuForwards = dMuForwards + dIncrement
    dMuBackwards = max(dMuInit - dIncrement,0+dIncrement)

    dCurSD = funcEstimateFeatureSD(dMuForwards,dFeatureBetaSD,dFeatureInterceptSD)
#    dCurSD = max(funcEstimateFeatureSD(dMuForwards,dFeatureBetaSD,dFeatureInterceptSD),0+dIncrement)
    dCurDepth = funcGetExp(dMuForwards, dCurSD)

    if(abs(dBestDepth-dTarget)>abs(dCurDepth-dTarget))
    {
      dBestMu = dMuForwards
      dBestDepth = dCurDepth
      dBestSD = dCurSD
    }

    if( 0 < (dMuBackwards) )
    {
      dCurSD = funcEstimateFeatureSD(dMuBackwards,dFeatureBetaSD,dFeatureInterceptSD)
#      dCurSD = max(funcEstimateFeatureSD(dMuBackwards,dFeatureBetaSD,dFeatureInterceptSD),0+dIncrement)
      dCurDepth = funcGetExp(dMuBackwards, dCurSD)

      if(abs(dBestDepth-dTarget)>abs(dCurDepth-dTarget))
      {
        dBestMu = dMuBackwards
        dBestDepth = dCurDepth
        dBestSD = dCurSD
      }
    }
    if(abs(dBestDepth-dTarget) < dError){break}

    fNoConverge = iIndex==iIterations
  }
  if(fNoConverge)
  {
    print(paste("funcGetParamsForReadDepth: Did not meet the optimization requirement (dError). Returning a potentially suboptimal result. This should give a read depth on average of ", dBestDepth, "but was targeting", dReadDepthTarget, ". Initial Mu=", dMuInit,"Sample Count=", dSampleCount, "SD Beta=", dFeatureBetaSD, "SD Intercept=",dFeatureInterceptSD))
  }
  return(list(dBestMu=dBestMu, dBestSD=dBestSD, dBestDepth=dBestDepth, dTargetDepth=dTarget))
}

# Tested 2
#!# Replace with a number function call
funcGetParamsForReadDepth <- function(
### Searches for the best Mu and SD combination to match the ReadDepth
### for the data set.
dMuInit,
### Initial starting guess for the Mu
dReadDepthTarget,
### Target read depth of interest
dFeatureCount,
### Describes the relationship between the Mu and SD (SD will be calculated). This is the coefficient of the relationship.
dGrandBetaSD,
### Describes the relationship between the Mu and SD (SD will be calculated). This is the intercept of the reationship.
dGrandInterceptSD,
### Number of features which will be generated
dError = .5,
### Minimal cost in the difference between the predicted ReadDepth and the target  read depth
iIterations = 10000,
### Maximal number of iterations before returning an answer even if the cost is less than the dError
dIncrement = .1
### The increment amount for searching for the optimal Mu
){
  # Initialize best to the first guess
  dBestMu = dMuInit
  dBestSD = funcEstimateGrandSD(dBestMu,dGrandBetaSD,dGrandInterceptSD)
  dBestDepth = funcGetExp(dBestMu,dBestSD)
  dTarget = dReadDepthTarget

  # Will search forwards and backwards for a mu and sd pair that meet the criteria for the read depth
  dMuForwards = dMuInit
  dMuBackwards = dMuInit

  # Will be flagged True if the interations run out without satistfying the dError criteria
  fNoConverge = FALSE

  # Iterate modifying the Mu less/more by the increment to find the optimal solution.
  for(iIndex in 1:iIterations)
  {
    dMuForwards = dMuForwards + dIncrement
    dMuBackwards = dMuInit - dIncrement

    dCurSD = funcEstimateGrandSD(dMuForwards,dGrandBetaSD,dGrandInterceptSD)
    dCurDepth = funcGetExp(dMuForwards, dCurSD)

    if(abs(dBestDepth-dTarget)>abs(dCurDepth-dTarget))
    {
      dBestMu = dMuForwards
      dBestDepth = dCurDepth
      dBestSD = dCurSD
    }

    if( 0 < (dMuBackwards) )
    {
      dCurSD = funcEstimateGrandSD(dMuBackwards,dGrandBetaSD,dGrandInterceptSD)
      dCurDepth = funcGetExp(dMuBackwards, dCurSD)

      if(abs(dBestDepth-dTarget)>abs(dCurDepth-dTarget))
      {
        dBestMu = dMuBackwards
        dBestDepth = dCurDepth
        dBestSD = dCurSD
      }
    }
    if(abs(dBestDepth-dTarget) < dError){break}

    fNoConverge = iIndex==iIterations
  }
  if(fNoConverge)
  {
    print(paste("funcGetParamsForReadDepth: Did not meet the optimization requirement (dError). Returning a potentially suboptimal result. This should give a read depth on average of ", dBestDepth*dFeatureCount, "but was targeting", dReadDepthTarget, "."))
  }
  return(list(dBestMu=dBestMu, dBestSD=dBestSD, dBestDepth=dBestDepth, dTargetDepth=dTarget))
}

# Reviewed but not tested
func_generate_mu_vector <- function(
### Generate the initial mu vector and associated SD and percent zeros if needed
# If a mu vector is given of the size of the int_number_samples then pass through
# Otherwise sample to that size with replacement.
int_number_features,
### Number of features
int_number_samples,
### Number of samples
iMinNumberSamples,
### Minimum number of samples
iReadDepth,
### Simulated read depth
vdMu = NA,
### Vector of Mu for the original mu distribution (means of features) if not supplied, one will be generated by rlnorm
vdSD = NA,
### The vector of SD matching the vdMus.
vdPercentZero = NA,
### The vector of percent zeros matching the vdMus
dBetaSD = c_dBetaSD,
### If vdSD is not given, SD will be generated by a relationship with Mu using this relationship between Mu and SD
dSDIntercept = c_dSDIntercept,
### If vdSD is not given, SD will be generated by a relationship with Mu using the intercept given here
dBetaZero = c_dBetaZero,
### If vdPercentZero is not given, vdPercentZero will be generated by a relationship with Mu using this relationship between Mu and percent zero with an intercept
dInterceptZero = c_dZeroIntercept,
### If vdPercentZero is not given, vdPercentZero will be generated by a relationship with Mu using this relationship between Mu and percent zero with an intercept
dGrandMu = c_dGrandMu,
### The guess for the starting point for the correct mu of mus value
dBetaGrandSD = c_dBetaGrandSD,
### The beta for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
dInterceptGrandSD = c_dInterceptGrandSD,
### The beta for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
fVerbose = FALSE
### Controls the plotting of graphic pdf logging (Default FALSE, TRUE indicates logging occurs)
){

  if(fVerbose & !is.na(vdMu))
  {
    hist(vdMu,main="func_generate_mu_vector: Original mu vector")
    hist(rlnorm(int_number_features, log(dGrandMu), log(dGrandMu)*dBetaGrandSD+dInterceptGrandSD), main=paste("func_generate_mu_vector: Predicting Original mu vector","dGrandMu", dGrandMu,"dBetaGrandSD", dBetaGrandSD))
  }

  if(is.na(vdMu))
  {
    print("func_generate_mu_vector: Generate vdMu Vector.")
    # Draw a vector of mu's
    # This will be the template distribution all bugs within a sample will be based on.
    # This allows to have structure within the data set of bugs more or less prevalent with a level of consistency
    # The mean of the distribution of bugs is derived from the max number of bug counts for a sample, making the
    # total number of bugs per sample the same within a small random margin
    lsParams = funcGetParamsForReadDepth(dMuInit=dGrandMu, dReadDepthTarget=iReadDepth/int_number_features, dFeatureCount=int_number_features, dGrandBetaSD=dBetaGrandSD, dGrandInterceptSD=dInterceptGrandSD)

    vdMu <- log(rlnorm(int_number_features, log(lsParams$dBestMu), log(lsParams$dBestSD)))
    # TODO Why?
    vdMu[vdMu<0]=0.00001
    if(fVerbose)
    {
      hist(vdMu,main=paste("Generated Mu vector. Actual:",round(mean(vdMu),2),"Real:",round(iReadDepth/int_number_features,2)))
      barplot(vdMu, main=paste("Generated actual average read depth",mean(vdMu)), xlab="Feature", ylab="log(Mu)")
      abline(mean(vdMu),0,col="violet")
    }
    # Since new mus have been created and SD and percent zero depend on mus, these are reset to zero so they will be regenerated.
    vdSD = NA
    vdPercentZero = NA
  }
  if(is.na(vdSD))
  {
    # Generate vector of SD based on mu since it is not known
    print("func_generate_mu_vector: Generate vdSD Vector.")
    vdSD = sapply(vdMu, function(x) funcEstimateFeatureSD(x,dBetaSD,dSDIntercept))
    #TODO why?
    vdSD[vdSD<0] = 0.00001
#    vdSD = sapply(vdMu, function(x) max(funcEstimateFeatureSD(x,dBetaSD,dSDIntercept), c_dMinimumSD))
    if(fVerbose){plot(vdMu, vdSD, main="Generated Relationship of Mu and SD", col="orange")}
  }

  # TODO not super excited about having this min here but it stays until the logging in the functions do not throw them off.
  # Min mu and sd requirements
  # Otherwise the logging of the parameter gives and error
#  viLowMu = union(which(vdMu<=1), which(vdSD<=1))
#  dLowestMu = funcEstimateMinimumFeatureSD(dBetaSD,dSDIntercept)
#  vdMu[viLowMu] = dLowestMu
#  vdSD[viLowMu] = 0
#  vdSD[viLowMu] = funcEstimateFeatureSD(dLowestMu, dBetaSD, dSDIntercept)

  if(is.na(vdPercentZero))
  {
    # Generate vector of percent zero based on mu since it is not known
    print("func_generate_mu_vector: Generate vdPercentZero Vector.")
    vdPercentZero = funcEstimatePercentZero(vdMu,dBetaZero,dInterceptZero)
    if(fVerbose){plot(exp(1)^vdMu, vdPercentZero, main="Generated Relationship of PercentZero and Mu", col="purple")}
  }
  if(!length(vdMu)==int_number_features)
  {
    # This is the scenario that the calibration file is used and the number of the samples needed are not equal
    # This correct number of are selected with replacement.
    print("func_generate_mu_vector: Reseting count of samples.")
    viWhich = sample(1:length(vdMu), size = int_number_features, replace = TRUE)
    vdMu = vdMu[viWhich]
    vdSD = vdSD[viWhich]
    vdPercentZero = vdPercentZero[viWhich]   
  }

  # QC and contraints
  #!# check min and max?
  # Make sure the percent zero passes the min
  # If there are not enough nonzeors, there is no signal to use.
  dMinPercent = iMinNumberSamples/int_number_samples
#  dMaxPercent = 1-dMinPercent
  print("dMinPercent")
  print(dMinPercent)
#  print("dMaxPercent")
#  print(dMaxPercent)
#  vdPercentZero[which(vdPercentZero>dMaxPercent)] = dMaxPercent
  vdPercentZero[which(vdPercentZero<dMinPercent)] = dMinPercent

  return(list(mu=vdMu,sd=vdSD,PercentZero=vdPercentZero))
}


funcMakeFeature <- function(
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
dTotalReadDepth,
### The max count a feature can have.
dMaxCount,
### Max count for an entry in the feature
fVerbose
### If pdf logging of feaure should occur
){
  # Generate feature
  dFeature = func_zero_inflate(log(dMu), dPercentZero, iNumberSamples, log(dSD))

  # Change out extreme outliers
  iOutliersLocations = which(dFeature>dMaxCount)
  iOutliers = dFeature[iOutliersLocations]
  if(length(iOutliers))
  {
    for(iIndexOutlier in 1:length(iOutliers))
    {
      dNewEntry = func_zero_inflate(log(dMu),0,1,log(dSD))
      iCount = 0 
      while((dNewEntry>dMaxCount)|iCount<1000)
      {
        dNewEntry = func_zero_inflate(log(dMu),0,1,log(dSD))
        iCount = iCount + 1
      }
      # Get the difference of the outlier and the replacement
      # and shuffle throughout the feature to try to keep average read depth
      # The extra shuffle is not added into back uniformly becaus that will create a band of values
      # or a large space between the 0 and other value entries. This makes more of a smooth
      # distribution.
      dDifference = dFeature[iOutliersLocations[iIndexOutlier]] - dNewEntry
      dFeature[iOutliersLocations[iIndexOutlier]] = dNewEntry
      if(dDifference>0)
      {
        viScatterLocations = sample(which(dFeature>0),dDifference,replace=TRUE)
        for(iScatterIndex in 1:length(viScatterLocations))
        {
          dFeature[viScatterLocations[iScatterIndex]] = dFeature[viScatterLocations[iScatterIndex]] + min(iScatterIndex,dDifference)
          dDifference = max(0,dDifference - iScatterIndex)
        }
      }
    }
  }
  # Min is placed here because if you ask for 2 samples but this is evaluating to find atleast 3 samples at a certain level you will inf loop
  # Look to see how many measurements are above the min.
  # If there are not enough counts more the the min number of counts
  # then make sure there are by drawing more samples (not zero inflated)
  # And changing any measurement less than the min to the min
  # Then randomly assigning this to a location that is zero
#  vdSignal = which(dFeature >= iMinNumberCounts)
#  if(length(vdSignal) < iMinNumberSamples)
#  {
#    viNonSignalLocations = which(dFeature<iMinNumberCounts)
#    print("which(dFeature<iMinNumberCounts)")
#    print(which(dFeature<iMinNumberCounts))
#    print("which(dFeature==0)")
#    print(which(dFeature==0))
#    print("viNonSignalLocations")
#    print(viNonSignalLocations)
#    # Randomly select locations to add additional signal to
#    vdNewSignalLocations = sample(x=viNonSignalLocations, iMinNumberSamples - length(vdSignal), replace=FALSE)
#
#    # Replace the feature with these nonzero measurements
#    print("vdSignal")
#    print(vdSignal)
#    dFeature[vdNewSignalLocations] = dFeature [sample(vdSignal,length(vdNewSignalLocations), replace=TRUE)] + iMinNumberCounts
#    print("dFeature")
#    print(dFeature)
#  }

  # If any count is bigger than the dTotalReadDepth it is reset to dTotalReadDepth
  dFeature[dFeature> dTotalReadDepth] = dTotalReadDepth

  # Extra useful measurements, the true and expected means
  dMean = mean(dFeature)
  dExpCal = funcGetExp(dMu,dSD)
  dExp = mean(dMu,dSD)

  if(fVerbose)
  {
    plot( dFeature, main = paste("funcMakeFeature","Mean",round(dMean,2),"SD",round(sd(log(dFeature)),2),"Exp",round(dExp,2)))
  }
  return(list(Feature = dFeature, Mean = dMean, Exp = dExp, ExpCal = dExpCal))
}

func_generate_random_lognormal_matrix <- function(
int_number_features,
### Number of features
int_number_samples,
### Number of samples,
iMinNumberCounts,
### Minimum number of counts for a feature to be considered in a smaple for the QC filtering
iMinNumberSamples,
### Created bugs must have a minimum number of samples (iMinNumberSamples) that have a minimum number of counts (iMinNumberSamples)
iReadDepth,
### Simulated read depth for sample creation
vdMu = NA,
### Vector of Mu for the original mu distribution (means of features) if not supplied, one will be generated by rlnorm
xPercentZero = c_dBetaZero,
### Either a vector of percent zero for each mu or a single parameter modeling the relationship between percent zeros and mus
vdSD = NA,
fZeroInflate = TRUE,
### Turns off Zero inflation if FALSE (default TRUE, zero inflation turned on)
dBetaSD = c_dBetaSD,
dSDIntercept = c_dInterceptSD,
dBetaZero = c_dZeroBeta,
dInterceptZero = c_dZeroIntercept,
dGrandMu = c_dGrandMu,
dBetaGrandSD = c_dBetaGrandSD,
dInterceptGrandSD = c_dInterceptGrandSD,
dMaxCount = c_dMaxCount,
fVerbose = FALSE
){
  print(paste("func_generate_random_lognormal_matrix::","dGrandMu", dGrandMu, "dGrandSD", funcEstimateGrandSD(dGrandMu, dBetaGrandSD, dInterceptGrandSD), "int_number_samples", int_number_samples, "int_number_features", int_number_features, "iReadDepth", iReadDepth, "dBetaSD", dBetaSD, "dSDIntercept", dSDIntercept, "fZeroInflate", fZeroInflate))

  # Preallocating for speed
  mat_bugs <- matrix(data=NA,nrow=int_number_features,ncol=int_number_samples)
  # Get the initial mu vector for generating features.
  lsInitialDistribution = func_generate_mu_vector(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberSamples=iMinNumberSamples, iReadDepth=iReadDepth, vdMu=vdMu, vdSD=vdSD, vdPercentZero=xPercentZero, dBetaSD=dBetaSD, dSDIntercept=dSDIntercept, dBetaZero=dBetaZero, dInterceptZero=dInterceptZero, dGrandMu = dGrandMu, dBetaGrandSD = dBetaGrandSD, dInterceptGrandSD = dInterceptGrandSD, fVerbose=fVerbose)

  print("lsInitialDistribution")
  print(lsInitialDistribution)

  # Need to set the distributions to the right maginitude, previously in the log form
  lsInitialDistribution$mu = exp(1)^lsInitialDistribution$mu
  lsInitialDistribution$sd = exp(1)^lsInitialDistribution$sd
  #!# Are the mu and the sd in the right range?

  # Update the Mu, SD and Percent zero bugs and report on distributions
  mu_vector = lsInitialDistribution[["mu"]]
  vdSD = lsInitialDistribution[["sd"]]
  xPercentZero = lsInitialDistribution[["PercentZero"]]
  if(fVerbose)
  {
    # Plot Distribution of mu vector
    hist(mu_vector, main=paste("Generated Initial Mus: Base Log Normal Matrix. RD=",round(sum(mu_vector)),"AveRD=",round(mean(mu_vector)),sep=" "))
    plot(mu_vector, main=paste("Generated Initial Mus: Base Log Normal Matrix. RD=",round(sum(mu_vector)),"AveRD=",round(mean(mu_vector)),sep=" "))
  }

  # Number of samples needed to have signal as a constraint
  iNumberSamples = min(int_number_samples, iMinNumberSamples)

  # Count the read depth left over
  dTotalReadDepth = iReadDepth

  # Measure how well the features were generated
  vdActual = c()
  vdPredicted = c()

  print("start bug generation")
  # Generate the bugs based on the previously defined distribution
  for (i in 1:int_number_features)
  {
    # If there is no more read depth for the new features
    # Take the most extreme feature and replace it with an estimate from a
    # mu / sd combination from the lower two quartiles of the mus
    print("dTotalReadDepth")
    print(dTotalReadDepth)
    if(dTotalReadDepth<0)
    {
      # Select feature of highest magnitude and reset.
      # Remove it's contribution to read depth
      iReset = which(vdActual==max(vdActual))[1]
      dTotalReadDepth = dTotalReadDepth+vdActual[iReset]
      print(paste("Removed sample ",iReset))
      # Select a new mu from values less than the median of what is already present
      # Get associated SD and percent zeros
      iNewMu = sample(which(vdActual<median(vdActual)),1)
      dCurMu = mu_vector[iNewMu]
      dCurSD = vdSD[iNewMu]
      dZeroInflate = 0
      if(fZeroInflate){dZeroInflate = xPercentZero[iNewMu]}

      # Create new feature
      lFeatureDetails = funcMakeFeature(dMu=dCurMu, dSD=dCurSD, dPercentZero=dZeroInflate, iNumberSamples=int_number_samples, iMinNumberCounts=iMinNumberCounts, iMinNumberSamples=iMinNumberSamples, dTotalReadDepth=dTotalReadDepth, dMaxCount = dMaxCount, fVerbose=fVerbose )
      print("lFeatureDetails")
      print(lFeatureDetails)

      # Update the matrix with the new feature
      mat_bugs[iReset,] = lFeatureDetails[["Feature"]]

      # Update the variables monitoring this step (read depth, actual mu, and predicted mu)
      dTotalReadDepth = dTotalReadDepth-lFeatureDetails[["Mean"]]
      vdActual[iReset] = lFeatureDetails[["Mean"]]
      vdPredicted[iReset] = lFeatureDetails[["Exp"]]

      # Update the global parameters for the features
      mu_vector[iReset] = mu_vector[iNewMu]
      vdSD[iReset] = vdSD[iNewMu]
      xPercentZero[iReset] = xPercentZero[iNewMu]
    }

    # Current parameters
    dCurMu = mu_vector[i]
    dCurSD = vdSD[i]
    # Turn off zero inflation if need be
    dZeroInflate = 0
    if(fZeroInflate){dZeroInflate = xPercentZero[i]}

    # Get zero-inflated feature
    lFeatureDetails = funcMakeFeature(dMu=dCurMu, dSD=dCurSD, dPercentZero=dZeroInflate, iNumberSamples=int_number_samples, iMinNumberCounts=iMinNumberCounts, iMinNumberSamples=iMinNumberSamples, dTotalReadDepth=dTotalReadDepth, dMaxCount = dMaxCount, fVerbose=fVerbose )

    # Update matrix
    mat_bugs[i,] <- lFeatureDetails[["Feature"]]

    # Update the read depth, predicted and actual expectations
    dTotalReadDepth = dTotalReadDepth-lFeatureDetails[["Mean"]]
    vdActual = c(vdActual,lFeatureDetails[["Mean"]])
    vdPredicted = c(vdPredicted,lFeatureDetails[["Exp"]])
  }
  print("Left over Read Depth")
  print(dTotalReadDepth)

  # Report on feature creation
  if(fVerbose)
  {
    plot(sapply(1:ncol(mat_bugs), function(x) sum(mat_bugs[,x])), main=paste("Read depth of Microbiome","( Exp ", iReadDepth,")") )
    abline(mean(sapply(1:ncol(mat_bugs), function(x) sum(mat_bugs[,x]))),0,col="red")
    abline(iReadDepth,0,col="grey")
  }

  ### Plot zero inflated model checks
  vdMeans = c()
  viZeroCounts = c()
  vdSD = c()

  for(iIndex in 1:dim(mat_bugs)[1])
  {
    vColumn = mat_bugs[iIndex,]
    ## Measure mean of each column
    vdMeans = c(vdMeans,mean(vColumn[which(vColumn!=0)]))

    ## Measure number of 0s in each column
    viZeroCounts = c(viZeroCounts, length(which(vColumn==0)))

    ## Measure SD in each column
    vdSD = c(vdSD,sd(vColumn[which(vColumn!=0)]))
  }
  viZeroCounts = viZeroCounts/length(viZeroCounts)
  vdSD[is.na(vdSD)]=0

  # Floor to Counts
  mat_bugs = floor(mat_bugs)

  ## Plot
  dActualSampleAverage = mean(colSums(mat_bugs))
  if(fVerbose)
  {
    barplot(colSums(mat_bugs),main=paste("Total Read Depth per Sample (Target:",round(iReadDepth),"Actual",round(dActualSampleAverage),")"), xlab="Sample", ylab="Read Depth")
    abline(iReadDepth,0,col="grey")
    abline(dActualSampleAverage,0,col="red")
    plot(vdMeans,viZeroCounts, main="GenRLNORM: Zero Occurence by Mean of Non-zero Measurements", xlab="Means (excluding 0s)", ylab="Zero Occurence")
    plot(vdMeans, vdSD, main="GenRLNORM: Standard Deviation by Mean of Non-zero Measurements", xlab="Means (excluding 0s)", ylab="Standard Deviation (excluding 0s)")
    if(sum(viZeroCounts>0)){
      plot(log(vdMeans,2),log(viZeroCounts,2), main="GenRLNORM: Log Zero Occurence by Mean of Non-zero Measurements", xlab=paste("Log Means (excluding 0s)",length(vdMeans),"Samples"), ylab="Log Zero Occurence")
    }
    plot(log(vdMeans), log(vdSD), main="GenRLNORM: Log Standard Deviation by Mean of Non-zero Measurements", xlab=paste("Log Means (excluding 0s)",length(vdMeans),"Samples"), ylab="Log Standard Deviation (excluding 0s)")
    ### END Plot zero inflated model checks
    hist(as.vector(mat_bugs), main="GenRLNORM: mat_bugs finish")
    hist(log(as.vector(mat_bugs)), main="GenRLNORM: logged final mat_bugs")
  }

  # Truth table for log normal data
  mtrxParameters <- matrix(data=NA, nrow=6, ncol=1)
  mtrxParameters[1,1] <- paste(c_strSyntheticMicrobiome, c_strRandom, sep='')
  mtrxParameters[2,1] <- paste(c_strNumberOfFeatures, int_number_features)
  mtrxParameters[3,1] <- paste(c_strNumberOfSamples, int_number_samples)
  mtrxParameters[4,1] <- paste(c_strTotalSampleBugOccurrence, iReadDepth)
  mtrxParameters[5,1] <- paste(c_strNumberCounts, iMinNumberCounts)
  mtrxParameters[6,1] <- paste(c_strNumberSamples, iMinNumberSamples)

  print("end func_generate_random_lognormal_matrix")
  return(list(mu_vector=mu_vector, mat_bugs=mat_bugs, mtrxParameters=mtrxParameters))
  ### Returns a row major matrix of log-normal data.
}

# Tested
funcIsFactorMetadataValid = function(
### Check to make sure a level of entrophy is in a discontinuous metadata so it has a minimal level to associate
vxMetadata,
### A metadata that can be factor data
iMin
### Minimum number of instances of a level
){
  vxFactorMetadata = as.factor(vxMetadata)
  lstrLevels = levels(vxFactorMetadata)
  for(strLevel in lstrLevels)
  {
    if(length(which(vxFactorMetadata==strLevel))<iMin)
    {
      return(FALSE)
    }
  }
  return(TRUE)
}

#TODO look over
func_generate_metadata <- function(
int_base_metadata_number,
int_number_samples,
### Number of samples
dMinLevelPercent
### The minimum percent of samples a level can have
){
  # Preallocate matrix
  mat_metadata <- matrix(data=NA,nrow=(int_base_metadata_number*c_iCountTypesOfMetadata),ncol=int_number_samples)

  # Used to report on metadata
  mtrxParameters = c(c_strMetadataDetails)

  # Continous metadata, generated means
  # Generating and padding the list of potential mean values
  li_mean_value_list = list(runif(1, c_dRunifMin, c_dRunifMax),1,100)
  if(length(li_mean_value_list) < (int_base_metadata_number*2))
  {
    for (k in (length(li_mean_value_list)+1):(int_base_metadata_number*2))
    {
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
    mat_metadata[i,] <- rnorm(int_number_samples,mean=mean_value,sd=mean_value/5)
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
      probability = runif(1, c_dRunifMin, c_dRunifMax)
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
    mat_metadata[i,] <- vsCurMetadata
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
    mat_metadata[i,] <- vsCurMetadata

    mtrxParameters = c(mtrxParameters, paste(c_strMetadata, i, " ", c_strFactor, " ", paste(levels(as.factor(mat_metadata[i,])), collapse = " "), sep=""))
    i = i+1
  }
  return(list(mat_metadata=mat_metadata, mtrxParameters=mtrxParameters))
}

funcBoxPlotOutliers = function(
### Generic plotter for a series of boxplots.
mtrxData,
strTitle,
lviHighlight = list(),
strColor = "red",
fWithZeros = TRUE,
fBySamples = TRUE
### Samples are expected to be columns
){
  strXLabel="Samples"
  if(!fBySamples)
  {
    strXLabel = "Features"
    mtrxData = t(mtrxData)
  }
  if(fWithZeros)
  {
    boxplot(mtrxData, main=paste(strTitle," - With Zeros - By ",strXLabel, sep=""), xlab = strXLabel, ylab="Distribution of Counts")
    lvdData = list()
    viPosition = c()
    for(iCol in 1:ncol(mtrxData))
    {
      lvdData[[length(lvdData)+1]] = mtrxData[,iCol]
      viPosition = c(viPosition,length(lvdData))
      lstrColoring = rep(c_strDefaultMarkerColor,length(mtrxData[,iCol]))

      if(!is.null(lviHighlight[iCol][[1]]))
      {
        lstrColoring[lviHighlight[[iCol]]] = strColor
      }

      # Highlight for the different levels (highlight and not)
      for(strLevel in levels(as.factor(lstrColoring)))
      {
        stripchart(mtrxData[,iCol][which(lstrColoring==strLevel)],at=iCol,add=TRUE,vertical=TRUE,method="jitter",jitter=0.3,col=strLevel)
      }
    }
 } else {
    liNoZeros = list()
    for(iColumn in 1:ncol(mtrxData))
    {
      viNotZeros = which(mtrxData[,iColumn]!=0)
      liNoZeros[[length(liNoZeros)+1]] = mtrxData[,iColumn][viNotZeros]

      fTrimHighlights = !is.null(lviHighlight[iColumn][[1]])
      # Trim down colors to the same indices as non zeros
      if(fTrimHighlights)
      {
        lviHighlight[[iColumn]] = intersect(lviHighlight[[iColumn]],viNotZeros)
        # Shift the indices down to account for removing zeros
        viShiftedIndices = c()
        for(iIndex in lviHighlight[[iColumn]])
        {
          viShiftedIndices = c(viShiftedIndices,1+length(which(viNotZeros<iIndex)))
        }
        if(is.null(viShiftedIndices))
        {
          lviHighlight[[iColumn]][1] = NULL
        } else {
          lviHighlight[[iColumn]] = viShiftedIndices
        }
      }
    }
    boxplot(liNoZeros, main=paste(strTitle," - Without Zeros - By ", strXLabel, sep=""), xlab = strXLabel, ylab="Distribution of Counts")
    for(iColumn in 1:ncol(mtrxData))
    {
      liNoZeros[[length(liNoZeros)+1]] = mtrxData[,iColumn][which(mtrxData[,iColumn]!=0)]
      lstrColoring = rep(c_strDefaultMarkerColor,length(mtrxData[,iColumn]))
      if(!is.null(lviHighlight[iColumn][[1]]))
      {
        lstrColoring[lviHighlight[[iColumn]]] = strColor
      }

      # Highlight for the different levels (highlight and not)
      for(strLevel in levels(as.factor(lstrColoring)))
      {
        stripchart(liNoZeros[[iColumn]][which(lstrColoring==strLevel)],at=iColumn,add=TRUE,vertical=TRUE,method="jitter",jitter=0.3,col=strLevel)
      }
    }
  }
}

func_generate_random_lognormal_with_outliers <- function(
### Generates a random log normal distribution of data as a null matrix
### The option of using a matrix passed in as a parameter as the null matrix is provided (mtrxBugs)
### A percent of samples are given outliers based on the percent parameter
int_number_features,
### Number of features
int_number_samples,
### Number of samples
iMinNumberCounts,
### Minimum number of counts for a feature to be considered in a smaple for the QC filtering
iMinNumberSamples,
### Created bugs must have a minimum number of samples (iMinNumberSamples) that have a minimum number of counts (iMinNumberSamples)
dMaxPercentOutliers,
### The maximum percent of outliers to create in each sample (0 =< dPercent =< 1.0)
dPercentSamples,
### Percent of samples given outliers (0 =< dPercent =< 1.0)
mtrxBugs,
### Precalcualted null matrix
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
        iMaxValue <- dfSorted[(int_number_features-(iOutlier-1)), iSample]
        iMinValue <- dfSorted[iOutlier+iBufferForZeros, iSample]
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

        mtrxBugs[c(iRowMin, iRowMax), iSample] <- mtrxBugs[c(iRowMax, iRowMin), iSample]
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
  for(iIndex in 1:length(lviSwapped))
  {
    for(iItemSwapped in lviSwapped[[iIndex]])
    {
      mtrxParameters = c(mtrxParameters, paste(c_strOutlierParameter,paste(c_strFeature,c_strOutlier, iItemSwapped, sep="_"),c_strSampleParameter,iIndex))
    }
  }

  print("Stop func_generate_random_lognormal_with_outliers")

  # And return
  return(list(mat_bugs=mtrxBugs, mtrxParameters=mtrxParameters))
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
  metadata_average <- funcGetRowMean(vdCurMetadata)
  metadata_sigma = funcGetSD(vdCurMetadata)

  # Get the average and SD of the feature (ignoring zeros if zero inflated)
  data_average <- mean(vdCurData[which(vdCurData>0)])
  data_sigma <- sd(vdCurData[which(vdCurData>0)])
  if(is.na(data_sigma)){return(NULL)}

  if(!fZeroInflated)
  {
    data_average <- mean(vdCurData)
    data_sigma <- sd(vdCurData)
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
  scaled_metadata <- (vdCurMetadata - metadata_average)/metadata_sigma
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
  if(is.null(vdSpikedBug)){return(list(PASS=FALSE, CommonCOunts=c()))}

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

func_generate_random_lognormal_with_multivariate_spikes <- function(
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
  # Tracks the bug of interest
  iIndexSpikedFeature = NA

  # Initialize froze levels if need be
  if(is.null(lsFrozeLevels))
  {
    lsFrozeLevels = list()
  }

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

  # For each spiked in bug
  for(iSpikedBug in 1:floor(int_number_features*percent_spikes))
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
    lxCurrentBestRun = list(Metadata = c(), MetadataNames = c(), SpikinBug = c(), BugIndex = 0, Count = -1)

    # Metadata indices that are selected for spikin
    viSelectedMetadata = c()

    # Find valid spike-in scenario or best choice
    while(fSpikeInFailed)
    {
      # Get the bug to attempt the association
      if(!is.null(liFrozeDataIndicies))
      {
        iIndexSpikedFeature <- liFrozeDataIndicies[[iSpikedBug]]
        iSpikeInLoopControl <- iLoopingControlIncrement + 1
      } else {
        iIndexSpikedFeature <- sample(viRemainingFeatures,1)
      }

      # Select which of the metadatum we will be using to scale
      if(!is.null(lviFrozeMetadataIndices))
      {
        viSelectedMetadata <- lviFrozeMetadataIndices[[iSpikedBug]]
      } else {
        viSelectedMetadata <- sample(1:nrow(metadata_matrix), multivariate_parameter, replace=TRUE)
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
        if(sum(lxQCInfo[["CommonCounts"]]>iMinSamples) > lxCurrentBestRun[["Count"]])
        {
          lxCurrentBestRun = list(Metadata = vdCurMetadata, MetadataNames = vstrSpikedMetadata, SpikinBug = vdSpikedBug, BugIndex = iIndexSpikedFeature,  Count = sum(lxQCInfo[["CommonCounts"]]>iMinSamples), MetadataIndices = viSelectedMetadata, Levels = vsCurFrozeLevels)
        }
        
        # If we have ran out of iteration, use the best failed scenario and indicate this.
        if(iSpikeInLoopControl > iLoopingControlIncrement)
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
      }
    }
    # Update metadata and data indices for the features used for spikins
    lviMetadata[[ length(lviMetadata)+1 ]] = viSelectedMetadata
    liData[[ length(liData)+1 ]] = iIndexSpikedFeature

    # Successful Spike-in so update
    mtrxBugs[iIndexSpikedFeature,] <- vdSpikedBug
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
  # Floor to counts (spike-ins will be real numbers)
  mtrxBugs <- floor(mtrxBugs)

  # Return normalized matrix and spike-in list
  return(list(mat_bugs=mtrxBugs, m_parameter_rec=m_parameter_rec, MetadataIndices=lviMetadata, DataIndices=liData, Levels=lsFrozeLevels))
}

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

# Starting the actual script
option_list = list(
make_option(c("-e","--read_depth"), type="integer", default=8030, help="Simulated read depth for counts"),
make_option(c("-s","--seed"), type="integer", default=NA, help="random seed"),
make_option(c("-a","--lefse_file"), type="character", default=NULL, help="folder for lefse inputs"),
make_option(c("-p","--number_metadata"), type="integer", default=5, dest='number_metadata',help="base number of metadata -- continuous = base*2, binary = base, quarternary = base"),
make_option(c("-f","--number_features"), type="integer", default=300, help="number of features"),
make_option(c("-n","--number_samples"), type="integer", default=50, help="number of samples"),
make_option(c("-o","--max_percent_outliers"), type="double", default=.05, help="When an outlier is spiked into a sample, the max percentage of outliers to spike (0.0 =< x =< 1.0)"),
make_option(c("-t","--percent_outlier_spikins"), type="double", default=.05, help="The percent of samples to spike in outliers (0.0 =< x =< 1.0)"),
make_option(c("-g","--int_multiplier_range"), type="double", default=3, help="Maximum value for the multiplier range."),
make_option(c("-m","--int_min_multiplier_range"), type="double", default=1, help="Minimum value for the multiplier range."),
make_option(c("-d","--int_multiplier_delta"), type="double", default=1, help="multiplier delta"),
make_option(c("-r","--collinear_range"), type="integer", default=5, help="collinear range for the spikes"),
make_option(c("-i","--collinear_increment"), type="integer", default=1, help="collinear delta for the spikes"),
make_option(c("-k","--percent_spiked"), type="double", default=.03, help="multiplier delta for the spikes"),
make_option(c("-c","--calibrate"), type="character", default=NA, help="Calibration file for generating the random log normal data. TSV file (column = feature)"),
make_option(c("-l","--minLevelPercent"), type="double", default=.1, help="Minimum number of instances a level can have in discontinuous metadata (Rounded up to the nearest count)."),
make_option(c("-u","--minOccurence"), type="integer", default=4, help="Minimum counts a bug can have for the ocurrence quality control filter used when creating bugs. ( Filtering minimum number of counts in a minimum number of samples)."),
make_option(c("-w","--minSample"), type="integer", default=4, help="Minimum samples a bug can be in for the ocurrence quality control filter used when creating bugs. ( Filtering minimum number of counts in a minimum number of samples)."),
make_option(c("-v","--verbose"), action="store_false", default = TRUE, help="If True logging and plotting is made by the underlying methodology."),
make_option(c("-z","--noZeroInflate"), action="store_true", default = FALSE, help="If given, zero inflation is not used and a lognormal distribution is generated for each feature")
)

main = function(
pArgs
){
lxArgs = parse_args(pArgs,positional_arguments = TRUE)

options <- lxArgs[['options']]
seed <- options[['seed']]
lefse_file <- options[['lefse_file']]
int_base_metadata_number <- options[['number_metadata']]
int_number_features <- options[['number_features']]
int_number_samples <- options[['number_samples']]
dPercentOutliers <- options[['max_percent_outliers']]
dPercentOutlierSpikins <- options[['percent_outlier_spikins']]
iReadDepth = options[['read_depth']]
int_multiplier_range <- options[['int_multiplier_range']]
int_min_multiplier_range <- options[['int_min_multiplier_range']]
int_multiplier_delta <- options[['int_multiplier_delta']]
collinear_range <- options[['collinear_range']]
collinear_increment <- options[['collinear_increment']]
dPercentMultSpiked <- options[['percent_spiked']]
strCalibrationFile = options[['calibrate']]
dMinLevelCountPercent = options[['minLevelPercent']]
dMinOccurenceCount = options[['minOccurence']]
dMinOccurenceSample = options[['minSample']]
fVerbose = options[['verbose']]
fZeroInflate = !options[['noZeroInflate']]

# locational arguments
file_names <- lxArgs[['args']]
strNormalizedFileName <- file_names[1]
strCountFileName = file_names[2]
parameter_filename <- file_names[3]
strBugBugInteractionFile <- file_names[4]
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
lsMetadataInfo <- func_generate_metadata(int_base_metadata_number,int_number_samples,dMinLevelCountPercent)
mat_metadata =  lsMetadataInfo[["mat_metadata"]]
metadata_parameters = lsMetadataInfo[["mtrxParameters"]]
vParametersAssociations = c(vParametersAssociations,lsMetadataInfo[["mtrxParameters"]])

# generate plain random lognormal bugs
pdf(file.path(dirname(strCountFileName),"FuncGenerateRLNorm.pdf"), useDingbats=FALSE)
# Get the fitted values for calibrating rlnorm
vdMu = NA
vdSD = NA
vdPercentZero = NA
dSDBeta = c_dSDBeta
dSDIntercept = c_dSDIntercept
dBetaZero = c_dBetaZero
dGrandMu = c_dGrandMu
dGrandSD = c_dGrandSD
dInterceptZero = c_dZeroIntercept
dGrandBeta = c_dBetaGrandSD
dGrandIntercept = c_dInterceptGrandSD
dMaxCount = iReadDepth

if(!is.na(strCalibrationFile) & (strCalibrationFile!="NA"))
{
  # Get the fit for the data
  print("Calibrating...")
  lsFit = funcCalibrateRLNormToMicrobiome(strCalibrationFile,fVerbose)

  vdMu = lsFit[["mu"]]
  vdSD = lsFit[["sd"]]
  vdPercentZero = lsFit[["percentZero"]]
  dInterceptZero = lsFit[["dZeroIntercept"]]
  dSDBeta  = lsFit[["dSDBeta"]]
  dSDIntercept  = lsFit[["dSDIntercept"]]
  dBetaZero = lsFit[["dZeroBeta"]]
  dGrandMu = lsFit[["dGrandMu"]]
  dGrandSD = lsFit[["dGrandSD"]]
  dGrandBeta = lsFit[["dGrandBeta"]]
  dGrandIntercept = lsFit[["dGrandIntercept"]]
  iReadDepth = lsFit[["dAverageReadDepth"]]
  dMaxCount = lsFit[["dMaxCount"]]
}
print("Parameters after Calibration File (if no calibration file is used, defaults are shown")
print(paste("Length vdMu", length(vdMu), "length vdSD", length(vdSD), "length vdPercentZero", length(vdPercentZero), "dSDBeta", dSDBeta,"dSDIntercept", dSDIntercept,"dBetaZero", dBetaZero, "dGrandMu", dGrandMu, "dGrandSD", dGrandSD, "dInterceptZero", dInterceptZero, "dGrandBeta", dGrandBeta, "dGrandIntercept", dGrandIntercept, "Read depth", iReadDepth))

mat_random_lognormal_bugs = func_generate_random_lognormal_matrix(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberCounts=dMinOccurenceCount, iMinNumberSamples=dMinOccurenceSample, iReadDepth=iReadDepth, vdMu=vdMu, xPercentZero=vdPercentZero, vdSD=vdSD, fZeroInflate=fZeroInflate, dBetaSD=dSDBeta, dSDIntercept=dSDIntercept, dBetaZero=dBetaZero, dInterceptZero=dInterceptZero, dGrandMu=dGrandMu, dBetaGrandSD=dGrandBeta, dInterceptGrandSD=dGrandIntercept, dMaxCount = dMaxCount, fVerbose=fVerbose)

lefse_lognormal = NULL
if(!is.null(lefse_file))
{
  lefse_file = dirname(lefse_file)
  lefse_lognormal <- func_generate_lefse_matrices(lefse_file, metadata_parameters, int_number_features, int_number_samples,mat_metadata, mat_random_lognormal_bugs[['mat_bugs']], 'lognormal')
}

vParametersAssociations = c(vParametersAssociations,mat_random_lognormal_bugs[["mtrxParameters"]])
list_of_bugs[[length(list_of_bugs) + 1]] = mat_random_lognormal_bugs[["mat_bugs"]]
hist(as.vector(mat_random_lognormal_bugs[["mat_bugs"]]), main="Final: Log normal matrix")
dev.off()

# generate random lognormal with outliers
pdf(file.path(dirname(strCountFileName),"LognormalWithOutliers.pdf"), useDingbats=FALSE)
mat_random_lognormal_outliers_bugs <- func_generate_random_lognormal_with_outliers(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberCounts=dMinOccurenceCount,iMinNumberSamples=dMinOccurenceSample, dMaxPercentOutliers=dPercentOutliers,dPercentSamples=dPercentOutlierSpikins,mtrxBugs=mat_random_lognormal_bugs[["mat_bugs"]],fVerbose=fVerbose)
lefse_outliers = NULL
if(!is.null(lefse_file))
{
  lefse_outliers <- func_generate_lefse_matrices(lefse_file, metadata_parameters,int_number_features,int_number_samples,mat_metadata, mat_random_lognormal_outliers_bugs[['mat_bugs']], 'outliers')
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

  for (mult in seq(int_min_multiplier_range, int_multiplier_range, int_multiplier_delta))
  {
    pdf(file.path(dirname(strCountFileName),paste("SpikeIn_n_", how_many_multivariates,"_m_", mult,".pdf",sep="")), useDingbats=FALSE)
    mat_random_lognormal_multivariate_spikes <- func_generate_random_lognormal_with_multivariate_spikes(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberCounts=dMinOccurenceCount, iMinNumberSamples=dMinOccurenceSample, percent_spikes=dPercentMultSpiked, multiplier=mult, metadata_matrix=mat_metadata, multivariate_parameter=how_many_multivariates, dMinLevelCountPercent=dMinLevelCountPercent, mtrxBugs=mat_random_lognormal_bugs[["mat_bugs"]],fZeroInflated=fZeroInflate, lviFrozeMetadataIndices=lviMetadata, liFrozeDataIndicies=liData, lsFrozeLevels=lsLevels, fVerbose=fVerbose)
    mat_random_lognormal_multivariate_spikes_bugs <- mat_random_lognormal_multivariate_spikes[["mat_bugs"]]
    lviMetadata <- mat_random_lognormal_multivariate_spikes[["MetadataIndices"]]
    liData <- mat_random_lognormal_multivariate_spikes[["DataIndices"]]
    lsLevels <- mat_random_lognormal_multivariate_spikes[["Levels"]]

    if(!is.null(lefse_file))
    {
      lefse_spike <- func_generate_lefse_matrices(lefse_file, metadata_parameters,int_number_features,int_number_samples,mat_metadata, mat_random_lognormal_multivariate_spikes_bugs, paste('multivariate_n_', how_many_multivariates, '_m_', mult,sep=""))
    }
    # generate known associations for random lognormal with spikes
    vParametersAssociations = c(vParametersAssociations,mat_random_lognormal_multivariate_spikes[["m_parameter_rec"]])
    list_of_bugs[[length(list_of_bugs)+1]] <- mat_random_lognormal_multivariate_spikes_bugs
    lsMicrobiomeKeys[[length(lsMicrobiomeKeys)+1]] = paste(c_strSpike,"n",how_many_multivariates,"m",sub(".","_",paste(mult),fixed=TRUE), sep="_")

    hist(as.vector(mat_random_lognormal_multivariate_spikes_bugs), main=paste("Final: Spiked matrix n_", how_many_multivariates,"_m_", mult,sep=""))
    dev.off()
  }
}

# Add bug associated bug microbiome
#lsBugBugInfo = func_generate_bug_bug_spiking_matrix(int_number_features=int_number_features, int_number_samples=int_number_samples, #iMinNumberSamples=dMinOccurenceSample, iReadDepth=iReadDepth, iMinNumberCounts=dMinOccurenceCount, vdMu=vdMu, vdSD=vdSD, dPercentZero=vdPercentZero, #fZeroInflate=fZeroInflate, dBetaSD=dSDBeta, dSDIntercept=dSDIntercept, dBetaZero=dBetaZero, fVerbose=fVerbose)
#list_of_bugs[[length(list_of_bugs)+1]] = lsBugBugInfo[["mtrxBugs"]]
#vParametersAssociations = c(vParametersAssociations,lsBugBugInfo[["vStrParameters"]])
#lsMicrobiomeKeys[[length(lsMicrobiomeKeys)+1]] = c_strBugBugAssocations

# preallocate final pcl  matrix
final_matrix <- matrix(data=NA,nrow=(number_metadata+int_number_features*length(list_of_bugs))+1, ncol=(int_number_samples+1))
final_matrix[1,1] = '#SampleID'
final_matrix[1,2:(int_number_samples + 1)] <- paste('Sample',1:int_number_samples,sep='')
final_matrix[2:(number_metadata+1),1] <- paste(c_strMetadata,1:number_metadata,sep='')
vdDim = dim(mat_metadata)
mat_metadata[(floor(vdDim[1]/2)+1):vdDim[1],] = paste("Group_",mat_metadata[(floor(vdDim[1]/2)+1):vdDim[1],],sep="")
final_matrix[2:(number_metadata+1),2:(int_number_samples+1)] <- mat_metadata

# Make a matrix for counts (use the other for normalized)
mtrxFinalCounts = final_matrix

start = 2 + number_metadata
end = (2+number_metadata) + (int_number_features-1)
iFirstData = start

for(iMatIndex in 1:length(list_of_bugs))
{
  final_matrix[start:end,1] <- paste(paste(c_strFeature,lsMicrobiomeKeys[iMatIndex],sep="_"),1:int_number_features,sep='_')
  # Normalize each column by thier sum and add to output
  final_matrix[start:end,2:(int_number_samples+1)] <- funcNormalizeMicrobiome(list_of_bugs[[iMatIndex]])
  mtrxFinalCounts[start:end,1] <- paste(paste(c_strFeature,lsMicrobiomeKeys[iMatIndex],sep="_"),1:int_number_features,sep='_')
  mtrxFinalCounts[start:end,2:(int_number_samples+1)] <- list_of_bugs[[iMatIndex]]
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
