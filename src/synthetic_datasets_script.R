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

source("synthetic_datasets_script_constants.R")
source("synthetic_datasets_script_helper_functions.R")
source("synthetic_datasets_script_bug_bug.R")
source("synthetic_datasets_script_options.R")

main = function(
pArgs
){
  lxArgs = parse_args(pArgs,positional_arguments = TRUE)
  print("lxArgs")
  print(lxArgs)
  options = lxArgs[['options']]

  # Get arguments and check defaults
  seed = options[[ 'seed' ]]

  lefse_file = options[[ 'lefse_file' ]]

  int_base_metadata_number = options[[ 'number_metadata' ]]
  if(int_base_metadata_number<1) stop("Please provide the base number for metadata generation as 1 or greater.")

  int_number_features = options[[ 'number_features' ]]
  if(int_number_features<1) stop("Please provide a number of features of atleast 1")

  strAssociationType = options[[ 'association_type' ]]
  if(strAssociationType != "linear") stop("Only linear associations supported at this point.")

  int_number_samples = options[[ 'number_samples' ]]
  if(int_number_samples<1) stop("Please provide a number of samples of atleast 1")

  dPercentOutliers = options[[ 'max_percent_outliers' ]]
  if( (dPercentOutliers>1) | (dPercentOutliers<0) ) stop("Please provide a percent outliers in the range of 0 to 1")

  dPercentOutlierSpikins = options[[ 'percent_outlier_spikins' ]]
  if( (dPercentOutlierSpikins>1) | (dPercentOutlierSpikins<0) ) stop("Please provide a percent spikins in the range of 0 to 1")

  iReadDepth = options[['read_depth']]
  if(iReadDepth < max(int_number_features, int_number_samples)) stop("Please provide a read depth of atleast equal to feature size or sample size (which ever is larger)")

  iNumAssociations = options[[ 'bugs_to_spike' ]]
  if(iNumAssociations<0) stop("Please provide a number of associations (bug-bug correlation) greater than or equal to 0")

  dVarScale = options[[ 'variance_scale' ]]
  if(dVarScale<0) stop("Please provide a variance scaling parameter greater than or equal to 0.")

  iMaxNumberCorrDomainBugs = options[[ 'max_domain_bugs' ]]
  if(iMaxNumberCorrDomainBugs<0) stop("Please provide a maximum number of domain features greater than or equal to 0")
  if(iMaxNumberCorrDomainBugs==0){
    iNumAssociations <-0
    warning("The maximum number of domain features is 0: setting the number of bug-bug associations to 0")
  }

  iNumberDatasets = options[[ "datasetCount" ]]
  if(iNumberDatasets<1) stop("Please provide a number of datasets which is at least 1.")

  dPercentMultSpiked = options[[ 'percent_spiked' ]]
  if( (dPercentMultSpiked>1) | (dPercentMultSpiked<0) ) stop("Please provide a percent multivariate spike in the range of 0 to 1")

  strCalibrationFile = options[[ 'calibrate' ]]

  dMinLevelCountPercent = options[[ 'minLevelPercent' ]]
  if( (dMinLevelCountPercent>1) | (dMinLevelCountPercent<0) ) stop("Please provide a min level percent in the range of 0 to 1")

  dMinOccurenceCount = options[[ 'minOccurence' ]]
  if(dMinOccurenceCount<0) stop("Please provide a min occurence greater than or equal to 0")

  dMinOccurenceSample = options[[ 'minSample' ]]
  if(dMinOccurenceSample<0) stop("Please provide a min sample greater than or equal to 0")
  if(dMinOccurenceSample>int_number_samples)
  {
    dMinOccurenceSample = int_number_samples
    print(paste("The min sample (for QC) was larger than the actual sample size, reset the min sample to the sample size, minSample is now equal to number_samples which is ",int_number_samples))
  }

  if ( options[[ 'scalePercentZeros' ]] < 0 ){ stop( "Please provide a scale percent zero greater than 0." ) }

  vdSpikeCount = as.integer(unlist( strsplit( options[[ 'spikeCount' ]], "," ) ) )

  vdSpikeStrength = as.double( unlist( strsplit( options[[ 'spikeStrength' ]], "," ) ) )

  if( length( vdSpikeCount ) != length( vdSpikeStrength ) )
  {
    if( length( vdSpikeCount ) == 1 )
    {
      vdSpikeCount = rep( vdSpikeCount, length( vdSpikeStrength ) )
    } else if( length( vdSpikeStrength ) == 1 )
    {
      vdSpikeStrength = rep( vdSpikeStrength, length( vdSpikeCount ) )
    } else {
      stop( "Please make sure to either provide the same length of spike-in counts and spike-n strengths or provide one that is 1 value. These values will be paired in order or one will be repeated to pair with the other values." )
    }
  }

  strBugBugCoef = options[[ 'bugBugCoef' ]]
  vecBugBugCoef = as.vector(
      sapply(
          strsplit(strBugBugCoef,',')[[1]],
          as.numeric
          )
      )
  if( length(vecBugBugCoef)<2 ){
      stop("Please provide both an integer and a slope value for the coefficients")
  }
  
  mtrxSpikeConfig = cbind( vdSpikeCount, vdSpikeStrength )

  fVerbose     = options[[ 'verbose' ]]
  fZeroInflate = !options[[ 'noZeroInflate' ]]
  fRunMetadata = !options[[ "noRunMetadata" ]]
  fRunBugBug   = options[[ "runBugBug" ]]

  # locational arguments
  file_names = lxArgs[[ 'args' ]]
  strNormalizedFileName = file_names[1]
  strCountFileName = file_names[2]
  parameter_filename = file_names[3]
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

  # Holds key words for the feature names of the microbiomes
  lsMicrobiomeKeys = c( )

  # Default number of metadata
  number_metadata = 0

  if(fRunMetadata){

    # generate the metadata
    lsMetadataInfo = func_generate_metadata( int_base_metadata_number, int_number_samples,dMinLevelCountPercent )
    mat_metadata =  lsMetadataInfo[[ "mat_metadata" ]]
    metadata_parameters = lsMetadataInfo[[ "mtrxParameters" ]]
    vParametersAssociations = c(vParametersAssociations,lsMetadataInfo[[ "mtrxParameters" ]])

    # generate plain random lognormal bugs
    pdf(file.path(dirname(strCountFileName),"FuncGenerateRLNorm.pdf"), useDingbats=FALSE)
    # Get the fitted values for calibrating rlnorm
    vdExp = NA
    vdMu = NA
    vdSD = NA
    vdPercentZero = NA
    #dSDBeta = c_dSDBeta
    #dBetaZero = c_dBetaZero
    #dGrandBeta = c_dBetaGrandSD

    print("Parameters BEFORE Calibration File")
    print(paste("Length exp",NA,"Length vdMu", NA, "length vdSD", NA, "length vdPercentZero", NA, "Read depth", iReadDepth))

    if(!is.na(strCalibrationFile) & (strCalibrationFile!="NA"))
    {
      # Get the fit for the data
      print("Calibrating...")
      lsFit = funcCalibrateRLNormToMicrobiome(strCalibrationFile, fVerbose)
      vdExp = lsFit[["exp"]]
      vdMu = lsFit[["mu"]]
      vdSD = lsFit[["sd"]]
      vdPercentZero = lsFit[["percentZero"]]
      iReadDepth = lsFit[["dAverageReadDepth"]]
      int_number_features = lsFit[["iFeatureCount"]]
    }
    print("Parameters AFTER Calibration File (if no calibration file is used, defaults are shown)")
    print( paste( "Length exp",           length(vdExp),
                  "Length vdMu",          length(vdMu), 
                  "length vdSD",          length(vdSD), 
                  "length vdPercentZero", length(vdPercentZero), 
                  "Read depth",           iReadDepth, 
                  "Feature Count",        int_number_features ) )

    mat_random_lognormal_bugs = func_generate_random_lognormal_matrix( int_number_features = int_number_features, 
                                                                       int_number_samples  = int_number_samples, 
                                                                       iMinNumberCounts    = dMinOccurenceCount, 
                                                                       iMinNumberSamples   = dMinOccurenceSample, 
                                                                       iReadDepth          = iReadDepth, 
                                                                       vdExp               = vdExp, 
                                                                       vdMu                = vdMu, 
                                                                       vdPercentZero       = vdPercentZero, 
                                                                       vdSD                = vdSD, 
                                                                       fZeroInflate        = fZeroInflate, 
                                                                       lSDRel              = list(BetaSD=c_dSDBeta, InterceptSD=c_dSDIntercept), 
                                                                       lPercentZeroRel     = list( InterceptZero = c_dInterceptZero, 
                                                                                                   BetaZero      = c_dBetaZero, 
                                                                                                   Beta2Zero     = c_dBeta2Zero, 
                                                                                                   Scale         = options[[ 'scalePercentZeros' ]]), 
                                                                       dBetaGrandSD        = c_dBetaGrandSD, 
                                                                       fVerbose            = fVerbose )

    lefse_lognormal = NULL
    if(!is.null(lefse_file))
    {
      lefse_file = dirname(lefse_file)
      lefse_lognormal = func_generate_lefse_matrices(lefse_file, metadata_parameters, int_number_features, int_number_samples, mat_metadata, mat_random_lognormal_bugs[['mat_bugs']], 'lognormal')
    }

    vParametersAssociations = c(vParametersAssociations,mat_random_lognormal_bugs[["mtrxParameters"]])
    list_of_bugs[[length(list_of_bugs) + 1]] = mat_random_lognormal_bugs[["mat_bugs"]]
    hist(as.vector(mat_random_lognormal_bugs[["mat_bugs"]]), main="Final: Log normal matrix")
    dev.off()

    # generate random lognormal with outliers
    pdf(file.path(dirname(strCountFileName),"LognormalWithOutliers.pdf"), useDingbats=FALSE)
    mat_random_lognormal_outliers_bugs = func_generate_random_lognormal_with_outliers( int_number_features = int_number_features, 
                                                                                       int_number_samples  = int_number_samples,
                                                                                       dMaxPercentOutliers = dPercentOutliers,
                                                                                       dPercentSamples     = dPercentOutlierSpikins,
                                                                                       mtrxBugs            = mat_random_lognormal_bugs[["mat_bugs"]],
                                                                                       fVerbose            = fVerbose )
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
    lsMicrobiomeKeys[[1]] = c_strRandom
    lsMicrobiomeKeys[[2]] = c_strOutlier

    # There are 4 groups of metadata (2 continuous, binary, and quarternery)
    number_metadata = c_iCountTypesOfMetadata * int_base_metadata_number

    # Hold the info to freeze the levels in the data
    llsLevels = list()
    lliMetadata = list()
    lliData = list()

    # Generate random lognormal with varying amounts of spikes
    for( iIndex in 1:nrow( mtrxSpikeConfig ) )
    {
      vMultConfig = mtrxSpikeConfig[ iIndex, ]
      iSpikeCount = vMultConfig[ 1 ]
      sKey = as.character( iSpikeCount )
      iSpikeStrength = vMultConfig[ 2 ]
      lsLevels = NULL
      liData = NULL
      lviMetadata = NULL

      if( sKey %in% names( llsLevels ) )
      { 
        lsLevels = llsLevels[[ sKey ]] 
        liData = lliData[[ sKey ]]
        lviMetadata = lliMetadata[[ sKey ]]
      }

      pdf(file.path(dirname(strCountFileName), paste("SpikeIn_n_", iSpikeCount,"_m_", iSpikeStrength,".pdf",sep="")), useDingbats=FALSE)
      mat_random_lognormal_multivariate_spikes = func_generate_random_lognormal_with_multivariate_spikes( int_number_features     = int_number_features, 
                                                                                                          int_number_samples      = int_number_samples,  
                                                                                                          percent_spikes          = dPercentMultSpiked, 
                                                                                                          multiplier              = iSpikeStrength, 
                                                                                                          metadata_matrix         = mat_metadata, 
                                                                                                          multivariate_parameter  = iSpikeCount, 
                                                                                                          dMinLevelCountPercent   = dMinLevelCountPercent, 
                                                                                                          mtrxBugs                = mat_random_lognormal_bugs[["mat_bugs"]], 
                                                                                                          fZeroInflated           = fZeroInflate, 
                                                                                                          lviFrozeMetadataIndices = lviMetadata, 
                                                                                                          liFrozeDataIndicies     = liData, 
                                                                                                          lsFrozeLevels           = lsLevels, 
                                                                                                          fVerbose                = fVerbose )
      mat_random_lognormal_multivariate_spikes_bugs = mat_random_lognormal_multivariate_spikes[["mat_bugs"]]

      lliMetadata[[ sKey ]] = mat_random_lognormal_multivariate_spikes[["MetadataIndices"]]
      lliData[[ sKey ]] = mat_random_lognormal_multivariate_spikes[["DataIndices"]]
      llsLevels[[ sKey ]] = mat_random_lognormal_multivariate_spikes[["Levels"]]

      if(!is.null(lefse_file))
      {
        lefse_spike = func_generate_lefse_matrices( lefse_file, 
                                                    metadata_parameters, 
                                                    int_number_features, 
                                                    int_number_samples,
                                                    mat_metadata, 
                                                    mat_random_lognormal_multivariate_spikes_bugs, 
                                                    paste('multivariate_n_', iSpikeCount, '_m_', mult,sep=""))
      }
      # generate known associations for random lognormal with spikes
      vParametersAssociations = c( vParametersAssociations, mat_random_lognormal_multivariate_spikes[["m_parameter_rec"]])
      list_of_bugs[[length(list_of_bugs)+1]] = mat_random_lognormal_multivariate_spikes_bugs
      lsMicrobiomeKeys[[length(lsMicrobiomeKeys)+1]] = paste(c_strSpike,"n", iSpikeCount,"m",sub(".","_",paste(iSpikeStrength),fixed=TRUE), sep="_")

      hist(as.vector(mat_random_lognormal_multivariate_spikes_bugs), main=paste("Final: Spiked matrix n_", iSpikeCount,"_m_", iSpikeStrength,sep=""))
      dev.off()
    }
  }

  if(fRunBugBug){

    # Get the association type and the coefficients
    if( strAssociationType=="linear" ){
        funcAssociation = func_linear_association
    }
    lAssociationParams = list(
        dIntercept = vecBugBugCoef[1],
        vdSlope     = vecBugBugCoef[-1]
        )

      
    # Get the fitted values for calibrating rlnorm
    vdExp = NA
    vdMu = NA
    vdSD = NA
    vdPercentZero = NA

    print("Parameters for Bug-Bug spikes BEFORE Calibration File")
    print(paste("Length exp",NA,"Length vdMu", NA, "length vdSD", NA, "length vdPercentZero", NA, "Read depth", iReadDepth))

    if(!is.na(strCalibrationFile) & (strCalibrationFile!="NA"))
    {
      # Get the fit for the data
      print("Calibrating...")
      lsFit = funcCalibrateRLNormToMicrobiome(strCalibrationFile, fVerbose)
      vdExp = lsFit[["exp"]]
      vdMu = lsFit[["mu"]]
      vdSD = lsFit[["sd"]]
      vdPercentZero = lsFit[["percentZero"]]
      iReadDepth = lsFit[["dAverageReadDepth"]]
      int_number_features = lsFit[["iFeatureCount"]]
    }
    print("Parameters for Bug-Bug spikes AFTER Calibration File (if no calibration file is used, defaults are shown)")
    print( paste( "Length exp",           length(vdExp),
                  "Length vdMu",          length(vdMu), 
                  "length vdSD",          length(vdSD), 
                  "length vdPercentZero", length(vdPercentZero), 
                  "Read depth",           iReadDepth, 
                  "Feature Count",        int_number_features ) )


    # Get the initial mu vector for generating features so that the datasets are iid   
    lsInitialDistribution = funcGenerateFeatureParameters(int_number_features = int_number_features,
                                                          int_number_samples  = int_number_samples,
                                                          iMinNumberSamples   = dMinOccurenceSample,
                                                          iReadDepth          = iReadDepth,
                                                          vdExp               = vdExp,
                                                          vdMu                = vdMu,
                                                          vdSD                = vdSD,
                                                          vdPercentZero       = vdPercentZero,
                                                          lSDRel              = list( BetaSD = c_dSDBeta, InterceptSD = c_dSDIntercept ),
                                                          lPercentZeroRel     = list( InterceptZero = c_dInterceptZero,
                                                                                      BetaZero      = c_dBetaZero,
                                                                                      Beta2Zero     = c_dBeta2Zero,
                                                                                      Scale         = options[[ 'scalePercentZeros' ]]), 
                                                          dBetaGrandSD        = c_dBetaGrandSD,
                                                          fVerbose            = fVerbose)


    # Update the Mu, SD and Percent zero bugs and report on distributions                                                                                                                                          
    vdMu = lsInitialDistribution[["mu"]]
    vdSD = lsInitialDistribution[["sd"]]
    vdPercentZero = lsInitialDistribution[["PercentZero"]]
    vdExp = lsInitialDistribution[["exp"]]

    mtrxDistributionParameters = matrix(NA, nrow=5, ncol=1)
    mtrxDistributionParameters[1,1] = paste(c_strSyntheticMicrobiome, c_strDistributionParameters, sep="")
    mtrxDistributionParameters[2,1] = paste(c_strMuVector, toString( round( vdMu, 2 ) ))
    mtrxDistributionParameters[3,1] = paste(c_strSDVector, toString( round( vdSD, 2 ) ))
    mtrxDistributionParameters[4,1] = paste(c_strPercentZeroVector, toString( round( vdPercentZero, 2 ) ))
    mtrxDistributionParameters[5,1] = paste(c_strExpVector, toString( round( vdExp, 2 ) ))
    vParametersAssociations = c( vParametersAssociations,  mtrxDistributionParameters)


    # Get the indices for the associations
    lsAssociationIdx = func_get_corr_indices( iNumAssociations = iNumAssociations,
                                              iMaxDomainNumber = iMaxNumberCorrDomainBugs,
                                              iNumFeatures     = int_number_features )

    # Update the parameters
    viIdxCorrRangeBugs  = lsAssociationIdx[["RangeBugs"]]
    liIdxCorrDomainBugs = lsAssociationIdx[["DomainBugs"]]

    # Flag to add the parameters only once, since the datasets are iid                                              
    fAddedParameters    = FALSE
    fAddBasisParameters = FALSE

    for( iDataset in 1:iNumberDatasets ){
      # generate plain random lognormal bugs
      pdf(file.path(dirname(strCountFileName),paste("FuncGenerateRLNormBugBug_d_",iDataset,".pdf",sep="")), useDingbats=FALSE)


      mat_random_lognormal_bugs = func_generate_random_lognormal_matrix( int_number_features = int_number_features, 
                                                                         int_number_samples  = int_number_samples, 
                                                                         iMinNumberCounts    = dMinOccurenceCount, 
                                                                         iMinNumberSamples   = dMinOccurenceSample, 
                                                                         iReadDepth          = iReadDepth, 
                                                                         vdExp               = vdExp, 
                                                                         vdMu                = vdMu, 
                                                                         vdPercentZero       = vdPercentZero, 
                                                                         vdSD                = vdSD, 
                                                                         fZeroInflate        = fZeroInflate, 
                                                                         lSDRel              = list(BetaSD=c_dSDBeta, InterceptSD=c_dSDIntercept), 
                                                                         lPercentZeroRel     = list( InterceptZero = c_dInterceptZero, 
                                                                                                     BetaZero      = c_dBetaZero, 
                                                                                                     Beta2Zero     = c_dBeta2Zero, 
                                                                                                     Scale         = options[[ 'scalePercentZeros' ]]), 
                                                                         dBetaGrandSD        = c_dBetaGrandSD, 
                                                                         fVerbose            = fVerbose )

      if( !fAddedParameters ){
        vParametersAssociations = c(vParametersAssociations,
                                    mat_random_lognormal_bugs[["mtrxBasisParameters"]],
                                    paste(c_strNumberDatasets,iNumberDatasets))
        fAddedParameters==TRUE
      }
      list_of_bugs[[length(list_of_bugs) + 1]] = mat_random_lognormal_bugs[["mat_basis"]]
      lsMicrobiomeKeys[[length(lsMicrobiomeKeys)+1]] = paste(c_strBugBugAssociations,c_strNull,"d",iDataset,sep="_")
      hist(as.vector(mat_random_lognormal_bugs[["mat_basis"]]), main="Final: Basis log normal matrix")
      dev.off()

      if(iNumAssociations > 0)  # Only add associations if associations are requested
      {
        pdf(file.path(dirname(strCountFileName), paste("BugBugSpikeIn_d_", iDataset,sep="")), useDingbats=FALSE)
        mat_random_lognormal_bugbug_spikes = func_generate_bug_bug_spikes( mtrxData            = mat_random_lognormal_bugs[["mat_basis"]],
                                                                           dVarScale           = dVarScale,
                                                                           funcAssociation     = funcAssociation,
                                                                           lAssociationParams  = lAssociationParams,
                                                                           viIdxCorrRangeBugs  = viIdxCorrRangeBugs,
                                                                           liIdxCorrDomainBugs = liIdxCorrDomainBugs,
                                                                           iMaxDomainNumber    = iMaxNumberCorrDomainBugs,
                                                                           vdPercentZero       = vdPercentZero,
                                                                           fZeroInflate        = TRUE,
                                                                           fVerbose            = fVerbose )
        mat_random_lognormal_bugbug_spikes_bugs = mat_random_lognormal_bugbug_spikes[["mat_bugs"]]

        if( !fAddedParameters ){
          vParametersAssociations = c( vParametersAssociations,
                                       mat_random_lognormal_bugbug_spikes[["mtrxAssnParameters"]],
                                       paste(c_strNumberDatasets,iNumberDatasets) )
          fAddedParameters = TRUE
        }
        list_of_bugs[[length(list_of_bugs)+1]] = mat_random_lognormal_bugbug_spikes_bugs
        lsMicrobiomeKeys[[length(lsMicrobiomeKeys)+1]] = paste(c_strBugBugAssociations,"a",iNumAssociations,"d",iDataset,sep="_")

        hist(as.vector(mat_random_lognormal_bugbug_spikes_bugs), main=paste("Final: Bug-Bug Spiked matrix d=",iDataset,sep=""))
        dev.off()

      } else if(!fAddBasisParameters){   # Run just to get truth file parameters for the null matrices

        pdf(file.path(dirname(strCountFileName), "BugBugSpikeIn"), useDingbats=FALSE)
        mat_random_lognormal_bugbug_spikes = func_generate_bug_bug_spikes( mtrxData            = mat_random_lognormal_bugs[["mat_basis"]],
                                                                           dVarScale           = dVarScale,
                                                                           funcAssociation     = funcAssociation,
                                                                           lAssociationParams  = lAssociationParams,
                                                                           viIdxCorrRangeBugs  = NULL,
                                                                           liIdxCorrDomainBugs = NULL,
                                                                           iMaxDomainNumber    = iMaxNumberCorrDomainBugs,
                                                                           vdPercentZero       = vdPercentZero,
                                                                           fZeroInflate        = TRUE,
                                                                           fVerbose            = fVerbose )
        mat_random_lognormal_bugbug_spikes_bugs = mat_random_lognormal_bugbug_spikes[["mat_bugs"]]

        if( !fAddedParameters ){
          vParametersAssociations = c( vParametersAssociations,
                                       mat_random_lognormal_bugbug_spikes[["mtrxAssnParameters"]],
                                       paste(c_strNumberDatasets,iNumberDatasets) )
          fAddedParameters = TRUE
        }
        list_of_bugs[[length(list_of_bugs)+1]] = mat_random_lognormal_bugbug_spikes_bugs
        lsMicrobiomeKeys[[length(lsMicrobiomeKeys)+1]] = paste(c_strBugBugAssociations,"a",iNumAssociations,"d",iDataset,sep="_")

        hist(as.vector(mat_random_lognormal_bugbug_spikes_bugs), main=paste("Final: Bug-Bug Spiked matrix d=",iDataset,sep=""))
        dev.off()

        fAddBasisParameters = TRUE
      }
    }
  }
  # preallocate final pcl  matrix
  final_matrix = matrix(data=NA,nrow=(number_metadata+int_number_features*length(list_of_bugs))+1, ncol=(int_number_samples+1))
  final_matrix[1,1] = '#SampleID'
  final_matrix[1,2:(int_number_samples + 1)] = paste('Sample',1:int_number_samples,sep='')
  if(number_metadata > 0){
    final_matrix[2:(number_metadata+1),1] = paste(c_strMetadata,1:number_metadata,sep='')
    vdDim = dim(mat_metadata)
    mat_metadata[(floor(vdDim[1]/2)+1):vdDim[1],] = paste("Group_",mat_metadata[(floor(vdDim[1]/2)+1):vdDim[1],],sep="")
    final_matrix[2:(number_metadata+1),2:(int_number_samples+1)] = mat_metadata
  }

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
