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
c_dRunifMin = .01
c_dRunifMax = .99
### Continuous metadata are draw from uniform distributions, these are the bounds
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

# Define the relationship between different feature properties
# These variables are associated with settings for
# Calculating the SD and percent Zero based on the mu or expectation.
# The constants for these variables have been estimated from
# a real (IBD) data set. These constants are used unless a
# calibration file is given which estimates the values in the
# same manner as the constants were estimated.
c_dSDBeta = 0.1251#0.242678502955986
c_dSDIntercept = 1.121
### The estimate for the relationship between SD and exp
c_dBetaZero = -0.05603027#-0.09194 # -0.119637999937095 # .86160
c_dBeta2Zero = -0.0111924
c_dInterceptZero = 0.9073969
### The estimate for the relationship between exp and zero percent
c_dBetaGrandSD = 0.04982219 #0.01699164
### The estimate for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)

####
c_vdExpLog = c( -0.877070018720874,1.24126858906963,-0.71743987312899,-0.579818495252942,-0.916290731874155,-1.42711635564015,-0.350976922824095,-2.68824757380603,4.45611317802658,-1.39432653281715,-0.601479992034121,-1.03282454813011,4.6737256284836,2.52348613175548,-0.867500567704723,-1.6094379124341,-2.57702193869581,0.799307376388336,1.99768903980758,-1.58963528513792,7.55066334615255,3.13948623719869,1.50318811259139,0.0353671438372913,2.43220886078755,4.85881966369566,1.67372640231646,1.92629036218566,3.4618533755355,1.93036131866568,4.68792920115586,-0.0921152889078056,1.84498423046535,-2.02495335639577,-1.83258146374831,2.14710019015365,-0.785262469467751,1.0123279200711,0.542324290825362,3.30630044979247,-1.65025990695436,0.525911261184032,2.43571640559723,-0.653926467406664,4.43585172136879,0.790273891290668,0.72270598280149,0.0506931143155182,-0.71743987312899,1.21669157673371,5.10315370174078,1.5295285292058,-1.20397280432594,3.88362353090645,-2.52572864430826,1.82841278687609,-1.20397280432594,4.13523055469444,-0.524248644098131,-0.551647618286246,-0.0790432073404529,-1.24479479884619,6.32837954283977,-1.10262031006565,-0.20334092401803,-0.415515443961666,0.5607579925142,-0.233193887167711,-2.52572864430826,5.19168938132325,5.3469170363841,5.12508240135524,-1.80788885115794,-1.65025990695436,5.59600422381589,6.02419307773006,0.845009529869191,-2.52572864430826,-0.317454230785451,2.86539577182599,2.59614982860398,-0.776528789498996,0.501986675098786,-0.742337424750717,0.249980205267769,3.63989927189222,-1.09064411901893,6.50169199200066,2.4308023907948,2.11529119457353,4.66317490827467,1.17803942229943,4.3994244118503,5.69332801675234,3.43501893013767,-3.03655426807425,-0.988861424708991,4.46558622777281,2.94528073005789,3.11688739511046,4.15211019907224,3.28929700916393,2.14006616349627,5.78226088252746,0.681074599325676,3.17588481306311,5.01359092192189,1.95103982687531,0.877134016672961,2.47754629538171,4.2285840374303,4.32772883126603,-1.57021719928082,-2.88240358824699,-1.96611285637283,-1.78379129957888,0.29266961396282,2.67662833109387,1.04802050255205,1.04942204447734,-1.58963528513792,0.0544881852840698,-0.0877389143080067,1.38729386145297,-2.52572864430826,-0.2484613592985,-1.80788885115794,-0.867500567704723,-3.12356564506388,-1.76026080216868,-0.988861424708991,-0.877070018720874,0.289680075114454,-0.53785429615391,-2.15416508787577,-1.99510039324608,0.666803205220343,1.70765295993106,-2.52572864430826,5.39355481643488,1.92861865194525,1.27759494419655,0.227932068046007,-2.12026353620009,4.45066611039057,0.0769610411361284,1.26186428274171,1.57608793275255,3.5184469417124,0.710987098688276,-2.95651156040071,0.4265740713184,1.86128553187667,-0.661648513500574,-2.12026353620009,2.77832225408754,5.86404018285232,5.91260530552877 )

c_vdExp = c(0.416,3.46,0.488,0.56,0.4,0.24,0.704,0.068,86.152,0.248,0.548,0.356,107.096,12.472,0.42,0.2,0.076,2.224,7.372,0.204,1902.004,23.092,4.496,1.036,11.384,128.872,5.332,6.864,31.876,6.892,108.628,0.912,6.328,0.132,0.16,8.56,0.456,2.752,1.72,27.284,0.192,1.692,11.424,0.52,84.424,2.204,2.06,1.052,0.488,3.376,164.54,4.616,0.3,48.6,0.08,6.224,0.3,62.504,0.592,0.576,0.924,0.288,560.248,0.332,0.816,0.66,1.752,0.792,0.08,179.772,209.96,168.188,0.164,0.192,269.348,413.308,2.328,0.08,0.728,17.556,13.412,0.46,1.652,0.476,1.284,38.088,0.336,666.268,11.368,8.292,105.972,3.248,81.404,296.88,31.032,0.048,0.372,86.972,19.016,22.576,63.568,26.824,8.5,324.492,1.976,23.948,150.444,7.036,2.404,11.912,68.62,75.772,0.208,0.056,0.14,0.168,1.34,14.536,2.852,2.856,0.204,1.056,0.916,4.004,0.08,0.78,0.164,0.42,0.044,0.172,0.372,0.416,1.336,0.584,0.116,0.136,1.948,5.516,0.08,219.984,6.88,3.588,1.256,0.12,85.684,1.08,3.532,4.836,33.732,2.036,0.052,1.532,6.432,0.516,0.12,16.092,352.144,369.668)

c_vdPercentZero = c( 0.964,0.672,0.988,0.92,0.988,0.952,0.88,0.984,0.428,0.936,0.864,0.916,0.472,0.508,0.972,0.956,0.98,0.92,0.672,0.952,0.064,0.72,0.832,0.984,0.66,0.332,0.864,0.908,0.728,0.884,0.392,0.888,0.748,0.96,0.976,0.84,0.952,0.972,0.788,0.948,0.956,0.784,0.848,0.944,0.7,0.972,0.9,0.884,0.924,0.776,0.096,0.696,0.936,0.58,0.984,0.772,0.936,0.352,0.96,0.88,0.932,0.936,0.072,0.952,0.92,0.912,0.876,0.92,0.968,0.136,0.052,0.12,0.988,0.944,0.2,0.02,0.88,0.98,0.916,0.464,0.452,0.892,0.8,0.932,0.78,0.22,0.936,0.212,0.504,0.608,0.484,0.756,0.256,0.064,0.752,0.988,0.944,0.54,0.948,0.824,0.556,0.492,0.548,0.016,0.788,0.932,0.092,0.48,0.784,0.396,0.152,0.732,0.972,0.984,0.98,0.98,0.932,0.7,0.884,0.972,0.924,0.916,0.936,0.94,0.984,0.944,0.964,0.936,0.984,0.988,0.964,0.98,0.912,0.952,0.972,0.972,0.928,0.728,0.976,0.416,0.932,0.952,0.952,0.984,0.416,0.876,0.968,0.772,0.608,0.82,0.988,0.948,0.892,0.932,0.98,0.84,0.216,0.612 )

c_vdSD = c(1.29288222850903,1.11656790017514,1.92641197859569,1.02105677634939,1.48524719110303,0.870378434477115,0.959635078912515,0.583945419522526,1.79348082582133,0.669655926119326,0.600623613084851,0.853648651939737,1.72677780263148,1.26147280133706,0.734982664869215,0.914826396438766,0.719267497333394,1.36445803823611,1.24480603198208,0.554831763342855,1.70747836197823,1.49834967389745,1.10454940596015,1.86473285616939,1.1385848782746,1.4834337653741,1.52109111483939,1.46795133122398,1.98185399587786,1.45950174080202,1.59366893620959,1.03774612854838,1.44800233026963,0.592107847055163,0.594358593206619,1.6320293719476,0.962434020985145,2.07513585839894,1.09851647985654,2.49835214683079,0.976161247538305,1.10000715772956,1.80999152731209,1.15455694103497,2.29447750308131,1.84856812708444,1.18034213829724,0.938036155588448,0.874004621310478,1.09917849384688,1.7434683468884,1.16454079297455,0.728859758118786,1.77099312319428,0.694335857646383,1.37928702782278,0.76009658704113,1.68201745534149,1.35054907979964,0.832670104071402,1.16865528891661,0.847430958060489,1.5524935357744,0.904423733977812,1.17860611568599,1.03253595130507,1.29763235709243,1.06536965945072,0.27424593232245,1.57376990950369,1.34803070031666,1.55956753125577,1.10470442660092,0.723280490439043,1.5977407349549,0.930850462498865,1.29438882062999,0.803822732571653,0.986243910474401,1.29553952694038,1.28823184910116,0.685283831641858,0.946823554040535,0.918659788437422,0.788068631230099,1.22250675713921,0.819552326795716,1.59822401385881,1.17461364613988,1.17170376599601,1.87008425727893,1.11396621015164,1.41133141289528,1.38435992039007,1.37755755926622,0.256091402042052,0.758162389852548,1.64572542671349,2.14130589750552,1.88179351330357,1.41089291208342,1.60224644297979,1.14667348092866,0.942793883548745,0.922341311891637,1.96208871494235,1.56299028564406,1.01027693685149,1.02567334310102,1.17471546316181,1.46981919823393,2.08186319151382,1.17306243431671,0.758804728095711,0.843696416029418,0.670981308304214,1.15201957724894,1.56952695767858,1.42621506263831,2.07361684190959,0.751160258829315,1.2533916605354,1.41373507230623,1.73889522835292,1.00489987867617,1.23624923115956,0.556618232166134,1.1257008349353,0.612456251969496,0.824857904639158,1.25061326841638,1.27907129252095,1.03099622126377,1.1743855423561,0.685788698475515,1.22474258547412,1.49553187222792,1.23068112783264,0.421237925499913,2.25088768513035,2.24179750266376,1.88939818126567,1.42078414973316,0.830682791171414,1.8296524259332,1.16414627194934,2.22431472689055,1.27214629360624,1.64269420055488,1.20601967125448,0.294925311390253,1.7163875280654,1.62333836869266,1.16870601697591,0.618214596868589,1.751366359945,2.40861367145808,1.4507511261923)

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
  print( "start funcCalibrateRLNormToMicrobiome" )
  # Read in file
  print( "Reading file." )
  dfData = read.table( sCalibrationFile )
  row.names( dfData ) = dfData[[1]]
  dfData = dfData[ -1, -1 ]

  # Get read depth of the samples (From a tsv file, samples = rows)
  ldReadDepths = sapply( 1:nrow( dfData ), function( x ) sum( as.numeric( as.matrix( dfData )[ x, ] ) ) )

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
  for( iIndex in 1:ncol( dfData ) ){
    # Get the percent zero before removing zeros for other measurements
    vdCur = as.numeric( as.vector( as.matrix( dfData[ iIndex ] ) ) )
    vdPercentZero = c( vdPercentZero, length( which( vdCur == 0 ) ) / length( vdCur ) )

    # Measure expectation of the feature with zeros
    vdExp = c( vdExp, mean( vdCur ) )

    # Remove zeros from data
    vdCur = vdCur[ which( vdCur != 0 ) ]
    
    #### Note
    #### rlnorm needs a mean and sd from a logged rlnorm distribution which would match this
    #### without further manipulation. The "mean" in the formula is actually not the expectation
    #### The expectation is e^mu+.5*sd^2 so this is always more than mu.

    # Log nonzero data
    vdLogCur = log( vdCur )
    vdLogSD = c( vdLogSD, sd( vdLogCur ) )
    vdMu = c( vdMu, funcGetMu( vdExp[iIndex], exp(sd( vdLogCur ) ) ) )
  }

  # Estimate the distribution parameters from the expectation vector
  # Includes the relationship between the grand mu and grand sd
  # The grand mu, grand expectation (logged) and the grand sd
  lParams = funcGenerateExpVectorParameters( vdExp, TRUE )

  ##### Get relationship between logSD and log(Exp)
  # Log to make a linear relationship
  # This is both values logged as SD is based on vdLogCur
  vdExpLog = log( vdExp )

  lmod = lm( vdLogSD ~ vdExpLog ) 
  dBetaSD = coef( lmod )[ "vdExpLog" ]
  dInterceptSD = coef(lmod)[ "(Intercept)" ]

  #### Percent Zero and Exp
  ### Estimated with a polynomial to the second degree
  dBeta2Zero = NA
  dBetaZero = NA
  dInterceptZero = NA
  if( sum( vdPercentZero ) > 0 ){
    # Fit with polynomial regression (degree 2)
    lmodPoly = lm( vdPercentZero ~ poly( vdExpLog, 2, raw=T ) )
    dBeta2Zero = coef( lmodPoly )[ "poly(vdExpLog, 2, raw = T)2" ]
    dBetaZero = coef( lmodPoly )[ "poly(vdExpLog, 2, raw = T)1" ]
    dInterceptZero = coef( lmodPoly )[ "(Intercept)" ]
  }

  if( fVerbose ){
    # Indicate the relationships found
    print("The following relationships were found.")
    print("***Grand distribution***")
    print(paste("Expectation of feature expectations vs SD of feature expectations:",lParams$GrandLogSDBeta))
    print("***Feature distributions***")
    print("Log Exp (with zeros) vs Log SD (without zeros):")
    print(paste("Intercept=", dInterceptSD, ", Beta=", dBetaSD))
    print("Log Exp (with zeros) vs Percent Zeros:")
    print(paste(paste("Intercept=", dInterceptZero, ", Beta(degree1)=", dBetaZero,", Beta (degree2)=", dBeta2Zero)))

    # Initial description of calibration data set
    barplot( ldReadDepths, main = paste( "Read depth of Template Microbiome, Average", mean( ldReadDepths ) ) )
    abline( mean( ldReadDepths ), 0, col = "red" )
    hist( as.numeric( as.vector( as.matrix( dfData ) ) ), main = "Histogram of Template Microbiome Counts" )
    hist( log( as.numeric( as.vector( as.matrix( dfData ) ) ) ), main = "Histogram of Logged Template Microbiome Counts" )
    plot( as.numeric( as.vector( as.matrix( dfData ) ) ), main = "Distribution of Template Microbiome" )

    # Get relationship between SD and Log Exp
    plot( vdExpLog, vdLogSD, main = "Estimating SD with Log Exp" )
    points( x = vdExpLog, y = funcEstimateFeatureSD( vdExpLog, dBetaSD, 0 ), col = "green" )
    plot( vdExpLog, vdLogSD, main = "Estimating SD with Log Exp with Intercept" )
    points( x = vdExpLog, y = funcEstimateFeatureSD( vdExpLog, dBetaSD, dInterceptSD ), col = "green" )

    # Percent zero
    if( !is.na( dBetaZero ) ){
      print(c(dBeta2Zero, dBetaZero,dInterceptZero))
      plot( vdExpLog, vdPercentZero, main = "Estimating Percent Zero" )
      points( x = vdExpLog, y = funcEstimatePercentZero( vdExpLog, dBetaZero, dBeta2Zero, dInterceptZero ), col = 'violet' ) 
      plot( vdExpLog, vdPercentZero, main = "Estimating Percent Zero no intercept" )
      points( x = vdExpLog, y = funcEstimatePercentZero( vdExpLog, dBetaZero, dBeta2Zero, 0 ), col = 'violet' ) 
    }
  }

  print( "stop funcCalibrateRLNormToMicrobiome" )
  return( list( exp = vdExp, mu = vdMu, sd = exp( vdLogSD ), percentZero = vdPercentZero,
              dAverageReadDepth = mean( ldReadDepths ), iFeatureCount = ncol( dfData ) ) )

### When returning the grand Mu remember that you are returning the Mu that gives the expectation for the Mus
### given the rlnorm function so this is different than the mus measured in the logged distribution (rlnorm)
### exp: Not Logged expectation of distribution ( mean(x) )
### mu: The exponentiated logMu (calculated not measured)
### sd: The exponentiated logSD (based on the exp calculated without zeros)
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
funcEstimateFeatureSD = function(
### Estimate the SD given the Log Exp and parameters modeling the relationship between the Log Exp and the SD
vdExpLog,
### The measured mean of feature values
dBetaSD,
### The beta for the relationship between the mu and the SD
dInterceptSD = 0
){
  return( dInterceptSD  + ( vdExpLog * dBetaSD ) )
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


funcEstimatePercentZero = function(
### Estimate the percent zero given the logged expectation and parameters modeling the relationship between the log Exp and the percent zero
vdExpLog,
### The measured mean of feature values
dBetaZero,
### The beta for the relationship between the logged exp and the percent zero (for degree 1)
dBeta2Zero,
### The beta for the relationship between the logged exp and the percent zero (for degree 2)
dInterceptZero = 0
### The intercept for the polynomial relationship
){
  return(dInterceptZero  + (vdExpLog * dBetaZero) + (vdExpLog^2 * dBeta2Zero))
  ### Returns the Percent Zero (between 0 and 1) for a feature
}


funcForceMinCountsInMinSamples = function(
### For a feature to pass the requirement of having a certain minimal count in a minimal number of samples
vdFeature,
### Vector of signal (integers). Same length as vdLeftOver.
vdLeftOver,
### Vector of left over signal to add in. Same length as vdLeftOver.
iMinNumberCounts = 0,
### Minimum number of counts for a sample to pass the filter
iMinNumberSamples = 0
### Min number of samples to have the minimum number of counts
){
  # Check to make sure there are enough none zeros to add signal to for the min number of samples
  iNonZeroSamplesNeeded = min( 0, iMinNumberSamples - length( which( vdFeature == 0 ) ) )
  if(iNonZeroSamplesNeeded > 0)
  {
    # If more samples are needed, add them back in as the mean and remove signal from the left over
    dSignalMean = round( mean( vdFeature[ which( vdFeature > 0 ) ] ) )
    viUpdate = sample( which( vdFeature == 0 ), iNonZeroSamplesNeeded )
    vdFeature[ viUpdate ] = dSignalMean

    # Try removing the signal from the appropriate LeftOver element
    # If there is not enough left over to remove, keep and then remove randomly
    dRandomRemove = 0
    for( iUpdatedSample in viUpdate )
    {
      vdLeftOver[ iUpdatedSample ] = vdLeftOver[ iUpdatedSample ] - dSignalMean
      if( vdLeftOver[ iUpdatedSample ] < 0 )
      {
        dRandomRemove = dRandomRemove + abs(vdLeftOver[ iUpdatedSample ])
        vdLeftOver[ iUpdatedSample ] = 0
      }
    }
    dRandomRemove = floor( dRandomRemove )
    # Remove Randomly if needed
    for( iCount in 1:dRandomRemove )
    {
      viNonZeroLeftOver = which( vdLeftOver > 0 )
      if( length( viNonZeroLeftOver ) )
      {
        iRemoveCount = sample( viNonZeroLeftOver, 1)
        vdLeftOver[ iRemoveCount ] = vdLeftOver[ iRemoveCount ] - 1
      }
    }
  }

  # Min is placed here because if you ask for 2 samples but this is evaluating to find atleast 3 samples at a certain level you will inf loop
  # Look to see how many measurements are above the min.
  # If there are not enough counts more the the min number of counts
  # then inflate the distribution up until there is enough counts.
  iNeededExtraValues = iMinNumberSamples - length( which( vdFeature >= iMinNumberCounts ) )

  # Probability for each sample to get a count is based on it's current percentage of samples (multinomial distribution)
  vdProbabilities = vdFeature / sum( vdFeature )

  # Keep adding left over to the distribution until you have samples above the min count.
  # Keep going even if you run out of LeftOver
  # Add and remove counts to indices that are not zero
  viAddIndices = which( vdFeature > 0 )
  viRemoveFromLeftOverIndices = which( vdLeftOver > 0 )

  # While we need to add counts
  while( iNeededExtraValues > 0 )
  {
    # Index to add to, if multiple possibilities, select using the current percentage counts.
    iIndexAdd = viAddIndices
    if( length( iIndexAdd ) > 1 )
    {
      iIndexAdd = sample( viAddIndices, 1, prob = vdProbabilities[ viAddIndices ] )
    }

    vdFeature[ iIndexAdd ] = vdFeature[ iIndexAdd ] + 1
    iNeededExtraValues = iMinNumberSamples - length( which( vdFeature >= iMinNumberCounts ) )

    # Remove a count from Left Over, select using the current percentage counts.
    if( length(viRemoveFromLeftOverIndices) > 1 )
    {
      viRemoveFromLeftOverIndices = sample( viRemoveFromLeftOverIndices, 1, prob = vdLeftOver[ viRemoveFromLeftOverIndices ] / sum( vdLeftOver[ viRemoveFromLeftOverIndices ] ) )
    }
    if( length( viRemoveFromLeftOverIndices == 1 ) )
    {
      vdLeftOver[ viRemoveFromLeftOverIndices ] = vdLeftOver[ viRemoveFromLeftOverIndices ] - 1
      vdLeftOver[ vdLeftOver < 0 ] = 0
      viRemoveFromLeftOverIndices = which( vdLeftOver > 0 )
    }
  }
  return( list( Feature = vdFeature, LeftOver = vdLeftOver) )
}


# 2 Tests 9/4/2013
funcGenerateExpVectorParameters = function(
### Get the point estimate for the relationship between the mu and sd
vdExpectations,
### These are the mus of each of the features from the calibration file.
### These are untransformed (not logged) and the expectation with zeros
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


#funcGenerateBugCorrelation = function(
#
#){
#
#
#
#}


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
### Vector of Mu parameters for the grand exp distribution (means of features) if not supplied, one will be generated by rlnorm
vdSD = NA,
### Vector of SD parameters for the grand exp distribution (means of features) if not supplied, one will be generated by rlnorm
vdPercentZero = NA,
### Vector of percent zero parameters for features if not supplied, one will be generated by rlnorm
fZeroInflate = TRUE,
### Controls if zero inflation is used.
lSDRel = NA,
### Parameters for the relatioship between the exp and sd in which generates the expectation vector
lPercentZeroRel = NA,
### Parameters for the relatioship between the percent zero and feature expectations
dBetaGrandSD = NA,
### The beta for the relationship between the grand Exp and the SD
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
                                       lSDRel = lSDRel,
                                       lPercentZeroRel = lPercentZeroRel,
                                       dBetaGrandSD = dBetaGrandSD,
                                       fVerbose = fVerbose)

  print("RRRR")
  print(sum(mu_vector[["exp"]]))

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
                                                    vdPercentZero = mu_vector[["PercentZero"]],
                                                    vdSD = mu_vector[["sd"]],
                                                    fZeroInflate = fZeroInflate,
                                                    lSDRel = lSDRel,
                                                    lPercentZeroRel = lPercentZeroRel,
                                                    dBetaGrandSD = dBetaGrandSD,
#                                                    funcUpdateData = funcGenerateBugCorrelation
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
lSDRel = NA,
### If vdSD is not given, SD will be generated by a relationship with Exp parameters using the parameters given here
lPercentZeroRel = NA,
### If vdPercentZero is not given, vdPercentZero will be generated by a relationship with Exp using the parameters given here
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
    print("funcGenerateFeatureParameters: Generating vdExp Vector.")

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

    # Remember this is a vector of expectation not of log mu parameters.
    lExpVectorReturn = funcTruncatedRLNorm(iNumberMeasurements=int_number_features, dLogMean=lsParams$dLogMu, dLogSD=lsParams$dLogSD, iThreshold=exp((c_iTimesSDIsOutlier*lsParams$dLogSD)+lsParams$dLogMu))
    vdExp = lExpVectorReturn$Feature

    ### Update the distribution to the sum (read depth) requested.
    ### Depending on how many features are requested this is more or less needed
    ### This is not needed at the limit with many features.
    vdExp = funcUpdateDistributionToSum( vdExp, iReadDepth )

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

  if(is.na(vdSD))
  {
    # Generate vector of SD based on mu since it is not known
    print("funcGenerateFeatureParameters: Generating vdSD Vector.")
    vdSD = exp(funcEstimateFeatureSD(log(vdExp), lSDRel$BetaSD, lSDRel$InterceptSD))

    # If there is only 1 feature then the feature can not vary because its expectation would
    # Also be the read depth that needs to be preserved. Given there is only 1 feature no other
    # Features will be available to compensate for fluctations.
    if(length(vdExp)==1){vdSD = 0}

    # Floor to close to 0 because less can not be logged
    viInvalid = which( vdSD <= 0 )
    if( length( viInvalid ) > 0 )
    {
      print( paste( "funcGenerateFeatureParameters: Changing low SDs to a little more than 0. # occurences = ", length( viInvalid ) ) )
      vdSD[ viInvalid ] = c_dALittleMoreThanZero
    }

    if(fVerbose)
    {
      vdVisExp = vdExp
      vdVisExp[which(vdVisExp==0)] = c_dALittleMoreThanZero
      vdVisSD = vdSD
      vdVisSD[which(vdVisSD==0)] = c_dALittleMoreThanZero
      plot(vdVisExp, vdVisSD, main="Generated Relationship of Exp and SD", col="orange")
    }
  }

  if(is.na(vdMu))
  {
    print("funcGenerateFeatureParameters: Generating vdMu Vector.")
    # We know the vdExp for each sample
    # We know the SD for each sample
    vdMu = sapply(1:length(vdExp), function(x) funcGetMu(vdExp[x],vdSD[x]))
  }

  if(is.na(vdPercentZero))
  {
    # Generate vector of percent zero based on exp since it is not known
    print("funcGenerateFeatureParameters: Generating vdPercentZero Vector.")
    vdPercentZero = funcEstimatePercentZero( log(vdExp), lPercentZeroRel$BetaZero, lPercentZeroRel$Beta2Zero, lPercentZeroRel$InterceptZero ) 

    viLessThanZero = which(vdPercentZero < 0)
    if(length(viLessThanZero>0))
    {
      print(paste("funcGenerateFeatureParameters: Changing negative Percent Zeros to 0. # occurences = ",length(viLessThanZero)))
      vdPercentZero[vdPercentZero < 0] = 0
    }
    if(fVerbose){plot(log(vdExp), vdPercentZero, main="Generated Relationship of PercentZero and Log exp", col="purple")}
  }

#!# Removed resampling of exp, if a calibrated file is used, the number of features are frozen.
#  if(!length(vdMu)==int_number_features)
#  {
#    # This is the scenario that the calibration file is used and the number of the samples needed are not equal
#    # This correct number of are selected with replacement.
#    print("funcGenerateFeatureParameters: Reseting count of samples.")
#    viWhich = sample(1:length(vdExp), size = int_number_features, replace = TRUE)
#    vdExp = vdExp[viWhich]
#    vdMu = vdMu[viWhich]
#    vdSD = vdSD[viWhich]
#    vdPercentZero = vdPercentZero[viWhich]   
#  }

  print("stop funcGenerateFeatureParameters")

  # QC and contraints for percent zero
  # Make sure the percent zero passes the max
  # If there are not enough nonzeros, there is no signal to use.
  # This should be the max. Given there is a certain number of samples that have to have signal
  # The percentage of zeros must allow for those samples not to be zero and so restricts
  # the max the percent zero can be.
  dMaxPercent = 1 - ( iMinNumberSamples / int_number_samples )
  vdPercentZero[ which( vdPercentZero > dMaxPercent ) ] = dMaxPercent
  return( list( exp = vdExp, mu = vdMu, sd = vdSD, PercentZero = vdPercentZero ) )

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
### The vector of expectations for each feature. If not provided one will be generated and vdMu, vdSD, and dPercentZero will be reset
vdMu = NA,
### Vector of Mu parameters for the original exp distribution (means of features) if not supplied, one will be generated
vdPercentZero = NA,
### Vector of percent zero parameters for the original exp distribution (means of features) if not supplied, one will be generated
vdSD = NA,
### Vector of SD parameters for the original exp distribution (means of features) if not supplied, one will be generated
fZeroInflate = TRUE,
### Turns off Zero inflation if FALSE (default TRUE, zero inflation turned on)
lSDRel = NA,
### If vdSD is not given, SD will be generated by a relationship with Exp parameters using the parameters given here
lPercentZeroRel = NA,
### If vdPercentZero is not given, vdPercentZero will be generated by a relationship with Exp using the parameters given here
dBetaGrandSD = c_dBetaGrandSD,
### The beta for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
funcUpdateData = NA,
### Function to update the bug matrix before read depth shuffling. If not provided, will not occur.
fVerbose = FALSE
### Turns on logging (typically generates pdfs)
){
  print(paste("func_generate_random_lognormal_matrix START"))

  # Preallocating for speed
  mat_bugs = matrix(data=NA,nrow=int_number_features,ncol=int_number_samples)
  mtrxWithLeftOver = matrix(data=NA,nrow=int_number_features,ncol=int_number_samples)
  # Get the initial mu vector for generating features.
  lsInitialDistribution = funcGenerateFeatureParameters(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberSamples=iMinNumberSamples, iReadDepth=iReadDepth, vdExp=vdExp, vdMu=vdMu, vdSD=vdSD, vdPercentZero=vdPercentZero, lSDRel, lPercentZeroRel, dBetaGrandSD = dBetaGrandSD, fVerbose=fVerbose)

  pdf("Check2.pdf")
  plot(sort(c_vdExpLog), sort(log(lsInitialDistribution$exp)))
  plot(sort(log(lsInitialDistribution$exp)),main="Exp")
  plot(sort(c_vdSD), sort(log(lsInitialDistribution$sd)), main="SD")
  plot(sort(c_vdPercentZero), sort(lsInitialDistribution$PercentZero),main="PercentZero")

  # Update the Mu, SD and Percent zero bugs and report on distributions
  vdMu = lsInitialDistribution[["mu"]]
  vdSD = lsInitialDistribution[["sd"]]
  vdPercentZero = lsInitialDistribution[["PercentZero"]]
  vdExp = lsInitialDistribution[["exp"]]

  plot(vdExp, funcGetExp(vdMu,vdSD), main="Check vdMu and vdSD pair")

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

  # Make features and assign feature samples to samples giving higher counts to lower read depth samples.
  print("func_generate_random_lognormal_matrix: START Making features")
  for(iReset in 1:int_number_features)
  {
    print(paste("feature",iReset))
    # Create new feature
    lFeatureDetails = funcMakeFeature(dMu=vdMu[iReset], dSD=vdSD[iReset], dPercentZero=vdPercentZero[iReset], iNumberSamples=int_number_samples, iMinNumberCounts=iMinNumberCounts, iMinNumberSamples=iMinNumberSamples, dTruncateThreshold=(c_iTimesSDIsOutlier*vdSD[iReset])+vdMu[iReset], fZeroInflate=fZeroInflate, fVerbose=fVerbose )

    #!# Remove
#    print("lFeatureDetails")
#    print(lFeatureDetails)
    ####

    # Store the left over counts per sample
    vdLeftOver = vdLeftOver + lFeatureDetails$LeftOver

    # Update the matrix with the new feature
    #!# remove leftover
    mat_bugs[iReset,] = lFeatureDetails$Feature
    mtrxWithLeftOver[iReset,] = lFeatureDetails$Feature + lFeatureDetails$LeftOver
  }
  print("func_generate_random_lognormal_matrix: Made features")
  print("vdLeftOver")
  print(vdLeftOver)

  #!# Remove
  plot(vdExp, funcGetRowMetric(mat_bugs,mean), main = "Expected vs Actual Read Depth: Initial")
  plot(vdExp, funcGetRowMetric(mtrxWithLeftOver,mean), main = "Expected vs Actual Read Depth: Initial mtrxWithLeftOver")
  plot(log(vdExp), log(funcGetRowMetric(mat_bugs,mean)), main = "Expected vs Actual Read Depth: Initial logged")
  plot(log(vdExp), log(funcGetRowMetric(mtrxWithLeftOver,mean)), main = "Expected vs Actual Read Depth: Initial mtrxWithLeftOver logged")
  ####

  # Remove any fully zero sample
  lZeroCorrectionResults = funcZeroCorrectMatrix(mtrxData=mat_bugs, vdFeatureMus=vdExp, vdLeftOver=vdLeftOver)
  mat_bugs = lZeroCorrectionResults[["Data"]]
  vdLeftOver = lZeroCorrectionResults[["LeftOver"]]
  print("func_generate_random_lognormal_matrix: Zero corrected")  

  plot(vdExp, funcGetRowMetric(mat_bugs,mean), main = "Expected vs Actual Read Depth: AfterZeroCorrect")

  # Allow a function to be called to manipulate the matrix before updating to the read depth
#  if( !is.na( funcUpdateData ) )
#  {
#    params = funcUpdateData( params )
#  }

  # Shuffle back in removed signal.
#  mat_bugs = funcShuffleMatrix(mtrxData=mat_bugs, vdFeatureMus=vdMu, vdShuffleIn=vdLeftOver)
#  print("func_generate_random_lognormal_matrix: Shuffled")

  plot(vdExp, funcGetRowMetric(mat_bugs,mean), main = "Expected vs Actual Read Depth: After Shuffle")

  # Round to counts
  mat_bugs = funcRoundMatrix(mtrxData=mat_bugs)

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
    print(iReadDepth)
    print("Feature mean summary")
    print(summary(funcGetRowMetric(mat_bugs,mean)))
  }
  dev.off()

  if(fVerbose)
  {
    ## Plot
    barplot(colSums(mat_bugs),main=paste("Read depth mean=",mean(colSums(mat_bugs))), xlab="Samples")
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
      # Controls breaking the following while loop incase a solution can not be found.
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
        print(iSpikeInLoopControl)
        # Get the bug to attempt the association
        # If previously associations have been made
        # liFrozenDataIndices makes the same associations happen here
        # This is so if multiple multipliers are given
        # The different matrices show the differences given increased size of effect
        # Not difference driven by selecting different bugs
        if(!is.null(liFrozeDataIndicies) & length(lviFrozeMetadataIndices)>0)
        {
          iIndexSpikedFeature = liFrozeDataIndicies[[iSpikedBug]]
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
        print(paste("Selected metadata", viSelectedMetadata))
        print(paste("Selected data", iIndexSpikedFeature))

        # lxQCInfo has the slots PASS = Boolean, CommonCounts = vector of integers
        fSpikeInFailed = !lxQCInfo[["PASS"]]

        # Check to see if the looping must end and no success is given
        # Also if the spikein failed, make sure to update the best failed scenario
        iSpikeInLoopControl = iSpikeInLoopControl + 1
        if(fSpikeInFailed)
        {
          # Keep the best failed scenario so far.
          if(!is.null(lxQCInfo[["CommonCounts"]]) && ( sum(lxQCInfo[["CommonCounts"]]>iMinSamples) > lxCurrentBestRun[["Count"]]))
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
            print(paste("While spiking in a relationship between metadata and bug was not able to meet the minimal percentage of spiked-in samples for the relationship. Min sample = ", iMinSamples))
            print(paste("The following spike-in does not pass quality control: Bug=", lxCurrentBestRun$BugIndex," Metadata=",paste(lxCurrentBestRun$MetadataNames,collapse=","),"Multiplier=", multiplier,"Count=", multivariate_parameter))
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
  # If not zero inflated
  if(!fZeroInflate){dPercentZero = 0}

  # Expectation of the feature
  dExpCal = funcGetExp(dMu,dSD)

  # Generate feature
  lFeatureData = func_zero_inflate(log(dMu), dPercentZero, iNumberSamples, log(dSD), dTruncateThreshold)
  dFeature = lFeatureData$Feature
  dLeftOver = lFeatureData$LeftOver

  # Update the distributions to the targeted expectations
  lsShiftResults = funcUpdateDistributionToExpectation( vdFeatures = dFeature, vdLeftOver = dLeftOver, dExp = dExpCal )
  dFeature = lsShiftResults$Feature

#!#  lsForceResults = funcForceMinCountsInMinSamples( vdFeature = dFeature, vdLeftOver = dLeftOver, iMinNumberCounts = iMinNumberCounts, iMinNumberSamples = iMinNumberSamples)
#!#Remove comment  dFeature = lsForceResults$Feature
#!#  dLeftOver = lsForceResults$LeftOver

  # Extra useful measurements, the true and expected means
  dMean = mean(dFeature)

  if(fVerbose)
  {
    plot( dFeature, main = paste("funcMakeFeature","Mean (",round(dMu,2),")",round(dMean,2),"SD(",round(dSD,2),")",round(sd(log(dFeature[which(dFeature>0)])),2),"Exp",round(dExpCal,2)))
  }
  return(list(Feature = dFeature, Exp = dMean, ExpCal = dExpCal, LeftOver = dLeftOver))
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
### Round a sparse matrix. If a value is greater than o but less than 1 then set to zero and then round.
### This keeps the current level of sparsity and does not add more zeros.
mtrxData,
### Matrix of data values to round
vdLeftOver
### Left over signal to use
){
  vdSubCount = intersect( which( mtrxData < 1 ), which( mtrxData > 0 ) )

  mtrxData[ vdSubCount ] = 1
  return( round( mtrxData ) )
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

  # Make left over an int
  viLeftOver = floor(vdShuffleIn)
  # Number of samples
  iNumberSamples = ncol(mtrxData)

  for(iShuffle in 1:iNumberSamples)
  {
    print(paste("Shuffling/Updating Sample", iShuffle))
    print(iShuffle)
    print("viLeftOver[iShuffle]")
    print(viLeftOver[iShuffle])
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
### The value used to define outliers. 
){
  # Get a feature measurement one at a time checking for outliers
  # and redrawing outliers
  vdFeature = rlnorm( iNumberMeasurements, dLogMean, dLogSD )
  vdDifference = rep( 0, iNumberMeasurements )

  # If a value is given to truncate outliers
  if( !is.na( iThreshold ) )
  {
    viOutliers = which( vdFeature > iThreshold )

    # If there are outliers
    # Change the outliers to the mean of the draws before it is drawn (indices less than it's own)
    # and optionally shuffle back in the differences between the outliers and the max
    if( length( viOutliers ) )
    {
      print("dMeanDraws")
      print(paste("dLogMean",round(exp(dLogMean),2)))
      print(paste("dLogSD",round(exp(dLogSD),2)))
      vdMeanDraws = c()
      for( iIndex in viOutliers )
      {
        dMeanDraws = mean( vdFeature[ sample( 1:iIndex, 2, replace = TRUE ) ] )
        vdMeanDraws= c(vdMeanDraws, dMeanDraws)
        vdDifference[ iIndex ] = vdFeature[ iIndex ] - dMeanDraws
        vdFeature[ iIndex ] = dMeanDraws
      }
      print(sort(vdMeanDraws))
    }
  }
  #Truncate negatives to zero
  vdFeature[ vdFeature < 0 ] = 0

  print("vdFeature end trunc")
  print(vdFeature)
  return(list(Feature=vdFeature, LeftOver=vdDifference))
### vdFeature: The signal (optionally truncated, log normal distribution)
### vdDifference: Left over signal that was removed from the original vdFeature
}


funcUpdateDistributionToExpectation = function(
### Updates a distribution to the mean while keeping the shape.
vdFeatures,
### Distribution of counts
vdLeftOver,
### Left over signal, this is not used but is updated. So as many counts as needed are used to update the distribution and are removed from vdLeftOver if there.
dExp
### The expectation to which to update the distribution
){
  print("funcUpdateDistributionToExpectation START")
  # Used to ignore zeros in these calculations
  viNonZeroIndices = which( vdFeatures > 0 )
  iLengthNoZeros = length( viNonZeroIndices )

  # Get the amount ofchange in signal needed
  dDifference = dExp - mean( vdFeatures[ viNonZeroIndices ] )

  # If signal is needed to be changed then update or remove signal by count
  if( abs( dDifference ) > 0 )
  {
    # Number of counts to shift
    dCounts = abs( floor( dDifference * iLengthNoZeros ) )
    # Add to the distributio and shift it up
    if(dDifference > 0)
    {
      vdUpdateIndices = sample( viNonZeroIndices, dCounts, replace = TRUE, prob = vdFeatures[ viNonZeroIndices ] / sum( vdFeatures[ viNonZeroIndices ] ) )
      for( iIndex in vdUpdateIndices )
      {
        vdFeatures[ iIndex ] = vdFeatures[ iIndex ] + 1
      }
    } else if( dDifference < 0 )
    {
      # Remove fom distribution
      for( iIndex in 1:dCounts )
      {
        viGreaterThan1 = which( vdFeatures > 1 )
        if( length( viGreaterThan1 ) > 1 )
        {
          iUpdateIndex = sample( viGreaterThan1, 1, prob = vdFeatures[ viGreaterThan1 ] / sum( vdFeatures[ viGreaterThan1 ] ) )
          vdFeatures[ iUpdateIndex ] = vdFeatures[ iUpdateIndex ] - 1
        }
      }
    }
    #!# update leftover
  }
  print("funcUpdateDistributionToExpectation STOP")
  return( list( Feature = vdFeatures, LeftOver = vdLeftOver ) )
}

funcUpdateDistributionToSum = function(
### Update the distribution to the sum (read depth) requested.
### If the current read depth is not as large as needed sample the difference using
### a multinomial model based on the current taxa means (values) (as proportions)
### If the current read depth is larger than needed
### Then remove values in the same manner as described above for adding values.
### Except that a count of 1 can not be removed. This would throw off the sparsity.
vdDistribution,
### A vector of numeric data that will be forced to sum to iTargetSum (+- less than 1)
iTargetSum
### The target sum to update the distribution to
){
  # Update vdExp to target
  dCurReadDepthDifference = iTargetSum - sum( vdDistribution )

  # Multinomial probability
  vdProbabilities = vdDistribution / sum( vdDistribution )

  # Need more counts for the sum
  if(dCurReadDepthDifference > 0)
  {
    # Use a multinomial model to select indices
    viSamplingIndices = sample( x = 1:length( vdDistribution ), size = abs( dCurReadDepthDifference ), replace = TRUE, prob = vdProbabilities )
    # Update indices
    for(iIndex in viSamplingIndices)
    {
      vdDistribution[iIndex] = vdDistribution[iIndex]+1
    }
  # Need less counts for sum
  } else if( dCurReadDepthDifference < 0 ) {
   for( iIndex in 1:round( abs( dCurReadDepthDifference ) ) )
    {
      # Find those entries that are greater than one (do not want to make zeros here)
      viGreaterThanOne = which( vdDistribution > 1 )
      if(length( viGreaterThanOne ) > 0)
      {
        # If there is only one entry that is not 1 then select it otherwise sample.
        iGreaterThanOne = viGreaterThanOne
        if( length( iGreaterThanOne ) > 1 ){ iGreaterThanOne = sample( x = iGreaterThanOne, size = 1, prob = vdProbabilities[viGreaterThanOne] ) }
        vdDistribution[ iGreaterThanOne ] = vdDistribution[ iGreaterThanOne ] - 1
      }
    }
  }
  return( vdDistribution )
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
dLogMean,
### Mean of the distribution (logged)
dPercentZeroInflated,
### Percentage of return which is zero
int_number_samples,
### The number of samples to create
dLogSD,
### The sd of the distribution (logged)
iThreshold = NA
### The threshold for outliers
){
  # Get feature given distribution parameters
  lFeatureData = funcTruncatedRLNorm( int_number_samples, dLogMean, dLogSD, iThreshold = iThreshold )
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
#dSDBeta = c_dSDBeta
#dBetaZero = c_dBetaZero
#dGrandBeta = c_dBetaGrandSD

print("Parameters BEFORE Calibration File")
print(paste("Length exp",NA,"Length vdMu", NA, "length vdSD", NA, "length vdPercentZero", NA, "Read depth", iReadDepth))

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
  iReadDepth = lsFit[["dAverageReadDepth"]]
  int_number_features = lsFit[["iFeatureCount"]]
}
print("Parameters AFTER Calibration File (if no calibration file is used, defaults are shown)")
print(paste("Length exp",length(vdExp),"Length vdMu", length(vdMu), "length vdSD", length(vdSD), "length vdPercentZero", length(vdPercentZero), "Read depth", iReadDepth, "Feature Count", int_number_features))

mat_random_lognormal_bugs = func_generate_random_lognormal_matrix(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberCounts=dMinOccurenceCount, iMinNumberSamples=dMinOccurenceSample, iReadDepth=iReadDepth, vdExp=vdExp, vdMu=vdMu, vdPercentZero=vdPercentZero, vdSD=vdSD, fZeroInflate=fZeroInflate, lSDRel=list(BetaSD=c_dSDBeta, InterceptSD=c_dSDIntercept), lPercentZeroRel = list(InterceptZero=c_dInterceptZero, BetaZero=c_dBetaZero, Beta2Zero=c_dBeta2Zero), dBetaGrandSD=c_dBetaGrandSD, fVerbose=fVerbose)

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
    print("how_many_multivariates")
    print(how_many_multivariates)
    print("mult")
    print(mult)
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
lsBugBugInfo = func_generate_bug_bug_spiking_matrix(int_number_features=int_number_features, int_number_samples=int_number_samples, iMinNumberCounts=dMinOccurenceCount, iMinNumberSamples=dMinOccurenceSample, iReadDepth=iReadDepth, vdExp=vdExp, vdMu=vdMu, vdPercentZero=vdPercentZero, vdSD=vdSD, fZeroInflate=fZeroInflate, lSDRel=list(BetaSD=c_dSDBeta, InterceptSD=c_dSDIntercept), lPercentZeroRel = list(InterceptZero=c_dInterceptZero, BetaZero=c_dBetaZero, Beta2Zero=c_dBeta2Zero), dBetaGrandSD=c_dBetaGrandSD, fVerbose=fVerbose, dVarScale=dVarScale,iNumAssociations=iNumAssociations,iMaxNumberCorrDomainBugs=iMaxNumberCorrDomainBugs)
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
