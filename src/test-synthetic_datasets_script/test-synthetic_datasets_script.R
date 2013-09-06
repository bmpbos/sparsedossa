c_strDir <- file.path(getwd( ),"..")

source(file.path(c_strDir,"synthetic_datasets_script.R"))

########## General Controls
# Control the percision used for rounding in the test results.
c_iRoundingPercision = 2
c_fVerbose = TRUE
sCalibrateFile = "Calibrate-RLNorm.tsv"
sCalibrateSparseFile = "Calibrate-Sparse.tsv"
##########

########## Start testing


### Test that the script runs to completion
#!#

### funcBoxPlotOutliers
### Plotting function.


### funcCalibrateRLNormToMicrobiome
### TODO(Timothy): Make a file and calculate the correct relationships and make sure that is found here.


### funcCycleVectorIndices
context( "Test funcCycleVectorIndices" )
viExpected1 <- c( 1, 2, 1, 2, 1, 2 )
viExpected2 <- c( 1, 2, 3, 1, 2, 3, 1, 2 )
viExpected3 <- 1:10
viExpected4 <- 1
viExpected5 <- rep( 1, 20 )
test_that( "Test that funcCycleVectorIndices cycles through indices correctly.", {
  expect_equal( viExpected1, funcCycleVectorIndices( 6, 2 ) )
  expect_equal( viExpected2, funcCycleVectorIndices( 8, 3 ) )
  expect_equal( viExpected3, funcCycleVectorIndices( 10, 20 ) )
  expect_equal( viExpected4, funcCycleVectorIndices( 1, 100 ) )
  expect_equal( viExpected5, funcCycleVectorIndices( 20, 1 ) )
})


### funcEstimateFeatureSD
context("Test funcEstimateFeatureSD")
test_that( "Test that funcEstimateFeatureSD gives the correct answer.", {
  expect_equal( 0, funcEstimateFeatureSD( 0, 12834 ) )
  expect_equal( 1, funcEstimateFeatureSD( 1, 1 ) )
  expect_equal( 1.5, funcEstimateFeatureSD( .5, 3 ) )
  expect_equal( 1010, funcEstimateFeatureSD( 101, 10 ) )
  expect_equal( -24, funcEstimateFeatureSD( -3, 8 ) )
})


### funcEstimateGrandSD
context("Test funcEstimateGrandSD")
test_that( "Test that funcEstimateGrandSD gives the correct answer.", {
  expect_equal( 0, funcEstimateGrandSD( 0, 12834 ) )
  expect_equal( 1, funcEstimateGrandSD( 1, 1 ) )
  expect_equal( 1.5, funcEstimateGrandSD( .5, 3 ) )
  expect_equal( 1010, funcEstimateGrandSD( 101, 10 ) )
  expect_equal( -24, funcEstimateGrandSD( -3, 8 ) )
})


### funcEstimatePercentZero
context("Test funcEstimatePercentZero")
test_that( "Test that funcEstimatePercentZero gives the correct answer.", {
  expect_equal( 0, funcEstimatePercentZero( 0, 12834 ) )
  expect_equal( 1, funcEstimatePercentZero( 1, 1 ) )
  expect_equal( 1.5, funcEstimatePercentZero( .5, 3 ) )
  expect_equal( 1010, funcEstimatePercentZero( 101, 10 ) )
  expect_equal( -24, funcEstimatePercentZero( -3, 8 ) )
})


### funcGenerateExpVectorParameters
context( "Test funcGenerateExpVectorParameters" )
dMu <- 2
dSD <- 3
iCount = 1000000
dLogMu <- log( dMu )
dLogSD <- log( dSD )
dExp <- funcGetExp( dMu, dSD )
dRelationship <- dLogSD / dExp
dReadDepth <- dExp * iCount
vdDistribution <- rlnorm( iCount, dLogMu, dLogSD )
lResult = funcGenerateExpVectorParameters(vdDistribution)
testthat( "Test that funcGenerateExpVectorParameters gives the correct answer for good case 1.", {
  expect_equal( round( dLogMu, 2 ), round( lResult$GrandLogMu, 2 ) )
  expect_equal( round( dExp, 2), round( lResult$GrandExp, 2 ) )
  expect_equal( round( dLogSD, 2), round( lResult$GrandLogSD, 2 ) )
  expect_true( dReadDepth - ( dReadDepth * 0.1 ) < lResult$dReadDepth )
  expect_true( dReadDepth + ( dReadDepth * 0.1 ) > lResult$dReadDepth )
  expect_equal( round( dRelationship, 2 ), round( lResult$GrandLogSDBeta, 2 ) )
})

dMu <- 5
dSD <- 4
iCount = 1000000
dLogMu <- log( dMu )
dLogSD <- log( dSD )
dExp <- funcGetExp( dMu, dSD )
dRelationship <- dLogSD / dExp
dReadDepth <- dExp * iCount
vdDistribution <- rlnorm( iCount, dLogMu, dLogSD )
lResult = funcGenerateExpVectorParameters(vdDistribution)
testthat( "Test that funcGenerateExpVectorParameters gives the correct answer for good case 2.", {
  expect_equal( round( dLogMu, 2 ), round( lResult$GrandLogMu, 2 ) )
  expect_equal( round( dExp, 2), round( lResult$GrandExp, 2 ) )
  expect_equal( round( dLogSD, 2), round( lResult$GrandLogSD, 2 ) )
  expect_true( dReadDepth - ( dReadDepth * 0.1 ) < lResult$dReadDepth )
  expect_true( dReadDepth + ( dReadDepth * 0.1 ) > lResult$dReadDepth )
  expect_equal( round( dRelationship, 2 ), round( lResult$GrandLogSDBeta, 2 ) )
})


### func_generate_bug_bug_spiking_matrix
#!#


### func_generate_lefse_matrices
#!#


### funcGenerateFeatureParameters
context( "Test funcGenerateFeatureParameters" )
iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- NA
vdMu <- NA
vdSD <- NA
vdPercentZero <- NA
dBetaSD <- 1
dBetaZero <- 0
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )

lResult = funcGenerateFeatureParameters( iNumberFeatures, iNumberSamples, iMinNumberSamples, iReadDepth, vdExp, vdMu, vdSD, vdPercentZero, dBetaSD, dBetaZero, dBetaGrandSD )
testthat( "funcGenerateFeatureParameters: Generating all vectors." {
  expect_equal( round( iReadDepth / iNumberFeatures, 2 ), round( mean( lResult$exp ), 2 ) )
  expect_equal( iReadDepth, sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) )
  expect_equal( round( lResult$exp, 2), round( lResult$sd, 2) )
  expect_equal( 0, sum( lResult$PercentZero ) )
})

iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdMu <- NA
vdSD <- NA
vdPercentZero <- NA
dBetaSD <- 1
dBetaZero <- 0
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )

lResult = funcGenerateFeatureParameters( iNumberFeatures, iNumberSamples, iMinNumberSamples, iReadDepth, vdExp, vdMu, vdSD, vdPercentZero, dBetaSD, dBetaZero, dBetaGrandSD )
testthat( "funcGenerateFeatureParameters: Starting with Mu and SD vector." {
  expect_equal( round( iReadDepth / iNumberFeatures, 2 ), round( mean( lResult$exp ), 2 ) )
  expect_equal( iReadDepth, sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) )
  expect_equal( round( lResult$exp, 2), round( lResult$sd, 2) )
  expect_equal( 0, sum( lResult$PercentZero ) )
  expect_equal( round( lResult$mu, 2), round( c(0.05285117, 0.04323700, 0.09918717, 0.04046334, 0.04762098, 0.07139393, 0.05309049, 0.08947141, 0.06587467, 0.04046334), 2 ) )
  expect_equal( round( lResult$sd, 2), round( c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829), 2 ) )
  expect_equal( round( lResult$PercentZero, 2), round( rep( 1, length( lResult$exp ) ), 2 ) )
})

iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdMu <- c(0.05285117, 0.04323700, 0.09918717, 0.04046334, 0.04762098, 0.07139393, 0.05309049, 0.08947141, 0.06587467, 0.04046334)
vdSD <- NA
vdPercentZero <- NA
dBetaSD <- 1
dBetaZero <- 0
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )

lResult = funcGenerateFeatureParameters( iNumberFeatures, iNumberSamples, iMinNumberSamples, iReadDepth, vdExp, vdMu, vdSD, vdPercentZero, dBetaSD, dBetaZero, dBetaGrandSD )
testthat( "funcGenerateFeatureParameters: Starting with Mu and SD and Mu vector." {
  expect_equal( round( iReadDepth / iNumberFeatures, 2 ), round( mean( lResult$exp ), 2 ) )
  expect_equal( iReadDepth, sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) )
  expect_equal( round( lResult$exp, 2), round( lResult$sd, 2) )
  expect_equal( 0, sum( lResult$PercentZero ) )
  expect_equal( round( lResult$sd, 2), round( c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829), 2 ) )
  expect_equal( round( lResult$PercentZero, 2), round( rep( 1, length( lResult$exp ) ), 2 ) )
})

iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdMu <- c(0.05285117, 0.04323700, 0.09918717, 0.04046334, 0.04762098, 0.07139393, 0.05309049, 0.08947141, 0.06587467, 0.04046334)
vdSD <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdPercentZero <- NA
dBetaSD <- 1
dBetaZero <- 0
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )

lResult = funcGenerateFeatureParameters( iNumberFeatures, iNumberSamples, iMinNumberSamples, iReadDepth, vdExp, vdMu, vdSD, vdPercentZero, dBetaSD, dBetaZero, dBetaGrandSD )
testthat( "funcGenerateFeatureParameters: Starting with Mu and SD and Mu vector." {
  expect_equal( round( iReadDepth / iNumberFeatures, 2 ), round( mean( lResult$exp ), 2 ) )
  expect_equal( iReadDepth, sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) )
  expect_equal( round( lResult$exp, 2), round( lResult$sd, 2) )
  expect_equal( 0, sum( lResult$PercentZero ) )
  expect_equal( round( lResult$sd, 2), round( c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829), 2 ) )
  expect_equal( round( lResult$PercentZero, 2), round( rep( 1, length( lResult$exp ) ), 2 ) )
})

iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdMu <- c(0.05285117, 0.04323700, 0.09918717, 0.04046334, 0.04762098, 0.07139393, 0.05309049, 0.08947141, 0.06587467, 0.04046334)
vdSD <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdPercentZero <- rep( 1, length( lResult$exp ) )
dBetaSD <- 1
dBetaZero <- 0
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )

lResult = funcGenerateFeatureParameters( iNumberFeatures, iNumberSamples, iMinNumberSamples, iReadDepth, vdExp, vdMu, vdSD, vdPercentZero, dBetaSD, dBetaZero, dBetaGrandSD )
testthat( "funcGenerateFeatureParameters: Starting with Mu and SD and Mu and PZ vector." {
  expect_equal( round( iReadDepth / iNumberFeatures, 2 ), round( mean( lResult$exp ), 2 ) )
  expect_equal( iReadDepth, sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) )
  expect_equal( round( lResult$exp, 2), round( lResult$sd, 2) )
  expect_equal( 0, sum( lResult$PercentZero ) )
})

iNumberFeatures <- 100
lResult = funcGenerateFeatureParameters( iNumberFeatures, iNumberSamples, iMinNumberSamples, iReadDepth, vdExp, vdMu, vdSD, vdPercentZero, dBetaSD, dBetaZero, dBetaGrandSD )
testthat( "funcGenerateFeatureParameters: Starting with all vectors but wrong size." {
  expect_equal(,)
})

testthat( "funcGenerateFeatureParameters: Starting with min number really high, min samples and sd = 1." {
  expect_equal(iNumberFeatures,length(lResult$mu))
  expect_equal(iNumberFeatures,length(lResult$exp))
  expect_equal(iNumberFeatures,length(lResult$sd))
  expect_equal(iNumberFeatures,length(lResult$PercentZero))
})

iNumberFeatures <- 1000
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 36500
vdExp <- NA
vdMu <- NA
vdSD <- NA
vdPercentZero <- NA
dBetaSD <- 1
dBetaZero <- 0
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )

testthat( "funcGenerateFeatureParameters: Read depth 36500." {
  expect_equal(36500,sum(lResult$exp))
})


### func_generate_metadata
#!#


### func_generate_random_lognormal_matrix
#!#


### func_generate_random_lognormal_with_outliers
context("Test func_generate_random_lognormal_with_outliers")
#!#


### func_generate_random_lognormal_with_multivariate_spikes
context("Test func_generate_random_lognormal_with_multivariate_spikes")
#!#


### funcGetExp
context("Test funcGetExp")
dSDOne = round(3,4)
dMuOne = round(4,4)
dExpOne = mean(rlnorm(1000000,log(dMuOne),log(dSDOne)))
dSDTwo = round(2,4)
dMuTwo = round(3,4)
dExpTwo = mean(rlnorm(1000000,log(dMuTwo),log(dSDTwo)))
dSDThree = round(1,4)
dMuThree = round(1,4)
dExpThree = mean(rlnorm(1000000,log(dMuThree),log(dSDThree)))
test_that("Test that the expected SD is calculated correctly.",{
  expect_equal(round(dExpOne,1),round(funcGetExp(dMuOne,dSDOne),1))
  expect_equal(round(dExpTwo,1),round(funcGetExp(dMuTwo,dSDTwo),1))
  expect_equal(round(dExpThree,1),round(funcGetExp(dMuThree,dSDThree),1))
})


### funcGetMu
context("Test funcGetMu")
dSDOne = round(3,4)
dMuOne = round(4,4)
dExpOne = mean(rlnorm(1000000,log(dMuOne),log(dSDOne)))
dSDTwo = round(2,4)
dMuTwo = round(3,4)
dExpTwo = mean(rlnorm(1000000,log(dMuTwo),log(dSDTwo)))
dSDThree = round(1,4)
dMuThree = round(1,4)
dExpThree = mean(rlnorm(1000000,log(dMuThree),log(dSDThree)))
test_that("Test that the expected SD is calculated correctly.",{
  expect_equal(dMuOne,round(funcGetMu(dExpOne,dSDOne),1))
  expect_equal(dMuTwo,round(funcGetMu(dExpTwo,dSDTwo),1))
  expect_equal(dMuThree,round(funcGetMu(dExpThree,dSDThree),1))
})


### funcGetParamsForReadDepth
context("Test funcGetParamsForReadDepth")
dMu <- 2
dSD <- 3
dExp <- funcGetExp( dMu, dSD )
iCount <- 1000000
vdDistribution <- rlnorm( iCount, log( dMu ), log( dSD ) )
lResult <- funcGetParamsForReadDepth( mean( vdDistribution ), log(dSD) / dExp )

test_that( "Test that the parameters measured by funcGetParamsForExpectation are off less than 1 decimal, test 1." ,{
  expect_equal( round( log(dMu) ), round( lResult$dLogMu ) )
  expect_equal( round( log(dSD) ), round( lResult$dLogSD ) )
})

dMu <- 4
dSD <- 2
dExp <- funcGetExp( dMu, dSD )
iCount <- 1000000
vdDistribution <- rlnorm( iCount, log( dMu ), log( dSD ) )
lResult <- funcGetParamsForReadDepth( mean( vdDistribution ), log(dSD) / dExp )

test_that( "Test that the parameters measured by funcGetParamsForExpectation are off less than 1 decimal, test 1." ,{
  expect_equal( round( log(dMu) ), round( lResult$dLogMu ) )
  expect_equal( round( log(dSD) ), round( lResult$dLogSD ) )
})


### funcGetRowMetric
viOne <- 1:10
viTwo <- c( 1, 5, 2, 8, 3, 4, 6, 9, 12, 2 )
viThree <- c( 1, 24, 32, 5, 21, 2, 4, 6, 8, 67 )
vdAnswer3 <- c( 5.5, 5.2, 17 )
test_that( "Test that the mean is calculated correctly for vectors or matrices." ,{
  expect_equal( funcGetRowMetric( viOne, mean ), 5.5 )
  expect_equal( funcGetRowMetric( matrix( viOne, ncol = length( viOne ), byrow = TRUE ), mean ), 5.5 )
  expect_equal( funcGetRowMetric( matrix( c( viOne, viTwo, viThree ), ncol = length( viOne ), byrow = TRUE ), mean), vdAnswer3 )
})


### Test funcGetSD
context("Test funcGetSD")
viOne = 1:10
viTwo = c(1,5,2,8,3,4,6,9,12,2)
viThree = c(1,24,32,5,21,2,4,6,8,67)
dSDOne = 3.02765
dSDTwo = 3.552777
dSDThree = 20.51016
test_that("Test that the SD is calculated correctly for vectors or matrices.",{
  expect_equal(round(funcGetSD(viOne),c_iRoundingPercision), round(dSDOne,c_iRoundingPercision))
  expect_equal(round(funcGetSD(matrix(c(viOne),ncol=length(viOne),byrow=TRUE)),c_iRoundingPercision), round(dSDOne,c_iRoundingPercision))
  expect_equal(round(funcGetSD(matrix(c(viOne, viTwo, viThree),ncol=length(viOne),byrow=TRUE)),c_iRoundingPercision), round(c(dSDOne ,dSDTwo, dSDThree),c_iRoundingPercision))
})


### funcIsFactorMetadataValid
context("Test funcIsFactorMetadataValid")
vxMetadataOne = c(1,2,2,1,1,1,2,2,1,2,1,2)
vxMetadataTwo = c(1,2,3,1,1,1,3,3,2,2,3,2,1,3)
test_that("Test that minimal levels of data are detected correctly, discontinuous data only.",{
  expect_equal(funcIsFactorMetadataValid(vxMetadataOne,1), TRUE)
  expect_equal(funcIsFactorMetadataValid(vxMetadataOne,8), FALSE)
  expect_equal(funcIsFactorMetadataValid(vxMetadataTwo,4), TRUE)
  expect_equal(funcIsFactorMetadataValid(vxMetadataTwo,5), FALSE)
})


#!#
### funcMakeFeature
context("Test funcMakeFeature for min counts")
# Test for minimum count in a sample
dMuOne = 1
dSDOne = 1
dPercentZeroOne = 0.0
iNumberSamplesOne = 1000
iMinNumberCountsOne = 1000
iMinNumberSamplesOne = iNumberSamplesOne
dTotalReadDepth = 1

iMinNumberCountsTwo = 50
iMinNumberCountsThree = 0

test_that("Test that the expected number of samples are given the minimum value.",{
  expect_equal(length(which(funcMakeFeature(dMuOne,dSDOne,dPercentZeroOne,iNumberSamplesOne,iMinNumberCountsOne,iMinNumberSamplesOne, dTotalReadDepth)$Feature==iMinNumberCountsOne)),iMinNumberSamplesOne)
  expect_equal(length(which(funcMakeFeature(dMuOne,dSDOne,dPercentZeroOne,iNumberSamplesOne,iMinNumberCountsTwo,iMinNumberSamplesOne, dTotalReadDepth)$Feature==iMinNumberCountsTwo)),iMinNumberSamplesOne)
  expect_equal(length(which(funcMakeFeature(dMuOne,dSDOne,dPercentZeroOne,iNumberSamplesOne,iMinNumberCountsThree,iMinNumberSamplesOne, dTotalReadDepth)$Feature==iMinNumberCountsThree)),0)
})

context("Test funcMakeFeature for parameters outside of min number counts")
# Test for minimum count in a sample
dMuOne = 10
dSDOne = 2
dPercentZeroOne = 0.4
iNumberSamplesOne = 1000
iMinNumberCounts = 0
iMinNumberSamples = 0

dMuTwo = 10
dSDTwo = 2
dPercentZeroTwo = 0.0
iNumberSamplesTwo = 1000

dTotalReadDepth = 1

vdFeatureOne = funcMakeFeature(dMuOne,dSDOne,dPercentZeroOne,iNumberSamplesOne,iMinNumberCounts,iMinNumberSamples)$Feature
viNotZeroOne = which(vdFeatureOne!=0)
vdFeatureTwo = funcMakeFeature(dMuTwo,dSDTwo,dPercentZeroTwo,iNumberSamplesTwo,iMinNumberCounts,iMinNumberSamples)$Feature
viNotZeroTwo = which(vdFeatureTwo!=0)

dRealMuOne = mean(log(vdFeatureOne[viNotZeroOne]))
dRealMuTwo = mean(log(vdFeatureTwo[viNotZeroTwo]))
dRealSDOne = sd(log(vdFeatureOne[viNotZeroOne]))
dRealSDTwo = sd(log(vdFeatureTwo[viNotZeroTwo]))
dRealPercentZeroOne = length(which(vdFeatureOne==0))/length(vdFeatureOne)
dRealPercentZeroTwo = length(which(vdFeatureTwo==0))/length(vdFeatureTwo)

test_that("Test that the expected number of samples are given the minimum value.",{
  expect_equal(round(dRealMuOne),round(log(dMuOne)))
  expect_equal(round(dRealMuTwo),round(log(dMuTwo)))
  expect_equal(round(dRealSDOne),round(log(dSDOne)))
  expect_equal(round(dRealSDTwo),round(log(dSDTwo)))
  expect_equal(round(dRealPercentZeroOne),round(dPercentZeroOne))
  expect_equal(round(dRealPercentZeroTwo),round(dPercentZeroTwo))
})


### funcNormalizeMicrobiome
context("Test funcNormalizeMicrobiome")
mtrxCountsOne = matrix(c(1,1,1,2,2,2,3,3,3,4,4,4), byrow=TRUE, nrow=4, ncol=3)
mtrxRelAbOne = matrix(c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4), byrow=TRUE, nrow=4, ncol=3)
mtrxCountsTwo = matrix(c(1,1,1,2,2,2,3,3,3,4,4,4), nrow=4, ncol=3)
mtrxRelAbTwo = matrix(c(1/5,1/5,1/5,2/5,2/10,2/10,3/10,3/10,3/15,4/15,4/15,4/15), nrow=4, ncol=3)
mtrxCountsZero = matrix(rep(0,4*3), byrow=TRUE, nrow=4, ncol=3)
mtrxRelAbZero = matrix(rep(0,4*3), byrow=TRUE, nrow=4, ncol=3)
mtrxCountsFour = matrix(c(1,0,1,0,2,0,3,0,3,0,4,0), byrow=TRUE, nrow=4, ncol=3)
mtrxRelAbFour = matrix(c(1/4,0/6,1/4,0/4,2/6,0/4,3/4,0/6,3/4,0/4,4/6,0/4), byrow=TRUE, nrow=4, ncol=3)
test_that("Test that normalization occurs by column and when zeros are present.",{
  expect_equal(mtrxRelAbOne,funcNormalizeMicrobiome(mtrxCountsOne))
  expect_equal(mtrxRelAbTwo,funcNormalizeMicrobiome(mtrxCountsTwo))
  expect_equal(mtrxRelAbZero,funcNormalizeMicrobiome(mtrxCountsZero))
  expect_equal(mtrxRelAbFour,funcNormalizeMicrobiome(mtrxCountsFour))
})


### funcNumericIsInt
context("Test funcNumericIsInt")
test_that("Test that integers are correctly identified.",{
  expect_equal(funcNumericIsInt(1), TRUE)
  expect_equal(funcNumericIsInt(2.0), TRUE)
  expect_equal(funcNumericIsInt(1.5), FALSE)
  expect_equal(funcNumericIsInt(0.00005), FALSE)
})


### funcPlotSpikeins
### Plotting function


#!#
### funcPrepareMetadata
context("Test funcPrepareMetadata with dummying")
vstrMetadataOne = c(1,2,2,1,1,1,2,2,1,2,1,2)
vstrMetadataTwo = c(1,2,3,1,1,3,2,3,1,2,1,2)
vstrMetadataThree = c(1,2,3,4,1,3,2,4,1,2,1,2)

lResult1 = funcPrepareMetadata(c(3),vstrMetadataOne,TRUE)
lResult2 = funcPrepareMetadata(c(8),vstrMetadataOne,TRUE)
lResult3 = funcPrepareMetadata(c(3,13,23),matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE),TRUE)
lResult4 = funcPrepareMetadata(c(1,2,3),matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE),TRUE)

strLevelSelected1 = strsplit(lResult1[["names"]],"_")[[1]]
strLevelSelected1 = strLevelSelected1[length(strLevelSelected1)]
strLevelSelected2 = strsplit(lResult2[["names"]],"_")[[1]]
strLevelSelected2 = strLevelSelected2[length(strLevelSelected2)]
strLevelSelected3 = unlist(lapply(strsplit(lResult3[["names"]],"_"), function(x) x[[length(x)]]))
strLevelSelected4 = unlist(lapply(strsplit(lResult4[["names"]],"_"), function(x) x[[length(x)]]))

vxMetadata1 = matrix(rep(1, length(vstrMetadataOne)),nrow=1)
vxMetadata2 = matrix(rep(1, length(vstrMetadataOne)),nrow=1)

vxMetadata1[which(vstrMetadataOne==strLevelSelected1)] = 2
vstrNames1 = c(paste(c_strMetadata,3,"_",c_strLevel,"_",strLevelSelected1,sep=""))
vxMetadata2[which(vstrMetadataOne==strLevelSelected2)] = 2
vstrNames2 = c(paste(c_strMetadata,8,"_",c_strLevel,"_",strLevelSelected2,sep=""))

vxMetadata3 = matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE)
iLevelIndices = which(vxMetadata3==strLevelSelected3)
liDim = dim(vxMetadata3)
vxMetadata3[setdiff(c(1:(liDim[1]*liDim[2])),iLevelIndices)] = 1
vxMetadata3[iLevelIndices] = 2
vstrNames3 = c(paste(c_strMetadata,3,"_",c_strLevel,"_", strLevelSelected3[[1]],sep=""),paste(c_strMetadata,13,"_",c_strLevel,"_", strLevelSelected3[[2]],sep=""),paste(c_strMetadata,23,"_",c_strLevel,"_", strLevelSelected3[[3]],sep=""))

vxMetadata4 = matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE)
iLevelIndices = which(vxMetadata4==strLevelSelected4)
liDim = dim(vxMetadata4)
vxMetadata4[setdiff(c(1:(liDim[1]*liDim[2])),iLevelIndices)] = 1
vxMetadata4[iLevelIndices] = 2
vstrNames4 = c(paste(c_strMetadata,1,"_",c_strLevel,"_", strLevelSelected4[[1]],sep=""),paste(c_strMetadata,2,"_",c_strLevel,"_",strLevelSelected4[[2]],sep=""),paste(c_strMetadata,3,"_",c_strLevel,"_",strLevelSelected4[[3]],sep=""))

test_that("Test that factor data of different varying levels are correctly dummied and their names with levels are recoded.",{
  expect_equal(sort(lResult1[["names"]]), sort(vstrNames1))
  expect_equal(lResult1[["metadata"]], vxMetadata1)
  expect_equal(sort(lResult2[["names"]]), sort(vstrNames2))
  expect_equal(lResult2[["metadata"]], vxMetadata2)
  expect_equal(sort(lResult3[["names"]]), sort(vstrNames3))
  expect_equal(lResult3[["metadata"]], vxMetadata3)
  expect_equal(sort(lResult4[["names"]]), sort(vstrNames4))
  expect_equal(lResult4[["metadata"]], vxMetadata4)
})

context("Test funcPrepareMetadata with OUT dummying")
vstrMetadataOne = c(1,2,2,1,1,1,2,2,1,2,1,2)
vstrMetadataTwo = c(1,2,3,1,1,3,2,3,1,2,1,2)
vstrMetadataThree = c(1,2,3,4,1,3,2,4,1,2,1,2)

lResult1 = funcPrepareMetadata(c(3),vstrMetadataOne,FALSE)
lResult2 = funcPrepareMetadata(c(8),vstrMetadataOne,FALSE)
lResult3 = funcPrepareMetadata(c(3,13,23),matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE),FALSE)
lResult4 = funcPrepareMetadata(c(1,2,3),matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE),FALSE)

vstrNames1 = paste(c_strMetadata,c(3),sep="")
vstrNames2 = paste(c_strMetadata,c(8),sep="")
vstrNames3 = paste(c_strMetadata,c(3,13,23),sep="")
vstrNames4 = paste(c_strMetadata,c(1,2,3),sep="")

vxMetadata1 = vstrMetadataOne
vxMetadata2 = vxMetadata1
vxMetadata3 = matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE)
vxMetadata4 = matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE)

test_that("Test that factor data of different varying levels are not dummied and their names with levels are not recorded.",{
  expect_equal(lResult1[["names"]], vstrNames1)
  expect_equal(lResult1[["metadata"]], vxMetadata1)
  expect_equal(lResult2[["names"]], vstrNames2)
  expect_equal(lResult2[["metadata"]], vxMetadata2)
  expect_equal(lResult3[["names"]], vstrNames3)
  expect_equal(lResult3[["metadata"]], vxMetadata3)
  expect_equal(lResult4[["names"]], vstrNames4)
  expect_equal(lResult4[["metadata"]], vxMetadata4)
})


### funcQCSpikin
context("Test funcQCSpikin")
viMetadata1 = c(2,1,2,2,1,1,2,2,1,1)
viMetadata2 = c(2,1,1,1,2,2,1,1,1,1)
vdBug1 = c(0.1,0.3,0.2,0.6,0.8,0.9,0.3,0.5,0.6,0.4)
vdBug2 = c(0.1,0.0,0.2,0.0,0.8,0.0,0.3,0.0,0.6,0.0)
vdBug3 = c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
test_that("Test that dummied data are correctly checked.",{
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 5)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 3)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 0)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 6)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 130)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug2, 2)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug2, 4)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug3, 1)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 3)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 2)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 0)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 4)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 100)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug2, 2)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug2, 3)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug3, 1)$PASS, FALSE)
})

viMetadata3 = c(2,1,2,2,2,2,2,2,1,2)
viMetadata4 = rep(1,length(viMetadata1))
mtrxMatrixOf2 = matrix(c(viMetadata1,viMetadata2), byrow=TRUE, nrow=2)
mtrxMatrixOf3 = matrix(c(viMetadata3, viMetadata1,viMetadata2), byrow=TRUE, nrow=3)
mtrxMatrixOf4 = matrix(c(viMetadata3, viMetadata4, viMetadata1,viMetadata2), byrow=TRUE, nrow=4)
test_that("Test that dummied data are correctly check when they are matrices (row major).",{
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1,3)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1,3)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1,3)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2,1)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2,2)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2,3)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2,6)$PASS, FALSE)
})

viMetadata1 = c(1,2,1,1,2,2,1,1,2,2)
viMetadata2 = c(1,2,3,3,1,1,2,3,3,2)
vdBug1 = c(0.1,0.3,0.2,0.6,0.8,0.9,0.3,0.5,0.6,0.4)
vdBug2 = c(0.1,0.0,0.2,0.0,0.8,0.0,0.3,0.0,0.6,0.0)
vdBug3 = c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)

test_that("Test that dummied data are correctly checked without dummying.",{
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 5, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 12, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 130, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug2, 2, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug2, 6, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug3, 1, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 9, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 2, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 12, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug2, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug2, 9, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug3, 1, FALSE)$PASS, FALSE)
})

viMetadata3 = c(1,1,2,2,3,3,4,4,5,5)
viMetadata4 = rep(1,length(viMetadata1))
mtrxMatrixOf2 = matrix(c(viMetadata1,viMetadata2), byrow=TRUE, nrow=2)
mtrxMatrixOf3 = matrix(c(viMetadata3, viMetadata1,viMetadata2), byrow=TRUE, nrow=3)
mtrxMatrixOf4 = matrix(c(viMetadata3, viMetadata4, viMetadata1,viMetadata2), byrow=TRUE, nrow=4)
test_that("Test that dummied data are correctly checked without dummying for matrices of data.",
{
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1, 6, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1, 12, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1, 9, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1, 8, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1, 2, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1, 8, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2, 12, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2, 0, FALSE)$PASS, TRUE)
})


### funcShuffleMatrix
#!#


### funcSpikeNewBug
context("Test funcSpikeNewBug")
vdContinuousMetadata1 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
vdContinuousMetadata2 = c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
vdContinuousMetadata3 = c(0.3,0.5,0.7,0.2,0.4,0.6,0.8,0.9,0.2,0.4)
viBinaryMetadata1 = c(1,1,1,1,1,2,2,2,2,2)
viBinaryMetadata2 = c(2,2,2,2,2,1,1,1,1,1)
viBinaryMetadata3 = c(2,2,1,1,2,2,1,1,2,1)
viFactorMetadata1 = c(1,1,1,2,2,2,3,3,3,4)
viFactorMetadata2 = c(1,2,3,4,1,2,3,4,1,4)
mtrxdMultContinuousMetadata1 = matrix(c(vdContinuousMetadata1, vdContinuousMetadata3),ncol=length(vdContinuousMetadata1), byrow=TRUE)
mtrxdMultContinuousMetadata2 = matrix(c(vdContinuousMetadata1, vdContinuousMetadata2),ncol=length(vdContinuousMetadata1), byrow=TRUE)
mtrxMultBinaryMetadata1 = matrix(c(viBinaryMetadata1, viBinaryMetadata3),ncol=length(viBinaryMetadata1), byrow=TRUE)
mtrxMultBinaryMetadata2 = matrix(c(viBinaryMetadata1, viBinaryMetadata2),ncol=length(viBinaryMetadata1), byrow=TRUE)
mtrxFactorMetadata1 = matrix(c(viFactorMetadata1, viFactorMetadata2),ncol=length(viFactorMetadata1), byrow=TRUE)
vdSparseBug1 = c(4,0,0,4,0,4,0,0,0,4)
vdSparseBug2 = c(3,0,0,4,0,8,0,0,0,4)
vdNotSparseBug1 = c(4,4,4,4,4,4,4,4,4,4)
vdNotSparseBug2 = c(3,6,5,4,7,8,9,2,3,4)
iMultiplier1 = 1
iMultiplier2 = 2
iMultiplier3 = 3

#test_that("funcSpikeNewBug: Test that relationships are made with simple inputs, not varying bug.",{
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata1, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(c(2.513699,2.843988,3.174277,3.504566,3.834855,4.165145,4.495434,4.825723,5.156012,5.486301),c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata2, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(c(5.486301,5.156012,4.825723,4.495434,4.165145,3.834855,3.504566,3.174277,2.843988,2.513699),c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viBinaryMetadata1, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viBinaryMetadata2, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viBinaryMetadata3, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viFactorMetadata1, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viFactorMetadata2, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#})

#test_that("funcSpikeNewBug: Test that relationships are made with simple inputs, varying bug.",{
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata1, vdNotSparseBug2, iMultiplier1),c_iRoundingPercision), round(c(1),c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata2, vdNotSparseBug2, iMultiplier1),c_iRoundingPercision), round(c(1),c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viBinaryMetadata1, vdNotSparseBug2, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viBinaryMetadata2, vdNotSparseBug2, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viBinaryMetadata3, vdNotSparseBug2, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viFactorMetadata1, vdNotSparseBug2, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viFactorMetadata2, vdNotSparseBug2, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#})

test_that("funcSpikeNewBug: Test that the multiplier increases the relationship.",{
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata1, vdNotSparseBug1, iMultiplier2),c_iRoundingPercision), round(c(1.027398,1.687976,2.348554,3.009133,3.669711,4.330289,4.990867,5.651446,6.312024,6.972602),c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata1, vdNotSparseBug1, iMultiplier3),c_iRoundingPercision), round(c(-0.4589032,0.5319641,1.5228315,2.5136989,3.5045663,4.4954337,5.4863011,6.4771685,7.4680359,8.4589032),c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata1, vdNotSparseBug2, iMultiplier2),c_iRoundingPercision), round(c(1),c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata1, vdNotSparseBug2, iMultiplier3),c_iRoundingPercision), round(c(1),c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata1, vdSparseBug1, iMultiplier2),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata1, vdSparseBug1, iMultiplier3),c_iRoundingPercision), round(1,c_iRoundingPercision))
})

#test_that("funcSpikeNewBug: Test that relationships are made with simple inputs and a sparse bug.",{
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata1, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(vdContinuousMetadata2, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viBinaryMetadata1, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viBinaryMetadata2, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viBinaryMetadata3, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viFactorMetadata1, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(viFactorMetadata2, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#})

#test_that("funcSpikeNewBug: Test that relationships are made with matrix inputs.",{
#  expect_equal(round(funcSpikeNewBug(mtrxdMultContinuousMetadata1, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(mtrxMultBinaryMetadata1, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(mtrxFactorMetadata1, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(mtrxdMultContinuousMetadata1, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(mtrxMultBinaryMetadata1, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(mtrxFactorMetadata1, vdSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#})

#test_that("funcSpikeNewBug: Test the effect of using metadata strongly indirectly associated.",{
#  expect_equal(round(funcSpikeNewBug(mtrxdMultContinuousMetadata2, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(mtrxMultBinaryMetadata2, vdNotSparseBug1, iMultiplier1),c_iRoundingPercision), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(mtrxdMultContinuousMetadata2, vdNotSparseBug1, iMultiplier3), round(1,c_iRoundingPercision))
#  expect_equal(round(funcSpikeNewBug(mtrxMultBinaryMetadata2, vdNotSparseBug1, iMultiplier3),c_iRoundingPercision), round(1,c_iRoundingPercision))
#})


### funcTruncatedRLNorm
iCount <- 1000000
dMu <- 2
dSD <- 3
dExp = funcGetExp( 2, 3 )
iThreshold <- NA
lResult <- funcTruncatedRLNorm( iCount, log( dMu ), log( dSD ), iThreshold )
test_that( "Test that funcTruncatedRLNorm creates an expected distribution without truncation." ,{
  expect_equal(round( dExp, 2), round( mean( lResult$Feature ), 2) )
  expect_equal( iCount, length( lResult$Feature ) )
  expect_equal( 0, length( lResult$LeftOver ) )
})

iCount <- 1000000
dMu <- 5
dSD <- 7
dExp = funcGetExp( 5, 7 )
iThreshold <- NA
lResult <- funcTruncatedRLNorm( iCount, log( dMu ), log( dSD ), iThreshold )
test_that( "Test that funcTruncatedRLNorm creates an expected distribution without truncation." ,{
  expect_equal(round( dExp ), round( mean( lResult$Feature ) ) )
  expect_equal( iCount, length( lResult$Feature ) )
  expect_equal( 0, length( lResult$LeftOver ) )
})

iCount <- 1000000
dMu <- 2
dSD <- 3
dExp = funcGetExp( 2, 3 )
iThreshold <- 3
lResult <- funcTruncatedRLNorm( iCount, log( dMu ), log( dSD ), iThreshold )
dSDThreeCount <- sum( lResult$LeftOver )
test_that( "Test that funcTruncatedRLNorm creates an expected distribution without truncation." ,{
  expect_equal(round( dExp, 2), round( mean( lResult$Feature + lResult$Feature ), 2) )
  expect_equal( iCount, length( lResult$Feature ) )
  expect_equal( 0, length( lResult$LeftOver ) )
})

iCount <- 1000000
dMu <- 2
dSD <- 7
dExp = funcGetExp( 2, 7 )
iThreshold <- 3
lResult <- funcTruncatedRLNorm( iCount, log( dMu ), log( dSD ), iThreshold )
test_that( "Test that funcTruncatedRLNorm creates an expected distribution without truncation." ,{
  expect_equal(round( dExp ), round( mean( lResult$Feature + lResult$LeftOver ) ) )
  expect_true( dSDThreeCount < sum( lResult$LeftOver ) )
  expect_equal( 0, length( which( lResult$Feature > iThreshold ) ) )
})

### funcZeroCorrectMatrix
mtrxShuffle = matrix( 1:10, ncol = 10, nrow = 10 )
mtrxResult = funcZeroCorrectMatrix( mtrxShuffle, funcGetRowMetric( mtrxShuffle, mean ), c() )
test_that( "funcZeroCorrectMatrix: Test that a matrx with no zeros is not shuffled." ,{
  expect_equal( mtrxResult, mtrxShuffle )
}

mtrxShuffle = matrix( 0, ncol = 10, nrow = 10 )
mtrxResult = funcZeroCorrectMatrix( mtrxShuffle, funcGetRowMetric( mtrxShuffle, mean ), c() )
test_that( "funcZeroCorrectMatrix: That a fully zero matrix is returned zero if no left over is available." ,{
  expect_equal( mtrxResult, mtrxShuffle )
}

viNumbers = 1:100
viNumbers[ c( 1, 12, 23, 34, 45, 56, 67, 78, 89, 100 ) ] = 0
mtrxShuffle = matrix( viNumbers, ncol = 10, nrow = 10 )
mtrxResult = funcZeroCorrectMatrix( mtrxShuffle, funcGetRowMetric( mtrxShuffle, mean ), c() )
test_that( "funcZeroCorrectMatrix: Zeros in samples that have other values are not shuffled." ,{
  expect_equal( mtrxResult, mtrxShuffle )
}

viNumbers = 1:100
viNumbers[ c( 1:10, 21:30, 41:50, 61:70, 81:90 ) ] = 0
mtrxShuffle = matrix( viNumbers, ncol = 10, nrow = 10 )
mtrxResult = funcZeroCorrectMatrix( mtrxShuffle, funcGetRowMetric( mtrxShuffle, mean ), c() )
test_that( "funcZeroCorrectMatrix: Zeros samples are shuffled out if there are samples with which to shuffle. Half of the samples are nonzero." ,{
  expect_equal( 0, length( which( colSums( mtrxResult$Data ) == 0 ) ) )
}

viNumbers = 1:100
viNumbers[ c( 1:90 ) ] = 0
mtrxShuffle = matrix( viNumbers, ncol = 10, nrow = 10 )
mtrxResult = funcZeroCorrectMatrix( mtrxShuffle, funcGetRowMetric( mtrxShuffle, mean ), c() )
test_that( "funcZeroCorrectMatrix: Zeros samples are shuffled out if there are samples with which to shuffle. One none zero sample." ,{
  expect_equal( 0, length( which( colSums( mtrxResult$Data ) == 0 ) ) )
}

## Makes negative numbers, not good
viNumbers = 1:100
viNumbers[ c( 1:95 ) ] = 0
mtrxShuffle = matrix( viNumbers, ncol = 10, nrow = 10 )
mtrxResult = funcZeroCorrectMatrix( mtrxShuffle, funcGetRowMetric( mtrxShuffle, mean ), c() )
test_that( "funcZeroCorrectMatrix: Zeros samples are shuffled out if there are samples with which to shuffle. One none zero sample, not enough signal to shuffle and no left over." ,{
  expect_equal( 0, length( which( colSums( mtrxResult$Data ) == 0 ) ) )
}

### Check Left Over use.
viNumbers = 1:100
viNumbers[ c( 1:95 ) ] = 0
mtrxShuffle = matrix( viNumbers, ncol = 10, nrow = 10 )
mtrxResult = funcZeroCorrectMatrix( mtrxShuffle, funcGetRowMetric( mtrxShuffle, mean ), rep( 100, 10 ) )
test_that( "funcZeroCorrectMatrix: Zeros samples are shuffled out if there are samples with which to shuffle. One none zero sample, not enough signal to shuffle and with left over." ,{
  expect_equal( 0, length( which( colSums( mtrxResult$Data ) == 0 ) ) )
}

### func_zero_inflate
context("Test func_zero_inflate checking that the made distributions are equal to the parameters.")
# No SD
dMeanOne = 5
dPercentZeroOne = .50
dSDOne = 1
iNumberSamples = 100000
vdTestOne = func_zero_inflate(dMeanOne,dPercentZeroOne,iNumberSamples,dSDOne)$Feature
viNotZeroOne = which(vdTestOne!=0)
dRealMeanOne = mean(log(vdTestOne[viNotZeroOne]))
dRealSDOne = sd(log(vdTestOne[viNotZeroOne]))
iRealLengthOne = length(vdTestOne)
dRealPercentZeroOne = length(which(vdTestOne==0))/iRealLengthOne
# SD
dMeanTwo = 10
dPercentZeroTwo = .3
dSDTwo = 2
iNumberSamples = 100000
vdTestTwo = func_zero_inflate(dMeanTwo,dPercentZeroTwo,iNumberSamples,dSDTwo)$Feature
viNotZeroTwo = which(vdTestTwo!=0)
dRealMeanTwo = mean(log(vdTestTwo[viNotZeroTwo]))
dRealSDTwo = sd(log(vdTestTwo[viNotZeroTwo]))
iRealLengthTwo = length(vdTestTwo)
dRealPercentZeroTwo = length(which(vdTestTwo==0))/iRealLengthTwo
# No negatives
dMeanThree = 1
dPercentZeroThree = .0
dSDThree = 5
iNumberSamples = 100000
vdTestThree = func_zero_inflate(dMeanThree,dPercentZeroThree,iNumberSamples,dSDThree)$Feature

test_that("Test that func_zero_inflate creates a distribution close to the given parameters.",{
  expect_equal(round(dRealMeanOne),round(dMeanOne))
  expect_equal(round(dRealSDOne),round(dSDOne))
  expect_equal(iRealLengthOne,iNumberSamples)
  expect_equal(round(dPercentZeroOne,1),round(dRealPercentZeroOne,1))

  expect_equal(round(dRealMeanTwo),round(dMeanTwo))
  expect_equal(round(dRealSDTwo),round(dSDTwo))
  expect_equal(iRealLengthTwo,iNumberSamples)
  expect_equal(round(dPercentZeroTwo,1),round(dRealPercentZeroTwo,1))

  expect_equal(length(which(vdTestThree<0)),0)
})