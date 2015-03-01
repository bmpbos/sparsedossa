\name{sparseDOSSA}
\alias{sparseDOSSA}
\title{
Sparse Data Observations for Simulating Synthetic Abundance}
\description{
Sparse Data Observations for Simulating Synthetic Abundance}
\usage{
sparseDOSSA(variance_scale = 1,
	 bugs_to_spike = 0,
	 calibrate = NA, 
	 datasetCount = 1, 
	 read_depth = 8030,
	 number_features = 300,
	 bugBugCoef = "0,0.5",
	 spikeCount = "1", 
	 percent_spiked = 0.03,
	 minLevelPercent = 0.1, 
	 max_domain_bugs = 2, 
	 number_samples = 50,
	 max_percent_outliers = 0.05, 
	 number_metadata = 5, 
	 spikeStrength = "1.0",
	 seed = NA, 
	 percent_outlier_spikins = 0.05,
	 minOccurence = 0,
	 verbose = TRUE, minSample = 0,
	 scalePercentZeros = 1, 
	 association_type = "linear",
	 noZeroInflate = FALSE, 
	 noRunMetadata = FALSE, 
	 runBugBug = FALSE, help = FALSE)
}
\arguments{
  \item{variance_scale}{
  Tuning parameters for noise in bug-bug associations, Non-negative values are expected,Multiple values should be comma-separated,Values will be recycled if the length doesn't match the number of associations }
  \item{bugs_to_spike}{
   Number of bugs to correlate with others.  A non-negative integer value is expected}
  \item{calibrate}{
  Calibration file for generating the random log normal data. TSV file (column = feature)
}
  \item{datasetCount}{
 The number of bug-bug spiked datasets to generate.  A positive integer value is expected
}
  \item{read_depth}{
  Simulated read depth for counts. A positive integer value is expected.
}
  \item{number_features}{
   The number of features per sample to create. A positive integer value is expected.
}
  \item{bugBugCoef}{
  A vector of string separated values for the association coefficients for the bug-bug associations - At least two values, an intercept and slope, must be given, Values are comma-separated -Example: 0,0.5.
}
  \item{spikeCount}{
  Counts of spiked metadata used in the spike-in dataset - These values should be comma delimited values, in the order of the spikeStrength values (if given),  Can be one value, in this case the value will be repeated to pair with the spikeCount values (if multiple are present),Example 1,2,3"}
  \item{percent_spiked}{
  The percent of features spiked-in. A real number between 0 and 1 is expected
}
  \item{minLevelPercent}{
  Minimum percent of measurements out of the total a level can have in a discontinuous metadata (Rounded up to the nearest count),A real number between 0 and 1 is expected}
  \item{max_domain_bugs}{
  Maximum number of bugs with which each correlated bug can be associated with, A positive integer greater than 0 is expected
}
  \item{number_samples}{
  The number of samples to generate. A positive integer greater than 0 is expected
 }
  \item{max_percent_outliers}{
  The maximum percent of outliers to spike into a sample. A real number between 0 and 1 is expected
 }
  \item{number_metadata}{
    Indicates how many metadata are created, number_metadata*2 = number continuous metadata, number_metadata = number binary metadata, number_metadata = number quaternary metadata, A positive integer greater than 0 is expected.
 }
  \item{spikeStrength}{
  Strength of the metadata association with the spiked-in feature, These values should be comma delimited and in the order of the spikeCount values (if given),Can be one value, in this case the value wil be repeated to pair with the spikeStrength values (if multiple are present)-Example 0.2,0.3,0.4.",}
  \item{seed}{
    "A seed to freeze the random generation of counts/relative abundance,If left as default (NA), generation is random - If seeded, data generation will be random within a run but identical if ran again under the same settings,an integer is expected. 
 }
  \item{percent_outlier_spikins}{
  The percent of samples to spike in outliers. A real number between 0 to 1 is expected.
}
  \item{minOccurence}{
   Minimum counts a bug can have for the ocurrence quality control filter used when creating bugs ,( Filtering minimum number of counts in a minimum number of samples)-  A positive integer is expected
 }
  \item{verbose}{
   If True logging and plotting is made by the underlying methodology,   This is a flag, it is either included or not included in the commandline, no value needed
}
  \item{minSample}{
  Minimum samples a bug can be in for the ocurrence quality control filter used when creating bugs,  ( Filtering minimum number of counts in a minimum number of samples) - A positive integer is expected. 
}
  \item{scalePercentZeros}{ 
  A scale used to multiply the percent zeros of all features across the sample after it is derived from the relatiohships with it and the feature abundance or calibration file, - Requires a number greater than 0 - A number greater than 1 increases sparsity, a number less than 1 decreases sparsity ,O removes sparsity, 1 (default) does not change the value and the value.",
 }
  \item{association_type}{
  The type of association to generate - Options are 'linear' or 'rounded_linear'
}
  \item{noZeroInflate}{
  If given, zero inflation is not used when generating a feature- This is a flag, it is either included or not included in the commandline, no value needed.
 }
  \item{noRunMetadata}{
   If given, no metadata files are generated, This is a flag, it is either included or not included in the commandline, no value needed.
 }
  \item{runBugBug}{
  If given, bug-bug interaction files are generated in addition to any metadata files- This is a flag, it is either included or not included in the commandline, no value needed}
  \item{help}{
  Display help
}
}
 
\value{
Need to add this
}
\author{
Timothy Tickle<ttickle@hsph.harvard.edu>,Curtis Huttenhower <chuttenh@hsph.harvard.edu>
}
\examples{
 
sparseDOSSA (
	variance_scale = 1,
	bugs_to_spike = 0,
	calibrate = NA,
	datasetCount = 1,
	read_depth = 8030,
	number_features = 300,
	bugBugCoef =  "0,0.5",
	spikeCount = "1",
	percent_spiked = 0.03,
	minLevelPercent =  0.1,
	max_domain_bugs = 2,
	number_samples = 50, 
	max_percent_outliers = 0.05,
	number_metadata = 5,
	spikeStrength =  "1.0",
	seed =  1,
	percent_outlier_spikins = 0.05,
	minOccurence =  0,
	verbose =  TRUE,
	minSample =  0,
	scalePercentZeros = 1,
	association_type =  "linear",
	noZeroInflate =  FALSE,
	noRunMetadata = FALSE,
	runBugBug =  FALSE,
	help =  FALSE
)
}
 
