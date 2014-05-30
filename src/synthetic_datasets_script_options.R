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

library(optparse)

option_list = list(
      make_option(
          c("-a","--variance_scale"),
          type="character",
          default = "1",
          help=paste(
              "Tuning parameters for noise in bug-bug associations",
              "Non-negative values are expected",
              "Multiple values should be comma-separated",
              "Values will be recycled if the length doesn't match the number of associations",
              sep=". "
              )
          ),
      make_option(
          c("-b","--bugs_to_spike"),
          type="integer",
          default=0,
          help="Number of bugs to correlate with others.  A non-negative integer value is expected."
          ),
      make_option(
          c("-c","--calibrate"),
          type="character",
          default=NA,
          help="Calibration file for generating the random log normal data. TSV file (column = feature)"
          ),
      make_option(
          c("-d", "--datasetCount"),
          type="integer",
          default = 1,
          help="The number of bug-bug spiked datasets to generate.  A positive integer value is expected."
          ),
      make_option(
          c("-e","--read_depth"),
          type="integer",
          default=8030,
          help="Simulated read depth for counts. A positive integer value is expected."
          ),
      make_option(
          c("-f","--number_features"),
          type="integer",
          default=300,
          help="The number of features per sample to create. A positive integer value is expected."
          ),
      make_option(
          c("-g","--bugBugCoef"),
          type="character",
          default = "0,0.5",
          help = paste(
              "A vector of string separated values for the association coefficients for the bug-bug associations",
              "At least two values, an intercept and slope, must be given",
              "Values are comma-separated",
              "Example: 0,0.5.",
              sep=". "
              )
          ),
      make_option(
          c("-i","--spikeCount"),
          type = "character",
          default = "1",
          help = paste(
              "Counts of spiked metadata used in the spike-in dataset",
              "These values should be comma delimited values, in the order of the spikeStrength values (if given)",
              "Can be one value, in this case the value will be repeated to pair with the spikeCount values (if multiple are present)",
              "Example 1,2,3",
              sep = ". "
              )
          ),
      make_option(
          c("-j","--lefse_file"),
          type="character",
          default=NULL,
          help="Folder containing lefSe inputs."
          ),
      make_option(
          c("-k","--percent_spiked"),
          type="double",
          default=.03,
          help="The percent of features spiked-in. A real number between 0 and 1 is expected."
          ),
      make_option(
          c("-l","--minLevelPercent"),
          type="double",
          default=.1,
          help=paste(
              "Minimum percent of measurements out of the total a level can have in a discontinuous metadata (Rounded up to the nearest count)",
              "A real number between 0 and 1 is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-m","--max_domain_bugs"),
          type="integer",
          default=2,
          help=paste(
              "Maximum number of bugs with which each correlated bug can be associated with",
              "A positive integer greater than 0 is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-n","--number_samples"),
          type="integer",
          default=50,
          help="The number of samples to generate. A positive integer greater than 0 is expected."
          ),
      make_option(
          c("-o","--max_percent_outliers"),
          type="double",
          default=.05,
          help="The maximum percent of outliers to spike into a sample. A real number between 0 and 1 is expected."
          ),
      make_option(
          c("-p","--number_metadata"),
          type="integer",
          default=5,
          dest='number_metadata',
          help=paste(
              "Indicates how many metadata are created",
              "number_metadata*2 = number continuous metadata, number_metadata = number binary metadata, number_metadata = number quaternary metadata",
              "A positive integer greater than 0 is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-r","--spikeStrength"),
          type = "character",
          default = "1.0",
          help = paste(
              "Strength of the metadata association with the spiked-in feature",
              "These values should be comma delimited and in the order of the spikeCount values (if given)",
              "Can be one value, in this case the value wil be repeated to pair with the spikeStrength values (if multiple are present)",
              "Example 0.2,0.3,0.4.",
              sep = ". "
              )
          ),
      make_option(
          c("-s","--seed"),
          type="integer",
          default=NA,
          help=paste(
              "A seed to freeze the random generation of counts/relative abundance",
              "If left as default (NA), generation is random",
              "If seeded, data generation will be random within a run but identical if ran again under the same settings",
              "An integer is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-t","--percent_outlier_spikins"),
          type="double",
          default=.05,
          help="The percent of samples to spike in outliers. A real number between 0 to 1 is expected."
          ),
      make_option(
          c("-u","--minOccurence"),
          type="integer",
          default=0,
          help=paste(
              "Minimum counts a bug can have for the ocurrence quality control filter used when creating bugs",
              "( Filtering minimum number of counts in a minimum number of samples)",
              "A positive integer is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-v","--verbose"),
          action="store_false",
          default = TRUE,
          help=paste(
              "If True logging and plotting is made by the underlying methodology",
              "This is a flag, it is either included or not included in the commandline, no value needed.",
              sep = ". "
              )
          ),
      make_option(
          c("-w","--minSample"),
          type="integer",
          default=0,
          help=paste(
              "Minimum samples a bug can be in for the ocurrence quality control filter used when creating bugs",
              "( Filtering minimum number of counts in a minimum number of samples)",
              "A positive integer is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-x","--scalePercentZeros"),
          type="double",
          default=1,
          help=paste(
              "A scale used to multiply the percent zeros of all features across the sample after it is derived from the relatiohships with it and the feature abundance or calibration file",
              "Requires a number greater than 0",
              "A number greater than 1 increases sparsity, a number less than 1 decreases sparsity",
              "O removes sparsity, 1 (default) does not change the value and the value.",
              sep = ". "
              )
          ),
      make_option(
          c("-y","--association_type"),
          type="character",
          default = "linear",
          help=paste(
              "The type of association to generate",
              "Options are 'linear' or 'rounded_linear'.",
              sep = ". "
              )
          ),
      make_option(
          c("-z","--noZeroInflate"),
          action="store_true",
          default = FALSE,
          help=paste(
              "If given, zero inflation is not used when generating a feature",
              "This is a flag, it is either included or not included in the commandline, no value needed.",
              sep = ". "
              )
          ),
      make_option(
          c("--noRunMetadata"),
          action="store_true",
          default = FALSE,
          help=paste(
              "If given, no metadata files are generated",
              "This is a flag, it is either included or not included in the commandline, no value needed.",
              sep = ". "
              )
          ),
      make_option(
          c("--runBugBug"),
          action="store_true",
          default = FALSE,
          help=paste(
              "If given, bug-bug interaction files are generated in addition to any metadata files",
              "This is a flag, it is either included or not included in the commandline, no value needed.",
              sep = ". "
              )
          )
    )
