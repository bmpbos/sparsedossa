[TOC]

#**SparseDOSSA: Sparse Data Observations for the Simulation of Synthetic Abundances**#

AUTHORS: Boyu Ren (bor158@mail.harvard.edu)

##**Description**##

SparseDOSSA introduces a hierarchical model of microbial ecological population structure. It is capable of simulating realistic metagenomic data with known correlation structures, and thus provides a gold standard to enable benchmarking of statistical metagenomics methods. SparseDOSSA's model captures the marginal distribution of each microbial feature as a truncated, zero-inflated log-normal distribution, with parameters derived in turn from a parent log-normal distribution. The model can be effectively fit to reference microbial datasets in order to parameterize their microbes and communities, or to simulate synthetic datasets of similar population structure. Most importantly, it allows users to include both known feature-feature and feature-metadata correlation structures.

If you use this software, please cite :
_Will add once the paper is out_

-------------

##**Pre-requisites**##

SparseDOSSA requires *R 3.0.2* or higher with [*optparse*](https://cran.r-project.org/web/packages/optparse/optparse.pdf) library installed.

----------------------

##**Installation**##

SparseDOSSA can be obtained by either [Download SparseDOSSA](https://bitbucket.org/biobakery/sparsedossa/get/default.zip)

**OR**

Cloning the repository via the following command
``
$ hg clone https://bitbucket.org/biobakery/sparsedossa
``

--------------------------

##**Basic Usage**##

This section presents some basic usages of SparseDOSSA.

SparseDOSSA's hierarchical model is calibrated using the [PRISM dataset](https://www.dropbox.com/s/akgv0bv8bbpzcqo/prism.tsv?dl=0) by default. If you have your own reference dataset and would like to simulate data based on it, please follow the example below. Your dataset must be in a QIIME OTU table format, that is taxonomic units in rows and samples in columns, with each cell indicates the observed counts. Assume the file is *reference_OTU.txt*, using the following command, we can simulate microbiome dataset that has the same dimension and follows similar patterns with *reference_OTU.txt*:

```
#!cmd
$ synthetic_datasets_script.R -c reference_OTU.txt
```

Here is a basic example of simulating dataset with 50 features (OTUs), 10 samples and 10 metadata for each type (binary, quaternary and continuous), without any correlation structure. We use the default model parameters:

```
#!cmd
$ synthetic_datasets_script.R -f 50 -k 0 -n 10
```

If we want to add feature-metadata correlation, with 2% of the features spiked and each spiked feature correlated with two randomly selected metadata, we can use:

```
#!cmd
$ synthetic_datasets_script.R -f 50 -i 2 -k 0.02 -n 10 -p 10
```

You can also simulate dataset with feature-feature correlation only. Assume each spiked feature is correlated with two other randomly selected features and 10 of the features are spiked:

```
#!cmd
$ synthetic_datasets_script.R -f 50 -b 10 -m 2 -n 10 -p 10 --runBugBug
```

##**Output of SparseDOSSA**##

For each simulation using default model parameters, SparseDOSSA will produce three txt files: *SyntheticMicrobiome.pcl*, *SyntheticMicrobiome-Counts.pcl*, *SyntheticMicrobiomeParameterFile.txt*. The first two files contain the actual microbiome abundance data and the third file records values of model parameters, diagnostic information and spike-in assignment.

###**SyntheticMicrobiome.pcl**###

This file records the synthetic microbiome data for null community (no spike-in and outliers), outlier-added community without spike-in and final spiked data. We put samples in columns and features in rows. The first chunk of the file is metadata, with row names **Metdata_***. The second chunk is for null community, with row names **Feature_Lognormal_***. The third chunk is for outlier-introduced community, with row names **Feature_Outlier_***. The last chunk is for spiked data, with row names **Feature_spike**. This file records relative abundance data.

###**SyntheticMicrobiome-Counts.pcl**###

This file has the same organization as *SyntheticMicrobiome.pcl* but records raw counts data.

###**SyntheticMicrobiomeParameterFile.txt**###

This file records diagnostic information and values of model paramters as well as the spike-in assignment. The most part of this file is used only for debugging. Users can focus on lines after **Minimum Spiked-in Samples**. Those lines record which metadata are correlated with which feature. The format is all metadata that are correlated with a specific features are listed under the name of the feature.


##**Full command-line options**##

The full command-line options of SparseDOSSA is listed below:

```
#!

$ Rscript R/synthetic_datasets_script.R --help
Usage: synthetic_datasets_script.R [options] NormalizedFile(Optional) CountFile(Optional) TrueFile(Optional)


Options:
        -a VARIANCE_SCALE, --variance_scale=VARIANCE_SCALE
                Tuning parameters for noise in bug-bug associations. Non-negative values are expected. Multiple values should be comma-separated. 
                Values will be recycled if the length does not match the number of associations

        -b BUGS_TO_SPIKE, --bugs_to_spike=BUGS_TO_SPIKE
                Number of bugs to correlate with others.  A non-negative integer value is expected.

        -c CALIBRATE, --calibrate=CALIBRATE
                Calibration file for generating the random log normal data. TSV file (column = feature)

        -d DATASETCOUNT, --datasetCount=DATASETCOUNT
                The number of bug-bug spiked datasets to generate.  A positive integer value is expected.

        -e READ_DEPTH, --read_depth=READ_DEPTH
                Simulated read depth for counts. A positive integer value is expected.

        -f NUMBER_FEATURES, --number_features=NUMBER_FEATURES
                The number of features per sample to create. A positive integer value is expected.

        -g BUGBUGCOEF, --bugBugCoef=BUGBUGCOEF
                A vector of string separated values for the association coefficients for the bug-bug associations. 
                At least two values, an intercept and slope, must be given. Values are comma-separated. Example: 0,0.5.

        -i SPIKECOUNT, --spikeCount=SPIKECOUNT
                Counts of spiked metadata used in the spike-in dataset. These values should be comma delimited values, in the order of the spikeStrength values (if given). 
                Can be one value, in this case the value will be repeated to pair with the spikeCount values (if multiple are present). Example 1,2,3

        -j LEFSE_FILE, --lefse_file=LEFSE_FILE
                Folder containing lefSe inputs.

        -k PERCENT_SPIKED, --percent_spiked=PERCENT_SPIKED
                The percent of features spiked-in. A real number between 0 and 1 is expected.

        -l MINLEVELPERCENT, --minLevelPercent=MINLEVELPERCENT
                Minimum percent of measurements out of the total a level can have in a discontinuous metadata (Rounded up to the nearest count). 
                A real number between 0 and 1 is expected.

        -m MAX_DOMAIN_BUGS, --max_domain_bugs=MAX_DOMAIN_BUGS
                Maximum number of bugs with which each correlated bug can be associated with. A positive integer greater than 0 is expected.

        -n NUMBER_SAMPLES, --number_samples=NUMBER_SAMPLES
                The number of samples to generate. A positive integer greater than 0 is expected.

        -o MAX_PERCENT_OUTLIERS, --max_percent_outliers=MAX_PERCENT_OUTLIERS
                The maximum percent of outliers to spike into a sample. A real number between 0 and 1 is expected.

        -p NUMBER_METADATA, --number_metadata=NUMBER_METADATA
                Indicates how many metadata are created. number_metadata*2 = number continuous metadata, number_metadata = number binary metadata, 
                number_metadata = number quaternary metadata. A positive integer greater than 0 is expected.

        -r SPIKESTRENGTH, --spikeStrength=SPIKESTRENGTH
                Strength of the metadata association with the spiked-in feature. These values should be comma delimited and in the order of the spikeCount values (if given). 
                Can be one value, in this case the value wil be repeated to pair with the spikeStrength values (if multiple are present). Example 0.2,0.3,0.4.

        -s SEED, --seed=SEED
                A seed to freeze the random generation of counts/relative abundance. If left as default (NA), generation is random. If seeded, 
                data generation will be random within a run but identical if ran again under the same settings. An integer is expected.

        -t PERCENT_OUTLIER_SPIKINS, --percent_outlier_spikins=PERCENT_OUTLIER_SPIKINS
                The percent of samples to spike in outliers. A real number between 0 to 1 is expected.

        -u MINOCCURENCE, --minOccurence=MINOCCURENCE
                Minimum counts a bug can have for the ocurrence quality control filter used when creating bugs. 
                (Filtering minimum number of counts in a minimum number of samples). A positive integer is expected.

        -v, --verbose
                If True logging and plotting is made by the underlying methodology. This is a flag, 
                it is either included or not included in the commandline, no value needed.

        -w MINSAMPLE, --minSample=MINSAMPLE
                Minimum samples a bug can be in for the ocurrence quality control filter used when creating bugs.
                (Filtering minimum number of counts in a minimum number of samples). A positive integer is expected.

        -x SCALEPERCENTZEROS, --scalePercentZeros=SCALEPERCENTZEROS
                A scale used to multiply the percent zeros of all features across the sample after it is derived from the relatiohships with it and 
                the feature abundance or calibration file. Requires a number greater than 0. A number greater than 1 increases sparsity, 
                a number less than 1 decreases sparsity. O removes sparsity, 1 (default) does not change the value and the value.

        -y ASSOCIATION_TYPE, --association_type=ASSOCIATION_TYPE
                The type of association to generate. Options are 'linear' or 'rounded_linear'.

        -z, --noZeroInflate
                If given, zero inflation is not used when generating a feature. This is a flag, it is either included or not included in the commandline, no value needed.

        --noRunMetadata
                If given, no metadata files are generated. This is a flag, it is either included or not included in the commandline, no value needed.

        --runBugBug
                If given, bug-bug interaction files are generated in addition to any metadata files. This is a flag, it is either included or 
                not included in the commandline, no value needed.

        -h, --help
                Show this help message and exit

```