[TOC]

#**SparseDOSSA: Sparse Data Observations for the Simulation of Synthetic Abundances**#

AUTHORS: Boyu Ren (bor158@mail.harvard.edu)

##**Description**##

SparseDOSSA introduces a hierarchical model of microbial ecological population structure. It is capable of simulating realistic metagenomic data with known correlation structures, and thus provides a gold standard to enable benchmarking of statistical metagenomics methods. SparseDOSSA's model captures the marginal distribution of each microbial feature as a truncated, zero-inflated log-normal distribution, with parameters derived in turn from a parent log-normal distribution. The model can be effectively fit to reference microbial datasets in order to parameterize their microbes and communities, or to simulate synthetic datasets of similar population structure. Most importantly, it allows users to include both known feature-feature and feature-metadata correlation structures.

-------------

##**Pre-requisites**##

SparseDOSSA requires R 3.0.2 or higher with argparse library installed.

----------------------

##**Installation**##

SparseDOSSA can be obtained by cloning the repository via the following commands
``$ hg clone https://bitbucket.org/biobakery/sparsedossa``

--------------------------

##**Basic Usage**##

This section presents some basic usages of SparseDOSSA.

SparseDOSSA's hierarchical model is calibrated using the PRISM dataset by default. If you have your own reference dataset and would like to simulate data based on it, here is the basic example of fitting a given dataset to SparseDOSSA's hierarchical model. Assume the dataset is in a QIIME OTU table format with name reference_OTU.txt. An output file named SyntheticMicrobiomeParameterFile.txt will record all the point estimates of the parameters:

```
#!cmd
$ synthetic_datasets_script.R -c reference_OTU.txt
```

Here is an example of simulating dataset with 50 features (microbes), 10 samples and 10 metadata for each category (binary, quarenary and continuous), without any correlation structure. We use the default estimates of model parameters:

```
#!cmd
$ synthetic_datasets_script.R -f 50 -k 0 -n 10
```

If we want to only add feature-metadata correlation, with 2\% of the features are spiked and each spiked feature is correlated with two metadata, we can use:

```
#!cmd
$ synthetic_datasets_script.R -f 50 -i 2 -k 0.02 -n 10 -p 10
```

You can also simulate dataset with feature-feature correlation only. Assume each spiked feature is correlated with two other features and 10 of the features are spiked:

```
#!cmd
$ synthetic_datasets_script.R -f 50 -b 10 -m 2 -n 10 -p 10 --runBugBug
```

For advanced options and other analysis types (such as strain tracking) please refer to the full command-line options.

##**Full command-line options**##
