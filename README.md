# Selection-Framework
Framework to analyse past selection using time-series data (more details in Loog et al. (2017) https://academic.oup.com/mbe/article-lookup/doi/10.1093/molbev/msx142) 

INPUT

The code takes input in the form of a comma separated variants (CSV) file. The two input files used in the analyses in Loog et al. 2017 are provided as examples. 

The example input files contain:

1) Information about the genotypes for each sample or pooled allele frequencies at a particular time point for a given locus. The genotype information used by the code is the form of counts (i.e. how many times have the ancestral and derived alleles been observed).

2) Sample dates (can be either a point estimates or a ranges). The dates in the example input file are in the form of years AD. Note that the code needs to be modified in order to accommodate any other date/age format (e.g. years before present). 

3) The input file also contains information about the geographical origins of the archaeological samples although this in not directly used in the analysis but can be used to subset the data should there be need. 

ANALYSIS

The code asks to specify fixed parameters and ranges for the three free parameters. Both currently set as specified in Loog et al. (2017).

The fixed parameters are:

1)	Date of start of the gene flow period (in years AD) (gf.period[1])
2)	Date of end of the gene flow period (in years AD) (gf.period[2])
3)	The amount/ fraction of gene flow from an external population to the population of interest per time unit of a year (total.gf)
4)	The frequency of the derived allele in the external population (f.external)

The three free parameters are: 

1)	Selection coefficient (s.values)
2)	Starting time of Selection (sel.start)
3)	Ancestral allele frequency (f.anc)

OUTPUT

As output the code produces two arrays. 

1)	A 3-dimentional array (results.cube) that contains the likelihoods of the data for each parameter combination considered in the parameters sweep. Each dimension represents a free parameter (dimension 1 = selection, dimension 2 = starting time of selection, dimension 3 = ancestral allele frequency). 
This array is used to calculate the marginal likelihoods for each free parameter (i.e. ML.s for the Selection Coefficients; ML.sel.start for the Starting Times of Selection; ML.f.anc for the Ancestral Allele Frequencies) or joint marginal likelihoods for a combination of free parameters.
Code is provided to calculate & plot marginal likelihoods for each free parameter.

2)	A matrix (ML.curve) that contains a likelihood value for each specified allele frequency at given, equally spaced time points. Code is provided to visualize the posterior distribution of allele frequency trajectories through time, which allows identifying not only the most likely allele frequency trajectory through time but also to quantify the uncertainty in these estimates.

Reference:
Loog et al. (2017) Inferring allele frequency trajectories from ancient DNA indicates that selection on a chicken gene coincided with changes in medieval husbandry practices. Molecular Biology and Evolution. doi:10.1093/molbev/msx142.
