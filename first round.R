# To run the workflow, user needs to enter the value of a number of variables, especially the location of the files and their names. 
# There is a parameter set for each part of the work flow
# /!\ Never change the *names* of the variables /!\
#######################################################################
# parameter for graphic outputs
# Set to 'FALSE' in case an error occurs when attempting to draw a plot 
# problem encountered when combining apple computers and unix workstation
draw_plot=TRUE
# Enter here free comments relative to the present dataset
comment0 <-paste("Data curation for the set of 240+2 individuals from plates 1 to 6. ",
				 "data are in", 
				 "\\scriptsize D://Mes donnees/sequencage/2017/article sequence/Tassel_to_J/data/",
                 "\\normalsize" ,sep="\n")				 
#suggestions for comment variables. 
# the following sample text shows how add some latex formatting code to a multiline comment texts. 
#•This is normal latex code except that the backslash needs to be doubled,
# otherwise it will be considered as part of a special character. 
# comment0<-paste("this is a sample text. \\textit{this is in italics},  this is not.   ",
#                 "\\textbf{some bold text}   ",
#                 "\\footnotesize a bit smaller    ",
#                 "\\scriptsize  even smaller    ",
#                 "\\tiny and this is very very small    ",
#                 "    ",
#                 "\\normalsize back to normal font",sep="\n")
# This allows the above to be printed as formatted text in the output document.
# 
##########################################################
# Parameters for step 1
##########################################################
# Tassel was executed and data were filtered. The FILTER.tab contains genotyping data with coverage information. 
# This step can take one or more input files (eg. because data were treated in parallel). 
# For each input file, a genotypic file and a marker info file are created.
##########################################################

# Enter here your free comments for step 1
comment1 <- ""
# place of initial data files
where_in <- "data/set1"
# name of input file.
# Note: to save memory, several files can be passed with syntax c(file1,file2, file3)
file.in <- "first_round.tab"
#place of files produced at step 1
where_out_1 <- "data/set1/step1"
##name of output genotype file
file.value <- "genot_first_round.txt"
# name of marker info file
# Note: to save memory, several files can be passed with syntax c(file1,file2, file3)
file.ident <- "markid_first_round.txt"
#p.value cut-off Suggested something like 1/(2* number of rows). Negative values cancel the filtef
P.value.min <- -1
# Cross type: (only BC1 is implemented currently)
type <- "BC1"
# Name the parents. Two values are required (In BC1, the heterzygous parent comes first).
# if one parent was not geotyped, set the correspondig flag to FALSE
parent <- c("P3637","P49")
# flags to activate/deactivate tests one one parent
# at least one of the following should be TRUE
has.hetero <- TRUE # heterozygous parent
has.homo <- TRUE #  homozygous parent
#list of  excluded genotypes
others <- NULL
#others<-c("DA10D","DegNour")
# names of the columns in the filter.tab file.
marker_column <- "marker"
scaffold_column <- "scaffold"
position_column <- "position"
# locus_column <- "marker" 
# scaffold_column<- "chomosome"
# position_column <- "position"

##########################################################
# Parameters for step 2
##########################################################
# At this stage, the genotypes are assigned to the parents. 
# Only loci whose segregation is compatible with the mating plan are retained. 
# Statistics are computed on the disjunction and the results are written into the locus info file.
# The columns of the genotype table are reordered. First the parents, then the PROGENIES in order of increasing number of missing data. 
# All genotyping files and all locus information file are merged at this step. 
##########################################################

# uses where_out_1, file.value, file.ident
# Enter here your free comments
comment2 <- ""
# name of input file
# place of files produced at step 2
where_out_2 <-"D:/Mes donnees/sequencage/2017/article sequence/Tassel_to_J/TtJGIT/data/set1/step2"
  
file.value.2 <- "genot_1to6_small.txt"
file.ident.2 <- "markid_1to6_small.txt"


# Thresholds
# discard individuals having less than 100*cut.ind % good data. 
# This parameter should be liberal since more tests will come later. 
cut.ind <- 0.5
# The following thresholds flag markers for elimination at a later stage
# minimum frequency of the minority allele
balanced <- .15
# minimum proportion of valid data (excluding missing and unexpected alleles)
few_good_data <- 0.7
# maximum proportion of unexpected alleles
many_wrong <- 0.1

##########################################################
# Parameters for step 3
##########################################################
# Codes the genotype file in a form usable by mapping software. Parents are removed from the output file.
##########################################################

# uses where_out_2, file.value.2 and file.ident.2 for data entry
# Enter here your free comments
comment3 <- ""
#place of final files
where_out_3 <- "D:/Mes donnees/sequencage/2017/article sequence/Tassel_to_J/TtJGIT/data/set1/step3"
# name of coded genotype file
file.value.3 <- "genot_1to6_small_3.txt"
#name of updated marker info file 
file.ident.3 <- "markid_1to6_small_3.txt"
# Use a four letter code to describe the genotypes  (preferred for Joinmap: "ABHU")
# letters are in upper case (at a letter stage, lower cases will be used to identify corrections)
# first letter: a homozygous parental genotype
# Second letter: the other homozygous parental genotype (not used in BC1)
# third letter: the heterozygous genotype
# fourth letter: unknown genotype
# Joinmap might require converting to lower case at the end of the process (to be checked)
genotype_code <- "ABHU"
# do not change the next two lines
code_U <- strsplit(genotype_code,split="")[[1]]
code_L <- tolower(code_U)
##########################################################
# Parameters for step 4
##########################################################
# the purpose of this treatment chain is to automatize the selection of useful markers. 
# There are two main issues: redundant markers and outliers ie markers that segregate
# independently from the other loci of the same scaffold (possibly chimeric scaffolds or wrong assignment of marker to scaffold).
##########################################################
# uses where_out_3, file.value.3 and file.ident.3 for data entry

# Enter here your free comments
comment4 <- ""
# place of final files
where_out_4 <- "D:/Mes donnees/sequencage/2017/article sequence/Tassel_to_J/TtJGIT/data/set1/step4"
# names of files
pair_outlier <- "pair_outlier.txt"
pair_redundant <- "pair_redundant.txt"
multiple_outlier <- "multiple_outlier.txt"
multiple_redundant <- "multiple_redundant.txt"
multiple_all <- "multiple_all.txt"
reduced_data <- "reduced_data.txt"
reduced_info <- "reduced_info.txt"
indiv_select <- "indiv_select.txt"
# two distances between loci of the same scaffold are computed
# the first one is the expectation of the *number* of recombinations
# under the hypothesis of no linkage (taking into account missing data)
# It is used to identify outliers
# suggested value:  a little less than N/2-sqrt(N) where N is the number of descendents
# This value may need to be adapted based on the histograms produced at this step.
t_outlier <- 60
# the second distance is the observed recombination rate, based only on individuals which have both loci 
# issues a warning for large gaps between loci of the same scaffold.
l_gap<-0.05
# in the second step of imputation, chosen individuals are assigned a value based on a majority rule. 
# Accept imputation only if the frequency of the most frequent allele is above 'majority'
majority <- 0.8
# Directories containing the functions
where_functions <- "functions"
#Parameters to select individuals and loci
stop_remove <- 20 # do not consider removing more than 'stop_remove' individuals
max_miss_loc <- 0.2 # retain only loci with 'min_miss_loc*100' % missing data. 
##########################################################
# Parameters for step 5
##########################################################
# Last step: Examine last figure of previous step and chose the number of individuals to eliminate.
# This step takes one data and one information file.  
# Another file set with less individual and less loci (but also less missing data) is produced.
##########################################################


# uses where_out_4,for data entry

# Enter here your free comments
comment5 <- ""
# place of final files
where_out_5 <- "D:/Mes donnees/sequencage/2017/article sequence/Tassel_to_J/TtJGIT/data/set1/step5"
#Number of individuals used for plotting

removed <- 18
