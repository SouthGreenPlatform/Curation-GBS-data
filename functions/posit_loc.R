#extracts the position of a locus on the scaffold.
#Here the label is the position is simply the number following "M"
#Obviously depends on the convention adopted to name loci. 
#The 3 dots '...' make it possible to have more parameters in modified function.

posit_loc<-function(locus,...)
  as.numeric(substring(locus,regexpr("M",locus)[1]+1))