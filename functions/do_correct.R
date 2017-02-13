#' Title correct errors on one scaffold
#'
#' @param scaffold character: name of the scaffold
#' @param info data.frame: information on the loci
#' @param genot data.frame: genoryping data
#' @param codeD character: genotype encoding in the data file
#' @param codeR character: genotype encoding in the resulting file
#'
#' @return corrected genotype in 
#' @export
#'
#' @examples
do_correct<-function(scaffold,info,genot,codeD="ABHU",codeR="ABHU"){
  #les locus
  select_loc<-as.character(info$scaffold) ==scaffold
  loc_scaff<-as.character(info$locus)[select_loc]
  #les parents
  p0<-as.character(info$parent)[select_loc] 
  #les statuts
  s0<-info$select[select_loc]
  #oter les outliers
  s1<-s0[s0 != "outlier"]
  loc_scaff<-loc_scaff[s0 != "outlier"]
  #nombre d'all?les retenus
  n_scaff<-length(loc_scaff)
  if(n_scaff==0) return(NULL)
  # qui est retenu?
  s2<-ifelse(s1 %in% c("large","small","chosen","single") , 1 ,0)
  #les g?notypes
  
  genot_scaff<-genot[loc_scaff,]
  
  # make  character strings
  u1<-apply(genot_scaff, 2, paste, collapse="")
  if(codeD == codeR) v1<-u1 else v1<-transcode(toupper(u1),codeD, codeR)
  w1<-correctGenotype(toupper(v1),scaffoldRules)
  
  if(codeD == codeR) x1<-w1 else x1<-transcode(toupper(w1),codeR, codeD)
  out<-sapply(x1,function(x) strsplit(x,"")[[1]])
  rownames(out)<-loc_scaff
  as.data.frame( out)
}
