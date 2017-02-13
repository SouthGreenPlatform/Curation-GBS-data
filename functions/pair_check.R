# working on scaffolds with 2 loci
#' Title
#'
#' @param scaff character: name of a scaffold with 2 loci 
#' @param genot_data  data.frame containing genotypic data
#' @param id data.frame containing locus identifiers
#' @param type character: type of cross
#' @param code character: code used to score data
#' @param mc character: name (in 'id') of the column containing the marker names
#' @param sc character: name (in 'id') of the column containing the scaffold names
#'
#' @return data_frame with proximity information (two lines) 
#' @export
#'
#' @examples
pair_check<-function(scaff,genot_data,id,type,code="abh-",mc=marker_column,sc=scaffold_column){
  setloc<-id[mc][id[sc]==scaff]
  g<-genot_data[setloc,]
  d<-distance(g[1,],g[2,],type,code)
  out<-data.frame(scaffold=scaff,locus = setloc[1],posit=posit_loc(setloc[1]),target= setloc[2],dist=d[[1]],f1=d[[2]],f2=d[[3]],signif=d[[4]])
  rbind(out,data.frame(scaffold=scaff,locus = setloc[2],posit=posit_loc(setloc[2]),target= setloc[1],dist=d[[1]],f1=d[[2]],f2=d[[3]],signif=d[[4]]))
}
#example
#pair_check("C69556603",genot,id,"BC1",code="DTHM")