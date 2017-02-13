#' Title repeatedly checking scaffolds with two locus
#'
#' @param scaff_list character vector containing the names of scaffolds with 2 loci 
#' @param genot_data data.frame containing genotypic data
#' @param id  data.frame containing locus identifiers
#' @param type character: type of cross 
#' @param code character: code used to score data
#' @param mc character: name (in 'id') of the column containing the marker names
#' @param sc character: name (in 'id') of the column containing the scaffold names
#'
#' @return data_frame with proximity information
#' @export
#'
#' @examples
pair_check_all<-function(scaff_list,genot_data,id,type,code="abh-",mc=marker_column,sc=scaffold_column){
  out<-NULL
  for(scaff in scaff_list)
    out<-rbind(out, pair_check(scaff,genot,id,"BC1",code=code,mc,sc) )
  out
}
