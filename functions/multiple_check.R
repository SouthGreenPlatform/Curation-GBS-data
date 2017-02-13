#working with multiple
#' Title
#'
#' @param scaff character: name of a scaffold with several loci 
#' @param genot_data data.frame containing genotypic data
#' @param id data.frame containing locus identifiers
#' @param type character: type of cross
#' @param code character: code used to score data
#'
#' @return data_frame with proximity information (one line per locus)
#' @export
#'
#' @examples
multiple_check<-function(scaff,genot_data,id,type,code="abh-",mc=marker_column,sc=scaffold_column){
  setloc<-id[marker_column][id[scaffold_column]==scaff]
  nloc<-length(setloc)
  g<-genot_data[setloc,]
  out<-NULL
  for(i in 1:nloc){
    dmin<-2*ncol(g)
    for(j in (1:nloc)[-i]){
      d<-distance(g[i,],g[j,],type,code)
      if(d[[1]]<dmin){
        dmin<-d[[1]]
        f1<-d[[2]]
        f2<-d[[3]]
        signif<-d[[4]]
        target=setloc[j]
      }
      #dmin<-min(dmin,distance(g[i,],g[j,],type,code)[[1]])
    }
    out<-rbind(out,data.frame(scaffold=scaff,locus=setloc[i],posit=posit_loc(setloc[i]),target,dist=dmin,f1,f2,signif))
  }
  out
}
#example
# multiple_check("S000008",genot,id,"BC1",code="DTHM")