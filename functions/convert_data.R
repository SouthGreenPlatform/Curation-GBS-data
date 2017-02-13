#' Title convert data
#'
#' @param file.value character vector: names of the genotype data files.
#' @param file.in character vector: names of the input files 
#' @param file.ident character vector: names of the locus identifier files
#' @param where_in character: pathway to the input file.
#' @param P.value.min numeric: minimum p.value (to filter out)
#' @param where_out character: pathway to the output file.
#'
#' @return information on the data
#' @export
#'
#' @examples
convert_data <- function (file.value, file.in, file.ident, where_in,  P.value.min, where_out,long_form) {
  if(length(unique(file.value))==length(file.in)&length(unique(file.ident))==length(file.in)){
    check<-NULL  
    for(i in 1:length(file.in)){
      ftab<-read.table(topath(where_in,file.in[i]),header = TRUE)
      ncol<-length(ftab)
      loc_num<-nrow(ftab)
      mark_id<-colnames(ftab)[1:3]
      stat_id<-colnames(ftab)[(ncol-1):ncol]
      ind_id<-colnames(ftab)[4:(ncol-2)]
      keep.mark<-ftab$P.value>P.value.min
      # Rewrite values if needed
      if(as.numeric(gregexpr(":",ftab[[ind_id[1]]][1])) <1)
        value<-sapply(ind_id, function(x) ftab[[x]][keep.mark])
      else
        value<-sapply(ind_id, function(x) 
          substr(ftab[[x]][keep.mark],1,as.numeric(gregexpr(":",ftab[[x]][keep.mark]))-1))
      
      
      rownames(value)<-ftab[[mark_id[1]]][keep.mark]
      colnames(value)<-ind_id
      ident<-ftab[keep.mark,c(mark_id,stat_id)]
      scaff_nb1<-length(unique(ftab[[2]]))
      scaff_nb2<-length(unique(ident[[2]]))
      write.table(cbind(rownames(value),value),topath(where_out_1,file.value[i]),sep="\t",row.names=FALSE,quote=FALSE)
      write.table(ident,topath(where_out_1,file.ident[i]),sep="\t",row.names=FALSE,quote=FALSE)
      check<-rbind(check,cbind( before_mark=nrow(ftab),before_scaff=scaff_nb1,after_mark=nrow(value),after_scaff=scaff_nb2, individuals=ncol(value)))
      
    }#1
  }else stop("There must be as many value and ident files as input file")
  check
}
# #old version
# convert_data <- function (file.value, file.in, file.ident, where_in,  P.value.min, where_out,long_form) {
#   if(length(unique(file.value))==length(file.in)&length(unique(file.ident))==length(file.in)){
#     check<-NULL  
#     for(i in 1:length(file.in)){
#       ftab<-read.table(topath(where_in,file.in[i]),header = TRUE)
#       ncol<-length(ftab)
#       loc_num<-nrow(ftab)
#       mark_id<-colnames(ftab)[1:3]
#       stat_id<-colnames(ftab)[(ncol-1):ncol]
#       ind_id<-colnames(ftab)[4:(ncol-2)]
#       keep.mark<-ftab$P.value>P.value.min
#       # Rewrite values
#       value<-sapply(ind_id, function(x) 
#         substr(ftab[[x]][keep.mark],1,as.numeric(gregexpr(":",ftab[[x]][keep.mark]))-1))
#       
#       rownames(value)<-ftab[[mark_id[1]]][keep.mark]
#       colnames(value)<-ind_id
#       ident<-ftab[keep.mark,c(mark_id,stat_id)]
#       scaff_nb1<-length(unique(ftab[[2]]))
#       scaff_nb2<-length(unique(ident[[2]]))
#       write.table(cbind(rownames(value),value),topath(where_out_1,file.value[i]),sep="\t",row.names=FALSE,quote=FALSE)
#       write.table(ident,topath(where_out_1,file.ident[i]),sep="\t",row.names=FALSE,quote=FALSE)
#       check<-rbind(check,cbind( before_mark=nrow(ftab),before_scaff=scaff_nb1,after_mark=nrow(value),after_scaff=scaff_nb2, individuals=ncol(value)))
#       
#     }#1
#   }else stop("There must be as many value and ident files as input file")
#   check
# }