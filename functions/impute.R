#' Title impute missing data in "ancestor" loci
#'
#' @param locus character: a locus name
#' @param genot_data data.frame: the genotypic data
#' @param select_file data.frame: information on the loci
#' @param type character: the type of cross (only "BC1" for the moment )
#' @param accept numeric: majority criterion if more than one value is proposed, the most frequent is accepted only if its frequency is above 'accept
#' @param code character: code used for data encoding. The default is ok for joinmap
#'
#' @return the genotypes at a locus, after imputation
#' @export
#'
#' @examples#' 
impute_x<-function(locus, genot_data,select_file,type="BC1",accept=0.8,code="ABHU"){
  codes<-strsplit(code,"")[[1]]
  this_data<-as.character(genot_data[locus,])
  #pattern<-rep("0",length(this_data))
  descend<-as.character(select_file$locus[select_file$parent==locus]  )
  descend_data<-genot_data[descend,]
  miss<-grep(codes[4],toupper(this_data))
  # next, no missing data
  if(length(miss)==0) return(this_data)
  # ?t least one missing value
  #pattern[miss]<-1
  out1<-this_data[miss]
  correct<-this_data
  if(length(miss)==1){
    out2<-impute_1(this_data[miss],toupper(descend_data[,miss]))
  }else{
    #Several missing values
    out2<-sapply(1:length(out1), function(i) impute_1(out1[i],toupper(descend_data[,miss][[i]])))
  }
  correct[miss]<-tolower(out2)
  return(correct)
}

#' Title corrects one value
#'
#' @param current character: missing genotypic value to be imputed (in the ancestor)
#' @param character vector: Source the genotypes of the descendent loci for the present individual
#' @param accept numeric: majority criterion (see assign_x)
#' @param code 
#'
#' @return character: genotypic value after imputation
#' @export
#'
#' @examples
impute_1<-function(current,Source,accept=0.8,code="ABHU"){
  # does just one correction
  # if descendents have just one value, return this value
  # if just one value except missing, return this value
  # if more than one value and one has more than 'accept' * number of (non missing) values, return it 
  # otherwise, return a missing value
  codes<-strsplit(code,"")[[1]]
  size<-length(Source)
  T<-table(Source)
  if(length(T)==1) return(names(T))
  T2<-T[names(T) !=codes[4]]
  if(length(T2)==1) return(names(T2))
  score=max(T2)/sum(T2)
  if(score>=accept) return(names(T)[which.max(T2)])
  return(current)
}
