#' Title select individuals to be removed and monitor the process 
#'
#' @param gdata data.frame: genotypic data    
#' @param idata data.frame: Information about the loci    
#' @param threshold numeric: the maximum acceptable proportion of missing data for a locus    
#' @param stop.at numeric: the number of steps to be executed    
#' @param code character: character: a 4 letter code for genotypes    
#' 
#' @details a recursive functions which helps selecting an
#'  optimal trade-off between number of individuals and number of loci.
#' @return a data.frame with 6 columns:     
#'  - size is the number of remaining individuals at the current step,     
#'  - ind1 is the maximum of the percentage of missing data for the individual with all loci,    
#'  - ind2 is the maximum of the percentage of missing data for the individual with the selected loci,   
#'  - loc is the number of loci that meet the 'threshold' criterion,    
#'  - scaff is the number of scaffolds represented by 'loc'    
#'  - remove is the name of the individual that will be removed at the *next* step    
#' 
#' @export
#'
#' @examples
select_ind <- function(gdata,idata,threshold = .2,stop.at = 3,code = "ABHU") {
  code_U <- strsplit(code,split = "")[[1]]
  code_L <- tolower(code_U)
  val1 <- sort(c(code_L, code_U))
  unknown <- c(code_L[4],code_U[4])
  zind <-
    t(apply(gdata,2,function(y)
      apply(sapply(val1, function(x)
        match(y, x,nomatch = 0)),2,sum)))
  zloc <-
    t(apply(gdata,1,function(y)
      apply(sapply(val1, function(x)
        match(y, x,nomatch = 0)),2,sum)))
  cmiss.ind <- apply(zind[,unknown],1,sum)
  cmiss.loc <- apply(zloc[,unknown],1,sum)
  pmiss.ind <- cmiss.ind / nrow(gdata)
  pmiss.loc <- cmiss.loc / ncol(gdata)
  p.keep.loc <- pmiss.loc[pmiss.loc < threshold]
  remove <- names(cmiss.ind)[which.max(cmiss.ind)]
  keep.ind <- setdiff(names(cmiss.ind),remove)
  keep.scaff <-
    unique(idata$scaffold[idata$locus %in% names(p.keep.loc)])
  zind2 <-
    t(apply(gdata[names(p.keep.loc),keep.ind],2,function(y)
      apply(sapply(val1, function(x)
        match(y, x,nomatch = 0)),2,sum)))
  cmiss2.ind <- apply(zind2[,c("u","U")],1,sum)
  pmiss2.ind <- cmiss2.ind / length(p.keep.loc)
  out <-
    data.frame(
      size = ncol(gdata),ind1 = round(max(pmiss.ind),2),
       ind2 = round(max(pmiss2.ind),2),loc = length(p.keep.loc),
        scaff = length(keep.scaff),remove)
  if (stop.at == 1)
    return(out)
  return(rbind(out,Recall(
    good_data[,keep.ind],idata,threshold, stop.at - 1,code
  )))
}
