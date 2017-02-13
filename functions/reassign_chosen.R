
#' Title Reassign the parents in a group of loci of the same scaffold
#' Previous functions created groups of markers brlonging to the same scaffold and differing only by missing data.
#' One of these loci was assigned the role of "parent" (the one which had the largest number of valid data)
#' If, after data correction and imputation, one of its descendents has more valid data, this function will assign it as the new parent 
#' 
#' @author Luc Baudouin
#'  
#'  
#' @param chosen_locus Character (atomic): a locus 
#' @param corrected_data  data.frame: the genotypic values 
#' @param locus_info  data.frame: tle locus information data
#'
#' @return data.frame: a new version of the locus_info data with two added columns for the locus group
#' @export
#'
#' @examples
reassign_chosen <- function(chosen_locus,corrected_data, locus_info,code="ABHU") {
  code_1<-strsplit(code,split="")[[1]]
  fam <-
    c(as.character(chosen_locus),as.character(selected$locus[selected$parent == chosen_locus]))
  data_fam <- as.matrix(corrected_data[fam,])
  info_fam <- locus_info[locus_info$locus %in% fam,]
  place <- match("chosen",info_fam$select)
  new_signif <- sapply(1:nrow(data_fam), function(k1) length(data_fam[k1,]) - sum(toupper(data_fam[k1,]) == code_1[4]))
  if (max(new_signif) > new_signif[place]) {
    info_fam$swap <- "yes"
    new_place <- which.max(new_signif)
    info_fam$select[new_place] <- "chosen"
    info_fam$select[place] <- "swaped"
    new_parent <- info_fam$locus[new_place]
    info_fam$parent <- new_parent
    info_fam$parent[new_place] <- NA
    
  }else
    info_fam$swap <- "no"
  cbind(info_fam,new_signif)
}
