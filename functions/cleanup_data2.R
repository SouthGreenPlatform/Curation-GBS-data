# eliminating outliers and redundent markers. 
# 
# the missing data in the selected loci will be imputed based on a majority vote among their descendents.


#' Title
#'
#' @param scaffold character: a scaffold
#' @param dist_data 
#' @param threshold 
#' @param large_gap 
#' 
#' @details  
#' Starting from the first marker.
#' 
#' * if it is outlier, mark it as discarded.
#' * if it is redundent and if the target is not discarded, mark it as discarded.
#' * if it is quasi-redundent, and if the target is not discarded, 
#' + if the target has more markers, mark the present marker as deleted.
#' + if not, retain the marker.
#' * if it is not redundent (*f*2>0), retain the marker. If *f*2 is above a certain value (eg 0.05), issue a warning...
#' 
#' Makes sure all loci point toward a selected marker by applyin an inheritence rule (if discarded marker 1 points toward discarded marker 2 pointing to selected marker 3, then marker 3 becomes the target of marker 1. Etc.) 
#' 
#' The missing data in the selected loci will be imputed based on a majority vote among their descendents.
#' 
#' @return
#' @export
#'
#' @examples cleanup_data("S000008",multiple_dist,5.5,0.05)  ???
cleanup_data<-function(scaffold,dist_data,threshold, large_gap){
  Subset<-dist_data[dist_data$scaffold==scaffold,]
  Subset<-Subset[order(Subset$posit),]
  Subset$select<-"_"
  Subset$parent<-""
  promise=NULL
  
  for(i in 1:nrow(Subset)){
    if(Subset$dist[i]>threshold) {
      # * if it is outlier, mark it as discarded.
      Subset$select[i]<-"outlier"
      next
    }
    if(Subset$f2[i]>0){
      # * if it is not redundent (*f*2>0), retain the marker. If *f*2 is above a certain value (eg 0.05), issue a warning...
      Subset$select[i]<-ifelse(Subset$f2[i]>large_gap,"large","small")
      next
    }
    # * if it is quasi-redundent, and if the target is not discarded, 
    # + if the target has more data, mark the present marker as ddiscarded.
    # + if not, retain the marker.
    # Ensure to have one and only one locus from a group
    value<-ifelse(Subset$dist[i]>0,"quasi-redundent","redundent")
    target_place<-as.character(Subset$locus)==as.character(Subset$target[i])
    # deal with symmetrical cases
     backward<-Subset$posit[i]>Subset$posit[target_place]
    if(Subset$signif[target_place]>=Subset$signif[i]) 
      promise<-c(promise,as.character(Subset$locus[target_place]))
    if(Subset$signif[target_place]>Subset$signif[i]) 
      promise<-setdiff(promise,as.character(Subset$locus[i]))
    cond1<-as.character(Subset$locus[i]) %in% promise
    cond2<-!as.character(Subset$locus[target_place]) %in% promise
    cond2<-cond2|backward
    if(cond1&cond2) {
      Subset$select[i]<-"chosen"
    }else{
      Subset$select[i]<-value
      Subset$parent[i]<-as.character(Subset$locus[target_place])
    }
     promise<-setdiff(promise,as.character(Subset$locus[i])) # already done?
  } # ENDFOR
  #Inheritence: all discarded loci point toward a chosen locus (except for outliers)
  chosen<-as.character(Subset$locus[Subset$select %in% c("chosen","large","small")])
  for(i in 1:nrow(Subset)){
    if(Subset$select[i] %in% c("redundent","quasi-redundent")){
      this_parent<-Subset$parent[i]
      #fun<-function(this_parent,chosen,Subset,i){
        repeat
          if(this_parent %in% chosen)
            break else 
              this_parent<- Subset$parent[as.character(Subset$locus)==this_parent] 
      #      this_parent
      #}
      #this_parent<-fun(this_parent,chosen,Subset,i)     
      #repeat if(this_parent %in% chosen) break else this_parent<- Subset$parent[as.numeric(Subset$locus[i])==this_parent]  
      Subset$parent[i]<-this_parent
    }
      
  }
  
  Subset
} #END
# example
# cleanup_data("S000008",multiple_dist,5.5,0.05)
# debug(cleanup_data)
