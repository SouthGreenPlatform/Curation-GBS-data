# extracts a specified subset of dataset data
# The condition applies to a locus but the whole scaffold is extracted.
subset_data<-function(loc_dist,genot_data,threshold,up=TRUE,distance="dist"){
  if(up) ld<-loc_dist[loc_dist[[distance]]>threshold,] else ld<-loc_dist[loc_dist[[distance]]<=threshold,]
  dist<-ld[[distance]]
  scaffold<-unique(ld$scaffold)
  ld1<-loc_dist[loc_dist$scaffold %in% scaffold,]
  out<-cbind(ld1[-2],genot_data[as.character(ld1$locus),])
  rownames(out)<-ld1[[2]]
  #sort(rownames(out))
  #out[sort(rownames(out)),]
  out
}