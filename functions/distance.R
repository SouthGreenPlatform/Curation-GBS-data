distance<-function(gmark1, gmark2,type,code="abh-"){
  # compute three kinds of recombination distances
  cv<- strsplit(code,split="")[[1]]
  tc<-grep(type,c("BC1","F2"))
  if( length(tc)==0) stop( "type must be in c('BC1','F2').")
  if(length(gmark1)!=length(gmark2)) stop("genotypic vectors should have the same length")
  v1<-factor(as.character(gmark1),cv)
  v2<-factor(as.character(gmark2),cv)
  if(type=="BC1"){
    enrec<-c(0,0,2,1,0,0,0,0,2,0,0,1,1,0,1,1)/2
    filter<-c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0)
  }
  if(type=="F2") {
    enrec<-c(0,2,1,1,2,0,1,1,1,1,1,1,1,1,1,1)
    filter<-c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
  }
  gg<-as.numeric(table(data.frame(v1,v2)))
  exp_n_rec<-sum(gg*enrec)
  f1<-exp_n_rec/sum(gg)
  signif<-sum(gg*filter)
  f2<-sum(gg*enrec*filter)/signif
  data.frame(exp_n_rec,f1,f2,signif)
}
