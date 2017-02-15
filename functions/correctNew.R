#passage à la F2
# ind<-pop[4:5,]
# dim(ind)
# toStrings(ind)
# 

#Je m'inspire de joinmap mais je prends des capitales. Les corrections seront en minuscule
#On peut toujours revenir aux minuscules... avec tolower.

#Voir aussi ceci qui rend la fonction inutile
#> chartr("abhu","DTHM",c("a","b","b","u","ab","z"))
#[1] "D"  "T"  "T"  "M"  "DT" "z" 

codePop<-function(p,code="ABHU"){
  # conversion of simulated data into string form (locus level) 
  # Applies to F2 populations
  code_1<-strsplit(code,split="")[[1]]
  ni<-floor(nrow(p)/2)
  sapply(1:ni, function(x)
    paste(
      sapply(1:ncol(ind), function(y)
        if(is.na(p[2*x-1,y])||is.na(p[2*x,y])) return(code_1[4])
        else if(p[2*x-1,y]&&p[2*x,y]) return(code_1[3]) 
        else if(!p[2*x-1,y]&&!p[2*x,y]) return(code_1[1]) 
        else return(code_1[2]))
      ,collapse="")
  )
}


mappingRules<-function(){
  #to be use for mapping
  # This is the F2 version (should work also for backcrosses!)
  # code: 
  #   A first parent (recurrent parent in a BC)
  #   B second parent (heterozygote in a BC)
  #   H heterozygote (not used in BC)
  #   U missing genotype
  #########################################################################
  # 'query' is the pattern collected in an individual genotype (a string)
  # the outcome is a string which may contain unwanted character
  # they are identified by 'problem' and replaced by 'correct' (see examples)
  #########################################################################
  # to produce alternative rule sets, just add/remove/alter rules
  # IMPORTANT: the length of the strings matched by 'problem' should 
  # have a fixed length, equal to 'length(correct)'
  # Remember that the order of rules counts!!!
  # To test new rules, use the following line as a template 
  # zzz<-zz<-"00090099999000990999999000";regmatches(zz,gregexpr("0990",zz,perl=TRUE)) <- "0000";rbind(zzz,zz)
  ########################################################################
  # new version with shorter missing sequences
  #########################################################################
  rule<-NULL

  # uniform sequences except sporadic errors or missing data. Problems are corrected.
  # Basically, this is the four-point method
  # example: "00010002000000990000" becomes "00000000000000000000"
  rule<-rbind(rule,cbind(query="A{3,}([BHU]{0,2}A{3,})*",
                         problem="[BHU]",
                         correct="a")) 
  rule<-rbind(rule,cbind(query="B{3,}([AHU]{0,2}B{3,})*",
                         problem="[AHU]",
                         correct="b"))
  rule<-rbind(rule,cbind(query="H{3,}([ABU]{0,2}H{3,})*",
						problem="[ABU]",
						correct="h"))
  
  # tandem error. Errors are corrected  
  # example: "000110900000022022000" becomes "000000000000000000000" 
  rule<-rbind(rule,cbind(query="A{3,}[BH]A{1,2}[BH]A{3,}" ,
                         problem="[BH]",
                         correct="a"))
  rule<-rbind(rule,cbind(query="B{3,}[AH]B{1,2}[AH]B{3,}" ,
                         problem="[AH]",
                         correct="b"))
  rule<-rbind(rule,cbind(query="H{3,}[AB]H{1,2}[AB]H{3,}" ,
                         problem="[AB]",
                         correct="h"))
  
  # only one missing value, not taken in charge by previous rules. isolated missing data corrected
  # example:  "000900999990000" becomes "000000999990000"
  
  rule<-rbind(rule,cbind(query="A{3,}(UA)*A{2,}" ,
                         problem="AUA",
                         correct="AaA"))
  rule<-rbind(rule,cbind(query="B{3,}(UB)*B{2,}" ,
                         problem="BUB",
                         correct="BbB"))
  rule<-rbind(rule,cbind(query="H{3,}(UH)*H{2,}" ,
                         problem="HUH",
                         correct="HhH"))

  # same as before but two missing values
  # example:  "000990999990000" becomes "000000999990000"
  rule<-rbind(rule,cbind(query="A{3,}(UUA)*A{2,}" ,
                         problem="AUUA",
                         correct="AaaA"))
  rule<-rbind(rule,cbind(query="B{3,}(UUB)*B{2,}" ,
                         problem="BUUB",
                         correct="BbbB"))
  rule<-rbind(rule,cbind(query="H{3,}(UUH)*H{2,}" ,
                         problem="HUUH",
                         correct="HhhH"))
  
  rule
}
# see independent file
# scaffoldRules<-function(){
#   #to be use for scaffolds
#   # This is the F2 version (should work also for backcrosses!)
#   # code: 
#   #   A first parent (recurrent parent in a BC)
#   #   B second parent (heterozygote in a BC)
#   #   H heterozygote (not used in BC)
#   #   U missing genotype
#   #########################################################################
#   # 'query' is the pattern collected in an individual genotype (a string)
#   # the outcome is a string which may contain unwanted character
#   # they are identified by 'problem' and replaced by 'correct' (see examples)
#   #########################################################################
#   # to produce alternative rule sets, just add/remove/alter rules
#   # IMPORTANT: the length of the strings matched by 'problem' should 
#   # have a fixed length, equal to 'length(correct)'
#   # Remember that the order of rules counts!!!
#   # To test new rules, use the following line as a template 
#   # zzz<-zz<-"00090099999000990999999000";regmatches(zz,gregexpr("0990",zz,perl=TRUE)) <- "0000";rbind(zzz,zz)
#   ########################################################################
#   # new version with shorter missing sequences
#   #########################################################################
#   rule<-NULL
# 
#   # uniform sequences except sporadic errors or missing data. Problems are corrected.
#   # At least two 
#   # example: "00010002000000990000" becomes "00000000000000000000"
#   rule<-rbind(rule,cbind(query="A{2,}([BHU]{0,2}A{2,})*",
#                          problem="[BHU]",
#                          correct="a")) 
#   rule<-rbind(rule,cbind(query="B{2,}([AHU]{0,2}B{2,})*",
#                          problem="[AHU]",
#                          correct="b"))
#   rule<-rbind(rule,cbind(query="H{2,}([ABU]{0,2}H{2,})*",
# 						problem="[ABU]",
# 						correct="h"))
#   
#   # tandem error. Errors are corrected  
#   # example: "000110900000022022000" becomes "000000000000000000000" 
#   rule<-rbind(rule,cbind(query="A{2,}[BH]A{1,2}[BH]A{2,}" ,
#                          problem="[BH]",
#                          correct="a"))
#   rule<-rbind(rule,cbind(query="B{2,}[AH]B{1,2}[AH]B{2,}" ,
#                          problem="[AH]",
#                          correct="b"))
#   rule<-rbind(rule,cbind(query="H{2,}[AB]H{1,2}[AB]H{2,}" ,
#                          problem="[AB]",
#                          correct="h"))
#   
#   # only one missing value, not taken in charge by previous rules. isolated missing data corrected
#   # example:  "000900999990000" becomes "000000999990000"
#   
#   rule<-rbind(rule,cbind(query="A{2,}(UA)*A{1,}" ,
#                          problem="AUA",
#                          correct="AaA"))
#   rule<-rbind(rule,cbind(query="B{2,}(UB)*B{1,}" ,
#                          problem="BUB",
#                          correct="BbB"))
#   rule<-rbind(rule,cbind(query="H{2,}(UH)*H{1,}" ,
#                          problem="HUH",
#                          correct="HhH"))
# 
#   # same as before but two missing values
#   # example:  "000990999990000" becomes "000000999990000"
#   rule<-rbind(rule,cbind(query="A{2,}(UUA)*A{1,}" ,
#                          problem="AUUA",
#                          correct="AaaA"))
#   rule<-rbind(rule,cbind(query="B{2,}(UUB)*B{1,}" ,
#                          problem="BUUB",
#                          correct="BbbB"))
#   rule<-rbind(rule,cbind(query="H{2,}(UUH)*H{1,}" ,
#                          problem="HUUH",
#                          correct="HhhH"))
#   
#   rule
# }


correctGenotype <- function (noisy,ruleSet,code) {
  # create a corrected version of noisy genotyping data. 
  # see 'correctRules()'
  # 
  # 
  if(!missing(code))  noisy<-transcode(noisy,code,"ABHU")
  rule<-ruleSet()
  # initialize with original data  
  corr<-noisy
  for(i in 1:nrow(rule)) {
    found <- gregexpr(rule[i,"query"],noisy,perl=TRUE)
    matches <- regmatches(noisy,found)   # identified patterns
    condition <- sapply(matches,length) != 0
    if (sum(condition) > 0) {
      for (x in 1:length(matches[condition]))
        # corrects problems
        regmatches(matches[condition][[x]],gregexpr(rule[i,"problem"],matches[condition][[x]],perl=TRUE)) <- rule[i,"correct"]
      # replace identified patterns with corrected version
      regmatches(corr,found) <- matches
    }
  }
  corr
}


diagnostic<-function(origin,noisy,corr){
  # used to test rules for correctGenotype
  # identifies remaining unsolved cases ("9")
  filter<-grepl("9+",corr)
  u<-lapply(1:length(origin), function(x) rbind(true=origin[x],noisy=noisy[x],corrected=corr[x]))
  u[filter]
}

unknown<-function(corr){
  # used to test rules for correctGenotype
  # calculates the percentage of unsolved cases ("9")
  unk<-regmatches(corr,gregexpr("9+",corr))
  sunk<-sum(sapply(unk,function(x) sum(nchar(x))))
  spos<-sum(nchar(corr))
  paste(sunk/spos*100,"% unknown genotypes")
}

# compareGenot<-function(n){
# # for debugging
#     cat("true     ",genot[n],"\n")
#   cat("noisy    ",observ[n],"\n")
#   cat("corrected",correct[n],"\n")
#   cat("\n")
# }


transcode<-function(Data,code1, code2){
  # traduit des données avec un codage différent.
  # les fonctions utilisent un codage "0129" 
  # code1 et code2 sont des chaines de caractères
  # ordre: parent 1, parent 2, heterozygote, manquant
  # cas des backcrosses, le 2e code est fourni mais  utilisé
  if(nchar(code1)!=nchar(code2)) stop("codes should have the same number of elements")
  code_1<-strsplit(code1,split="")[[1]]
  code_2<-strsplit(code2,split="")[[1]]
  for(i in 1:length(code_1)){
    regmatches(Data,gregexpr(code_1[i],Data))<-code_2[i]
  }
  Data
}

findCO2 <-function (corr) {
  # locates crossing overs (CO) to intervals with missing data (code "9")
  # in case of missing data, intervals are assigned equally distributed CO probability.
  # This is the F2 version (should work also for backcrosses!)
  pat<-"09*1|19*0|19*2|29*1"
  a<-gregexpr(pat,corr,perl=TRUE)
  out1<-NULL
  for(i in 1:length(corr)){
    zz<-rep(0,nchar(corr[1])-1)
    if(a[[i]][1]!=-1){
      beg<-as.numeric(a[[i]])
      end<-beg+attr(a[[i]], "match.length")-2
      for(j in 1:length(a[[i]])) zz[beg[j]:end[j]]<-1/(end[j]-beg[j]+1)
    }
    out1<-rbind(out1,zz)
    if(ncol(out1)!=length(zz)) cat(c(i,length(zz),ncol(out1),"\n"))
  }
  out1
}


#faire une seule fonction avec les 2 suivantes
testPhase<-function(l1,l2,code="0129"){
  # Base function for testing phases
  # l1 and l2 is the score at consecutive loci for an individual
  # code is four character: Parent 1, parent 2, heterozygote, missing data  
  # output:1 = undetermined, 2=identical, 3=permitted, 4=forbidden
  code_1<-strsplit(code,split="")[[1]]
  if(!((l1 %in% code_1)&(l1 %in% code_1))) stop("values should be in '", code,"'")
  if(code_1[4] %in% c(l1,l2)) return(1)
  if(l1==l2) return(2)
  if(code_1[3] %in% c(l1,l2)) return(3)
  return(4)
}


testPhasePop<-function(genot,tol=0.02,...){
  # population level function for testing phases
  # individuals as rows, loci as columns
  # tol is there to allow F<P  if F<= (P+I)/tol
  out0<-NULL
  for(j in 2:ncol(genot)){
    ##########################  
    out<-c(U=0,I=0,P=0,F=0,S=0,G=0)
    for(i in 1:nrow(genot)) {
      n<-testPhase(genot[i,j-1],genot[i,j],...)
      out[n]<-out[n]+1
    }
    maxF<-(out[3]+out[2])*tol
    out[5]<-ifelse(out[4]>max(out[3],maxF),1,0)
    ###############################
    out0<-rbind(out0,out)
  }
  for(i in 2:nrow(out0)) out0[i,6]<-ifelse(out0[i-1,6]==out0[i,5],0,1)
  rownames(out0)<-NULL
  out0
}


plotPop2 <- function (pop,code="abhm",pal=c("magenta","cyan","plum","peachpuff4"),...) {
# New version of plotPop
# the matrix ix individual*locus  
# Use a four letter code (parent 1, parent2, heterozygote and missing)
    code_1<-strsplit(code,split="")[[1]]
  
  gr<-data.frame(ind=rep(1:nrow(pop),times=ncol(pop)),
                 loc=rep(1:ncol(pop),each=nrow(pop)),
                 gen=as.vector(pop),
                 h=rep(1,nrow(pop)*ncol(pop)),
                 v=rep(1,nrow(pop)*ncol(pop)))
  gra<-gr[gr$gen==code_1[1]&!is.na(gr$gen),]
  grb<-gr[gr$gen==code_1[2]&!is.na(gr$gen),]
  grh<-gr[gr$gen==code_1[3]&!is.na(gr$gen),]
  grm<-gr[gr$gen==code_1[4]&!is.na(gr$gen),]
  nrow(gra)+nrow(grb)+nrow(grh)
  plot(0,type="n",ylim=c(0,ncol(pop)),xlim=c(0,nrow(pop)),xlab="individual",ylab="locus",frame.plot=FALSE,...)
  rect(gra$ind-1,gra$loc-1,gra$ind,gra$loc,border=NA,col=pal[1]) 
  rect(grb$ind-1,grb$loc-1,grb$ind,grb$loc,border=NA,col=pal[2]) 
  rect(grh$ind-1,grh$loc-1,grh$ind,grh$loc,border=NA,col=pal[3]) 
  rect(grm$ind-1,grm$loc-1,grm$ind,grm$loc,border=NA,col=pal[4]) 
}

locStat<-function(genotM,loc) {
  ###########################################
  # statistiques élémentaires sur les données
  code="abhm"
  code_1<-strsplit(code,split="")[[1]]
  z<-apply(genotM,2, table)
  if(class(z)=="list"){
  na<-sapply(z,function(x) x[code_1[1]])
  nb<-sapply(z,function(x) x[code_1[2]])
  nh<-sapply(z,function(x) x[code_1[3]])
  nm<-sapply(z,function(x) x[code_1[4]])
  }else {
    # z is a matrix
    na=z[code_1[1],]
    nb=z[code_1[2],]
    nh=z[code_1[3],]
    nm=z[code_1[4],]
  }
  out<-data.frame(loc,na,nb,nh,nm)
  out[2:4][is.na(out[2:4])]<-0
  out$obs<-out$na+out$nb+out$nh
  out$miss<-out$nm/(out$obs+out$nm)
  out$fa<-out$na/out$obs
  out$fb<-out$nb/out$obs
  out$fh<-out$nh/out$obs
  out$p.value<-sapply(1:nrow(out), function(x)
    chisq.test(as.numeric(out[x,2:4]),p=c(1/4,1/4,1/2))$p.value)
  out
}


findModif<-function(first,second,code="abhm"){
# Changes made on the genotype data
  code_1<-strsplit(code,split="")[[1]]
  sapply(1:nrow(first), function(x)
    sapply(1:ncol(first), function(y){
      if(code_1[4]==second[x,y]) return("U") 
      if(code_1[4]==first[x,y]) return("M")
      if(first[x,y]==second[x,y]) return("I")
      return("C")
    }))
}


summaryModif<-function (genotM, correctg,...) {
  #summary, locus by locus of changes 
  # will serve as a basis to evaluate error and missing data rate.
  #needs improving
  # I: no change
  # U: genotype set to missing 
  # M genotype was missing
  # C: the genotype changed
  diag<-findModif(t(genotM),t(correctg))
  diagsum<-t(apply(diag,2,function(x) c(I=length(x[x=="I"]),U=length(x[x=="U"]),M=length(x[x=="M"]),C=length(x[x=="C"]))))
  diagsum<-data.frame(diagsum)
  rownames(diagsum)<-colnames(correctg)
  diagsum
  
}
# testPhaseInt(genotM[1:2,],"abhm")

# xx<-testPhasePop(genotM,"abhm")

# testPhase("a","a","abhm")
# testPhase("a","b","abhm")
# testPhase("a","h","abhm")
# testPhase("a","m","abhm")
# testPhase("b","a","abhm")
# testPhase("b","b","abhm")
# testPhase("b","h","abhm")
# testPhase("b","m","abhm")
# testPhase("h","a","abhm")
# testPhase("h","b","abhm")
# testPhase("h","h","abhm")
# testPhase("h","m","abhm")
# testPhase("m","a","abhm")
# testPhase("m","b","abhm")
# testPhase("m","h","abhm")
# testPhase("m","m","abhm")
# 
# testPhase("0","0")

#  codeLoc(TRUE,TRUE)
#  codeLoc(TRUE,FALSE)
#  codeLoc(FALSE,TRUE)
#  codeLoc(FALSE,FALSE)
#  codeLoc(TRUE,NA)
#  codeInd<-function(indDat)
#  sapply(1:ncol(ind), function(x) codeLoc(indDat[1,x],indDat[2,x]))
#  
#  codeInd(ind)
#  
# genot<-codePop(pop)
# observ<-codePop(popN)
# correct<-correctGenotype2(observ) 
# correct
# compareGenot(1)
# bid<-sapply(1:10,compareGenot)
# observ3<-transcode(observ,"0129","abcd") 
# 
# rm("observ2")
# debug(transcode)
# test<-"0000111122229111999000222"
# pat<-"09*1|19*0|19*2|29*1"
# regmatches(test,gregexpr(pat,test))
# 
# 
# before<-findCO2(genot)
# xb<-apply(before,2,sum)
# after<-findCO2(correct)
# xa<-apply(after,2,sum)
# noisy<-findCO2(observ)
# xn<-apply(noisy,2,sum)
# 
# loc2<-(loc-loc[1])[-1]
# loc2<-loc2*L/L0 # les abscisses "physiques" des locus
# plot(loc2,type="l",ylim=c(0,150),xlab="loci",ylab="distance(cM)", main ="actual and estimated distances",frame.plot=FALSE)
# lines(cumsum(xb/N*100) ,col="orange")# crossing-overs sans bruit
# lines(cumsum(xn/N*100),col="red") # crossing-overs avec bruit
# lines(cumsum(xa/N*100),col="cyan") # données corrigées
# legend("bottomright",c("actual","estimated without errors", "estimated without corrections","estimated with corrections" ), lty=rep(1,4),col=c("black","orange","red","cyan"),inset=0.02,cex=.7)
# abline(h=d*(L-1),lty=3)
