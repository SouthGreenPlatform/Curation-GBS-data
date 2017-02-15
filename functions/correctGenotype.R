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

