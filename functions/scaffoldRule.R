#' Title Defines correction rules for scaffolds
#'
#' @return character matrix: the rules
#' @export
#'
#' @examples

scaffoldRules<-function(){
  
  rule<-NULL
  
  # uniform sequences except sporadic errors or missing data. Problems are corrected.
  # At most two 
  # example: "AAAHAAABAAAAAABBAAAAA" becomes "AAAaAAAaAAAAAAaaAAAAA"
  rule<-rbind(rule,cbind(query="A{2,}([BH]{1,2}A{2,})*",
                         problem="[BH]",
                         correct="a")) 
  rule<-rbind(rule,cbind(query="B{2,}([AH]{1,2}B{2,})*",
                         problem="[AH]",
                         correct="b"))
  rule<-rbind(rule,cbind(query="H{2,}([AB]{1,2}H{2,})*",
                         problem="[AB]",
                         correct="h"))
  # Easier rules for missing data
  # example: "AAAUAAAUAAAAAAUUA" becomes "AAAaAAAaAAAAAAaaA"
  
  rule<-rbind(rule,cbind(query="A{1,}(U{1,3}A+)+",
                         problem="U",
                         correct="a")) 
  rule<-rbind(rule,cbind(query="B{1,}(U{1,3}B+)+",
                         problem="U",
                         correct="b"))
  rule<-rbind(rule,cbind(query="H{1,}(U{1,3}H+)+",
                         problem="U",
                         correct="h"))
  # interferences between errors and missing data case 1
  # example: "AAUHAA" becomes "AAaaAA
  
  rule<-rbind(rule,cbind(query="A{2,}(U[BH]A{2,})+",
                         problem="U[BH]",
                         correct="aa")) 
  rule<-rbind(rule,cbind(query="B{2,}(U[AH]B{2,})+",
                         problem="U[AH]",
                         correct="bb")) 
  rule<-rbind(rule,cbind(query="H{2,}(U[AB]A{2,})+",
                         problem="U[AB]",
                         correct="hh")) 

  # interferences between errors and missing data case 2
  # example: "AAHUAA" becomes "AAaaAA
  
  rule<-rbind(rule,cbind(query="A{2,}([BH]UA{2,})+",
                         problem="[BH]U",
                         correct="aa")) 
  rule<-rbind(rule,cbind(query="B{2,}([AH]UB{2,})+",
                         problem="[AH]U",
                         correct="bb")) 
  rule<-rbind(rule,cbind(query="H{2,}([AB]UA{2,})+",
                         problem="[AB]U",
                         correct="hh")) 
  
  # interferences between errors and missing data case 3
  # example: "AAUAHAA" becomes "AAaAaAA
  
  rule<-rbind(rule,cbind(query="A{2,}(UA[BH]A{2,})+",
                         problem="UA[BH]",
                         correct="aAa")) 
  rule<-rbind(rule,cbind(query="B{2,}(UB[AH]B{2,})+",
                         problem="UB[AH]",
                         correct="bBb")) 
  rule<-rbind(rule,cbind(query="H{2,}(UH[AB]A{2,})+",
                         problem="UH[AB]",
                         correct="hHh")) 
  
  # interferences between errors and missing data case 4
  # example: "AAUAHAA" becomes "AAaAaAA
  
  rule<-rbind(rule,cbind(query="A{2,}([BH]AUA{2,})+",
                         problem="[BH]AU",
                         correct="aAa")) 
  rule<-rbind(rule,cbind(query="B{2,}([AH]BUB{2,})+",
                         problem="[AH]BU",
                         correct="bBb")) 
  rule<-rbind(rule,cbind(query="H{2,}([AB]HUA{2,})+",
                         problem="[AB]HU",
                         correct="hHh")) 
  # Ambiguous sequence
  # Rationale "AABABB" could be "AAaABB" or"AABbBB". Therefore, it is converted to "AAuuBB" 
  rule<-rbind(rule,cbind(query="ABAB",
                         problem="BA",
                         correct="uu")) 
  rule<-rbind(rule,cbind(query="AHAB",
                         problem="HA",
                         correct="uu")) 
  rule<-rbind(rule,cbind(query="AHAH",
                         problem="HA",
                         correct="uu")) 
  rule<-rbind(rule,cbind(query="ABAH",
                         problem="BA",
                         correct="uu"))
  
  rule<-rbind(rule,cbind(query="BABA",
                         problem="AB",
                         correct="uu")) 
  rule<-rbind(rule,cbind(query="BHBA",
                         problem="HB",
                         correct="uu")) 
  rule<-rbind(rule,cbind(query="BABH",
                         problem="AB",
                         correct="uu")) 
  rule<-rbind(rule,cbind(query="BHBH",
                         problem="HB",
                         correct="uu")) 
  
  rule<-rbind(rule,cbind(query="HABA",
                         problem="AB",
                         correct="uu")) 
  rule<-rbind(rule,cbind(query="HBHA",
                         problem="BH",
                         correct="uu")) 
  rule<-rbind(rule,cbind(query="BABH",
                         problem="AB",
                         correct="uu")) 
  rule<-rbind(rule,cbind(query="BHBH",
                         problem="HB",
                         correct="uu")) 
  rule
}

bid<-function(){
  #Pour tester les règles...
  # tandem error. Errors are corrected  
  # example: "000110900000022022000" becomes "000000000000000000000" 
  rule<-rbind(rule,cbind(query="A{2,}[BH]A{1,2}[BH]A{2,}" ,
                         problem="[BH]",
                         correct="a"))
  rule<-rbind(rule,cbind(query="B{2,}[AH]B{1,2}[AH]B{2,}" ,
                         problem="[AH]",
                         correct="b"))
  rule<-rbind(rule,cbind(query="H{2,}[AB]H{1,2}[AB]H{2,}" ,
                         problem="[AB]",
                         correct="h"))  
  # only one missing value, not taken in charge by previous rules. isolated missing data corrected
  # example:  "000900999990000" becomes "000000999990000"
  
  rule<-rbind(rule,cbind(query="A{2,}(UA)*A{1,}" ,
                         problem="AUA",
                         correct="AaA"))
  rule<-rbind(rule,cbind(query="B{2,}(UB)*B{1,}" ,
                         problem="BUB",
                         correct="BbB"))
  rule<-rbind(rule,cbind(query="H{2,}(UH)*H{1,}" ,
                         problem="HUH",
                         correct="HhH"))
  
  # same as before but two missing values
  # example:  "000990999990000" becomes "000000999990000"
  rule<-rbind(rule,cbind(query="A{2,}(UUA)*A{1,}" ,
                         problem="AUUA",
                         correct="AaaA"))
  rule<-rbind(rule,cbind(query="B{2,}(UUB)*B{1,}" ,
                         problem="BUUB",
                         correct="BbbB"))
  rule<-rbind(rule,cbind(query="H{2,}(UUH)*H{1,}" ,
                         problem="HUUH",
                         correct="HhhH"))
}