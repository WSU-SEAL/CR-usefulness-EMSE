source("common.R")

pValue_extract <- function(x){
  z <- summary(x)$coefficients/summary(x)$standard.errors
  # 2-tailed Wald z tests to test significance of coefficients
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  p
}

#print(datafile)
  datasetCR = read.csv("./usefulness_dataset.csv", header = TRUE)
  
  
 
  datasetCR$log_RI =log(datasetCR$RI)
  datasetCR$log_CCR =log(datasetCR$CCR)
  datasetCR$log_RCE =log(datasetCR$RCE+1)
  datasetCR$log_MR =log(datasetCR$MR+1)
  datasetCR$log_PRC =log(datasetCR$PRC+1)
  datasetCR$log_PFRC =log(datasetCR$PFRC+1)
  datasetCR$log_CL =log(datasetCR$CL)
  datasetCR$log_COMC =log(datasetCR$reviewer_made_commit+1)
  
  
  dynamic_vars = c(# participant attributes
    'WRC',
    'log_MR',
    'RPT',
    'log_PFRC',
    'log_COMC',
    'RCS',
    'RCSH',
    'log_RCE',
    'log_PRC',
    # contextual attributes
    'PN',
    'CV',
    'FUR',
    'DUR',
    'log_RI',
    'log_CCR',
    'IBF',
    'INF',
    'CCY',
    'log_CL',
    'PDR',
    'PDL' )
  

 

  survived_vars = AutoSpearman(dataset = datasetCR,
                                   metrics = dynamic_vars,verbose = T)
  print(survived_vars)
  
  
  #Plot hierarchical clusters and the spearman's correlation threshold of 0.7
  vc <- varclus(~ ., data=datasetCR[,dynamic_vars], trans="abs")
   plot(vc)
  threshold <- 0.7
  abline(h=1-threshold, col = "red", lty = 2)
  

  #Create formula string
  formula_string ="sample_rating ~  IBF "
    
  for (variable in survived_vars){
      if (variable !="IBF")
        formula_string =paste(formula_string, "+ ", variable,  sep=" ")
  
  }
    
  print (formula_string )
  
  model_fit <- glm(as.formula(formula_string) , data=datasetCR,  x=T, y=T  )
  model_summary=summary(model_fit)

 
print(model_summary)

beta_values =coef(model_fit)
conf_int =confint(model_fit)

variable_power= Anova(model_fit, test="Wald", type ="II")
print(variable_power)

MLR_LATEX =""


for (variable in survived_vars){
  
  beta =formatC(beta_values[variable],format = "e", digits = 2)
  confint_start =formatC( conf_int[variable,]["2.5 %"],format = "e", digits = 2)
  confint_end =formatC(conf_int[variable,]["97.5 %"], format = "e", digits = 2)
  p_value = round(model_summary$coefficients[variable, "Pr(>|t|)"],4)
  p_code =get_p_code(p_value)
  chisq = round(variable_power[variable,]$Chisq,3)
  MLR_LATEX =paste(MLR_LATEX, variable, " &  ", beta, " &  [", confint_start, " , ",
                   confint_end, "] &", chisq, "& ",
                   p_value, p_code, "\\ \\hline NEWLINE",  sep=" ")
  
}

print(MLR_LATEX)

print(get.r.squared(model_fit))

#log likelihood test
 lrtest(model_fit)
 
 
 #----------multinomial model-------------------------------
 
#Exclude others category
 excludeOtherCR =datasetCR[datasetCR$comment_group!='OTHER',]
 
 excludeOtherCR$CommentG <-factor(excludeOtherCR$comment_group)
 print(levels(excludeOtherCR$CommentG))
 

 #Relevel to make functional as reference
 excludeOtherCR$CommentG <- relevel(excludeOtherCR$CommentG, ref=4)
 print(levels(excludeOtherCR$CommentG))

 MLR_formula_string ="CommentG ~ IBF "
 
 for (variable in survived_vars){
   if (variable !="IBF")
     MLR_formula_string =paste(MLR_formula_string, "+ ", variable,  sep=" ")
   
 }
 
 print(MLR_formula_string)


 multinom_model <- multinom(as.formula(MLR_formula_string), data = excludeOtherCR, model = TRUE)
 summary(multinom_model)
 

 # example from https://stats.oarc.ucla.edu/r/dae/multinomial-logistic-regression/
 
 z <- summary(multinom_model)$coefficients/summary(multinom_model)$standard.errors
 p <- (1 - pnorm(abs(z), 0, 1)) * 2
 print(p)
 # Checking the model
 odds_ratios =exp(coef(multinom_model))
 print(odds_ratios)
 
 
 MLR_table=""
 
 for (variable in survived_vars){
   
   #beta =formatC(beta_values[variable],format = "e", digits = 2)
   
   MLR_table=paste(MLR_table,variable," & ",format(round(odds_ratios[1,variable], 2), nsmall = 2)," & ",format(round(p[1,variable], 4), nsmall = 2),
                   get_p_code(p[1,variable])," & ",format(round(odds_ratios[2,variable], 2), nsmall = 2),
                   " & ",format(round(p[2,variable], 4), nsmall = 2),get_p_code(p[2,variable]),
                   " & ",format(round(odds_ratios[3,variable], 2), nsmall = 2),
                   " & ",format(round(p[3,variable], 4), nsmall = 2),get_p_code(p[3,variable])," & ",
                   format(round(odds_ratios[4,variable], 2), nsmall = 2)," & ",
                   format(round(p[4,variable], 4), nsmall = 2),get_p_code(p[4,variable]),"\\ \\hline NEWLINE " , sep="")
   
 }
 
 print(MLR_table)
 
 PseudoR2(multinom_model, which = "Nagelkerke")
 lrtest(multinom_model)
 
