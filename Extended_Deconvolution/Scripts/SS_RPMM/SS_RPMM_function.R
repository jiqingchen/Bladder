################################################################################
################### GENERAL PACKAGES REQUIRED FOR SS-RPMM   ####################
###################          By: Devin C. Koestler          ####################
################################################################################
# load in the necessary packages.
library(RPMM)
library(nlme)
library(survival)
################################################################################
####################### GENERAL FUNCTIONS FOR SS-RPMM   ########################
################################################################################
################################################################################
# FUNCTION(S):  llikeHMObject, ebayes, and ClassAssign
# PURPOSE:      Used together to Predict methylation class for future subjects
#
################################################################################
llikeHMObject = function(o, x, type){
  J = dim(x)[2]
  L = rep(0, dim(x)[1])
  for(j in 1:J){
    if(type=="blc") {
      L = L + dbeta(x[,j], o$unsplit$a[j], o$unsplit$b[j], log=TRUE)
    }
    else if(type=="glc"){
      L = L + dnorm(x[,j], o$unsplit$mu[j], o$unsplit$sigma[j], log=TRUE)
    }
  }
  L }
ebayes = function(rpmm,x,type,nodelist=NULL){
  ApplyFun = get(paste(type,"TreeApply",sep=""))
  if(is.null(nodelist)){
    nodelist = unlist(ApplyFun(rpmm,
                               function(u,tree) u,
                               terminalOnly=TRUE, asObject=FALSE))
  }
  NN = length(nodelist)
  nx = dim(x)[1]
  L = matrix(NA, nx, NN)
  eta = rep(NA, NN)
  for(i in 1:NN){
    L[,i] = llikeHMObject(rpmm[[nodelist[i]]],x,type)
    eta[i] = sum(rpmm[[nodelist[i]]]$weight)
  }
  eta = eta/sum(eta)
  L[is.na(L)] = min(L,na.rm=TRUE)
  for(i in 1:nx){
    L[i,] = exp(L[i,] - max(L[i,])) * eta
    L[i,] = L[i,]/sum(L[i,])
  }
  
  L }
ClassAssign = function(x, d) {
  classval = NULL;
  for(i in 1:dim(x)[1]) {
    t0 = 0
    tmp1 = NULL
    for(j in 1:dim(x)[2]) {
      if (x[i,j] >= t0) {
        classval[i] = d[j]
        tmp1 = classval[i]
        t0 = x[i,j] }
      else classval[i] = tmp1
    }
  }
  return(classval)
}
################################################################################
# FUNCTION:  TrainTestSplit
#    Stratified random split of the data into training and testing sets
#
# ARGUMENTS:
# Y:
#    Covariates:
# Strat: #
# seed:
#    propTrain:
#
# RETURNS:   A list whose first item represents the training data and second
#            item represents the testing data
#
################################################################################
TrainTestSplit = function(Y, Covariates, Strat = NULL, seed = NULL, propTrain = 1/2) {
  fullData = data.frame(Covariates, Y)
  strat = with(fullData, paste(fullData[,Strat]))
  sflag = T
  nMan = dim(fullData)[1]
  stratL = split((1:nMan)[sflag], strat[sflag])
  nStrat = length(stratL)
  fullData$TrainingSet = rep(-1, nMan)
  if(is.null(seed)){
    SEED = abs(round(rnorm(1, mean = 1000, sd = 500)))
  }
  else SEED = seed
  set.seed(SEED)
  xntr = 0
  for(s in 1:nStrat){
    
    ixs = stratL[[s]]
    nixs = length(ixs)
    if(nixs %% 2 == 1){
      ntr = propTrain*(nixs-1) + xntr
      xntr = 1-xntr
    }
    else ntr = propTrain*nixs
    if(ntr>0){
      #hack for the n=1 problem
      if(nixs==1) fullData$TrainingSet[ixs] = 1
      else{
        trs = sample(ixs, ntr)
        fullData$TrainingSet[trs] = 1
      }
    } }
  TrainingData = subset(fullData, TrainingSet == 1)[,1:(dim(fullData)[2]-1)]
  TestingData = subset(fullData, TrainingSet == -1)[,1:(dim(fullData)[2]-1)]
  list(TrainingData, TestingData)
}
################################################################################
# FUNCTION:  DummyVars
#    Function that generates a matrix of dummy variables from a factor
#
# ARGUMENTS:
#    dat:    data.frame, matrix, or vector of factor(s)
#
# RETURNS:  Matrix n x (L-1) where L is the number of levels for the factor(s)
#
################################################################################
DummyVars = function(dat) {
  if(dim(as.matrix(dat))[2] == 1) {
    dat.fact = as.factor(dat)
    mat = matrix(nrow = length(dat), ncol = length(levels(dat.fact))-1)
    for( i in 1:(length(levels(dat.fact))-1)) {
      mat[,i] = ifelse(dat.fact == levels(dat.fact)[i], 1, 0)
    }
  } else {
    numLevels = NumLevels(dat) - dim(dat)[2]
    mat = matrix(nrow = dim(dat)[1], ncol = numLevels)
    index = 1
    for(j in 1:dim(dat)[2]) {
      dat.fact = as.factor(dat[,j])
      for(i in 1:(length(levels(dat.fact))-1)) {
        mat[,index] = ifelse(dat.fact == levels(dat.fact)[i], 1, 0)
        index = index+1
      }
    } }
  mat
}
################################################################################
# FUNCTION:  NumLevels
#    Function that computes the number of total levels for factor(s)
#
# ARGUMENTS:
#    dat:    data.frame, matrix, or vector of factor(s)
#
# RETURNS:  Numeric value indicating the number of factor levels
#
################################################################################
NumLevels = function(dat) {
  numlevel = 0
  dat1 = data.frame(dat)
  for(i in 1:dim(as.matrix(dat))[2]) {
    numlevel = numlevel + length(levels(as.factor(dat1[,i])))
  }
  numlevel }
################################################################################
# FUNCTION:  MostImpCpGsSurvival
#    Function that computes the raw Cox-scores for each of the CpG loci
#
# ARGUMENTS:
# Y:
#    covariates:
#
# times:
# censor:
# terms: #
#
# factors: #
#    strat:
#
#
#
# RETURNS:
#
################################################################################
MostImpCpGsSurvival = function(Y, covariates, times, censor, terms = NULL,
                               factors = NULL, strat = NULL)  {
  terms1 = terms[!(terms %in% c(factors, strat))]
  J = dim(Y)[2]
  coxscores = NULL
  if(is.null(terms)& is.null(factors) & is.null(strat)) {
    for(i in 1:J) {
      surv = coxph(Surv(covariates[,times], covariates[,censor]) ~ Y[,i])
      coxscores[i] = summary(surv)[[7]][[4]]
    } }
  
  else if (!is.null(terms) & is.null(factors) & is.null(strat)) {
    for(i in 1:J) {
      surv = coxph(Surv(covariates[,times], covariates[,censor]) ~ Y[,i] +
                     as.matrix(covariates[,as.character(terms1)]))
      coxscores[i] = summary(surv)[[7]][[((length(terms1)+1)*3 + 1)]]
    } }
  else if (!is.null(terms) & !is.null(factors) & is.null(strat)) {
    factMat = DummyVars(covariates[,factors])
    for(i in 1:J) {
      surv = coxph(Surv(covariates[,times], covariates[,censor]) ~ Y[,i] +
                     as.matrix(cbind(covariates[,as.character(terms1)],factMat)))
      coxscores[i] = summary(surv)[[7]][[(length(terms1)+
                                            NumLevels(covariates[,factors])-length(factors)+1)*3+1]]
    } }
  else if (!is.null(terms) & !is.null(factors) & !is.null(strat)) {
    factMat = DummyVars(covariates[,factors])
    for( i in 1:J) {
      surv = coxph(Surv(covariates[,times], covariates[,censor]) ~ Y[,i] +
                     as.matrix(cbind(covariates[,as.character(terms1)],factMat))+
                     strata(covariates[,as.character(strat)]))
      coxscores[i] = summary(surv)[[7]][[(length(terms1)+
                                            NumLevels(covariates[,factors])-length(factors)+1)*3+1]]
    } }
  else if (!is.null(terms) & is.null(factors) & !is.null(strat)) {
    for( i in 1:J) {
      surv = coxph(Surv(covariates[,times], covariates[,censor]) ~ Y[,i] +
                     as.matrix(covariates[,as.character(terms1)])+
                     strata(covariates[,as.character(strat)]))
      coxscores[i] = summary(surv)[[7]][[(length(terms1)+1)*3+1]]
    } }
  abscoxscores = abs(coxscores)
  CpGnames = colnames(Y)
  AbsCox = data.frame(CpGnames, abscoxscores)
  CoxOrd = AbsCox[order(abscoxscores, decreasing = T),]
  rownames(CoxOrd) = CoxOrd[,1]
  CoxOrd1 = data.frame(CoxOrd[,-1])
  rownames(CoxOrd1) = CoxOrd[,1]
  colnames(CoxOrd1) = "AbsCoxScore"
  CoxOrd1
}
################################################################################
# FUNCTION:  MostImpCpGs
#    Function that computes the raw Tscores for each of the CpG loci
#
# ARGUMENTS:
# 
#   Y: Data frame/matrix (n x J) of beta values for dataset.
#   covariates: Data frame (n x P) of covariates for dataset, including
#               the clinical variable of interest and random effect term
#   clinvar: Clinical variable of interest (i.e. case status)
#   terms: Terms used in the model given as a character vector.  Leave
#          this blank if you don't want to control for any other covariates.
#   factors: Terms in the model to be treated as factors (i.e. used for 
#            categorical covariates).  Given as a character vector.
#   randomEffect: Name of term to be included as a random effect
#   is.factor: is the clinical variable of interest categorial - TRUE if yes.
#
# RETURNS:   The raw TScores for each of the J CpG loci
#
################################################################################

MostImpCpGs = function(Y, covariates, clinvar, terms = NULL, factors = NULL,
                       randomEffect = NULL, is.factor = FALSE)  {
  terms1 = terms[!(terms %in% c(factors, clinvar, randomEffect))]
  J = dim(Y)[2]
  TScores = NULL
  if(is.null(randomEffect)) {
    if(is.null(terms1) & is.null(factors)) {
      for(i in 1:J) {
        asinSqrt = Y[,i]
        reg = lm(asinSqrt ~ covariates[,clinvar])
        TScores[i] = summary(reg)[[4]][[6]]
      } }
    else if (!is.null(terms1) & is.null(factors)) {
      for( i in 1:J) {
        asinSqrt = Y[,i]
        reg = lm(asinSqrt ~ covariates[,clinvar] + as.matrix(covariates[,terms1]))
        TScores[i] = summary(reg)[[4]][[(length(terms1)+2)*3-length(terms1)]]
      } }
    else if (!is.null(terms1) & !is.null(factors)) {
      factMat = DummyVars(covariates[,factors])
      for( i in 1:J) {
        asinSqrt = Y[,i]
        reg = lm(asinSqrt ~ covariates[,clinvar] +
                   as.matrix(cbind(covariates[,terms1],factMat)))
        TScores[i] = summary(reg)[[4]][[(length(terms1)+NumLevels(covariates[,factors])-
                                           length(factors)+2)*3 -
                                          (length(terms1)+NumLevels(covariates[,factors])-
                                             length(factors))]]
      } }
    else if (is.null(terms1) & !is.null(factors)) {
      factMat = DummyVars(covariates[,factors])
      for( i in 1:J) {
        asinSqrt = Y[,i]
        reg = lm(asinSqrt ~ covariates[,clinvar] + as.matrix(factMat))
        TScores[i] = summary(reg)[[4]][[(NumLevels(covariates[,factors]) -
                                           length(factors)+2)*3 - (NumLevels(covariates[,factors]) -
                                                                     length(factors))]]
      } }
  }
  else if (!is.null(randomEffect)) {
    if(is.null(terms1) & is.null(factors)) {
      cv = as.matrix(covariates[,clinvar])
      re = covariates[,randomEffect]
      for(i in 1:J) {
        asinSqrt = Y[,i]
        reg = try(lme(asinSqrt ~cv,random = ~1|as.factor(re),
                      na.action = na.omit, control = lmeControl(opt = "optim")), silent = T)
        if(is.list(summary(reg)) == T) TScores[i] = summary(reg)[[19]][[8]]
        else TScores[i] = NA
      }
    }
    else if (!is.null(terms1) & is.null(factors)) {
      cv = as.matrix(covariates[,clinvar])
      re = covariates[,randomEffect]
      ot = as.matrix(covariates[,terms1])
      for( i in 1:J) {
        asinSqrt = Y[,i]
        reg = try(lme(asinSqrt ~cv+ot,random = ~1|as.factor(re),
                      na.action = na.omit, control = lmeControl(opt = "optim")), silent = T)
        if(is.list(summary(reg)) == T) {
          TScores[i] = summary(reg)[[19]][[(length(terms1)+2)*4-(length(terms1))]]}
        else TScores[i] = NA
      } }
    else if (!is.null(terms1) & !is.null(factors)) {
      factMat = DummyVars(covariates[,factors])
      cv = as.matrix(covariates[,clinvar])
      re = covariates[,randomEffect]
      ot = as.matrix(cbind(covariates[,terms1],factMat))
      numVar = dim(ot)[2]
      for( i in 1:J) {
        asinSqrt = Y[,i]
        reg = try(lme(asinSqrt ~cv+ot,random = ~1|as.factor(re),
                      na.action = na.omit, control = lmeControl(opt = "optim")), silent = T)
        if(is.list(summary(reg)) == T) {
          TScores[i] = summary(reg)[[19]][[(numVar+2)*4-(numVar)]]}
        else TScores[i] = NA
      } }
    else if (is.null(terms1) & !is.null(factors)) {
      factMat = DummyVars(covariates[,factors])
      cv = as.matrix(covariates[,clinvar])
      re = covariates[,randomEffect]
      ot = as.matrix(factMat)
      numVar = dim(ot)[2]
      for( i in 1:J) {
        asinSqrt = Y[,i]
        reg = try(lme(asinSqrt ~cv+ot,random = ~1|as.factor(re),
                      na.action = na.omit, control = lmeControl(opt = "optim")), silent = T)
        
        if(is.list(summary(reg)) == T) {
          TScores[i] = summary(reg)[[19]][[(numVar+2)*4-(numVar)]]}
        else TScores[i] = NA
      } }
  }
  absTScores = abs(TScores)
  CpGnames = colnames(Y)
  AbsT = data.frame(CpGnames, absTScores)
  TOrd = AbsT[order(absTScores, decreasing = T),]
  rownames(TOrd) = TOrd[,1]
  TOrd1 = data.frame(TOrd[,-1])
  rownames(TOrd1) = TOrd[,1]
  colnames(TOrd1) = "AbsTScore"
  TOrd1
}
################################################################################
# FUNCTION:  NestedXValidation
#    Function that performs a Nested X-Validation to determine the number of
#    CpG loci to be used in fitting RPMM
#
# ARGUMENTS:
#     Y: Data frame/matrix (n x J) of beta values for training dataset.
#    covariates: Data frame of the covariates for the training dataset
#    TScores: TScores object from MostImpGpGs function
#    clinvar: Clinical variable of interest (i.e. case status). Supplied as a vector
#    vartype: Variable type of the clinical variable of interest.  Options
#             are: "binary" and "continuous"
#    mrange: Minimum and maximum number of CpGs to be used in fitting
#            RPMM.  Supplied as vector (ie. c(5, 50), where 5 is the
#    method: Fits a guassian- or beta-distributed RPMM. Defaults
#            to a gaussian-distributed RPMM minimum number of CpGs considered and 50 is the max).
#            Defaults to c(5,50).
#    L: Number of nested splits of the Training data. Defaults to 10 splits
#    seeds: 
#
#
#
# RETURNS:
#
################################################################################
NestedXValidation = function(Y, covariates, TScores, clinvar, vartype, 
                             method = "gaussian", mrange = c(5,50), 
                             L = 10,seeds = NULL) {
M = mrange[2]
m = mrange[1]
if(is.null(seeds)) SEEDS = round(abs(rnorm(L, mean = 100, sd = 25)),0)
else SEEDS = seeds

k = round(dim(Y)[1]/2)
index = 1:dim(Y)[1]
p = dim(covariates)[2]
Results = matrix(nrow = (M-m+1), ncol = L)
for(j in 1:L) {
  print(paste("L = ", j, sep = ""))
  splitDat = TrainTestSplit(Y[,as.character(rownames(TScores)[1:M])],
                            covariates, Strat = clinvar, seed = SEEDS[j])
  Learning = splitDat[[1]]
  Validation = splitDat[[2]]
  LearningData = Learning[-(1:p)]
  LearningDataPheno = Learning[,clinvar]
  ValidationData = Validation[-(1:p)]
  ValidationDataPheno = Validation[,clinvar]
  
  if(method == "gaussian") {
    for(i in 1:(M-m+1)) {
        print(paste("M = ", i+(m-1), sep = ""))
        ClassesValidation = NULL
    RPMMobject = try(glcTree(as.matrix(LearningData[,as.character(rownames(TScores)[1:(i+(m-1))])]),
                                 verbose = 0),silent = FALSE)
    ClassesRPMM = levels(glcTreeLeafClasses(RPMMobject))
    if(length(ClassesRPMM) > 1 ) {
      PosteriorValidation = ebayes(RPMMobject,
                             as.matrix(ValidationData[,as.character(rownames(TScores)[1:(i+(m-1))])]),,"glc")
                             ClassesValidation = ClassAssign(PosteriorValidation, ClassesRPMM)
                             if(length(levels(as.factor(ClassesValidation))) > 1) {
                               if(vartype == "binary") {
                                 Results[i,j] = summary(table(ClassesValidation,
                                                              ValidationDataPheno))$p.value}
                               else if (vartype == "continuous") {
                                 Results[i,j] =
                                   kruskal.test(ValidationDataPheno,factor(ClassesValidation))$p.value}
                             }
                             else Results[i,j] = NA
        }
        else  Results[i,j] = NA
      } }
    else if(method == "beta") {
      for(i in 1:(M-m+1)) {
        print(paste("M = ", i+(m-1), sep = ""))
        ClassesValidation = NULL
        RPMMobject = try(blcTree(as.matrix(LearningData[,as.character(rownames(TScores)[1:(i+(m-1))])]),
                                 verbose = 0),silent = FALSE)
        ClassesRPMM = levels(blcTreeLeafClasses(RPMMobject))
        if(length(ClassesRPMM) > 1 ) {
          PosteriorValidation = ebayes(RPMMobject,
                                       as.matrix(ValidationData[,as.character(rownames(TScores)[1:(i+(m-1))])]),
                                       "blc")
          ClassesValidation = ClassAssign(PosteriorValidation, ClassesRPMM)
          if(length(levels(as.factor(ClassesValidation))) > 1) {
            if(vartype == "binary") {
              Results[i,j] = summary(table(ClassesValidation,
                                           ValidationDataPheno))$p.value }
            else if (vartype == "continuous") {
              Results[i,j] =
                kruskal.test(ValidationDataPheno,factor(ClassesValidation))$p.value}
          }
          else Results[i,j] = NA
        }
        else  Results[i,j] = NA
      } }
  }
  rownames(Results) = paste("M = ", m:M)
  Results
}

################################################################################
# FUNCTION:  NestedXValidationSurvival
#    Function that performs a Nested X-Validation to determine the number of
#    CpG loci to be used in fitting RPMM.  To be used when the clinical variable
#    of interest is survival
#
# ARGUMENTS:
# Y:
#    covariates:
# CoxScores:
# times:
# censor: #
# mrange: #
#
#
# method: #
# L: #
#    seeds:
#
#
#
# RETURNS:
#
################################################################################

NestedXValidationSurvival = function(Y, covariates, CoxScores, times, censor, method =
                                       "gaussian", mrange = c(5,50), L = 10,
                                     seeds = NULL) {
  M = mrange[2]
  m = mrange[1]
  if(is.null(seeds)) SEEDS = round(abs(rnorm(L, mean = 100, sd = 25)),0)
  else SEEDS = seeds
  survDat = c(times, censor)
  k = round(dim(Y)[1]/2)
  index = 1:dim(Y)[1]
  p = dim(covariates)[2]
  Results = matrix(nrow = (M-m+1), ncol = L)
  for(j in 1:L) {
    print(paste("L = ", j, sep = ""))
    splitDat = TrainTestSplit(Y[,as.character(rownames(CoxScores)[1:M])],
                              covariates, Strat = censor, seed = SEEDS[j])
    Learning = splitDat[[1]]
    Validation = splitDat[[2]]
    LearningData = Learning[-(1:p)]
    LearningDataPheno = Learning[,survDat]
    ValidationData = Validation[-(1:p)]
    ValidationDataPheno = Validation[,survDat]
    if(method == "gaussian") {
      for(i in 1:(M-m+1)) {
        print(paste("M = ", i+(m-1), sep = ""))
        ClassesValidation = NULL
        RPMMobject = try(glcTree(as.matrix(LearningData[,as.character(rownames(CoxScores)[1:(i
                                                                                             +(m-1))])]),
                                 verbose = 0),silent = FALSE)
        ClassesRPMM = levels(glcTreeLeafClasses(RPMMobject))
        if(length(ClassesRPMM) > 1 ) {
          PosteriorValidation = ebayes(RPMMobject,
                                       as.matrix(ValidationData[,as.character(rownames(CoxScores)[1:(i+(m-1))])]), "glc")
          ClassesValidation = ClassAssign(PosteriorValidation, ClassesRPMM)
          if(length(levels(as.factor(ClassesValidation))) > 1) {
            lr.chisq = survdiff(Surv(ValidationDataPheno[,1],ValidationDataPheno[,2])~
                                  as.factor(ClassesValidation))[[5]]
            Results[i,j] = pchisq(lr.chisq, df =
                                    (length(levels(as.factor(ClassesValidation)))-1),
                                  lower.tail = F)
          }
          else Results[i,j] = NA
        }
        else  Results[i,j] = NA
      }
    }
    else if(method == "beta") {
      for(i in 1:(M-m+1)) {
        print(paste("M = ", i+(m-1), sep = ""))
        ClassesValidation = NULL
        RPMMobject = try(blcTree(as.matrix(LearningData[,as.character(rownames(CoxScores)[1:(i
                                                                                             +(m-1))])]),
                                 verbose = 0),silent = FALSE)
        ClassesRPMM = levels(blcTreeLeafClasses(RPMMobject))
        if(length(ClassesRPMM) > 1 ) {
          PosteriorValidation = ebayes(RPMMobject,
                                       as.matrix(ValidationData[,as.character(rownames(CoxScores)[1:(i+(m-1))])]), "blc")
          ClassesValidation = ClassAssign(PosteriorValidation, ClassesRPMM)
          if(length(levels(as.factor(ClassesValidation))) > 1) {
            lr.chisq = survdiff(Surv(ValidationDataPheno[,1],ValidationDataPheno[,2])~
                                  as.factor(ClassesValidation))[[5]]
            Results[i,j] = pchisq(lr.chisq, df =
                                    (length(levels(as.factor(ClassesValidation)))-1),
                                  lower.tail = F)
          }
          else Results[i,j] = NA
        }
        else  Results[i,j] = NA
      }
    } }
  rownames(Results) = paste("M = ", m:M)
  apply(Results, 1, median, na.rm = T)
}

################################################################################
# FUNCTION:  PredMethClasses
#    Function that fits an RPMM to a Training data and predicts methylation
#    class membership for observations in a Testing data
#
# ARGUMENTS:
# Ytrain:
# Ytest
# Scores: #
# M: #
#    method:
#
#
# RETURNS:
#
################################################################################

PredMethClasses = function(Ytrain, Ytest, Scores, M, method = "gaussian") {
  if(method == "gaussian") {
    RPMMTrain = glcTree(as.matrix(Ytrain[,as.character(rownames(Scores)[1:M])]),
                        verbose = 0)
    RPMMTrainClasses = levels(glcTreeLeafClasses(RPMMTrain))
    PosteriorTesting = ebayes(RPMMTrain,
                              as.matrix(Ytest[,as.character(rownames(Scores)[1:M])]), "glc")
    TestingClasses = ClassAssign(PosteriorTesting, RPMMTrainClasses)
    as.factor(TestingClasses)
  }
  else if(method == "beta") {
    RPMMTrain = blcTree(as.matrix(Ytrain[,as.character(rownames(Scores)[1:M])]),
                        verbose = 0)
    RPMMTrainClasses = levels(blcTreeLeafClasses(RPMMTrain))
    PosteriorTesting = ebayes(RPMMTrain,
                              as.matrix(Ytest[,as.character(rownames(Scores)[1:M])]), "blc")
    TestingClasses = ClassAssign(PosteriorTesting, RPMMTrainClasses)
    as.factor(TestingClasses)
  }
}

################################################################################
################################################################################
#################    FUNCTIONS FOR POST SS-RPMM ANALYSES    ####################
################################################################################
################################################################################

permTestChiSquare <- function(x,y,R=10000,report=1000){
  obs <- !is.na(x) & !is.na(y)
  x <- x[obs]
  y <- y[obs]
  stat0 <- summary(table(x, y))$stat
  stats <- rep(NA,R)
  n <- length(y)
  for(r in 1:R){
    if(r %% report==0) cat(r,"\n")
    stats[r] <- summary(table(x, y[sample(1:n,n,replace=FALSE)]))$stat
  }
  mean(stats>stat0)
}
permTestKruskal <- function(x,y,R=10000,report=1000){
  obs <- !is.na(x) & !is.na(y)
  x <- x[obs]
  y <- y[obs]
  stat0 <- kruskal.test(y~x)$stat
  stats <- rep(NA,R)
  n <- length(y)
  for(r in 1:R){
    if(r %% report==0) cat(r,"\n")
    yy <- y[sample(1:n,n,replace=FALSE)]
    stats[r] <- kruskal.test(yy~x)$stat
  }
  mean(stats>stat0)
}
plotMethByClass = function(BETA, CLASS, sep="red", ordSubj=NULL,
                           col=colorRampPalette(c("yellow","black","blue"),space="Lab")(64)){
  maxBETA = max(BETA, na.rm = T)
  minBETA = min(BETA, na.rm = T)
  nColor = length(col)
  nCpG = dim(BETA)[2]
  nSubj = dim(BETA)[1]
  index = split(1:nSubj, CLASS)
  nClass = length(index)
  plot(c(0,nCpG), c(0,nSubj), type="n", xlab="", ylab="",
       xaxt="n", yaxt="n", bty="n")
  if(is.null(ordSubj)){
    ordSubj = hclust(dist(t(BETA)), method="ward")$ord
  }
  BETA = BETA[,ordSubj]
  k=0
  ordRet = NULL
  for(i in 1:nClass){
    ii = index[[i]]
    nii = length(ii)
    ord = hclust(dist(BETA[ii,]), method="ward")$ord
    for(j in 1:nii){
      colori = ceiling((BETA[ii[ord[j]],] - minBETA)/(maxBETA - minBETA)*nColor)
      rect(1:nCpG-1,k,1:nCpG,k+1,col=col[as.numeric(colori)],density=-1,border=NA)
      k = k+1
    }
    if(i<nClass) lines(c(0,nCpG), c(k,k), col=sep, lwd = 4)
    ordRet = c(ordRet,ii[ord])
  }
  nn = cumsum(c(0,unlist(lapply(index,length))))
  axp = 0.5*(nn[-1]+nn[-nClass-1])
  axis(2, axp, names(index), las=2, cex.axis = 1.5)
  invisible(ordRet)
}
