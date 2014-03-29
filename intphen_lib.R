# Library of scripts for Phenotypic Integration Analysis


###############################################################
# To cite these functions, use: ###############################
# Torices R & Méndez M. Resource allocation to inflorescence components is highly integrated despite differences between allocation currencies and sites. International Journal of Plant Sciences (In press)

###############################################################
# INDEX #######################################################
# The functions included in this script are:

# 1. intphen ###################################################
# It allows to estimate:
#    INT: The variance of the eigenvalues (Wagner 1984; Cheverud et al. 1989; Pavlicev et al 2009)
#    INTc: A correction of INT based on the number of traits (p) and the number of individuals (n): INTc= INT - (p-1/n)
#    RelINT: The Relative Phenotypic Integration index: INT index by the maximum possible integration value for a given number of traits (INT/p-1). RelINT can be estimated
#    for INT or for INTc

# 2. intphen.boot ##############################################
# This function estimates the confidence intervals of INTc (the corrected INT) using randomizations

# 3. cor.par ###################################################
# This function calculates the partial correlation matrix between a set of traits and a third control variable. It is utilised internally
# by 'intphen.par' and 'intphen.par.boot' to estimate the phenotypic integration index from partial correlations.

# 4. intphen.par ############################################
# It allows to estimate the same set of indices than intphen but from partial correlation matrix instead of the correlation matrix:
#    INTpc: The variance of the eigenvalues using the partial correlation matrix instead of a correlation matrix. We proposed this
#    index (Torices & Méndez 2014) for analysing phenotypic integration of resource allocation to different components in an individual
#    or a given organ when resource allocation data come from observational studies in which resource availability is not controlled and
#    therefore could lead to components correlations only by the fact of different resource availability. Thus, INT can be estimated using
#    the matrix of partial correlations in which size of the organ or individual is used as the third control variable.
#    INTpc.c: The INTpc corrected by the number of traits and individuals of each population
#    RelINTpc: The same as RelINT but estimated with the partial correlation matrix

# 5. intphen.par.boot ##############################################
# This function estimates the confidence intervals of INTc (the corrected INT) using randomizations

# 6. IIsdeCor ######################################################
# This function estimates the standard deviation of eigenvalues of a correlation matrix (Based on Haber (2011)). 
# It is a small modification of that script provided by Haber (2011) to estimate this index using correlation matrices



########################################################################
# INT, INTc, RelINT, RelINTc ----> intphen()
# To use: intphen(yourtraits)
# You must provide a dataframe or a matrix object with all your traits in which columns are the traits and rows are the individuals
# IMPORTANT: Do not include missing values in your data set. Otherwise corrected indices are wrong!!!!!!
# At the moment the functions return the next values:
# INT (Phenotypic integration index) =   
# RelINT (% of maximum possible integration) = 
# INT.c (Corrected phenotypic integration index) = 
# RelINT.c (% of maximum possible integration) = 

intphen<-function(traits){
  X<-na.exclude(traits)
  cor_X <-cor(X)
  eig_X<-eigen(cor_X, only.values=TRUE)$values
  d <- eig_X
  p <- length (d)
  n <- nrow(X)
  INT<-sum((d-1)^2)/(p-1)
  INT.c<-(INT-((p-1)/n))
  pref="INT (Phenotypic integration index) =  "
  print(paste(pref, round(INT, 3)))
  pref2="RelINT (% of maximum possible integration) = "
  perc<-(INT/(p-1))*100
  print(paste(pref2, round(perc, 3)))
  pref3="INT.c (Corrected phenotypic integration index) = "
  print(paste(pref3, round(INT.c, 3)))
  pref4="RelINT.c (% of maximum possible integration) = "
  perc.c<-(INT.c/(p-1))*100
  print(paste(pref4, round(perc.c, 3)))
  pref5="Number of observations used = "
  print(paste(pref5, n))
}



########################################################################
# Randomization test of INT significance ----> intphen.boot()
# To use: intphen.boot (yourtraits, the_number_of_replicates)
# 'yourtraits' must be a dataframe or a matrix object with all your traits in which columns are the traits and rows are the individuals
# 'the_number_of_replicates' indicates the number of randomizations that will be performed (1000 by default)
# IMPORTANT: Do not include missing values in your data set. Otherwise corrected indices are wrong!!!!!!


intphen.boot<-function (X,Y=1000){
  X<-na.exclude(X)
  INT = list()
  length (INT) = Y
  for (i in 1:Y){
    cor_X<-cor(X[sample(nrow(X), replace=TRUE),])
    d<-eigen(cor_X, only.values=TRUE)$values
    p <- length (d)
    n <- nrow(X)
    Int<-sum((d-1)^2)/(p-1)
    Int.c<-(Int-((p-1)/n))
    INT[i]=Int.c
  }
  Intphen <-as.numeric(INT)
  pref="Mean = "
  print(paste(pref, round(mean(Intphen), 3)))
  pref="Median ="
  print(paste(pref, round(median(Intphen), 3)))
  pref2="SD = "
  print(paste(pref2, round(sd(Intphen), 3)))
  pref3="SE = "
  se<-(sd(Intphen)/sqrt(nrow(X)))
  print (paste(pref3, round(se, 3)))
  pref4="Lower IC 99% = "
  print(paste(pref4,round(quantile(Intphen, probs=0.5/100), 3)))
  pref5="Higher IC 99% = " 
  print(paste(pref5,round(quantile(Intphen, probs=99.5/100), 3)))  
  pref6="Lower IC 95% = "
  print(paste(pref6,round(quantile(Intphen, probs=2.5/100), 3)))
  pref7="Higher IC 95% = " 
  print(paste(pref7,round(quantile(Intphen, probs=97.5/100), 3)))  
  pref8="Number of replicates = "
  print (paste(pref8, length(INT)))
}



########################################################################
# Partial correlation matrix -----> cor.par()
# To use: cor.par (traits, c.trait, trait.names=FALSE)
# 'traits': a data frame or a matrix with your traits. (the same format than 'intphen')
# 'c.trait': the third control variable to estimate partial correlations, for instance: the size of the organ or the individual
# 'trait.names'(logical, default=FALSE): 'TRUE' to get the names of the traits in the output; 'FALSE': correlation matrix without row and column names.
# IMPORTANT: Do not include missing values in your data set. 
# Dependence: library(ppcor)

library(ppcor)

cor.par<-function(traits, z, trait.names=FALSE) {
  ntraits<-ncol(traits)
  r <- matrix(0, nrow = ntraits, ncol = ntraits)
  for (i in seq_len(ntraits)) {
    for (j in seq_len(i)) {
      if (i==j) r[i, j] <- 1
      else {
        x2 <- traits[, i]
        y2 <- traits[, j]
        z2 <- z
        r[i, j] <-pcor.test(x2, y2, z2, method="pearson")$estimate
      }}}
  r <- r + t(r) - diag(diag(r))
  if (trait.names==FALSE) return(r)
  else{
    rownames(r) <- colnames(traits)
    colnames(r) <- colnames(traits)
    r
  }}  



#########################################################################
# INTpc, INTpc.c, RelINTpc, RelINTpc.c  ------->  intphen.par()
# To use: intphen.par (yourtraits)
# 'yourtraits' must be a dataframe or a matrix object with all your traits in which columns are the traits and rows are the individuals.
#             The third control variable (for instance the size) must always be the last column in the dataframe
# At the moment the functions return the next values:
# INTpc (Phenotypic integration index of partial correlations) =   
# RelINTpc (% of maximum possible integration of partial correlations) = 
# INTpc.c (Corrected phenotypic integration index of partial correlations) = 
# RelINTpc.c (% of maximum possible integration of partial correlations) = 
# Dependence: libary (ppcor), and the function 'cor.par'
# This function calls internally the function 'cor.par'. Then you need to load it first.


intphen.par<-function(X){
  X<-na.exclude(X)
  N<-ncol(X)
  traits<-X[,1:(N-1)]
  c.trait<-X[,N]
  cor_X<-cor.par(traits, c.trait)
  eig_X<-eigen(cor_X, only.values=TRUE)$values
  d <- eig_X
  p <- length (d)
  n <- nrow(X)
  INT<-sum((d-1)^2)/(p-1)
  INT.c<-(INT-((p-1)/n))
  pref="INTpc (Phenotypic integration index of partial correlations) =  "
  print(paste(pref, round(INT, 3)))
  pref2="RelINTpc (% of maximum possible integration of partial correlations) = "
  perc<-(INT/(p-1))*100
  print(paste(pref2, round(perc, 3)))
  pref3="INTpc.c (Corrected phenotypic integration index of partial correlations) = "
  print(paste(pref3, round(INT.c, 3)))
  pref4="RelINTpc.c (% of maximum possible integration of partial correlations) = "
  perc.c<-(INT.c/(p-1))*100
  print(paste(pref4, round(perc.c, 3)))
  pref5="Number of observations used = "
  print(paste(pref5, n))
}



##########################################################################
# Randomization test of INTpc significance ----> intphen.boot.par()
# To use: intphen.boot.par(yourtraits, the_number_of_replicates)
# 'yourtraits' must be a dataframe or a matrix object with all your traits in which columns are the traits and rows are the individuals.
#             The third control variable (for instance the size) must always be the last column in the dataframe
# 'the_number_of_replicates' indicates the number of randomizations that will be performed (1000 by default)
# IMPORTANT: Do not include missing values in your data set. Otherwise corrected indices are wrong!!!!!!
# Dependence: libary (ppcor), and the function 'cor.par'
# This function calls internally the function 'cor.par'. Then you need to load it first.

library(ppcor)

intphen.boot.par<-function (X,Y=1000){
  INT = list()
  length (INT) = Y
  X<-na.exclude(X)
  for (i in 1:Y){
    t.sample<-X[sample(nrow(X), replace=TRUE),]
    Z<-ncol(t.sample)
    traits<-t.sample[,1:(Z-1)]
    c.trait<-t.sample[,Z]
    cor_X<-cor.par(traits, c.trait)
    d <-eigen(cor_X, only.values=TRUE)$values
    p <- length (d)
    n <- nrow(t.sample)
    Int<-sum((d-1)^2)/(p-1)
    Int.c<-(Int-((p-1)/n))
    INT[i]=Int.c
  }
  Intphen <-as.numeric(INT)
  pref="Mean = "
  print(paste(pref, round(mean(Intphen), 3)))
  pref="Median ="
  print(paste(pref, round(median(Intphen), 3)))
  pref2="SD = "
  print(paste(pref2, round(sd(Intphen), 3)))
  pref3="SE = "
  se<-(sd(Intphen)/sqrt(nrow(X)))
  print (paste(pref3, round(se, 3)))
  pref4="Lower IC 99% = "
  print(paste(pref4,round(quantile(Intphen, probs=0.5/100), 3)))
  pref5="Higher IC 99% = " 
  print(paste(pref5,round(quantile(Intphen, probs=99.5/100), 3)))  
  pref6="Lower IC 95% = "
  print(paste(pref6,round(quantile(Intphen, probs=2.5/100), 3)))
  pref7="Higher IC 95% = " 
  print(paste(pref7,round(quantile(Intphen, probs=97.5/100), 3)))  
  pref8="Number of replicates = "
  print (paste(pref8, length(INT)))
}



############################################################################
# IIsdeCor ------> IIsdeCor ()
# To use: IIsdeCor (the_matrix_of_correlation)
# You can use either the matrix of correlation or the matrix of partial correlations as is calculated by 'par.cor'
# A modification of 'IIsde' From Haber (2011)

IIsdeCor <- function(X) {
  d <- eigen(X)$values
  p <- length(d)
  sqrt(sum((d-1)^2)/(p-1))
}
