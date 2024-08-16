#Util functions
#Univariate methods

#'This function computes the pairwise spearman correlation between two data.frames
#'@param X a data.frame
#'@param Y a data.frame
#'@return a list of lists of cor.test objects
#'
compute.pairwise.spearman = function(X,Y){
  
  ncol.X = ncol(X)
  ncol.Y = ncol(Y)
  
  lapply(seq_len(ncol.X), function(x){
    lapply(seq_len(ncol.Y), function(y){
      cor.test(X[,x], Y[,y], method="spearman")
    })
  })
  
}

#This function assesses the global association, while adjusting for multiplicity (ACAT)
#'@param cor.test.by.element a list returned by compute.pairwise.spearman
#'@param X a data.frame
#'@param Y a data.frame
#'@return a list of ACAT-combined p-values

ACAT.by.element = function(cor.test.by.element, X, Y){
  
  ncol.X = ncol(X)
  ncol.Y = ncol(Y)
  
  sapply(seq_len(ncol.X), function(x){
    ACAT::ACAT(sapply(seq_len(ncol.Y), function(y){
      cor.test.by.element[[x]][[y]]$p.value
    }))
  })
}

#This function computes the Jaccard index 
#' @param p.values a vector of p-values
#' @param index.true.associated a vector of index for true associations
#' @param threshold.signi a scalar
#' @return the Jaccard index 

compute.jaccard = function(p.values, index.true.associated, threshold.signi=0.05){
  
  index.signi.p.values = which(p.values<=threshold.signi)
  
  sum(index.signi.p.values%in%index.true.associated)/
    (length(index.signi.p.values)+length(index.true.associated)-sum(index.signi.p.values%in%index.true.associated))
}

PoTA = function(spearman.correlation, Associations, threshold=0.05){
  
  associations = lapply(1:length(names(Associations)), function(x){
    r = as.numeric(names(Associations)[x])
    sapply(Associations[[x]], function(y){
      c(r, y)
    })
  })
  
  df.associations = do.call("cbind", associations)
  
  titi = lapply(1:length(spearman.correlation), function(x) which(spearman.correlation[[x]]<=threshold))
  i.titi = which(sapply(titi, length)>0)
  
  ss = c()
  for(i in i.titi){
    
    for(j in titi[[i]]){
      fd = c(i, j)
      s = c()
      for(p in 1:ncol(df.associations)){
        s = c(s,all(fd==df.associations[,p]))
      }
      
      ss = c(ss,sum(s))
    }
    
  }
  
  return(sum(ss)/(length(ss)+ncol(dfp)-sum(ss)))
  
}

#This function computes the proportion of true associations between pairs of elements
#'@param cor.test.by.element a object returned by cor.test 
#'@param threshold.signi a scalar
#'@param index.X a vector of index for associated elements in X
#''@param index.Y a vector of index for associated elements in Y 
#'@return the proportion of true associations

find.true.associations = function(cor.test.by.element,threshold.signi=0.05, index.X, index.Y){
  titi = lapply(cor.test.by.element, function(x) {
    which(sapply(x, function(y){
      y$p.value<=threshold.signi
    }))
  })
  
  pairs = list()
  
  for(i in 1:length(index.X)){
    a = c(index.X[i],
          index.Y[i])
    pairs[[i]] = a
  }
  
  index.non.null = which(sapply(titi, function(x) length(x))>0)
  
  l=c()
  
  for(i in index.non.null){
    
    u = c()
    for(j in titi[[i]]){
      c = c(i,j)
      n = c()
      for(p in 1:length(index.X)){
        n = c(n,all(c  == pairs[[p]]))
      }
      
      u = c(u,sum(n))
      
    }
    
    l = c(l,sum(n))
  }
  
  sum(l)/(length(index.X)+length(l)-sum(l))
}

#This function builds the confusion matrix for a given significance threshold
#'@param p.values a vector of p-values
#'@param true.associated a vector of index of true associated
#'@param true.non.associated a vector of index of true non-associated
#'@param signi.threshold a scalar
#'

confusion.matrix = function(pvalues, true.associated, true.non.associated ,signi.threshold=0.05){
  
  significant.associations = which(pvalues<=signi.threshold)
  non.significant.associations = which(pvalues>signi.threshold)
  
  TP = sum(significant.associations%in%true.associated)
  TN = sum(non.significant.associations%in%true.non.associated)
  FP = sum(significant.associations%in%true.non.associated)
  FN = sum(non.significant.associations%in%true.associated)
  
  return(matrix(c(TP, FN,FP, TN), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg"))))
}

#This function computes the sensivity, specificity from a confusion matrix
#'@param confusion.matrix a matrix built by the confusion.matrix function 
#'@return a vector with sensitivity (1) and specificity (2) 
#'
sensivity.specificity.from.confusion.matrix = function(confusion.matrix){
  
  sensitivity = confusion.matrix[1,1]/(confusion.matrix[1,1]+confusion.matrix[2,1])
  specificity = confusion.matrix[2,2]/(confusion.matrix[2,2]+confusion.matrix[1,2])
  
  return(c(sensitivity,specificity))
}
       
#This function computes the F1-score (harmonic mean of sensitivity and specificity) from a confusion matrix
#'@param confusion.matrix a matrix built by the confusion.matrix function 
#'@return a scalar corresponding to the F1-score
#'
f1.score = function(confusion.matrix){
  m = sensivity.specificity.from.confusion.matrix(confusion.matrix)
  
  return(2/((1/m[1])+(1/m[2])))
}       


#This function builds the ROC curve 
#'@param pvalues a vector of p-values
#'@param true.associated a vector of index for true associated elements
#'@param true.non.associated a vector of index for true non-associated elements
#'@param signi.levels a vector with all the significance threshold levels
#'@return a plot with ROC

build.ROC = function(pvalues, true.associated, true.non.associated, signi.levels = seq(0,1, 0.0001)){
  
  confusion.matrix.for.all.signi = lapply(signi.levels, function(rho) confusion.matrix(pvalues,true.associated,true.non.associated,rho))
  
  specificity.sensitivity.for.all.confusion.matrix = sapply(confusion.matrix.for.all.signi, function(x){
    sensivity.specificity.from.confusion.matrix(x)
  })
  
  specificity.sensitivity.for.all.confusion.matrix = data.frame(t(specificity.sensitivity.for.all.confusion.matrix))
  
  colnames(specificity.sensitivity.for.all.confusion.matrix) = c("Sensitivity","Specificity")
  specificity.sensitivity.for.all.confusion.matrix$One_Minus_Specificity = 1 - specificity.sensitivity.for.all.confusion.matrix$Specificity
  
  ggpar(ggline(specificity.sensitivity.for.all.confusion.matrix, x="One_Minus_Specificity", y="Sensitivity", numeric.x.axis = TRUE), xlim = c(0,1), ylim=c(0,1), xlab = "1-Specificity", ylab="Sensitivity")+geom_abline(intercept = 0, slope=1, linetype=2)+
    annotate("text",x=0.8,y=0.1, label=paste0("AUC=", round(compute.auc(pvalues, true.associated),2)))
}


#This function computes the area under the curve
#'
#'@param pvalues a vector of p-values
#'@param true.associated a vector of true associations
#'@return a numeric corresponding to the area under the curve
#'
compute.auc = function(pvalues,true.associated){
  ert = 1:length(pvalues)
  
  ert = as.numeric(ert%in%true.associated)
  auc(roc(ert,pvalues))
}

#Lasso-ZINB 
find.true.associations.lasso = function(lasso.output, index.X, index.Y){
  
  names(lasso.output) = 1:length(lasso.output)
  non.null.coeffs.lasso = lapply(lasso.output, function(x) which(x[-1]!=0))
  non.null.coeffs.lasso = non.null.coeffs.lasso[which(sapply(non.null.coeffs.lasso, length)>0)]
  
  pairs = list()
  
  for(i in 1:length(index.Y)){
    a = c(index.Y[i],
          index.X[i])
    pairs[[i]] = a
  }
  
  
  l=c()
  
  for(i in names(non.null.coeffs.lasso)){
    u = c()
    for(j in non.null.coeffs.lasso[[i]]){
      c = c(i,j)
      n = c()
      for(p in 1:length(index.Y)){
        n = c(n,all(c  == pairs[[p]]))
      }
      
      u = c(u,sum(n))
    }
    l = c(l,sum(n))
  }
  
  sum(l)/(length(index.Y)+length(l)-sum(l))
}


#CCA
#This function returns different types of redundancy from a CCA object
#'@param X a data.frame
#'@param Y a data.frame
#'@param cca.object a object generated by CCA::cc
#'@return a list with all types of redundancy
#'
compute.redundancy.cca = function(X,Y,cca.object){
  
  #Proportion of variance by each variate
  variance.explained.X = colSums(cca.object$scores$corr.X.xscores^2)/ncol(X)
  variance.explained.Y = colSums(cca.object$scores$corr.Y.yscores^2)/ncol(Y)
  
  #Proportion of variance in the first pair explained by other member of the pair
  proportion.variance.canonical.cor = cca.object$cor^2
  
  #Proportion of X explained by the canonical variates 
  cond.redundancy.X = variance.explained.X*proportion.variance.canonical.cor
  cond.redundancy.Y = variance.explained.Y*proportion.variance.canonical.cor
  
  
  return(list("Redundancy" = list("X"=variance.explained.X, "Y"=variance.explained.Y),"VarianceCanonicalCor" = proportion.variance.canonical.cor,"ConditionalRedundancy" = list("X" = cond.redundancy.X, "Y"=cond.redundancy.Y)))
}


#PLS
#'@param X a data.frame
#'@param Y a data.frame
#'@param pls.object an object return by mixOmics::pls
#'@return a list with all types of redundancy
#'
compute.redundancy.pls = function(X,Y,pls.object){
  
  
  X.scaled = apply(X, 2, function(col) (col-mean(col))/sd(col))
  Y.scaled = apply(Y, 2, function(col) (col-mean(col))/sd(col))
  
  m = min(c(ncol(X), ncol(Y)))
  #Proportion of variance by each variate
  variance.explained.X = sapply(1:m, function(y) sum(sapply(1:m, function(x) cor(X.scaled[,x], pls.object$variates$X[,y])^2))/ncol(X))
  variance.explained.Y = sapply(1:m, function(y) sum(sapply(1:m, function(x) cor(Y.scaled[,x], pls.object$variates$Y[,y])^2))/ncol(Y))
  
  #Proportion of variance in the first pair explained by other member of the pair
  proportion.variance.canonical.cor = sapply(1:m, function(x) cor(pls.object$variates$X[,x], pls.object$variates$Y[,x])^2)
  
  #Proportion of variance explained by the canonical variates
  cond.redundancy.X = variance.explained.X * proportion.variance.canonical.cor
  cond.redundancy.Y = variance.explained.Y * proportion.variance.canonical.cor
  
  return(list("Redundancy" = list("X"=variance.explained.X, "Y"=variance.explained.Y),"VarianceCanonicalCor" = proportion.variance.canonical.cor,"ConditionalRedundancy" = list("X" = cond.redundancy.X, "Y"=cond.redundancy.Y)))
}

#This function computes the Jaccard index for a cca object
#' @param cca.object a pls object returned by CCA::cc
#' @param index.true.associated a vector of index for true associations
#' @return the Jaccard index 

compute.jaccard.cca = function(cca.object, index.true.associated, type=c("X","Y") ,index.component = 1:2){
  
  keep = length(index.true.associated)
  
  if(type=="X"){
    order.index.by.component = lapply(index.component, function(n){
      order(abs(cca.object$scores$corr.X.xscores[,n]),decreasing = T)[1:keep]})
  }
  
  else if(type=="Y"){
    order.index.by.component = lapply(index.component, function(n){
      order(abs(cca.object$scores$corr.Y.yscores[,n]),decreasing = T)[1:keep]})
  }
  
  sapply(order.index.by.component, function(x) {
    sum(x%in%index.true.associated)/
      (length(x)+length(index.true.associated)-sum(x%in%index.true.associated))
  })
}


#This function computes the Jaccard index for a pls object
#' @param pls.object a pls object returned by mixOmics::pls
#' @param index.true.associated a vector of index for true associations
#' @return the Jaccard index 

compute.jaccard.pls = function(pls.object, index.true.associated, type=c("X","Y") ,index.component = 1:2){
  
  keep = length(index.true.associated)
  
  if(type=="X"){
    order.index.by.component = lapply(index.component, function(n){
      order(abs(pls.object$loadings$X[,n]),decreasing = T)[1:keep]})
  }
  
  else if(type=="Y"){
    order.index.by.component = lapply(index.component, function(n){
      order(abs(pls.object$loadings$Y[,n]),decreasing = T)[1:keep]})
  }
  
  sapply(order.index.by.component, function(x) {
    sum(x%in%index.true.associated)/
    (length(x)+length(index.true.associated)-sum(x%in%index.true.associated))
  })
}


#RDA

#This function computes the proportion of explained variance from Redundancy analysis
#'@param rda.object a rda object from vegan rda function
#'@return the proportion of explained variance 
#'
compute.explained.variance.rda = function(rda.object){
  
  return(rda.object$CCA$tot.chi/ rda.object$tot.chi)
}
