#This function simulates data for benchmark analysis under the NORtA algorithm
#'@param n.replicates is an integer corresponding to the number of replicates generated
#'@param n.indiv is an integer corresponding to the number of individuals
#'@param n.Microbiotes is an integer corresponding to the number of Microbiotes
#'@param n.Metabolites is an integer corresponding to the number of Metabolites
#'@param n.associated.Microbiotes is an integer corresponding to the number of associated Microbiotes
#'@param n.associated.Metabolites is an integer corresponding to the number of associated Metabolites
#'@param prop.associated.Microbiotes is a vector corresponding to the minimal/maximal proportion of associated Microbiotes
#'@param prop.associated.Metabolites is a vector corresponding to the minimal/maximal proportion of associated Metabolites
#'@param n.impacted.Microbiotes is an integer corresponding to the number of impacted Microbiotes
#'@param n.impacted.Metabolites is an integer corresponding to the number of impacted Metabolites
#'@param prop.impacted.Microbiotes is a vector corresponding to the minimal/maximal proportion of impacted Microbiotes
#'@param prop.impacted.Metabolites is a vector corresponding to the minimal/maximal proportion of impacted Metabolites
#'@param fixed a boolean whether the considered effect is fixed or drawn for a Gaussian N(0, 0.2)
#'@param effect a numeric corresponding to the effect size
#'@param output.filename the filename corresponding to the file where are stored


simulate_data = function(n.replicates=1, n.indiv=100, n.Microbiotes=25, n.Metabolites=25, 
                         n.associated.Microbiotes=3, n.associated.Metabolites=3, prop.associated.Microbiotes=NULL, prop.associated.Metabolites=NULL, 
                         n.impacted.Microbiotes=3, n.impacted.Metabolites=3, prop.impacted.Microbiotes=NULL, prop.impacted.Metabolites=NULL, fixed=FALSE, effect=NULL, output.filename=NULL){
  
  l = list()
  
  list.index = list()
  list.associations = list()
  list.Microbiotes = list()
  list.Metabolites = list()
  list.effects = list()
  
  if(!is.null(prop.associated.Microbiotes)&!is.null(n.associated.Microbiotes)) stop("Please chose either a fixed number of associated microbiotes or a range of proportions, not both...")
  if(!is.null(prop.associated.Metabolites)&!is.null(n.associated.Metabolites)) stop("Please chose either a fixed number of associated metabolites or a range of proportions, not both...")

  if(!is.null(prop.impacted.Microbiotes)&!is.null(n.impacted.Microbiotes)) stop("Please chose either a fixed number of impacted microbiotes or a range of proportions, not both...")
  if(!is.null(prop.impacted.Metabolites)&!is.null(n.impacted.Metabolites)) stop("Please chose either a fixed number of impacted metabolites or a range of proportions, not both...")
  
  if(fixed==TRUE&is.null(effect)) stop("No value for the effect parameter, please provide a correct value")
  if(fixed==FALSE&!is.null(effect)) stop("Effect parameter passed, while fixed = F...")
  
  for(rep in 1:n.replicates){
    L0 = matrix(0, ncol=n.Microbiotes, nrow=n.Microbiotes)
    #Variance for microbiotes throughout samples
    #We can play on this parameter in order to assess biological variability
    diag(L0) = runif(n.Microbiotes,1.5,2.5)
  
    #Off-diagonal elements are randomly selected to have either 0 covariance or a positive or negative covariance based on uniform distribution
    L0[lower.tri(L0)] = sapply(1:length(L0[lower.tri(L0)]), function(x) sample(c(0,runif(1,-1.5,1.5)),1, prob = c(0.7,0.3)))
  
    Precision0 = L0%*%t(L0)
  
    #We obtain Covariance matrix based on Cholesky decomposition of lower triangular matrix
    Sigma0 = solve(Precision0)
    Cor0 = cov2cor(Sigma0)
  
    multi.norm = MASS::mvrnorm(n.indiv, rep(0,n.Microbiotes), Cor0)
    
    L1 = matrix(0, ncol=n.Metabolites, nrow=n.Metabolites)
    #Variance for microbiotes throughout samples
    #We can play on this parameter in order to assess biological variability
    diag(L1) = runif(n.Metabolites,1.5,2.5)
  
    #Off-diagonal elements are randomly selected to have either 0 covariance or a positive or negative covariance based on uniform distribution
    L1[lower.tri(L1)] = sapply(1:length(L1[lower.tri(L1)]), function(x) sample(c(0,runif(1,1,1.5)),1, prob = c(0.7,0.3)))
  
    Precision1 = L1%*%t(L1)
  
    #We obtain Covariance matrix based on Cholesky decomposition of lower triangular matrix
    Sigma1 = solve(Precision1)
    Cor1 = cov2cor(Sigma1)
  
    #The multivariate normal distribution is generated for 100 individuals with mean 0 and the
    #Correlation structure
    multi.norm1 = MASS::mvrnorm(n.indiv, rep(0,n.Metabolites), Cor1)
  
  
    if(!is.null(prop.associated.Microbiotes) & !is.null(prop.associated.Metabolites) & !is.null(prop.impacted.Microbiotes)&  !is.null(prop.impacted.Metabolites)){
          #print("Scenario1")
          n1 = ceiling(n.Microbiotes*round(runif(1,prop.associated.Microbiotes[1],prop.associated.Microbiotes[2]),2))
          n2 = ceiling(n.Metabolites*round(runif(1,prop.impacted.Metabolites[1],prop.impacted.Metabolites[2]),2))
          n3 = ceiling(n.Metabolites*round(runif(1,prop.associated.Metabolites[1],prop.associated.Metabolites[2]),2))
          n4 = ceiling(n.Microbiotes*round(runif(1,prop.impacted.Microbiotes[1],prop.impacted.Microbiotes[2]),2))
          
          
          
          index.associated.microbiote = sample(1:n.Microbiotes, n1)
          index.impacted.metabolite = sample(1:n.Metabolites, n2)
  
  
          index.associated.metabolite = sample(c(1:n.Metabolites)[-index.impacted.metabolite], n3)
          index.impacted.microbiote = sample(c(1:n.Microbiotes)[-index.associated.microbiote], n4)
    }
  
    else if(!is.null(n.associated.Microbiotes) & !is.null(prop.associated.Metabolites) & !is.null(prop.impacted.Microbiotes)&  !is.null(prop.impacted.Metabolites) ){
      #print("Scenario2")
      
      n1 = n.associated.Microbiotes
      n2 = ceiling(n.Metabolites*round(runif(1,prop.impacted.Metabolites[1],prop.impacted.Metabolites[2]),2))
      n3 = ceiling(n.Metabolites*round(runif(1,prop.associated.Metabolites[1],prop.associated.Metabolites[2]),2))
      n4 = ceiling(n.Microbiotes*round(runif(1,prop.impacted.Microbiotes[1],prop.impacted.Microbiotes[2]),2))
      
      index.associated.microbiote = sample(1:n.Microbiotes, n.associated.Microbiotes)
      index.impacted.metabolite = sample(1:n.Metabolites, n2)
      
      
      index.associated.metabolite = sample(c(1:n.Metabolites)[-index.impacted.metabolite], n3)
      index.impacted.microbiote = sample(c(1:n.Microbiotes)[-index.associated.microbiote], n4)
    }
    
    else if(!is.null(n.associated.Microbiotes) & !is.null(n.associated.Metabolites) & !is.null(prop.impacted.Microbiotes)&  !is.null(prop.impacted.Metabolites) ){
      #print("Scenario3")
      
      n1 = n.associated.Microbiotes
      n2 = ceiling(n.Metabolites*round(runif(1,prop.impacted.Metabolites[1],prop.impacted.Metabolites[2]),2))
      n3 = n.associated.Metabolites
      n4 = ceiling(n.Microbiotes*round(runif(1,prop.impacted.Microbiotes[1],prop.impacted.Microbiotes[2]),2))
      
      index.associated.microbiote = sample(1:n.Microbiotes, n1)
      index.impacted.metabolite = sample(1:n.Metabolites, n2)
      
      
      index.associated.metabolite = sample(c(1:n.Metabolites)[-index.impacted.metabolite], n3)
      index.impacted.microbiote = sample(c(1:n.Microbiotes)[-index.associated.microbiote], n4)
    }
    
    else if(!is.null(n.associated.Microbiotes) & !is.null(n.associated.Metabolites) & !is.null(n.impacted.Microbiotes)&  !is.null(prop.impacted.Metabolites) ){
      #print("Scenario4")
      
      n1 = n.associated.Microbiotes
      n2 = ceiling(n.Metabolites*round(runif(1,prop.impacted.Metabolites[1],prop.impacted.Metabolites[2]),2))
      n3 = n.associated.Metabolites
      n4 = n.impacted.Microbiotes
      
      index.associated.microbiote = sample(1:n.Microbiotes, n.associated.Microbiotes)
      index.impacted.metabolite = sample(1:n.Metabolites, round(n.Metabolites*round(runif(1,prop.impacted.Metabolites[1],prop.impacted.Metabolites[2]), 2)))
      
      
      index.associated.metabolite = sample(c(1:n.Metabolites)[-index.impacted.metabolite], n.associated.Metabolites)
      index.impacted.microbiote = sample(c(1:n.Microbiotes)[-index.associated.microbiote], n.impacted.Microbiotes)
    }
    else if(!is.null(n.associated.Microbiotes) & !is.null(n.associated.Metabolites) & !is.null(n.impacted.Microbiotes)&  !is.null(n.impacted.Metabolites) ){
      #print("Scenario5")
      
      n1 = n.associated.Microbiotes
      n2 = n.impacted.Metabolites
      n3 = n.associated.Metabolites
      n4 = n.associated.Metabolites
      
      index.associated.microbiote = sample(1:n.Microbiotes, n1)
      index.impacted.metabolite = sample(1:n.Metabolites, n2)
      
      
      index.associated.metabolite = sample(c(1:n.Metabolites)[-index.impacted.metabolite], n3)
      index.impacted.microbiote = sample(c(1:n.Microbiotes)[-index.associated.microbiote], n4)
    }
    
    else if(!is.null(n.associated.Microbiotes) & !is.null(n.associated.Metabolites) & !is.null(prop.impacted.Microbiotes)&  !is.null(n.impacted.Metabolites) ){
      #print("Scenario6")
      
      n1 = n.associated.Microbiotes
      n2 = n.impacted.Metabolites
      n3 = n.associated.Metabolites
      n4 = ceiling(n.Microbiotes*round(runif(1,prop.impacted.Microbiotes[1],prop.impacted.Microbiotes[2]),2))
      
      index.associated.microbiote = sample(1:n.Microbiotes, n1)
      index.impacted.metabolite = sample(1:n.Metabolites, n2)
      
      
      index.associated.metabolite = sample(c(1:n.Metabolites)[-index.impacted.metabolite], n3)
      index.impacted.microbiote = sample(c(1:n.Microbiotes)[-index.associated.microbiote], n4)
    }
    
    else if(!is.null(n.associated.Microbiotes) & !is.null(prop.associated.Metabolites) & !is.null(n.impacted.Microbiotes)&  !is.null(prop.impacted.Metabolites) ){
      #print("Scenario7")
      
      n1 = n.associated.Microbiotes
      n2 = ceiling(n.Metabolites*round(runif(1,prop.impacted.Metabolites[1],prop.impacted.Metabolites[2]),2))
      n3 = ceiling(n.Metabolites*round(runif(1,prop.associated.Metabolites[1],prop.associated.Metabolites[2]),2))
      n4 = ceiling(n.Microbiotes*round(runif(1,prop.impacted.Microbiotes[1],prop.impacted.Microbiotes[2]),2))
      
      index.associated.microbiote = sample(1:n.Microbiotes, n1)
      index.impacted.metabolite = sample(1:n.Metabolites, n2)
      
      
      index.associated.metabolite = sample(c(1:n.Metabolites)[-index.impacted.metabolite], n3)
      index.impacted.microbiote = sample(c(1:n.Microbiotes)[-index.associated.microbiote], n4)
    }
    
    else if(!is.null(prop.associated.Microbiotes) & !is.null(prop.associated.Metabolites) & !is.null(n.impacted.Microbiotes)&  !is.null(prop.impacted.Metabolites) ){
      #print("Scenario8")
      
      n1 = ceiling(n.Microbiotes*round(runif(1,prop.associated.Microbiotes[1],prop.associated.Microbiotes[2]),2))
      n2 = ceiling(n.Metabolites*round(runif(1,prop.impacted.Metabolites[1],prop.impacted.Metabolites[2]),2))
      n3 = ceiling(n.Metabolites*round(runif(1,prop.associated.Metabolites[1],prop.associated.Metabolites[2]),2))
      n4 = n.impacted.Microbiotes
      
      index.associated.microbiote = sample(1:n.Microbiotes, n1)
      index.impacted.metabolite = sample(1:n.Metabolites,n2)
      
      
      index.associated.metabolite = sample(c(1:n.Metabolites)[-index.impacted.metabolite], n3)
      index.impacted.microbiote = sample(c(1:n.Microbiotes)[-index.associated.microbiote], n4)
    }
    
    else if(!is.null(prop.associated.Microbiotes) & !is.null(prop.associated.Metabolites) & !is.null(prop.impacted.Microbiotes)&  !is.null(n.impacted.Metabolites) ){
      #print("Scenario9")
      
      n1 = ceiling(n.Microbiotes*round(runif(1,prop.associated.Microbiotes[1],prop.associated.Microbiotes[2]),2))
      n2 = n.impacted.Metabolites
      n3 = ceiling(n.Metabolites*round(runif(1,prop.associated.Metabolites[1],prop.associated.Metabolites[2]),2))
      n4 = ceiling(n.Microbiotes*round(runif(1,prop.impacted.Microbiotes[1],prop.impacted.Microbiotes[2]),2))
      
      
      index.associated.microbiote = sample(1:n.Microbiotes, n1)
      index.impacted.metabolite = sample(1:n.Metabolites, n2)
      
      
      index.associated.metabolite = sample(c(1:n.Metabolites)[-index.impacted.metabolite], n3)
      index.impacted.microbiote = sample(c(1:n.Microbiotes)[-index.associated.microbiote], n4)
    }
    
    if(length(index.associated.microbiote)>1 & length(index.associated.metabolite)>1){
      n.associated.by.metabolite = sapply(index.impacted.metabolite, function(x) sample(1:length(index.associated.microbiote),1))
      lhg = lapply(1:length(index.impacted.metabolite), function(x) sample(index.associated.microbiote,n.associated.by.metabolite[x]))
      
      n.associated.by.microbiote = sapply(index.impacted.microbiote, function(x) sample(1:length(index.associated.metabolite),1))
      lhh = lapply(1:length(index.impacted.microbiote), function(x) sample(index.associated.metabolite,n.associated.by.microbiote[x]))
    }
    
    else if(length(index.associated.microbiote)==1){
      #n.associated.by.metabolite = sapply(index.impacted.metabolite, function(x) sample(1:length(index.associated.microbiote),1))
      lhg = lapply(1:length(index.impacted.metabolite), function(x) index.associated.microbiote)
      
      n.associated.by.microbiote = sapply(index.impacted.microbiote, function(x) sample(1:length(index.associated.metabolite),1))
      lhh = lapply(1:length(index.impacted.microbiote), function(x) sample(index.associated.metabolite,n.associated.by.microbiote[x]))
    }
    
    else if(length(index.associated.metabolite)==1){
      n.associated.by.metabolite = sapply(index.impacted.metabolite, function(x) sample(1:length(index.associated.microbiote),1))
      lhh = lapply(1:length(index.impacted.microbiote), function(x) index.associated.metabolite)
      
      n.associated.by.metabolite = sapply(index.impacted.metabolite, function(x) sample(1:length(index.associated.microbiote),1))
      lhg = lapply(1:length(index.impacted.metabolite), function(x) sample(index.associated.microbiote,n.associated.by.metabolite[x]))
      
    }

    if(!fixed & is.null(effect)){
      lhg.beta = lapply(lhg, function(x) rnorm(length(x), 0,0.2))
      lhh.beta = lapply(lhh, function(x) rnorm(length(x), 0,0.2))
    }
    
    if(fixed & !is.null(effect)){
      lhg.beta = lapply(lhg, function(x) rep(effect,length(x)))
      lhh.beta = lapply(lhh, function(x) rep(effect,length(x)))
    }
  
    mu.microbiotes = replicate(n.Microbiotes,runif(n.indiv,1,5))
    mu.microbiotes[,index.associated.microbiote] = replicate(length(index.associated.microbiote), runif(n.indiv,1.5,2.5))
  
    mu.metabolites = replicate(n.Metabolites, runif(n.indiv,5,10))
    mu.metabolites[,index.associated.metabolite] = replicate(length(index.associated.metabolite), runif(n.indiv,1.5,2.5))
    
    size.microbiotes = runif(n.Microbiotes,0.3,1)
    size.metabolites = runif(n.Metabolites,1,2)
    
    simulated.microbiotes = matrix(qnbinom(pnorm(multi.norm), mu=mu.microbiotes, size=size.microbiotes), ncol=n.Microbiotes,nrow=n.indiv)
    simulated.metabolites = matrix(qnbinom(pnorm(multi.norm1), mu=mu.metabolites, size=size.metabolites), ncol=n.Metabolites,nrow=n.indiv)
  
    mu.microbiotes[,index.impacted.microbiote]=sapply(1:length(index.impacted.microbiote) , function(x) exp(simulated.metabolites[,lhh[[x]], drop=F]%*%lhh.beta[[x]]))
    mu.metabolites[,index.impacted.metabolite]=sapply(1:length(index.impacted.metabolite) , function(x) exp(simulated.microbiotes[,lhg[[x]], drop=F]%*%lhg.beta[[x]]))
  
    simulated.microbiotes = matrix(qnbinom(pnorm(multi.norm), mu=mu.microbiotes, size=size.microbiotes), ncol=n.Microbiotes,nrow=n.indiv)
    simulated.metabolites = matrix(qnbinom(pnorm(multi.norm1), mu=mu.metabolites, size=size.metabolites), ncol=n.Metabolites,nrow=n.indiv)
  
  
    names(lhg) = index.impacted.metabolite
    names(lhg.beta) = index.impacted.metabolite
  
    names(lhh) = index.impacted.microbiote
    names(lhh.beta) = index.impacted.microbiote
    
    list.index[[rep]] = list("Associated"= list("Microbiotes" = index.associated.microbiote, 
                                                "Metabolites" = index.associated.metabolite),
                             "Impacted" = list("Microbiotes" = index.impacted.microbiote, 
                                               "Metabolites" = index.impacted.metabolite))

    list.associations[[rep]] = list("Microbiote" = lhh, "Metabolite" = lhg)
    list.effects[[rep]] = list("Microbiote" = lhh.beta, "Metabolite" = lhg.beta)
    list.Microbiotes[[rep]] = simulated.microbiotes
    list.Metabolites[[rep]] = simulated.metabolites
  
  }
  
  l = list("Index" = list.index,
           "Associations" = list.associations,
           "Effects" = list.effects,
           "Simulated.Microbiotes" = list.Microbiotes,
           "Simulated.Metabolites" = list.Metabolites)
  
  if(is.null(output.filename)){
    return(l)
  }
  else{
    saveRDS(l, output.filename)
  }
  


}

simulate_data_null = function(nrep=1, nindiv=100, n.Microbiotes=25, n.Metabolites=25){
  
  list.Microbiotes = list()
  list.Metabolites = list()
  
  for(rep in 1:nrep){
    L0 = matrix(0, ncol=n.Microbiotes, nrow=n.Microbiotes)
    #Variance for microbiotes throughout samples
    #We can play on this parameter in order to assess biological variability
    diag(L0) = runif(n.Microbiotes,1.5,2.5)
    
    #Off-diagonal elements are randomly selected to have either 0 covariance or a positive or negative covariance based on uniform distribution
    L0[lower.tri(L0)] = sapply(1:length(L0[lower.tri(L0)]), function(x) sample(c(0,runif(1,-1.5,1.5)),1, prob = c(0.7,0.3)))
    
    Precision0 = L0%*%t(L0)
    
    #We obtain Covariance matrix based on Cholesky decomposition of lower triangular matrix
    Sigma0 = solve(Precision0)
    Cor0 = cov2cor(Sigma0)
    
    multi.norm = MASS::mvrnorm(nindiv, rep(0,n.Microbiotes), Cor0)
    
    L1 = matrix(0, ncol=n.Metabolites, nrow=n.Metabolites)
    #Variance for microbiotes throughout samples
    #We can play on this parameter in order to assess biological variability
    diag(L1) = runif(n.Metabolites,1.5,2.5)
    
    #Off-diagonal elements are randomly selected to have either 0 covariance or a positive or negative covariance based on uniform distribution
    L1[lower.tri(L1)] = sapply(1:length(L1[lower.tri(L1)]), function(x) sample(c(0,runif(1,1,1.5)),1, prob = c(0.7,0.3)))
    
    Precision1 = L1%*%t(L1)
    
    #We obtain Covariance matrix based on Cholesky decomposition of lower triangular matrix
    Sigma1 = solve(Precision1)
    Cor1 = cov2cor(Sigma1)
    
    #The multivariate normal distribution is generated for 100 individuals with mean 0 and the
    #Correlation structure
    multi.norm1 = MASS::mvrnorm(nindiv, rep(0,n.Metabolites), Cor1)
    
    mu.microbiotes = replicate(n.Microbiotes,runif(nindiv,1,5))
    mu.metabolites = replicate(n.Metabolites, runif(nindiv,5,10))
    
    size.microbiotes = runif(n.Microbiotes,0.3,1)
    size.metabolites = runif(n.Metabolites,1,2)
    
    simulated.microbiotes = matrix(qnbinom(pnorm(multi.norm), mu=mu.microbiotes, size=size.microbiotes), ncol=n.Microbiotes,nrow=nindiv)
    simulated.metabolites = matrix(qnbinom(pnorm(multi.norm1), mu=mu.metabolites, size=size.metabolites), ncol=n.Metabolites,nrow=nindiv)
    
    list.Microbiotes[[rep]] = simulated.microbiotes
    list.Metabolites[[rep]] = simulated.metabolites
  }
  
  return(list("Simulated.Microbiotes"=list.Microbiotes,"Simulated.Metabolites"=list.Metabolites))
}
