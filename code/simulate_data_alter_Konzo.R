#This function simulates data based under the alternative for our Konzo Scenario

#scenario == Konzo: list.parameters should have 5 elements
#2 sparse correlation matrices
#mean for Poisson
#and the two parameters for the NB (mu, size)

#list.parameters.Konzo can be found in the input\Konzo directory

list.parameters.Konzo = readRDS("input\\Konzo\\list_parameters_Konzo.RDS")


simulate.realistic.data.alter.Konzo = function(list.parameters, nrep=1000, nindiv=c(100,100), random.seed=1234){
  
  #nindiv is a vector where the first element is the number of individual in species, the second in metabolites
  
  Cov.Species.Konzo = list.parameters$Cov_species_Konzo
  Cov.Metabolites.Konzo = list.parameters$Cov_metabolites_Konzo
    
  #Mean estimated by MLE 
  lambda.metabolites.Konzo = list.parameters$lambda_metabolites
  
  #NB parameters estimated by MLE
  size.species.Konzo = list.parameters$size_species
  mu.species.Konzo = list.parameters$mu_species
  
  set.seed(random.seed) #For reproducibility
  
  alter.rep.Konzo = lapply(1: nrep, function(rep){
    print(rep)
    
    M1 = MASS::mvrnorm(nindiv[1], rep(0,ncol(Cov.Species.Konzo)), Cov.Species.Konzo)
    M2 = MASS::mvrnorm(nindiv[2], rep(0,ncol(Cov.Metabolites.Konzo)), Cov.Metabolites.Konzo)
    
    
    prop.associated.metabolites = round(runif(1,0.01,0.1),2)
    prop.associated.species = round(runif(1,0.01,0.1),2)
    
    prop.impacted.metabolites = round(runif(1,0.01,0.1),2)
    prop.impacted.species = round(runif(1,0.01,0.1),2)
    
    n.impacted.metabolites = round(ncol(Cov.Metabolites.Konzo)*prop.impacted.metabolites)
    n.impacted.species = round(ncol(Cov.Species.Konzo)*prop.impacted.species)
    
    n.associated.metabolites = ceiling(ncol(Cov.Metabolites.Konzo)*prop.associated.metabolites)
    n.associated.species = ceiling(ncol(Cov.Species.Konzo)*prop.associated.species)
    
    index.associated.species = sample(1:ncol(Cov.Species.Konzo), n.associated.species)
    index.associated.metabolites = sample(1:ncol(Cov.Metabolites.Konzo), n.associated.metabolites )
    
    index.impacted.species = sample(1:ncol(Cov.Species.Konzo), n.impacted.species )
    index.impacted.metabolites = sample(1:ncol(Cov.Metabolites.Konzo),  n.impacted.metabolites)
    
    vb=lapply(index.impacted.metabolites, function(x) sample(index.associated.species, sample(1:n.associated.species,1)))
    wb = lapply(vb, function(x) rnorm(length(x), 0,0.2))
    vn=lapply(index.impacted.species, function(x) sample(index.associated.metabolites, sample(1:n.associated.metabolites,1)))
    wn = lapply(vn, function(x) rnorm(length(x), 0,0.2))
    
    
    simulated.microbiotes = matrix(qnbinom(p = pnorm(M1), size = size.species.Konzo, mu = mu.species.Konzo), ncol=ncol(Cov.Species.Konzo),nrow=nindiv[1])
    simulated.metabolites = matrix(qpois(p = pnorm(M2), lambda = lambda.metabolites.Konzo), ncol=ncol(Cov.Metabolites.Konzo),nrow=nindiv[2])
    
    generate.assos.metabolites = sapply(1:length(index.impacted.metabolites), function(x) {
      mu = exp(simulated.microbiotes[,vb[[x]],drop=FALSE]%*%wb[[x]])
      mu[mu>1e+7] = 1e+7
      rpois(nrow(simulated.metabolites),mu)
    })
    
    generate.assos.species = sapply(1:length(index.impacted.species), function(x) {
      mu = exp(simulated.metabolites[,vn[[x]],drop=FALSE]%*%wn[[x]])
      mu[is.infinite(mu)] = max(mu[!is.infinite(mu)])
      
      mu[mu>1e+7] = 1e+7
      rnbinom(nrow(simulated.microbiotes), mu = mu, size=size.species.Konzo[index.impacted.species[x]])
    })
    
    simulated.metabolites[, index.impacted.metabolites] = generate.assos.metabolites
    simulated.microbiotes[, index.impacted.species] = generate.assos.species
    
    colnames(simulated.microbiotes) = paste0("Micro", 1:ncol(simulated.microbiotes))
    colnames(simulated.metabolites) = paste0("Metabo", 1:ncol(simulated.metabolites))
    
    list("Simulated.Microbiotes"=simulated.microbiotes,"Simulated.Metabolites"=simulated.metabolites, "index.associated.species"=index.associated.species, "index.associated.metabolites"=index.associated.metabolites, "index.impacted.species"=index.impacted.species, "index.impacted.metabolites"=index.impacted.metabolites, "associations.metabolites"=vb, "associations.species"=vn)
    })
  
  return(alter.rep.Konzo)
}
