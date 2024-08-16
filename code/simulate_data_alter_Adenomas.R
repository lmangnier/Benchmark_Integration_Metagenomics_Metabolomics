#This function simulates data based under the alternative for our Adenomas Scenario

#scenario == Adenomas: list.parameters should have 6 elements
#2 sparse correlation matrices
#mean for Gaussian 
#and the three parameters for the ZINB (mu, size and pstr0)

#list.parameters.Adenomas can be found in the input\Adenomas directory

list.parameters.Adenomas = readRDS("input\\Adenomas\\list_parameters_Adenomas.RDS")

simulate.realistic.data.alter.Adenomas(list.parameters.Adenomas, nrep = 10)
simulate.realistic.data.alter.Adenomas = function(list.parameters, nrep=1000, nindiv=c(100,100), random.seed=1234){
  
  #nindiv is a vector where the first element is the number of individual in species, the second in metabolites
  
  #Sparse correlation matrices estimated by Spiec-Easi
  Cov.Species.Adenomas = list.parameters$Cov_species_Adenomas
  Cov.Metabolites.Adenomas = list.parameters$Cov_metabolites_Adenomas
  
  #Mean estimated by MLE 
  mean.model.metabolites.Adenomas = list.parameters$mean_metabolites
  
  #ZINB parameters estimated by MLE
  size.species.Adenomas = list.parameters$size_species
  prop.zeros.species.Adenomas = list.parameters$propzeros_species
  mu.species.Adenomas = list.parameters$mu_species
  
  set.seed(random.seed) #For reproducibility
  
  alter.rep.Adenomas = lapply(1: nrep, function(rep){
    print(rep)
    
    M1 = MASS::mvrnorm(nindiv[1], rep(0,ncol(Cov.Species.Adenomas)), Cov.Species.Adenomas)
    simulated.metabolites = MASS::mvrnorm(nindiv[2], mean.model.metabolites.Adenomas, Cov.Metabolites.Adenomas)
    
    
    prop.associated.metabolites = round(runif(1,0.01,0.1),2)
    prop.associated.species = round(runif(1,0.01,0.1),2)
    
    prop.impacted.metabolites = round(runif(1,0.01,0.1),2)
    prop.impacted.species = round(runif(1,0.01,0.1),2)
    
    n.impacted.metabolites = round(ncol(Cov.Metabolites.Adenomas)*prop.impacted.metabolites)
    n.impacted.species = round(ncol(Cov.Species.Adenomas)*prop.impacted.species)
    
    n.associated.metabolites = ceiling(ncol(Cov.Metabolites.Adenomas)*prop.associated.metabolites)
    n.associated.species = ceiling(ncol(Cov.Species.Adenomas)*prop.associated.species)
    
    index.associated.species = sample(1:ncol(Cov.Species.Adenomas), n.associated.species)
    index.associated.metabolites = sample(1:ncol(Cov.Metabolites.Adenomas), n.associated.metabolites )
    
    index.impacted.species = sample(1:ncol(Cov.Species.Adenomas), n.impacted.species )
    index.impacted.metabolites = sample(1:ncol(Cov.Metabolites.Adenomas),  n.impacted.metabolites)
    
    vb=lapply(index.impacted.metabolites, function(x) sample(index.associated.species, sample(1:n.associated.species,1)))
    wb = lapply(vb, function(x) rnorm(length(x), 0,0.2))
    vn=lapply(index.impacted.species, function(x) sample(index.associated.metabolites, sample(1:n.associated.metabolites,1)))
    wn = lapply(vn, function(x) rnorm(length(x), 0,0.2))
    
    
    simulated.microbiotes = matrix(VGAM::qzinegbin(p = pnorm(M1), size = size.species.Adenomas, munb  = mu.species.Adenomas, pstr0 = prop.zeros.species.Adenomas), ncol=ncol(Cov.Species.Adenomas),nrow=nindiv[1])
    
    generate.assos.metabolites = sapply(1:length(index.impacted.metabolites), function(x) {
      mu = simulated.microbiotes[,vb[[x]],drop=FALSE]%*%wb[[x]]
      rnorm(nindiv[2],mu, sqrt(Cov.Metabolites.Adenomas[index.impacted.metabolites[x],index.impacted.metabolites[x]]))
    })
    
    generate.assos.species = sapply(1:length(index.impacted.species), function(x) {
      mu = exp(simulated.metabolites[,vn[[x]],drop=FALSE]%*%wn[[x]])
      mu[is.infinite(mu)] = max(mu[!is.infinite(mu)])
      
      mu[mu>1e+7] = 1e+7
      VGAM::rzinegbin(nindiv[1], munb = mu, size=size.species.Adenomas[index.impacted.species[x]], pstr0 = prop.zeros.species.Adenomas[index.impacted.species[x]])
    })
    
    simulated.metabolites[, index.impacted.metabolites] = generate.assos.metabolites
    simulated.microbiotes[, index.impacted.species] = generate.assos.species
    
    colnames(simulated.microbiotes) = paste0("Micro", 1:ncol(simulated.microbiotes))
    colnames(simulated.metabolites) = paste0("Metabo", 1:ncol(simulated.metabolites))
    
    list("Simulated.Microbiotes"=simulated.microbiotes,"Simulated.Metabolites"=simulated.metabolites, "index.associated.species"=index.associated.species, "index.associated.metabolites"=index.associated.metabolites, "index.impacted.species"=index.impacted.species, "index.impacted.metabolites"=index.impacted.metabolites, "associations.metabolites"=vb, "associations.species"=vn)
  })
  
  return(alter.rep.Adenomas)
}
