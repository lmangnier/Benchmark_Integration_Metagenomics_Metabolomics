#This script reproduces all the analyses shown in the wiki section
library(data.table)
library(dplyr)
library(coda4microbiome)
library(ggplot2)
library(patchwork)
library(vegan)
library(compositions)
library(MiRKAT)

#Functions 

plot_loadings = function(spls_res, titre){
  #This function plots loadings for PLS
  loadings_X <- spls_res$loadings$X %>% as.data.frame() %>% dplyr::mutate(block = "OTU")
  loadings_Y <- spls_res$loadings$Y %>% as.data.frame() %>% dplyr::mutate(block = "Metabolites")
  loadings <- dplyr::bind_rows(loadings_X, loadings_Y) %>%
    tibble::rownames_to_column("feature") %>%
    dplyr::select(-comp2) %>%
    dplyr::mutate(abs_comp1 = abs(comp1),
                  sign = sign(comp1)) %>%
    dplyr::mutate(SIGN = ifelse(sign == "-1", "-", "+")) %>%
    dplyr::arrange(desc(abs_comp1)) %>% .[1:20,] %>%
    dplyr::mutate(feature = factor(feature, levels = feature)) %>%
    dplyr::filter(comp1 != 0)
  
  ggplot(loadings, aes(y = abs_comp1, x = feature)) +
    geom_col(aes(fill = block), width = 0.5) +
    geom_label(aes(label = SIGN)) +
    coord_flip() +
    # facet_wrap(score ~ comp, scales = "free") +
    theme_bw() +
    ylab("Comp1") + xlab("Feature") +
    ggtitle(titre)
}


apply_var_coef <- function(x){   
  #This function keeps only half the number of individuals most variable features
  cv_abs <- lapply(as.data.frame(x), function(i){     
    sd(i)/mean(i)   }) %>% unlist() %>% abs() %>% sort(decreasing = TRUE)   
  most_var <- head(cv_abs, n = nrow(x)/2)   
  return(x[, names(most_var)]) 
}


apply_nzv_filter = function(x, freqCut = 95/5, uniqueCut = 10){
  #This function filters out features with few variability
  nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
  if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
  return(x)}




networks_from_CODA_LASSO = function(l){
  
  #This function plot a network from a list of CODA-LASSO objects
  library(igraph)
  df_metabo_species_CODA = do.call("rbind",lapply(1:length(l), function(x){
    if(length(l[[x]]$taxa.name)>0){
      data.frame("from"=names(l)[x], "to" = l[[x]]$taxa.name, "Weight" = l[[x]]$`log-contrast coefficients` )
      
    }
  }))
  
  graph_metabo_species_CODA = igraph::graph_from_data_frame(df_metabo_species_CODA, vertices = c(unique(df_metabo_species_CODA$from), unique(df_metabo_species_CODA$to)), directed = FALSE)
  E(graph_metabo_species_CODA)$width = abs(df_metabo_species_CODA$Weight)
  E(graph_metabo_species_CODA)$Color = ifelse(df_metabo_species_CODA$Weight<0, "green", "red")
  V(graph_metabo_species_CODA)$Type = ifelse(names(V(graph_metabo_species_CODA))%in%unique(df_metabo_species_CODA $from), "Metabolites", "Species" )
  
  col = c("Metabolites" = "#BC3C29FF", "Species" = "#0072B5FF") #same colors as in the paper
  ggnet::ggnet2(graph_metabo_species_CODA, color= "Type", color.palette =col ,size=4,edge.size="width", edge.color = "Color", color.legend = "Type")
  
}

metabolites_Adenomas = fread("metabolites_Adenomas.tsv", sep="\t", header = TRUE)
species_Adenomas = fread("genera_counts_Adenomas.tsv", sep="\t", header=TRUE)
metadata_Adenomas = fread("metadata_Adenomas.tsv", sep="\t", header=TRUE)

#Here we are removing Dataset, Publication.Name, Subject, DOI, and Fastq.Sample.ID columns in the metadata object
metadata_Adenomas = metadata_Adenomas %>% select(-c("Dataset", "Publication.Name", "DOI" ,"Subject", "Fastq.Sample.ID"))
#We created a new column grouping both Adenomas and Carcinoma in the same group
metadata_Adenomas = metadata_Adenomas %>% mutate(Phenotype = if_else(Study.Group=="Control", 0, 1))

#We proceed to some joins for facilitating analyses
metabolites_Adenomas_annotated = left_join(x=metabolites_Adenomas, y =metadata_Adenomas, by="Sample")
species_Adenomas_annotated = left_join(x=species_Adenomas, y =metadata_Adenomas, by="Sample")


#We split data based on phenotype
species_Adenomas_affected = species_Adenomas_annotated %>% filter(Phenotype == 1) %>% select(-c("Sample", colnames(metadata_Adenomas)))
species_Adenomas_unaffected  = species_Adenomas_annotated %>% filter(Phenotype == 0) %>% select(-c("Sample", colnames(metadata_Adenomas)))

metabolites_Adenomas_affected = metabolites_Adenomas_annotated %>% filter(Phenotype == 1) %>% select(-c("Sample", colnames(metadata_Adenomas)))
metabolites_Adenomas_unaffected  = metabolites_Adenomas_annotated %>% filter(Phenotype == 0) %>% select(-c("Sample", colnames(metadata_Adenomas)))


#We now apply the Mantel

metabolites_Adenomas_affected = metabolites_Adenomas_affected %>% apply_nzv_filter()
metabolites_Adenomas_unaffected = metabolites_Adenomas_unaffected %>% apply_nzv_filter()

ilr_species_Adenomas_affected = species_Adenomas_affected %>%  ilr %>% apply_nzv_filter()
ilr_species_Adenomas_unaffected = species_Adenomas_unaffected %>% ilr %>% apply_nzv_filter()


Euclidean_ilr_species_Adenomas_affected = dist(ilr_species_Adenomas_affected)
Euclidean_ilr_species_Adenomas_unaffected = dist(ilr_species_Adenomas_unaffected)


Euclidean_metabolites_Adenomas_affected = dist(metabolites_Adenomas_affected)
Euclidean_metabolites_Adenomas_unaffected = dist(metabolites_Adenomas_unaffected)


mantel(Euclidean_ilr_species_Adenomas_affected, Euclidean_metabolites_Adenomas_affected, method="spearman")
mantel(Euclidean_ilr_species_Adenomas_unaffected, Euclidean_metabolites_Adenomas_unaffected, method = "spearman")


#MMiRKAT 
covar_affected_Adenomas = metadata_Adenomas %>% filter(Phenotype==1) %>% select(c("Gender", "Chem ID / Smoking history", "Race","Age"))
covar_unaffected_Adenomas = metadata_Adenomas %>% filter(Phenotype==0) %>% select(c("Gender", "Chem ID / Smoking history", "Race", "Age"))

ilr_species_Adenomas_affected = as.data.frame(ilr_species_Adenomas_affected)
ilr_species_Adenomas_unaffected = as.data.frame(ilr_species_Adenomas_unaffected)


ilr_species_Adenomas_affected_top = ilr_species_Adenomas_affected %>% apply_var_coef()
ilr_species_Adenomas_unaffected_top = ilr_species_Adenomas_unaffected %>% apply_var_coef()

metabolites_Adenomas_affected_top = as.data.frame(metabolites_Adenomas_affected) %>% apply_var_coef()
metabolites_Adenomas_unaffected_top = as.data.frame(metabolites_Adenomas_unaffected) %>% apply_var_coef()


Euclidean_ilr_species_Adenomas_affected_top = dist(ilr_species_Adenomas_affected_top)
Euclidean_ilr_species_Adenomas_unaffected_top = dist(ilr_species_Adenomas_unaffected_top)

D2K_affected_top = D2K(as.matrix(Euclidean_ilr_species_Adenomas_affected_top))
D2K_unaffected_top = D2K(as.matrix(Euclidean_ilr_species_Adenomas_unaffected_top))

MMiRKAT(Y=as.matrix(metabolites_Adenomas_affected_top), Ks = D2K_affected_top, X = covar_affected_Adenomas, returnR2 = T,returnKRV = T)
MMiRKAT(Y=as.matrix(metabolites_Adenomas_unaffected_top), Ks = D2K_unaffected_top, X = covar_unaffected_Adenomas,returnKRV = T, returnR2 = T)


MMiRKAT(Y=as.matrix(metabolites_Adenomas_affected_top), Ks = D2K_affected_top, returnR2 = T,returnKRV = T)
MMiRKAT(Y=as.matrix(metabolites_Adenomas_unaffected_top), Ks = D2K_unaffected_top,returnKRV = T, returnR2 = T)


#MiRKAT
MiRKAT_affected = list()
MiRKAT_unaffected = list()

for(col in colnames(metabolites_Adenomas_affected)){
 MiRKAT_affected[[col]] = MiRKAT(y = metabolites_Adenomas_affected %>% select(all_of(col)) %>% unlist(),Ks=D2K(as.matrix(Euclidean_ilr_species_Adenomas_affected)), returnKRV = T, returnR2 = T)
  
}

for(col in colnames(metabolites_Adenomas_unaffected)){
  MiRKAT_unaffected[[col]] = MiRKAT(y = metabolites_Adenomas_unaffected %>% select(all_of(col)) %>% unlist(),Ks=D2K(as.matrix(Euclidean_ilr_species_Adenomas_unaffected)), returnKRV = T, returnR2 = T)
  
}

df_results_MiRKAT_affected = data.frame(do.call("rbind", lapply(1:ncol(metabolites_Adenomas_affected), function(x) c(MiRKAT_affected[[x]]$p_values, MiRKAT_affected[[x]]$R2))))
colnames(df_results_MiRKAT_affected) = c("pvalues", "R2")
df_results_MiRKAT_affected$Metabo = colnames(metabolites_Adenomas_affected)
df_results_MiRKAT_unaffected = data.frame(do.call("rbind", lapply(1:ncol(metabolites_Adenomas_unaffected), function(x) c(MiRKAT_unaffected[[x]]$p_values, MiRKAT_unaffected[[x]]$R2))))
colnames(df_results_MiRKAT_unaffected) = c("pvalues", "R2")
df_results_MiRKAT_unaffected$Metabo = colnames(metabolites_Adenomas_unaffected)


signif_metabolites_unaffected = df_results_MiRKAT_unaffected$Metabo[which(df_results_MiRKAT_unaffected$pvalues<=(0.05/nrow(df_results_MiRKAT_unaffected)))]
signif_metabolites_affected = df_results_MiRKAT_affected$Metabo[which(df_results_MiRKAT_affected$pvalues<=(0.05/nrow(df_results_MiRKAT_affected)))]

signif_metabolites_unaffected[signif_metabolites_unaffected%in%signif_metabolites_affected]

#RDA

clr_species_Adenomas = species_Adenomas %>% clr() %>% apply_nzv_filter() 
clr_species_Adenomas_affected = species_Adenomas_affected %>% clr() %>% apply_nzv_filter()
clr_species_Adenomas_unaffected = species_Adenomas_unaffected %>% clr() %>% apply_nzv_filter()


spe.list = list( "affected" = clr_species_Adenomas_affected %>% apply_var_coef(), "unaffected" = clr_species_Adenomas_unaffected %>%  apply_var_coef())
met.list = list( "affected" = metabolites_Adenomas_affected %>% as.data.frame() %>% apply_var_coef(), "unaffected" = metabolites_Adenomas_unaffected %>% as.data.frame() %>% apply_var_coef())

spe.rda = purrr::map2(spe.list, met.list, ~rda(.x ~ ., data = .y))


var_prop.rda = purrr::map(spe.rda, ~ summary(eigenvals(.x)) %>% t %>% as.data.frame() %>% tibble::rownames_to_column("RDA") %>%
                             dplyr::filter(stringr::str_detect(RDA, "RDA")) %>%
                             dplyr::select("Cumulative Proportion") %>%
                             dplyr::mutate(index = 1:nrow(.)))

var_prop_plot.rda = purrr::imap(var_prop.rda, ~ ggplot(.x, aes(x = index, y = `Cumulative Proportion`)) +
                                   geom_line() +
                                   geom_point() +
                                   theme_minimal() +
                                   xlab("Number of components") + ylab("Cummulative % of explained var.") +
                                   coord_cartesian(ylim = c(0,1)) +
                                   scale_x_continuous(breaks = c(1,seq(10,170, by = 10)), minor_breaks = NULL) +
                                   scale_y_continuous(minor_breaks = NULL) +
                                   ggtitle(.y))

var_prop_plot.rda$affected + var_prop_plot.rda$unaffected

#PLS 
pls.reg.res = purrr::map2(spe.list, met.list, ~{
  mixOmics::pls(X = .x, Y = .y, ncomp = 2, mode = 'regression')
})

po = purrr::imap(pls.reg.res, ~plot_loadings(.x,  titre = .y))
 po$affected / po$unaffected

#clr-lm 

results_clr_lm_affected = list()
results_clr_lm_unaffected = list()

for(c0 in signif_metabolites_unaffected){
  print(c0)
  for(c1 in 1:ncol(clr_species_Adenomas_unaffected)){
    results_clr_lm_unaffected[[c0]][[c1]] = summary(lm(as.data.frame(metabolites_Adenomas_unaffected)[,c0]~as.data.frame(clr_species_Adenomas_unaffected)[,c1]))$coefficients
  }
  
}

for(c0 in signif_metabolites_affected){
  print(c0)
  for(c1 in 1:ncol(clr_species_Adenomas_affected)){
    results_clr_lm_affected[[c0]][[c1]] = summary(lm(as.data.frame(metabolites_Adenomas_affected)[,c0]~as.data.frame(clr_species_Adenomas_affected)[,c1]))$coefficients
  }
  
}

df_results_ex_affected = data.frame(t(sapply(results_clr_lm_affected$`N-acetyl-cadaverine`, function(x) x[2,c(1,4)])))
colnames(df_results_ex_affected) = c("Estimate", "P")
df_results_ex_affected$signi = df_results_ex_affected$P<=0.05/nrow(df_results_ex_affected)

df_results_ex_unaffected = data.frame(t(sapply(results_clr_lm_unaffected$`N-acetyl-cadaverine`, function(x) x[2,c(1,4)])))
colnames(df_results_ex_unaffected) = c("Estimate", "P")
df_results_ex_unaffected$signi = df_results_ex_unaffected$P<=0.05/nrow(df_results_ex_unaffected)

df_results_ex_affected$Species = ifelse(df_results_ex_affected$signi, colnames(clr_species_Adenomas_affected), NA)
df_results_ex_unaffected$Species = colnames(clr_species_Adenomas_unaffected)
df_results_ex_unaffected$Species = ifelse(df_results_ex_unaffected$signi, colnames(clr_species_Adenomas_unaffected), NA)

a = ggplot(df_results_ex_affected, aes(x=Estimate, y=-log10(P), colour=signi, label=Species))+geom_point()+geom_hline(yintercept = -log10(0.05/nrow(df_results_ex_affected)))+ggtitle("Affected")+theme_bw()
b = ggplot(df_results_ex_unaffected, aes(x=Estimate, y=-log10(P), colour=signi, label=Species))+geom_point()+geom_hline(yintercept = -log10(0.05/nrow(df_results_ex_unaffected)))+ggtitle("Unaffected")+theme_bw()

a+b


lapply(1:length(results_clr_lm_unaffected), function(y) results_clr_lm_unaffected[[y]][sapply(results_clr_lm_unaffected[[y]], function(x) x[2,4])<= (0.05/179)])
lapply(1:length(results_clr_lm_unaffected), function(y) colnames(species_Adenomas_unaffected)[which(sapply(results_clr_lm_unaffected[[y]], function(x) x[2,4])<= (0.05/179))])

#CODA-LASSO

species_Adenomas_affected_filtered = species_Adenomas_affected %>% as.data.frame() %>% apply_nzv_filter()
species_Adenomas_unaffected_filtered = species_Adenomas_unaffected %>% as.data.frame() %>% apply_nzv_filter()

results_CODA_affected = list()

for(col in signif_metabolites_affected){
  print(col)
  results_CODA_affected[[col]] = coda_glmnet(species_Adenomas_affected_filtered , metabolites_Adenomas_affected %>% select(all_of(col)) %>% unlist(), showPlots = F)
}


results_CODA_unaffected = list()

for(col in signif_metabolites_unaffected){
  print(col)
  results_CODA_unaffected[[col]] = coda_glmnet(species_Adenomas_unaffected_filtered, metabolites_Adenomas_unaffected %>% select(all_of(col)) %>% unlist(), showPlots = F)
}

(networks_from_CODA_LASSO(results_CODA_affected)+ggtitle("Affected")) + (networks_from_CODA_LASSO(results_CODA_unaffected)+ggtitle("Unaffected"))

#sPLS-Reg

spe.list = list( "affected" = clr_species_Adenomas_affected , "unaffected" = clr_species_Adenomas_unaffected )
met.list = list( "affected" = metabolites_Adenomas_affected, "unaffected" = metabolites_Adenomas_unaffected)

spls.tune = purrr::map2(spe.list, met.list, ~{
  
  set.seed(123)
  try(mixOmics::tune.spls(X = .x, Y = .y, ncomp = 2,
                          test.keepX = round(seq(from = ncol(.x)*0.02, to = ncol(.x)*0.2, length.out = 5)),
                          test.keepY = round(seq(from = ncol(.y)*0.02, to = ncol(.y)*0.2, length.out = 5)),
                          nrepeat = 1, folds = 10, # use 10 folds
                          mode = 'regression', measure="cor"))
})


spls.reg.final <- purrr::imap(pls.tune, ~mixOmics::spls(X = spe.list[[.y]],
                                                         Y = met.list[[.y]],
                                                         ncomp = 2,
                                                         keepX = spls.tune[[.y]]$choice.keepX,
                                                         keepY = spls.tune[[.y]]$choice.keepY,
                                                         mode = "regression"))

# circle correlation plot
purrr::imap(spls.reg.final, ~mixOmics::plotVar(.x, var.names = F, pch = c(19,19), title = .y))
