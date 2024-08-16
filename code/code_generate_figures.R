library(ggplot2)
library(ggpubr)
library(wesanderson)
library(igraph)

gg_qqplot <- function(ps, ci = 0.95, colour="blue", title="") {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3, colour=colour) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)+theme_bw()+ggtitle(title)#+annotate(geom="text", x=x,y=y,label=text)
}

gg_qqplot_facet_grid <- function(list_pvalues,ci = 0.95, title="") {
  #n  <- length(ps)
  
  dfs <- lapply(1:length(list_pvalues), function(i) {
    n = length(list_pvalues[[i]])
    data.frame(
      Method = names(list_pvalues)[i],
      observed = -log10(sort(list_pvalues[[i]])),
      expected = -log10(ppoints(n)),
      clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
    )})
  
  all.dfs = do.call("rbind", dfs)
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  ggplot(all.dfs) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed, color=Method), size = 2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5)+xlab(log10Pe) +
    ylab(log10Po)+theme_bw()+scale_color_npg()#+theme(legend.position = "none")
  
}


mantel.ilr.spearman.100.100.500 = list.files("data\\Mantel\\ilr\\100_100_500\\spearman", full.names = TRUE)
mantel.ilr.spearman.25.25.100 = list.files("data\\Mantel\\ilr\\25_25_100\\spearman", full.names = TRUE)


pvalues.null.mantel.alpha.log.euclidean.100.100.500 = readRDS(mantel.alpha.spearman.100.100.500[8])

pvalues.null.mantel.ilr.log.canberra.100.100.500 = readRDS(mantel.ilr.spearman.100.100.500[7])
pvalues.null.mantel.ilr.log.euclidean.100.100.500 = readRDS(mantel.ilr.spearman.100.100.500[8])
pvalues.null.mantel.ilr.log.manhattan.100.100.500 = readRDS(mantel.ilr.spearman.100.100.500[9])


pvalues.alternative.mantel.ilr.log.canberra.100.100.500 = readRDS(mantel.ilr.spearman.100.100.500[1])
pvalues.alternative.mantel.ilr.log.euclidean.100.100.500 = readRDS(mantel.ilr.spearman.100.100.500[2])
pvalues.alternative.mantel.ilr.log.manhattan.100.100.500 = readRDS(mantel.ilr.spearman.100.100.500[3])



df.power.Mantel.ilr.log.100.100.500 = data.frame("Power" = c(sum(pvalues.alternative.mantel.ilr.log.euclidean.100.100.500<=0.05)/1000,
                                                         sum(pvalues.alternative.mantel.ilr.log.canberra.100.100.500<= 0.05)/1000,
                                                         sum(pvalues.alternative.mantel.ilr.log.manhattan.100.100.500<= 0.05)/1000), "Distance.Kernel"= c("Euclidean","Canberra","Manhattan"))


pvalues.null.mantel.ilr.log.euclidean.100.100.500 = readRDS(mantel.ilr.spearman.100.100.500[8])

pvalues.null.MMiRKAT.100.100.500 = readRDS("data\\MMiRKAT\\list.results.null.MMiRKAT.compositional.100.100.500.RDS")
pvalues.alter.MMiRKAT.100.100.500 = readRDS("data\\MMiRKAT\\list.results.power.MMiRKAT.compositional.100.100.500.RDS")

pvalues.null.MMiRKAT.25.25.100 = readRDS("data\\MMiRKAT\\list.results.null.MMiRKAT.compositional.25.25.100.RDS")
pvalues.alter.MMiRKAT.25.25.100 = readRDS("data\\MMiRKAT\\list.results.power.MMiRKAT.compositional.25.25.100.RDS")




df.power.MMiRKAT.ilr.log.100.100.500 = data.frame("Power" = c(sum(pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.ilr.log$power.p.values.100.100.500indiv.MMiRKAT.Canberra.ilr.log<=0.05)/1000,
                                                              sum(pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.ilr.log$power.p.values.100.100.500indiv.MMiRKAT.Euclidean.ilr.log<=0.05)/1000,
                                                              sum(pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.ilr.log$power.p.values.100.100.500indiv.MMiRKAT.Manhattan.ilr.log<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))


MMiRKAT.alpha.100.100.500 = list.files("data\\MMiRKAT\\alpha\\100_100_500", full.names=TRUE)

pvalues.null.MMiRKAT.log.Canberra.100.100.500 = readRDS(MMiRKAT.alpha.100.100.500[1])
pvalues.null.MMiRKAT.log.Euclidean.100.100.500 = readRDS(MMiRKAT.alpha.100.100.500[3])
pvalues.null.MMiRKAT.log.Manhattan.100.100.500 = readRDS(MMiRKAT.alpha.100.100.500[5])

pvalues.alter.MMiRKAT.log.Canberra.100.100.500 = readRDS(MMiRKAT.alpha.100.100.500[7])
pvalues.alter.MMiRKAT.log.Euclidean.100.100.500 = readRDS(MMiRKAT.alpha.100.100.500[9])
pvalues.alter.MMiRKAT.log.Manhattan.100.100.500 = readRDS(MMiRKAT.alpha.100.100.500[11])

df.power.MMiRKAT.alpha.log.100.100.500 = data.frame("Power"=c(sum(pvalues.alter.MMiRKAT.log.Canberra.100.100.500<=0.05)/1000,
                                                              sum(pvalues.alter.MMiRKAT.log.Euclidean.100.100.500<=0.05)/1000,
                                                              sum(pvalues.alter.MMiRKAT.log.Manhattan.100.100.500<=0.05)/1000), 
                                                    "Distance-Kernel"=c("Canberra", "Euclidean", "Manhattan"))


ggplot(df.power.MMiRKAT.alpha.log.100.100.500, aes(x=Distance.Kernel, y=Power, fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", colour="black")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("#DF8F44FF","#00A1D5FF" ,"#B24745FF"))+ xlab("Distance Kernel")
df.power.alpha.log.MMiRKAT.Mantel.euclidean.100.100.500 = data.frame("Power" = c(sum(pvalues.alter.mantel.alpha.log.euclidean.100.100.500<=0.05)/1000,
                                                                     sum(pvalues.alter.MMiRKAT.log.Euclidean.100.100.500<=0.05)/1000), "Method" = c("Mantel", "MMiRKAT") )

ggplot(df.power.alpha.log.MMiRKAT.Mantel.euclidean.100.100.500, aes(x=Method, y=Power, fill=Method)) + geom_bar(position="dodge", stat="identity", colour="black")+
  theme_bw()+ylab("Power")+theme(legend.position = "none")+scale_fill_manual(values=wes_palette(n=2, name="FrenchDispatch"))



#Multivariate factorization-based methods

CCA.alpha.log.100.100.500 = readRDS("data\\results_CCA\\alpha\\100_100_500\\redundancy_cca_alpha_microbiome_log_metabolome_100_100_500indiv.RDS")
PLS.Can.alpha.log.100.100.500 = readRDS("data\\results_PLS\\alpha\\100_100_500\\canonical\\redundancy_pls_alpha_microbiome_log_metabolome_100_100_500indiv_canonical.RDS")
PLS.Reg.alpha.log.100.100.500 = readRDS("data\\results_PLS\\alpha\\100_100_500\\regression\\redundancy_pls_alpha_microbiome_log_metabolome_100_100_500indiv_regression.RDS")
RDA.alpha.log.100.100.500 = readRDS("data\\results_RDA\\alpha\\100_100_500\\explained_variance_alpha_microbiome_log_metabolome_rda_100_100_500indiv.RDS")
MOFA.alpha.log.100.100.500 = readRDS("data\\results_MOFA2\\alpha\\100_100_500\\MOFA2_alpha_microbiome_log_metabolome_Scenario_100_100_500indiv.RDS")

df.results.alpha.log.multivariate.factorization.based.100.100.500 = data.frame("Explained-Variance"=c(CCA.alpha.log.100.100.500*100,
                                                                               PLS.Reg.alpha.log.100.100.500*100,
                                                                               PLS.Can.alpha.log.100.100.500*100,
                                                                               RDA.alpha.log.100.100.500*100,
                                                                               sapply(1:1000, function(X) MOFA.alpha.log.100.100.500[[X]]$r2_total$group1[1])),
                                                                               "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))

plot.redundancy.alpha = ggplot(df.results.alpha.log.multivariate.factorization.based.100.100.500, aes(x=Method, y=Explained.Variance))+geom_boxplot(fill="grey")+
  theme_bw()+ylab("% Explained Variance")+theme(legend.position = "none")+ylim(c(0,100))

CCA.ilr.log.100.100.500 = readRDS("data\\results_CCA\\ilr\\100_100_500\\redundancy_cca_ilr_microbiome_log_metabolome_100_100_500indiv.RDS")
PLS.Can.ilr.log.100.100.500 = readRDS("data\\results_PLS\\ilr\\100_100_500\\canonical\\redundancy_pls_ilr_microbiome_log_metabolome_100_100_500indiv_canonical.RDS")
PLS.Reg.ilr.log.100.100.500 = readRDS("data\\results_PLS\\ilr\\100_100_500\\regression\\redundancy_pls_ilr_microbiome_log_metabolome_100_100_500indiv_regression.RDS")
RDA.ilr.log.100.100.500 = readRDS("data\\results_RDA\\ilr\\100_100_500\\explained_variance_ilr_microbiome_log_metabolome_rda_100_100_500indiv.RDS")
MOFA.ilr.log.100.100.500 = readRDS("data\\results_MOFA2\\ilr\\100_100_500\\MOFA2_ilr_microbiome_log_metabolome_Scenario_100_100_500indiv.RDS")

df.results.ilr.log.multivariate.factorization.based.100.100.500 = data.frame("Explained-Variance"=c(CCA.ilr.log.100.100.500*100,
                                                                                                      PLS.Reg.ilr.log.100.100.500*100,
                                                                                                      PLS.Can.ilr.log.100.100.500*100,
                                                                                                      RDA.ilr.log.100.100.500*100,
                                                                                                      sapply(1:1000, function(X) MOFA.ilr.log.100.100.500[[X]]$r2_total$group1[1])),
                                                                               "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))




plot.redundancy.ilr = ggplot(df.results.ilr.log.multivariate.factorization.based.100.100.500, aes(x=Method, y=Explained.Variance))+geom_boxplot(fill="grey")+theme_bw()+ylab("% Explained Variance")+theme(legend.position = "none")+ylim(c(0,100))



CCA.clr.log.100.100.500 = readRDS("data\\results_CCA\\clr\\100_100_500\\redundancy_cca_clr_microbiome_log_metabolome_100_100_500indiv.RDS")
PLS.Can.clr.log.100.100.500 = readRDS("data\\results_PLS\\clr\\100_100_500\\canonical\\redundancy_pls_clr_microbiome_log_metabolome_100_100_500indiv_canonical.RDS")
PLS.Reg.clr.log.100.100.500 = readRDS("data\\results_PLS\\clr\\100_100_500\\regression\\redundancy_pls_clr_microbiome_log_metabolome_100_100_500indiv_regression.RDS")
RDA.clr.log.100.100.500 = readRDS("data\\results_RDA\\clr\\100_100_500\\explained_variance_clr_microbiome_log_metabolome_rda_100_100_500indiv.RDS")
MOFA.clr.log.100.100.500 = readRDS("data\\results_MOFA2\\clr\\100_100_500\\MOFA2_clr_microbiome_log_metabolome_Scenario_100_100_500indiv.RDS")


df.results.clr.log.multivariate.factorization.based.100.100.500 = data.frame("Explained-Variance"=c(sapply(CCA.clr.log.100.100.500, as.numeric)*100,
                                                                                                  PLS.Reg.clr.log.100.100.500*100,
                                                                                                  PLS.Can.clr.log.100.100.500*100,
                                                                                                  RDA.clr.log.100.100.500*100,
                                                                                                  sapply(1:1000, function(X) MOFA.clr.log.100.100.500[[X]]$r2_total$group1[1])),
                                                                           "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA"), each=1000))

plot.redundancy.clr.100.100.500 = ggplot(df.results.clr.log.multivariate.factorization.based.100.100.500, aes(x=Method, y=Explained.Variance, fill=Method))+geom_boxplot()+theme_bw()+scale_fill_manual(values=wes_palette(n=5, name="Royal2"))+ylab("% Explained Variance")+theme(legend.position = "none")+ylim(c(0,100))



df.power.Mantel.ilr.log.100.100.500$Method = "Mantel Test"
df.power.MMiRKAT.ilr.log.100.100.500$Method = "MMiRKAT"

df.power.GA = rbind(df.power.Mantel.ilr.log.100.100.500,df.power.MMiRKAT.ilr.log.100.100.500)
df.results.alpha.log.multivariate.factorization.based.100.100.500$Normalization = "Alpha"
df.results.ilr.log.multivariate.factorization.based.100.100.500$Normalization = "ILR"
df.results.multivariate.factorization.based.100.100.500 = rbind(df.results.alpha.log.multivariate.factorization.based.100.100.500,df.results.ilr.log.multivariate.factorization.based.100.100.500)


ggarrange(ggarrange(ggarrange(gg_qqplot_facet_grid(list("Canberra"=pvalues.null.mantel.ilr.log.canberra.100.100.500,"Euclidean"=pvalues.null.mantel.ilr.log.euclidean.100.100.500,"Manhattan"=pvalues.null.mantel.ilr.log.manhattan.100.100.500))+ggtitle("Mantel Test")+theme(legend.position = c(0.35,0.75))+labs(color="Distance Kernel"),
          gg_qqplot_facet_grid(list("Canberra"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.ilr.log$null.p.values.100.100.500indiv.MMiRKAT.Canberra.ilr.log,"Euclidean"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.ilr.log$null.p.values.100.100.500indiv.MMiRKAT.Euclidean.ilr.log,"Manhattan"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.ilr.log$null.p.values.100.100.500indiv.MMiRKAT.Manhattan.ilr.log))+ggtitle("MMiRKAT")+theme(legend.position = c(0.35,0.75))+labs(color="Distance Kernel"),
          labels = c("A","B")),
ggarrange(ggplot(df.power.GA, aes(x=Distance.Kernel,y=Power, fill=Distance.Kernel)) + geom_bar(position="dodge", stat = "identity",colour="black", width=0.7)+xlab("Distance Kernel")+theme_bw()+ theme(legend.position = "none")+ylim(c(0,1))+facet_grid(.~Method)+scale_fill_npg(), labels = "C"),ncol=2),
ggarrange(ggplot(df.results.multivariate.factorization.based.100.100.500, aes(x=Method, y=Explained.Variance, fill=Method))+geom_violin()+scale_fill_nejm()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+ylim(c(0,100)), labels="D"),nrow=2)

#Multivariate penalized factorization-based methods

simulated.data.100.100.500.alter = readRDS("data\\data_alternative_100_100_500indiv.RDS")
simulated.data.25.25.100.alter = readRDS("data\\data_alternative_25_25_100indiv.RDS")


sCCA.clr.log.25.25.100 = readRDS("data\\results_sCCA\\clr\\25_25_100\\sCCA_clrMicrobiome_logMetabolome_25_25_100ind.RDS")
sPLS.canonical.clr.log.25.25.100 = readRDS("data\\results_sPLS\\clr\\25_25_100\\canonical\\loadings_sPLScanonical_clr_Microbiome_log_Metabolome_25_25_100indiv.RDS")
sPLS.regression.clr.log.25.25.100 = readRDS("data\\results_sPLS\\clr\\25_25_100\\regression\\loadings_sPLSregression_clr_Microbiome_log_Metabolome_25_25_100indiv.RDS")

sCCA.clr.log.100.100.500 = readRDS("data\\results_sCCA\\clr\\100_100_500\\sCCA_clrMicrobiome_logMetabolome_100_100_500ind.RDS")
sPLS.canonical.clr.log.100.100.500 = readRDS("data\\results_sPLS\\clr\\100_100_500\\canonical\\loadings_sPLScanonical_clr_Microbiome_log_Metabolome_100_100_500indiv.RDS")
sPLS.regression.clr.log.100.100.500 = readRDS("data\\results_sPLS\\clr\\100_100_500\\regression\\loadings_sPLSregression_clr_Microbiome_log_Metabolome_100_100_500indiv.RDS")


sCCA.original.log.25.25.100 = readRDS("data\\results_sCCA\\original\\25_25_100\\sCCA_noneMicrobiome_logMetabolome_25_25_100ind.RDS")
sPLS.canonical.original.log.25.25.100 = readRDS("data\\results_sPLS\\original\\25_25_100\\canonical\\loadings_sPLScanonical_Microbiome_log_Metabolome_25_25_100indiv.RDS")
sPLS.regression.original.log.25.25.100 = readRDS("data\\results_sPLS\\original\\25_25_100\\regression\\loadings_sPLSregression_Microbiome_log_Metabolome_25_25_100indiv.RDS")

sCCA.original.log.100.100.500 = readRDS("data\\results_sCCA\\original\\100_100_500\\sCCA_noneMicrobiome_logMetabolome_100_100_500ind.RDS")
sPLS.canonical.original.log.100.100.500 = readRDS("data\\results_sPLS\\original\\100_100_500\\canonical\\loadings_sPLScanonical_Microbiome_log_Metabolome_100_100_500indiv.RDS")
sPLS.regression.original.log.100.100.500 = readRDS("data\\results_sPLS\\original\\100_100_500\\regression\\loadings_sPLSregression_Microbiome_log_Metabolome_100_100_500indiv.RDS")

sCCA.original.original.25.25.100 = readRDS("data\\results_sCCA\\original\\25_25_100\\sCCA_noneMicrobiome_noneMetabolome_25_25_100ind.RDS")
sPLS.canonical.original.original.25.25.100 = readRDS("data\\results_sPLS\\original\\25_25_100\\canonical\\loadings_sPLScanonical_Microbiome_Metabolome_25_25_100indiv.RDS")
sPLS.regression.original.original.25.25.100 = readRDS("data\\results_sPLS\\original\\25_25_100\\regression\\loadings_sPLSregression_Microbiome_Metabolome_25_25_100indiv.RDS")

sCCA.original.original.100.100.500 = readRDS("data\\results_sCCA\\original\\100_100_500\\sCCA_noneMicrobiome_noneMetabolome_100_100_500ind.RDS")
sPLS.canonical.original.original.100.100.500 = readRDS("data\\results_sPLS\\original\\100_100_500\\canonical\\loadings_sPLScanonical_Microbiome_Metabolome_100_100_500indiv.RDS")
sPLS.regression.original.original.100.100.500 = readRDS("data\\results_sPLS\\original\\100_100_500\\regression\\loadings_sPLSregression_Microbiome_Metabolome_100_100_500indiv.RDS")


index.nonnull.sPLS.regression.clr.log.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.regression.clr.log.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.regression.clr.log.100.100.500[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.clr.log.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.canonical.clr.log.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.canonical.clr.log.100.100.500[[X]]$Y)!=0))
})

index.nonnull.sPLS.regression.clr.log.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.regression.clr.log.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.regression.clr.log.25.25.100[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.clr.log.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.canonical.clr.log.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.canonical.clr.log.25.25.100[[X]]$Y)!=0))
})


index.nonnull.sPLS.regression.original.log.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.regression.original.log.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.regression.original.log.100.100.500[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.original.log.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.canonical.original.log.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.canonical.original.log.100.100.500[[X]]$Y)!=0))
})

index.nonnull.sPLS.regression.original.original.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.regression.original.original.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.regression.original.original.100.100.500[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.original.original.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.canonical.original.original.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.canonical.original.original.100.100.500[[X]]$Y)!=0))
})

index.nonnull.sPLS.regression.original.log.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.regression.original.log.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.regression.original.log.25.25.100[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.original.log.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.canonical.original.log.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.canonical.original.log.25.25.100[[X]]$Y)!=0))
})


index.nonnull.sPLS.regression.original.original.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.regression.original.original.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.regression.original.original.25.25.100[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.original.original.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.canonical.original.original.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.canonical.original.original.25.25.100[[X]]$Y)!=0))
})

confusion.matrix = function(index.associations.Microbiome,index.non.associations.Microbiome,index.associations.Metabolome, index.non.associations.Metabolome,true.associated.Microbiome, true.non.associated.Microbiome,
                            true.associated.Metabolome, true.non.associated.Metabolome){
  
  TP = sum(index.associations.Microbiome%in%true.associated.Microbiome) + sum(index.associations.Metabolome%in%true.associated.Metabolome)
  TN = sum(index.non.associations.Microbiome%in%true.non.associated.Microbiome) + sum(index.non.associations.Metabolome%in%true.non.associated.Metabolome)
  FP = sum(index.associations.Microbiome%in%true.non.associated.Microbiome) + sum(index.associations.Metabolome%in%true.non.associated.Metabolome)
  FN = sum(index.non.associations.Microbiome%in%true.associated.Microbiome) + sum(index.non.associations.Metabolome%in%true.associated.Metabolome)
  
  return(matrix(c(TP, FN,FP, TN), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg"))))
}

sensivity.specificity.from.confusion.matrix = function(confusion.matrix){
  
  sensitivity = confusion.matrix[1,1]/(confusion.matrix[1,1]+confusion.matrix[2,1])
  specificity = confusion.matrix[2,2]/(confusion.matrix[2,2]+confusion.matrix[1,2])
  
  return(c(sensitivity,specificity))
}

f1.score = function(confusion.matrix){
  m = sensivity.specificity.from.confusion.matrix(confusion.matrix)
  
  return(2/((1/m[1])+(1/m[2])))
}    


f1.Score.sCCA.clr.log.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.clr.log.100.100.500[[X]]$index$Microbiome,c(1:100)[-sCCA.clr.log.100.100.500[[X]]$index$Microbiome],
                 sCCA.clr.log.100.100.500[[X]]$index$Metabolome,c(1:100)[-sCCA.clr.log.100.100.500[[X]]$index$Metabolome],
                 c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                 c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                 c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                 c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.clr.log.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.clr.log.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.clr.log.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.clr.log.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.clr.log.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})



f1.Score.sPLS.canonical.clr.log.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.clr.log.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.clr.log.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.clr.log.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.clr.log.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})


F1Score.all.methods.clr.log.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.clr.log.100.100.500,f1.Score.sPLS.regression.clr.log.100.100.500,f1.Score.sPLS.canonical.clr.log.100.100.500), "Scenario"="High-Dimensional", "Metric"="F1-Score")

sparsity.all.methods.clr.log.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.clr.log.100.100.500[[X]]$index$Microbiome) +length(sCCA.clr.log.100.100.500[[X]]$index$Metabolome))/200 )
,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.clr.log.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.clr.log.100.100.500[[X]]$Metabolome))/200)
,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.clr.log.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.clr.log.100.100.500[[X]]$Metabolome))/200)), "Scenario" = "High-Dimensional", "Metric" = "Sparsity")

df.results.clr.log.multi.penalized.100.100.500 = rbind(F1Score.all.methods.clr.log.100.100.500,sparsity.all.methods.clr.log.100.100.500)
df.results.clr.log.multi.penalized.100.100.500$Metric = factor(df.results.clr.log.multi.penalized.100.100.500$Metric, levels = c("Sparsity", "F1-Score"))



f1.Score.sCCA.clr.log.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.clr.log.25.25.100[[X]]$index$Microbiome,c(1:100)[-sCCA.clr.log.25.25.100[[X]]$index$Microbiome],
                            sCCA.clr.log.25.25.100[[X]]$index$Metabolome,c(1:100)[-sCCA.clr.log.25.25.100[[X]]$index$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.clr.log.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.clr.log.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.clr.log.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.clr.log.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.clr.log.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.canonical.clr.log.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.clr.log.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.clr.log.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.clr.log.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.clr.log.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})

F1Score.all.methods.clr.log.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.clr.log.25.25.100,f1.Score.sPLS.regression.clr.log.25.25.100,f1.Score.sPLS.canonical.clr.log.25.25.100), "Scenario"="Low-Dimensional", "Metric"="F1-Score")

sparsity.all.methods.clr.log.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.clr.log.25.25.100[[X]]$index$Microbiome) +length(sCCA.clr.log.25.25.100[[X]]$index$Metabolome))/200 )
                                                                                                                               ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.clr.log.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.clr.log.25.25.100[[X]]$Metabolome))/200)
                                                                                                                               ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.clr.log.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.clr.log.25.25.100[[X]]$Metabolome))/200)), "Scenario" = "Low-Dimensional", "Metric" = "Sparsity")

df.results.clr.log.multi.penalized.25.25.100 = rbind(F1Score.all.methods.clr.log.25.25.100,sparsity.all.methods.clr.log.25.25.100)


df.results.clr.log.multi.penalized.25.25.100$Metric = factor(df.results.clr.log.multi.penalized.25.25.100$Metric, levels=c("Sparsity","F1-Score"))

plot.results.penalized.multivariate = ggarrange(ggplot(df.results.clr.log.multi.penalized.25.25.100, aes(x=Method, y=Value*100, fill=Method))+geom_violin()+theme_bw()+facet_grid(~Metric)+scale_fill_futurama()+ylab("%")+theme(legend.position = "none")+ylim(c(0,100)),
          ggplot(df.results.clr.log.multi.penalized.100.100.500, aes(x=Method, y=Value*100, fill=Method))+geom_violin()+theme_bw()+facet_grid(~Metric)+scale_fill_futurama()+ylab("%")+theme(legend.position = "none")+ylim(c(0,100)),labels = c("C","F"), ncol=2)



f1.Score.sCCA.original.log.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.original.log.100.100.500[[X]]$index$Microbiome,c(1:100)[-sCCA.original.log.100.100.500[[X]]$index$Microbiome],
                            sCCA.original.log.100.100.500[[X]]$index$Metabolome,c(1:100)[-sCCA.original.log.100.100.500[[X]]$index$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.original.log.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.original.log.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.original.log.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.original.log.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.original.log.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})


f1.Score.sPLS.canonical.original.log.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.original.log.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.original.log.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.original.log.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.original.log.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sCCA.original.original.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.original.log.100.100.500[[X]]$index$Microbiome,c(1:100)[-sCCA.original.original.100.100.500[[X]]$index$Microbiome],
                            sCCA.original.log.100.100.500[[X]]$index$Metabolome,c(1:100)[-sCCA.original.original.100.100.500[[X]]$index$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.original.original.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.original.original.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.original.original.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.original.original.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.original.original.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})


f1.Score.sPLS.canonical.original.original.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.original.original.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.original.original.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.original.original.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.original.original.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})


F1Score.all.methods.original.log.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.original.log.100.100.500,f1.Score.sPLS.regression.original.log.100.100.500,f1.Score.sPLS.canonical.original.log.100.100.500), "Scenario" = "High-Dimensional", "Metric" = "F1-Score")

sparsity.all.methods.original.log.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.original.log.100.100.500[[X]]$index$Microbiome) +length(sCCA.original.log.100.100.500[[X]]$index$Metabolome))/200 )
                                                                                                                                  ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.original.log.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.original.log.100.100.500[[X]]$Metabolome))/200)
                                                                                                                                  ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.original.log.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.original.log.100.100.500[[X]]$Metabolome))/200)), "Scenario"="High-Dimensional", "Metric"="Sparsity")
df.results.original.log.multi.penalized.100.100.500 = rbind(F1Score.all.methods.original.log.100.100.500,sparsity.all.methods.original.log.100.100.500)

plot.F1Score.original.log.100.100.500 = ggplot(F1Score.all.methods.original.log.100.100.500, aes(x=Method,y=F1.Score, fill=Method))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("F1-Score")+theme(legend.position = "none")+ylim(c(0,1))

plot.sparsity.original.log.100.100.500 = ggplot(sparsity.all.methods.original.log.100.100.500, aes(x=Method,y=Sparsity, fill=Method))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("% Selected Features")+theme(legend.position = "none")+ylim(c(0,1))


f1.Score.sCCA.original.log.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.original.log.25.25.100[[X]]$index$Microbiome,c(1:100)[-sCCA.original.log.25.25.100[[X]]$index$Microbiome],
                            sCCA.original.log.25.25.100[[X]]$index$Metabolome,c(1:100)[-sCCA.original.log.25.25.100[[X]]$index$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.original.log.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.original.log.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.original.log.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.original.log.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.original.log.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})


f1.Score.sPLS.canonical.original.log.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.original.log.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.original.log.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.original.log.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.original.log.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})


f1.Score.sCCA.original.original.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.original.original.25.25.100[[X]]$index$Microbiome,c(1:100)[-sCCA.original.original.25.25.100[[X]]$index$Microbiome],
                            sCCA.original.original.25.25.100[[X]]$index$Metabolome,c(1:100)[-sCCA.original.original.25.25.100[[X]]$index$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.original.original.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.original.original.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.original.original.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.original.original.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.original.original.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})


f1.Score.sPLS.canonical.original.original.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.original.original.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.original.original.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.original.original.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.original.original.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})

sapply(1:1000, function(X) (length(sCCA.original.log.25.25.100[[X]]$index$Microbiome) +length(sCCA.original.log.25.25.100[[X]]$index$Metabolome))/50 )
sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.original.log.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.original.log.25.25.100[[X]]$Metabolome))/50)
sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.original.log.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.original.log.25.25.100[[X]]$Metabolome))/50)


F1Score.all.methods.original.log.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.original.log.25.25.100,f1.Score.sPLS.regression.original.log.25.25.100,f1.Score.sPLS.canonical.original.log.25.25.100), "Scenario" = "Low-Dimensional", "Metric" = "F1-Score")

sparsity.all.methods.original.log.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.original.log.25.25.100[[X]]$index$Microbiome) +length(sCCA.original.log.25.25.100[[X]]$index$Metabolome))/50 )
                                                                                                                                    ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.original.log.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.original.log.25.25.100[[X]]$Metabolome))/50)
                                                                                                                                    ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.original.log.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.original.log.25.25.100[[X]]$Metabolome))/50)), "Scenario"="Low-Dimensional", "Metric"="Sparsity")
df.results.original.log.multi.penalized.25.25.100 = rbind(F1Score.all.methods.original.log.25.25.100,sparsity.all.methods.original.log.25.25.100)

ggarrange(ggplot(df.results.original.log.multi.penalized.25.25.100, aes(x=Method,y=Value,fill=Method))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("%")+theme(legend.position = "none")+ylim(c(0,1))+facet_grid(~Metric),
  ggplot(df.results.original.log.multi.penalized.100.100.500, aes(x=Method,y=Value,fill=Method))+geom_boxplot()+theme_bw()+
    scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("%")+theme(legend.position = "none")+ylim(c(0,1))+facet_grid(~Metric), labels=c("A","B"), ncol=2)


F1Score.all.methods.original.original.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.original.original.100.100.500,f1.Score.sPLS.regression.original.original.100.100.500,f1.Score.sPLS.canonical.original.original.100.100.500), "Scenario" = "High-Dimensional", "Metric" = "F1-Score")

sparsity.all.methods.original.original.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.original.original.100.100.500[[X]]$index$Microbiome) +length(sCCA.original.original.100.100.500[[X]]$index$Metabolome))/200 )
                                                                                                                                    ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.original.original.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.original.original.100.100.500[[X]]$Metabolome))/200)
                                                                                                                                    ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.original.original.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.original.original.100.100.500[[X]]$Metabolome))/200)), "Scenario"="High-Dimensional", "Metric"="Sparsity")
df.results.original.original.multi.penalized.100.100.500 = rbind(F1Score.all.methods.original.original.100.100.500,sparsity.all.methods.original.original.100.100.500)


F1Score.all.methods.original.original.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.original.original.25.25.100,f1.Score.sPLS.regression.original.original.25.25.100,f1.Score.sPLS.canonical.original.original.25.25.100), "Scenario" = "Low-Dimensional", "Metric" = "F1-Score")

sparsity.all.methods.original.original.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.original.original.25.25.100[[X]]$index$Microbiome) +length(sCCA.original.original.25.25.100[[X]]$index$Metabolome))/50 )
                                                                                                                                  ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.original.original.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.original.original.25.25.100[[X]]$Metabolome))/50)
                                                                                                                                  ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.original.original.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.original.original.25.25.100[[X]]$Metabolome))/50)), "Scenario"="Low-Dimensional", "Metric"="Sparsity")
df.results.original.original.multi.penalized.25.25.100 = rbind(F1Score.all.methods.original.original.25.25.100,sparsity.all.methods.original.original.25.25.100)

plot.F1Score.original.original.25.25.100 = ggplot(F1Score.all.methods.original.original.25.25.100, aes(x=Method,y=F1.Score, fill=Method))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("F1-Score")+theme(legend.position = "none")+ylim(c(0,1))

plot.sparsity.original.original.25.25.100 = ggplot(sparsity.all.methods.original.original.25.25.100, aes(x=Method,y=Sparsity, fill=Method))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("% Selected Features")+theme(legend.position = "none")+ylim(c(0,1))

ggarrange(ggplot(df.results.original.log.multi.penalized.25.25.100, aes(x=Method,y=Value,fill=Method))+geom_boxplot()+theme_bw()+
            scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("%")+theme(legend.position = "none")+ylim(c(0,1))+facet_grid(~Metric),
          ggplot(df.results.original.log.multi.penalized.100.100.500, aes(x=Method,y=Value,fill=Method))+geom_boxplot()+theme_bw()+
            scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("%")+theme(legend.position = "none")+ylim(c(0,1))+facet_grid(~Metric),ggplot(df.results.original.original.multi.penalized.25.25.100, aes(x=Method,y=Value,fill=Method))+geom_boxplot()+theme_bw()+
            scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("%")+theme(legend.position = "none")+ylim(c(0,1))+facet_grid(~Metric),
          ggplot(df.results.original.original.multi.penalized.100.100.500, aes(x=Method,y=Value,fill=Method))+geom_boxplot()+theme_bw()+
            scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+ylab("%")+theme(legend.position = "none")+ylim(c(0,1))+facet_grid(~Metric), labels=c("A","B","C","D"), ncol=2, nrow=2)

#Log-Contrast 
pvalues.null.LCLR.25.25.100 = readRDS("data\\results_LogContrast\\pvalues_null_LCLR_25_25_100indiv.RDS")
pvalues.null.LCLR.100.100.500 = readRDS("data\\results_LogContrast\\pvalues_null_LCLR_100_100_500indiv.RDS")


pvalues.alternative.LCLR.25.25.100 = readRDS("data\\results_LogContrast\\pvalues_alternative_LCLR_25_25_100indiv.RDS")
pvalues.alternative.LCLR.100.100.500 = readRDS("data\\results_LogContrast\\pvalues_alternative_LCLR_100_100_500indiv.RDS")


df.power.LCLR = data.frame("Power" = c(sum(unlist(pvalues.alternative.LCLR.100.100.500)<=0.05)/(1000*100),
                                       sum(unlist(pvalues.alternative.LCLR.25.25.100)<=0.05)/(1000*25)),
                           "Scenario" = c("High Dimensional", "Low Dimensional"), "Method" = "Log-Contrast")

pvalues.null.MiRKAT.100.100.500 = readRDS("data\\list.results.null.MiRKAT.100.100.500.RDS")
pvalues.null.MiRKAT.25.25.100 = readRDS("data\\list.results.null.MiRKAT.25.25.100.RDS")

pvalues.alter.ilr.log.MiRKAT.25.25.100 = readRDS("data\\power.p.values.25.25.100indiv.MiRKAT.ilr.log.unlist.RDS")
pvalues.alter.ilr.log.MiRKAT.100.100.500 = readRDS("data\\power.p.values.100.100.500indiv.MiRKAT.irl.log.RDS")

df.power.MiRKAT = data.frame("Power" = c(sum(unlist(pvalues.alter.ilr.log.MiRKAT.100.100.500)<=0.05)/(1000*100),
                                       sum(pvalues.alter.ilr.log.MiRKAT.25.25.100<=0.05)/(1000*25)),
                           "Scenario" = c("High Dimensional", "Low Dimensional"), "Method" = "MiRKAT")

df.power.compositional.predictors = rbind(df.power.LCLR,df.power.MiRKAT)
df.power.compositional.predictors$Scenario = factor(df.power.compositional.predictors$Scenario, levels = c("Low Dimensional","High Dimensional"))
df.power.compositional.predictors$Method = factor(df.power.compositional.predictors$Method, levels = c("MiRKAT","Log-Contrast"))
df.power.compositional.predictors$Scenario_Method = paste0(df.power.compositional.predictors$Scenario,df.power.compositional.predictors$Method)


#CODA-LASSO

coda.LASSO.25.25.100 = readRDS("data\\CODA_LASSO\\coda_LASSO_25_25_100indiv.RDS")
coda.LASSO.100.100.500 = readRDS("data\\CODA_LASSO\\coda_LASSO_100_100_500indiv.RDS")



true.associations.25.25.100 = lapply(1:1000, function(y){
  d1 = do.call("rbind", lapply(names(simulated.data.25.25.100.alter$Associations[[y]]$Microbiote), function(x){
    expand.grid(as.numeric(x), simulated.data.25.25.100.alter$Associations[[y]]$Microbiote[[x]])
  }))
  
  colnames(d1) = c("Microbiote","Metabolite")
  d1 = d1[,c("Metabolite","Microbiote")]
  
  d2 = do.call("rbind", lapply(names(simulated.data.25.25.100.alter$Associations[[y]]$Metabolite), function(x){
    expand.grid(as.numeric(x), simulated.data.25.25.100.alter$Associations[[y]]$Metabolite[[x]])
  }))
  
  colnames(d2) = c("Metabolite", "Microbiote")
  
  d = rbind(d1, d2)
})

true.associations.100.100.500 = lapply(1:1000, function(y){
  d1 = do.call("rbind", lapply(names(simulated.data.100.100.500.alter$Associations[[y]]$Microbiote), function(x){
    expand.grid(as.numeric(x), simulated.data.100.100.500.alter$Associations[[y]]$Microbiote[[x]])
  }))
  
  colnames(d1) = c("Microbiote","Metabolite")
  d1 = d1[,c("Metabolite","Microbiote")]
  d2 = do.call("rbind", lapply(names(simulated.data.100.100.500.alter$Associations[[y]]$Metabolite), function(x){
    expand.grid(as.numeric(x), simulated.data.100.100.500.alter$Associations[[y]]$Metabolite[[x]])
  }))
  
  colnames(d2) = c("Metabolite", "Microbiote")
  
  d = rbind(d1, d2)
})





associations.coda.LASSO.25.25.100 = lapply(1:1000, function(x) {
  
  d = do.call("rbind", lapply(which(sapply(coda.LASSO.25.25.100[[x]], length)>0), function(y) {
  expand.grid(y, coda.LASSO.25.25.100[[x]][[y]] )}))
  
  if(!is.null(d)){
    colnames(d) = c("Metabolite", "Microbiote")
  } else {
    d = NULL
  }
  
  d
})

associations.coda.LASSO.100.100.500 = lapply(1:100, function(x) {
  
  d = do.call("rbind", lapply(which(sapply(coda.LASSO.100.100.500[[x]], length)>0), function(y) {
    expand.grid(y, coda.LASSO.100.100.500[[x]][[y]] )}))
  
  if(!is.null(d)){
    colnames(d) = c("Metabolite", "Microbiote")
  } else {
    d = NULL
  }
  
  d
})



all.possible.pairs.25.25.100 = 625
all.possible.pairs.100.100.500 = 10000

F1Score.CODA.LASSO.25.25.100 = sapply(1:1000, function(rep){
  
  if(is.null(associations.coda.LASSO.25.25.100[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.25.25.100[[rep]]),0, all.possible.pairs.25.25.100 - nrow(true.associations.25.25.100[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.25.25.100[[rep]], 1, function(row1) {
      apply(associations.coda.LASSO.25.25.100[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.coda.LASSO.25.25.100[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.25.25.100[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.25.25.100 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    
    f1.score(confusion.matrix)
  }

})

F1Score.CODA.LASSO.100.100.500 = sapply(1:100, function(rep){
  
  if(is.null(associations.coda.LASSO.100.100.500[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.100.100.500[[rep]]),0, all.possible.pairs.100.100.500 - nrow(true.associations.100.100.500[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.100.100.500[[rep]], 1, function(row1) {
      apply(associations.coda.LASSO.100.100.500[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.coda.LASSO.100.100.500[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.100.100.500[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.100.100.500 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    
    f1.score(confusion.matrix)
  }
  
})

sparsity.CODA.LASSO.25.25.100 = sapply(1:1000, function(x) {
  if(length(associations.coda.LASSO.25.25.100[[x]]) > 0){
    nrow(associations.coda.LASSO.25.25.100[[x]])/625
  } else {
    0
  }
})

sparsity.CODA.LASSO.100.100.500 = sapply(1:100, function(x) {
  if(length(associations.coda.LASSO.100.100.500[[x]]) > 0){
    nrow(associations.coda.LASSO.100.100.500[[x]])/10000
  } else {
    0
  }
})

df.results.CODA.LASSO.25.25.100 = data.frame("Value" = c(F1Score.CODA.LASSO.25.25.100, sparsity.CODA.LASSO.25.25.100), "Scenario" = "Low-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"= "CODA-LASSO")
df.results.CODA.LASSO.25.25.100$Metric = factor(df.results.CODA.LASSO.25.25.100$Metric, levels = c("Sparsity","F1-Score"))
df.results.CODA.LASSO.100.100.500 = data.frame("Value" = c(F1Score.CODA.LASSO.100.100.500, sparsity.CODA.LASSO.100.100.500), "Scenario" = "High-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=100), "Method"="CODA-LASSO")
df.results.CODA.LASSO.100.100.500$Metric = factor(df.results.CODA.LASSO.100.100.500$Metric, levels = c("Sparsity","F1-Score"))




#LASSO
LASSO.log.clr.25.25.100 = readRDS("data\\LASSO\\coefs_LASSO_log_Metabolome_clr_Microbiome_25_25_100indiv.RDS")
LASSO.log.clr.100.100.500 = readRDS("data\\LASSO\\coefs_LASSO_log_Metabolome_clr_Microbiome_100_100_500indiv.RDS")


associations.LASSO.log.clr.25.25.100 = lapply(1:1000, function(x) {
  m = as.matrix(do.call("cbind",LASSO.log.clr.25.25.100[[x]]))
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:25)
  do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
  expand.grid(which(colSums(m)!=0)[y], apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0))[[y]])
}))
})

F1Score.LASSO.log.clr.25.25.100 = sapply(1:1000, function(rep){
  
  if(is.null(associations.LASSO.log.clr.25.25.100[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.25.25.100[[rep]]),0, all.possible.pairs.25.25.100 - nrow(true.associations.25.25.100[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.25.25.100[[rep]], 1, function(row1) {
      apply(associations.LASSO.log.clr.25.25.100[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.LASSO.log.clr.25.25.100[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.25.25.100[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.25.25.100 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    
    f1.score(confusion.matrix)
  }
  
})
associations.LASSO.log.clr.100.100.500 = lapply(1:1000, function(x) {
  m = as.matrix(do.call("cbind",LASSO.log.clr.100.100.500[[x]]))
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:100)
  do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
    expand.grid(which(colSums(m)!=0)[y], apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0))[[y]])
  }))
})

F1Score.LASSO.log.clr.100.100.500 = sapply(1:1000, function(rep){
  
  if(is.null(associations.LASSO.log.clr.100.100.500[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.100.100.500[[rep]]),0, all.possible.pairs.100.100.500 - nrow(true.associations.100.100.500[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.100.100.500[[rep]], 1, function(row1) {
      apply(associations.LASSO.log.clr.100.100.500[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.LASSO.log.clr.100.100.500[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.100.100.500[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.100.100.500 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    
    f1.score(confusion.matrix)
  }
  
})


sparsity.LASSO.log.clr.25.25.100 = sapply(1:1000, function(x) {
  if(length(associations.LASSO.log.clr.25.25.100[[x]]) > 0){
    nrow(associations.LASSO.log.clr.25.25.100[[x]])/625
  } else {
    0
  }
})

sparsity.LASSO.log.clr.100.100.500 = sapply(1:1000, function(x) {
  if(length(associations.LASSO.log.clr.100.100.500[[x]]) > 0){
    nrow(associations.LASSO.log.clr.100.100.500[[x]])/10000
  } else {
    0
  }
})

df.results.LASSO.log.clr.25.25.100 = data.frame("Value" = c(F1Score.LASSO.log.clr.25.25.100, sparsity.LASSO.log.clr.25.25.100), "Scenario" = "Low-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"= "clr-LASSO")
df.results.LASSO.log.clr.25.25.100$Metric = factor(df.results.LASSO.log.clr.25.25.100$Metric, levels=c("Sparsity","F1-Score"))
df.results.LASSO.log.clr.100.100.500 = data.frame("Value" = c(F1Score.LASSO.log.clr.100.100.500, sparsity.LASSO.log.clr.100.100.500), "Scenario" = "High-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"="clr-LASSO")
df.results.LASSO.log.clr.100.100.500$Metric = factor(df.results.LASSO.log.clr.100.100.500$Metric, levels=c("Sparsity","F1-Score"))


MLASSO.log.clr.25.25.100 = readRDS("data\\LASSO\\coefs_MLASSO_log_Metabolome_clr_Microbiome_25_25_100indiv.RDS")
MLASSO.log.clr.100.100.500 = readRDS("data\\LASSO\\coefs_MLASSO_log_Metabolome_clr_Microbiome_100_100_500indiv.RDS")

associations.MLASSO.log.clr.25.25.100 = lapply(1:1000, function(x) {
  print(x)
  m = as.matrix(do.call("cbind",MLASSO.log.clr.25.25.100[[x]]))
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:25)
  if(length(which(colSums(m)!=0)) >0) {
    do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
      expand.grid(which(colSums(m)!=0)[y], matrix(apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0)), ncol=length(which(colSums(m)!=0)))[,y,drop=F])
    }))
  } else {
    NULL
  }
})

F1Score.MLASSO.log.clr.25.25.100 = sapply(1:1000, function(rep){
  
  if(is.null(associations.MLASSO.log.clr.25.25.100[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.25.25.100[[rep]]),0, all.possible.pairs.25.25.100 - nrow(true.associations.25.25.100[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.25.25.100[[rep]], 1, function(row1) {
      apply(associations.MLASSO.log.clr.25.25.100[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.MLASSO.log.clr.25.25.100[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.25.25.100[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.25.25.100 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  }
  
})


associations.MLASSO.log.clr.100.100.500 = lapply(1:1000, function(x) {
  print(x)
  m = as.matrix(do.call("cbind",MLASSO.log.clr.100.100.500[[x]]))
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:100)
  if(length(which(colSums(m)!=0)) >0) {
    do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
      expand.grid(which(colSums(m)!=0)[y], matrix(apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0)), ncol=length(which(colSums(m)!=0)))[,y,drop=F])
    }))
  } else {
    NULL
  }
})

F1Score.MLASSO.log.clr.100.100.500 = sapply(1:1000, function(rep){
  
  if(is.null(associations.MLASSO.log.clr.100.100.500[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.100.100.500[[rep]]),0, all.possible.pairs.100.100.500 - nrow(true.associations.100.100.500[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.100.100.500[[rep]], 1, function(row1) {
      apply(associations.MLASSO.log.clr.100.100.500[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.MLASSO.log.clr.100.100.500[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.100.100.500[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.100.100.500 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  }
  
})


sparsity.MLASSO.log.clr.25.25.100 = sapply(1:1000, function(x) {
  if(length(associations.MLASSO.log.clr.25.25.100[[x]]) > 0){
    nrow(associations.MLASSO.log.clr.25.25.100[[x]])/625
  } else {
    0
  }
})

sparsity.MLASSO.log.clr.100.100.500 = sapply(1:1000, function(x) {
  if(length(associations.MLASSO.log.clr.100.100.500[[x]]) > 0){
    nrow(associations.MLASSO.log.clr.100.100.500[[x]])/10000
  } else {
    0
  }
})

df.results.MLASSO.log.clr.25.25.100 = data.frame("Value" = c(F1Score.MLASSO.log.clr.25.25.100, sparsity.MLASSO.log.clr.25.25.100), "Scenario" = "Low-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"= "clr-MLASSO")
df.results.MLASSO.log.clr.25.25.100$Metric = factor(df.results.MLASSO.log.clr.25.25.100$Metric, levels=c("Sparsity", "F1-Score"))
df.results.MLASSO.log.clr.100.100.500 = data.frame("Value" = c(F1Score.MLASSO.log.clr.100.100.500, sparsity.MLASSO.log.clr.100.100.500), "Scenario" = "High-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"="clr-MLASSO")
df.results.MLASSO.log.clr.100.100.500$Metric = factor(df.results.MLASSO.log.clr.100.100.500$Metric, levels=c("Sparsity", "F1-Score"))

df.results.penalized.univariate.log.metabolome.25.25.100 = rbind(df.results.CODA.LASSO.25.25.100,df.results.LASSO.log.clr.25.25.100,df.results.MLASSO.log.clr.25.25.100,df.results.MLASSO.log.clr.100.100.500)
df.results.penalized.univariate.log.metabolome.100.100.500 = rbind(df.results.CODA.LASSO.100.100.500,df.results.LASSO.log.clr.100.100.500,df.results.MLASSO.log.clr.100.100.500)


ggarrange(ggplot(df.results.penalized.univariate.log.metabolome.25.25.100, aes(x=Method, y=Value, fill=Method))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=wes_palette(n=3, name="AsteroidCity3"))+ylab("%")+theme(legend.position = "none")+ylim(c(0,1))+facet_grid(~Metric),
ggplot(df.results.penalized.univariate.log.metabolome.100.100.500, aes(x=Method, y=Value, fill=Method))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=wes_palette(n=3, name="AsteroidCity3"))+ylab("%")+theme(legend.position = "none")+ylim(c(0,1))+facet_grid(~Metric), labels = c("A","B"))


LASSO.clr.log.25.25.100 = readRDS("data\\LASSO\\coefs_LASSO_clr_Microbiome_log_Metabolome_25_25_100indiv.RDS")
LASSO.clr.log.100.100.500 = readRDS("data\\LASSO\\coefs_LASSO_clr_Microbiome_log_Metabolome_100_100_500indiv.RDS")




associations.LASSO.clr.log.25.25.100 = lapply(1:1000, function(x) {
  m = as.matrix(do.call("cbind",LASSO.clr.log.25.25.100[[x]]))
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:25)
  do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
    expand.grid(which(colSums(m)!=0)[y], apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0))[[y]])
  }))
})

F1Score.LASSO.clr.log.25.25.100 = sapply(1:1000, function(rep){
  
  if(is.null(associations.LASSO.clr.log.25.25.100[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.25.25.100[[rep]]),0, all.possible.pairs.25.25.100 - nrow(true.associations.25.25.100[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.25.25.100[[rep]], 1, function(row1) {
      apply(associations.LASSO.clr.log.25.25.100[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.LASSO.clr.log.25.25.100[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.25.25.100[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.25.25.100 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    
    f1.score(confusion.matrix)
  }
  
})
associations.LASSO.clr.log.100.100.500 = lapply(1:1000, function(x) {
  m = as.matrix(do.call("cbind",LASSO.clr.log.100.100.500[[x]]))
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:100)
  do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
    expand.grid(which(colSums(m)!=0)[y], apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0))[[y]])
  }))
})

F1Score.LASSO.clr.log.100.100.500 = sapply(1:1000, function(rep){
  
  if(is.null(associations.LASSO.clr.log.100.100.500[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.100.100.500[[rep]]),0, all.possible.pairs.100.100.500 - nrow(true.associations.100.100.500[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.100.100.500[[rep]], 1, function(row1) {
      apply(associations.LASSO.clr.log.100.100.500[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.LASSO.clr.log.100.100.500[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.100.100.500[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.100.100.500 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    
    f1.score(confusion.matrix)
  }
  
})

sparsity.LASSO.clr.log.25.25.100 = sapply(1:1000, function(x) {
  if(length(associations.LASSO.clr.log.25.25.100[[x]]) > 0){
    nrow(associations.LASSO.clr.log.25.25.100[[x]])/625
  } else {
    0
  }
})

sparsity.LASSO.clr.log.100.100.500 = sapply(1:1000, function(x) {
  if(length(associations.LASSO.clr.log.100.100.500[[x]]) > 0){
    nrow(associations.LASSO.clr.log.100.100.500[[x]])/10000
  } else {
    0
  }
})

df.results.LASSO.clr.log.25.25.100 = data.frame("Value" = c(F1Score.LASSO.clr.log.25.25.100, sparsity.LASSO.clr.log.25.25.100), "Scenario" = "Low-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"= "LASSO")
df.results.LASSO.clr.log.25.25.100$Metric = factor(df.results.LASSO.clr.log.25.25.100$Metric, levels = c("Sparsity", "F1-Score"))
df.results.LASSO.clr.log.100.100.500 = data.frame("Value" = c(F1Score.LASSO.clr.log.100.100.500, sparsity.LASSO.clr.log.100.100.500), "Scenario" = "High-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"="LASSO")
df.results.LASSO.clr.log.100.100.500$Metric = factor(df.results.LASSO.clr.log.100.100.500$Metric, levels = c("Sparsity", "F1-Score"))

#MLASSO
MLASSO.clr.log.25.25.100 = readRDS("data\\LASSO\\coefs_MLASSO_clr_Microbiome_log_Metabolome_25_25_100indiv.RDS")
MLASSO.clr.log.100.100.500 = readRDS("data\\LASSO\\coefs_MLASSO_clr_Microbiome_log_Metabolome_100_100_500indiv.RDS")


associations.MLASSO.clr.log.25.25.100 = lapply(1:1000, function(x) {
  print(x)
  m = as.matrix(do.call("cbind",MLASSO.clr.log.25.25.100[[x]]))
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:25)
  if(length(which(colSums(m)!=0)) >0) {
    do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
      expand.grid(which(colSums(m)!=0)[y], matrix(apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0)), ncol=length(which(colSums(m)!=0)))[,y,drop=F])
    }))
  } else {
    NULL
  }
})

F1Score.MLASSO.clr.log.25.25.100 = sapply(1:1000, function(rep){
  
  if(is.null(associations.MLASSO.clr.log.25.25.100[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.25.25.100[[rep]]),0, all.possible.pairs.25.25.100 - nrow(true.associations.25.25.100[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.25.25.100[[rep]], 1, function(row1) {
      apply(associations.MLASSO.clr.log.25.25.100[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.MLASSO.clr.log.25.25.100[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.25.25.100[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.25.25.100 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  }
  
})


associations.MLASSO.clr.log.100.100.500 = lapply(1:1000, function(x) {
  print(x)
  m = as.matrix(do.call("cbind",MLASSO.clr.log.100.100.500[[x]]))
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:100)
  if(length(which(colSums(m)!=0)) >0) {
    do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
      expand.grid(which(colSums(m)!=0)[y], matrix(apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0)), ncol=length(which(colSums(m)!=0)))[,y,drop=F])
    }))
  } else {
    NULL
  }
})

F1Score.MLASSO.clr.log.100.100.500 = sapply(1:1000, function(rep){
  
  if(is.null(associations.MLASSO.clr.log.100.100.500[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.100.100.500[[rep]]),0, all.possible.pairs.100.100.500 - nrow(true.associations.100.100.500[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.100.100.500[[rep]], 1, function(row1) {
      apply(associations.MLASSO.clr.log.100.100.500[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.MLASSO.clr.log.100.100.500[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.100.100.500[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.100.100.500 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  }
  
})


sparsity.MLASSO.clr.log.25.25.100 = sapply(1:1000, function(x) {
  if(length(associations.MLASSO.clr.log.25.25.100[[x]]) > 0){
    nrow(associations.MLASSO.clr.log.25.25.100[[x]])/625
  } else {
    0
  }
})

sparsity.MLASSO.clr.log.100.100.500 = sapply(1:1000, function(x) {
  if(length(associations.MLASSO.clr.log.100.100.500[[x]]) > 0){
    nrow(associations.MLASSO.clr.log.100.100.500[[x]])/10000
  } else {
    0
  }
})

df.results.MLASSO.clr.log.25.25.100 = data.frame("Value" = c(F1Score.MLASSO.clr.log.25.25.100, sparsity.MLASSO.clr.log.25.25.100), "Scenario" = "Low-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"= "MLASSO")
df.results.MLASSO.clr.log.25.25.100$Metric = factor(df.results.MLASSO.clr.log.25.25.100$Metric,levels=c("Sparsity", "F1-Score"))
df.results.MLASSO.clr.log.100.100.500 = data.frame("Value" = c(F1Score.MLASSO.clr.log.100.100.500, sparsity.MLASSO.clr.log.100.100.500), "Scenario" = "High-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"="MLASSO")
df.results.MLASSO.clr.log.100.100.500$Metric = factor(df.results.MLASSO.clr.log.100.100.500$Metric, levels=c("Sparsity", "F1-Score"))

#Sparse Dirichlet
sparse.DM.log.25.25.100 = readRDS("data\\LASSO\\sparse_DM_log_Metabolome_25_25_100indiv.RDS")
sparse.DM.log.100.100.500 = readRDS("data\\LASSO\\sparse_DM_log_Metabolome_100_100_500indiv.RDS")

sparse.GDM.log.25.25.100 = readRDS("data\\LASSO\\sparse_GDM_log_Metabolome_25_25_100indiv.RDS")
sparse.GDM.log.100.100.500 = readRDS("data\\LASSO\\sparse_GDM_log_Metabolome_100_100_500indiv.RDS")


associations.sparse.DM.log.25.25.100 = lapply(1:1000, function(x) {
  print(x)
  m = sparse.DM.log.25.25.100[[x]]
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:25)
  
  if(length(which(colSums(m)!=0)) >0) {
    do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
      expand.grid(which(colSums(m)!=0)[y], matrix(apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0)), ncol=length(which(colSums(m)!=0)))[[y]])
    }))
  } else {
    NULL
  }
})

F1Score.sparse.DM.log.25.25.100 = sapply(1:1000, function(rep){
  
  if(is.null(associations.sparse.DM.log.25.25.100[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.25.25.100[[rep]]),0, all.possible.pairs.25.25.100 - nrow(true.associations.25.25.100[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.25.25.100[[rep]], 1, function(row1) {
      apply(associations.sparse.DM.log.25.25.100[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.sparse.DM.log.25.25.100[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.25.25.100[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.25.25.100 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  }
  
})


associations.sparse.DM.log.100.100.500 = lapply(1:100, function(x) {
  print(x)
  m = sparse.DM.log.100.100.500[[x]]
  m = m[-1,]
  colnames(m) = paste0("Metabolite",1:100)
  
  if(length(which(colSums(m)!=0)) >0) {
    do.call("rbind",lapply(1:length(which(colSums(m)!=0)), function(y){
      expand.grid(which(colSums(m)!=0)[y], matrix(apply(m[,which(colSums(m)!=0)],2,function(col) which(col!=0)), ncol=length(which(colSums(m)!=0)))[[y]])
    }))
  } else {
    NULL
  }
})

F1Score.sparse.DM.log.100.100.500 = sapply(1:100, function(rep){
  
  if(is.null(associations.sparse.DM.log.100.100.500[[rep]])){
    confusion.matrix = matrix(c(0, nrow(true.associations.100.100.500[[rep]]),0, all.possible.pairs.100.100.500 - nrow(true.associations.100.100.500[[rep]])), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  } else {
    tre = apply(true.associations.100.100.500[[rep]], 1, function(row1) {
      apply(associations.sparse.DM.log.100.100.500[[rep]], 1, function(row2){
        sum(all(row1 == row2))
      })
    })
    
    false.pos = nrow(associations.sparse.DM.log.100.100.500[[rep]]) - sum(tre ==1 )
    true.pos = sum(tre ==1 )
    false.neg = nrow(true.associations.100.100.500[[rep]]) - sum(tre ==1 )
    true.neg = all.possible.pairs.100.100.500 - false.pos - true.pos- false.neg
    
    confusion.matrix = matrix(c(true.pos, false.neg,false.pos, true.neg), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg")))
    #print(confusion.matrix)
    f1.score(confusion.matrix)
  }
  
})

sparsity.sparse.DM.log.25.25.100 = sapply(1:1000, function(x) {
  if(length(associations.sparse.DM.log.25.25.100[[x]]) > 0){
    nrow(associations.sparse.DM.log.25.25.100[[x]])/625
  } else {
    0
  }
})

sparsity.sparse.DM.log.100.100.500 = sapply(1:100, function(x) {
  if(length(associations.sparse.DM.log.100.100.500[[x]]) > 0){
    nrow(associations.sparse.DM.log.100.100.500[[x]])/10000
  } else {
    0
  }
})

df.results.sparse.DM.log.25.25.100 = data.frame("Value" = c(F1Score.sparse.DM.log.25.25.100, sparsity.sparse.DM.log.25.25.100), "Scenario" = "Low-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=1000), "Method"= "Sparse Dirichlet")
df.results.sparse.DM.log.25.25.100$Metric = factor(df.results.sparse.DM.log.25.25.100$Metric, levels = c("Sparsity", "F1-Score"))
df.results.sparse.DM.log.100.100.500 = data.frame("Value" = c(F1Score.sparse.DM.log.100.100.500, sparsity.sparse.DM.log.100.100.500), "Scenario" = "High-Dimensional", "Metric" = rep(c("F1-Score", "Sparsity"), each=100), "Method"= "Sparse Dirichlet")
df.results.sparse.DM.log.100.100.500$Metric = factor(df.results.sparse.DM.log.100.100.500$Metric,levels = c("Sparsity", "F1-Score"))

df.results.penalized.univariate.clr.microbiome.25.25.100 = rbind(df.results.sparse.DM.log.25.25.100,df.results.LASSO.clr.log.25.25.100, df.results.MLASSO.clr.log.25.25.100)
df.results.penalized.univariate.clr.microbiome.25.25.100$Method = factor(df.results.penalized.univariate.clr.microbiome.25.25.100$Method, levels = c("LASSO", "MLASSO","Sparse Dirichlet"))

df.results.penalized.univariate.clr.microbiome.100.100.500 = rbind(df.results.sparse.DM.log.100.100.500,df.results.LASSO.clr.log.100.100.500, df.results.MLASSO.clr.log.100.100.500)
df.results.penalized.univariate.clr.microbiome.100.100.500$Method = factor(df.results.penalized.univariate.clr.microbiome.100.100.500$Method, levels = c("LASSO", "MLASSO","Sparse Dirichlet"))

ggarrange(ggarrange(ggplot(df.results.penalized.univariate.log.metabolome.25.25.100, aes(x=Method, y=Value*100, fill=Method))+geom_violin()+theme_bw()+
            scale_fill_simpsons()+ylab("%")+theme(legend.position = "none")+ylim(c(0,100))+facet_grid(~Metric)+ggtitle("Low-Dimensional"),
          ggplot(df.results.penalized.univariate.log.metabolome.100.100.500, aes(x=Method, y=Value*100, fill=Method))+geom_violin()+theme_bw()+
            scale_fill_simpsons()+ylab("%")+theme(legend.position = "none")+ylim(c(0,100))+facet_grid(~Metric)+ggtitle("High-Dimensional"), labels = c("A","D")),

ggarrange(ggplot(df.results.penalized.univariate.clr.microbiome.25.25.100, aes(x=Method, y=Value*100, fill=Method))+geom_violin()+theme_bw()+
            scale_fill_uchicago()+ylab("%")+theme(legend.position = "none")+ylim(c(0,100))+facet_grid(~Metric),
          ggplot(df.results.penalized.univariate.clr.microbiome.100.100.500, aes(x=Method, y=Value*100, fill=Method))+geom_violin()+theme_bw()+
            scale_fill_uchicago()+ylab("%")+theme(legend.position = "none")+ylim(c(0,100))+facet_grid(~Metric), labels = c("B","E")), 
plot.results.penalized.multivariate, nrow=3)


#Regression
pvalues.null.lm.clr.log.25.25.100 = readRDS("data\\Regression\\25_25_100\\null\\pvalues_null_lm_clr_microbiome_log_metabolome_25_25_100indiv.RDS")
pvalues.null.Dir.log.25.25.100 = readRDS("data\\Regression\\25_25_100\\null\\pvalues_null_Dir_log_metabolome_25_25_100indiv.RDS")
pvalues.null.log.clr.25.25.100 = readRDS("data\\Regression\\25_25_100\\null\\pvalues_null_loglinear_clr_microbiome_25_25_100indiv.RDS")


pvalues.alternative.lm.clr.log.25.25.100 = readRDS("data\\Regression\\25_25_100\\alternative\\pvalues_alternative_lm_clr_microbiome_log_metabolome_25_25_100indiv.RDS")
pvalues.alternative.Dir.log.25.25.100 = readRDS("data\\Regression\\25_25_100\\alternative\\pvalues_alternative_Dir_log_metabolome_25_25_100indiv.RDS")
pvalues.alternative.log.clr.25.25.100 = readRDS("data\\Regression\\25_25_100\\alternative\\pvalues_alter_loglinear_clr_microbiome_25_25_100indiv.RDS")

pvalues.null.log.clr.100.100.500 = readRDS("data\\Regression\\100_100_500\\null\\pvalues_null_loglinear_clr_microbiome_100_100_500indiv.RDS")
pvalues.alternative.log.clr.100.100.500 = readRDS("data\\Regression\\100_100_500\\alternative\\pvalues_alter_loglinear_clr_microbiome_100_100_500indiv.RDS")

#Correlation
spearman.null.clr.log.25.25.100 = readRDS("data\\Correlation\\25_25_100\\Spearman\\null\\pvalues_null_spearman_clr_microbiome_metabolome_25_25_100indiv.RDS")
pearson.null.clr.log.25.25.100 = readRDS("data\\Correlation\\25_25_100\\Pearson\\null\\pvalues_null_pearson_clr_microbiome_log_metabolome_25_25_100indiv.RDS")

spearman.alternative.clr.log.25.25.100 = readRDS("data\\Correlation\\25_25_100\\Spearman\\alternative\\pvalues_alternative_spearman_clr_microbiome_metabolome_25_25_100indiv.RDS")
pearson.alternative.clr.log.25.25.100 = readRDS("data\\Correlation\\25_25_100\\Pearson\\alternative\\pvalues_alternative_pearson_clr_microbiome_log_metabolome_25_25_100indiv.RDS")

spearman.alternative.clr.log.100.100.500 = readRDS("data\\Correlation\\100_100_500\\Spearman\\alternative\\pvalues_alternative_spearman_clr_microbiome_metabolome_100_100_500indiv.RDS")
pearson.alternative.clr.log.100.100.500 = readRDS("data\\Correlation\\100_100_500\\Pearson\\alternative\\pvalues_alternative_pearson_clr_microbiome_log_metabolome_100_100_500indiv.RDS")


df.power.compositional.predictors.low.dimensional = data.frame("Scenario"="Low-Dimensional"
  ,"Method"=c("Log-Contrast", "MiRKAT", "clr-lm"), 
                                               "Power" = c(sum(unlist(pvalues.alternative.LCLR.25.25.100) <= 0.05)/25000,
                                                           sum(pvalues.alter.ilr.log.MiRKAT.25.25.100<=0.05)/25000,
                                                           sum(unlist(lapply(1:1000, function(x) sapply(pvalues.alternative.log.clr.25.25.100[[x]], ACAT::ACAT)))<=0.05)/25000))

#df.power.compositional.predictors.low.dimensional$Method = factor(df.power.compositional.predictors.low.dimensional$Method, levels=c("Log-Contrast", "MiRKAT", "clr-lm"))

df.power.compositional.predictors.high.dimensional = data.frame("Scenario"="High-Dimensional",
  "Method"=c("Log-Contrast", "MiRKAT", "clr-lm"), 
                                                               "Power" = c(sum(unlist(pvalues.alternative.LCLR.100.100.500) <= 0.05)/100000,
                                                                           sum(unlist(pvalues.alter.ilr.log.MiRKAT.100.100.500)<=0.05)/100000,
                                                                           sum(unlist(lapply(1:1000, function(x) sapply(pvalues.alternative.log.clr.100.100.500[[x]], ACAT::ACAT)))<=0.05)/100000))

#df.power.compositional.predictors.high.dimensional$Method = factor(df.power.compositional.predictors.high.dimensional$Method, levels=c("Log-Contrast", "MiRKAT", "clr-lm"))


df.power.compositional.predictors = rbind(df.power.compositional.predictors.low.dimensional,df.power.compositional.predictors.high.dimensional)
df.power.compositional.predictors$Scenario = factor(df.power.compositional.predictors$Scenario, levels=c("Low-Dimensional", "High-Dimensional"))


ggarrange(ggarrange(gg_qqplot_facet_grid(list("Log-Contrast"= unlist(pvalues.null.LCLR.25.25.100), "MiRKAT"=unlist(pvalues.null.MiRKAT.25.25.100$null.p.values.25.25.100indiv.MiRKAT.ilr.log),
                                              "clr-lm"=unlist(lapply(1:1000, function(x) sapply(pvalues.null.log.clr.25.25.100[[x]], ACAT::ACAT)))))+ggtitle("Low-Dimensional")+theme(legend.position = c(0.7, 0.15))
                    ,gg_qqplot_facet_grid(list("Log-Contrast"= unlist(pvalues.null.LCLR.100.100.500), "MiRKAT"=unlist(pvalues.null.MiRKAT.100.100.500$null.p.values.100.100.500indiv.MiRKAT.ilr.log),
                                               "clr-lm"=unlist(lapply(1:1000, function(x) sapply(pvalues.null.log.clr.100.100.500[[x]], ACAT::ACAT)))))+ggtitle("High-Dimensional")+theme(legend.position = c(0.7, 0.15))
                    , ncol=2, labels = c("A","B")),ggarrange(ggplot(df.power.compositional.predictors, aes(x=Method,y=Power, fill=Method))+geom_bar(stat="identity", position = "dodge", colour="black")+facet_grid(~Scenario)+
          theme_bw()+scale_fill_npg()+ylab("Power")+theme(legend.position = "none")+ylim(c(0,1)), labels = "C"), ncol=2)



df.power.compositional.outcomes = data.frame("Method"=c("Dirichlet", "lm"), 
                                               "Power" = c(sum(unlist(lapply(1:1000, function(x) sapply(pvalues.alternative.Dir.log.25.25.100[[x]], ACAT::ACAT)))<=0.05)/25000,
                                             sum(unlist(lapply(1:1000, function(x) sapply(pvalues.alternative.lm.clr.log.25.25.100[[x]], ACAT::ACAT)))<=0.05)/25000))

ggarrange(ggarrange(gg_qqplot(unlist(pvalues.null.LCLR.25.25.100), colour="#3da38d", title = "Log-Contrast"),
                    gg_qqplot(unlist(pvalues.null.MiRKAT.25.25.100), colour="#b17340", title = "MiRKAT"),
                    gg_qqplot(unlist(lapply(1:1000, function(x) sapply(pvalues.null.log.clr.25.25.100[[x]], ACAT::ACAT))), colour = "chocolate1", title="clr-lm"),
                    ncol=3),
          ggarrange(ggplot(df.power.compositional.predictors, aes(x=Method,y=Power, fill=Method))+geom_bar(stat="identity", position = "dodge", colour="black") +theme_bw()+
  scale_fill_manual(values=c("#3da38d", "#b17340", "chocolate1" ))+ylab("Power")+theme(legend.position = "none")+ylim(c(0,1))),
  ggarrange(gg_qqplot(unlist(lapply(1:1000, function(x) sapply(pvalues.null.Dir.log.25.25.100[[x]], ACAT::ACAT))), colour = "deepskyblue", title = "Dirichlet"),
            gg_qqplot(unlist(lapply(1:1000, function(x) sapply(pvalues.null.lm.clr.log.25.25.100[[x]], ACAT::ACAT))), colour = "darkgoldenrod1", title="lm")),
  ggplot(df.power.compositional.outcomes, aes(x=Method,y=Power, fill=Method))+geom_bar(stat="identity", position = "dodge", colour="black") +theme_bw()+
    scale_fill_manual(values=c("deepskyblue", "darkgoldenrod1" ))+ylab("Power")+theme(legend.position = "none")+ylim(c(0,1)),
  labels=c("A","B","C","D"),ncol=2, nrow=2)

ggplot(data=data.frame("pvalues"=unlist(pvalues.null.LCLR.25.25.100)), aes(sample=pvalues)) + geom_qq(distribution = stats::qunif)+geom_qq_line(distribution = stats::qunif)

#Konzo application 
#Some data are not available on the repo
#Please contact Matthew Bramble for data access

library(MOFA2)
#MOFA2
Microbiome.data = read.table("data\\Konzo\\Microbiome_unnormalized.tsv", header=TRUE, sep="\t")
head(Microbiome.data[,1:10])
Metabolome.data = read.table("data\\Konzo\\Metabolome_unnormalized.tsv", header=TRUE, sep="\t")
head(Metabolome.data[,1:10])
metadata.Konzo = read.csv2("data\\Konzo\\Metadata_Konzo2021.csv", header=TRUE, sep=",")

sample.in.Microbiome = which(colnames(Microbiome.data)=="Sample")
sample.in.Metabolome = which(colnames(Metabolome.data)=="Sample")

sample = Microbiome.data$Sample
metadata.Konzo = metadata.Konzo[metadata.Konzo$Sample%in%sample,]

Metabolome.data = Metabolome.data[,-sample.in.Metabolome]
Microbiome.data = Microbiome.data[,-sample.in.Microbiome]

colnames.Microbiome = colnames(Microbiome.data)
colnames.Metabolome = colnames(Metabolome.data)

metabo.chemical.annotations = read.csv2("data\\Konzo\\Chemical_annot_Metabo.csv", sep=",", header=TRUE)
MOFA2.Konzo.output = readRDS("data\\Konzo\\MOFA_Konzo.RDS")
plot_variance_explained(MOFA2.Konzo.output)                                                                                                          
get_variance_explained(MOFA2.Konzo.output)

MOFA2.Weights.Micro.Metabo = get_weights(MOFA2.Konzo.output)

Weights.Metabo.first.factor = as.data.frame(MOFA2.Weights.Micro.Metabo$Metabo[,1,drop=FALSE])
Weights.Metabo.first.factor$sign = ifelse(Weights.Metabo.first.factor[,1] <0, "-", "+")

Weights.Metabo.first.factor$Factor1.scaled = Weights.Metabo.first.factor$Factor1/max(abs(Weights.Metabo.first.factor$Factor1))
order.Weights.Metabo.first.factor = Weights.Metabo.first.factor[order(abs(Weights.Metabo.first.factor$Factor1.scaled), decreasing = TRUE),]
order.Weights.Metabo.first.factor$Metabolite = sub(".","",rownames(order.Weights.Metabo.first.factor))
top10.Metabo = order.Weights.Metabo.first.factor[1:10,]

top10.Metabo$Chemical = NA
for(i in 1:nrow(top10.Metabo)){
  m=as.numeric(top10.Metabo[i,"Metabolite"])
  chemical = metabo.chemical.annotations[metabo.chemical.annotations$CHEM_ID==m,"CHEMICAL_NAME"]
  top10.Metabo$Chemical[i] = chemical
}

ggplot(top10.Metabo, aes(x = reorder(Chemical,-abs(Factor1.scaled)), y = abs(Factor1.scaled), fill=sign, label=sign)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Metabolite")+geom_label(size=3.5)+ylim(c(0,1))


Weights.Microbio.first.factor = as.data.frame(MOFA2.Weights.Micro.Metabo$Micro[,1,drop=FALSE])
Weights.Microbio.first.factor$sign = ifelse(Weights.Microbio.first.factor[,1] <0, "-", "+")

Weights.Microbio.first.factor$Factor1.scaled = Weights.Microbio.first.factor$Factor1/max(abs(Weights.Microbio.first.factor$Factor1))
order.Weights.Microbio.first.factor = Weights.Microbio.first.factor[order(abs(Weights.Microbio.first.factor$Factor1.scaled), decreasing = TRUE),]
order.Weights.Microbio.first.factor$Species = colnames.Microbiome[as.numeric(sub(".*feature_ *(.*?) *_v2.*", "\\1", rownames(order.Weights.Microbio.first.factor)))]

order.Weights.Microbio.first.factor$Species = stringr::str_replace_all(order.Weights.Microbio.first.factor$Species, "\\.", " ")
top10.Microbio = order.Weights.Microbio.first.factor[1:10,]


ggplot(top10.Microbio, aes(x = reorder(Species,-abs(Factor1.scaled)), y = abs(Factor1.scaled), fill=sign, label=sign)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Species")+geom_label(size=3.5)+ylim(c(0,1))

sPLS.sCCA.Konzo.output = readRDS("data\\Konzo\\results_sCCA_sPLS_Konzo.RDS")
non.null.loadings.species.KONZO.sPLS = sPLS.sCCA.Konzo.output$sPLS$loadings$X[rowSums(sPLS.sCCA.Konzo.output$sPLS$loadings$X)!=0,]
non.null.loadings.metabo.KONZO.sPLS = sPLS.sCCA.Konzo.output$sPLS$loadings$Y[rowSums(sPLS.sCCA.Konzo.output$sPLS$loadings$Y)!=0,]

plot.t = mixOmics::plotVar(sPLS.sCCA.Konzo.output$sPLS)
df.sPLS.KONZO.loadings = data.frame("Type" = c(rep("Microorganisms", nrow(non.null.loadings.species.KONZO.sPLS)), rep("Metabolites",nrow(non.null.loadings.metabo.KONZO.sPLS))),
                                    "1st Component" = plot.t$x, "2nd Component"=plot.t$y)




l = list()
for(m in rownames(sPLS.sCCA.Konzo.output$sPLS$loadings$Y[rowSums(sPLS.sCCA.Konzo.output$sPLS$loadings$Y)!=0,])){
  l[[m]] = coda4microbiome::coda_glmnet(Microbiome.data[,rownames(sPLS.sCCA.Konzo.output$sPLS$loadings$X[rowSums(sPLS.sCCA.Konzo.output$sPLS$loadings$X)!=0,])], log(Metabolome.data[,m]))
  
}



library(dplyr)
df.CODA2.nasso = reshape2::melt(table(sapply(l, function(x) length(x$taxa.name))))
df.CODA2.effect.var = reshape2::melt(tapply(unlist(lapply(l, function(x){
  p = x$`log-contrast coefficients`
  names(p) = x$taxa.name
  p
})),sapply(strsplit(names(unlist(lapply(l, function(x){
  p = x$`log-contrast coefficients`
  names(p) = x$taxa.name
  p
}))), ".", fixed=T), function(y) paste0(y[2:length(y)], collapse = "_")), var))
findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}
df.CODA2.effect.var$Var1 = gsub("_"," ", df.CODA2.effect.var$Var1)
df.CODA2.effect.var = df.CODA2.effect.var %>% mutate(outlier = if_else(findoutlier(value), Var1, NA))
df.CODA2.effect.var$Species = "Microorganisms"

df.CODA.example = data.frame("Coefficient" = l$X100010869$`log-contrast coefficients`,"Species"=l$X100010869$taxa.name)
df.CODA.example$Species = gsub("\\."," ", df.CODA.example$Species)
df.CODA.example$sign = ifelse(df.CODA.example$Coefficient<0, "-", "+")


ggarrange(ggarrange(ggplot(top10.Metabo, aes(x = reorder(Chemical,-abs(Factor1.scaled)), y = abs(Factor1.scaled), fill=sign, label=sign)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
                      ylab("Weights")+xlab("Metabolite")+geom_label(size=3.5)+ylim(c(0,1))
                    ,ggplot(top10.Microbio, aes(x = reorder(Species,-abs(Factor1.scaled)), y = abs(Factor1.scaled), fill=sign, label=sign)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
                      ylab("Weights")+xlab("Microorganism")+geom_label(size=3.5)+ylim(c(0,1))
,labels = c("A","B")),ggarrange(ggplot(df.sPLS.KONZO.loadings, aes(x=X1st.Component, y=X2nd.Component, colour=Type))+geom_point()+theme_bw()+scale_color_nejm()+xlab("1st Component")+ylab("2nd Component"),ggplot(df.CODA2.nasso, aes(x=as.factor(Var1), y=value))+geom_bar(stat="identity", position="dodge",colour="black",fill="#2b8cbe")+theme_bw()+xlab("#Associated Microorganisms")+ylab("Frequency"), labels=c("C","D"),ncol=2),
ggarrange(ggplot(df.CODA.example, aes(x = reorder(Species,-abs(Coefficient)), y = abs(Coefficient), fill=sign, label=sign)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
            ylab("Coefficients")+xlab("Microorganism")+geom_label(size=3.5)+ylim(c(0,1)),
          ggplot(df.CODA2.effect.var, aes(x=Species, y=value)) +
            geom_violin(width=0.3,fill="#9ecae1",color="black") + theme_bw() +xlab("")+ ylab("Variance of Coefficients")+
            ggrepel::geom_text_repel(data = df.CODA2.effect.var[!is.na(df.CODA2.effect.var$outlier), ], aes(y=value,label=outlier, colour="red"), show.legend = F), ncol=2, labels=c("E","F")),nrow=3)


models_lc_subset_sPLS = readRDS("data\\Konzo\\results_LCLR_Konzo_subset_sPLS.RDS")


pvalues.effects.lc.subset.sPLS = sapply(1:length(models_lc_subset_sPLS), function(x) models_lc_subset_sPLS[[x]]$pvalues)
signi.models.lc.subset.sPLS = which(pvalues.effects.lc.subset.sPLS<=0.05/249)
subset.signi.models_lc_subset_sPLS = models_lc_subset_sPLS[signi.models.lc.subset.sPLS]


df.metabo.species.LC = do.call("rbind",lapply(1:3, function(x) {
  t = data.frame("from"=rownames(non.null.loadings.metabo.KONZO.sPLS)[signi.models.lc.subset.sPLS[x]],"Weight"=subset.signi.models_lc_subset_sPLS[[x]]$coefs[-1])
  t$to = rownames(t)
  rownames(t) =NULL
  t = t[,c("from", "to", "Weight")]
  t
}))

graph.metabo.species.LC = igraph::graph_from_data_frame(df.metabo.species.LC, vertices = c(unique(df.metabo.species.LC$from), unique(df.metabo.species.LC$to)), directed = FALSE)
E(graph.metabo.species.LC)$width = abs(df.metabo.species.LC$Weight)
E(graph.metabo.species.LC)$Color = ifelse(df.metabo.species.LC$Weight<0, "green", "red")
V(graph.metabo.species.LC)$Type = ifelse(names(V(graph.metabo.species.LC))%in%unique(df.metabo.species.LC$from), "Metabolites", "Species" )
V(graph.metabo.species.LC)$Color = ifelse(V(graph.metabo.species.LC)$Type=="Metabolites", "pink", "brown" )

col = c("Metabolites" = "pink", "Species" = "gold")
ggnet::ggnet2(graph.metabo.species.LC, color= "Type", color.palette =col ,size=4,edge.size="width", edge.color = "Color", color.legend = "Type")


plot(graph.metabo.species.LC, vertex.label=NA, vertex.size=7, vertex.color= V(graph.metabo.species.LC)$Color, edge.width = E(graph.metabo.species.LC)$width*1.5, edge.color=E(graph.metabo.species.LC)$Color)
legend(x=-1, y=-1.1, c("Metabolites","Species"), pch=21,
       
       col="#777777", pt.bg=c("pink", "brown"), pt.cex=2, cex=.8, bty="n", ncol=1)

df.metabo.species.CODA = do.call("rbind",lapply(1:length(l), function(x){
  if(length(l[[x]]$taxa.name)>0){
    data.frame("from"=names(l)[x], "to" = l[[x]]$taxa.name, "Weight" = l[[x]]$`log-contrast coefficients` )
    
  }
}))

graph.metabo.species.CODA = igraph::graph_from_data_frame(df.metabo.species.CODA, vertices = c(unique(df.metabo.species.CODA$from), unique(df.metabo.species.CODA$to)), directed = FALSE)
E(graph.metabo.species.CODA)$width = abs(df.metabo.species.CODA$Weight)
E(graph.metabo.species.CODA)$Color = ifelse(df.metabo.species.CODA$Weight<0, "green", "red")
V(graph.metabo.species.CODA)$Type = ifelse(names(V(graph.metabo.species.CODA))%in%unique(df.metabo.species.CODA$from), "Metabolites", "Species" )
V(graph.metabo.species.CODA)$VColor = ifelse(V(graph.metabo.species.CODA)$Type=="Metabolites", "pink", "brown" )

col = c("Metabolites" = "pink", "Species" = "brown")
ggnet::ggnet2(graph.metabo.species.CODA, color= "Type", color.palette =col ,size=4,edge.size="width", edge.color = "Color", color.legend = "Type")

igraph::degree(graph.metabo.species.CODA)[names(igraph::degree(graph.metabo.species.CODA)) == "Faecalibacterium.prausnitzii"]
sort(igraph::degree(graph.metabo.species.CODA),decreasing=T)[1:100]

cluster_graph_metabo_species_CODA = igraph::cluster_louvain(graph.metabo.species.CODA)

par(mfrow=c(1,1))
plot(graph.metabo.species.CODA,vertex.label=NA, vertex.size=5, vertex.color= V(graph.metabo.species.CODA)$VColor, edge.width = E(graph.metabo.species.CODA)$width*1.5, edge.color=E(graph.metabo.species.CODA)$Color)
 legend(x=-1, y=-0.8, c("Metabolites","Species"), pch=21,
       
       col="#777777", pt.bg=c("pink", "brown"), pt.cex=2, cex=.8, bty="n", ncol=1)

plot(cluster_graph_metabo_species_CODA, graph.metabo.species.CODA, vertex.label=NA, vertex.size=4)
