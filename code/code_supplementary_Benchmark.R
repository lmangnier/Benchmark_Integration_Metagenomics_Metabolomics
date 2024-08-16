#Script for generating supplementary figures
library(ggplot2)
library(ggpubr)
library(wesanderson)

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

gg_qqplot_facet_grid <- function(list_pvalues,ci = 0.95, title="") {
  #n  <- length(ps)
  
  dfs <- lapply(1:length(list_pvalues), function(i) {
    n = length(list_pvalues[[i]])
    data.frame(
      Distance_Kernel = names(list_pvalues)[i],
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
    geom_point(aes(expected, observed, color=Distance_Kernel), size = 2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5)+xlab(log10Pe) +
    ylab(log10Po)+theme_bw()#+theme(legend.position = "none")
  
}

power_figure = function(df){
  
  ggplot(df, aes(x=Distance_Kernel, y=Power, fill=Distance_Kernel)) + geom_bar(stat="identity", position="dodge", colour="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))
}
#Global Association methods
mantel.ilr.spearman.100.100.500 = list.files("data\\Mantel\\ilr\\100_100_500\\spearman", full.names = TRUE)
mantel.clr.spearman.100.100.500 = list.files("data\\Mantel\\clr\\100_100_500\\spearman", full.names = TRUE)
mantel.alpha.spearman.100.100.500 = list.files("data\\Mantel\\alpha\\100_100_500\\spearman", full.names = TRUE)

mantel.ilr.spearman.25.25.100 = list.files("data\\Mantel\\ilr\\25_25_100", full.names = TRUE)
mantel.clr.spearman.25.25.100 = list.files("data\\Mantel\\clr\\25_25_100\\spearman", full.names = TRUE)
mantel.alpha.spearman.25.25.100 = list.files("data\\Mantel\\alpha\\25_25_100\\spearman", full.names = TRUE)


mantel.ilr.pearson.100.100.500 = list.files("data\\Mantel\\ilr\\100_100_500\\pearson", full.names = TRUE)
mantel.clr.pearson.100.100.500 = list.files("data\\Mantel\\clr\\100_100_500\\pearson", full.names = TRUE)
mantel.alpha.pearson.100.100.500 = list.files("data\\Mantel\\alpha\\100_100_500\\pearson", full.names = TRUE)

mantel.ilr.pearson.25.25.100 = list.files("data\\Mantel\\ilr\\25_25_100\\pearson", full.names = TRUE)
mantel.clr.pearson.25.25.100 = list.files("data\\Mantel\\clr\\25_25_100\\pearson", full.names = TRUE)
mantel.alpha.pearson.25.25.100 = list.files("data\\Mantel\\alpha\\25_25_100\\pearson", full.names = TRUE)


list.pvalues.null.mantel.spearman.clr.log.100.100.500 = list("Canberra"=readRDS(mantel.clr.spearman.100.100.500[7]),
                                                              "Euclidean"=readRDS(mantel.clr.spearman.100.100.500[8]),
                                                              "Manhattan"=readRDS(mantel.clr.spearman.100.100.500[9]))



list.pvalues.null.mantel.spearman.alpha.log.100.100.500 = list("Canberra"=readRDS(mantel.alpha.spearman.100.100.500[7]),
                                                              "Euclidean"=readRDS(mantel.alpha.spearman.100.100.500[8]),
                                                              "Manhattan"=readRDS(mantel.alpha.spearman.100.100.500[9]))

df.power.Mantel.spearman.clr.log.100.100.500 = data.frame("Power" = c(sum(readRDS(mantel.clr.spearman.100.100.500[1])<=0.05)/1000,
                                                              sum(readRDS(mantel.clr.spearman.100.100.500[2])<=0.05)/1000,
                                                              sum(readRDS(mantel.clr.spearman.100.100.500[3])<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.Mantel.spearman.alpha.log.100.100.500 = data.frame("Power" = c(sum(readRDS(mantel.alpha.spearman.100.100.500[1])<=0.05)/1000,
                                                                        sum(readRDS(mantel.alpha.spearman.100.100.500[2])<=0.05)/1000,
                                                                        sum(readRDS(mantel.alpha.spearman.100.100.500[3])<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))



ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.clr.log.100.100.500)+ggtitle("CLR")+labs(color="Distance Kernel"),
          gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.alpha.log.100.100.500)+ggtitle("Alpha")+labs(color="Distance Kernel"),
          ggplot(df.power.Mantel.spearman.clr.log.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.Mantel.spearman.alpha.log.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),labels = c("A","B","C","D"), nrow=2, ncol=2)


list.pvalues.null.mantel.pearson.clr.log.100.100.500 = list("Canberra"=readRDS(mantel.clr.pearson.100.100.500[7]),
                                                             "Euclidean"=readRDS(mantel.clr.pearson.100.100.500[8]),
                                                             "Manhattan"=readRDS(mantel.clr.pearson.100.100.500[9]))

list.pvalues.null.mantel.pearson.ilr.log.100.100.500 = list("Canberra"=readRDS(mantel.ilr.pearson.100.100.500[7]),
                                                              "Euclidean"=readRDS(mantel.ilr.pearson.100.100.500[8]),
                                                              "Manhattan"=readRDS(mantel.ilr.pearson.100.100.500[9]))

list.pvalues.null.mantel.pearson.alpha.log.100.100.500 = list("Canberra"=readRDS(mantel.alpha.pearson.100.100.500[7]),
                                                               "Euclidean"=readRDS(mantel.alpha.pearson.100.100.500[8]),
                                                               "Manhattan"=readRDS(mantel.alpha.pearson.100.100.500[9]))

df.power.Mantel.pearson.clr.log.100.100.500 = data.frame("Power" = c(sum(readRDS(mantel.clr.pearson.100.100.500[1])<=0.05)/1000,
                                                                      sum(readRDS(mantel.clr.pearson.100.100.500[2])<=0.05)/1000,
                                                                      sum(readRDS(mantel.clr.pearson.100.100.500[3])<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.Mantel.pearson.ilr.log.100.100.500 = data.frame("Power" = c(sum(readRDS(mantel.ilr.pearson.100.100.500[1])<=0.05)/1000,
                                                                     sum(readRDS(mantel.ilr.pearson.100.100.500[2])<=0.05)/1000,
                                                                     sum(readRDS(mantel.ilr.pearson.100.100.500[3])<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.Mantel.pearson.alpha.log.100.100.500 = data.frame("Power" = c(sum(readRDS(mantel.alpha.pearson.100.100.500[1])<=0.05)/1000,
                                                                        sum(readRDS(mantel.alpha.pearson.100.100.500[2])<=0.05)/1000,
                                                                        sum(readRDS(mantel.alpha.pearson.100.100.500[3])<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))



ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.clr.log.100.100.500)+ggtitle("CLR")+labs(color="Distance Kernel"),
          gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.ilr.log.100.100.500)+ggtitle("ILR")+labs(color="Distance Kernel"),
          gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.alpha.log.100.100.500)+ggtitle("Alpha")+labs(color="Distance Kernel"),
          ggplot(df.power.Mantel.pearson.clr.log.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.Mantel.pearson.ilr.log.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.Mantel.pearson.alpha.log.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),labels = c("A","B","C","D","E","F"), nrow=2, ncol=3)

list.pvalues.null.mantel.spearman.ilr.log.25.25.100 = list("Canberra"=readRDS(mantel.ilr.spearman.25.25.100[7]),
                                                             "Euclidean"=readRDS(mantel.ilr.spearman.25.25.100[8]),
                                                             "Manhattan"=readRDS(mantel.ilr.spearman.25.25.100[9]))




list.pvalues.null.mantel.spearman.clr.log.100.100.500 = list("Canberra"=readRDS(mantel.clr.spearman.100.100.500[7]),
                                                             "Euclidean"=readRDS(mantel.clr.spearman.100.100.500[8]),
                                                             "Manhattan"=readRDS(mantel.clr.spearman.100.100.500[9]))



list.pvalues.null.mantel.spearman.alpha.log.100.100.500 = list("Canberra"=readRDS(mantel.alpha.spearman.100.100.500[7]),
                                                               "Euclidean"=readRDS(mantel.alpha.spearman.100.100.500[8]),
                                                               "Manhattan"=readRDS(mantel.alpha.spearman.100.100.500[9]))

df.power.Mantel.spearman.clr.log.100.100.500 = data.frame("Power" = c(sum(readRDS(mantel.clr.spearman.100.100.500[1])<=0.05)/1000,
                                                                      sum(readRDS(mantel.clr.spearman.100.100.500[2])<=0.05)/1000,
                                                                      sum(readRDS(mantel.clr.spearman.100.100.500[3])<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.Mantel.spearman.alpha.log.100.100.500 = data.frame("Power" = c(sum(readRDS(mantel.alpha.spearman.100.100.500[1])<=0.05)/1000,
                                                                        sum(readRDS(mantel.alpha.spearman.100.100.500[2])<=0.05)/1000,
                                                                        sum(readRDS(mantel.alpha.spearman.100.100.500[3])<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))



ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.clr.log.100.100.500)+ggtitle("CLR")+labs(color="Distance Kernel"),
          gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.alpha.log.100.100.500)+ggtitle("Alpha")+labs(color="Distance Kernel"),
          ggplot(df.power.Mantel.spearman.clr.log.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.Mantel.spearman.alpha.log.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),labels = c("A","B","C","D"), nrow=2, ncol=2)
list.pvalues.null.mantel.spearman.ilr.none.100.100.500 = list("Canberra"=readRDS(mantel.ilr.spearman.100.100.500[10]),
"Euclidean"=readRDS(mantel.ilr.spearman.100.100.500[11]),
"Manhattan"=readRDS(mantel.ilr.spearman.100.100.500[12]))


gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.ilr.none.100.100.500)


list.pvalues.null.mantel.spearman.clr.100.100.500 = list("Canberra"=readRDS(mantel.clr.spearman.100.100.500[10]),
                                                         "Euclidean"=readRDS(mantel.clr.spearman.100.100.500[11]),
                                                         "Manhattan"=readRDS(mantel.clr.spearman.100.100.500[12]))


list.pvalues.null.mantel.spearman.ilr.100.100.500 = list("Canberra"=readRDS(mantel.ilr.spearman.100.100.500[10]),
                                                             "Euclidean"=readRDS(mantel.ilr.spearman.100.100.500[11]),
                                                             "Manhattan"=readRDS(mantel.ilr.spearman.100.100.500[12]))



list.pvalues.null.mantel.spearman.alpha.100.100.500 = list("Canberra"=readRDS(mantel.alpha.spearman.100.100.500[10]),
                                                               "Euclidean"=readRDS(mantel.alpha.spearman.100.100.500[11]),
                                                               "Manhattan"=readRDS(mantel.alpha.spearman.100.100.500[12]))


df.power.Mantel.spearman.ilr.100.100.500 = data.frame("Power"=c(sum(readRDS(mantel.ilr.spearman.100.100.500[4])<=0.05)/1000,
                                                                sum(readRDS(mantel.ilr.spearman.100.100.500[5])<=0.05)/1000,
                                                                sum(readRDS(mantel.ilr.spearman.100.100.500[6])<=0.05)/1000),
                                                      "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

df.power.Mantel.spearman.alpha.100.100.500 = data.frame("Power"=c(sum(readRDS(mantel.alpha.spearman.100.100.500[4])<=0.05)/1000,
                                                                sum(readRDS(mantel.alpha.spearman.100.100.500[5])<=0.05)/1000,
                                                                sum(readRDS(mantel.alpha.spearman.100.100.500[6])<=0.05)/1000),
                                                      "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="Alpha")
df.power.Mantel.spearman.clr.100.100.500 = data.frame("Power"=c(sum(readRDS(mantel.clr.spearman.100.100.500[4])<=0.05)/1000,
                                                                sum(readRDS(mantel.clr.spearman.100.100.500[5])<=0.05)/1000,
                                                                sum(readRDS(mantel.clr.spearman.100.100.500[6])<=0.05)/1000),
                                                      "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="CLR")

ggarrange(ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.clr.100.100.500)+labs(colour="Distance Kernel")+ggtitle("CLR"),
          gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.ilr.100.100.500)+labs(colour="Distance Kernel")+ggtitle("ILR"),
          gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.alpha.100.100.500)+labs(colour="Distance Kernel")+ggtitle("Alpha"), ncol=3, labels = c("A","B","C")), 
          ggarrange(power_figure(df.power.Mantel.spearman.clr.100.100.500), power_figure(df.power.Mantel.spearman.ilr.100.100.500), power_figure(df.power.Mantel.spearman.alpha.100.100.500), ncol=3,labels=c("D", "E","F")),nrow=2)

list.pvalues.null.mantel.pearson.clr.100.100.500 = list("Canberra"=readRDS(mantel.clr.pearson.100.100.500[10]),
                                                         "Euclidean"=readRDS(mantel.clr.pearson.100.100.500[11]),
                                                         "Manhattan"=readRDS(mantel.clr.pearson.100.100.500[12]))


list.pvalues.null.mantel.pearson.ilr.100.100.500 = list("Canberra"=readRDS(mantel.ilr.pearson.100.100.500[10]),
                                                         "Euclidean"=readRDS(mantel.ilr.pearson.100.100.500[11]),
                                                         "Manhattan"=readRDS(mantel.ilr.pearson.100.100.500[12]))



list.pvalues.null.mantel.pearson.alpha.100.100.500 = list("Canberra"=readRDS(mantel.alpha.pearson.100.100.500[10]),
                                                           "Euclidean"=readRDS(mantel.alpha.pearson.100.100.500[11]),
                                                           "Manhattan"=readRDS(mantel.alpha.pearson.100.100.500[12]))


df.power.Mantel.pearson.ilr.100.100.500 = data.frame("Power"=c(sum(readRDS(mantel.ilr.pearson.100.100.500[4])<=0.05)/1000,
                                                                sum(readRDS(mantel.ilr.pearson.100.100.500[5])<=0.05)/1000,
                                                                sum(readRDS(mantel.ilr.pearson.100.100.500[6])<=0.05)/1000),
                                                      "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

df.power.Mantel.pearson.alpha.100.100.500 = data.frame("Power"=c(sum(readRDS(mantel.alpha.pearson.100.100.500[4])<=0.05)/1000,
                                                                  sum(readRDS(mantel.alpha.pearson.100.100.500[5])<=0.05)/1000,
                                                                  sum(readRDS(mantel.alpha.pearson.100.100.500[6])<=0.05)/1000),
                                                        "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="Alpha")
df.power.Mantel.pearson.clr.100.100.500 = data.frame("Power"=c(sum(readRDS(mantel.clr.pearson.100.100.500[4])<=0.05)/1000,
                                                                sum(readRDS(mantel.clr.pearson.100.100.500[5])<=0.05)/1000,
                                                                sum(readRDS(mantel.clr.pearson.100.100.500[6])<=0.05)/1000),
                                                      "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="CLR")

ggarrange(ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.clr.100.100.500)+labs(colour="Distance Kernel")+ggtitle("CLR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.ilr.100.100.500)+labs(colour="Distance Kernel")+ggtitle("ILR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.alpha.100.100.500)+labs(colour="Distance Kernel")+ggtitle("Alpha"), ncol=3, labels = c("A","B","C")), 
          ggarrange(power_figure(df.power.Mantel.pearson.clr.100.100.500), power_figure(df.power.Mantel.pearson.ilr.100.100.500), power_figure(df.power.Mantel.pearson.alpha.100.100.500), ncol=3,labels=c("D", "E","F")),nrow=2)


list.pvalues.null.mantel.spearman.clr.log.25.25.100 = list("Canberra"=readRDS(mantel.clr.spearman.25.25.100[7]),
                                                         "Euclidean"=readRDS(mantel.clr.spearman.25.25.100[8]),
                                                         "Manhattan"=readRDS(mantel.clr.spearman.25.25.100[9]))


list.pvalues.null.mantel.spearman.ilr.log.25.25.100 = list("Canberra"=readRDS(mantel.ilr.spearman.25.25.100[7]),
                                                         "Euclidean"=readRDS(mantel.ilr.spearman.25.25.100[8]),
                                                         "Manhattan"=readRDS(mantel.ilr.spearman.25.25.100[9]))



list.pvalues.null.mantel.spearman.alpha.log.25.25.100 = list("Canberra"=readRDS(mantel.alpha.spearman.25.25.100[7]),
                                                           "Euclidean"=readRDS(mantel.alpha.spearman.25.25.100[8]),
                                                           "Manhattan"=readRDS(mantel.alpha.spearman.25.25.100[9]))


df.power.Mantel.spearman.ilr.log.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.ilr.spearman.25.25.100[1])<=0.05)/1000,
                                                                sum(readRDS(mantel.ilr.spearman.25.25.100[2])<=0.05)/1000,
                                                                sum(readRDS(mantel.ilr.spearman.25.25.100[3])<=0.05)/1000),
                                                      "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

df.power.Mantel.spearman.alpha.log.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.alpha.spearman.25.25.100[1])<=0.05)/1000,
                                                                  sum(readRDS(mantel.alpha.spearman.25.25.100[2])<=0.05)/1000,
                                                                  sum(readRDS(mantel.alpha.spearman.25.25.100[3])<=0.05)/1000),
                                                        "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="Alpha")
df.power.Mantel.spearman.clr.log.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.clr.spearman.25.25.100[1])<=0.05)/1000,
                                                                sum(readRDS(mantel.clr.spearman.25.25.100[2])<=0.05)/1000,
                                                                sum(readRDS(mantel.clr.spearman.25.25.100[3])<=0.05)/1000),
                                                      "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="CLR")

ggarrange(ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.clr.log.25.25.100)+labs(colour="Distance Kernel")+ggtitle("CLR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.ilr.log.25.25.100)+labs(colour="Distance Kernel")+ggtitle("ILR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.alpha.log.25.25.100)+labs(colour="Distance Kernel")+ggtitle("Alpha"), ncol=3, labels = c("A","B","C")), 
          ggarrange(power_figure(df.power.Mantel.spearman.clr.log.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.spearman.ilr.log.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.spearman.alpha.log.25.25.100)+xlab("Distance Kernel"), ncol=3,labels=c("D", "E","F")),nrow=2)



list.pvalues.null.mantel.pearson.clr.log.25.25.100 = list("Canberra"=readRDS(mantel.clr.pearson.25.25.100[7]),
                                                           "Euclidean"=readRDS(mantel.clr.pearson.25.25.100[8]),
                                                           "Manhattan"=readRDS(mantel.clr.pearson.25.25.100[9]))


list.pvalues.null.mantel.pearson.ilr.log.25.25.100 = list("Canberra"=readRDS(mantel.ilr.pearson.25.25.100[9]),
                                                           "Euclidean"=readRDS(mantel.ilr.pearson.25.25.100[10]),
                                                           "Manhattan"=readRDS(mantel.ilr.pearson.25.25.100[11]))



list.pvalues.null.mantel.pearson.alpha.log.25.25.100 = list("Canberra"=readRDS(mantel.alpha.pearson.25.25.100[7]),
                                                             "Euclidean"=readRDS(mantel.alpha.pearson.25.25.100[8]),
                                                             "Manhattan"=readRDS(mantel.alpha.pearson.25.25.100[9]))


df.power.Mantel.pearson.ilr.log.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.ilr.pearson.25.25.100[3])<=0.05)/1000,
                                                                  sum(readRDS(mantel.ilr.pearson.25.25.100[4])<=0.05)/1000,
                                                                  sum(readRDS(mantel.ilr.pearson.25.25.100[5])<=0.05)/1000),
                                                        "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

df.power.Mantel.pearson.alpha.log.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.alpha.pearson.25.25.100[1])<=0.05)/1000,
                                                                    sum(readRDS(mantel.alpha.pearson.25.25.100[2])<=0.05)/1000,
                                                                    sum(readRDS(mantel.alpha.pearson.25.25.100[3])<=0.05)/1000),
                                                          "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="Alpha")
df.power.Mantel.pearson.clr.log.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.clr.pearson.25.25.100[1])<=0.05)/1000,
                                                                  sum(readRDS(mantel.clr.pearson.25.25.100[2])<=0.05)/1000,
                                                                  sum(readRDS(mantel.clr.pearson.25.25.100[3])<=0.05)/1000),
                                                        "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="CLR")

ggarrange(ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.clr.log.25.25.100)+labs(colour="Distance Kernel")+ggtitle("CLR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.ilr.log.25.25.100)+labs(colour="Distance Kernel")+ggtitle("ILR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.alpha.log.25.25.100)+labs(colour="Distance Kernel")+ggtitle("Alpha"), ncol=3, labels = c("A","B","C")), 
          ggarrange(power_figure(df.power.Mantel.pearson.clr.log.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.pearson.ilr.log.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.pearson.alpha.log.25.25.100)+xlab("Distance Kernel"), ncol=3,labels=c("D", "E","F")),nrow=2)


list.pvalues.null.mantel.spearman.clr.25.25.100 = list("Canberra"=readRDS(mantel.clr.spearman.25.25.100[10]),
                                                           "Euclidean"=readRDS(mantel.clr.spearman.25.25.100[11]),
                                                           "Manhattan"=readRDS(mantel.clr.spearman.25.25.100[12]))


list.pvalues.null.mantel.spearman.ilr.25.25.100 = list("Canberra"=readRDS(mantel.ilr.spearman.25.25.100[10]),
                                                           "Euclidean"=readRDS(mantel.ilr.spearman.25.25.100[11]),
                                                           "Manhattan"=readRDS(mantel.ilr.spearman.25.25.100[12]))



list.pvalues.null.mantel.spearman.alpha.25.25.100 = list("Canberra"=readRDS(mantel.alpha.spearman.25.25.100[10]),
                                                             "Euclidean"=readRDS(mantel.alpha.spearman.25.25.100[11]),
                                                             "Manhattan"=readRDS(mantel.alpha.spearman.25.25.100[12]))


df.power.Mantel.spearman.ilr.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.ilr.spearman.25.25.100[4])<=0.05)/1000,
                                                                  sum(readRDS(mantel.ilr.spearman.25.25.100[5])<=0.05)/1000,
                                                                  sum(readRDS(mantel.ilr.spearman.25.25.100[6])<=0.05)/1000),
                                                        "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

df.power.Mantel.spearman.alpha.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.alpha.spearman.25.25.100[4])<=0.05)/1000,
                                                                    sum(readRDS(mantel.alpha.spearman.25.25.100[5])<=0.05)/1000,
                                                                    sum(readRDS(mantel.alpha.spearman.25.25.100[6])<=0.05)/1000),
                                                          "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="Alpha")
df.power.Mantel.spearman.clr.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.clr.spearman.25.25.100[4])<=0.05)/1000,
                                                                  sum(readRDS(mantel.clr.spearman.25.25.100[5])<=0.05)/1000,
                                                                  sum(readRDS(mantel.clr.spearman.25.25.100[6])<=0.05)/1000),
                                                        "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="CLR")

ggarrange(ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.clr.25.25.100)+labs(colour="Distance Kernel")+ggtitle("CLR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.ilr.25.25.100)+labs(colour="Distance Kernel")+ggtitle("ILR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.spearman.alpha.25.25.100)+labs(colour="Distance Kernel")+ggtitle("Alpha"), ncol=3, labels = c("A","B","C")), 
          ggarrange(power_figure(df.power.Mantel.spearman.clr.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.spearman.ilr.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.spearman.alpha.25.25.100)+xlab("Distance Kernel"), ncol=3,labels=c("D", "E","F")),nrow=2)



list.pvalues.null.mantel.pearson.clr.25.25.100 = list("Canberra"=readRDS(mantel.clr.pearson.25.25.100[10]),
                                                       "Euclidean"=readRDS(mantel.clr.pearson.25.25.100[11]),
                                                       "Manhattan"=readRDS(mantel.clr.pearson.25.25.100[12]))


list.pvalues.null.mantel.pearson.ilr.25.25.100 = list("Canberra"=readRDS(mantel.ilr.pearson.25.25.100[12]),
                                                       "Euclidean"=readRDS(mantel.ilr.pearson.25.25.100[13]),
                                                       "Manhattan"=readRDS(mantel.ilr.pearson.25.25.100[14]))



list.pvalues.null.mantel.pearson.alpha.25.25.100 = list("Canberra"=readRDS(mantel.alpha.pearson.25.25.100[11]),
                                                         "Euclidean"=readRDS(mantel.alpha.pearson.25.25.100[12]),
                                                         "Manhattan"=readRDS(mantel.alpha.pearson.25.25.100[13]))


df.power.Mantel.pearson.ilr.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.ilr.pearson.25.25.100[6])<=0.05)/1000,
                                                              sum(readRDS(mantel.ilr.pearson.25.25.100[7])<=0.05)/1000,
                                                              sum(readRDS(mantel.ilr.pearson.25.25.100[8])<=0.05)/1000),
                                                    "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

df.power.Mantel.pearson.alpha.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.alpha.pearson.25.25.100[4])<=0.05)/1000,
                                                                sum(readRDS(mantel.alpha.pearson.25.25.100[5])<=0.05)/1000,
                                                                sum(readRDS(mantel.alpha.pearson.25.25.100[6])<=0.05)/1000),
                                                      "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="Alpha")

df.power.Mantel.pearson.clr.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.clr.pearson.25.25.100[4])<=0.05)/1000,
                                                              sum(readRDS(mantel.clr.pearson.25.25.100[5])<=0.05)/1000,
                                                              sum(readRDS(mantel.clr.pearson.25.25.100[6])<=0.05)/1000),
                                                    "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"),"Normalization"="CLR")

ggarrange(ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.clr.25.25.100)+labs(colour="Distance Kernel")+ggtitle("CLR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.ilr.25.25.100)+labs(colour="Distance Kernel")+ggtitle("ILR"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.alpha.25.25.100)+labs(colour="Distance Kernel")+ggtitle("Alpha"), ncol=3, labels = c("A","B","C")), 
          ggarrange(power_figure(df.power.Mantel.pearson.clr.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.pearson.ilr.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.pearson.alpha.25.25.100)+xlab("Distance Kernel"), ncol=3,labels=c("D", "E","F")),nrow=2)

mantel.original.25.25.100 = list.files("C:/Users/loicm/Desktop/Projets/Benchmark_Microbiome_Metabolome/results/Mantel/original/25_25_100/", full.names = TRUE)
mantel.original.100.100.500 = list.files("C:/Users/loicm/Desktop/Projets/Benchmark_Microbiome_Metabolome/results/Mantel/original/100_100_500/", full.names = TRUE)

list.pvalues.null.mantel.original.original.25.25.100 = list("Canberra"=readRDS("C:/Users/loicm/Desktop/Projets/Benchmark_Microbiome_Metabolome/results/Mantel/original/25_25_100/pvalues_null_mantel_spearman_microbiome_metabolome_euclidean_canberra_25_25_100indiv.RDS"),
                                                       "Euclidean"=readRDS("C:/Users/loicm/Desktop/Projets/Benchmark_Microbiome_Metabolome/results/Mantel/original/25_25_100/pvalues_null_mantel_spearman_microbiome_metabolome_euclidean_euclidean_25_25_100indiv.RDS"),
                                                       "Manhattan"=readRDS("C:/Users/loicm/Desktop/Projets/Benchmark_Microbiome_Metabolome/results/Mantel/original/25_25_100/pvalues_null_mantel_spearman_microbiome_metabolome_euclidean_manhattan_25_25_100indiv.RDS"))


list.pvalues.null.mantel.original.original.100.100.100 = list("Canberra"=readRDS(mantel.original.100.100.500[22]),
                                                            "Euclidean"=readRDS(mantel.original.100.100.500[23]),
                                                            "Manhattan"=readRDS(mantel.original.100.100.500[24]))


df.power.Mantel.spearman.original.original.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.original.25.25.100[10])<=0.05)/1000,
                                                              sum(readRDS(mantel.original.25.25.100[11])<=0.05)/1000,
                                                              sum(readRDS(mantel.original.25.25.100[12])<=0.05)/1000),
                                                    "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

df.power.Mantel.spearman.original.original.100.100.500 = data.frame("Power"=c(sum(readRDS(mantel.original.100.100.500[10])<=0.05)/1000,
                                                                            sum(readRDS(mantel.original.100.100.500[11])<=0.05)/1000,
                                                                            sum(readRDS(mantel.original.100.100.500[12])<=0.05)/1000),
                                                                  "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

ggarrange(ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.original.original.25.25.100)+labs(colour="Distance Kernel")+ggtitle("Low-Dimensional"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.original.original.100.100.100)+labs(colour="Distance Kernel")+ggtitle("High-Dimensional"), ncol=2, labels = c("A","B")), 
          ggarrange(power_figure(df.power.Mantel.spearman.original.original.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.spearman.original.original.100.100.500)+xlab("Distance Kernel"), ncol=2,labels=c("C","D")),nrow=2)


mantel.pearson.original.25.25.100 = list.files("C:/Users/loicm/Desktop/Projets/Benchmark_Microbiome_Metabolome/results/Mantel/original/25_25_100/pearson/", full.names = TRUE)
mantel.pearson.original.100.100.500 = list.files("C:/Users/loicm/Desktop/Projets/Benchmark_Microbiome_Metabolome/results/Mantel/original/100_100_500/pearson/", full.names = TRUE)

list.pvalues.null.mantel.pearson.original.original.25.25.100 = list("Canberra"=readRDS(mantel.pearson.original.25.25.100[4]),
                                                            "Euclidean"=readRDS(mantel.pearson.original.25.25.100[5]),
                                                            "Manhattan"=readRDS(mantel.pearson.original.25.25.100[6]))


list.pvalues.null.mantel.pearson.original.original.100.100.100 = list("Canberra"=readRDS(mantel.pearson.original.100.100.500[10]),
                                                              "Euclidean"=readRDS(mantel.pearson.original.100.100.500[11]),
                                                              "Manhattan"=readRDS(mantel.pearson.original.100.100.500[12]))


df.power.Mantel.pearson.original.original.25.25.100 = data.frame("Power"=c(sum(readRDS(mantel.pearson.original.25.25.100[1])<=0.05)/1000,
                                                                            sum(readRDS(mantel.pearson.original.25.25.100[2])<=0.05)/1000,
                                                                            sum(readRDS(mantel.pearson.original.25.25.100[3])<=0.05)/1000),
                                                                  "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

df.power.Mantel.pearson.original.original.100.100.500 = data.frame("Power"=c(sum(readRDS(mantel.pearson.original.100.100.500[4])<=0.05)/1000,
                                                                              sum(readRDS(mantel.pearson.original.100.100.500[5])<=0.05)/1000,
                                                                              sum(readRDS(mantel.pearson.original.100.100.500[6])<=0.05)/1000),
                                                                    "Distance_Kernel"=c("Canberra","Euclidean","Manhattan"), "Normalization"="ILR")

ggarrange(ggarrange(gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.original.original.25.25.100)+labs(colour="Distance Kernel")+ggtitle("Low-Dimensional"),
                    gg_qqplot_facet_grid(list.pvalues.null.mantel.pearson.original.original.100.100.100)+labs(colour="Distance Kernel")+ggtitle("High-Dimensional"), ncol=2, labels = c("A","B")), 
          ggarrange(power_figure(df.power.Mantel.pearson.original.original.25.25.100)+xlab("Distance Kernel"), power_figure(df.power.Mantel.pearson.original.original.100.100.500)+xlab("Distance Kernel"), ncol=2,labels=c("C","D")),nrow=2)

MMiRKAT.alpha.25.25.100 = list.files("data\\MMiRKAT\\alpha\\25.25.100", full.names=TRUE)
MMiRKAT.alpha.100.100.500 = list.files("data\\MMiRKAT\\alpha\\100_100_500", full.names=TRUE)

pvalues.null.MMiRKAT.100.100.500 = readRDS("data\\MMiRKAT\\list.results.null.MMiRKAT.compositional.100.100.500.RDS")
pvalues.alter.MMiRKAT.100.100.500 = readRDS("data\\MMiRKAT\\list.results.power.MMiRKAT.compositional.100.100.500.RDS")

pvalues.null.MMiRKAT.25.25.100 = readRDS("data\\MMiRKAT\\list.results.null.MMiRKAT.compositional.25.25.100.RDS")
pvalues.alter.MMiRKAT.25.25.100 = readRDS("data\\MMiRKAT\\list.results.power.MMiRKAT.compositional.25.25.100.RDS")


list.pvalues.null.MMiRKAT.clr.log.25.25.100 = list("Euclidean"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.clr.log$null.p.values.25.25.100indiv.MMiRKAT.Euclidean.clr.log,
                                                   "Canberra"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.clr.log$null.p.values.25.25.100indiv.MMiRKAT.Canberra.clr.log,
                                                   "Manhattan"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.clr.log$null.p.values.25.25.100indiv.MMiRKAT.Manhattan.clr.log)

list.pvalues.null.MMiRKAT.clr.log.100.100.500 = list("Euclidean"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.clr.log$null.p.values.100.100.500indiv.MMiRKAT.Euclidean.clr.log,
                                                   "Canberra"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.clr.log$null.p.values.100.100.500indiv.MMiRKAT.Canberra.clr.log,
                                                   "Manhattan"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.clr.log$null.p.values.100.100.500indiv.MMiRKAT.Manhattan.clr.log)

list.pvalues.null.MMiRKAT.ilr.log.25.25.100 = list("Euclidean"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.ilr.log$null.p.values.25.25.100indiv.MMiRKAT.Euclidean.ilr.log,
                                                   "Canberra"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.ilr.log$null.p.values.25.25.100indiv.MMiRKAT.Canberra.ilr.log,
                                                   "Manhattan"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.ilr.log$null.p.values.25.25.100indiv.MMiRKAT.Manhattan.ilr.log)

list.pvalues.null.MMiRKAT.alpha.log.25.25.100 = list("Euclidean"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\null.p.values.25.25.100indiv.MMiRKAT.Euclidean.alpha.log.RDS"),
                                                   "Canberra"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\null.p.values.25.25.100indiv.MMiRKAT.Canberra.alpha.log.RDS"),
                                                   "Manhattan"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\null.p.values.25.25.100indiv.MMiRKAT.Manhattan.alpha.log.RDS"))

list.pvalues.null.MMiRKAT.alpha.log.100.100.500 = list("Euclidean"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\null.p.values.100.100.500indiv.MMiRKAT.Euclidean.alpha.log.RDS"),
                                                     "Canberra"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\null.p.values.100.100.500indiv.MMiRKAT.Canberra.alpha.log.RDS"),
                                                     "Manhattan"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\null.p.values.100.100.500indiv.MMiRKAT.Manhattan.alpha.log.RDS"))

list.pvalues.power.MMiRKAT.clr.log.25.25.100 = list("Euclidean"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.clr.log$power.p.values.25.25.100indiv.MMiRKAT.Euclidean.clr.log,
                                                  "Canberra"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.clr.log$power.p.values.25.25.100indiv.MMiRKAT.Canberra.clr.log,
                                                  "Manhattan"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.clr.log$power.p.values.25.25.100indiv.MMiRKAT.Manhattan.clr.log)

list.pvalues.power.MMiRKAT.clr.log.100.100.500 = list("Euclidean"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.clr.log$power.p.values.100.100.500indiv.MMiRKAT.Euclidean.clr.log,
                                                    "Canberra"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.clr.log$power.p.values.100.100.500indiv.MMiRKAT.Canberra.clr.log,
                                                    "Manhattan"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.clr.log$power.p.values.100.100.500indiv.MMiRKAT.Manhattan.clr.log)


list.pvalues.power.MMiRKAT.ilr.log.25.25.100 = list("Euclidean"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.ilr.log$power.p.values.25.25.100indiv.MMiRKAT.Euclidean.ilr.log,
                                                    "Canberra"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.ilr.log$power.p.values.25.25.100indiv.MMiRKAT.Canberra.ilr.log,
                                                    "Manhattan"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.ilr.log$power.p.values.25.25.100indiv.MMiRKAT.Manhattan.ilr.log)

list.pvalues.power.MMiRKAT.alpha.log.100.100.500 = list("Euclidean"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\power.p.values.100.100.500indiv.MMiRKAT.Euclidean.alpha.log.RDS"),
                                                     "Canberra"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\power.p.values.100.100.500indiv.MMiRKAT.Canberra.alpha.log.RDS"),
                                                     "Manhattan"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\power.p.values.100.100.500indiv.MMiRKAT.Manhattan.alpha.log.RDS"))

list.pvalues.power.MMiRKAT.alpha.log.25.25.100 = list("Euclidean"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\power.p.values.25.25.100indiv.MMiRKAT.Euclidean.alpha.log.RDS"),
                                                      "Canberra"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\power.p.values.25.25.100indiv.MMiRKAT.Canberra.alpha.log.RDS"),
                                                      "Manhattan"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\power.p.values.25.25.100indiv.MMiRKAT.Manhattan.alpha.log.RDS"))

df.power.MMiRKAT.clr.log.25.25.100 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.clr.log.25.25.100$Canberra<=0.05)/1000,
                                                                      sum(list.pvalues.power.MMiRKAT.clr.log.25.25.100$Euclidean<=0.05)/1000,
                                                                      sum(list.pvalues.power.MMiRKAT.clr.log.25.25.100$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.MMiRKAT.clr.log.100.100.500 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.clr.log.100.100.500$Canberra<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.clr.log.100.100.500$Euclidean<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.clr.log.100.100.500$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))


df.power.MMiRKAT.ilr.log.25.25.100 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.ilr.log.25.25.100$Canberra<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.ilr.log.25.25.100$Euclidean<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.ilr.log.25.25.100$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))


df.power.MMiRKAT.alpha.log.25.25.100 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.alpha.log.25.25.100$Canberra<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.alpha.log.25.25.100$Euclidean<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.alpha.log.25.25.100$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.MMiRKAT.alpha.log.100.100.500 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.alpha.log.100.100.500$Canberra<=0.05)/1000,
                                                              sum(list.pvalues.power.MMiRKAT.alpha.log.100.100.500$Euclidean<=0.05)/1000,
                                                              sum(list.pvalues.power.MMiRKAT.alpha.log.100.100.500$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

ggarrange(gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.clr.log.25.25.100)+labs(colour="Distance Kernel")+ggtitle("CLR"),
          gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.ilr.log.25.25.100)+labs(colour="Distance Kernel")+ggtitle("ILR"),
          gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.alpha.log.25.25.100)+labs(colour="Distance Kernel")+ggtitle("Alpha"), 
          ggplot(df.power.MMiRKAT.clr.log.25.25.100, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.MMiRKAT.ilr.log.25.25.100, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.MMiRKAT.alpha.log.25.25.100, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          labels = c("A","B","C","D","E","F"), ncol=3, nrow=2)

ggarrange(gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.clr.log.100.100.500)+labs(colour="Distance Kernel")+ggtitle("CLR"),
          gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.alpha.log.100.100.500)+labs(colour="Distance Kernel")+ggtitle("Alpha"),
          ggplot(df.power.MMiRKAT.clr.log.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.MMiRKAT.alpha.log.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          labels = c("A","B","C","D"), ncol=2, nrow=2)



list.pvalues.null.MMiRKAT.clr.original.25.25.100 = list("Euclidean"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.clr.original$null.p.values.25.25.100indiv.MMiRKAT.Euclidean.compositionnal,
                                                   "Canberra"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.clr.original$null.p.values.25.25.100indiv.MMiRKAT.Canberra.compositionnal,
                                                   "Manhattan"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.clr.original$null.p.values.25.25.100indiv.MMiRKAT.Manhattan.compositionnal)

list.pvalues.null.MMiRKAT.clr.original.100.100.500 = list("Euclidean"=pvalues.null.MMiRKAT.100.100.500$`list.results.null.MMiRKAT.100.100.500.clr.original `$null.p.values.100.100.500indiv.MMiRKAT.Euclidean.clr.original,
                                                     "Canberra"=pvalues.null.MMiRKAT.100.100.500$`list.results.null.MMiRKAT.100.100.500.clr.original `$null.p.values.100.100.500indiv.MMiRKAT.Canberra.clr.original,
                                                     "Manhattan"=pvalues.null.MMiRKAT.100.100.500$`list.results.null.MMiRKAT.100.100.500.clr.original `$null.p.values.100.100.500indiv.MMiRKAT.Manhattan.clr.original)

list.pvalues.null.MMiRKAT.ilr.original.25.25.100 = list("Euclidean"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.ilr.original$null.p.values.25.25.100indiv.MMiRKAT.Euclidean.ilr.original,
                                                   "Canberra"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.ilr.original$null.p.values.25.25.100indiv.MMiRKAT.Canberra.ilr.original,
                                                   "Manhattan"=pvalues.null.MMiRKAT.25.25.100$list.results.null.MMiRKAT.25.25.100.ilr.original$null.p.values.25.25.100indiv.MMiRKAT.Manhattan.ilr.original)

list.pvalues.null.MMiRKAT.ilr.original.100.100.500 = list("Euclidean"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.ilr.original$null.p.values.100.100.500indiv.MMiRKAT.Euclidean.ilr.original,
                                                        "Canberra"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.ilr.original$null.p.values.100.100.500indiv.MMiRKAT.Canberra.ilr.original,
                                                        "Manhattan"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.ilr.original$null.p.values.100.100.500indiv.MMiRKAT.Manhattan.ilr.original)

list.pvalues.null.MMiRKAT.alpha.original.25.25.100 = list("Euclidean"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\null.p.values.25.25.100indiv.MMiRKAT.Euclidean.alpha.original.RDS"),
                                                     "Canberra"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\null.p.values.25.25.100indiv.MMiRKAT.Canberra.alpha.original.RDS"),
                                                     "Manhattan"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\null.p.values.25.25.100indiv.MMiRKAT.Manhattan.alpha.original.RDS"))

list.pvalues.null.MMiRKAT.alpha.original.100.100.500 = list("Euclidean"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\null.p.values.100.100.500indiv.MMiRKAT.Euclidean.alpha.original.RDS"),
                                                       "Canberra"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\null.p.values.100.100.500indiv.MMiRKAT.Canberra.alpha.original.RDS"),
                                                       "Manhattan"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\null.p.values.100.100.500indiv.MMiRKAT.Manhattan.alpha.original.RDS"))

list.pvalues.power.MMiRKAT.clr.original.25.25.100 = list("Euclidean"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.clr.original$power.p.values.25.25.100indiv.MMiRKAT.Euclidean.compositionnal,
                                                    "Canberra"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.clr.original$power.p.values.25.25.100indiv.MMiRKAT.Canberra.compositionnal,
                                                    "Manhattan"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.clr.original$power.p.values.25.25.100indiv.MMiRKAT.Manhattan.compositionnal)

list.pvalues.power.MMiRKAT.clr.original.100.100.500 = list("Euclidean"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.clr.original$power.p.values.100.100.500indiv.MMiRKAT.Euclidean.clr.original,
                                                      "Canberra"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.clr.original$power.p.values.100.100.500indiv.MMiRKAT.Canberra.clr.original,
                                                      "Manhattan"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.clr.original$power.p.values.100.100.500indiv.MMiRKAT.Manhattan.clr.original)


list.pvalues.power.MMiRKAT.ilr.original.25.25.100 = list("Euclidean"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.ilr.original$power.p.values.25.25.100indiv.MMiRKAT.Euclidean.ilr.original,
                                                    "Canberra"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.ilr.original$power.p.values.25.25.100indiv.MMiRKAT.Canberra.ilr.original,
                                                    "Manhattan"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.ilr.original$power.p.values.25.25.100indiv.MMiRKAT.Manhattan.ilr.original)

list.pvalues.power.MMiRKAT.ilr.original.100.100.500 = list("Euclidean"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.ilr.original$power.p.values.100.100.500indiv.MMiRKAT.Euclidean.ilr.original,
                                                         "Canberra"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.ilr.original$power.p.values.100.100.500indiv.MMiRKAT.Canberra.ilr.original,
                                                         "Manhattan"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.ilr.original$power.p.values.100.100.500indiv.MMiRKAT.Manhattan.ilr.original)

list.pvalues.power.MMiRKAT.alpha.original.100.100.500 = list("Euclidean"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\power.p.values.100.100.500indiv.MMiRKAT.Euclidean.alpha.original.RDS"),
                                                        "Canberra"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\power.p.values.100.100.500indiv.MMiRKAT.Canberra.alpha.original.RDS"),
                                                        "Manhattan"=readRDS("data\\MMiRKAT\\alpha\\100_100_500\\power.p.values.100.100.500indiv.MMiRKAT.Manhattan.alpha.original.RDS"))

list.pvalues.power.MMiRKAT.alpha.original.25.25.100 = list("Euclidean"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\power.p.values.25.25.100indiv.MMiRKAT.Euclidean.alpha.original.RDS"),
                                                      "Canberra"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\power.p.values.25.25.100indiv.MMiRKAT.Canberra.alpha.original.RDS"),
                                                      "Manhattan"=readRDS("data\\MMiRKAT\\alpha\\25_25_100\\power.p.values.25.25.100indiv.MMiRKAT.Manhattan.alpha.original.RDS"))

df.power.MMiRKAT.clr.original.25.25.100 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.clr.original.25.25.100$Canberra<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.clr.original.25.25.100$Euclidean<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.clr.original.25.25.100$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.MMiRKAT.clr.original.100.100.500 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.clr.original.100.100.500$Canberra<=0.05)/1000,
                                                              sum(list.pvalues.power.MMiRKAT.clr.original.100.100.500$Euclidean<=0.05)/1000,
                                                              sum(list.pvalues.power.MMiRKAT.clr.original.100.100.500$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))


df.power.MMiRKAT.ilr.original.25.25.100 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.ilr.original.25.25.100$Canberra<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.ilr.original.25.25.100$Euclidean<=0.05)/1000,
                                                            sum(list.pvalues.power.MMiRKAT.ilr.original.25.25.100$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.MMiRKAT.ilr.original.100.100.500 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.ilr.original.100.100.500$Canberra<=0.05)/1000,
                                                                 sum(list.pvalues.power.MMiRKAT.ilr.original.100.100.500$Euclidean<=0.05)/1000,
                                                                 sum(list.pvalues.power.MMiRKAT.ilr.original.100.100.500$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.MMiRKAT.alpha.original.25.25.100 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.alpha.original.25.25.100$Canberra<=0.05)/1000,
                                                              sum(list.pvalues.power.MMiRKAT.alpha.original.25.25.100$Euclidean<=0.05)/1000,
                                                              sum(list.pvalues.power.MMiRKAT.alpha.original.25.25.100$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.MMiRKAT.alpha.original.100.100.500 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.alpha.original.100.100.500$Canberra<=0.05)/1000,
                                                                sum(list.pvalues.power.MMiRKAT.alpha.original.100.100.500$Euclidean<=0.05)/1000,
                                                                sum(list.pvalues.power.MMiRKAT.alpha.original.100.100.500$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

ggarrange(gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.clr.original.25.25.100)+labs(colour="Distance Kernel")+ggtitle("CLR"),
          gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.ilr.original.25.25.100)+labs(colour="Distance Kernel")+ggtitle("ILR"),
          gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.alpha.original.25.25.100)+labs(colour="Distance Kernel")+ggtitle("Alpha"), 
          ggplot(df.power.MMiRKAT.clr.original.25.25.100, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.MMiRKAT.ilr.original.25.25.100, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.MMiRKAT.alpha.original.25.25.100, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          labels = c("A","B","C","D","E","F"), ncol=3, nrow=2)

ggarrange(gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.clr.original.100.100.500)+labs(colour="Distance Kernel")+ggtitle("CLR"),
          gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.ilr.original.100.100.500)+labs(colour="Distance Kernel")+ggtitle("ILR"),
          gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.alpha.original.100.100.500)+labs(colour="Distance Kernel")+ggtitle("Alpha"), 
          ggplot(df.power.MMiRKAT.clr.original.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.MMiRKAT.ilr.original.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.MMiRKAT.alpha.original.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          labels = c("A","B","C","D","E","F"), ncol=3, nrow=2)

list.pvalues.null.MMiRKAT.original.original.100.100.500 = list("Euclidean"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.original.original$null.p.values.100.100.500indiv.MMiRKAT.Euclidean.original.original,
                                                          "Canberra"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.original.original$null.p.values.100.100.500indiv.MMiRKAT.Canberra.original.original,
                                                          "Manhattan"=pvalues.null.MMiRKAT.100.100.500$list.results.null.MMiRKAT.100.100.500.original.original$null.p.values.100.100.500indiv.MMiRKAT.Manhattan.original.original)

list.pvalues.null.MMiRKAT.original.original.25.25.100 = list("Euclidean"=pvalues.null.MMiRKAT.25.25.100$`list.results.null.MMiRKAT.25.25.100.original.original `$null.p.values.25.25.100indiv.MMiRKAT.Euclidean.original.original,
                                                        "Canberra"=pvalues.null.MMiRKAT.25.25.100$`list.results.null.MMiRKAT.25.25.100.original.original `$null.p.values.25.25.100indiv.MMiRKAT.Canberra.original.original,
                                                        "Manhattan"=pvalues.null.MMiRKAT.25.25.100$`list.results.null.MMiRKAT.25.25.100.original.original `$null.p.values.25.25.100indiv.MMiRKAT.Manhattan.original.original)


list.pvalues.power.MMiRKAT.original.original.25.25.100 = list("Euclidean"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.original.original$power.p.values.25.25.100indiv.MMiRKAT.Euclidean.original.original,
                                                         "Canberra"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.original.original$power.p.values.25.25.100indiv.MMiRKAT.Canberra.original.original,
                                                         "Manhattan"=pvalues.alter.MMiRKAT.25.25.100$list.results.power.MMiRKAT.25.25.100.original.original$power.p.values.25.25.100indiv.MMiRKAT.Canberra.original.original)

list.pvalues.power.MMiRKAT.original.original.100.100.500 = list("Euclidean"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.original.original$power.p.values.100.100.500indiv.MMiRKAT.Euclidean.original.original,
                                                              "Canberra"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.original.original$power.p.values.100.100.500indiv.MMiRKAT.Canberra.original.original,
                                                              "Manhattan"=pvalues.alter.MMiRKAT.100.100.500$list.results.power.MMiRKAT.100.100.500.original.original$power.p.values.100.100.500indiv.MMiRKAT.Canberra.original.original)

df.power.MMiRKAT.original.original.25.25.100 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.original.original.25.25.100$Canberra<=0.05)/1000,
                                                                 sum(list.pvalues.power.MMiRKAT.original.original.25.25.100$Euclidean<=0.05)/1000,
                                                                 sum(list.pvalues.power.MMiRKAT.original.original.25.25.100$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))

df.power.MMiRKAT.original.original.100.100.500 = data.frame("Power" = c(sum(list.pvalues.power.MMiRKAT.original.original.100.100.500$Canberra<=0.05)/1000,
                                                                   sum(list.pvalues.power.MMiRKAT.original.original.100.100.500$Euclidean<=0.05)/1000,
                                                                   sum(list.pvalues.power.MMiRKAT.original.original.100.100.500$Manhattan<=0.05)/1000), "Distance.Kernel" = c("Canberra", "Euclidean", "Manhattan"))


ggarrange(gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.original.original.25.25.100)+labs(colour="Distance Kernel")+ggtitle("Low Dimensional"),
          gg_qqplot_facet_grid(list.pvalues.null.MMiRKAT.original.original.100.100.500)+labs(colour="Distance Kernel")+ggtitle("High Dimensional"),
          ggplot(df.power.MMiRKAT.original.original.25.25.100, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          ggplot(df.power.MMiRKAT.original.original.100.100.500, aes(x=Distance.Kernel,y=Power,fill=Distance.Kernel))+geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+theme(legend.position = "none")+ylim(c(0,1))+xlab("Distance Kernel"),
          labels = c("A","B","C","D"), ncol=2, nrow=2)
#Data Summarization

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
                                                                             "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))

ggplot(df.results.clr.log.multivariate.factorization.based.100.100.500, aes(x=Method, y=Explained.Variance, fill=Method))+geom_boxplot()+theme_bw()+ylab("% Explained Variance")+theme(legend.position = "none")+ylim(c(0,100))



CCA.alpha.log.25.25.100 = readRDS("data\\results_CCA\\alpha\\25_25_100\\redundancy_cca_alpha_microbiome_log_metabolome_25_25_100indiv.RDS")
PLS.Can.alpha.log.25.25.100 = readRDS("data\\results_PLS\\alpha\\25_25_100\\canonical\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_canonical.RDS")
PLS.Reg.alpha.log.25.25.100 = readRDS("data\\results_PLS\\alpha\\25_25_100\\regression\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_regression.RDS")
RDA.alpha.log.25.25.100 = readRDS("data\\results_RDA\\alpha\\25_25_100\\explained_variance_alpha_microbiome_log_metabolome_rda_25_25_100indiv.RDS")
MOFA.alpha.log.25.25.100 = readRDS("data\\results_MOFA2\\alpha\\25_25_100\\MOFA2_alpha_microbiome_log_metabolome_Scenario_25_25_100indiv.RDS")

CCA.ilr.log.25.25.100 = readRDS("data\\results_CCA\\ilr\\25_25_100\\redundancy_cca_ilr_microbiome_log_metabolome_25_25_100indiv.RDS")
PLS.Can.ilr.log.25.25.100 = readRDS("data\\results_PLS\\ilr\\25_25_100\\canonical\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_canonical.RDS")
PLS.Reg.ilr.log.25.25.100 = readRDS("data\\results_PLS\\ilr\\25_25_100\\regression\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_regression.RDS")
RDA.ilr.log.25.25.100 = readRDS("data\\results_RDA\\ilr\\25_25_100\\explained_variance_ilr_microbiome_log_metabolome_rda_25_25_100indiv.RDS")
MOFA.ilr.log.25.25.100 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv.RDS")


CCA.clr.log.25.25.100 = readRDS("data\\results_CCA\\clr\\25_25_100\\redundancy_cca_clr_microbiome_log_metabolome_25_25_100indiv.RDS")
PLS.Can.clr.log.25.25.100 = readRDS("data\\results_PLS\\clr\\25_25_100\\canonical\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_canonical.RDS")
PLS.Reg.clr.log.25.25.100 = readRDS("data\\results_PLS\\clr\\25_25_100\\regression\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_regression.RDS")
RDA.clr.log.25.25.100 = readRDS("data\\results_RDA\\clr\\25_25_100\\explained_variance_clr_microbiome_log_metabolome_rda_25_25_100indiv.RDS")
MOFA.clr.log.25.25.100 = readRDS("data\\results_MOFA2\\clr\\25_25_100\\MOFA2_clr_microbiome_log_metabolome_Scenario_25_25_100indiv.RDS")


df.results.alpha.log.multivariate.factorization.based.25.25.100 = data.frame("Explained-Variance"=c(CCA.alpha.log.25.25.100*100,
                                                                                                      PLS.Reg.alpha.log.25.25.100*100,
                                                                                                      PLS.Can.alpha.log.25.25.100*100,
                                                                                                      RDA.alpha.log.25.25.100*100,
                                                                                                      sapply(1:1000, function(X) MOFA.alpha.log.25.25.100[[X]]$r2_total$group1[1])),
                                                                               "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))


df.results.ilr.log.multivariate.factorization.based.25.25.100 = data.frame("Explained-Variance"=c(CCA.ilr.log.25.25.100*100,
                                                                                                      PLS.Reg.ilr.log.25.25.100*100,
                                                                                                      PLS.Can.ilr.log.25.25.100*100,
                                                                                                      RDA.ilr.log.25.25.100*100,
                                                                                                      sapply(1:1000, function(X) MOFA.ilr.log.25.25.100[[X]]$r2_total$group1[1])),
                                                                               "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))


df.results.clr.log.multivariate.factorization.based.25.25.100 = data.frame("Explained-Variance"=c(as.numeric(CCA.clr.log.25.25.100)*100,
                                                                                                    PLS.Reg.clr.log.25.25.100*100,
                                                                                                    PLS.Can.clr.log.25.25.100*100,
                                                                                                    RDA.clr.log.25.25.100*100,
                                                                                                    sapply(1:1000, function(X) MOFA.clr.log.25.25.100[[X]]$r2_total$group1[1])),
                                                                             "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))



df.results.alpha.log.multivariate.factorization.based.25.25.100$Normalization = "Alpha"
df.results.ilr.log.multivariate.factorization.based.25.25.100$Normalization = "ILR"
df.results.clr.log.multivariate.factorization.based.25.25.100$Normalization = "CLR"

df.results.multivariate.factorization.based.25.25.100 = rbind(df.results.alpha.log.multivariate.factorization.based.25.25.100,df.results.ilr.log.multivariate.factorization.based.25.25.100,
                                                              df.results.clr.log.multivariate.factorization.based.25.25.100)


ggplot(df.results.multivariate.factorization.based.25.25.100, aes(x=Method, y=Explained.Variance, fill=Method))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")


CCA.alpha.25.25.100 = readRDS("data\\results_CCA\\alpha\\25_25_100\\redundancy_cca_alpha_microbiome_metabolome_25_25_100indiv.RDS")
PLS.Can.alpha.25.25.100 = readRDS("data\\results_PLS\\alpha\\25_25_100\\canonical\\redundancy_pls_alpha_microbiome_metabolome_25_25_100indiv_canonical.RDS")
PLS.Reg.alpha.25.25.100 = readRDS("data\\results_PLS\\alpha\\25_25_100\\regression\\redundancy_pls_alpha_microbiome_metabolome_25_25_100indiv_regression.RDS")
RDA.alpha.25.25.100 = readRDS("data\\results_RDA\\alpha\\25_25_100\\explained_variance_alpha_microbiome_metabolome_rda_25_25_100indiv.RDS")
MOFA.alpha.25.25.100 = readRDS("data\\results_MOFA2\\alpha\\25_25_100\\MOFA2_alpha_microbiome_metabolome_Scenario_25_25_100indiv.RDS")

CCA.ilr.25.25.100 = readRDS("data\\results_CCA\\ilr\\25_25_100\\redundancy_cca_ilr_microbiome_metabolome_25_25_100indiv.RDS")
PLS.Can.ilr.25.25.100 = readRDS("data\\results_PLS\\ilr\\25_25_100\\canonical\\redundancy_pls_ilr_microbiome_metabolome_25_25_100indiv_canonical.RDS")
PLS.Reg.ilr.25.25.100 = readRDS("data\\results_PLS\\ilr\\25_25_100\\regression\\redundancy_pls_ilr_microbiome_metabolome_25_25_100indiv_regression.RDS")
RDA.ilr.25.25.100 = readRDS("data\\results_RDA\\ilr\\25_25_100\\explained_variance_ilr_microbiome_metabolome_rda_25_25_100indiv.RDS")
MOFA.ilr.25.25.100 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_metabolome_Scenario_25_25_100indiv.RDS")


CCA.clr.25.25.100 = readRDS("data\\results_CCA\\clr\\25_25_100\\redundancy_cca_clr_microbiome_metabolome_25_25_100indiv.RDS")
PLS.Can.clr.25.25.100 = readRDS("data\\results_PLS\\clr\\25_25_100\\canonical\\redundancy_pls_clr_microbiome_metabolome_25_25_100indiv_canonical.RDS")
PLS.Reg.clr.25.25.100 = readRDS("data\\results_PLS\\clr\\25_25_100\\regression\\redundancy_pls_clr_microbiome_metabolome_25_25_100indiv_regression.RDS")
RDA.clr.25.25.100 = readRDS("data\\results_RDA\\clr\\25_25_100\\explained_variance_clr_microbiome_metabolome_rda_25_25_100indiv.RDS")
MOFA.clr.25.25.100 = readRDS("data\\results_MOFA2\\clr\\25_25_100\\MOFA2_clr_microbiome_metabolome_Scenario_25_25_100indiv.RDS")


df.results.alpha.multivariate.factorization.based.25.25.100 = data.frame("Explained-Variance"=c(CCA.alpha.25.25.100*100,
                                                                                                    PLS.Reg.alpha.25.25.100*100,
                                                                                                    PLS.Can.alpha.25.25.100*100,
                                                                                                    RDA.alpha.25.25.100*100,
                                                                                                    sapply(1:1000, function(X) MOFA.alpha.25.25.100[[X]]$r2_total$group1[1])),
                                                                             "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))


df.results.ilr.multivariate.factorization.based.25.25.100 = data.frame("Explained-Variance"=c(CCA.ilr.25.25.100*100,
                                                                                                  PLS.Reg.ilr.25.25.100*100,
                                                                                                  PLS.Can.ilr.25.25.100*100,
                                                                                                  RDA.ilr.25.25.100*100,
                                                                                                  sapply(1:1000, function(X) MOFA.ilr.25.25.100[[X]]$r2_total$group1[1])),
                                                                           "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))


df.results.clr.multivariate.factorization.based.25.25.100 = data.frame("Explained-Variance"=c(as.numeric(CCA.clr.25.25.100)*100,
                                                                                                  PLS.Reg.clr.25.25.100*100,
                                                                                                  PLS.Can.clr.25.25.100*100,
                                                                                                  RDA.clr.25.25.100*100,
                                                                                                  sapply(1:1000, function(X) MOFA.clr.25.25.100[[X]]$r2_total$group1[1])),
                                                                           "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))



df.results.alpha.multivariate.factorization.based.25.25.100$Normalization = "Alpha"
df.results.ilr.multivariate.factorization.based.25.25.100$Normalization = "ILR"
df.results.clr.multivariate.factorization.based.25.25.100$Normalization = "CLR"

df.results.multivariate.factorization.based.original.metabolome.25.25.100 = rbind(df.results.alpha.multivariate.factorization.based.25.25.100,df.results.ilr.multivariate.factorization.based.25.25.100,
                                                              df.results.clr.multivariate.factorization.based.25.25.100)


ggplot(df.results.multivariate.factorization.based.original.metabolome.25.25.100, aes(x=Method, y=Explained.Variance, fill=Method))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")

CCA.alpha.100.100.500 = readRDS("data\\results_CCA\\alpha\\100_100_500\\redundancy_cca_alpha_microbiome_metabolome_100_100_500indiv.RDS")
PLS.Can.alpha.100.100.500 = readRDS("data\\results_PLS\\alpha\\100_100_500\\canonical\\redundancy_pls_alpha_microbiome_metabolome_100_100_500indiv_canonical.RDS")
PLS.Reg.alpha.100.100.500 = readRDS("data\\results_PLS\\alpha\\100_100_500\\regression\\redundancy_pls_alpha_microbiome_metabolome_100_100_500indiv_regression.RDS")
RDA.alpha.100.100.500 = readRDS("data\\results_RDA\\alpha\\100_100_500\\explained_variance_alpha_microbiome_metabolome_rda_100_100_500indiv.RDS")
MOFA.alpha.100.100.500 = readRDS("data\\results_MOFA2\\alpha\\100_100_500\\MOFA2_alpha_microbiome_metabolome_Scenario_100_100_500indiv.RDS")

CCA.ilr.100.100.500 = readRDS("data\\results_CCA\\ilr\\100_100_500\\redundancy_cca_ilr_microbiome_metabolome_100_100_500indiv.RDS")
PLS.Can.ilr.100.100.500 = readRDS("data\\results_PLS\\ilr\\100_100_500\\canonical\\redundancy_pls_ilr_microbiome_metabolome_100_100_500indiv_canonical.RDS")
PLS.Reg.ilr.100.100.500 = readRDS("data\\results_PLS\\ilr\\100_100_500\\regression\\redundancy_pls_ilr_microbiome_metabolome_100_100_500indiv_regression.RDS")
RDA.ilr.100.100.500 = readRDS("data\\results_RDA\\ilr\\100_100_500\\explained_variance_ilr_microbiome_metabolome_rda_100_100_500indiv.RDS")
MOFA.ilr.100.100.500 = readRDS("data\\results_MOFA2\\ilr\\100_100_500\\MOFA2_ilr_microbiome_metabolome_Scenario_100_100_500indiv.RDS")


CCA.clr.100.100.500 = readRDS("data\\results_CCA\\clr\\100_100_500\\redundancy_cca_clr_microbiome_metabolome_100_100_500indiv.RDS")
PLS.Can.clr.100.100.500 = readRDS("data\\results_PLS\\clr\\100_100_500\\canonical\\redundancy_pls_clr_microbiome_metabolome_100_100_500indiv_canonical.RDS")
PLS.Reg.clr.100.100.500 = readRDS("data\\results_PLS\\clr\\100_100_500\\regression\\redundancy_pls_clr_microbiome_metabolome_100_100_500indiv_regression.RDS")
RDA.clr.100.100.500 = readRDS("data\\results_RDA\\clr\\100_100_500\\explained_variance_clr_microbiome_metabolome_rda_100_100_500indiv.RDS")
MOFA.clr.100.100.500 = readRDS("data\\results_MOFA2\\clr\\100_100_500\\MOFA2_clr_microbiome_metabolome_Scenario_100_100_500indiv.RDS")


df.results.alpha.multivariate.factorization.based.100.100.500 = data.frame("Explained-Variance"=c(CCA.alpha.100.100.500*100,
                                                                                                PLS.Reg.alpha.100.100.500*100,
                                                                                                PLS.Can.alpha.100.100.500*100,
                                                                                                RDA.alpha.100.100.500*100,
                                                                                                sapply(1:1000, function(X) MOFA.alpha.100.100.500[[X]]$r2_total$group1[1])),
                                                                         "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))


df.results.ilr.multivariate.factorization.based.100.100.500 = data.frame("Explained-Variance"=c(CCA.ilr.100.100.500*100,
                                                                                              PLS.Reg.ilr.100.100.500*100,
                                                                                              PLS.Can.ilr.100.100.500*100,
                                                                                              RDA.ilr.100.100.500*100,
                                                                                              sapply(1:1000, function(X) MOFA.ilr.100.100.500[[X]]$r2_total$group1[1])),
                                                                       "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))


df.results.clr.multivariate.factorization.based.100.100.500 = data.frame("Explained-Variance"=c(as.numeric(CCA.clr.100.100.500)*100,
                                                                                              PLS.Reg.clr.100.100.500*100,
                                                                                              PLS.Can.clr.100.100.500*100,
                                                                                              RDA.clr.100.100.500*100,
                                                                                              sapply(1:1000, function(X) MOFA.clr.100.100.500[[X]]$r2_total$group1[1])),
                                                                       "Method" = rep(c("CCA","PLS-Reg","PLS-Can","RDA","MOFA2"), each=1000))



df.results.alpha.multivariate.factorization.based.100.100.500$Normalization = "Alpha"
df.results.ilr.multivariate.factorization.based.100.100.500$Normalization = "ILR"
df.results.clr.multivariate.factorization.based.100.100.500$Normalization = "CLR"

df.results.multivariate.factorization.based.original.metabolome.100.100.500 = rbind(df.results.alpha.multivariate.factorization.based.100.100.500,df.results.ilr.multivariate.factorization.based.100.100.500,
                                                                                  df.results.clr.multivariate.factorization.based.100.100.500)


ggplot(df.results.multivariate.factorization.based.original.metabolome.100.100.500, aes(x=Method, y=Explained.Variance, fill=Method))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")

CCA.clr.log.25.25.100.2.2 = readRDS("data\\results_CCA\\clr\\25_25_100\\redundancy_cca_clr_microbiome_log_metabolome_25_25_100indiv_2_2.RDS")
CCA.clr.log.25.25.100.4.4 = readRDS("data\\results_CCA\\clr\\25_25_100\\redundancy_cca_clr_microbiome_log_metabolome_25_25_100indiv_4_4.RDS")
CCA.clr.log.25.25.100.6.6 = readRDS("data\\results_CCA\\clr\\25_25_100\\redundancy_cca_clr_microbiome_log_metabolome_25_25_100indiv_6_6.RDS")

CCA.ilr.log.25.25.100.2.2 = readRDS("data\\results_CCA\\ilr\\25_25_100\\redundancy_cca_ilr_microbiome_log_metabolome_25_25_100indiv_2_2.RDS")
CCA.ilr.log.25.25.100.4.4 = readRDS("data\\results_CCA\\ilr\\25_25_100\\redundancy_cca_ilr_microbiome_log_metabolome_25_25_100indiv_4_4.RDS")
CCA.ilr.log.25.25.100.6.6 = readRDS("data\\results_CCA\\ilr\\25_25_100\\redundancy_cca_ilr_microbiome_log_metabolome_25_25_100indiv_6_6.RDS")

CCA.alpha.log.25.25.100.2.2 = readRDS("data\\results_CCA\\alpha\\25_25_100\\redundancy_cca_alpha_microbiome_log_metabolome_25_25_100indiv_2_2.RDS")
CCA.alpha.log.25.25.100.4.4 = readRDS("data\\results_CCA\\alpha\\25_25_100\\redundancy_cca_alpha_microbiome_log_metabolome_25_25_100indiv_4_4.RDS")
CCA.alpha.log.25.25.100.6.6 = readRDS("data\\results_CCA\\alpha\\25_25_100\\redundancy_cca_alpha_microbiome_log_metabolome_25_25_100indiv_6_6.RDS")


PLS.Reg.clr.log.25.25.100.2.2 = readRDS("data\\results_PLS\\clr\\25_25_100\\regression\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_regression_2_2_effect_0_2.RDS")
PLS.Reg.clr.log.25.25.100.4.4 = readRDS("data\\results_PLS\\clr\\25_25_100\\regression\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effect_0_2.RDS")
PLS.Reg.clr.log.25.25.100.6.6 = readRDS("data\\results_PLS\\clr\\25_25_100\\regression\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_regression_6_6_effect_0_2.RDS")

PLS.Reg.ilr.log.25.25.100.2.2 = readRDS("data\\results_PLS\\ilr\\25_25_100\\regression\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_regression_2_2_effect0_2.RDS")
PLS.Reg.ilr.log.25.25.100.4.4 = readRDS("data\\results_PLS\\ilr\\25_25_100\\regression\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effect0_2.RDS")
PLS.Reg.ilr.log.25.25.100.6.6 = readRDS("data\\results_PLS\\ilr\\25_25_100\\regression\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_regression_6_6_effect0_2.RDS")

PLS.Reg.alpha.log.25.25.100.2.2 = readRDS("data\\results_PLS\\alpha\\25_25_100\\regression\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_regression_2_2_effect_0_2.RDS")
PLS.Reg.alpha.log.25.25.100.4.4 = readRDS("data\\results_PLS\\alpha\\25_25_100\\regression\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effect_0_2.RDS")
PLS.Reg.alpha.log.25.25.100.6.6 = readRDS("data\\results_PLS\\alpha\\25_25_100\\regression\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_regression_6_6_effect_0_2.RDS")


PLS.Canonical.clr.log.25.25.100.2.2 = readRDS("data\\results_PLS\\clr\\25_25_100\\canonical\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_canonical_2_2_effect_0_2.RDS")
PLS.Canonical.clr.log.25.25.100.4.4 = readRDS("data\\results_PLS\\clr\\25_25_100\\canonical\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effect_0_2.RDS")
PLS.Canonical.clr.log.25.25.100.6.6 = readRDS("data\\results_PLS\\clr\\25_25_100\\canonical\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_canonical_6_6_effect_0_2.RDS")

PLS.Canonical.ilr.log.25.25.100.2.2 = readRDS("data\\results_PLS\\ilr\\25_25_100\\canonical\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_canonical_2_2_effect0_2.RDS")
PLS.Canonical.ilr.log.25.25.100.4.4 = readRDS("data\\results_PLS\\ilr\\25_25_100\\canonical\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effect0_2.RDS")
PLS.Canonical.ilr.log.25.25.100.6.6 = readRDS("data\\results_PLS\\ilr\\25_25_100\\canonical\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_canonical_6_6_effect0_2.RDS")

PLS.Canonical.alpha.log.25.25.100.2.2 = readRDS("data\\results_PLS\\alpha\\25_25_100\\canonical\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_canonical_2_2_effect_0_2.RDS")
PLS.Canonical.alpha.log.25.25.100.4.4 = readRDS("data\\results_PLS\\alpha\\25_25_100\\canonical\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effect_0_2.RDS")
PLS.Canonical.alpha.log.25.25.100.6.6 = readRDS("data\\results_PLS\\alpha\\25_25_100\\canonical\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_canonical_6_6_effect_0_2.RDS")


RDA.clr.log.25.25.100.2.2 = readRDS("data\\results_RDA\\clr\\25_25_100\\explained_variance_clr_microbiome_log_metabolome_rda_25_25_100indiv_2_2_effect0_2.RDS")
RDA.clr.log.25.25.100.4.4 = readRDS("data\\results_RDA\\clr\\25_25_100\\explained_variance_clr_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effect0_2.RDS")
RDA.clr.log.25.25.100.6.6 = readRDS("data\\results_RDA\\clr\\25_25_100\\explained_variance_clr_microbiome_log_metabolome_rda_25_25_100indiv_6_6_effect0_2.RDS")

RDA.ilr.log.25.25.100.2.2 = readRDS("data\\results_RDA\\ilr\\25_25_100\\explained_variance_ilr_microbiome_log_metabolome_rda_25_25_100indiv_2_2_effect0_2.RDS")
RDA.ilr.log.25.25.100.4.4 = readRDS("data\\results_RDA\\ilr\\25_25_100\\explained_variance_ilr_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effect0_2.RDS")
RDA.ilr.log.25.25.100.6.6 = readRDS("data\\results_RDA\\ilr\\25_25_100\\explained_variance_ilr_microbiome_log_metabolome_rda_25_25_100indiv_6_6_effect0_2.RDS")

RDA.alpha.log.25.25.100.2.2 = readRDS("data\\results_RDA\\alpha\\25_25_100\\explained_variance_alpha_microbiome_log_metabolome_rda_25_25_100indiv_2_2_effect0_2.RDS")
RDA.alpha.log.25.25.100.4.4 = readRDS("data\\results_RDA\\alpha\\25_25_100\\explained_variance_alpha_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effect0_2.RDS")
RDA.alpha.log.25.25.100.6.6 = readRDS("data\\results_RDA\\alpha\\25_25_100\\explained_variance_alpha_microbiome_log_metabolome_rda_25_25_100indiv_6_6_effect0_2.RDS")


MOFA.clr.log.25.25.100.2.2 = readRDS("data\\results_MOFA2\\clr\\25_25_100\\MOFA2_clr_microbiome_log_metabolome_Scenario_25_25_100indiv_2_2.RDS")
MOFA.clr.log.25.25.100.4.4 = readRDS("data\\results_MOFA2\\clr\\25_25_100\\MOFA2_clr_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4.RDS")
MOFA.clr.log.25.25.100.6.6 = readRDS("data\\results_MOFA2\\clr\\25_25_100\\MOFA2_clr_microbiome_log_metabolome_Scenario_25_25_100indiv_6_6.RDS")

MOFA.ilr.log.25.25.100.2.2 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv_2_2.RDS")
MOFA.ilr.log.25.25.100.4.4 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4.RDS")
MOFA.ilr.log.25.25.100.6.6 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv_6_6.RDS")

MOFA.alpha.log.25.25.100.2.2 = readRDS("data\\results_MOFA2\\alpha\\25_25_100\\MOFA2_alpha_microbiome_log_metabolome_Scenario_25_25_100indiv_2_2.RDS")
MOFA.alpha.log.25.25.100.4.4 = readRDS("data\\results_MOFA2\\alpha\\25_25_100\\MOFA2_alpha_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4.RDS")
MOFA.alpha.log.25.25.100.6.6 = readRDS("data\\results_MOFA2\\alpha\\25_25_100\\MOFA2_alpha_microbiome_log_metabolome_Scenario_25_25_100indiv_6_6.RDS")

df.results.cca.multivariate.factorization.based.25.25.100.various.associations = data.frame("Explained-Variance"=c(as.numeric(CCA.clr.log.25.25.100.2.2)*100,
                                                                                                as.numeric(CCA.clr.log.25.25.100.4.4)*100,
                                                                                                as.numeric(CCA.clr.log.25.25.100.6.6)*100,
                                                                                            as.numeric(CCA.ilr.log.25.25.100.2.2)*100,
                                                                                            as.numeric(CCA.ilr.log.25.25.100.4.4)*100,
                                                                                            as.numeric(CCA.ilr.log.25.25.100.6.6)*100,
                                                                                            as.numeric(CCA.alpha.log.25.25.100.2.2)*100,
                                                                                            as.numeric(CCA.alpha.log.25.25.100.4.4)*100,
                                                                                            as.numeric(CCA.alpha.log.25.25.100.6.6)*100),
                                                                                            "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=3000),
                                                                         "No_Associations" = as.factor(rep(rep(c(2,4,6), each=1000),3)))



df.results.PLS.Reg.multivariate.factorization.based.25.25.100.various.associations = data.frame("Explained-Variance"=c(as.numeric(PLS.Reg.clr.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(PLS.Reg.clr.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(PLS.Reg.clr.log.25.25.100.6.6)*100,
                                                                                                                   as.numeric(PLS.Reg.ilr.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(PLS.Reg.ilr.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(PLS.Reg.ilr.log.25.25.100.6.6)*100,
                                                                                                                   as.numeric(PLS.Reg.alpha.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(PLS.Reg.alpha.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(PLS.Reg.alpha.log.25.25.100.6.6)*100),
                                                                                            "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=3000),
                                                                                            "No_Associations" = as.factor(rep(rep(c(2,4,6), each=1000),3)))




df.results.PLS.Canonical.multivariate.factorization.based.25.25.100.various.associations = data.frame("Explained-Variance"=c(as.numeric(PLS.Canonical.clr.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(PLS.Canonical.clr.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(PLS.Canonical.clr.log.25.25.100.6.6)*100,
                                                                                                                   as.numeric(PLS.Canonical.ilr.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(PLS.Canonical.ilr.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(PLS.Canonical.ilr.log.25.25.100.6.6)*100,
                                                                                                                   as.numeric(PLS.Canonical.alpha.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(PLS.Canonical.alpha.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(PLS.Canonical.alpha.log.25.25.100.6.6)*100),
                                                                                            "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=3000),
                                                                                            "No_Associations" = as.factor(rep(rep(c(2,4,6), each=1000),3)))



df.results.RDA.multivariate.factorization.based.25.25.100.various.associations = data.frame("Explained-Variance"=c(as.numeric(RDA.clr.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(RDA.clr.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(RDA.clr.log.25.25.100.6.6)*100,
                                                                                                                   as.numeric(RDA.ilr.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(RDA.ilr.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(RDA.ilr.log.25.25.100.6.6)*100,
                                                                                                                   as.numeric(RDA.alpha.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(RDA.alpha.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(RDA.alpha.log.25.25.100.6.6)*100),
                                                                                            "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=3000),
                                                                                            "No_Associations" = as.factor(rep(rep(c(2,4,6), each=1000),3)))


df.results.cca.multivariate.factorization.based.25.25.100.various.associations = data.frame("Explained-Variance"=c(as.numeric(CCA.clr.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(CCA.clr.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(CCA.clr.log.25.25.100.6.6)*100,
                                                                                                                   as.numeric(CCA.ilr.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(CCA.ilr.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(CCA.ilr.log.25.25.100.6.6)*100,
                                                                                                                   as.numeric(CCA.alpha.log.25.25.100.2.2)*100,
                                                                                                                   as.numeric(CCA.alpha.log.25.25.100.4.4)*100,
                                                                                                                   as.numeric(CCA.alpha.log.25.25.100.6.6)*100),
                                                                                            "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=3000),
                                                                                            "No_Associations" = as.factor(rep(rep(c(2,4,6), each=1000),3)))



df.results.MOFA.multivariate.factorization.based.25.25.100.various.associations = data.frame("Explained-Variance"=c(sapply(1:1000, function(X) MOFA.clr.log.25.25.100.2.2[[X]]$r2_total$group1[1]),
                                                                                                                    sapply(1:1000, function(X) MOFA.clr.log.25.25.100.4.4[[X]]$r2_total$group1[1]),
                                                                                                                    sapply(1:1000, function(X) MOFA.clr.log.25.25.100.6.6[[X]]$r2_total$group1[1]),
                                                                                                                    sapply(1:1000, function(X) MOFA.ilr.log.25.25.100.2.2[[X]]$r2_total$group1[1]),
                                                                                                                    sapply(1:1000, function(X) MOFA.ilr.log.25.25.100.4.4[[X]]$r2_total$group1[1]),
                                                                                                                    sapply(1:1000, function(X) MOFA.ilr.log.25.25.100.6.6[[X]]$r2_total$group1[1]),
                                                                                                                    sapply(1:1000, function(X) MOFA.alpha.log.25.25.100.2.2[[X]]$r2_total$group1[1]),
                                                                                                                    sapply(1:1000, function(X) MOFA.alpha.log.25.25.100.4.4[[X]]$r2_total$group1[1]),
                                                                                                                    sapply(1:1000, function(X) MOFA.alpha.log.25.25.100.6.6[[X]]$r2_total$group1[1])),
                                                                                                "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=3000),
                                                                                                "No_Associations" = as.factor(rep(rep(c(2,4,6), each=1000),3)))



ggarrange(ggplot(df.results.cca.multivariate.factorization.based.25.25.100.various.associations, aes(x=No_Associations, y=Explained.Variance, fill=No_Associations))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("#Associations")+ylim(c(0,100)),
ggplot(df.results.PLS.Reg.multivariate.factorization.based.25.25.100.various.associations, aes(x=No_Associations, y=Explained.Variance, fill=No_Associations))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("#Associations")+ylim(c(0,100)),
ggplot(df.results.PLS.Canonical.multivariate.factorization.based.25.25.100.various.associations, aes(x=No_Associations, y=Explained.Variance, fill=No_Associations))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("#Associations")+ylim(c(0,100)),
ggplot(df.results.RDA.multivariate.factorization.based.25.25.100.various.associations, aes(x=No_Associations, y=Explained.Variance, fill=No_Associations))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("#Associations")+ylim(c(0,100)),
ggplot(df.results.MOFA.multivariate.factorization.based.25.25.100.various.associations, aes(x=No_Associations, y=Explained.Variance, fill=No_Associations))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("#Associations")+ylim(c(0,100)),nrow=5,labels = c("A","B","C","D","E"))


CCA.clr.log.25.25.100.4.4.neg02 = readRDS("data\\results_CCA\\clr\\25_25_100\\redundancy_cca_clr_microbiome_log_metabolome_25_25_100indiv_4_4_neg0.2.RDS")
CCA.clr.log.25.25.100.4.4.neg01 = readRDS("data\\results_CCA\\clr\\25_25_100\\redundancy_cca_clr_microbiome_log_metabolome_25_25_100indiv_4_4_neg0.1.RDS")
CCA.clr.log.25.25.100.4.4.01 = readRDS("data\\results_CCA\\clr\\25_25_100\\redundancy_cca_clr_microbiome_log_metabolome_25_25_100indiv_4_4_0.1.RDS")

CCA.ilr.log.25.25.100.4.4.neg02 = readRDS("data\\results_CCA\\ilr\\25_25_100\\redundancy_cca_ilr_microbiome_log_metabolome_25_25_100indiv_4_4_effectneg0.2.RDS")
CCA.ilr.log.25.25.100.4.4.neg01 = readRDS("data\\results_CCA\\ilr\\25_25_100\\redundancy_cca_ilr_microbiome_log_metabolome_25_25_100indiv_4_4_effectneg0.1.RDS")
CCA.ilr.log.25.25.100.4.4.01 = readRDS("data\\results_CCA\\ilr\\25_25_100\\redundancy_cca_ilr_microbiome_log_metabolome_25_25_100indiv_4_4_effect0.1.RDS")

CCA.alpha.log.25.25.100.4.4.neg02 = readRDS("data\\results_CCA\\alpha\\25_25_100\\redundancy_cca_alpha_microbiome_log_metabolome_25_25_100indiv_4_4_effectneg0.2.RDS")
CCA.alpha.log.25.25.100.4.4.neg01 = readRDS("data\\results_CCA\\alpha\\25_25_100\\redundancy_cca_alpha_microbiome_log_metabolome_25_25_100indiv_4_4_effectneg0.1.RDS")
CCA.alpha.log.25.25.100.4.4.01 = readRDS("data\\results_CCA\\alpha\\25_25_100\\redundancy_cca_alpha_microbiome_log_metabolome_25_25_100indiv_4_4_effect0.1.RDS")


df.results.cca.multivariate.factorization.based.25.25.100.various.effects = data.frame("Explained-Variance"=c(as.numeric(CCA.clr.log.25.25.100.4.4.neg02)*100,
                                                                                                                   as.numeric(CCA.clr.log.25.25.100.4.4.neg01)*100,
                                                                                                                   as.numeric(CCA.clr.log.25.25.100.4.4.01)*100,
                                                                                                                   as.numeric(CCA.clr.log.25.25.100.4.4)*100,
                                                                                                              as.numeric(CCA.ilr.log.25.25.100.4.4.neg02)*100,
                                                                                                              as.numeric(CCA.ilr.log.25.25.100.4.4.neg01)*100,
                                                                                                              as.numeric(CCA.ilr.log.25.25.100.4.4.01)*100,
                                                                                                              as.numeric(CCA.ilr.log.25.25.100.4.4)*100,
                                                                                                              as.numeric(CCA.alpha.log.25.25.100.4.4.neg02)*100,
                                                                                                              as.numeric(CCA.alpha.log.25.25.100.4.4.neg01)*100,
                                                                                                              as.numeric(CCA.alpha.log.25.25.100.4.4.01)*100,
                                                                                                              as.numeric(CCA.alpha.log.25.25.100.4.4)*100),
                                                                                            "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=4000),
                                                                                            "Association_Strength" = factor(rep(rep(c("-0.2","-0.1","0.1","0.2"), each=1000),3), levels = c("-0.2","-0.1","0.1","0.2")))


PLS.Reg.clr.log.25.25.100.4.4.neg02 = readRDS("data\\results_PLS\\clr\\25_25_100\\regression\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effectneg0_2.RDS")
PLS.Reg.clr.log.25.25.100.4.4.neg01 = readRDS("data\\results_PLS\\clr\\25_25_100\\regression\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effectneg0_1.RDS")
PLS.Reg.clr.log.25.25.100.4.4.01 = readRDS("data\\results_PLS\\clr\\25_25_100\\regression\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effect0_1.RDS")

PLS.Reg.ilr.log.25.25.100.4.4.neg02 = readRDS("data\\results_PLS\\ilr\\25_25_100\\regression\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effectneg0_2.RDS")
PLS.Reg.ilr.log.25.25.100.4.4.neg01 = readRDS("data\\results_PLS\\ilr\\25_25_100\\regression\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effectneg0_1.RDS")
PLS.Reg.ilr.log.25.25.100.4.4.01 = readRDS("data\\results_PLS\\ilr\\25_25_100\\regression\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effect0_1.RDS")

PLS.Reg.alpha.log.25.25.100.4.4.neg02 = readRDS("data\\results_PLS\\alpha\\25_25_100\\regression\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effectneg0_2.RDS")
PLS.Reg.alpha.log.25.25.100.4.4.neg01 = readRDS("data\\results_PLS\\alpha\\25_25_100\\regression\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effectneg0_1.RDS")
PLS.Reg.alpha.log.25.25.100.4.4.01 = readRDS("data\\results_PLS\\alpha\\25_25_100\\regression\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_regression_4_4_effect0_1.RDS")


df.results.PLS.Reg.multivariate.factorization.based.25.25.100.various.effects = data.frame("Explained-Variance"=c(as.numeric(PLS.Reg.clr.log.25.25.100.4.4.neg02)*100,
                                                                                                              as.numeric(PLS.Reg.clr.log.25.25.100.4.4.neg01)*100,
                                                                                                              as.numeric(PLS.Reg.clr.log.25.25.100.4.4.01)*100,
                                                                                                              as.numeric(PLS.Reg.clr.log.25.25.100.4.4)*100,
                                                                                                              as.numeric(PLS.Reg.ilr.log.25.25.100.4.4.neg02)*100,
                                                                                                              as.numeric(PLS.Reg.ilr.log.25.25.100.4.4.neg01)*100,
                                                                                                              as.numeric(PLS.Reg.ilr.log.25.25.100.4.4.01)*100,
                                                                                                              as.numeric(PLS.Reg.ilr.log.25.25.100.4.4)*100,
                                                                                                              as.numeric(PLS.Reg.alpha.log.25.25.100.4.4.neg02)*100,
                                                                                                              as.numeric(PLS.Reg.alpha.log.25.25.100.4.4.neg02)*100,
                                                                                                              as.numeric(PLS.Reg.alpha.log.25.25.100.4.4.01)*100,
                                                                                                              as.numeric(PLS.Reg.alpha.log.25.25.100.4.4)*100),
                                                                                       "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=4000),
                                                                                       "Association_Strength" = factor(rep(rep(c("-0.2","-0.1","0.1","0.2"), each=1000),3), levels = c("-0.2","-0.1","0.1","0.2")))

PLS.Canonical.clr.log.25.25.100.4.4.neg02 = readRDS("data\\results_PLS\\clr\\25_25_100\\canonical\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effectneg0_2.RDS")
PLS.Canonical.clr.log.25.25.100.4.4.neg01 = readRDS("data\\results_PLS\\clr\\25_25_100\\canonical\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effectneg0_1.RDS")
PLS.Canonical.clr.log.25.25.100.4.4.01 = readRDS("data\\results_PLS\\clr\\25_25_100\\canonical\\redundancy_pls_clr_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effect0_1.RDS")

PLS.Canonical.ilr.log.25.25.100.4.4.neg02 = readRDS("data\\results_PLS\\ilr\\25_25_100\\canonical\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effectneg0_2.RDS")
PLS.Canonical.ilr.log.25.25.100.4.4.neg01 = readRDS("data\\results_PLS\\ilr\\25_25_100\\canonical\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effectneg0_1.RDS")
PLS.Canonical.ilr.log.25.25.100.4.4.01 = readRDS("data\\results_PLS\\ilr\\25_25_100\\canonical\\redundancy_pls_ilr_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effect0_1.RDS")

PLS.Canonical.alpha.log.25.25.100.4.4.neg02 = readRDS("data\\results_PLS\\alpha\\25_25_100\\canonical\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effectneg0_2.RDS")
PLS.Canonical.alpha.log.25.25.100.4.4.neg01 = readRDS("data\\results_PLS\\alpha\\25_25_100\\canonical\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effectneg0_1.RDS")
PLS.Canonical.alpha.log.25.25.100.4.4.01 = readRDS("data\\results_PLS\\alpha\\25_25_100\\canonical\\redundancy_pls_alpha_microbiome_log_metabolome_25_25_100indiv_canonical_4_4_effect0_1.RDS")


df.results.PLS.Canonical.multivariate.factorization.based.25.25.100.various.effects = data.frame("Explained-Variance"=c(as.numeric(PLS.Canonical.clr.log.25.25.100.4.4.neg02)*100,
                                                                                                                  as.numeric(PLS.Canonical.clr.log.25.25.100.4.4.neg01)*100,
                                                                                                                  as.numeric(PLS.Canonical.clr.log.25.25.100.4.4.01)*100,
                                                                                                                  as.numeric(PLS.Canonical.clr.log.25.25.100.4.4)*100,
                                                                                                                  as.numeric(PLS.Canonical.ilr.log.25.25.100.4.4.neg02)*100,
                                                                                                                  as.numeric(PLS.Canonical.ilr.log.25.25.100.4.4.neg01)*100,
                                                                                                                  as.numeric(PLS.Canonical.ilr.log.25.25.100.4.4.01)*100,
                                                                                                                  as.numeric(PLS.Canonical.ilr.log.25.25.100.4.4)*100,
                                                                                                                  as.numeric(PLS.Canonical.alpha.log.25.25.100.4.4.neg02)*100,
                                                                                                                  as.numeric(PLS.Canonical.alpha.log.25.25.100.4.4.neg02)*100,
                                                                                                                  as.numeric(PLS.Canonical.alpha.log.25.25.100.4.4.01)*100,
                                                                                                                  as.numeric(PLS.Canonical.alpha.log.25.25.100.4.4)*100),
                                                                                           "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=4000),
                                                                                           "Association_Strength" = factor(rep(rep(c("-0.2","-0.1","0.1","0.2"), each=1000),3), levels = c("-0.2","-0.1","0.1","0.2")))

RDA.clr.log.25.25.100.4.4.neg02 = readRDS("data\\results_RDA\\clr\\25_25_100\\explained_variance_clr_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effectneg0_2.RDS")
RDA.clr.log.25.25.100.4.4.neg01 = readRDS("data\\results_RDA\\clr\\25_25_100\\explained_variance_clr_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effectneg0_1.RDS")
RDA.clr.log.25.25.100.4.4.01 = readRDS("data\\results_RDA\\clr\\25_25_100\\explained_variance_clr_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effect0_1.RDS")

RDA.ilr.log.25.25.100.4.4.neg02 = readRDS("data\\results_RDA\\ilr\\25_25_100\\explained_variance_ilr_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effectneg0_2.RDS")
RDA.ilr.log.25.25.100.4.4.neg01 = readRDS("data\\results_RDA\\ilr\\25_25_100\\explained_variance_ilr_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effectneg0_1.RDS")
RDA.ilr.log.25.25.100.4.4.01 = readRDS("data\\results_RDA\\ilr\\25_25_100\\explained_variance_ilr_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effect0_1.RDS")

RDA.alpha.log.25.25.100.4.4.neg02 = readRDS("data\\results_RDA\\alpha\\25_25_100\\explained_variance_alpha_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effectneg0_2.RDS")
RDA.alpha.log.25.25.100.4.4.neg01 = readRDS("data\\results_RDA\\alpha\\25_25_100\\explained_variance_alpha_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effectneg0_1.RDS")
RDA.alpha.log.25.25.100.4.4.01 = readRDS("data\\results_RDA\\alpha\\25_25_100\\explained_variance_alpha_microbiome_log_metabolome_rda_25_25_100indiv_4_4_effect0_1.RDS")


df.results.RDA.multivariate.factorization.based.25.25.100.various.effects = data.frame("Explained-Variance"=c(as.numeric(RDA.clr.log.25.25.100.4.4.neg02)*100,
                                                                                                              as.numeric(RDA.clr.log.25.25.100.4.4.neg01)*100,
                                                                                                              as.numeric(RDA.clr.log.25.25.100.4.4.01)*100,
                                                                                                              as.numeric(RDA.clr.log.25.25.100.4.4)*100,
                                                                                                              as.numeric(RDA.ilr.log.25.25.100.4.4.neg02)*100,
                                                                                                              as.numeric(RDA.ilr.log.25.25.100.4.4.neg01)*100,
                                                                                                              as.numeric(RDA.ilr.log.25.25.100.4.4.01)*100,
                                                                                                              as.numeric(RDA.ilr.log.25.25.100.4.4)*100,
                                                                                                              as.numeric(RDA.alpha.log.25.25.100.4.4.neg02)*100,
                                                                                                              as.numeric(RDA.alpha.log.25.25.100.4.4.neg01)*100,
                                                                                                              as.numeric(RDA.alpha.log.25.25.100.4.4.01)*100,
                                                                                                              as.numeric(RDA.alpha.log.25.25.100.4.4)*100),
                                                                                       "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=4000),
                                                                                       "Association_Strength" = factor(rep(rep(c("-0.2","-0.1","0.1","0.2"), each=1000),3), levels = c("-0.2","-0.1","0.1","0.2")))

MOFA.clr.log.25.25.100.4.4.neg02 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4_effectneg0.2.RDS")
MOFA.clr.log.25.25.100.4.4.neg01 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4_effectneg0.1.RDS")
MOFA.clr.log.25.25.100.4.4.01 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4_effect0.1.RDS")

MOFA.ilr.log.25.25.100.4.4.neg02 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4_effectneg0.2.RDS")
MOFA.ilr.log.25.25.100.4.4.neg01 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4_effectneg0.1.RDS")
MOFA.ilr.log.25.25.100.4.4.01 = readRDS("data\\results_MOFA2\\ilr\\25_25_100\\MOFA2_ilr_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4_effect0.1.RDS")

MOFA.alpha.log.25.25.100.4.4.neg02 = readRDS("data\\results_MOFA2\\alpha\\25_25_100\\MOFA2_alpha_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4_effectneg0.2.RDS")
MOFA.alpha.log.25.25.100.4.4.neg01 = readRDS("data\\results_MOFA2\\alpha\\25_25_100\\MOFA2_alpha_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4_effectneg0.1.RDS")
MOFA.alpha.log.25.25.100.4.4.01 = readRDS("data\\results_MOFA2\\alpha\\25_25_100\\MOFA2_alpha_microbiome_log_metabolome_Scenario_25_25_100indiv_4_4_effect0.1.RDS")


df.results.MOFA.multivariate.factorization.based.25.25.100.various.effects = data.frame("Explained-Variance"=c(sapply(1:1000, function(X) MOFA.clr.log.25.25.100.4.4.neg02[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.clr.log.25.25.100.4.4.neg01[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.clr.log.25.25.100.4.4.01[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.clr.log.25.25.100.4.4[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.ilr.log.25.25.100.4.4.neg02[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.ilr.log.25.25.100.4.4.neg01[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.ilr.log.25.25.100.4.4.01[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.ilr.log.25.25.100.4.4[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.alpha.log.25.25.100.4.4.neg02[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.alpha.log.25.25.100.4.4.neg01[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.alpha.log.25.25.100.4.4.01[[X]]$r2_total$group1[1]),
                                                                                                               sapply(1:1000, function(X) MOFA.alpha.log.25.25.100.4.4[[X]]$r2_total$group1[1])),
                                                                                       "Normalization" = rep(c("CLR", "ILR", "Alpha"), each=4000),
                                                                                       "Association_Strength" = factor(rep(rep(c("-0.2","-0.1","0.1","0.2"), each=1000),3), levels = c("-0.2","-0.1","0.1","0.2")))

df.results.PLS.Reg.multivariate.factorization.based.25.25.100.various.effects

ggarrange(ggplot(df.results.cca.multivariate.factorization.based.25.25.100.various.effects, aes(x=Association_Strength, y=Explained.Variance, fill=Association_Strength))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("Association Strength")+ylim(c(0,100)),
          ggplot(df.results.PLS.Reg.multivariate.factorization.based.25.25.100.various.effects, aes(x=Association_Strength, y=Explained.Variance, fill=Association_Strength))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("Association Strength")+ylim(c(0,100)),
          ggplot(df.results.PLS.Canonical.multivariate.factorization.based.25.25.100.various.effects, aes(x=Association_Strength, y=Explained.Variance, fill=Association_Strength))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("Association Strength")+ylim(c(0,100)),
          ggplot(df.results.RDA.multivariate.factorization.based.25.25.100.various.effects, aes(x=Association_Strength, y=Explained.Variance, fill=Association_Strength))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("Association Strength")+ylim(c(0,100)),
          ggplot(df.results.MOFA.multivariate.factorization.based.25.25.100.various.effects, aes(x=Association_Strength, y=Explained.Variance, fill=Association_Strength))+geom_boxplot()+facet_grid(.~Normalization)+theme_bw()+theme(legend.position = "none")+ylab("% Explained Variance")+xlab("Association Strength")+ylim(c(0,100)), labels = c("A","B","C","D","E"), nrow=5)


#Individual Associations

spearman.null.clr.log.25.25.100 = readRDS("data\\Correlation\\25_25_100\\Spearman\\null\\pvalues_null_spearman_clr_microbiome_metabolome_25_25_100indiv.RDS")
pearson.null.clr.log.25.25.100 = readRDS("data\\Correlation\\25_25_100\\Pearson\\null\\pvalues_null_pearson_clr_microbiome_log_metabolome_25_25_100indiv.RDS")

spearman.alternative.clr.log.25.25.100 = readRDS("data\\Correlation\\25_25_100\\Spearman\\alternative\\pvalues_alternative_spearman_clr_microbiome_metabolome_25_25_100indiv.RDS")
pearson.alternative.clr.log.25.25.100 = readRDS("data\\Correlation\\25_25_100\\Pearson\\alternative\\pvalues_alternative_pearson_clr_microbiome_log_metabolome_25_25_100indiv.RDS")

spearman.null.clr.log.100.100 = readRDS("data\\Correlation\\100_100_500\\Spearman\\null\\pvalues_null_spearman_clr_microbiome_metabolome_100_100_500indiv.RDS")
pearson.null.clr.log.100.100.500 = readRDS("data\\Correlation\\100_100_500\\Pearson\\null\\pvalues_null_pearson_clr_microbiome_log_metabolome_100_100_500indiv.RDS")

spearman.alternative.clr.log.100.100.500 = readRDS("data\\Correlation\\100_100_500\\Spearman\\alternative\\pvalues_alternative_spearman_clr_microbiome_metabolome_100_100_500indiv.RDS")
pearson.alternative.clr.log.100.100.500 = readRDS("data\\Correlation\\100_100_500\\Pearson\\alternative\\pvalues_alternative_pearson_clr_microbiome_log_metabolome_100_100_500indiv.RDS")


df.power.correlation = data.frame("Scenario"=factor(rep(c("Low-Dimensional", "High-Dimensional"), each=2), levels = c("Low-Dimensional", "High-Dimensional")), 
                                  "Correlation"=rep(c("Pearson","Spearman"),2),
                                  "Power" = c(sum(unlist(lapply(1:1000, function(x) sapply(pearson.alternative.clr.log.25.25.100[[x]], ACAT::ACAT))) <= 0.05)/25000,
                                              sum(unlist(lapply(1:1000, function(x) sapply(spearman.alternative.clr.log.25.25.100[[x]], ACAT::ACAT))) <= 0.05)/25000,
                                              sum(unlist(lapply(1:1000, function(x) sapply(pearson.alternative.clr.log.100.100.500[[x]], ACAT::ACAT))) <= 0.05)/100000,
                                              sum(unlist(lapply(1:1000, function(x) sapply(spearman.alternative.clr.log.100.100.500[[x]], ACAT::ACAT))) <= 0.05)/100000))

ggarrange(gg_qqplot_facet_grid(list("Pearson"=unlist(lapply(1:1000, function(x) sapply(pearson.null.clr.log.25.25.100[[x]], ACAT::ACAT))), 
                                           "Spearman"=unlist(lapply(1:1000, function(x) sapply(spearman.null.clr.log.25.25.100[[x]], ACAT::ACAT)))))+labs(colour="Correlation")+ggtitle("Low-Dimensional"),
          gg_qqplot_facet_grid(list("Pearson"=unlist(lapply(1:1000, function(x) sapply(pearson.null.clr.log.100.100.500[[x]], ACAT::ACAT))), 
                                    "Spearman"=unlist(lapply(1:1000, function(x) sapply(spearman.null.clr.log.100.100[[x]], ACAT::ACAT)))))+labs(colour="Correlation")+ggtitle("High-Dimensional"),
          ggplot(df.power.correlation, aes(x=Correlation, y=Power, fill=Correlation)) + geom_bar(stat="identity", position="dodge",color="black")+theme_bw()+theme(legend.position = "none")+facet_grid(.~Scenario)+ylim(c(0,1)), labels=c("A","B","C"), ncol=3)

"data\\power.p.values.100.100.500indiv.MiRKAT.alpha.original.RDS"
pvalues.null.alpha.log.MiRKAT.25.25.100 = readRDS("data\\null.p.values.25.25.100indiv.MiRKAT.alpha.log.RDS")
pvalues.null.alpha.original.MiRKAT.25.25.100 = readRDS("data\\null.p.values.25.25.100indiv.MiRKAT.alpha.original.RDS")

pvalues.null.alpha.log.MiRKAT.100.100.500 = readRDS("data\\null.p.values.100.100.500indiv.MiRKAT.alpha.log.RDS")
pvalues.null.alpha.original.MiRKAT.100.100.500 = readRDS("data\\null.p.values.100.100.500indiv.MiRKAT.alpha.original.RDS")

pvalues.null.MiRKAT.100.100.500 = readRDS("data\\list.results.null.MiRKAT.100.100.500.RDS")
pvalues.null.MiRKAT.25.25.100 = readRDS("data\\list.results.null.MiRKAT.25.25.100.RDS")

ggarrange(gg_qqplot_facet_grid(list("CLR-Log"=unlist(pvalues.null.MiRKAT.25.25.100$null.p.values.25.25.100indiv.MiRKAT.clr.log),"Alpha-Log"=unlist(pvalues.null.alpha.log.MiRKAT.25.25.100), "ILR-Original"=unlist(pvalues.null.MiRKAT.25.25.100$null.p.values.25.25.100indiv.MiRKAT.ilr.original),
                          "CLR-Original"=unlist(pvalues.null.MiRKAT.25.25.100$null.p.values.25.25.100indiv.MiRKAT.clr.original),"Alpha-Original"=unlist(pvalues.null.alpha.original.MiRKAT.25.25.100),"Original-Log"=unlist(pvalues.null.MiRKAT.25.25.100$list.results.null.MiRKAT.25.25.100.original.log)))+ggtitle("Low-Dimensional")+labs(colour="Normalizations"),
          gg_qqplot_facet_grid(list("CLR-Log"=unlist(pvalues.null.MiRKAT.100.100.500$null.p.values.100.100.500indiv.MiRKAT.clr.log),"Alpha-Log"=unlist(pvalues.null.alpha.log.MiRKAT.100.100.500), "ILR-Original"=unlist(pvalues.null.MiRKAT.100.100.500$null.p.values.100.100.500indiv.MiRKAT.ilr.original),
                                      "CLR-Original"=unlist(pvalues.null.MiRKAT.100.100.500$null.p.values.100.100.500indiv.MiRKAT.clr.original),"Alpha-Original"=unlist(pvalues.null.alpha.original.MiRKAT.100.100.500),"Original-Log"=unlist(pvalues.null.MiRKAT.100.100.500$null.p.values.100.100.500indiv.MiRKAT.original.log)))+ggtitle("High-Dimensional")+labs(colour="Normalizations"), labels=c("A","B"))


pvalues.alter.alpha.log.MiRKAT.25.25.100 = unlist(readRDS("data\\power.p.values.25.25.100indiv.MiRKAT.alpha.log.RDS"))
pvalues.alter.alpha.original.MiRKAT.25.25.100 = unlist(readRDS("data\\power.p.values.25.25.100indiv.MiRKAT.alpha.original.RDS"))

pvalues.alter.clr.log.MiRKAT.25.25.100 = unlist(readRDS("data\\power.p.values.25.25.100indiv.MiRKAT.clr.log.RDS"))
pvalues.alter.clr.original.MiRKAT.25.25.100 = unlist(readRDS("data\\power.p.values.25.25.100indiv.MiRKAT.clr.original.RDS"))

pvalues.alter.ilr.original.MiRKAT.25.25.100 = unlist(readRDS("data\\power.p.values.25.25.100indiv.MiRKAT.ilr.original.RDS"))

pvalues.alter.original.log.MiRKAT.25.25.100 = unlist(readRDS("data\\power.p.values.25.25.100indiv.MiRKAT.original.log.RDS"))


pvalues.alter.alpha.log.MiRKAT.100.100.500 = unlist(readRDS("data\\power.p.values.100.100.500indiv.MiRKAT.alpha.log.RDS"))
pvalues.alter.alpha.original.MiRKAT.100.100.500 = unlist(readRDS("data\\power.p.values.100.100.500indiv.MiRKAT.alpha.original.RDS"))

pvalues.alter.clr.log.MiRKAT.100.100.500 = unlist(readRDS("data\\power.p.values.100.100.500indiv.MiRKAT.clr.log.RDS"))
pvalues.alter.clr.original.MiRKAT.100.100.500 = unlist(readRDS("data\\power.p.values.100.100.500indiv.MiRKAT.clr.original.RDS"))

pvalues.alter.ilr.original.MiRKAT.100.100.500 = unlist(readRDS("data\\power.p.values.100.100.500indiv.MiRKAT.ilr.original.RDS"))

pvalues.alter.original.log.MiRKAT.100.100.500 = unlist(readRDS("data\\power.p.values.100.100.500indiv.MiRKAT.original.log.RDS"))

df.power.MiRKAT.25.25.100 = data.frame("Power" = c(sum(pvalues.alter.clr.log.MiRKAT.25.25.100<=0.05)/25000,
                                                   sum(pvalues.alter.alpha.log.MiRKAT.25.25.100<=0.05)/25000,
                                                   sum(pvalues.alter.ilr.original.MiRKAT.25.25.100<=0.05)/25000,
                                                   sum(pvalues.alter.clr.original.MiRKAT.25.25.100<=0.05)/25000,
                                                   sum(pvalues.alter.alpha.original.MiRKAT.25.25.100<=0.05)/25000,
                                                   sum(pvalues.alter.original.log.MiRKAT.25.25.100<=0.05)/25000), 
                                       "Normalizations"= c("CLR-Log", "Alpha-Log", "ILR-Original","CLR-Original", "Alpha-Original", "Original-Log"),
                                       "Scenario"= "Low-Dimensional")

df.power.MiRKAT.100.100.500 = data.frame("Power" = c(sum(pvalues.alter.clr.log.MiRKAT.100.100.500<=0.05)/100000,
                                                     sum(pvalues.alter.alpha.log.MiRKAT.100.100.500<=0.05)/100000,
                                                     sum(pvalues.alter.ilr.original.MiRKAT.100.100.500<=0.05)/100000,
                                                     sum(pvalues.alter.clr.original.MiRKAT.100.100.500<=0.05)/100000,
                                                     sum(pvalues.alter.alpha.original.MiRKAT.100.100.500<=0.05)/100000,
                                                     sum(pvalues.alter.original.log.MiRKAT.100.100.500<=0.05)/100000), 
                                        "Normalizations"= c("CLR-Log", "Alpha-Log", "ILR-Original","CLR-Original", "Alpha-Original", "Original-Log"),
                                        "Scenario"= "High-Dimensional")

ggarrange(ggarrange(gg_qqplot_facet_grid(list("CLR-Log"=unlist(pvalues.null.MiRKAT.25.25.100$null.p.values.25.25.100indiv.MiRKAT.clr.log),"Alpha-Log"=unlist(pvalues.null.alpha.log.MiRKAT.25.25.100), "ILR-Original"=unlist(pvalues.null.MiRKAT.25.25.100$null.p.values.25.25.100indiv.MiRKAT.ilr.original),
                                    "CLR-Original"=unlist(pvalues.null.MiRKAT.25.25.100$null.p.values.25.25.100indiv.MiRKAT.clr.original),"Alpha-Original"=unlist(pvalues.null.alpha.original.MiRKAT.25.25.100),"Original-Log"=unlist(pvalues.null.MiRKAT.25.25.100$list.results.null.MiRKAT.25.25.100.original.log)))+ggtitle("Low-Dimensional")+labs(colour="Normalizations"),
          gg_qqplot_facet_grid(list("CLR-Log"=unlist(pvalues.null.MiRKAT.100.100.500$null.p.values.100.100.500indiv.MiRKAT.clr.log),"Alpha-Log"=unlist(pvalues.null.alpha.log.MiRKAT.100.100.500), "ILR-Original"=unlist(pvalues.null.MiRKAT.100.100.500$null.p.values.100.100.500indiv.MiRKAT.ilr.original),
                                    "CLR-Original"=unlist(pvalues.null.MiRKAT.100.100.500$null.p.values.100.100.500indiv.MiRKAT.clr.original),"Alpha-Original"=unlist(pvalues.null.alpha.original.MiRKAT.100.100.500),"Original-Log"=unlist(pvalues.null.MiRKAT.100.100.500$null.p.values.100.100.500indiv.MiRKAT.original.log)))+ggtitle("High-Dimensional")+labs(colour="Normalizations"),labels=c("A","B"), ncol=2),
          ggarrange(ggplot(df.power.MiRKAT, aes(x=Normalizations, y=Power, fill=Normalizations)) + geom_bar(stat="identity", position="dodge", color="black") + theme_bw() + theme(legend.position = "none")+ylim(c(0,1))+facet_grid(.~Scenario), labels = c("C")),nrow=2)
df.power.MiRKAT = rbind(df.power.MiRKAT.25.25.100,df.power.MiRKAT.100.100.500)
df.power.MiRKAT$Scenario = factor(df.power.MiRKAT$Scenario, levels = c("Low-Dimensional", "High-Dimensional"))
ggplot(df.power.MiRKAT, aes(x=Normalizations, y=Power, fill=Normalizations)) + geom_bar(stat="identity", position="dodge", color="black") + theme_bw() + theme(legend.position = "none")+ylim(c(0,1))+facet_grid(.~Scenario)


pvalues.null.Dir.log.25.25.100 = readRDS("data\\Regression\\25_25_100\\null\\pvalues_null_Dir_log_metabolome_25_25_100indiv.RDS")
pvalues.alternative.Dir.log.25.25.100 = readRDS("data\\Regression\\25_25_100\\alternative\\pvalues_alternative_Dir_log_metabolome_25_25_100indiv.RDS")


pvalues.null.log.clr.25.25.100 = readRDS("data\\Regression\\25_25_100\\null\\pvalues_null_loglinear_clr_microbiome_25_25_100indiv.RDS")
pvalues.alternative.lm.clr.log.25.25.100 = readRDS("data\\Regression\\25_25_100\\alternative\\pvalues_alternative_lm_clr_microbiome_log_metabolome_25_25_100indiv.RDS")



df.power.compositional.outcome = data.frame("Method"=c("Dirichlet","LM"),
                                            "Power"=c(sum(unlist(lapply(1:1000, function(x) sapply(pvalues.alternative.Dir.log.25.25.100[[x]], ACAT::ACAT))) <= 0.05)/20000,
                                                      sum(unlist(lapply(1:1000, function(x) sapply(pvalues.alternative.lm.clr.log.25.25.100[[x]], ACAT::ACAT))) <= 0.05)/25000))


ggarrange(gg_qqplot_facet_grid(list("LM"=unlist(lapply(1:1000, function(x) sapply(pvalues.null.log.clr.25.25.100[[x]], ACAT::ACAT))),"Dirichlet"=unlist(lapply(1:1000, function(x) sapply(pvalues.null.Dir.log.25.25.100[[x]], ACAT::ACAT))))),
          ggplot(df.power.compositional.outcome, aes(x=Method, y=Power, fill=Method)) + geom_bar(stat="identity", position="dodge", color="black") + theme_bw() + theme(legend.position = "none")+ylim(c(0,1)),labels = c("A","B"), ncol=2)

#Univariate feature selection methods

#Multivariate feature selection methods

simulated.data.100.100.500.alter = readRDS("\data\\data_alternative_100_100_500indiv.RDS")
simulated.data.25.25.100.alter = readRDS("data\\data_alternative_25_25_100indiv.RDS")


sCCA.none.log.25.25.100 = readRDS("data\\results_sCCA\\original\\25_25_100\\sCCA_noneMicrobiome_logMetabolome_25_25_100ind.RDS")
sCCA.none.none.25.25.100 = readRDS("data\\results_sCCA\\original\\25_25_100\\sCCA_noneMicrobiome_logMetabolome_25_25_100ind.RDS")

sCCA.none.log.100.100.500 = readRDS("data\\results_sCCA\\original\\100_100_500\\sCCA_noneMicrobiome_logMetabolome_100_100_500ind.RDS")
sCCA.none.none.100.100.500 = readRDS("data\\results_sCCA\\original\\100_100_500\\sCCA_noneMicrobiome_logMetabolome_100_100_500ind.RDS")


sPLS.Reg.none.log.25.25.100 = readRDS("data\\results_sPLS\\original\\25_25_100\\regression\\loadings_sPLSregression_Microbiome_log_Metabolome_25_25_100indiv.RDS")
sPLS.Reg.none.none.25.25.100 = readRDS("data\\results_sPLS\\original\\25_25_100\\regression\\loadings_sPLSregression_Microbiome_Metabolome_25_25_100indiv.RDS")

sPLS.Reg.none.log.100.100.500 = readRDS("data\\results_sPLS\\original\\100_100_500\\regression\\loadings_sPLSregression_Microbiome_log_Metabolome_100_100_500indiv.RDS")
sPLS.Reg.none.none.100.100.500 = readRDS("data\\results_sPLS\\original\\100_100_500\\regression\\loadings_sPLSregression_Microbiome_Metabolome_100_100_500indiv.RDS")

sPLS.Canonical.none.log.25.25.100 = readRDS("data\\results_sPLS\\original\\25_25_100\\canonical\\loadings_sPLScanonical_Microbiome_log_Metabolome_25_25_100indiv.RDS")
sPLS.Canonical.none.none.25.25.100 = readRDS("data\\results_sPLS\\original\\25_25_100\\canonical\\loadings_sPLScanonical_Microbiome_Metabolome_25_25_100indiv.RDS")

sPLS.Canonical.none.log.100.100.500 = readRDS("data\\results_sPLS\\original\\100_100_500\\canonical\\loadings_sPLScanonical_Microbiome_log_Metabolome_100_100_500indiv.RDS")
sPLS.Canonical.none.none.100.100.500 = readRDS("data\\results_sPLS\\original\\100_100_500\\canonical\\loadings_sPLScanonical_Microbiome_Metabolome_100_100_500indiv.RDS")

index.nonnull.sPLS.regression.none.log.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.Reg.none.log.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.Reg.none.log.25.25.100[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.none.log.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.Canonical.none.log.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.Canonical.none.log.25.25.100[[X]]$Y)!=0))
})


index.nonnull.sPLS.regression.none.none.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.Reg.none.none.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.Reg.none.none.25.25.100[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.none.none.25.25.100 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.Canonical.none.none.25.25.100[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.Canonical.none.none.25.25.100[[X]]$Y)!=0))
})



index.nonnull.sPLS.regression.none.log.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.Reg.none.log.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.Reg.none.log.100.100.500[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.none.log.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.Canonical.none.log.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.Canonical.none.log.100.100.500[[X]]$Y)!=0))
})


index.nonnull.sPLS.regression.none.none.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.Reg.none.none.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.Reg.none.none.100.100.500[[X]]$Y)!=0))
})

index.nonnull.sPLS.canonical.none.none.100.100.500 = lapply(1:1000, function(X){
  list("Microbiome"= which(rowSums(sPLS.Canonical.none.none.100.100.500[[X]]$X)!=0), 
       "Metabolome"=which(rowSums(sPLS.Canonical.none.none.100.100.500[[X]]$Y)!=0))
})


f1.Score.sCCA.none.log.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.none.log.100.100.500[[X]]$index$Microbiome,c(1:100)[-sCCA.none.log.100.100.500[[X]]$index$Microbiome],
                            sCCA.none.log.100.100.500[[X]]$index$Metabolome,c(1:100)[-sCCA.none.log.100.100.500[[X]]$index$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.none.log.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.none.log.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.none.log.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.none.log.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.none.log.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})



f1.Score.sPLS.canonical.none.log.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.none.log.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.none.log.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.none.log.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.none.log.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})


F1Score.all.methods.none.log.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.none.log.100.100.500,f1.Score.sPLS.regression.none.log.100.100.500,f1.Score.sPLS.canonical.none.log.100.100.500), "Scenario"="High-Dimensional", "Metric"="F1-Score")

f1.Score.sCCA.none.none.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.none.none.100.100.500[[X]]$index$Microbiome,c(1:100)[-sCCA.none.none.100.100.500[[X]]$index$Microbiome],
                            sCCA.none.none.100.100.500[[X]]$index$Metabolome,c(1:100)[-sCCA.none.none.100.100.500[[X]]$index$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.none.none.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.none.none.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.none.none.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.none.none.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.none.none.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})



f1.Score.sPLS.canonical.none.none.100.100.500 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.none.none.100.100.500[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.none.none.100.100.500[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.none.none.100.100.500[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.none.none.100.100.500[[X]]$Metabolome],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Microbiotes,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.100.100.500.alter$Index[[X]]$Associated$Metabolites,simulated.data.100.100.500.alter$Index[[X]]$Impacted$Metabolites )]))
})


F1Score.all.methods.none.none.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.none.none.100.100.500,f1.Score.sPLS.regression.none.none.100.100.500,f1.Score.sPLS.canonical.none.none.100.100.500), "Scenario"="High-Dimensional", "Metric"="F1-Score")


sparsity.all.methods.none.log.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.none.log.100.100.500[[X]]$index$Microbiome) +length(sCCA.none.log.100.100.500[[X]]$index$Metabolome))/200 )
                                                                                                                                         ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.none.log.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.none.log.100.100.500[[X]]$Metabolome))/200)
                                                                                                                                                                                                                                                                         ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.none.log.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.none.log.100.100.500[[X]]$Metabolome))/200)), "Scenario"="High-Dimensional", "Metric"="Sparsity")

sparsity.all.methods.none.none.100.100.500 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.none.none.100.100.500[[X]]$index$Microbiome) +length(sCCA.none.none.100.100.500[[X]]$index$Metabolome))/200 )
                                                                                                                                ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.none.none.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.none.none.100.100.500[[X]]$Metabolome))/200)
                                                                                                                                ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.none.none.100.100.500[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.none.none.100.100.500[[X]]$Metabolome))/200)), "Scenario"="High-Dimensional", "Metric"="Sparsity")

df.metrics.multivariate.log.none.100.100.500 = rbind(F1Score.all.methods.none.log.100.100.500,sparsity.all.methods.none.log.100.100.500)
df.metrics.multivariate.log.none.100.100.500$Metric = factor(df.metrics.multivariate.log.none.100.100.500$Metric, levels=c("Sparsity","F1-Score"))


df.metrics.multivariate.none.none.100.100.500 = rbind(F1Score.all.methods.none.none.100.100.500,sparsity.all.methods.none.none.100.100.500)
df.metrics.multivariate.none.none.100.100.500$Metric = factor(df.metrics.multivariate.none.none.100.100.500$Metric, levels=c("Sparsity","F1-Score"))


f1.Score.sCCA.none.log.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.none.log.25.25.100[[X]]$index$Microbiome,c(1:100)[-sCCA.none.log.25.25.100[[X]]$index$Microbiome],
                            sCCA.none.log.25.25.100[[X]]$index$Metabolome,c(1:100)[-sCCA.none.log.25.25.100[[X]]$index$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.none.log.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.none.log.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.none.log.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.none.log.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.none.log.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})



f1.Score.sPLS.canonical.none.log.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.none.log.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.none.log.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.none.log.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.none.log.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})


F1Score.all.methods.none.log.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.none.log.25.25.100,f1.Score.sPLS.regression.none.log.25.25.100,f1.Score.sPLS.canonical.none.log.25.25.100), "Scenario"="High-Dimensional", "Metric"="F1-Score")

f1.Score.sCCA.none.none.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(sCCA.none.none.25.25.100[[X]]$index$Microbiome,c(1:100)[-sCCA.none.none.25.25.100[[X]]$index$Microbiome],
                            sCCA.none.none.25.25.100[[X]]$index$Metabolome,c(1:100)[-sCCA.none.none.25.25.100[[X]]$index$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})

f1.Score.sPLS.regression.none.none.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.regression.none.none.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.regression.none.none.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.regression.none.none.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.regression.none.none.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})



f1.Score.sPLS.canonical.none.none.25.25.100 = sapply(1:1000, function(X) {
  f1.score(confusion.matrix(index.nonnull.sPLS.canonical.none.none.25.25.100[[X]]$Microbiome,c(1:100)[-index.nonnull.sPLS.canonical.none.none.25.25.100[[X]]$Microbiome],
                            index.nonnull.sPLS.canonical.none.none.25.25.100[[X]]$Metabolome,c(1:100)[-index.nonnull.sPLS.canonical.none.none.25.25.100[[X]]$Metabolome],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Microbiotes,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Microbiotes )],
                            c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites ),
                            c(1:100)[-c(simulated.data.25.25.100.alter$Index[[X]]$Associated$Metabolites,simulated.data.25.25.100.alter$Index[[X]]$Impacted$Metabolites )]))
})


F1Score.all.methods.none.none.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(f1.Score.sCCA.none.none.25.25.100,f1.Score.sPLS.regression.none.none.25.25.100,f1.Score.sPLS.canonical.none.none.25.25.100), "Scenario"="High-Dimensional", "Metric"="F1-Score")


sparsity.all.methods.none.log.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.none.log.25.25.100[[X]]$index$Microbiome) +length(sCCA.none.log.25.25.100[[X]]$index$Metabolome))/50 )
                                                                                                                                ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.none.log.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.none.log.25.25.100[[X]]$Metabolome))/50)
                                                                                                                                ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.none.log.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.none.log.25.25.100[[X]]$Metabolome))/50)), "Scenario"="High-Dimensional", "Metric"="Sparsity")

sparsity.all.methods.none.none.25.25.100 = data.frame("Method" = rep(c("sCCA", "sPLS-Reg", "sPLS-Can"), each=1000),"Value" = c(sapply(1:1000, function(X) (length(sCCA.none.none.25.25.100[[X]]$index$Microbiome) +length(sCCA.none.none.25.25.100[[X]]$index$Metabolome))/50 )
                                                                                                                                 ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.regression.none.none.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.regression.none.none.25.25.100[[X]]$Metabolome))/50)
                                                                                                                                 ,sapply(1:1000, function(X) (length(index.nonnull.sPLS.canonical.none.none.25.25.100[[X]]$Microbiome) + length(index.nonnull.sPLS.canonical.none.none.25.25.100[[X]]$Metabolome))/50)), "Scenario"="High-Dimensional", "Metric"="Sparsity")

df.metrics.multivariate.log.none.25.25.100 = rbind(F1Score.all.methods.none.log.25.25.100,sparsity.all.methods.none.log.25.25.100)
df.metrics.multivariate.log.none.25.25.100$Metric = factor(df.metrics.multivariate.log.none.25.25.100$Metric, levels=c("Sparsity","F1-Score"))


df.metrics.multivariate.none.none.25.25.100 = rbind(F1Score.all.methods.none.none.25.25.100,sparsity.all.methods.none.none.25.25.100)
df.metrics.multivariate.none.none.25.25.100$Metric = factor(df.metrics.multivariate.none.none.25.25.100$Metric, levels=c("Sparsity","F1-Score"))

ggarrange(ggarrange(ggplot(df.metrics.multivariate.log.none.25.25.100, aes(x=Method, fill=Method, y=Value*100))+geom_boxplot()+facet_grid(.~Metric)+theme_bw()+theme(legend.position = "none")+ylab("%")+ylim(c(0,100))+ggtitle("Original Microbiome-Log Metabolome"),
          ggplot(df.metrics.multivariate.none.none.25.25.100, aes(x=Method, fill=Method, y=Value*100))+geom_boxplot()+facet_grid(.~Metric)+theme_bw()+theme(legend.position = "none")+ylab("%")+ylim(c(0,100))+ggtitle("Original Microbiome-Original Metabolome"),labels = c("A","B"), ncol=2),

ggarrange(ggplot(df.metrics.multivariate.log.none.100.100.500, aes(x=Method, fill=Method, y=Value*100))+geom_boxplot()+facet_grid(.~Metric)+theme_bw()+theme(legend.position = "none")+ylab("%")+ylim(c(0,100)),
ggplot(df.metrics.multivariate.none.none.100.100.500, aes(x=Method, fill=Method, y=Value*100))+geom_boxplot()+facet_grid(.~Metric)+theme_bw()+theme(legend.position = "none")+ylab("%")+ylim(c(0,100)),labels = c("C","D"), ncol=2),nrow=2)
