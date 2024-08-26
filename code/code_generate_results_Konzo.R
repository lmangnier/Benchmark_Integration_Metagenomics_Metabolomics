#This code reproduces the results from the real-data application on Konzo
#Metagenomics and Metabolomics data are available upon request

library(compositions)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggsci)
library(igraph)
library(GGally)
library(vegan)

Microbiome.data.aff = read.table("Microbiome_unnormalized_affected.tsv", header=TRUE, sep="\t")
Microbiome.data.unaff = read.table("Microbiome_unnormalized_unaffected.tsv", header=TRUE, sep="\t")

Metabolome.data.aff = read.table("Metabolome_unnormalized_affected.tsv", header=TRUE, sep="\t")
Metabolome.data.unaff = read.table("Metabolome_unnormalized_unaffected.tsv", header=TRUE, sep="\t")

sample.in.Microbiome = which(colnames(Microbiome.data.aff)=="Sample")
sample.in.Metabolome = which(colnames(Metabolome.data.aff)=="Sample")

Metabolome.data.aff = Metabolome.data.aff[,-sample.in.Metabolome]
Microbiome.data.aff = Microbiome.data.aff[,-sample.in.Microbiome]

Metabolome.data.unaff = Metabolome.data.unaff[,-sample.in.Metabolome]
Microbiome.data.unaff = Microbiome.data.unaff[,-sample.in.Microbiome]

colnames.Microbiome = colnames(Microbiome.data.aff)
colnames.Metabolome = colnames(Metabolome.data.aff)

clr_microbiome_affected = clr(Microbiome.data.aff)
ilr_microbiome_affected = ilr(Microbiome.data.aff)
log_metabolome_affected = log(Metabolome.data.aff+1)

clr_microbiome_unaffected = clr(Microbiome.data.unaff)
ilr_microbiome_unaffected = ilr(Microbiome.data.unaff)
log_metabolome_unaffected = log(Metabolome.data.unaff+1)

sPLS.sCCA.Konzo.output.aff = readRDS("results_sCCA_sPLS_Konzo_affected.RDS")
non.null.loadings.species.KONZO.sPLS.aff = sPLS.sCCA.Konzo.output.aff$sPLS$loadings$X[rowSums(sPLS.sCCA.Konzo.output.aff$sPLS$loadings$X)!=0,]
non.null.loadings.metabo.KONZO.sPLS.aff = sPLS.sCCA.Konzo.output.aff$sPLS$loadings$Y[rowSums(sPLS.sCCA.Konzo.output.aff$sPLS$loadings$Y)!=0,]


plot.t = mixOmics::plotVar(sPLS.sCCA.Konzo.output.aff$sPLS)
df.sPLS.KONZO.loadings = data.frame("Type" = c(rep("Microorganisms", nrow(non.null.loadings.species.KONZO.sPLS.aff)), rep("Metabolites",nrow(non.null.loadings.metabo.KONZO.sPLS.aff))),
                                    "1st Component" = plot.t$x, "2nd Component"=plot.t$y)

Figure6CA = ggplot(df.sPLS.KONZO.loadings, aes(x=X1st.Component, y=X2nd.Component, colour=Type))+geom_point()+theme_bw()+scale_color_nejm()+xlab("1st Component")+ylab("2nd Component")+ggtitle("Affected")


sPLS.sCCA.Konzo.output.unaff = readRDS("results_sCCA_sPLS_Konzo_unaffected.RDS")
non.null.loadings.species.KONZO.sPLS.unaff = sPLS.sCCA.Konzo.output.unaff$sPLS$loadings$X[rowSums(sPLS.sCCA.Konzo.output.unaff$sPLS$loadings$X)!=0,]
non.null.loadings.metabo.KONZO.sPLS.unaff = sPLS.sCCA.Konzo.output.unaff$sPLS$loadings$Y[rowSums(sPLS.sCCA.Konzo.output.unaff$sPLS$loadings$Y)!=0,]

plot.tt = mixOmics::plotVar(sPLS.sCCA.Konzo.output.unaff$sPLS)
df.sPLS.KONZO.loadings.2 = data.frame("Type" = c(rep("Microorganisms", nrow(non.null.loadings.species.KONZO.sPLS.unaff)), rep("Metabolites",nrow(non.null.loadings.metabo.KONZO.sPLS.unaff))),
                                    "1st Component" = plot.tt$x, "2nd Component"=plot.tt$y)

Figure6CB = ggplot(df.sPLS.KONZO.loadings.2, aes(x=X1st.Component, y=X2nd.Component, colour=Type))+geom_point()+theme_bw()+scale_color_nejm()+xlab("1st Component")+ylab("2nd Component")+ggtitle("Unaffected")

Figure6C = Figure6CA + Figure6CB


sum(rownames(non.null.loadings.species.KONZO.sPLS.unaff)%in%rownames(non.null.loadings.species.KONZO.sPLS.aff))
sum(rownames(non.null.loadings.metabo.KONZO.sPLS.unaff)%in%rownames(non.null.loadings.metabo.KONZO.sPLS.aff))

l.aff = list()
for(m in rownames(sPLS.sCCA.Konzo.output.aff$sPLS$loadings$Y[rowSums(sPLS.sCCA.Konzo.output.aff$sPLS$loadings$Y)!=0,])){
  print(m)
  l.aff[[m]] = coda4microbiome::coda_glmnet(Microbiome.data.aff[,rownames(sPLS.sCCA.Konzo.output.aff$sPLS$loadings$X[rowSums(sPLS.sCCA.Konzo.output.aff$sPLS$loadings$X)!=0,])], log(Metabolome.data.aff[,m]))
  
}

l.unaff = list()
for(m in rownames(sPLS.sCCA.Konzo.output.unaff$sPLS$loadings$Y[rowSums(sPLS.sCCA.Konzo.output.unaff$sPLS$loadings$Y)!=0,])){
  print(m)
  l.unaff[[m]] = coda4microbiome::coda_glmnet(Microbiome.data.unaff[,rownames(sPLS.sCCA.Konzo.output.unaff$sPLS$loadings$X[rowSums(sPLS.sCCA.Konzo.output.unaff$sPLS$loadings$X)!=0,])], log(Metabolome.data.unaff[,m]))
  
}
pty = c()
for(p in names(l.unaff[names(l.unaff)%in%names(l.aff)])){
  pty = c(pty,length(l.unaff[[p]]$taxa.name[l.unaff[[p]]$taxa.name %in% l.aff[[p]]$taxa.name])/length(l.unaff[[p]]$taxa.name))
}

pty[is.na(pty)]= 0
summary(pty)

sapply(l.unaff[names(l.unaff)%in%names(l.aff)], function(x) length(x$taxa.name))

boxplot(sapply(l.aff[names(l.aff)%in%names(l.unaff)], function(x) length(x$taxa.name)), 
            sapply(l.aff[!names(l.aff)%in%names(l.unaff)], function(x) length(x$taxa.name)))

names(l.aff[!names(l.aff)%in%names(l.unaff)])
names(l.unaff[!names(l.unaff)%in%names(l.aff)])

l.aff$X30$taxa.name[l.aff$X30$taxa.name%in%l.aff$X111$taxa.name]

l.aff$X30$`signature plot`
l.aff$X111$`signature plot`

df.CODA.example.A = data.frame("Coefficient" = l.aff$X30$`log-contrast coefficients`,"Species"=l.aff$X30$taxa.name)
df.CODA.example.A$Species = gsub("\\."," ", df.CODA.example.A$Species)
df.CODA.example.A$sign = ifelse(df.CODA.example.A$Coefficient<0, "-", "+")

a=ggplot(df.CODA.example.A, aes(x = reorder(Species,-abs(Coefficient)), y = abs(Coefficient), fill=sign, label=sign)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Coefficients")+xlab("Microorganism")+geom_label(size=3.5)+ylim(c(0,1))

df.CODA.example.B = data.frame("Coefficient" = l.aff$X111$`log-contrast coefficients`,"Species"=l.aff$X111$taxa.name)
df.CODA.example.B$Species = gsub("\\."," ", df.CODA.example.B$Species)
df.CODA.example.B$sign = ifelse(df.CODA.example.B$Coefficient<0, "-", "+")

b=ggplot(df.CODA.example.B, aes(x = reorder(Species,-abs(Coefficient)), y = abs(Coefficient), fill=sign, label=sign)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Coefficients")+xlab("Microorganism")+geom_label(size=3.5)+ylim(c(0,1))

Figure6D = (a + ggtitle("mevalonate")+theme(axis.title.y=element_text(vjust=5)))/(b+ggtitle("3-hydroxyisobutyrate")+theme(axis.title.y=element_text(vjust=5)))


gr = igraph::graph_from_data_frame(rbind(data.frame("from"="mevalonate", "to"=l.aff$X30$taxa.name, "weight"=l.aff$X30$`log-contrast coefficients`),
                                         data.frame("from"="3-hydroxyisobutyrate", "to"=l.aff$X111$taxa.name, "weight"=l.aff$X111$`log-contrast coefficients`),
                                         data.frame("from"="1-palmitoyl-2-oleoyl-GPE (16:0/18:1)", "to"=l.aff$X1526$taxa.name, "weight"=l.aff$X1526$`log-contrast coefficients`),
                                         
                                         data.frame("from"="2-hydroxybutyrate/2-hydroxyisobutyrate", "to"=l.aff$X100008928$taxa.name, "weight"=l.aff$X100008928$`log-contrast coefficients`),
                                         data.frame("from"="3-hydroxystearate", "to"=l.aff$X100008933$taxa.name, "weight"=l.aff$X100008933$`log-contrast coefficients`),
                                         data.frame("from"="1,2-dipalmitoyl-GPE (16:0/16:0)", "to"=l.aff$X100009204$taxa.name, "weight"=l.aff$X100009204$`log-contrast coefficients`),
                                         data.frame("from"="2-hydroxynervonate", "to"=l.aff$X100015638$taxa.name, "weight"=l.aff$X100015638$`log-contrast coefficients`),
                                         data.frame("from"="1-deoxysphinganine (m18:0)", "to"=l.aff$X100021873$taxa.name, "weight"=l.aff$X100021873$`log-contrast coefficients`)), directed = FALSE)

V(gr)$Type = ifelse(V(gr)$name %in% c("mevalonate", "3-hydroxyisobutyrate", "1-palmitoyl-2-oleoyl-GPE (16:0/18:1)",
                                       "2-hydroxybutyrate/2-hydroxyisobutyrate", "3-hydroxystearate", "1,2-dipalmitoyl-GPE (16:0/16:0)",
                                       "2-hydroxynervonate", "1-deoxysphinganine (m18:0)"), "Microorganisms", "Metabolites")
V(gr)$col = ifelse(V(gr)$Type =="Microorganisms", "brown", "pink")
V(gr)$size = ifelse(V(gr)$name %in% c("mevalonate", "3-hydroxyisobutyrate"), 6, 1)
V(gr)$label.cex = ifelse(V(gr)$name %in% c("mevalonate", "3-hydroxyisobutyrate"), 2, 0.5)
V(gr)$type = bipartite.mapping(gr)$type
E(gr)$color = ifelse(E(gr)$weight>0, "lightgreen", "pink")
E(gr)$weight = abs(E(gr)$weight)*2

LO = layout_as_bipartite(gr)
LO = LO[,c(2,1)]

Figure6E = ggnet2(gr, color = "Type", edge.size = "weight", edge.color="color", size=4, label=ifelse(V(gr)$name%in%c("mevalonate", "3-hydroxyisobutyrate"), V(gr)$name, NA), palette = c("Metabolites"="#BC3C29FF", "Microorganisms"="#0072B5FF"))

load("RDA.rda")
load("RDA_inverted.rda")


df.explained.RDA.metabolome.affected = summary(eigenvals(spe.rda.affected.inverted)) %>% t %>%  as.data.frame() %>% tibble::rownames_to_column("RDA") %>%   dplyr::select("Cumulative Proportion") %>%   mutate(index = 1:nrow(.))
df.explained.RDA.metabolome.affected$Group = "Affected"
df.explained.RDA.metabolome.unaffected = summary(eigenvals(spe.rda.unnaffected.inverted)) %>% t %>%  as.data.frame() %>% tibble::rownames_to_column("RDA") %>%   dplyr::select("Cumulative Proportion") %>%   mutate(index = 1:nrow(.))
df.explained.RDA.metabolome.unaffected$Group = "Unaffected"

df.explained.RDA.metabolome = rbind(df.explained.RDA.metabolome.affected, df.explained.RDA.metabolome.unaffected)
df.explained.RDA.metabolome$Omics = "Metabolome"

df.explained.RDA.microbiome.affected = summary(eigenvals(spe.rda.affected)) %>% t %>%  as.data.frame() %>% tibble::rownames_to_column("RDA") %>%   dplyr::select("Cumulative Proportion") %>%   mutate(index = 1:nrow(.))
df.explained.RDA.microbiome.unaffected = summary(eigenvals(spe.rda.unnaffected)) %>% t %>%  as.data.frame() %>% tibble::rownames_to_column("RDA") %>%   dplyr::select("Cumulative Proportion") %>%   mutate(index = 1:nrow(.))
df.explained.RDA.microbiome.affected$Group = "Affected"
df.explained.RDA.microbiome.unaffected$Group = "Unaffected"

df.explained.RDA.microbiome = rbind(df.explained.RDA.microbiome.unaffected, df.explained.RDA.microbiome.affected)
df.explained.RDA.microbiome$Omics = "Microbiome"

df.explained.RDS.all = rbind(df.explained.RDA.metabolome,df.explained.RDA.microbiome)

Figure6A = ggplot(df.explained.RDS.all, aes(x=index,y=`Cumulative Proportion`, color=Group))+geom_point()+geom_line()+facet_grid(.~Omics)+xlab("RDA Components")+ylab("% Cumulative Explained Variance") + theme_bw()

Chem_info = read.csv2("Chemical_annot_Metabo.csv", header=T, sep=",")
Chem_info$CHEM_ID = paste0("X",Chem_info$CHEM_ID)
library(tidyr)

df.contrib.unaffected.species = data.frame(scores(spe.rda.unnaffected)$species)
df.contrib.unaffected.species$species = rownames(df.contrib.unaffected.species)
df.contrib.unaffected.species$sign_RDA1 = ifelse(df.contrib.unaffected.species$RDA1<0,"-","+")
df.contrib.unaffected.species$sign_RDA2 = ifelse(df.contrib.unaffected.species$RDA2<0,"-","+")
df.contrib.unaffected.species$abs_RDA1 = abs(df.contrib.unaffected.species$RDA1)
df.contrib.unaffected.species$abs_RDA2 = abs(df.contrib.unaffected.species$RDA2)

df.contrib.unaffected.species.top20.RDA1 = df.contrib.unaffected.species[order(df.contrib.unaffected.species$abs_RDA1, decreasing = T)[1:20],]

Figure6AA = ggplot(df.contrib.unaffected.species.top20.RDA1, aes(x = reorder(species,-abs(abs_RDA1)), y = abs(abs_RDA1), fill=sign_RDA1, label=sign_RDA1)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Species")+geom_label(size=3.5)+ylim(c(0,1))


df.contrib.affected.species = data.frame(scores(spe.rda.affected)$species)
df.contrib.affected.species$species = rownames(df.contrib.affected.species)
df.contrib.affected.species$sign_RDA1 = ifelse(df.contrib.affected.species$RDA1<0,"-","+")
df.contrib.affected.species$sign_RDA2 = ifelse(df.contrib.affected.species$RDA2<0,"-","+")
df.contrib.affected.species$abs_RDA1 = abs(df.contrib.affected.species$RDA1)
df.contrib.affected.species$abs_RDA2 = abs(df.contrib.affected.species$RDA2)

df.contrib.affected.species.top20.RDA1 = df.contrib.affected.species[order(df.contrib.affected.species$abs_RDA1, decreasing = T)[1:20],]

Figure6AB = ggplot(df.contrib.affected.species.top20.RDA1, aes(x = reorder(species,-abs(abs_RDA1)), y = abs(abs_RDA1), fill=sign_RDA1, label=sign_RDA1)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Species")+geom_label(size=3.5)+ylim(c(0,1))


df.contrib.unaffected.metabolites = data.frame(scores(spe.rda.unnaffected.inverted)$species)
df.contrib.unaffected.metabolites$CHEM_ID = rownames(df.contrib.unaffected.metabolites)
df.contrib.unaffected.metabolites = merge(df.contrib.unaffected.metabolites, Chem_info[,c("CHEM_ID","CHEMICAL_NAME")], by="CHEM_ID")
df.contrib.unaffected.metabolites$sign_RDA1 = ifelse(df.contrib.unaffected.metabolites$RDA1<0,"-","+")
df.contrib.unaffected.metabolites$sign_RDA2 = ifelse(df.contrib.unaffected.metabolites$RDA2<0,"-","+")
df.contrib.unaffected.metabolites$abs_RDA1 = abs(df.contrib.unaffected.metabolites$RDA1)
df.contrib.unaffected.metabolites$abs_RDA2 = abs(df.contrib.unaffected.metabolites$RDA2)

df.contrib.unaffected.metabolites.top20.RDA1 = df.contrib.unaffected.metabolites[order(df.contrib.unaffected.metabolites$abs_RDA1, decreasing = T)[1:20],]

Figure6BA = ggplot(df.contrib.unaffected.metabolites.top20.RDA1, aes(x = reorder(CHEMICAL_NAME,-abs(abs_RDA1)), y = abs(abs_RDA1), fill=sign_RDA1, label=sign_RDA1)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Metabolites")+geom_label(size=3.5)+ylim(c(0,1))


df.contrib.affected.metabolites = data.frame(scores(spe.rda.affected.inverted)$species)
df.contrib.affected.metabolites$CHEM_ID = rownames(df.contrib.affected.metabolites)
df.contrib.affected.metabolites = merge(df.contrib.affected.metabolites, Chem_info[,c("CHEM_ID","CHEMICAL_NAME")], by="CHEM_ID")
df.contrib.affected.metabolites$sign_RDA1 = ifelse(df.contrib.affected.metabolites$RDA1<0,"-","+")
df.contrib.affected.metabolites$sign_RDA2 = ifelse(df.contrib.affected.metabolites$RDA2<0,"-","+")
df.contrib.affected.metabolites$abs_RDA1 = abs(df.contrib.affected.metabolites$RDA1)
df.contrib.affected.metabolites$abs_RDA2 = abs(df.contrib.affected.metabolites$RDA2)

df.contrib.affected.metabolites.top20.RDA1 = df.contrib.affected.metabolites[order(df.contrib.affected.metabolites$abs_RDA1, decreasing = T)[1:20],]

Figure6BB = ggplot(df.contrib.affected.metabolites.top20.RDA1, aes(x = reorder(CHEMICAL_NAME,-abs(abs_RDA1)), y = abs(abs_RDA1), fill=sign_RDA1, label=sign_RDA1)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Metabolites")+geom_label(size=3.5)+ylim(c(0,1))

Figure6B = ((Figure6AB+ggtitle("Affected"))+Figure6BB ) /((Figure6AA+ggtitle("Unaffected"))+Figure6BA )


library(ggpubr)

ggarrange(ggarrange(Figure6A, Figure6CA, Figure6D, nrow=3,labels=c("A.","C.","E.")), ggarrange(Figure6B, Figure6CB, Figure6E, labels=c("B.", "D","F." ), nrow=3), ncol=2)



library(MiRKAT)


D1= as.matrix(dist(ilr_microbiome_affected, "euclidean"))
D2=as.matrix(dist(ilr_microbiome_affected, "manhattan"))
D3= as.matrix(dist(ilr_microbiome_affected, "canberra"))

D1[is.na(D1)]=0
D2[is.na(D2)]=0
D3[is.na(D3)]=0


l = list("E"=D2K(D1),"M"=D2K(D2),"C"=D2K(D3))

qq = rownames(sPLS.sCCA.Konzo.output.aff$sPLS$loadings$Y[rowSums(sPLS.sCCA.Konzo.output.aff$sPLS$loadings$Y)!=0,])


for(q in qq) {
  print(q)
  print(MiRKAT(y=log_metabolome_affected[,q],Ks=l, omnibus="cauchy")$omnibus_p <= 0.05/30)
}



df.contrib.unaffected.species.top20.RDA2 = df.contrib.unaffected.species[order(df.contrib.unaffected.species$abs_RDA2, decreasing = T)[1:20],]

FigureS10A = ggplot(df.contrib.unaffected.species.top20.RDA2, aes(x = reorder(species,-abs(abs_RDA2)), y = abs(abs_RDA2), fill=sign_RDA2, label=sign_RDA2)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Species")+geom_label(size=3.5)+ylim(c(0,1))


df.contrib.affected.species.top20.RDA2 = df.contrib.affected.species[order(df.contrib.affected.species$abs_RDA2, decreasing = T)[1:20],]

FigureS10B = ggplot(df.contrib.affected.species.top20.RDA2, aes(x = reorder(species,-abs(abs_RDA2)), y = abs(abs_RDA2), fill=sign_RDA2, label=sign_RDA2)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Species")+geom_label(size=3.5)+ylim(c(0,1))


df.contrib.unaffected.metabolites.top20.RDA2 = df.contrib.unaffected.metabolites[order(df.contrib.unaffected.metabolites$abs_RDA2, decreasing = T)[1:20],]

FigureS10C = ggplot(df.contrib.unaffected.metabolites.top20.RDA2, aes(x = reorder(CHEMICAL_NAME,-abs(abs_RDA2)), y = abs(abs_RDA2), fill=sign_RDA2, label=sign_RDA2)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Metabolites")+geom_label(size=3.5)+ylim(c(0,1))

df.contrib.affected.metabolites.top20.RDA2 = df.contrib.affected.metabolites[order(df.contrib.affected.metabolites$abs_RDA2, decreasing = T)[1:20],]

FigureS10D = ggplot(df.contrib.affected.metabolites.top20.RDA2, aes(x = reorder(CHEMICAL_NAME,-abs(abs_RDA2)), y = abs(abs_RDA2), fill=sign_RDA2, label=sign_RDA2)) + geom_bar(stat="identity", position="dodge", color="black")+coord_flip()+theme_bw()+theme(legend.position = "none")+
  ylab("Weights")+xlab("Metabolites")+geom_label(size=3.5)+ylim(c(0,1))

FigureS10 = ((FigureS10C+ggtitle("Unaffected"))+(FigureS10D+xlab("")+ggtitle("Affected") )) /((FigureS10A)+(FigureS10B+xlab("")))

sapply(l.unaff[!names(l.unaff)%in%names(l.aff)], function(x) x$`log-contrast coefficients`)


gr2 = igraph::graph_from_data_frame(do.call("rbind",lapply(names(l.unaff[!names(l.unaff)%in%names(l.aff)]), function(x){
  if(length(l.unaff[[x]]$taxa.name)>0) data.frame("from"=x, "to"=l.unaff[[x]]$taxa.name, "weight"= l.unaff[[x]]$`log-contrast coefficients`)
  
})) , directed = FALSE)

V(gr2)$Type = ifelse(V(gr2)$name %in% names(l.unaff), "Metabolites", "Microorganisms")
V(gr2)$col = ifelse(V(gr2)$Type =="Microorganisms", "brown", "pink")
V(gr2)$type = bipartite.mapping(gr2)$type
E(gr2)$color = ifelse(E(gr2)$weight>0, "lightgreen", "pink")
E(gr2)$weight = abs(E(gr2)$weight)*2

LO = layout_as_bipartite(gr2)
LO = LO[,c(2,1)]

ggnet2(gr2, color = "Type", edge.size = "weight", edge.color="color", size=4, label=NA, palette = c("Metabolites"="#BC3C29FF", "Microorganisms"="#0072B5FF"))
