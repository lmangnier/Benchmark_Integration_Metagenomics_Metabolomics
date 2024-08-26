#This code produces panel 3E and Figure 4 presented in the paper https://doi.org/10.1101/2024.01.26.577441
#Data can be found there: 10.6084/m9.figshare.25234915
#Also supplementary figures are produced related to individual associations and data summarization

library(ggh4x)
library(dplyr)
library(stringr)
library(patchwork)

data_HALLA_Autism = list.files("results_HALLA_Autism", full.names = T)
data_HALLA_Adenomas = list.files("results_HALLA_Adenomas", full.names = T)
data_HALLA_Konzo = list.files("results_HALLA_Konzo", full.names = T)

#Null: Main scenario clr-log
HALLA_pvalues_null_clr_log_Autism = read.table(data_HALLA_Autism[2], header=TRUE)
HALLA_pvalues_null_clr_original_Adenomas = read.table(data_HALLA_Adenomas[2], header=TRUE)
HALLA_pvalues_null_clr_log_Konzo = read.table(data_HALLA_Konzo[2], header=TRUE)

#Alternative: Main scenario clr-log
HALLA_pvalues_alter_clr_log_Autism = read.table(data_HALLA_Autism[1], header=TRUE)
HALLA_pvalues_alter_clr_original_Adenomas = read.table(data_HALLA_Adenomas[1], header=TRUE)
HALLA_pvalues_alter_clr_log_Konzo = read.table(data_HALLA_Konzo[1], header=TRUE)

#Correlation
null_Cor_autism_spearman = readRDS("ACAT_combined_cor_null_Autism.RDS")
null_Cor_adenomas_spearman = readRDS("ACAT_combined_cor_null_Adenomas.RDS")
null_Cor_konzo_spearman = readRDS("ACAT_combined_null_cor_Konzo.RDS")

null_Cor_autism_pearson = readRDS("ACAT_cor_pearson_null_Autism")
null_Cor_adenomas_pearson = readRDS("ACAT_combined_cor_null_Adenomas.RDS")
null_Cor_konzo_pearson = readRDS("ACAT_cor_pearson_null_Konzo")


df_qqplot_spearman = data.frame("Dataset" = c(rep("Adenomas", 120000), rep("Autism", 22000), rep("Konzo", 85000)), 
                                "observed"= -log10(c(sort(c(null_Cor_adenomas_spearman)),sort(c(null_Cor_autism_spearman)),sort(c(null_Cor_konzo_spearman)))), 
                                  "expected"= c(-log10(ppoints(120000)), -log10(ppoints(22000)), -log10(ppoints(85000))),
                                  "clower"   = c(-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:120000, shape2 = 12000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))),
                                  "cupper"   = c(-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:120000, shape2 = 12000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))))

df_qqplot_pearson = data.frame("Dataset" = c(rep("Adenomas", 120000), rep("Autism", 22000), rep("Konzo", 85000)), 
                                "observed"= -log10(c(sort(c(null_Cor_adenomas_spearman)),sort(c(null_Cor_autism_spearman)),sort(c(null_Cor_konzo_spearman)))), 
                                "expected"= c(-log10(ppoints(120000)), -log10(ppoints(22000)), -log10(ppoints(85000))),
                                "clower"   = c(-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:120000, shape2 = 12000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))),
                                "cupper"   = c(-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:120000, shape2 = 12000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))))

ggplot(df_qqplot_spearman) +
  geom_ribbon(
    mapping = aes(x = expected, ymin = clower, ymax = cupper),
    alpha = 0.1
  ) +
  geom_point(aes(expected, observed, col = Dataset), size = 2) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  xlab(expression(paste("Expected -log"[10], plain(P)))) +
  ylab(expression(paste("Observed -log"[10], plain(P)))) +
  theme_bw()#+theme(legend.position = "none")

ggplot(df_qqplot_pearson) +
  geom_ribbon(
    mapping = aes(x = expected, ymin = clower, ymax = cupper),
    alpha = 0.1
  ) +
  geom_point(aes(expected, observed, col = Dataset), size = 2) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  xlab(expression(paste("Expected -log"[10], plain(P)))) +
  ylab(expression(paste("Observed -log"[10], plain(P)))) +
  theme_bw()

#QQplot Cor
pvalues_null_Cor = c(c(sapply(1:1000, function(x) sapply(null_Cor_ade[[x]], function(y) ACAT::ACAT(sapply( y, function(z) z$spearman))))),
                    c(sapply(1:1000, function(x) sapply(null_Cor_autism[[x]], function(y) ACAT::ACAT(sapply( y, function(z) z$spearman))))),
                     c(sapply(1:1000, function(x) sapply(null_Cor_autism[[x]], function(y) ACAT::ACAT(sapply( y, function(z) z$spearman))))))

#Power Cor
alter_Cor_autism = readRDS("Cor_results\\Autism\\results_clr_cor_null_subset_Autism.RDS")
alter_Cor_adenomas = readRDS("Cor_results\\Adenomas\\results_clr_cor_null_subset_Adenomas.RDS")
alter_Cor_konzo = readRDS("Cor_results\\Konzo\\results_clr_cor_null_subset_Konzo.RDS")




univariate_autism = readRDS("list_results_univariate_Autism.RDS")
univariate_konzo = readRDS("list_results_univariate_Konzo.RDS")
univariate_adenomas = readRDS("list_results_univariate_Adenomas.RDS")

MiRKAT_autism = readRDS("results_MiRKAT_Autism.RDS")
MiRKAT_konzo = readRDS("results_MiRKAT_Konzo.RDS")
MiRKAT_adenomas = readRDS("results_MiRKAT_Adenomas.RDS")



pvalues_univariate_null= c(sort(c(univariate_adenomas$`Log-Contrast`$null)),
sort(c(univariate_autism$`Log-contrast`$null)),
sort(c(univariate_konzo$`Log-Contrast`$null)),
sort(c(sapply(MiRKAT_adenomas$null$ILR, function(x) sapply(x, function(y) y$omnibus_p)))),
sort(c(sapply(MiRKAT_autism$null$ILR, function(x) sapply(x, function(y) y$omnibus_p)))),
sort(c(sapply(MiRKAT_konzo$null$ILR, function(x) sapply(x, function(y) y$omnibus_p)))),
sort(c(univariate_adenomas$`clr-lm`$null)),
sort(c(univariate_autism$`clr-lm`$null)),
sort(c(univariate_konzo$`clr-lm`$null)),
sort(HALLA_pvalues_null_clr_original_Adenomas$p.values),
sort(HALLA_pvalues_null_clr_log_Autism$p.values),
sort(HALLA_pvalues_null_clr_log_Konzo$p.values))


df_qqplot_univariate = data.frame("Dataset" = rep(c(rep("Adenomas", 120000), rep("Autism", 22000), rep("Konzo", 85000)),4), "Method" = rep(c("Log-contrast", "MiRKAT", "clr-lm", "HALLA"), each=227000 ),"observed"= -log10(pvalues_univariate_null), 
                                  "expected"= rep(c(-log10(ppoints(120000)), -log10(ppoints(22000)), -log10(ppoints(85000))),4),
                                  "clower"   = rep(c(-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:120000, shape2 = 120000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))),4),
                                  "cupper"   = rep(c(-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:120000, shape2 = 120000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))),4))

a = ggplot(df_qqplot_univariate) +
  geom_ribbon(
    mapping = aes(x = expected, ymin = clower, ymax = cupper),
    alpha = 0.1
  ) +
  geom_point(aes(expected, observed, col = Method), size = 2) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  xlab(expression(paste("Expected -log"[10], plain(P)))) +
  ylab(expression(paste("Observed -log"[10], plain(P)))) +
  theme_bw() +#+theme(legend.position = "none")
  facet_grid(.~Dataset) +
  ggtitle("A.")

pvalues_univariate_alter= c(mean(c(univariate_adenomas$`Log-Contrast`$alter)<=0.05),
                           mean(c(univariate_autism$`Log-contrast`$alter)<=0.05),
                           mean(unlist(univariate_konzo$`Log-Contrast`$alter)<=0.05),
                           mean(c(sapply(MiRKAT_adenomas$alter$ILR, function(x) sapply(x, function(y) y$omnibus_p)))<=0.05),
                           mean(c(sapply(MiRKAT_autism$alter$ILR, function(x) sapply(x, function(y) y$omnibus_p)))<=0.05),
                           mean(unlist(sapply(MiRKAT_konzo$alter$ILR, function(x) sapply(x, function(y) y$omnibus_p)))<=0.05),
                           mean(c(univariate_adenomas$`clr-lm`$alter)<=0.05),
                           mean(c(univariate_autism$`clr-lm`$alter)<=0.05),
                           mean(unlist(univariate_konzo$`clr-lm`$alter)<=0.05),
                           mean(HALLA_pvalues_alter_clr_original_Adenomas$p.values<=0.05),
                           mean(HALLA_pvalues_alter_clr_log_Autism$p.values<=0.05),
                           mean(HALLA_pvalues_alter_clr_log_Konzo$p.values<=0.05))




df_power_univariate = data.frame("Dataset" = rep(c("Adenomas","Autism","Konzo"), 4),"Method" = rep(c("Log-contrast", "MiRKAT", "clr-lm", "HALLA"), each=3),"power"= pvalues_univariate_alter)

df_mean_power_spearman = data.frame("Dataset"=c("Adenomas","Autism","Konzo"), "Power"=c(0.1048,0.058,0.108))


b= ggplot(df_power_univariate, aes(x=Method, y=power, fill=Method))+geom_bar(colour="black",stat="identity", position="dodge")+geom_hline(data=df_mean_power_spearman, aes(yintercept=Power,linetype="awesome"))+
facet_grid(.~Dataset)+ylim(c(0,1))+ylab("Power")+theme_bw()+theme(legend.position = "none")+ggtitle("B.")

a + b

#Multivariate
#CCA
cca.ilr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_CCA_ilr.RDS")
cca.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_CCA_ilr_log.RDS")
cca.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_CCA_ilr_log.RDS")

cca.clr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_CCA_clr.RDS")
cca.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_CCA_clr_log.RDS")
cca.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_CCA_clr_log.RDS")

cca.alpha.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_CCA_alpha.RDS")
cca.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_CCA_alpha_log.RDS")
cca.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_CCA_alpha_log.RDS")

#RDA
RDA.ilr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\RDA1_ilr_Adenomas.RDS")
RDA.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA1_ilr_log_Konzo.RDS")
RDA.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA1_ilr_log_Autism.RDS")

RDA.clr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\RDA1_clr_Adenomas.RDS")
RDA.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA1_clr_log_Konzo.RDS")
RDA.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA1_clr_log_Autism.RDS")

RDA.alpha.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\RDA1_alpha_Adenomas.RDS")
RDA.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA1_alpha_log_Konzo.RDS")
RDA.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA1_alpha_log_Autism.RDS")

RDA2.ilr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\RDA2_ilr_Adenomas.RDS")
RDA2.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA2_ilr_log_Konzo.RDS")
RDA2.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA2_ilr_log_Autism.RDS")

RDA2.clr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\RDA2_clr_Adenomas.RDS")
RDA2.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA2_clr_log_Konzo.RDS")
RDA2.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA2_clr_log_Autism.RDS")

RDA2.alpha.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\RDA2_alpha_Adenomas.RDS")
RDA2.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA2_alpha_log_Konzo.RDS")
RDA2.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA2_alpha_log_Autism.RDS")


#MOFA2

MOFA2.ilr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\MOFA2_ilr_Adenomas.RDS")
MOFA2.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\MOFA2_ilr_log_Konzo.RDS")
MOFA2.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\MOFA2_ilr_log_Autism.RDS")

MOFA2.clr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\MOFA2_clr_Adenomas.RDS")
MOFA2.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\MOFA2_clr_log_Konzo.RDS")
MOFA2.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\MOFA2_clr_log_Autism.RDS")

MOFA2.alpha.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\MOFA2_alpha_Adenomas.RDS")
MOFA2.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\MOFA2_alpha_log_Konzo.RDS")
MOFA2.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\MOFA2_alpha_log_Autism.RDS")

#PLS-Can
PLS.Can.ilr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_PLS_canonical_ilr.RDS")
PLS.Can.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_canonical_ilr_log.RDS")
PLS.Can.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_canonical_ilr_log.RDS")

PLS.Can.clr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_PLS_canonical_clr.RDS")
PLS.Can.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_canonical_clr_log.RDS")
PLS.Can.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_canonical_clr_log.RDS")

PLS.Can.alpha.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_PLS_canonical_alpha.RDS")
PLS.Can.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_canonical_alpha_log.RDS")
PLS.Can.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_canonical_alpha_log.RDS")


#PLS-Reg
PLS.Reg1.ilr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_PLS_regression1_ilr.RDS")
PLS.Reg1.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression1_ilr_log.RDS")
PLS.Reg1.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression1_ilr_log.RDS")

PLS.Reg1.clr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_PLS_regression1_clr.RDS")
PLS.Reg1.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression1_clr.RDS")
PLS.Reg1.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression1_clr.RDS")

PLS.Reg1.alpha.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_PLS_regression1_alpha.RDS")
PLS.Reg1.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression1_alpha.RDS")
PLS.Reg1.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression1_alpha.RDS")

PLS.Reg2.ilr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_PLS_regression2_ilr.RDS")
PLS.Reg2.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression2_ilr.RDS")
PLS.Reg2.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression2_ilr.RDS")

PLS.Reg2.clr.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_PLS_regression2_clr.RDS")
PLS.Reg2.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression2_clr.RDS")
PLS.Reg2.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression2_clr.RDS")

PLS.Reg2.alpha.Adenomas = readRDS("results_multivariate_Adenomas\\results_multivariate_Adenomas\\results_PLS_regression2_alpha.RDS")
PLS.Reg2.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression2_alpha.RDS")
PLS.Reg2.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression2_alpha.RDS")


explained.variance.ilr.PLS.Can1.Adenomas = sapply(PLS.Can.ilr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Can1.Konzo = sapply(PLS.Can.ilr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Can1.Autism = sapply(PLS.Can.ilr.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$Y)
  else NA
})

explained.variance.ilr.PLS.Can2.Adenomas = sapply(PLS.Can.ilr.Adenomas, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.ilr.PLS.Can2.Konzo = sapply(PLS.Can.ilr.Konzo, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.ilr.PLS.Can2.Autism = sapply(PLS.Can.ilr.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$X)
  else NA
})

explained.variance.ilr.PLS.Reg1.Adenomas = sapply(PLS.Reg1.ilr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Reg1.Konzo = sapply(PLS.Reg1.ilr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Reg1.Autism = sapply(PLS.Reg1.ilr.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.ilr.PLS.Reg2.Adenomas = sapply(PLS.Reg2.ilr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Reg2.Konzo = sapply(PLS.Reg2.ilr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Reg2.Autism = sapply(PLS.Reg2.ilr.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.ilr.CCA1.Adenomas = sapply(cca.ilr.Adenomas, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))
explained.variance.ilr.CCA1.Konzo = sapply(cca.ilr.Konzo, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))
explained.variance.ilr.CCA1.Autism = sapply(cca.ilr.Autism, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))

explained.variance.ilr.CCA2.Adenomas = sapply(cca.ilr.Adenomas, function(x) sum(x$Redundancy$ConditionalRedundancy$X))
explained.variance.ilr.CCA2.Konzo = sapply(cca.ilr.Konzo, function(x) sum(x$Redundancy$ConditionalRedundancy$X))
explained.variance.ilr.CCA2.Autism = sapply(cca.ilr.Autism, function(x) sum(x$Redundancy$ConditionalRedundancy$X))


explained.variance.ilr.RDA1.Adenomas = unlist(RDA.ilr.Adenomas)
explained.variance.ilr.RDA1.Konzo = unlist(RDA.ilr.Konzo)
explained.variance.ilr.RDA1.Autism = unlist(RDA.ilr.Autism)

explained.variance.ilr.RDA2.Adenomas = unlist(RDA2.ilr.Adenomas)
explained.variance.ilr.RDA2.Konzo = unlist(RDA2.ilr.Konzo)
explained.variance.ilr.RDA2.Autism = unlist(RDA2.ilr.Autism)


explained.variance.ilr.MOFA.Metabo.Adenomas = sapply(MOFA2.ilr.Adenomas, function(x) x$r2_total$group1[1])
explained.variance.ilr.MOFA.Metabo.Konzo = sapply(MOFA2.ilr.Konzo, function(x) x$r2_total$group1[1])
explained.variance.ilr.MOFA.Metabo.Autism = sapply(MOFA2.ilr.Autism, function(x) x$r2_total$group1[1])

explained.variance.ilr.MOFA.Metabo.Adenomas[explained.variance.ilr.MOFA.Metabo.Adenomas<0] = 0
explained.variance.ilr.MOFA.Metabo.Konzo[explained.variance.ilr.MOFA.Metabo.Konzo<0] = 0
explained.variance.ilr.MOFA.Metabo.Autism[explained.variance.ilr.MOFA.Metabo.Autism<0] = 0

explained.variance.ilr.MOFA.Micro.Adenomas = sapply(MOFA2.ilr.Adenomas, function(x) x$r2_total$group1[2])
explained.variance.ilr.MOFA.Micro.Konzo = sapply(MOFA2.ilr.Konzo, function(x) x$r2_total$group1[2])
explained.variance.ilr.MOFA.Micro.Autism = sapply(MOFA2.ilr.Autism, function(x) x$r2_total$group1[2])

explained.variance.ilr.MOFA.Micro.Adenomas[explained.variance.ilr.MOFA.Micro.Adenomas<0] = 0
explained.variance.ilr.MOFA.Micro.Konzo[explained.variance.ilr.MOFA.Micro.Konzo<0] = 0
explained.variance.ilr.MOFA.Micro.Autism[explained.variance.ilr.MOFA.Micro.Autism<0] = 0

summary(explained.variance.ilr.MOFA.Micro.Konzo)
df.explained.variance.ilr = data.frame("Normalization"="ilr","Method"=rep(c("PLS-Can1", "PLS-Can2", "PLS-Reg1", "PLS-Reg2", "CCA1", "CCA2", "RDA1","RDA2", "MOFA2_1", "MOFA2_2"), each="3000"), "Dataset"=rep(rep(c("Adenomas", "Konzo", "Autism"), each=1000),10),
                                       "Explained_Variance"=c(explained.variance.ilr.PLS.Can1.Adenomas,explained.variance.ilr.PLS.Can1.Konzo,explained.variance.ilr.PLS.Can1.Autism,
                                                              explained.variance.ilr.PLS.Can2.Adenomas,explained.variance.ilr.PLS.Can2.Konzo,explained.variance.ilr.PLS.Can2.Autism,
                                                              explained.variance.ilr.PLS.Reg1.Adenomas,explained.variance.ilr.PLS.Reg1.Konzo,explained.variance.ilr.PLS.Reg1.Autism,
                                                              explained.variance.ilr.PLS.Reg2.Adenomas,explained.variance.ilr.PLS.Reg2.Konzo,explained.variance.ilr.PLS.Reg2.Autism,
                                                              explained.variance.ilr.CCA1.Adenomas,explained.variance.ilr.CCA1.Konzo,explained.variance.ilr.CCA1.Autism,
                                                              explained.variance.ilr.CCA2.Adenomas,explained.variance.ilr.CCA2.Konzo,explained.variance.ilr.CCA2.Autism,
                                                              explained.variance.ilr.RDA1.Adenomas,explained.variance.ilr.RDA1.Konzo,explained.variance.ilr.RDA1.Autism,
                                                              explained.variance.ilr.RDA2.Adenomas,explained.variance.ilr.RDA2.Konzo,explained.variance.ilr.RDA2.Autism,
                                                              explained.variance.ilr.MOFA.Metabo.Adenomas/100,explained.variance.ilr.MOFA.Metabo.Konzo/100,explained.variance.ilr.MOFA.Metabo.Autism/100,
                                                              explained.variance.ilr.MOFA.Micro.Adenomas/100,explained.variance.ilr.MOFA.Micro.Konzo/100,explained.variance.ilr.MOFA.Micro.Autism/100))
library(ggplot2)
ggplot(df.explained.variance.ilr, aes(x=Method,y=Explained_Variance))+geom_boxplot()+facet_grid(.~Dataset)


explained.variance.clr.PLS.Can1.Adenomas = sapply(PLS.Can.clr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Can1.Konzo = sapply(PLS.Can.clr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Can1.Autism = sapply(PLS.Can.clr.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$Y)
  else NA
})

explained.variance.clr.PLS.Can2.Adenomas = sapply(PLS.Can.clr.Adenomas, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.clr.PLS.Can2.Konzo = sapply(PLS.Can.clr.Konzo, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.clr.PLS.Can2.Autism = sapply(PLS.Can.clr.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$X)
  else NA
})

explained.variance.clr.PLS.Reg1.Adenomas = sapply(PLS.Reg1.clr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Reg1.Konzo = sapply(PLS.Reg1.clr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Reg1.Autism = sapply(PLS.Reg1.clr.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.clr.PLS.Reg2.Adenomas = sapply(PLS.Reg2.clr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Reg2.Konzo = sapply(PLS.Reg2.clr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Reg2.Autism = sapply(PLS.Reg2.clr.Autism, function(x) sum(x$ConditionalRedundancy$Y))


explained.variance.clr.CCA1.Adenomas = sapply(cca.clr.Adenomas, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$Y)
  else NA
})
explained.variance.clr.CCA1.Konzo = sapply(cca.clr.Konzo, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$Y)
  else NA
})
explained.variance.clr.CCA1.Autism = sapply(cca.clr.Autism, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$Y)
  else NA
})

explained.variance.clr.CCA2.Adenomas = sapply(cca.clr.Adenomas, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$X)
  else NA
})
explained.variance.clr.CCA2.Konzo = sapply(cca.clr.Konzo, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$X)
  else NA
})
explained.variance.clr.CCA2.Autism = sapply(cca.clr.Autism, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$X)
  else NA
})

explained.variance.clr.RDA1.Adenomas = unlist(RDA.clr.Adenomas)
explained.variance.clr.RDA1.Konzo = unlist(RDA.clr.Konzo)
explained.variance.clr.RDA1.Autism = unlist(RDA.clr.Autism)

explained.variance.clr.RDA2.Adenomas = unlist(RDA2.clr.Adenomas)
explained.variance.clr.RDA2.Konzo = unlist(RDA2.clr.Konzo)
explained.variance.clr.RDA2.Autism = unlist(RDA2.clr.Autism)


explained.variance.clr.MOFA.Metabo.Adenomas = sapply(MOFA2.clr.Adenomas, function(x) x$r2_total$group1[1])
explained.variance.clr.MOFA.Metabo.Konzo = sapply(MOFA2.clr.Konzo, function(x) x$r2_total$group1[1])
explained.variance.clr.MOFA.Metabo.Autism = sapply(MOFA2.clr.Autism, function(x) x$r2_total$group1[1])

explained.variance.clr.MOFA.Metabo.Adenomas[explained.variance.clr.MOFA.Metabo.Adenomas<0] = 0
explained.variance.clr.MOFA.Metabo.Konzo[explained.variance.clr.MOFA.Metabo.Konzo<0] = 0
explained.variance.clr.MOFA.Metabo.Autism[explained.variance.clr.MOFA.Metabo.Autism<0] = 0

explained.variance.clr.MOFA.Micro.Adenomas = sapply(MOFA2.clr.Adenomas, function(x) x$r2_total$group1[2])
explained.variance.clr.MOFA.Micro.Konzo = sapply(MOFA2.clr.Konzo, function(x) x$r2_total$group1[2])
explained.variance.clr.MOFA.Micro.Autism = sapply(MOFA2.clr.Autism, function(x) x$r2_total$group1[2])

explained.variance.clr.MOFA.Micro.Adenomas[explained.variance.clr.MOFA.Micro.Adenomas<0] = 0
explained.variance.clr.MOFA.Micro.Konzo[explained.variance.clr.MOFA.Micro.Konzo<0] = 0
explained.variance.clr.MOFA.Micro.Autism[explained.variance.clr.MOFA.Micro.Autism<0] = 0


df.explained.variance.clr = data.frame("Normalization"="clr", "Method"=rep(c("PLS-Can1", "PLS-Can2", "PLS-Reg1", "PLS-Reg2", "CCA1", "CCA2", "RDA1","RDA2", "MOFA2_1", "MOFA2_2"), each="3000"), "Dataset"=rep(rep(c("Adenomas", "Konzo", "Autism"), each=1000),10),
                                       "Explained_Variance"=c(explained.variance.clr.PLS.Can1.Adenomas,explained.variance.clr.PLS.Can1.Konzo,explained.variance.clr.PLS.Can1.Autism,
                                                              explained.variance.clr.PLS.Can2.Adenomas,explained.variance.clr.PLS.Can2.Konzo,explained.variance.clr.PLS.Can2.Autism,
                                                              explained.variance.clr.PLS.Reg1.Adenomas,explained.variance.clr.PLS.Reg1.Konzo,explained.variance.clr.PLS.Reg1.Autism,
                                                              explained.variance.clr.PLS.Reg2.Adenomas,explained.variance.clr.PLS.Reg2.Konzo,explained.variance.clr.PLS.Reg2.Autism,
                                                              explained.variance.clr.CCA1.Adenomas,explained.variance.clr.CCA1.Konzo,explained.variance.clr.CCA1.Autism,
                                                              explained.variance.clr.CCA2.Adenomas,explained.variance.clr.CCA2.Konzo,explained.variance.clr.CCA2.Autism,
                                                              explained.variance.clr.RDA1.Adenomas,explained.variance.clr.RDA1.Konzo,explained.variance.clr.RDA1.Autism,
                                                              explained.variance.clr.RDA2.Adenomas,explained.variance.clr.RDA2.Konzo,explained.variance.clr.RDA2.Autism,
                                                              explained.variance.clr.MOFA.Metabo.Adenomas/100,explained.variance.clr.MOFA.Metabo.Konzo/100,explained.variance.clr.MOFA.Metabo.Autism/100,
                                                              explained.variance.clr.MOFA.Micro.Adenomas/100,explained.variance.clr.MOFA.Micro.Konzo/100,explained.variance.clr.MOFA.Micro.Autism/100))

ggplot(df.explained.variance.clr, aes(x=Method,y=Explained_Variance))+geom_boxplot()+facet_grid(.~Dataset)


explained.variance.alpha.PLS.Can1.Adenomas = sapply(PLS.Can.alpha.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Can1.Konzo = sapply(PLS.Can.alpha.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Can1.Autism = sapply(PLS.Can.alpha.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$Y)
  else NA
})

explained.variance.alpha.PLS.Can2.Adenomas = sapply(PLS.Can.alpha.Adenomas, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.alpha.PLS.Can2.Konzo = sapply(PLS.Can.alpha.Konzo, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.alpha.PLS.Can2.Autism = sapply(PLS.Can.alpha.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$X)
  else NA
})

explained.variance.alpha.PLS.Reg1.Adenomas = sapply(PLS.Reg1.alpha.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Reg1.Konzo = sapply(PLS.Reg1.alpha.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Reg1.Autism = sapply(PLS.Reg1.alpha.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.alpha.PLS.Reg2.Adenomas = sapply(PLS.Reg2.alpha.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Reg2.Konzo = sapply(PLS.Reg2.alpha.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Reg2.Autism = sapply(PLS.Reg2.alpha.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.alpha.CCA1.Adenomas = sapply(cca.alpha.Adenomas, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))
explained.variance.alpha.CCA1.Konzo = sapply(cca.alpha.Konzo, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))
explained.variance.alpha.CCA1.Autism = sapply(cca.alpha.Autism, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))

explained.variance.alpha.CCA2.Adenomas = sapply(cca.alpha.Adenomas, function(x) sum(x$Redundancy$ConditionalRedundancy$X))
explained.variance.alpha.CCA2.Konzo = sapply(cca.alpha.Konzo, function(x) sum(x$Redundancy$ConditionalRedundancy$X))
explained.variance.alpha.CCA2.Autism = sapply(cca.alpha.Autism, function(x) sum(x$Redundancy$ConditionalRedundancy$X))


explained.variance.alpha.RDA1.Adenomas = unlist(RDA.alpha.Adenomas)
explained.variance.alpha.RDA1.Konzo = unlist(RDA.alpha.Konzo)
explained.variance.alpha.RDA1.Autism = unlist(RDA.alpha.Autism)

explained.variance.alpha.RDA2.Adenomas = unlist(RDA2.alpha.Adenomas)
explained.variance.alpha.RDA2.Konzo = unlist(RDA2.alpha.Konzo)
explained.variance.alpha.RDA2.Autism = unlist(RDA2.alpha.Autism)


explained.variance.alpha.MOFA.Metabo.Adenomas = sapply(MOFA2.alpha.Adenomas, function(x) x$r2_total$group1[1])
explained.variance.alpha.MOFA.Metabo.Konzo = sapply(MOFA2.alpha.Konzo, function(x) x$r2_total$group1[1])
explained.variance.alpha.MOFA.Metabo.Autism = sapply(MOFA2.alpha.Autism, function(x) x$r2_total$group1[1])

explained.variance.alpha.MOFA.Metabo.Adenomas[explained.variance.alpha.MOFA.Metabo.Adenomas<0] = 0
explained.variance.alpha.MOFA.Metabo.Konzo[explained.variance.alpha.MOFA.Metabo.Konzo<0] = 0
explained.variance.alpha.MOFA.Metabo.Autism[explained.variance.alpha.MOFA.Metabo.Autism<0] = 0

explained.variance.alpha.MOFA.Micro.Adenomas = sapply(MOFA2.alpha.Adenomas, function(x) x$r2_total$group1[2])
explained.variance.alpha.MOFA.Micro.Konzo = sapply(MOFA2.alpha.Konzo, function(x) x$r2_total$group1[2])
explained.variance.alpha.MOFA.Micro.Autism = sapply(MOFA2.alpha.Autism, function(x) x$r2_total$group1[2])

explained.variance.alpha.MOFA.Micro.Adenomas[explained.variance.alpha.MOFA.Micro.Adenomas<0] = 0
explained.variance.alpha.MOFA.Micro.Konzo[explained.variance.alpha.MOFA.Micro.Konzo<0] = 0
explained.variance.alpha.MOFA.Micro.Autism[explained.variance.alpha.MOFA.Micro.Autism<0] = 0

summary(explained.variance.alpha.MOFA.Micro.Konzo)
df.explained.variance.alpha = data.frame("Normalization"="alpha","Method"=rep(c("PLS-Can1", "PLS-Can2", "PLS-Reg1", "PLS-Reg2", "CCA1", "CCA2", "RDA1","RDA2", "MOFA2_1", "MOFA2_2"), each="3000"), "Dataset"=rep(rep(c("Adenomas", "Konzo", "Autism"), each=1000),10),
                                         "Explained_Variance"=c(explained.variance.alpha.PLS.Can1.Adenomas,explained.variance.alpha.PLS.Can1.Konzo,explained.variance.alpha.PLS.Can1.Autism,
                                                                explained.variance.alpha.PLS.Can2.Adenomas,explained.variance.alpha.PLS.Can2.Konzo,explained.variance.alpha.PLS.Can2.Autism,
                                                                explained.variance.alpha.PLS.Reg1.Adenomas,explained.variance.alpha.PLS.Reg1.Konzo,explained.variance.alpha.PLS.Reg1.Autism,
                                                                explained.variance.alpha.PLS.Reg2.Adenomas,explained.variance.alpha.PLS.Reg2.Konzo,explained.variance.alpha.PLS.Reg2.Autism,
                                                                explained.variance.alpha.CCA1.Adenomas,explained.variance.alpha.CCA1.Konzo,explained.variance.alpha.CCA1.Autism,
                                                                explained.variance.alpha.CCA2.Adenomas,explained.variance.alpha.CCA2.Konzo,explained.variance.alpha.CCA2.Autism,
                                                                explained.variance.alpha.RDA1.Adenomas,explained.variance.alpha.RDA1.Konzo,explained.variance.alpha.RDA1.Autism,
                                                                explained.variance.alpha.RDA2.Adenomas,explained.variance.alpha.RDA2.Konzo,explained.variance.alpha.RDA2.Autism,
                                                                explained.variance.alpha.MOFA.Metabo.Adenomas/100,explained.variance.alpha.MOFA.Metabo.Konzo/100,explained.variance.alpha.MOFA.Metabo.Autism/100,
                                                                explained.variance.alpha.MOFA.Micro.Adenomas/100,explained.variance.alpha.MOFA.Micro.Konzo/100,explained.variance.alpha.MOFA.Micro.Autism/100))

ggplot(df.explained.variance.alpha, aes(x=Method,y=Explained_Variance))+geom_boxplot()+facet_grid(.~Dataset)


df.explained.variance = rbind(df.explained.variance.alpha, df.explained.variance.clr, df.explained.variance.ilr)
unique(df.explained.variance$Method)
pipoo = df.explained.variance %>% mutate(Methode = str_remove(Method,"_"))%>% mutate(Methode = str_remove(Methode,".$"))%>% mutate(mode = str_extract(Method,".$"))

pipoo$Omics = ifelse(pipoo$mode ==1,"Species", "Metabolites")

ggplot(pipoo, aes(x=Methode,y=Explained_Variance, fill=Omics))+geom_boxplot()+scale_fill_manual(values=c("red", "white"))+facet_nested(Normalization~Dataset+Methode, scales = "free_x")+ylab("% Explained Variance")+ggtitle("E.")+theme_bw()+theme(axis.text.x = element_blank(), legend.position = "bottom")+xlab("")

write.table(pipoo, "C:\\Users\\loicm\\Downloads\\data_for_Figure3E.txt", col.names = T, row.names = F, quote=F)

unique(pipoo$Method)


pipa = do.call("rbind",lapply(unique(pipoo$Method), function(x){
  p = do.call("rbind",lapply(unique(pipoo$Dataset), function(y){
    do.call("rbind",lapply(unique(pipoo$Normalization), function(norm){
      data.frame("Method" = x, "Dataset"= y, "Normalization"=norm,"Mean_Explained_Variance"=mean(pipoo[pipoo$Method==x&pipoo$Dataset==y&pipoo$Normalization==norm, "Explained_Variance"], na.rm=T))
      
    }))
  }))
  p
}))

pipa_agg = aggregate(Mean_Explained_Variance~Method+Dataset,pipa, sd)
pipa_agg = pipa_agg %>% mutate(Methode = str_remove(Method,"_"))%>% mutate(Methode = str_remove(Methode,".$"))%>% mutate(mode = str_extract(Method,".$"))

pipa_agg$Omics = ifelse(pipa_agg$mode ==1,"Species", "Metabolites")

ggplot(pipa_agg, aes(x=Methode,y=Mean_Explained_Variance, fill=Omics))+geom_bar(position="dodge", stat="identity",color="black")+scale_fill_manual(values=c("#BC3C29FF", "#0072B5FF"))+facet_grid(.~Dataset)+ylab("STD % Explained Variance")+theme_bw()+theme(legend.position = "bottom")+xlab("")

mean(pipoo[pipoo$Method=="CCA1"&pipoo$Dataset=="Adenomas"&pipoo$Normalization=="alpha", "Explained_Variance"], na.rm=T)
mean(pipoo[pipoo$Method=="CCA1"&pipoo$Dataset=="Adenomas"&pipoo$Normalization=="ilr", "Explained_Variance"], na.rm=T)
mean(pipoo[pipoo$Method=="CCA1"&pipoo$Dataset=="Adenomas"&pipoo$Normalization=="clr", "Explained_Variance"], na.rm=T)

#Multivariate
#CCA
cca.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_CCA_ilr.RDS")
cca.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_CCA_ilr.RDS")

cca.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_CCA_clr_log.RDS")
cca.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_CCA_clr_log.RDS")

cca.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_CCA_alpha_log.RDS")
cca.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_CCA_alpha_log.RDS")

#RDA
RDA.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA1_ilr_Konzo.RDS")
RDA.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA1_ilr_Autism.RDS")

RDA.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA1_clr_Konzo.RDS")
RDA.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA1_clr_Autism.RDS")

RDA.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA1_alpha_Konzo.RDS")
RDA.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA1_alpha_Autism.RDS")

RDA2.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA2_ilr_Konzo.RDS")
RDA2.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA2_ilr_Autism.RDS")

RDA2.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA2_clr_Konzo.RDS")
RDA2.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA2_clr_Autism.RDS")

RDA2.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\RDA2_alpha_Konzo.RDS")
RDA2.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\RDA2_alpha_Autism.RDS")


#MOFA2

MOFA2.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\MOFA2_ilr_Konzo.RDS")
MOFA2.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\MOFA2_ilr_Autism.RDS")

MOFA2.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\MOFA2_clr_Konzo.RDS")
MOFA2.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\MOFA2_clr_Autism.RDS")

MOFA2.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\MOFA2_alpha_Konzo.RDS")
MOFA2.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\MOFA2_alpha_Autism.RDS")

#PLS-Can
PLS.Can.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_canonical_ilr_log.RDS")
PLS.Can.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_canonical_ilr_log.RDS")

PLS.Can.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_canonical_clr_log.RDS")
PLS.Can.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_canonical_clr_log.RDS")

PLS.Can.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_canonical_alpha_log.RDS")
PLS.Can.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_canonical_alpha_log.RDS")


#PLS-Reg
PLS.Reg1.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression1_ilr_log.RDS")
PLS.Reg1.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression1_ilr_log.RDS")

PLS.Reg1.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression1_clr.RDS")
PLS.Reg1.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression1_clr.RDS")

PLS.Reg1.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression1_alpha.RDS")
PLS.Reg1.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression1_alpha.RDS")

PLS.Reg2.ilr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression2_ilr.RDS")
PLS.Reg2.ilr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression2_ilr.RDS")

PLS.Reg2.clr.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression2_clr.RDS")
PLS.Reg2.clr.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression2_clr.RDS")

PLS.Reg2.alpha.Konzo = readRDS("results_multivariate_Konzo\\results_multivariate_Konzo\\results_PLS_regression2_alpha.RDS")
PLS.Reg2.alpha.Autism = readRDS("results_multivariate_Autism\\results_multivariate_Autism\\results_PLS_regression2_alpha.RDS")


explained.variance.ilr.PLS.Can1.Adenomas = sapply(PLS.Can.ilr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Can1.Konzo = sapply(PLS.Can.ilr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Can1.Autism = sapply(PLS.Can.ilr.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$Y)
  else NA
})

explained.variance.ilr.PLS.Can2.Adenomas = sapply(PLS.Can.ilr.Adenomas, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.ilr.PLS.Can2.Konzo = sapply(PLS.Can.ilr.Konzo, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.ilr.PLS.Can2.Autism = sapply(PLS.Can.ilr.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$X)
  else NA
})

explained.variance.ilr.PLS.Reg1.Adenomas = sapply(PLS.Reg1.ilr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Reg1.Konzo = sapply(PLS.Reg1.ilr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Reg1.Autism = sapply(PLS.Reg1.ilr.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.ilr.PLS.Reg2.Adenomas = sapply(PLS.Reg2.ilr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Reg2.Konzo = sapply(PLS.Reg2.ilr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.ilr.PLS.Reg2.Autism = sapply(PLS.Reg2.ilr.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.ilr.CCA1.Adenomas = sapply(cca.ilr.Adenomas, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))
explained.variance.ilr.CCA1.Konzo = sapply(cca.ilr.Konzo, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))
explained.variance.ilr.CCA1.Autism = sapply(cca.ilr.Autism, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))

explained.variance.ilr.CCA2.Adenomas = sapply(cca.ilr.Adenomas, function(x) sum(x$Redundancy$ConditionalRedundancy$X))
explained.variance.ilr.CCA2.Konzo = sapply(cca.ilr.Konzo, function(x) sum(x$Redundancy$ConditionalRedundancy$X))
explained.variance.ilr.CCA2.Autism = sapply(cca.ilr.Autism, function(x) sum(x$Redundancy$ConditionalRedundancy$X))


explained.variance.ilr.RDA1.Adenomas = unlist(RDA.ilr.Adenomas)
explained.variance.ilr.RDA1.Konzo = unlist(RDA.ilr.Konzo)
explained.variance.ilr.RDA1.Autism = unlist(RDA.ilr.Autism)

explained.variance.ilr.RDA2.Adenomas = unlist(RDA2.ilr.Adenomas)
explained.variance.ilr.RDA2.Konzo = unlist(RDA2.ilr.Konzo)
explained.variance.ilr.RDA2.Autism = unlist(RDA2.ilr.Autism)


explained.variance.ilr.MOFA.Metabo.Adenomas = sapply(MOFA2.ilr.Adenomas, function(x) x$r2_total$group1[1])
explained.variance.ilr.MOFA.Metabo.Konzo = sapply(MOFA2.ilr.Konzo, function(x) x$r2_total$group1[1])
explained.variance.ilr.MOFA.Metabo.Autism = sapply(MOFA2.ilr.Autism, function(x) x$r2_total$group1[1])

explained.variance.ilr.MOFA.Metabo.Adenomas[explained.variance.ilr.MOFA.Metabo.Adenomas<0] = 0
explained.variance.ilr.MOFA.Metabo.Konzo[explained.variance.ilr.MOFA.Metabo.Konzo<0] = 0
explained.variance.ilr.MOFA.Metabo.Autism[explained.variance.ilr.MOFA.Metabo.Autism<0] = 0

explained.variance.ilr.MOFA.Micro.Adenomas = sapply(MOFA2.ilr.Adenomas, function(x) x$r2_total$group1[2])
explained.variance.ilr.MOFA.Micro.Konzo = sapply(MOFA2.ilr.Konzo, function(x) x$r2_total$group1[2])
explained.variance.ilr.MOFA.Micro.Autism = sapply(MOFA2.ilr.Autism, function(x) x$r2_total$group1[2])

explained.variance.ilr.MOFA.Micro.Adenomas[explained.variance.ilr.MOFA.Micro.Adenomas<0] = 0
explained.variance.ilr.MOFA.Micro.Konzo[explained.variance.ilr.MOFA.Micro.Konzo<0] = 0
explained.variance.ilr.MOFA.Micro.Autism[explained.variance.ilr.MOFA.Micro.Autism<0] = 0

summary(explained.variance.ilr.MOFA.Micro.Konzo)
df.explained.variance.ilr = data.frame("Normalization"="ilr","Method"=rep(c("PLS-Can1", "PLS-Can2", "PLS-Reg1", "PLS-Reg2", "CCA1", "CCA2", "RDA1","RDA2", "MOFA2_1", "MOFA2_2"), each="3000"), "Dataset"=rep(rep(c("Adenomas", "Konzo", "Autism"), each=1000),10),
                                       "Explained_Variance"=c(explained.variance.ilr.PLS.Can1.Adenomas,explained.variance.ilr.PLS.Can1.Konzo,explained.variance.ilr.PLS.Can1.Autism,
                                                              explained.variance.ilr.PLS.Can2.Adenomas,explained.variance.ilr.PLS.Can2.Konzo,explained.variance.ilr.PLS.Can2.Autism,
                                                              explained.variance.ilr.PLS.Reg1.Adenomas,explained.variance.ilr.PLS.Reg1.Konzo,explained.variance.ilr.PLS.Reg1.Autism,
                                                              explained.variance.ilr.PLS.Reg2.Adenomas,explained.variance.ilr.PLS.Reg2.Konzo,explained.variance.ilr.PLS.Reg2.Autism,
                                                              explained.variance.ilr.CCA1.Adenomas,explained.variance.ilr.CCA1.Konzo,explained.variance.ilr.CCA1.Autism,
                                                              explained.variance.ilr.CCA2.Adenomas,explained.variance.ilr.CCA2.Konzo,explained.variance.ilr.CCA2.Autism,
                                                              explained.variance.ilr.RDA1.Adenomas,explained.variance.ilr.RDA1.Konzo,explained.variance.ilr.RDA1.Autism,
                                                              explained.variance.ilr.RDA2.Adenomas,explained.variance.ilr.RDA2.Konzo,explained.variance.ilr.RDA2.Autism,
                                                              explained.variance.ilr.MOFA.Metabo.Adenomas/100,explained.variance.ilr.MOFA.Metabo.Konzo/100,explained.variance.ilr.MOFA.Metabo.Autism/100,
                                                              explained.variance.ilr.MOFA.Micro.Adenomas/100,explained.variance.ilr.MOFA.Micro.Konzo/100,explained.variance.ilr.MOFA.Micro.Autism/100))
library(ggplot2)
ggplot(df.explained.variance.ilr, aes(x=Method,y=Explained_Variance))+geom_boxplot()+facet_grid(.~Dataset)


explained.variance.clr.PLS.Can1.Adenomas = sapply(PLS.Can.clr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Can1.Konzo = sapply(PLS.Can.clr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Can1.Autism = sapply(PLS.Can.clr.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$Y)
  else NA
})

explained.variance.clr.PLS.Can2.Adenomas = sapply(PLS.Can.clr.Adenomas, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.clr.PLS.Can2.Konzo = sapply(PLS.Can.clr.Konzo, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.clr.PLS.Can2.Autism = sapply(PLS.Can.clr.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$X)
  else NA
})

explained.variance.clr.PLS.Reg1.Adenomas = sapply(PLS.Reg1.clr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Reg1.Konzo = sapply(PLS.Reg1.clr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Reg1.Autism = sapply(PLS.Reg1.clr.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.clr.PLS.Reg2.Adenomas = sapply(PLS.Reg2.clr.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Reg2.Konzo = sapply(PLS.Reg2.clr.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.clr.PLS.Reg2.Autism = sapply(PLS.Reg2.clr.Autism, function(x) sum(x$ConditionalRedundancy$Y))


explained.variance.clr.CCA1.Adenomas = sapply(cca.clr.Adenomas, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$Y)
  else NA
})
explained.variance.clr.CCA1.Konzo = sapply(cca.clr.Konzo, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$Y)
  else NA
})
explained.variance.clr.CCA1.Autism = sapply(cca.clr.Autism, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$Y)
  else NA
})

explained.variance.clr.CCA2.Adenomas = sapply(cca.clr.Adenomas, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$X)
  else NA
})
explained.variance.clr.CCA2.Konzo = sapply(cca.clr.Konzo, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$X)
  else NA
})
explained.variance.clr.CCA2.Autism = sapply(cca.clr.Autism, function(x) {
  if(!all(is.na(x$Redundancy))) sum(x$Redundancy$ConditionalRedundancy$X)
  else NA
})

explained.variance.clr.RDA1.Adenomas = unlist(RDA.clr.Adenomas)
explained.variance.clr.RDA1.Konzo = unlist(RDA.clr.Konzo)
explained.variance.clr.RDA1.Autism = unlist(RDA.clr.Autism)

explained.variance.clr.RDA2.Adenomas = unlist(RDA2.clr.Adenomas)
explained.variance.clr.RDA2.Konzo = unlist(RDA2.clr.Konzo)
explained.variance.clr.RDA2.Autism = unlist(RDA2.clr.Autism)


explained.variance.clr.MOFA.Metabo.Adenomas = sapply(MOFA2.clr.Adenomas, function(x) x$r2_total$group1[1])
explained.variance.clr.MOFA.Metabo.Konzo = sapply(MOFA2.clr.Konzo, function(x) x$r2_total$group1[1])
explained.variance.clr.MOFA.Metabo.Autism = sapply(MOFA2.clr.Autism, function(x) x$r2_total$group1[1])

explained.variance.clr.MOFA.Metabo.Adenomas[explained.variance.clr.MOFA.Metabo.Adenomas<0] = 0
explained.variance.clr.MOFA.Metabo.Konzo[explained.variance.clr.MOFA.Metabo.Konzo<0] = 0
explained.variance.clr.MOFA.Metabo.Autism[explained.variance.clr.MOFA.Metabo.Autism<0] = 0

explained.variance.clr.MOFA.Micro.Adenomas = sapply(MOFA2.clr.Adenomas, function(x) x$r2_total$group1[2])
explained.variance.clr.MOFA.Micro.Konzo = sapply(MOFA2.clr.Konzo, function(x) x$r2_total$group1[2])
explained.variance.clr.MOFA.Micro.Autism = sapply(MOFA2.clr.Autism, function(x) x$r2_total$group1[2])

explained.variance.clr.MOFA.Micro.Adenomas[explained.variance.clr.MOFA.Micro.Adenomas<0] = 0
explained.variance.clr.MOFA.Micro.Konzo[explained.variance.clr.MOFA.Micro.Konzo<0] = 0
explained.variance.clr.MOFA.Micro.Autism[explained.variance.clr.MOFA.Micro.Autism<0] = 0


df.explained.variance.clr = data.frame("Normalization"="clr", "Method"=rep(c("PLS-Can1", "PLS-Can2", "PLS-Reg1", "PLS-Reg2", "CCA1", "CCA2", "RDA1","RDA2", "MOFA2_1", "MOFA2_2"), each="3000"), "Dataset"=rep(rep(c("Adenomas", "Konzo", "Autism"), each=1000),10),
                                       "Explained_Variance"=c(explained.variance.clr.PLS.Can1.Adenomas,explained.variance.clr.PLS.Can1.Konzo,explained.variance.clr.PLS.Can1.Autism,
                                                              explained.variance.clr.PLS.Can2.Adenomas,explained.variance.clr.PLS.Can2.Konzo,explained.variance.clr.PLS.Can2.Autism,
                                                              explained.variance.clr.PLS.Reg1.Adenomas,explained.variance.clr.PLS.Reg1.Konzo,explained.variance.clr.PLS.Reg1.Autism,
                                                              explained.variance.clr.PLS.Reg2.Adenomas,explained.variance.clr.PLS.Reg2.Konzo,explained.variance.clr.PLS.Reg2.Autism,
                                                              explained.variance.clr.CCA1.Adenomas,explained.variance.clr.CCA1.Konzo,explained.variance.clr.CCA1.Autism,
                                                              explained.variance.clr.CCA2.Adenomas,explained.variance.clr.CCA2.Konzo,explained.variance.clr.CCA2.Autism,
                                                              explained.variance.clr.RDA1.Adenomas,explained.variance.clr.RDA1.Konzo,explained.variance.clr.RDA1.Autism,
                                                              explained.variance.clr.RDA2.Adenomas,explained.variance.clr.RDA2.Konzo,explained.variance.clr.RDA2.Autism,
                                                              explained.variance.clr.MOFA.Metabo.Adenomas/100,explained.variance.clr.MOFA.Metabo.Konzo/100,explained.variance.clr.MOFA.Metabo.Autism/100,
                                                              explained.variance.clr.MOFA.Micro.Adenomas/100,explained.variance.clr.MOFA.Micro.Konzo/100,explained.variance.clr.MOFA.Micro.Autism/100))

ggplot(df.explained.variance.clr, aes(x=Method,y=Explained_Variance))+geom_boxplot()+facet_grid(.~Dataset)


explained.variance.alpha.PLS.Can1.Adenomas = sapply(PLS.Can.alpha.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Can1.Konzo = sapply(PLS.Can.alpha.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Can1.Autism = sapply(PLS.Can.alpha.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$Y)
  else NA
})

explained.variance.alpha.PLS.Can2.Adenomas = sapply(PLS.Can.alpha.Adenomas, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.alpha.PLS.Can2.Konzo = sapply(PLS.Can.alpha.Konzo, function(x) sum(x$ConditionalRedundancy$X))
explained.variance.alpha.PLS.Can2.Autism = sapply(PLS.Can.alpha.Autism, function(x) {
  if(!all(is.na(x))) sum(x$ConditionalRedundancy$X)
  else NA
})

explained.variance.alpha.PLS.Reg1.Adenomas = sapply(PLS.Reg1.alpha.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Reg1.Konzo = sapply(PLS.Reg1.alpha.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Reg1.Autism = sapply(PLS.Reg1.alpha.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.alpha.PLS.Reg2.Adenomas = sapply(PLS.Reg2.alpha.Adenomas, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Reg2.Konzo = sapply(PLS.Reg2.alpha.Konzo, function(x) sum(x$ConditionalRedundancy$Y))
explained.variance.alpha.PLS.Reg2.Autism = sapply(PLS.Reg2.alpha.Autism, function(x) sum(x$ConditionalRedundancy$Y))

explained.variance.alpha.CCA1.Adenomas = sapply(cca.alpha.Adenomas, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))
explained.variance.alpha.CCA1.Konzo = sapply(cca.alpha.Konzo, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))
explained.variance.alpha.CCA1.Autism = sapply(cca.alpha.Autism, function(x) sum(x$Redundancy$ConditionalRedundancy$Y))

explained.variance.alpha.CCA2.Adenomas = sapply(cca.alpha.Adenomas, function(x) sum(x$Redundancy$ConditionalRedundancy$X))
explained.variance.alpha.CCA2.Konzo = sapply(cca.alpha.Konzo, function(x) sum(x$Redundancy$ConditionalRedundancy$X))
explained.variance.alpha.CCA2.Autism = sapply(cca.alpha.Autism, function(x) sum(x$Redundancy$ConditionalRedundancy$X))


explained.variance.alpha.RDA1.Adenomas = unlist(RDA.alpha.Adenomas)
explained.variance.alpha.RDA1.Konzo = unlist(RDA.alpha.Konzo)
explained.variance.alpha.RDA1.Autism = unlist(RDA.alpha.Autism)

explained.variance.alpha.RDA2.Adenomas = unlist(RDA2.alpha.Adenomas)
explained.variance.alpha.RDA2.Konzo = unlist(RDA2.alpha.Konzo)
explained.variance.alpha.RDA2.Autism = unlist(RDA2.alpha.Autism)


explained.variance.alpha.MOFA.Metabo.Adenomas = sapply(MOFA2.alpha.Adenomas, function(x) x$r2_total$group1[1])
explained.variance.alpha.MOFA.Metabo.Konzo = sapply(MOFA2.alpha.Konzo, function(x) x$r2_total$group1[1])
explained.variance.alpha.MOFA.Metabo.Autism = sapply(MOFA2.alpha.Autism, function(x) x$r2_total$group1[1])

explained.variance.alpha.MOFA.Metabo.Adenomas[explained.variance.alpha.MOFA.Metabo.Adenomas<0] = 0
explained.variance.alpha.MOFA.Metabo.Konzo[explained.variance.alpha.MOFA.Metabo.Konzo<0] = 0
explained.variance.alpha.MOFA.Metabo.Autism[explained.variance.alpha.MOFA.Metabo.Autism<0] = 0

explained.variance.alpha.MOFA.Micro.Adenomas = sapply(MOFA2.alpha.Adenomas, function(x) x$r2_total$group1[2])
explained.variance.alpha.MOFA.Micro.Konzo = sapply(MOFA2.alpha.Konzo, function(x) x$r2_total$group1[2])
explained.variance.alpha.MOFA.Micro.Autism = sapply(MOFA2.alpha.Autism, function(x) x$r2_total$group1[2])

explained.variance.alpha.MOFA.Micro.Adenomas[explained.variance.alpha.MOFA.Micro.Adenomas<0] = 0
explained.variance.alpha.MOFA.Micro.Konzo[explained.variance.alpha.MOFA.Micro.Konzo<0] = 0
explained.variance.alpha.MOFA.Micro.Autism[explained.variance.alpha.MOFA.Micro.Autism<0] = 0

summary(explained.variance.alpha.MOFA.Micro.Konzo)
df.explained.variance.alpha = data.frame("Normalization"="alpha","Method"=rep(c("PLS-Can1", "PLS-Can2", "PLS-Reg1", "PLS-Reg2", "CCA1", "CCA2", "RDA1","RDA2", "MOFA2_1", "MOFA2_2"), each="3000"), "Dataset"=rep(rep(c("Adenomas", "Konzo", "Autism"), each=1000),10),
                                         "Explained_Variance"=c(explained.variance.alpha.PLS.Can1.Adenomas,explained.variance.alpha.PLS.Can1.Konzo,explained.variance.alpha.PLS.Can1.Autism,
                                                                explained.variance.alpha.PLS.Can2.Adenomas,explained.variance.alpha.PLS.Can2.Konzo,explained.variance.alpha.PLS.Can2.Autism,
                                                                explained.variance.alpha.PLS.Reg1.Adenomas,explained.variance.alpha.PLS.Reg1.Konzo,explained.variance.alpha.PLS.Reg1.Autism,
                                                                explained.variance.alpha.PLS.Reg2.Adenomas,explained.variance.alpha.PLS.Reg2.Konzo,explained.variance.alpha.PLS.Reg2.Autism,
                                                                explained.variance.alpha.CCA1.Adenomas,explained.variance.alpha.CCA1.Konzo,explained.variance.alpha.CCA1.Autism,
                                                                explained.variance.alpha.CCA2.Adenomas,explained.variance.alpha.CCA2.Konzo,explained.variance.alpha.CCA2.Autism,
                                                                explained.variance.alpha.RDA1.Adenomas,explained.variance.alpha.RDA1.Konzo,explained.variance.alpha.RDA1.Autism,
                                                                explained.variance.alpha.RDA2.Adenomas,explained.variance.alpha.RDA2.Konzo,explained.variance.alpha.RDA2.Autism,
                                                                explained.variance.alpha.MOFA.Metabo.Adenomas/100,explained.variance.alpha.MOFA.Metabo.Konzo/100,explained.variance.alpha.MOFA.Metabo.Autism/100,
                                                                explained.variance.alpha.MOFA.Micro.Adenomas/100,explained.variance.alpha.MOFA.Micro.Konzo/100,explained.variance.alpha.MOFA.Micro.Autism/100))

ggplot(df.explained.variance.alpha, aes(x=Method,y=Explained_Variance))+geom_boxplot()+facet_grid(.~Dataset)


df.explained.variance = rbind(df.explained.variance.alpha, df.explained.variance.clr, df.explained.variance.ilr)
unique(df.explained.variance$Method)
pipoo = df.explained.variance %>% mutate(Methode = str_remove(Method,"_"))%>% mutate(Methode = str_remove(Methode,".$"))%>% mutate(mode = str_extract(Method,".$"))

pipoo$Omics = ifelse(pipoo$mode ==1,"Species", "Metabolites")

ggplot(pipoo[pipoo$Dataset!="Adenomas",], aes(x=Methode,y=Explained_Variance, fill=Omics))+geom_boxplot()+scale_fill_manual(values=c("red", "white"))+facet_nested(Normalization~Dataset+Methode, scales = "free_x")+ylab("% Explained Variance")+theme_bw()+theme(axis.text.x = element_blank(), legend.position = "bottom")+xlab("")

write.table(pipoo[pipoo$Dataset!="Adenomas",], "C:\\Users\\loicm\\Downloads\\data_for_Figure3E_no_log.txt", col.names = T, row.names = F, quote=F)


df.MiRKAT.null = data.frame(""sort(c(sapply(MiRKAT_adenomas$null$CLR, function(x) sapply(x, function(y) y$omnibus_p)))),
sort(c(sapply(MiRKAT_autism$null$CLR, function(x) sapply(x, function(y) y$omnibus_p)))),
sort(c(sapply(MiRKAT_konzo$null$CLR, function(x) sapply(x, function(y) y$omnibus_p))))

df.MiRKAT.null.clr = data.frame("Dataset" = c(rep("Adenomas", 120000), rep("Autism", 22000), rep("Konzo", 85000)), 
           "observed"= -log10(c(sort(c(sapply(MiRKAT_adenomas$null$CLR, function(x) sapply(x, function(y) y$omnibus_p)))),sort(c(sapply(MiRKAT_autism$null$CLR, function(x) sapply(x, function(y) y$omnibus_p)))),sort(c(sapply(MiRKAT_konzo$null$CLR, function(x) sapply(x, function(y) y$omnibus_p)))))), 
           "expected"= c(-log10(ppoints(120000)), -log10(ppoints(22000)), -log10(ppoints(85000))),
           "clower"   = c(-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:120000, shape2 = 12000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))),
           "cupper"   = c(-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:120000, shape2 = 12000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))), "Normalization"="CLR")

df.MiRKAT.null.alpha = data.frame("Dataset" = c(rep("Adenomas", 120000), rep("Autism", 22000), rep("Konzo", 85000)), 
                                "observed"= -log10(c(sort(c(sapply(MiRKAT_adenomas$null$Alpha, function(x) sapply(x, function(y) y$omnibus_p)))),sort(c(sapply(MiRKAT_autism$null$Alpha, function(x) sapply(x, function(y) y$omnibus_p)))),sort(c(sapply(MiRKAT_konzo$null$Alpha, function(x) sapply(x, function(y) y$omnibus_p)))))), 
                                "expected"= c(-log10(ppoints(120000)), -log10(ppoints(22000)), -log10(ppoints(85000))),
                                "clower"   = c(-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:120000, shape2 = 12000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))),
                                "cupper"   = c(-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:120000, shape2 = 12000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:22000, shape2 = 22000:1)),-log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:85000, shape2 = 85000:1))),"Normalization"="Alpha")


df.MiRKAT.null = rbind(df.MiRKAT.null.clr, df.MiRKAT.null.alpha)

qrt = ggplot(df.MiRKAT.null) +
  geom_ribbon(
    mapping = aes(x = expected, ymin = clower, ymax = cupper),
    alpha = 0.1
  ) +
  geom_point(aes(expected, observed, col = Normalization), size = 2) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  xlab(expression(paste("Expected -log"[10], plain(P)))) +
  ylab(expression(paste("Observed -log"[10], plain(P)))) +
  theme_bw() +#+theme(legend.position = "none")
  facet_grid(.~Dataset)+ggtitle("A.")

power_MiRKAT_CLR = c(mean(c(sapply(MiRKAT_adenomas$alter$CLR, function(x) sapply(x, function(y) y$omnibus_p)))<=0.05),
mean(c(sapply(MiRKAT_autism$alter$CLR, function(x) sapply(x, function(y) y$omnibus_p)))<=0.05),
mean(unlist(sapply(MiRKAT_konzo$alter$CLR, function(x) sapply(x, function(y) y$omnibus_p)))<=0.05))

df_power_MiRKAT_CLR = data.frame("Dataset" = c("Adenomas","Autism","Konzo"),"power"= power_MiRKAT_CLR, "Normalization"="CLR")

power_MiRKAT_alpha = c(mean(c(sapply(MiRKAT_adenomas$alter$Alpha, function(x) sapply(x, function(y) y$omnibus_p)))<=0.05),
                     mean(c(sapply(MiRKAT_autism$alter$Alpha, function(x) sapply(x, function(y) y$omnibus_p)))<=0.05),
                     mean(unlist(sapply(MiRKAT_konzo$alter$Alpha, function(x) sapply(x, function(y) y$omnibus_p)))<=0.05))

df_power_MiRKAT_alpha = data.frame("Dataset" = c("Adenomas","Autism","Konzo"),"power"= power_MiRKAT_alpha, "Normalization"="Alpha")

df_power_MiRKAT = rbind(df_power_MiRKAT_CLR, df_power_MiRKAT_alpha)

qru = ggplot(df_power_MiRKAT, aes(x=Normalization, y=power, fill=Normalization))+geom_bar(colour="black",stat="identity", position="dodge")+
  facet_grid(.~Dataset)+ylim(c(0,1))+ylab("Power")+theme_bw()+theme(legend.position = "none")+ggtitle("B.")

qrt / qru


alter_Cor_autism_pearson = readRDS("ACAT_cor_pearson_alter_Autism")
alter_Cor_adenomas_pearson = readRDS("ACAT_combined_cor_alter_Adenomas.RDS")
alter_Cor_konzo_pearson = readRDS("ACAT_cor_pearson_alter_Konzo")

alter_Cor_autism_spearman = readRDS("ACAT_combined_cor_alter_Autism.RDS")
alter_Cor_adenomas_spearman = readRDS("ACAT_combined_cor_alter_Adenomas.RDS")
alter_Cor_konzo_spearman = readRDS("ACAT_combined_alter_cor_Konzo.RDS")


pvalues_univariate_alter_FDR= c(mean(c(sapply(1:1000, function(x) p.adjust(univariate_adenomas$`Log-Contrast`$alter[,x], method = "fdr")))<=0.05),
                                mean(c(sapply(1:1000, function(x) p.adjust(univariate_autism$`Log-contrast`$alter[,x], method = "fdr")))<=0.05),
                                mean(unlist(sapply(1:1000, function(x) p.adjust(univariate_konzo$`Log-Contrast`$alter[[x]], method = "fdr")))<=0.05),
                                
                                mean(c(sapply(1:1000, function(y) p.adjust(sapply(MiRKAT_adenomas$alter$ILR[[y]], function(x) x$omnibus_p), method="fdr")))<=0.05),
                                mean(c(sapply(1:1000, function(y) p.adjust(sapply(MiRKAT_autism$alter$ILR[[y]], function(x) x$omnibus_p), method="fdr")))<=0.05),
                                mean(unlist(sapply(1:1000, function(y) p.adjust(sapply(MiRKAT_konzo$alter$ILR[[y]], function(x) x$omnibus_p), method="fdr")))<=0.05),
                                
                                mean(c(sapply(1:1000, function(x) p.adjust(univariate_adenomas$`clr-lm`$alter[,x], method = "fdr")))<=0.05),
                                mean(c(sapply(1:1000, function(x) p.adjust(univariate_autism$`clr-lm`$alter[,x], method = "fdr")))<=0.05),
                                mean(unlist(sapply(1:1000, function(x) p.adjust(univariate_konzo$`clr-lm`$alter[[x]], method = "fdr")))<=0.05),
                                
                                mean(c(apply(aggregate(p.values ~ i, HALLA_pvalues_alter_clr_original_Adenomas, function(x) p.adjust(x, method="fdr"))[,-1], 1, function(row) row<=0.05))),
                                mean(c(apply(aggregate(p.values ~ i, HALLA_pvalues_alter_clr_log_Autism, function(x) p.adjust(x, method="fdr"))[,-1], 1, function(row) row<=0.05))),
                                mean(unlist(apply(aggregate(p.values ~ i, HALLA_pvalues_alter_clr_log_Konzo, function(x) p.adjust(x, method="fdr")), 1, function(row) row$p.values<=0.05))),
                                
                                
                                mean(c(apply(alter_Cor_adenomas_spearman, 2, function(col) p.adjust(col, method="fdr"))<=0.05)),
                                mean(c(apply(alter_Cor_autism_spearman, 2, function(col) p.adjust(col, method="fdr"))<=0.05)),
                                mean(unlist(lapply(alter_Cor_konzo_spearman, function(rep) p.adjust(rep, method="fdr")))<=0.05),
                                
                                mean(c(apply(alter_Cor_adenomas_pearson, 2, function(col) p.adjust(col, method="fdr"))<=0.05)),
                                mean(c(apply(alter_Cor_autism_pearson, 2, function(col) p.adjust(col, method="fdr"))<=0.05)),
                                mean(unlist(lapply(alter_Cor_konzo_pearson, function(rep) p.adjust(rep, method="fdr")))<=0.05)
)


df_power_univariate_FDR = data.frame("Dataset" = rep(c("Adenomas","Autism","Konzo"), 6),"Method" = rep(c("Log-contrast", "MiRKAT", "clr-lm", "HALLA", "Spearman", "Pearson"), each=3),"power"= pvalues_univariate_alter_FDR)



ggplot(df_power_univariate_FDR, aes(x=Method, y=power, fill=Method))+geom_bar(colour="black",stat="identity", position="dodge")+
  facet_grid(.~Dataset)+ylim(c(0,1))+ylab("Power")+theme_bw()+theme(legend.position = "none")



pvalues_univariate_alter_init = c(c(sapply(1:1000, function(x) mean(univariate_adenomas$`Log-Contrast`$alter[,x]<=0.05))),
                                c(sapply(1:1000, function(x) mean(univariate_autism$`Log-contrast`$alter[,x]<=0.05))),
                                unlist(sapply(1:1000, function(x) mean(univariate_konzo$`Log-Contrast`$alter[[x]]<=0.05))),
                                
                                c(sapply(1:1000, function(y) mean(sapply(MiRKAT_adenomas$alter$ILR[[y]], function(x) x$omnibus_p)<=0.05))),
                                c(sapply(1:1000, function(y) mean(sapply(MiRKAT_autism$alter$ILR[[y]], function(x) x$omnibus_p)<=0.05))),
                                unlist(sapply(1:1000, function(y) mean(sapply(MiRKAT_konzo$alter$ILR[[y]], function(x) x$omnibus_p)<=0.05))),
                                
                                c(sapply(1:1000, function(x) mean(univariate_adenomas$`clr-lm`$alter[,x]<=0.05))),
                                c(sapply(1:1000, function(x) mean(univariate_autism$`clr-lm`$alter[,x]<=0.05))),
                                unlist(sapply(1:1000, function(x) mean(univariate_konzo$`clr-lm`$alter[[x]]<=0.05))),
                                
                                c(apply(aggregate(p.values ~ i, HALLA_pvalues_alter_clr_original_Adenomas, function(x) x)[,-1], 1, function(row) mean(row<=0.05))),
                                c(apply(aggregate(p.values ~ i, HALLA_pvalues_alter_clr_log_Autism, function(x) x)[,-1], 1, function(row) mean(row<=0.05))),
                                unlist(apply(aggregate(p.values ~ i, HALLA_pvalues_alter_clr_log_Konzo, function(x) x), 1, function(row) mean(row$p.values<=0.05))),
                                
                                
                                c(apply(alter_Cor_adenomas_spearman, 2, function(col) mean(col<=0.05))),
                                c(apply(alter_Cor_autism_spearman, 2, function(col) mean(col<=0.05))),
                                unlist(lapply(alter_Cor_konzo_spearman, function(rep) mean(rep<=0.05))),
                                
                                c(apply(alter_Cor_adenomas_pearson, 2, function(col) mean(col<=0.05))),
                                c(apply(alter_Cor_autism_pearson, 2, function(col) mean(col<=0.05))),
                                unlist(lapply(alter_Cor_konzo_pearson, function(rep) mean(rep<=0.05))))


df_power_univariate_FDR = data.frame("Dataset" = rep(rep(c("Adenomas","Autism","Konzo"), each=1000),6),"Method" = rep(c("Log-contrast", "MiRKAT", "clr-lm", "HALLA", "Spearman", "Pearson"), each=3000),"power"= pvalues_univariate_alter_init)



ggplot(df_power_univariate_FDR, aes(x=Method, y=power, fill=Method))+geom_boxplot()+
  facet_grid(.~Dataset)+ylim(c(0,1))+ylab("Power")+theme_bw()+theme(legend.position = "none")+geom_hline(yintercept = 0.10)

