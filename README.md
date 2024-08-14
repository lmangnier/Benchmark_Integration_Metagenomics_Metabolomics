In this repo you can find the code for reproducing results and figures from our paper "*Decoding the Microbiome-Metabolome Nexus: A Systematic Benchmark of Integrative Strategies*". We summarize important results and findings below. We also provide a complete tutorial in the Wiki section on combining the best methods to identify complementary biological processes in one real dataset. 


## Work Summary
In our work we compared 19 statistical methods across 3 synthetic datasets to identify the most promising methods across a variety of scientific questions, such as **global associations**, **data summarization**, **individual associations**, and **feature selection**. Datasets considered exhibit various structures. We illustrated best methods in a real-world data application for Konzo. Below we summarize the most important findings, providing practical guidelines for helping the community.

## Global Associations

Global associations refers to the analyze of statistical relationships occuring globally. Compared to indivual associations where we fitted a model between each pair of microbiota-metabolite, global associations provides a global measure of the relationship occuring between two datasets. We evaluated three methods: the Mantel test, MMiRKAT and the Procrustes Analysis. 
Methods were compared both on the Type-I error rate and Power. 


![Figure3_GA_github](https://github.com/user-attachments/assets/ba97f958-1005-42c0-9ec7-1fbe8692a744)

*Important Findings*:
- the Mantel test is the most powerful method controlling-well the Type-I error rate across our scenarios and microbiome normalizations. That is providing the accurate amount of false positives while detecting significant signals when an association is expected. 
- to improve result transferability ILR transformation should represent the default choice when performing global associations with the Mantel test
- MMiRKAT could be considered if correction for bias is expected, the method cannot deal with P>>N scenarios, however. Procrustes Analysis offers graphical vizualization which is interesting in certain scenarios.

## Data Summarization

Data Summarization refers to the process of recapitulating data variability accounting for the correlation between the two omics. In our benchmark we considered five different methods: CCA, PLS-Reg, PLS-Can, RDA, and MOFA2. Methods were compared based on the data variability contained in the latent factors. 

![Figure3_DS_github](https://github.com/user-attachments/assets/552bf666-d982-43f5-94b7-f31ae2dd7fde)

*Important Findings*:
- RDA captures most variability across latent factors, while being robust to microbiome normalizations.
- Directionality is a key driver of method performance depending on the underlying data structure, particularly for regression-based approaches (PLS-Reg or RDA). 

## Individual Associations

Unlike Global Associations, Individual Associations provide a measure of association between one metabolite and one or more species. We exploited different strategies, for evaluating individual relationships: clr-lm, MiRKAT, Log-contrast and HALLA. Methods were compared against Spearman's or Pearson's correlations regarding on the Type-I error rate and power.

![Rplot67](https://github.com/user-attachments/assets/ef0c8107-bc57-40fe-8b4c-7d4d82b4a5fe)

*Important Findings*:
- MiRKAT is the best method to study the impact of species on metabolite, providing a collective framework robust to data normalizations. Post-hoc analyses have to be applied for detecting the microbial drivers. clr-lm is our default choice in this case.
- Unlike other methods, MiRKAT is able to provide a framework to adjust for confounding.
- Cautious is important here, since systematic applications at large scale may suffer from extreme running times and difficut interpretations.  

## Feature Selection 

Feature selection refers to the statistical practice of subsetting important features associated with a variable. In our context, we leverage on two different strategies: univariate and multivariate feature selection methods. We compared three different methods for the two strategies: CODA-LASSO, clr-LASSO and clr-MLASSO, and sCCA, sPLS-Can and sPLS-Reg, for univariate and multivariate approaches, respectively. We considered evaluating methods for compositional predictors solely.

![figure_5 (1)](https://github.com/user-attachments/assets/f04ab27d-6e91-40da-9517-28f0d9b5c2c8)



*Important Findings*:
- **Univariate**: CODA-LASSO is the best trade-off between sparsity and reliability to select the core species associated with individual metabolites. However the method is strongly impacted by the underlying data structure such as the proportion of zeros.
- **Multivariate**: sPLS-Reg is the best method to identify core species when accounting for the between-omics correlation.



## Reproducibility 
Data for the different methods benchmarked in the paper can be found at: https://doi.org/10.6084/m9.figshare.25234915. 
Details on results and methods can be found at: https://doi.org/10.1101/2024.01.26.577441

## Contact

Questions and comments should be adressed to: loic.mangnier@gmail.com

