In this repo you can find the code for reproducing results and figures from our paper "Decoding the Microbiome-Metabolome Nexus: A Systematic Benchmark of Integrative Strategies".

## Work Summary
In our work we compared 19 statistical methods across 3 synthetic datasets to identify the most promising methods across a variety of scientific questions, such as **global associations**, **data summarization**, **individual associations**, and **feature selection**. We illustrated best methods in a real-world data application for Konzo. Below we summarize the most important findings, providing practical guidelines for helping the community.

## Global Associations

Global associations refers to the analyze of statistical relationships occuring globally. Compared to indivual association where we fitted a model between each pair of microbiota-metabolite, global associations provides a global measure between two datasets. We evaluated three methods: the Mantel test, MMiRKAT and the Procrustes Analysis. 
Methods were compared both on the Type-I error rate and Power. 


![Figure3_GA_github](https://github.com/user-attachments/assets/ba97f958-1005-42c0-9ec7-1fbe8692a744)

*Important Findings*:
- the Mantel test is the most powerful method controlling-well the Type-I error rate across our scenarios and microbiome normalizations. That is providing the accurate amount of false positives while detecting significant signals when an association is expected. 
- to improve result transferability ILR transformation should represent the default use when performing global associations with the Mantel test
- MMiRKAT could be considered if correction for bias is expected, the method cannot deal with P>>N scenarios. Procrustes Analysis offers graphical vizualization which is interesting in certain scenarios.

## Data Summarization

Data Summarization refers to the process of recapitulating data variability accounting for the correlation between the two omics. In our benchmark we considered five different methods: CCA, PLS-Reg, PLS-Can, RDA, and MOFA2. Methods were compared based on the variability contained in the latent factors. 

![Figure3_DS_github](https://github.com/user-attachments/assets/552bf666-d982-43f5-94b7-f31ae2dd7fde)



Data for the different methods benchmarked in the paper can be found at: https://doi.org/10.6084/m9.figshare.25234915. 
Details on results and methods can be found at: https://doi.org/10.1101/2024.01.26.577441
We provide a complete tutorial in the Wiki section on how combining best methods to identify complementary biological processes on one real dataset:



