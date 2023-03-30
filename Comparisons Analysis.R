#Programmatic access to EWAS Catalogue datasets
setwd("/Users/claudiabredemeyer/Desktop/Us vs Them")
library(readr)
TheirAD<- read_tsv("alzheimers_disease.tsv")
TheirAmyloidosis <- read_tsv("burden_of_neuritic_amyloid_plaques.tsv")
TheirBraak <- read_tsv("alzheimers_disease_braak_stage.tsv")
TheirPsychosis <- read_tsv("psychotic_experiences.tsv")

#We also took only values < e-4

#Accessing Their Data
ad_cpg <- TheirAD$cpg  #De Jager et al 
amy_cpg <- TheirAmyloidosis$cpg #De Jager et al 
braak_cpg <- TheirBraak$cpg #De Jager, Lunnon 2014, L Zhang 2020 
psych_cpg <- TheirPsychosis$cpg #S Roberts et al, Longitudinal adolescent psychosis 

ad_genes <- unique(unlist(str_split(TheirAD$gene, ";")))
amy_genes <- unique(unlist(str_split(TheirAmyloidosis$gene, ";")))
braak_genes <- unique(unlist(str_split(TheirBraak$gene, ";")))
psych_genes <- unique(unlist(str_split(TheirPsychosis$gene, ";")))

#Accessing our data
OurDMPs <- read.csv("SigDMP.csv")
our_cpg <- OurDMPs$Name
dmp_freq <- read.csv("DMP_Freq.csv")
our_genes <- dmp_freq$gene_vec

#Testing for common CpGs 
intersect(our_cpg, ad_cpg)
intersect(our_cpg, amy_cpg)
intersect(our_cpg, psych_cpg)
intersect(our_cpg, braak_cpg)

#Testing for common genes
intersect(our_genes, ad_genes)
intersect(our_genes, amy_genes)
intersect(our_genes, psych_genes)
intersect(our_genes, braak_genes)
#No common genes except in Braak staging, wherein we found ANKRD35; CCDC42B; 
# FBXL18; NAV1; NEURL1B; SEMA4B; ZFAND3

#Accessing Pishva et al data 
pishva <- read.csv("Pishva.csv")
pishva <- subset(pishva, pishva$p<0.0001)
pishva_cpg <- pishva$CpGs
pishva_genes <- pishva$Genes

intersect(our_genes, pishva_genes)
intersect(our_cpg, pishva_cpg)
# No common genes or cpgs :c 
