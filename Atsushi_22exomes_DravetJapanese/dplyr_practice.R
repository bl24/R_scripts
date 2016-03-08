# dplyr practice
# Used for 22 Japanese exome cohort in Dropbox/JapanEpilepsy/Branden_Ranalysis/2016_Jan_(10mild:12severe)
# 2016-Feb-01

install.packages('dplyr')
library(dplyr)

# For Test1
setwd("~/Dropbox/JapanEpilepsy/Branden_Ranalysis/2016_Jan_(10mild:12severe)/Test1")
GenoTest=read.table("GenotypeTest.snps.txt", header=T, sep="\t")

# To filter only exonic variants
ExonicTest1 = GenoTest %>% filter(Func.refGene=="exonic")
manhattan(ExonicTest1,chr="Chr",bp="Start",p="pval.geno",logp=T, main="Exonic Variants (pval.geno)")

# To filter only UTR variants
UTRTest1 = GenoTest %>% filter(Func.refGene=="UTR3"|Func.refGene=="UTR5")
manhattan(UTRTest1,chr="Chr",bp="Start",p="pval.geno",logp=T, main="UTR Variants (pval.geno)")

ExonUTRTest1 = GenoTest %>% filter(Func.refGene=="UTR3"|Func.refGene=="UTR5"|Func.refGene=="exonic")

