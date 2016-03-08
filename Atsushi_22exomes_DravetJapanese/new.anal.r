library(data.table)
library(plyr)
library(qqman)

#### CHANGE THESE FOR YOUR SYSTEM ####
setwd("~/Dropbox/JapanEpilepsy")
epi=read.table("japan.txt",header=T,sep="\t")
#HGVD data w/ genotypes
jap=read.table("DBexome20131010.tab",sep="\t")

#HGVD data in ANNOVAR, no genotypes
jap2=read.table("DBexome20131010.tab.hgvd2annovar-all.txt")
names(jap)=c("Chr","Start","rsID_freq","Ref.HGVD","Alt.HGVD","nSample","Filter","Mean_depth","SD_depth","RR","RA","AA","NR","NA","Gene")
names(jap2)=c("Chr","Start","End","Ref.HGVD","Alt.HGVD","Alt.freq")

############################################################################################
#First Analysis - Fisher exact test on genotype counts for all mild vs severe - no HGVD

#compute the fisher exact test on genotype counts between mild and severe for all sites in epi
temp=numeric(nrow(epi))
temp2=numeric(nrow(epi))
temp3=numeric(nrow(epi))
for (i in 1:nrow(epi)) {
  if (i%%10000==0) { print(i)}
  temp[i]=fisher.test(rbind(c(epi$mild00[i],epi$mild01[i],epi$mild11[i]),c(epi$severe00[i],epi$severe01[i],epi$severe11[i])))$p.value
  temp2[i]=fisher.test(rbind(c(epi$mild00[i],epi$mild01[i]+epi$mild11[i]),c(epi$severe00[i],epi$severe01[i]+epi$severe11[i])))$p.value
  temp3[i]=fisher.test(rbind(c(epi$mild00[i]+epi$mild01[i],epi$mild11[i]),c(epi$severe00[i]+epi$severe01[i],epi$severe11[i])))$p.value
}
epi$pval.geno=temp
epi$pval.dom=temp2
epi$pval.rec=temp3

#to remove pseudogenes
epi.t=epi[ -grep("GL", epi$Chr) , ]

#rename 'X', 'Y', 'MT' in "Chr" column to numbers 23, 24, 25
epi.t$Chr <- gsub('X', '23', epi.t$Chr)
epi.t$Chr <- gsub('Y', '24', epi.t$Chr)
epi.t$Chr <- gsub('MT', '25', epi.t$Chr)

write.table(epi.t,"GenotypeTest.snps.txt",quote=F,sep="\t",row.names=F)

#manhattan plot
#t=epi;t$chr=as.numeric(t$Chr), adds new column "chr" with numeric values for "Chr", 
#commented out because numeric conversion changes values i.e. 2 in "Chr" becomes 12 in "chr"

Test1=read.table("GenotypeTest.snps.txt", header=T, sep = "\t")
manhattan(Test1,chr="Chr",bp="Start",p="pval.geno",logp=T)

############################################################################################
#Second Analysis - PBS for mild and severe

#Use inner join to get sites in common for HGVD and epilepsy patients
#note that indels and ref/alt issues lead to many dropped sites
d=join(epi,jap,type="inner",by=c("Chr","Start")) 
#66111 sites after join

#check ref/alt between HGVD and epi
t=d[which(as.character(d$Ref)!=as.character(d$Ref.HGVD)),c(1,2,3,4,5,17,18)]
t2=d[which(as.character(d$Alt)!=as.character(d$Alt.HGVD)),c(1,2,3,4,5,17,18)]

#remove sites w/ Ref/Alt diffs
d=d[-which(as.character(d$Ref)!=as.character(d$Ref.HGVD)),]
d=d[-which(as.character(d$Alt)!=as.character(d$Alt.HGVD)),]
#54293 sites after above filter

#remove sites from HGVD with <100 individuals
d=d[-which(d$nSample<100),]
#54020 sites after above filter

#compute allele frequencies
#p is freq of Ref allele, q is freq of Alt allele
d$p.mild=(2*d$mild00+d$mild01)/(2*(d$mild00+d$mild01+d$mild11))
d$q.mild=1-d$p.mild
d$p.severe=(2*d$severe00+d$severe01)/(2*(d$severe00+d$severe01+d$severe11))
d$q.severe=1-d$p.severe
d$p.jap=(2*d$RR+d$RA)/(2*(d$RR+d$RA+d$AA))
d$q.jap=1-d$p.jap

#compute FST from the genotype frequencies
#compute sample sizes at each site
d$n.mild=d$mild00+d$mild01+d$mild11
d$n.severe=d$severe00+d$severe01+d$severe11
d$n.jap=d$RR+d$RA+d$AA

#compute allele freqs for groups
d$p.mild.severe=(d$p.mild*d$n.mild + d$p.severe*d$n.severe)/(d$n.mild+d$n.severe)
d$q.mild.severe=1-d$p.mild.severe

d$p.mild.jap=(d$p.mild*d$n.mild + d$p.jap*d$n.jap)/(d$n.mild+d$n.jap)
d$q.mild.jap=1-d$p.mild.jap

d$p.severe.jap=(d$p.severe*d$n.severe + d$p.jap*d$n.jap)/(d$n.severe+d$n.jap)
d$q.severe.jap=1-d$p.severe.jap

#compute expected het
d$mild.exp.het=2*d$p.mild*d$q.mild
d$severe.exp.het=2*d$p.severe*d$q.severe
d$jap.exp.het=2*d$p.jap*d$q.jap

#compute Hs and Ht
d$Hs.mild.severe=(d$mild.exp.het * d$n.mild + d$severe.exp.het * d$n.severe)/(d$n.mild + d$n.severe)
d$Ht.mild.severe=2*d$p.mild.severe*d$q.mild.severe

d$Hs.mild.jap=(d$mild.exp.het*d$n.mild + d$jap.exp.het*d$n.jap)/(d$n.mild+d$n.jap)
d$Ht.mild.jap=2*d$p.mild.jap*d$q.mild.jap

d$Hs.severe.jap=(d$jap.exp.het*d$n.jap + d$severe.exp.het*d$n.severe)/(d$n.jap+d$n.severe)
d$Ht.severe.jap=2*d$p.severe.jap*d$q.severe.jap

#compute FST (really Nei's GST=(Ht-Hs)/Ht )
d$FST.mild.severe=(d$Ht.mild.severe-d$Hs.mild.severe)/d$Ht.mild.severe
d$FST.mild.jap=(d$Ht.mild.jap-d$Hs.mild.jap)/d$Ht.mild.jap
d$FST.severe.jap=(d$Ht.severe.jap-d$Hs.severe.jap)/d$Ht.severe.jap

#diff=alt.freq.severe-alt.freq.mild 
d$diff.severe.mild=d$q.severe-d$q.mild
d$pbs.severe=(1/2)*(-log(1-d$FST.mild.severe)-log(1-d$FST.severe.jap)+log(1-d$FST.mild.jap))
d$pbs.mild=(1/2)*(-log(1-d$FST.mild.severe)-log(1-d$FST.mild.jap)+log(1-d$FST.severe.jap))

#by gene computation
d.dt=data.table(d)
d.gene=d.dt[,list(chr=unique(Chr),
                   posStart=min(Start),
                   posEnd=max(End),
                   pbs.mild.mean=mean(pbs.mild),
                   pbs.mild.max=max(pbs.mild),
                   pbs.severe.mean=mean(pbs.severe),
                   pbs.severe.max=max(pbs.severe),
                   diff.mean=mean(diff.severe.mild),
                   diff.max=max(diff.severe.mild),
                   fst.mean=mean(FST.mild.severe),
                   fst.max=max(FST.mild.severe),
                   pval.mean=mean(pval.geno),
                   pval.min=min(pval.geno),
                   nSNPs=length(pbs.mild)),
             by=Gene.refGene]

#to remove NA values (remember to insert/edit column numbers when appropriate)
d.r=d[rowSums(is.na(d[,c(59,60)]))==0,] 

#rename 'X', 'Y', 'MT' in "Chr" column to numbers 23, 24, 25
d.r$Chr <- gsub('X', '23', d.r$Chr)
d.r$Chr <- gsub('Y', '24', d.r$Chr)
d.r$Chr <- gsub('MT', '25', d.r$Chr)

write.table(d.r,"PBS.snps.txt",quote=F,sep="\t",row.names=F)
write.table(d.gene,"PBS.byGene.txt",quote=F,sep="\t",row.names=F)

#manhattan plot
#d.t=d;d.t$chr=as.numeric(d.t$Chr), adds new column "chr" with numeric values for "Chr", 
#commented out because numeric conversion changes values i.e. 2 in "Chr" becomes 12 in "chr"

Test2=read.table("PBS.snps.txt", header=T, sep = "\t")
manhattan(Test2,chr="Chr",bp="Start",p="pbs.mild",logp=F)
manhattan(Test2,chr="Chr",bp="Start",p="pbs.severe",logp=F)

############################################################################################
####THIRD ANALYSIS
#use inner join to get sites in common for HGVD and epilepsy patients
#d2 is second round using ANNOVAR formated HGVD from Laurel - lost genotypes but still have alt allele freq
d2=join(epi,jap2,type="inner",by=c("Chr","Start")) 
#80858 sites after join

#check ref/alt between HGVD and epi
u=d2[which(as.character(d2$Ref)!=as.character(d2$Ref.HGVD)),c(1,2,3,4,5,17,18)]
u2=d2[which(as.character(d2$Alt)!=as.character(d2$Alt.HGVD)),c(1,2,3,4,5,17,18)]

#remove sites w/ Ref/Alt diffs
d2=d2[-which(as.character(d2$Ref)!=as.character(d2$Ref.HGVD)),]
d2=d2[-which(as.character(d2$Alt)!=as.character(d2$Alt.HGVD)),]
#66425 sites after above filter

d2$n.mild=d2$mild00+d2$mild01+d2$mild11
d2$n.severe=d2$severe00+d2$severe01+d2$severe11

#add in allele count diff and prob of allele count diff cols
#alt allele count diff
d2$alt.allele.count.mild=2*d2$mild11+d2$mild01
d2$alt.allele.count.severe=2*d2$severe11+d2$severe01
d2$alt.allele.count.diff=d2$alt.allele.count.mild-d2$alt.allele.count.severe

#prob of alt allele count diff
temp=numeric(nrow(d2))
for (i in 1:nrow(d2)) {
  if (i%%10000==0) { print(i)}
  temp[i]=diffBin(d2$alt.allele.count.diff[i], d2$n.mild[i]*2, d2$n.severe[i]*2, d2$Alt.freq[i])
}
d2$diff.pvalue=temp

d2$diff.pvalue.corrected=d2$diff.pvalue*nrow(d2)

d2$Chr <- gsub('X', '23', d2$Chr)
d2$Chr <- gsub('Y', '24', d2$Chr)
d2$Chr <- gsub('MT', '25', d2$Chr)

write.table(d2,"AlleleFreqDiffTest.snps.txt",quote=F,sep="\t",row.names=F)

Test3=read.table("AlleleFreqDiffTest.snps.txt", header=T, sep = "\t")

manhattan(Test3,chr="Chr",bp="Start",p="diff.pvalue",logp=T)

#function to compute prob of diff between two binomials both with/ prob p
diffBin<-function(diff,n1,n2,p){
  prob=0
  if (diff>=0){  
    for (i in 0:n1){     
      prob=prob+dbinom(i+diff,n1,p)*dbinom(i,n2,p)
    }
  }
  else
  {
    for (i in 0:n2){     
      prob=prob+dbinom(i+diff,n1,p)*dbinom(i,n2,p)
    }
  }
  return(prob)
}


