library(methylKit)
# Load methylated sites covered in all samples in three generations
setwd("/Volumes/Epiguru/Thesis chapters/final version of each chapter with analysis/sticklbeack/submission to Genetics/heritable cpg in f1f2")
load("./meth.RData")
# Subset F1 and F2 samples
f1f2_subset=read.csv("./stickleback_metadata.csv")
f1f2_subset=f1f2_subset[(f1f2_subset$Generation=="F1" | f1f2_subset$Generation=="F2"), ]
f1f2_subset$Label=as.character(f1f2_subset$Label)
f1f2_subset$Lane=as.factor(f1f2_subset$Lane)
f1f2_subset$Family=as.factor(f1f2_subset$Family)
f1f2_subset$Generation=as.factor(f1f2_subset$Generation)
# levels(f1f2_subset$Generation)=c("1", "0", NA)
f1f2_subset=f1f2_subset[,c(1,3,4,5,7,8)]

# Order is important, adjust order of sample id based on orginal sample id in meth, which is the same of snp individual id
ind=read.table("~/Desktop/heritable cpg in f1f2/f1f2.012.indv")
ind=as.character(ind$V1)

library(dplyr)
f1f2_subset=f1f2_subset %>%
  slice(match(ind, Label))

# Subset methylation levels
meth=reorganize(meth.all,
                sample.ids = as.character(f1f2_subset[,1]),
                treatment = as.numeric(f1f2_subset[,3])) # HL is 1, KL is 2

PCASamples(meth) # check if parent assignment is correct
# KL1M and KL1F are parents of KL_F2
# KL_F1_1 is the family for KL_F2

##########################################
#Exclude SNPs for methylation analyses####
##########################################
bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}

# The C/T and A/G location is produced before purging high LD sites, using parental samples
# load CT SNPs and duplicate column and save as bedfile, result_pat_C_T_maf.bed
CT_maf= read.csv(file="./heritable cpg in f1f2/result_pat_C_T.txt", sep="\t", header=FALSE)
CT_maf = cbind(CT_maf,V3=rep(CT_maf$V2))
write.table(CT_maf,file="./heritable cpg in f1f2/result_pat_C_T_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_CT_bed<- bed_to_granges(file="./heritable cpg in f1f2/result_pat_C_T_mafSTARTEND.bed")

# load GA SNPs and duplicate column and save as bedfile, result_pat_G_A_maf.bed is created in step 5.2
AG_maf= read.csv(file="./heritable cpg in f1f2/result_pat_G_A.txt", sep="\t", header=FALSE)
AG_maf = cbind(AG_maf,V3=rep(AG_maf$V2))
write.table(AG_maf,file="./heritable cpg in f1f2/result_pat_G_A_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_GA_bed<- bed_to_granges(file="./heritable cpg in f1f2/result_pat_G_A_mafSTARTEND.bed")

# interesect bedfile and unite-file
#### create overlap --> these positions are corrected for CT SNPs
unite_norm_10x_GRanges <- as(meth, "GRanges")
Overlap_CT=unite_norm_10x_GRanges[countOverlaps(unite_norm_10x_GRanges, blacklist_CT_bed) == 0L]
meth.1=makeMethylDB(meth,"methylBaseDB")
unite_norm_10x_CT <- selectByOverlap(meth.1, granges(Overlap_CT))

# make methylkit database
objDB_CT=makeMethylDB(unite_norm_10x_CT,"methylBaseDB")

##### create overlap --> these positions are corrected for GA SNPs
unite_norm_10x_CT_GRanges <- as(unite_norm_10x_CT,"GRanges")
Overlap_GA=unite_norm_10x_CT_GRanges[countOverlaps(unite_norm_10x_CT_GRanges, blacklist_GA_bed) == 0L]
unite_norm_10x_CT_GA <- selectByOverlap(objDB_CT, Overlap_GA)

#remove sex chr.
unite_norm_10x_CT_GA_Granges <- as(unite_norm_10x_CT_GA,"GRanges")
unite_norm_10x_CT_GA_Granges_nosex <- unite_norm_10x_CT_GA_Granges[!seqnames(unite_norm_10x_CT_GA_Granges) == 'groupXIX']

unite_norm_10x_CT_GA_DB=makeMethylDB(unite_norm_10x_CT_GA,"methylBaseDB")
unite_norm_10x_CT_GA_nosex_DB <- selectByOverlap(unite_norm_10x_CT_GA_DB,unite_norm_10x_CT_GA_Granges_nosex)

# Check sample id and adjust it if needed
unite_norm_10x_CT_GA_nosex_DB@sample.ids

# Find DMCs between the same strain in F1 and F2 and exclude them
# HL
hl_subset=subset(cov, Origin=="HL")
hl_subset=as.matrix(hl_subset)
hl_subset=data.frame(hl_subset)
unite_norm_10x_CT_GA_nosex_DB_hl=reorganize(unite_norm_10x_CT_GA_nosex_DB, 
                                            sample.ids = as.character(hl_subset[,1]),
                                            treatment = as.numeric(hl_subset[,3]))

myDiff.hl=calculateDiffMeth(unite_norm_10x_CT_GA_nosex_DB_hl,
                            covariates = hl_subset[,2,drop=F], # lane as a covariate
                            mc.cores = 2)

myDiff.sig.hl=getMethylDiff(myDiff.hl, difference=15, qvalue=0.01)
nrow(myDiff.sig.hl) # 137

# KL
kl_subset=subset(cov, Origin=="KL")
kl_subset=as.matrix(kl_subset)
kl_subset=data.frame(kl_subset)
kl_subset=kl_subset[kl_subset$Family=="KL_F1_1" | kl_subset$Family=="KL_F2_1",]
unite_norm_10x_CT_GA_nosex_DB_kl=reorganize(unite_norm_10x_CT_GA_nosex_DB, 
                                            sample.ids = as.character(kl_subset[,1]),
                                            treatment = as.numeric(kl_subset[,3]))

myDiff.kl=calculateDiffMeth(unite_norm_10x_CT_GA_nosex_DB_kl,
                            covariates = kl_subset[,2,drop=F], # lane and family as covariates
                            mc.cores = 2)

myDiff.sig.kl=getMethylDiff(myDiff.kl, difference=15, qvalue=0.01)
nrow(myDiff.sig.kl) # 82

# exclude DMCs to only keep stable (heritable) CpGs
Overlap_DMC.hl=unite_norm_10x_CT_GA_nosex_DB[countOverlaps(as(unite_norm_10x_CT_GA_nosex_DB, "GRanges"), as(myDiff.sig.hl, "GRanges")) == 0L]
unite_norm_10x_CT_GA_nosex_DB.1=makeMethylDB(unite_norm_10x_CT_GA_nosex_DB, "methylBaseDB")
heritable_CpG_no_hl=selectByOverlap(unite_norm_10x_CT_GA_nosex_DB.1, as(Overlap_DMC.hl, "GRanges"))

Overlap_DMC.kl=as(heritable_CpG_no_hl, "GRanges")[countOverlaps(as(heritable_CpG_no_hl, "GRanges"), as(myDiff.sig.kl, "GRanges")) == 0L]
heritable_CpG_no_hl.1=makeMethylDB(heritable_CpG_no_hl, "methylBaseDB")
heritable_CpG_no_hl_kl_DMC=selectByOverlap(heritable_CpG_no_hl.1, as(Overlap_DMC.kl, "GRanges"))

# For ploting purpose
f1=reorganize(heritable_CpG_no_hl_kl_DMC,
              sample.ids = cov[cov$Generation=="F1",]$Label,
              treatment = cov[cov$Generation=="F1",]$Origin)
clusterSamples(f1, plot = T)

f2=reorganize(heritable_CpG_no_hl_kl_DMC,
              sample.ids = cov[cov$Generation=="F2",]$Label,
              treatment = cov[cov$Generation=="F2",]$Origin)
clusterSamples(f2, plot = T)

clusterSamples(heritable_CpG_no_hl_kl_DMC, plot = T)

# Heritable CpGs location----------
heritable_grange=as(heritable_CpG_no_hl_kl_DMC, "GRanges")

# Check if DMCs in parental generation between marine and freshwater ecotypes also heritable
dmcs_p_ecotype=read.table("./dmcs_parental.txt")
colnames(dmcs_p_ecotype)=c("chr", "start", "end")
dmcs_p_ecotype_grange=as(dmcs_p_ecotype, "GRanges")
count=countOverlaps(dmcs_p_ecotype_grange, heritable_grange)
length(count[count>0]) # 845 DMCs between ecotypes are heritable across F1 and F2

# check if the percentage of DMCs that are also stable sites is deviated from the percentage of sites that are stable among all filtered sites
# test if hyper vs. hypo DMCs is biased
library(DescTools)
expected=c(0.996,1-0.996) # 52,729 out of 52,940 are stable sites
observed=c(845/891,1-845/891)

GTest(x=observed,
      p=expected,
      correct="none")

# G = 0.17116, X-squared df = 1, p-value = 0.6791

# Constitutive loci----------
meth.regions=regionCounts(meth.all, regions = as(meth.all, "GRanges"))
data1=data.frame(percMethylation(meth.regions))
data1$chr=as.character(unlist(meth.regions@.Data[1]))
data1$start=as.numeric(unlist(meth.regions@.Data[2]))
data1$end=as.numeric(unlist(meth.regions@.Data[3]))
data1=data1[,c(95:97,1:94)]

# hyper and hypo postion 
hyper=apply(data1[,4:ncol(data1)], 1, function(x)(mean(x, na.rm = T)>90))
cpg_hyper=data1[hyper,]

hyper_position=cpg_hyper[,1:3]
hyper_grange=as(hyper_position, "GRanges")

hypo=apply(data1[,4:ncol(data1)], 1, function(x)(mean(x, na.rm = T)<10))
cpg_hypo=data1[hypo,]

hypo_position=cpg_hypo[,1:3]
hypo_grange=as(hypo_position, "GRanges")

## Count the number of sites that are heritable and are also constitutive hyer/hypo
heritable_hyper_grange=heritable_grange[countOverlaps(heritable_grange, hyper_grange)==1L] # 6462
heritable_hypo_grange=heritable_grange[countOverlaps(heritable_grange, hypo_grange)==1L] # 28005

heritable_variable_grange=heritable_grange[countOverlaps(heritable_grange, hyper_grange)==0L]
heritable_variable_grange=heritable_variable_grange[countOverlaps(heritable_variable_grange, hypo_grange)==0L]

# there are 6462, 28005, and 18262 heritable CpGs that are hyper, hypo, and variable

# Next check the overlap with genomic elements
# First build a reference for all genomic elements in stickleback
# read the gene BED file, prodduced by gfftToGenePred and genePredtoBed, use stickleback gff3 file as input
library(genomation)
gene.obj=readTranscriptFeatures(location = "~/Desktop/heritable cpg in f1f2/stickle.bed")
# replace seqlevels in gene.obj
seqlevels(gene.obj)=gsub("chr", "group", seqlevels(gene.obj))

# General distribution of all heritable sites
Ann.heritable.all=annotateWithGeneParts(heritable_grange, gene.obj)
Ann.heritable.all
# 32.93, 8.88, 13.17, 45.02
# promoter, exons, introns, intergenic regions:
# G = 0.00045324, X-squared df = 1, p-value = 0.983
# G = 4.9368e-05, X-squared df = 1, p-value = 0.9944
# G = 3.4949e-05, X-squared df = 1, p-value = 0.9953
# G = 0.00010099, X-squared df = 1, p-value = 0.992

Ann.heritable.hyper=annotateWithGeneParts(heritable_hyper_grange, gene.obj)
Ann.heritable.hyper
# 11.92, 27.62, 18.86, 41.60
# promoter, exons, introns, intergenic regions:
# G = 23.59, X-squared df = 1, p-value = 1.192e-06 less
# G = 29.26, X-squared df = 1, p-value = 6.329e-08 more
# G = 2.5268, X-squared df = 1, p-value = 0.1119
# G = 0.48902, X-squared df = 1, p-value = 0.4844

Ann.heritable.hypo=annotateWithGeneParts(heritable_hypo_grange, gene.obj)
Ann.heritable.hypo
# 53.37, 2.39, 4.43, 39.81
# promoter, exons, introns, intergenic regions:
# G = 17.827, X-squared df = 1, p-value = 2.419e-05 more
# G = 7.19, X-squared df = 1, p-value = 0.007331 less
# G = 8.7089, X-squared df = 1, p-value = 0.003167 less
# G = 1.1276, X-squared df = 1, p-value = 0.2883

Ann.heritable.variable=annotateWithGeneParts(heritable_variable_grange, gene.obj)
Ann.heritable.variable
# 9.03, 12.19, 24.56, 54.22
# promoter, exons, introns, intergenic regions:
# G = 31.871, X-squared df = 1, p-value = 1.647e-08 less
# G = 1.2094, X-squared df = 1, p-value = 0.2715
# G = 9.3548, X-squared df = 1, p-value = 0.002224 more
# G = 3.3599, X-squared df = 1, p-value = 0.0668

# Null distribution of all filtered CpG in the genome
gCpG=read.table("./all_filtered_CpGs.txt", header = F) # 52,940 CpGs in total in the filtered CpG dataset
gCpG$end=gCpG$V2
colnames(gCpG)[1:2]=c("seqnames", "start")

g.CpG=annotateWithGeneParts(as(gCpG, "GRanges"), gene.obj)
g.CpG
# 32.83, 8.90, 13.19, 45.07
# promoter, exons, introns, intergenic regions

# G-test
library(DescTools)
expected=c(0.4507,1-0.4507)
observed=c(54.22,100-54.22)

GTest(x=observed,
      p=expected,
      correct="none")

# Genes overlapped with heritable loci-----
# extract genomic information from reference genome
library(GenomicFeatures)
library(biomaRt)
library(ChIPpeakAnno)

stickle=makeTxDbFromEnsembl(organism = "Gasterosteus aculeatus")
genes=genes(stickle)

gene_in_heritable_cpg=annotatePeakInBatch(heritable_grange, AnnotationData = genes, output = "overlapping")
gene_in_heritable_cpg=data.frame(gene_in_heritable_cpg)
gene_in_heritable_cpg=gene_in_heritable_cpg[!is.na(gene_in_heritable_cpg$fromOverlappingOrNearest),]
length(unique(gene_in_heritable_cpg$feature)) # 3128

gene_in_heritable_cpg_hyper=annotatePeakInBatch(heritable_hyper_grange, AnnotationData = genes, output = "overlapping")
gene_in_heritable_cpg_hyper=data.frame(gene_in_heritable_cpg_hyper)
gene_in_heritable_cpg_hyper=gene_in_heritable_cpg_hyper[!is.na(gene_in_heritable_cpg_hyper$fromOverlappingOrNearest),]
length(unique(gene_in_heritable_cpg_hyper$feature)) # 1016

gene_in_heritable_cpg_hypo=annotatePeakInBatch(heritable_hypo_grange, AnnotationData = genes, output = "overlapping")
gene_in_heritable_cpg_hypo=data.frame(gene_in_heritable_cpg_hypo)
gene_in_heritable_cpg_hypo=gene_in_heritable_cpg_hypo[!is.na(gene_in_heritable_cpg_hypo$fromOverlappingOrNearest),]
length(unique(gene_in_heritable_cpg_hypo$feature)) # 1450

gene_in_heritable_cpg_variable=annotatePeakInBatch(heritable_variable_grange, AnnotationData = genes, output = "overlapping")
gene_in_heritable_cpg_variable=data.frame(gene_in_heritable_cpg_variable)
gene_in_heritable_cpg_variable=gene_in_heritable_cpg_variable[!is.na(gene_in_heritable_cpg_variable$fromOverlappingOrNearest),]
length(unique(gene_in_heritable_cpg_variable$feature)) # 1905

# Gene ontology test
mart=useDataset("gaculeatus_gene_ensembl", useMart("ensembl"))
# Pool genes are genes associated with filtered CpGs
pool_gene=annotatePeakInBatch(as(unite_norm_10x_CT_GA_nosex_DB, "GRanges"), AnnotationData = genes, output = "overlapping")
pool_gene=data.frame(pool_gene)
pool_gene=pool_gene[!is.na(pool_gene$fromOverlappingOrNearest),]
gene_pool=getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "external_gene_name", "go_id"), 
                values = unique(pool_gene$feature), 
                mart = mart)
heritable_gene=getBM(filters = "ensembl_gene_id", 
                           attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "description"),  
                           values = unique(gene_in_heritable_cpg$feature), 
                           mart = mart)
heritable_hyper_gene=getBM(filters = "ensembl_gene_id", 
                     attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "description"),  
                     values = unique(gene_in_heritable_cpg_hyper$feature), 
                     mart = mart)
heritable_hypo_gene=getBM(filters = "ensembl_gene_id", 
                           attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "description"),  
                           values = unique(gene_in_heritable_cpg_hypo$feature), 
                           mart = mart)
heritable_variable_gene=getBM(filters = "ensembl_gene_id", 
                          attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "description"),  
                          values = unique(gene_in_heritable_cpg_variable$feature), 
                          mart = mart)
# Remove blank entries
gene_pool <- gene_pool[gene_pool$go_id != '',]
heritable_gene <- heritable_hyper_gene[heritable_gene$go_id != '',]
heritable_gene_hyper <- heritable_hyper_gene[heritable_hyper_gene$go_id != '',]
heritable_gene_hypo <- heritable_hypo_gene[heritable_hypo_gene$go_id != '',]
heritable_gene_variable <- heritable_variable_gene[heritable_variable_gene$go_id != '',]

# convert from table format to list format
geneID2GO <- by(gene_pool$go_id,
                gene_pool$ensembl_gene_id,
                function(x) as.character(x))

heriGenes=by(heritable_gene$go_id,
                  heritable_gene$ensembl_gene_id,
                  function(x) as.character(x))

heriGenes.hyper=by(heritable_gene_hyper$go_id,
             heritable_gene_hyper$ensembl_gene_id,
             function(x) as.character(x))

heriGenes.hypo=by(heritable_gene_hypo$go_id,
                   heritable_gene_hypo$ensembl_gene_id,
                   function(x) as.character(x))

heriGenes.var=by(heritable_gene_variable$go_id,
                 heritable_gene_variable$ensembl_gene_id,
                  function(x) as.character(x))
correction<-"fdr"
geneNames = names(geneID2GO)

myInterestingGenesNames=names(heriGenes)
geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))
names(geneList) <- geneNames

setwd("~/Desktop/")

library(topGO)
ontology=c("MF","BP","CC")

# GO analysis for all heritable-------
for (i in 1:length(ontology)) {
  tgData = new("topGOdata", 
               ontology = ontology[i], 
               allGenes = geneList, 
               annot = annFUN.gene2GO, 
               gene2GO = geneID2GO)
  fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
  fisherResCor = p.adjust(score(fisherRes), method=correction)
  weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
  weightResCor = p.adjust(score(weightRes), method=correction)
  allRes = GenTable(tgData, 
                    classic=fisherRes, 
                    weight=weightRes, 
                    orderBy="weight", 
                    ranksOf="classic", 
                    topNodes=200)
  allRes$fisher.COR = fisherResCor[allRes$GO.ID]
  allRes$weight.COR = weightResCor[allRes$GO.ID]
  write.csv(allRes, paste("heritable_genes",ontology[i],"csv",sep="."))
} 

# GO analysis for hyper-------
myInterestingGenesNames.hyper=names(heriGenes.hyper)
geneList.hyper = factor(as.integer(geneNames %in% myInterestingGenesNames.hyper))
names(geneList.hyper) <- geneNames

setwd("~/Desktop/")

library(topGO)
ontology=c("MF","BP","CC")

for (i in 1:length(ontology)) {
  tgData = new("topGOdata", 
               ontology = ontology[i], 
               allGenes = geneList.hyper, 
               annot = annFUN.gene2GO, 
               gene2GO = geneID2GO)
  fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
  fisherResCor = p.adjust(score(fisherRes), method=correction)
  weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
  weightResCor = p.adjust(score(weightRes), method=correction)
  allRes = GenTable(tgData, 
                    classic=fisherRes, 
                    weight=weightRes, 
                    orderBy="weight", 
                    ranksOf="classic", 
                    topNodes=200)
  allRes$fisher.COR = fisherResCor[allRes$GO.ID]
  allRes$weight.COR = weightResCor[allRes$GO.ID]
  write.csv(allRes, paste("heritable_genes.hyper",ontology[i],"csv",sep="."))
} 

myInterestingGenesNames.hypo=names(heriGenes.hypo)
geneList.hypo = factor(as.integer(geneNames %in% myInterestingGenesNames.hypo))
names(geneList.hypo) <- geneNames

setwd("~/Desktop/")

library(topGO)
ontology=c("MF","BP","CC")

# GO analysis for hypo-------
for (i in 1:length(ontology)) {
  tgData = new("topGOdata", 
               ontology = ontology[i], 
               allGenes = geneList.hypo, 
               annot = annFUN.gene2GO, 
               gene2GO = geneID2GO)
  fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
  fisherResCor = p.adjust(score(fisherRes), method=correction)
  weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
  weightResCor = p.adjust(score(weightRes), method=correction)
  allRes = GenTable(tgData, 
                    classic=fisherRes, 
                    weight=weightRes, 
                    orderBy="weight", 
                    ranksOf="classic", 
                    topNodes=200)
  allRes$fisher.COR = fisherResCor[allRes$GO.ID]
  allRes$weight.COR = weightResCor[allRes$GO.ID]
  write.csv(allRes, paste("heritable_genes.hypo",ontology[i],"csv",sep="."))
} 

myInterestingGenesNames.var=names(heriGenes.var)
geneList.var = factor(as.integer(geneNames %in% myInterestingGenesNames.var))
names(geneList.var) <- geneNames

setwd("~/Desktop/")

library(topGO)
ontology=c("MF","BP","CC")

# GO analysis for var-------
for (i in 1:length(ontology)) {
  tgData = new("topGOdata", 
               ontology = ontology[i], 
               allGenes = geneList.var, 
               annot = annFUN.gene2GO, 
               gene2GO = geneID2GO)
  fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
  fisherResCor = p.adjust(score(fisherRes), method=correction)
  weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
  weightResCor = p.adjust(score(weightRes), method=correction)
  allRes = GenTable(tgData, 
                    classic=fisherRes, 
                    weight=weightRes, 
                    orderBy="weight", 
                    ranksOf="classic", 
                    topNodes=200)
  allRes$fisher.COR = fisherResCor[allRes$GO.ID]
  allRes$weight.COR = weightResCor[allRes$GO.ID]
  write.csv(allRes, paste("heritable_genes.var",ontology[i],"csv",sep="."))
} 

# Calculate methylation percentage of DMC in each sample, for the use of ewas-------
meth.heritable.site=regionCounts(unite_norm_10x_CT_GA_nosex_DB, regions = as(heritable_CpG_no_hl_kl_DMC, "GRanges")) # Removed 54044-52729=1315 CpGs due to SNP, sex, or DMCs between generations
perc.meth=percMethylation(meth.heritable.site)
perc.meth.1=as.data.frame(perc.meth)
perc.meth.1=perc.meth.1/100

# Check order of samples, should be the same to Label in covariate and individual in SNPs
colnames(perc.meth.1)

write.table(perc.meth.1, file="./heritable cpg in f1f2/cpg_ewas.pylmm", col.names=F, row.names=F, quote=F, sep = "\t")

####################################
### Heat map for heritable CpGs ####
####################################
library(pheatmap)
sampleinfo=data.frame(Generation=cov$Generation, Strain=cov$Origin)
sampleinfo$Generation=factor(sampleinfo$Generation, levels = c(NA, "F1", "F2"))
sampleinfo$Strain=factor(sampleinfo$Strain, levels = c(NA, "HL", "KL"))
ann_colors=list(Generation=c(F1="#0066FF", F2="black"), Strain=c(HL="#009E73", KL="#D55E00"))
rownames(sampleinfo)=colnames(perc.meth.1)
heatmap=pheatmap(perc.meth.1, 
                 cluster_rows=TRUE,
                 kmeans_k = 5, 
                 cluster_cols=F, 
                 show_rownames=FALSE,
                 show_colnames=TRUE,
                 border_color = NA,
                 fontsize = 8,
                 annotation_col = sampleinfo,
                 annotation_colors = ann_colors,
                 annotation_names_col = F)
heatmap


### SNPs----
raw=read.table("~/Desktop/heritable cpg in f1f2/f1f2.012", row.names = 1) # 362 SNP
raw[raw==-1]="nan"

pos=read.table("~/Desktop/heritable cpg in f1f2/f1f2.012.pos")
sexchr=grep("XIX", pos$V1)
pos=pos[-sexchr,]
colnames(pos)=c("chr", "start")

ind=read.table("~/Desktop/heritable cpg in f1f2/f1f2.012.indv")
ind=ind$V1
ind

raw=raw[,-sexchr]

final=data.frame(t(raw)) # 350 SNP

write.table(final, file = "~/Desktop/heritable cpg in f1f2/snp_gwas.pylmm", row.names = F, col.names = F, quote = F, sep = "\t")

# Output covariate table
cov=read.csv("./stickleback_metadata.csv")
cov=cov[(cov$Generation=="F1" | cov$Generation=="F2"), ]
cov$Label=as.character(cov$Label)
cov$Lane=as.factor(cov$Lane)
cov$Family=as.factor(cov$Family)
cov$Generation=as.factor(cov$Generation)
# levels(cov$Generation)=c("1", "0", NA)
cov=cov[,c(1,3,4,5,7,8)]

# Adjust order of samples based on SNPs in covariate table
library(dplyr)
cov=cov %>%
  slice(match(ind, Label))

cov$Label

write.table(t(cov[,2,drop=F]), file="./heritable cpg in f1f2/cov_ewas.txt", col.names=F, row.names=F, quote=F, sep = "\t")
# No family as covariate because kinship matrix includes such information

# Plot results of significantly associated SNPS---------
ewas=read.table("./heritable cpg in f1f2/CG.ewas", skip = 1)
colnames(ewas)=c("SNP_ID",	"BETA",	"BETA_SD",	"F_STAT",	"P_VALUE")

ewas=cbind(pos, ewas)
ewas$logp=-log10(ewas$P_VALUE)
ewas$fdr=p.adjust(ewas$P_VALUE, method = "BH", n=length(ewas$P_VALUE))
ewas$sig=ifelse(ewas$fdr<0.1, "sig", "nonsig")
# remove 43 scaffold
scaffold=grep("scaffold", ewas$chr)
ewas=ewas[-scaffold,]
ewas$chr=factor(ewas$chr, levels = paste("group", as.roman(c(1:18,20:21)), sep = ""), ordered = T)

library(ggplot2)
ewas.plot=ggplot(ewas, aes(colour=sig, x=start, y=logp))+
  geom_point(size=2)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position = "none",
        axis.text.x=element_text(angle = 90, size = 14),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        panel.border = element_blank())+
  # geom_hline(yintercept = c(4), linetype="dashed")+
  scale_color_manual(values = c("#CCCCCC", "#000000"))+
  labs(x="",y=expression(paste("-", log["10"], italic("P"))))

ewas.plot

ewas.plot1=ewas.plot+
  facet_grid(.~chr, scales = 'free_x', switch = 'x')+
  theme(axis.text.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        panel.spacing = unit(0.1,"lines"))

ewas.plot1

# extract genes associated with sig SNPs
sig_snp=ewas[ewas$sig=="sig",]
sig_snp=sig_snp[,1:2]
sig_snp$end=sig_snp$start

sig_snp_gene=annotatePeakInBatch(as(sig_snp, "GRanges"), AnnotationData = genes, output = "shortestDistance")
sig_snp_gene=data.frame(sig_snp_gene)

sig_snp_gene_name=getBM(filters = "ensembl_gene_id", 
                        attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006", "description"), 
                        values = sig_snp_gene$feature, 
                        mart = mart)
sig_snp_gene_name=merge(sig_snp_gene, sig_snp_gene_name, by.x="feature", by.y="ensembl_gene_id")
sig_snp_gene_name=unique(sig_snp_gene_name)

# names for the six genes are: ENSGACG00000019612, GALNT10, taok3a, erbin, klf5b, epha7

# Distance between sig snps and CpGs in stickleback genome
dist=annotatePeakInBatch(as(sig_snp, "GRanges"), AnnotationData = as(gCpG, "GRanges"), output = "shortestDistance")
dist=data.frame(dist)
