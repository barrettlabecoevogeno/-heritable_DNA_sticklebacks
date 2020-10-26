library(methylKit)
load("./meth.RData")
# Subset HL_F2 samples
f2_subset=read.csv("./stickleback_metadata.csv")
f2_subset=f2_subset[f2_subset$Generation=="F2" & f2_subset$Origin=="HL", ]
f2_subset$Label=as.character(f2_subset$Label)
f2_subset$Lane=as.factor(f2_subset$Lane)
f2_subset=f2_subset[,c(1,3,7)]

# Adjust order of samples based on SNP individual order
library(dplyr)
snps_ind=read.table("./hl_f2.012.indv")
snps_ind$V1=as.character(snps_ind$V1)
f2_subset=f2_subset %>%
  slice(match(snps_ind$V1, Label))

# Subset methylation levels
meth=reorganize(meth.all,
                sample.ids = as.character(f2_subset[,1]),
                treatment = as.numeric(f2_subset[,3])) # HL is 2

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

# load CT SNPs and duplicate column and save as bedfile, result_pat_C_T_maf.bed
CT_maf= read.csv(file="./result_pat_C_T.txt", sep="\t", header=FALSE)
CT_maf = cbind(CT_maf,V3=rep(CT_maf$V2))
write.table(CT_maf,file="./result_pat_C_T_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blachlist_CT_bed<- bed_to_granges(file="./result_pat_C_T_mafSTARTEND.bed")

# load GA SNPs and duplicate column and save as bedfile, result_pat_G_A_maf.bed is created in step 5.2
AG_maf= read.csv(file="./result_pat_G_A.txt", sep="\t", header=FALSE)
AG_maf = cbind(AG_maf,V3=rep(AG_maf$V2))
write.table(AG_maf,file="./result_pat_G_A_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blachlist_GA_bed<- bed_to_granges(file="./result_pat_G_A_mafSTARTEND.bed")

# interesect bedfile and unite-file
#### create overlap --> these positions are corrected for CT SNPs
unite_norm_10x_GRanges <- as(meth, "GRanges")
Overlap_CT=unite_norm_10x_GRanges[countOverlaps(unite_norm_10x_GRanges, blachlist_CT_bed) == 0L]
meth.1=makeMethylDB(meth,"methylBaseDB")
unite_norm_10x_CT <- selectByOverlap(meth.1, granges(Overlap_CT))

# make methylkit database
objDB_CT=makeMethylDB(unite_norm_10x_CT,"methylBaseDB")

##### create overlap --> these positions are corrected for GA SNPs
unite_norm_10x_CT_GRanges <- as(unite_norm_10x_CT,"GRanges")
Overlap_GA=unite_norm_10x_CT_GRanges[countOverlaps(unite_norm_10x_CT_GRanges, blachlist_GA_bed) == 0L]
unite_norm_10x_CT_GA <- selectByOverlap(objDB_CT, Overlap_GA)

#remove sex chr.
unite_norm_10x_CT_GA_Granges <- as(unite_norm_10x_CT_GA,"GRanges")
unite_norm_10x_CT_GA_Granges_nosex <- unite_norm_10x_CT_GA_Granges[!seqnames(unite_norm_10x_CT_GA_Granges) == 'groupXIX']

unite_norm_10x_CT_GA_DB=makeMethylDB(unite_norm_10x_CT_GA,"methylBaseDB")
unite_norm_10x_CT_GA_nosex_DB <- selectByOverlap(unite_norm_10x_CT_GA_DB,unite_norm_10x_CT_GA_Granges_nosex)

# Totally removed 54044-52940=1104 sites: C/T, G/A SNPs or CpGs in LG19

# Calculate methylation levels
meth.count=regionCounts(unite_norm_10x_CT_GA_nosex_DB, regions = as(unite_norm_10x_CT_GA_nosex_DB, "GRanges"))
perc.meth.all=percMethylation(meth.count)

colnames(perc.meth.all) == snps_ind$V1 # double check order of samples

cpg=data.frame(perc.meth.all)

cpg$var=apply(cpg, 1, function(x)(max(x)-min(x)))
nrow(cpg[cpg$var>10,])
keep=as.numeric(rownames(cpg[cpg$var>10,]))
cpg=cpg[cpg$var>10,]
cpg$var=NULL

cpg=cpg[,snps_ind$V1]/100 # make sure order of samples in cpg is the same to snp
rownames(cpg)=paste(meth.count$chr[keep], meth.count$start[keep], sep = ",")
cpg$CpGid=rownames(cpg)
cpg=cpg[,c(29,1:28)]

write.table(cpg, "./cpg.hl.f2.txt", row.names = F, sep = "\t", quote = F)

# exclude constitutive loci across marine and freshwater
# cpg$mean=apply(cpg, 1, mean)
# constitutive=cpg$mean>0.9 | cpg$mean<0.1
# cpg_filter=cpg[!constitutive,]
# cpg_filter$mean=NULL

# cpg_filter$CpGid=rownames(cpg_filter)
# cpg_filter=cpg_filter[,c(65, 1:64)] # adjust columns

# write.table(cpg_filter, "./f2/cpg.f2.txt", row.names = F, sep = "\t", quote = F)

#### Second, make a matrix of snp-------
snps=read.table("./hl_f2.012", row.names = 1) # 546 SNPs

# exclude snps on sex chromosome, LG19
snps_pos=read.table("./hl_f2.012.pos")
sexchr=grep("XIX", snps_pos$V1)
snps=snps[,-sexchr]

snps_pos=snps_pos[-sexchr,]
snps_pos=paste(snps_pos$V1, snps_pos$V2, sep = ",")

snps_ind=read.table("./hl_f2.012.indv")
snps_ind$V1=as.character(snps_ind$V1)
rownames(snps)=snps_ind$V1
colnames(snps)=snps_pos

rownames(snps)

snps=data.frame(t(snps))
snps$snpid=rownames(snps)
snps=snps[,c(29, 1:28)]
snps[snps==-1]=NA

write.table(snps, "./snps.hl.f2.txt", row.names = F, sep = "\t", quote = F) # 2692 SNPs

#### Last make a matrix of covariate(s)-------
covariate=read.csv("./stickleback_metadata.csv")
covariate=covariate[covariate$Generation=="F2" & covariate$Origin=="HL",]
covariate=covariate[,c(1,3)]
covariate$Lane=factor(covariate$Lane)
covariate$Label=as.character(covariate$Label)

library(dplyr)
covariate=covariate %>%
  slice(match(snps_ind$V1, Label))
covariate=data.frame(t(covariate))

write.table(covariate, "./covariate.hl.f2.txt", col.names = F, sep = "\t", quote = F)

library(MatrixEQTL)
base.dir = "./"
useModel = modelLINEAR
SNP_file_name = paste(base.dir, "snps.hl.f2.txt", sep="")
meth_file_name = paste(base.dir, "cpg.hl.f2.txt", sep="")
covariates_file_name = paste(base.dir, "covariate.hl.f2.txt", sep="") 
output_file_name = paste(base.dir, "output.hl.f2.txt", sep="") 
pvOutputThreshold = 1e-5
errorCovariance = as.numeric()

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10000;     # read file in pieces of 10,000 rows
snps$LoadFile(SNP_file_name)    # 525 SNPs

meth.qtl = SlicedData$new();
meth.qtl$fileDelimiter = "\t";      # the TAB character
meth.qtl$fileOmitCharacters = "NA"; # denote missing values;
meth.qtl$fileSkipRows = 1;          # one row of column labels
meth.qtl$fileSkipColumns = 1;       # one column of row labels
meth.qtl$fileSliceSize = 100000;   # read file in pieces of 1,000,000 rows
meth.qtl$LoadFile(meth_file_name)   # 27614 CpGs

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 10;        # read file in pieces of 10 rows
cvrt$LoadFile(covariates_file_name)

me = Matrix_eQTL_engine(
  snps = snps,
  gene = meth.qtl,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  # errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = 100,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE) 

plot(me)

output=read.table("./output.hl.f2.txt", sep = "\t", skip = 1)
colnames(output)=c("SNP", "CpG", "beta", "t-stat", "p-value", "FDR")
output$BH=p.adjust(output$`p-value`, method = "BH", n = length(output$`p-value`))
p.adjust=5e-8/27614
output.sig=output[output$BH<p.adjust,] 
nrow(output.sig) # 968 sig mQTL

write.table(output.sig, "sig.meqtl_hl_f2.txt", sep = "\t", quote = F)

# Check number of unique SNPs and CpGs
output.sig$SNP=as.character(output.sig$SNP) 
output.sig$CpG=as.character(output.sig$CpG)

length(unique(output.sig$SNP)) # 335 SNPS
length(unique(output.sig$CpG)) # 75 CpGs, 75/27614=0.002716014

# Functional analysis of SNPs and CpGs in mQTL
library(GenomicFeatures)
library(biomaRt)
library(ChIPpeakAnno)

# extract genomic information from reference genome--------
stickle=makeTxDbFromEnsembl(organism = "Gasterosteus aculeatus")

genes=genes(stickle)

# Convert meQTL to GRanges
library(tidyr)
output.sig=separate(data = output.sig, col = SNP, into = c("snp_chr", "snp_start"), sep = ",")
output.sig=separate(data = output.sig, col = CpG, into = c("cpg_chr", "cpg_start"), sep = ",")

# cis and tran mQTL------
output.sig$cis_trans=ifelse(output.sig$snp_chr==output.sig$cpg_chr & abs(as.numeric(output.sig$snp_start)-as.numeric(output.sig$cpg_start))<1000000, "cis", "trans")
nrow(output.sig[output.sig$cis_trans=="cis",]) # 3 cis
nrow(output.sig[output.sig$cis_trans=="trans",]) # 965 trans
output.sig$site=paste(output.sig$snp_chr, output.sig$snp_start, sep = " ")

# G-test
observed = c(3,965)
expected = c(0.5, 0.5)
library(DescTools)
GTest(x=observed,
      p=expected,
      correct="none") 
# G = 1301.3, X-squared df = 1, p-value < 2.2e-16

# make sig mQTL to GRange
mqtl=output.sig[,1:2]
mqtl$end=mqtl$snp_start
colnames(mqtl)[1:2]=c("chr", "start")
mqtl=unique(mqtl)
mqtl_grange=as(mqtl, "GRanges")

# Test overlap between meQTLs and eQTLs in Ishikawa et al. 2017. Only shared and 10% SW is considered.
eqtl=read.csv("./eqtl_location.csv")
eqtl_grange=as(eqtl, "GRanges")
overlaps_eqtl_mqtl=annotatePeakInBatch(mqtl_grange, AnnotationData = eqtl_grange, output = "overlapping")
overlaps_eqtl_mqtl=data.frame(overlaps_eqtl_mqtl)
overlaps_eqtl_mqtl=overlaps_eqtl_mqtl[!is.na(overlaps_eqtl_mqtl$fromOverlappingOrNearest),]
nrow(unique(overlaps_eqtl_mqtl[,1:3])) # 24 overlap, 335-24=311 not overlap

# how many input SNPs overlap with eqtl regions
library(tidyr)
input_snp=data.frame(snps_pos)
input_snp=separate(data = input_snp, col = snps_pos, into = c("chr", "start"), sep = ",")
input_snp$end=input_snp$start
input_snp_grange=as(input_snp, "GRanges")
overlaps_eqtl_input_snp=annotatePeakInBatch(input_snp_grange, AnnotationData = eqtl_grange, output = "overlapping")
overlaps_eqtl_input_snp=data.frame(overlaps_eqtl_input_snp)
overlaps_eqtl_input_snp=overlaps_eqtl_input_snp[!is.na(overlaps_eqtl_input_snp$fromOverlappingOrNearest),]
nrow(unique(overlaps_eqtl_input_snp[,1:3])) # 39 overlap, 525-39=486 not overlap

# test whether meqtl overlap with eqtl is significant
# G-test
observed = c(24, 311)
expected = c(39/525, 486/525)
library(DescTools)
GTest(x=observed,
      p=expected,
      correct="none") 
# G = 0.034432, X-squared df = 1, p-value = 0.8528

# Genes in meqtl
overlaps.gene_mqtl=annotatePeakInBatch(mqtl_grange, AnnotationData = genes, output = "overlapping")
overlaps.gene_mqtl=data.frame(overlaps.gene_mqtl)
overlaps.gene_mqtl=overlaps.gene_mqtl[!is.na(overlaps.gene_mqtl$fromOverlappingOrNearest),]

# GO analysis for all mQTLs
# Convert Ensembl gene names to real gene names
mart=useDataset("gaculeatus_gene_ensembl", useMart("ensembl"))
gene_mqtl=getBM(filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id", "external_gene_name","go_id", "name_1006", "description"), 
                 values = overlaps.gene_mqtl$feature, 
                 mart = mart)

# Pool genes are genes overlapped with SNPs passed the filter
# Genes in meqtl
pool_gene_pos=read.table("./hl_f2.012.pos")
pool_gene_pos$end=pool_gene_pos$V2
colnames(pool_gene_pos)[1:2]=c("chr", "start")
pool.gene=annotatePeakInBatch(as(pool_gene_pos, "GRanges"), AnnotationData = genes, output = "overlapping")
pool.gene=data.frame(pool.gene)
pool.gene=pool.gene[!is.na(pool.gene$fromOverlappingOrNearest),]

gene_pool=getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "external_gene_name","go_id", "name_1006", "description"), 
                values = pool.gene$feature, 
                mart = mart)

# Remove blank entries
gene_pool <- gene_pool[gene_pool$go_id != '',]
gene_mqtl <- gene_mqtl[gene_mqtl$go_id != '',]

# convert from table format to list format
geneID2GO <- by(gene_pool$go_id,
                gene_pool$ensembl_gene_id,
                function(x) as.character(x))

mqtlGenes=by(gene_mqtl$go_id,
              gene_mqtl$ensembl_gene_id,
              function(x) as.character(x))

correction<-"fdr"
geneNames = names(geneID2GO)

myInterestingGenesNames=names(mqtlGenes)
geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))
names(geneList) <- geneNames

setwd("./")

library(topGO)
ontology=c("MF","BP","CC")

# GO analysis for share-------
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
  write.csv(allRes, paste("hl_mqtl_genes",ontology[i],"csv",sep="."))
} 

## Test if mQTL falls into differential regions between marine and freshwater
## Table comes from Pritchard et al. Regulatory Architecture of Gene Expression Variation in the Threespine Stichleback Gasterosteus aculeatus
## Recompile regions in Table S6 from the original three studies
library(tidyr)
jones=read.csv("./Jones.csv")
jones=separate(data = jones, col = region, into = c("chr", "regions"), sep = ":")
jones=separate(data = jones, col = regions, into = c("start", "end"), sep = "-")
jones$chr=sub("chr", "group", jones$chr)

hohenlohe=read.csv("./hohenlohe.csv")
terekhanova=read.csv("./terekhanova.csv")

divergence_region=rbind(jones, hohenlohe, terekhanova)
divergence_region=unique(divergence_region) # excluding duplicated regions
divergence_region_granges=as(divergence_region, "GRanges")

# Find overlap genes with mQTLs
mqtl_divergence_overlap=annotatePeakInBatch(mqtl_grange, AnnotationData =  divergence_region_granges, output = "overlapping")
mqtl_divergence_overlap=data.frame(mqtl_divergence_overlap)
mqtl_divergence_overlap=mqtl_divergence_overlap[!is.na(mqtl_divergence_overlap$fromOverlappingOrNearest),] 
# These are the overlapping mQTL in divergence regions

# Count how many unique mqtl overlap with divergence region
unique_mqtl_divergence_overlap=unique(mqtl_divergence_overlap[,1:3])
nrow(unique_mqtl_divergence_overlap) # 4 unique mQTL overlap with divergence region
colnames(unique_mqtl_divergence_overlap)[1]="chr"

# Check if there is overrepresentation of cis/trans in divergence-overlapping mQTL
unique_mqtl_divergence_overlap$site=paste(unique_mqtl_divergence_overlap$chr, unique_mqtl_divergence_overlap$start, sep = " ")

unique_mqtl_cis_trans=merge(unique_mqtl_divergence_overlap, output.sig, by="site", all=T)
unique_mqtl_cis_trans=unique_mqtl_cis_trans[!is.na(unique_mqtl_cis_trans$chr),]
unique_mqtl_cis_trans=unique_mqtl_cis_trans[,c(1,14)]
# 4 unique meQTL associated with 11 SNP-CpG associations
nrow(unique_mqtl_cis_trans[unique_mqtl_cis_trans$cis_trans=="cis",]) 
# 0
nrow(unique_mqtl_cis_trans[unique_mqtl_cis_trans$cis_trans=="trans",]) 
# 11

# G-test
observed = c(0, 11)
expected = c(3/968, 965/968)
library(DescTools)
GTest(x=observed,
      p=expected,
      correct="none") 
# G = 0.068288, X-squared df = 1, p-value = 0.7933

## effect size of the 11 meQTLs
effect_size_unique_mqtl_divergence_overlap=merge(unique_mqtl_divergence_overlap, output.sig, by="site", all=T)
effect_size_unique_mqtl_divergence_overlap=effect_size_unique_mqtl_divergence_overlap[!is.na(effect_size_unique_mqtl_divergence_overlap$chr),]
effect_size_unique_mqtl_divergence_overlap=effect_size_unique_mqtl_divergence_overlap[,c(1,5:9)]
summary(effect_size_unique_mqtl_divergence_overlap$beta)

write.table(effect_size_unique_mqtl_divergence_overlap, "effect_size_unique_mqtl_divergence_overlap.txt", quote = F, row.names = F)

# Genes associated with divergence-overlapped mQTLs
unique_mqtl_divergence_overlap_grange=as(unique_mqtl_divergence_overlap, "GRanges")
unique_mqtl_divergence_overlap_gene=annotatePeakInBatch(unique_mqtl_divergence_overlap_grange, AnnotationData =  genes, output = "overlapping")
unique_mqtl_divergence_overlap_gene=data.frame(unique_mqtl_divergence_overlap_gene)
# find out common gene names and description
unique_mqtl_divergence_overlap_gene_name=getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "external_gene_name", "description"), 
                values = unique_mqtl_divergence_overlap_gene$feature, 
                mart = mart)

metadata_unique_mqtl_divergence_overlap_gene=merge(unique_mqtl_divergence_overlap_gene, unique_mqtl_divergence_overlap_gene_name,
                                                   by.x="feature", by.y="ensembl_gene_id", all=T)
metadata_unique_mqtl_divergence_overlap_gene=metadata_unique_mqtl_divergence_overlap_gene[,c(1,16:17)] 
metadata_unique_mqtl_divergence_overlap_gene=unique(metadata_unique_mqtl_divergence_overlap_gene)
colnames(metadata_unique_mqtl_divergence_overlap_gene)=c("Ensembl Gene ID", "Symbol", "Description")
# 4 mQTL overlapped 2 genes
write.csv(metadata_unique_mqtl_divergence_overlap_gene, "./hl_mqtl_divergence.csv", row.names = F, quote = F)

# Check what SNPs underlying DMCs between marine and freshwater ecotypes in parental generation
dmcs_parental=read.table("./dmcs_parental.txt")
dmcs_parental$end=dmcs_parental$V2
colnames(dmcs_parental)[1:2]=c("chr", "start")
dmcs_parental_granges=as(dmcs_parental, "GRanges")
# identify overlaps between dmcs in parental generation and cpgs in significant meqtls
cpg_mqtl=output.sig[,3:4]
cpg_mqtl$end=cpg_mqtl$cpg_start
colnames(cpg_mqtl)[1:2]=c("chr", "start")
cpg_mqtl_granges=as(cpg_mqtl, "GRanges")

overlaps_dmcs_parental_cpg_in_mqtl=annotatePeakInBatch(dmcs_parental_granges, AnnotationData = cpg_mqtl_granges, output = "overlapping")
overlaps_dmcs_parental_cpg_in_mqtl=data.frame(overlaps_dmcs_parental_cpg_in_mqtl)
overlaps_dmcs_parental_cpg_in_mqtl=overlaps_dmcs_parental_cpg_in_mqtl[!is.na(overlaps_dmcs_parental_cpg_in_mqtl$fromOverlappingOrNearest),]

# extract SNP position from sig meqtls using above result as an index
# overlaps_dmcs_parental_cpg_in_mqtl$seqnames=as.character(overlaps_dmcs_parental_cpg_in_mqtl$seqnames)
# output.sig$cpg_start=as.numeric(output.sig$cpg_start)
# snp_underlying_dmcs_in_parental=output.sig[output.sig$cpg_chr %in% overlaps_dmcs_parental_cpg_in_mqtl$seqnames & output.sig$cpg_start %in% overlaps_dmcs_parental_cpg_in_mqtl$start, ]
# No dmc in parental generation is found in meQTL result
# nrow(snp_underlying_dmcs_in_parental) # 0
# nrow(snp_underlying_dmcs_in_parental[snp_underlying_dmcs_in_parental$cis_trans=="cis",]) # 0
# nrow(snp_underlying_dmcs_in_parental[snp_underlying_dmcs_in_parental$cis_trans=="trans",]) # 0

# check genes associated with the 0 SNPs
# snp_underlying_dmcs_in_parental_loc=snp_underlying_dmcs_in_parental[,1:2]
# colnames(snp_underlying_dmcs_in_parental_loc)=c("chr", "start")
# snp_underlying_dmcs_in_parental_loc$end=snp_underlying_dmcs_in_parental_loc$start
# snp_underlying_dmcs_in_parental_loc_granges=as(snp_underlying_dmcs_in_parental_loc, "GRanges")
# genes_in_snp_underlying_dmcs_in_parental=annotatePeakInBatch(snp_underlying_dmcs_in_parental_loc_granges, AnnotationData = genes, output = "overlapping")
# genes_in_snp_underlying_dmcs_in_parental=data.frame(genes_in_snp_underlying_dmcs_in_parental)
# genes_in_snp_underlying_dmcs_in_parental=genes_in_snp_underlying_dmcs_in_parental[!is.na(genes_in_snp_underlying_dmcs_in_parental$fromOverlappingOrNearest),]
# genes_in_snp_underlying_dmcs_in_parental$feature %in% unique_mqtl_divergence_overlap_gene_name$ensembl_gene_id # No genes in the 9 SNPs are in divergence regions of marine and freshwater ecotypes

# gene_names_in_snp_underlying_dmcs_in_parental=getBM(filters = "ensembl_gene_id", 
#                 attributes = c("ensembl_gene_id", "external_gene_name","go_id", "name_1006", "description"), 
#                 values = genes_in_snp_underlying_dmcs_in_parental$feature, 
#                 mart = mart)

# merge all information for the 51 SNPs (locations, names, go id, go terms, etc...)
# special_cpg_associated_snp=merge(gene_names_in_snp_underlying_dmcs_in_parental, genes_in_snp_underlying_dmcs_in_parental, by.x = "ensembl_gene_id", by.y="feature")
# length(unique(special_cpg_associated_snp$ensembl_gene_id)) # 30 genes annotated with the 51 SNPs
# for export
# special_cpg_associated_snp.1=special_cpg_associated_snp[,1:8]
# special_cpg_associated_snp.1=special_cpg_associated_snp.1[,c(6:8, 1:5)]
# special_cpg_associated_snp.1$end=NULL
# colnames(special_cpg_associated_snp.1)[2]="location"

# write.csv(special_cpg_associated_snp.1, "./hl_snps_underlying_special_dmc_full_info.csv", row.names = F)

# Check # of SNPs associated with the special CpG locate within genomic regions of high divergence between ecotypes

# snp_underlying_dmcs_in_parental_overlap_divergence=annotatePeakInBatch(snp_underlying_dmcs_in_parental_loc_granges, AnnotationData = divergence_region_granges, output = "overlapping")
# snp_underlying_dmcs_in_parental_overlap_divergence=data.frame(snp_underlying_dmcs_in_parental_overlap_divergence)
# snp_underlying_dmcs_in_parental_overlap_divergence 
# no SNPs associated with the special CpG locates within genomic regions of high divergence between ecotypes

# what is the gene associated with this special dmc, because only one dmc is associated with sig meQTL, so we just use the first row as an presentation
# special_dmc_loc=overlaps_dmcs_parental_cpg_in_mqtl[1,1:3]
# special_dmc_granges=as(special_dmc_loc, "GRanges")
# gene_special_dmc=annotatePeakInBatch(special_dmc_granges, AnnotationData = genes, output = "overlapping")
# gene_special_dmc=as.data.frame(gene_special_dmc)

# full information of this gene associated with this special dmc
# gene_special_dmc_info=getBM(filters = "ensembl_gene_id", 
#                             attributes = c("ensembl_gene_id", "external_gene_name","go_id", "name_1006", "description"), 
#                             values = gene_special_dmc$feature, 
#                             mart = mart)

# write.csv(snp_underlying_dmcs_in_parental, "./hl_snp_underlying_dmcs_in_parental.csv")
# write.csv(gene_special_dmc_info, "./hl_gene_special_dmc_info.csv")
