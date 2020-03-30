library(methylKit)
# Load methylated sites covered in all samples in three generations
setwd("/Volumes/Epiguru/Thesis chapters/final version of each chapter with analysis/sticklbeack/submission to Genetics/p_ecotypes/")
load("./meth.RData")
# Subset F1 and F2 samples
p_subset=read.csv("./stickleback_metadata.csv")
p_subset=p_subset[p_subset$Generation=="P", ]
p_subset$Label=as.character(p_subset$Label)
p_subset$Lane=as.factor(p_subset$Lane)
p_subset$Generation=as.factor(p_subset$Generation)
# levels(f1f2_subset$Generation)=c("1", "0", NA)
p_subset=p_subset[,c(1,2,6,7)]

# Order is important, adjust order of sample id based on orginal sample id in meth
ind=meth.all@sample.ids[84:94]

library(dplyr)
p_subset=p_subset %>%
  slice(match(ind, Label))

# Subset methylation levels
meth=reorganize(meth.all,
                sample.ids = as.character(p_subset[,1]),
                treatment = as.numeric(p_subset[,3])) # HL is 1, KL is 2

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
CT_maf= read.csv(file="./result_pat_C_T.txt", sep="\t", header=FALSE)
CT_maf = cbind(CT_maf,V3=rep(CT_maf$V2))
write.table(CT_maf,file="./result_pat_C_T_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_CT_bed<- bed_to_granges(file="./result_pat_C_T_mafSTARTEND.bed")

# load GA SNPs and duplicate column and save as bedfile, result_pat_G_A_maf.bed is created in step 5.2
AG_maf= read.csv(file="./result_pat_G_A.txt", sep="\t", header=FALSE)
AG_maf = cbind(AG_maf,V3=rep(AG_maf$V2))
write.table(AG_maf,file="./result_pat_G_A_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_GA_bed<- bed_to_granges(file="./result_pat_G_A_mafSTARTEND.bed")

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
unite_norm_10x_CT_GA_nosex_DB@treatment # marine is 3, freshwater is 2

# Find DMCs between the marine and freshwater ecotypes
myDiff.marine_vs_fresh=calculateDiffMeth(unite_norm_10x_CT_GA_nosex_DB,
                            covariates = p_subset[,3,drop=F], # Origin (sampling site) as a covariate
                            mc.cores = 2)

myDiff.sig.marine_vs_fresh=getMethylDiff(myDiff.marine_vs_fresh, difference=15, qvalue=0.01)
nrow(myDiff.sig.marine_vs_fresh) # 891

# number of hyper vs. hypo DMCs
myDiff.sig.marine_vs_fresh.hyper=getMethylDiff(myDiff.marine_vs_fresh, difference=15, qvalue=0.01,type = "hyper")
nrow(myDiff.sig.marine_vs_fresh.hyper) # 430
myDiff.sig.marine_vs_fresh.hypo=getMethylDiff(myDiff.marine_vs_fresh, difference=15, qvalue=0.01,type = "hypo")
nrow(myDiff.sig.marine_vs_fresh.hypo) # 461

# test if hyper vs. hypo DMCs is biased
library(DescTools)
expected=c(0.5,1-0.5)
observed=c(430/(430+461),461/(430+461))

GTest(x=observed,
      p=expected,
      correct="none")

# G = 0.0012108, X-squared df = 1, p-value = 0.9722

# For ploting purpose
# get methylation percentage for all CpG sites
DMCs=regionCounts(meth, regions = as(myDiff.sig.marine_vs_fresh, "GRanges"))
perc.DMCs=percMethylation(DMCs)

# heatmap based on DMCs
library(pheatmap)
sampleinfo=data.frame(Line=p_subset$Origin,
  Habitat=p_subset$Habitat)
rownames(sampleinfo)=colnames(perc.DMCs)
# levels(sampleinfo$Habitat)[2]="Freshwater"
ann_colors=list(Habitat=c(Marine="#00BFC4", Fresh="#F8766D"),
                Line=c(BI="#000000", HL="#CCCCCC", KL="#999999"))
pheatmap(perc.DMCs, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         show_colnames = FALSE,
         cluster_cols=TRUE,
         border_color = NA,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = sampleinfo,
         annotation_colors = ann_colors,
         annotation_names_col = F)

# Export coordinates of DMCs in parental generation
dmcs_loc=getData(DMCs)
dmcs_loc$chr=as.character(dmcs_loc$chr)
write.table(dmcs_loc[,1:3], "./dmcs_parental.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# Export all filtered sites to build null distribution
all_filtered_cpgs=getData(unite_norm_10x_CT_GA_nosex_DB)
all_filtered_cpgs$chr=as.character(all_filtered_cpgs$chr)
write.table(all_filtered_cpgs[,1:2], "/Volumes/Epiguru/Thesis chapters/final version of each chapter with analysis/sticklbeack/submission to Genetics/all_filtered_CpGs.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# identify genes and GO terms associated with DMCs
library(GenomicFeatures)
library(biomaRt)
library(ChIPpeakAnno)

stickle=makeTxDbFromEnsembl(organism = "Gasterosteus aculeatus")
genes=genes(stickle)

dmcs_loc_granges=as(dmcs_loc, "GRanges")
gene_in_DMCs=annotatePeakInBatch(dmcs_loc_granges, AnnotationData = genes, output = "overlapping")
gene_in_DMCs=data.frame(gene_in_DMCs)
gene_in_DMCs=gene_in_DMCs[!is.na(gene_in_DMCs$fromOverlappingOrNearest),]
length(unique(gene_in_DMCs$feature)) # 228 genes

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
DMCs_gene=getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006",  "description"),  
                values = unique(gene_in_DMCs$feature), 
                mart = mart)

# Merge information of DMCs (ensembl gene id, symbol, description)
DMCs_full_info=merge(DMCs_gene, gene_in_DMCs, by.x="ensembl_gene_id", by.y="feature")
DMCs_full_info=DMCs_full_info[,c(1:8)]

# Representatives of GO terms in methylation divergence between marine and freshwater
library(GO.db)
immune.offspring=GOBPOFFSPRING[["GO:0002376"]]
immune.go.id=c(immune.offspring, "GO:0002376")

metabolic.offspring=GOBPOFFSPRING[["GO:0008152"]]
metabolic.go.id=c(metabolic.offspring, "GO:0008152")

catalytic.offspring=GOBPOFFSPRING[["GO:0003824"]]
catalytic.go.id=c(catalytic.offspring, "GO:0003824")

DMCs_immune_gene=DMCs_full_info[DMCs_full_info$go_id %in% immune.go.id,]
DMCs_metabolic_gene=DMCs_full_info[DMCs_full_info$go_id %in% metabolic.go.id,]
DMCs_catalytic_gene=DMCs_full_info[DMCs_full_info$go_id %in% catalytic.go.id,]


# Check if dmcs in this study overlapped with dmcs in Smith et al. 2015. Mol. Biol. Evol.
dmcs_smith=read.csv("./dmcs_smith.csv") # DMCs are extracted from suppl. material Table S5 from Smith et al. 2015. Mol. Biol. Evol.
dmcs_smith$ID=as.character(dmcs_smith$ID)
overlap_dmcs=DMCs_gene[DMCs_gene$ensembl_gene_id %in% dmcs_smith$ID,]
# Check if dmcs in this study overlapped with DMGs in Artemov et al. 2017. Mol. Biol. Evol.
degs_artemov=read.csv("./DEGs_Artemov.csv") # DMCs are extracted from suppl. material Table S5C from Smith et al. 2015. Mol. Biol. Evol.
degs_artemov$ID=as.character(degs_artemov$ID)
degs_artemov=degs_artemov[degs_artemov$FDR<0.05,]
overlap_degs=DMCs_gene[DMCs_gene$ensembl_gene_id %in% degs_artemov$ID,]

# Remove blank entries
gene_pool <- gene_pool[gene_pool$go_id != '',]
gene_dmcs <- DMCs_gene[DMCs_gene$go_id != '',]

# convert from table format to list format
geneID2GO <- by(gene_pool$go_id,
                gene_pool$ensembl_gene_id,
                function(x) as.character(x))

dmcGenes=by(gene_dmcs$go_id,
             gene_dmcs$ensembl_gene_id,
             function(x) as.character(x))

correction<-"fdr"
geneNames = names(geneID2GO)

myInterestingGenesNames=names(dmcGenes)
geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))
names(geneList) <- geneNames

setwd("~/Desktop/")

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
  write.csv(allRes, paste("parental_dmcs_genes",ontology[i],"csv",sep="."))
} 
