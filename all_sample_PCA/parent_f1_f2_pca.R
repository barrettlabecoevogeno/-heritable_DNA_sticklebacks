library(methylKit)
setwd("/Volumes/Epiguru/Thesis chapters/final version of each chapter with analysis/sticklbeack/submission to Genetics/PCA")
# Metadata for all samples
info=read.csv("./stickleback_metadata.csv")

directory="/Volumes/Epiguru/Thesis chapters/raw_data/stickleback/stickle_R/combine/"

filenames=list.files(path=directory, full.names=TRUE)
filenames=as.list(filenames)
names=list.files(path=directory)
names=gsub(".bismark.cov", "", names)
names=as.list(names)

my.methRaw=methRead(location = filenames,
                    sample.id = names,
                    assembly = "stickleback",
                    pipeline = 'bismarkCoverage',
                    context = "CpG",
                    treatment = c(rep(1,11), # This is a fake treatment lable, no real meaning at all
                                  rep(0,83)),
                    mincov = 10)

filtered.my.methRaw <- filterByCoverage(my.methRaw, 
                                        lo.count=10,
                                        lo.perc = NULL,
                                        hi.count = NULL,
                                        hi.perc = 99.9)

# normalize read coverages between samples to avoid bias introduced by systematically more sequnenced sameples

normalized.myobj=normalizeCoverage(filtered.my.methRaw, method="median")

# merging samples for DMCs (step should be saved)------

meth.all=unite(normalized.myobj, destrand = F, mc.cores = 2)
save(meth.all, file = "meth.RData")

setwd("/Volumes/Epiguru/Thesis chapters/final version of each chapter with analysis/sticklbeack/submission to Genetics/PCA")
load("./meth.RData")
info=read.csv("./stickleback_metadata.csv")
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
CT_maf= read.csv(file="./PCA/result_pat_C_T.txt", sep="\t", header=FALSE)
CT_maf = cbind(CT_maf,V3=rep(CT_maf$V2))
write.table(CT_maf,file="./PCA/result_pat_C_T_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_CT_bed<- bed_to_granges(file="./PCA/result_pat_C_T_mafSTARTEND.bed")

# load GA SNPs and duplicate column and save as bedfile, result_pat_G_A_maf.bed is created in step 5.2
AG_maf= read.csv(file="./PCA/result_pat_G_A.txt", sep="\t", header=FALSE)
AG_maf = cbind(AG_maf,V3=rep(AG_maf$V2))
write.table(AG_maf,file="./PCA/result_pat_G_A_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_GA_bed<- bed_to_granges(file="./PCA/result_pat_G_A_mafSTARTEND.bed")

# interesect bedfile and unite-file
#### create overlap --> these positions are corrected for CT SNPs
unite_norm_10x_GRanges <- as(meth.all, "GRanges")
Overlap_CT=unite_norm_10x_GRanges[countOverlaps(unite_norm_10x_GRanges, blacklist_CT_bed) == 0L]
meth.1=makeMethylDB(meth.all,"methylBaseDB")
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

# PCA----
library(ggbiplot)

meth.count=regionCounts(unite_norm_10x_CT_GA_nosex_DB, regions = as(unite_norm_10x_CT_GA_nosex_DB, "GRanges")) # removed 54044-52940=1104 sites
perc.meth.all=percMethylation(meth.count)

pca=prcomp(t(perc.meth.all), center = T)
summary(pca)

# Add population info to pca matrix
df=as.data.frame(pca$x)
library(dplyr)
info$Label=as.character(info$Label)
info=info %>%
  slice(match(rownames(df), Label))

# PCA in all individuals-------
# calculate variance explained by each PC
var.all=as.data.frame(pca$x)
var.all=apply(var.all, 2, var)
var.all=var.all/sum(var.all)
var.all=var.all[2:3]
# PC2        PC3 
# 0.05146615 0.03178719

axes=as.data.frame(pca$x)
# only use PC2 and PC3
axes=axes[,c(2,3)]
axes.all=axes
axes.all$Label=rownames(axes.all)

axes.all=merge(axes.all, info, by="Label")

axes.all$Generation=factor(axes.all$Generation, levels = c("P", "F1", "F2"), ordered = T)
axes.all$Line=axes.all$Origin

g1=ggplot(axes.all, aes(PC2, PC3))
g2=g1+geom_point(aes(color=Generation, fill=Generation, shape=Line), size=3)+
  scale_color_manual(values = c("#999999", "#009E73", "#D55E00"))+
  scale_fill_manual(values = c("#999999", "#009E73", "#D55E00"))+
  scale_shape_manual(values = c(21:22, 24))+
  labs(x="PC2", y="PC3")
g3=g2 + theme_bw()+theme(panel.grid=element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank())
g3


# Extract top CpGs with highest or lowest loadings on PC1
loading_name=paste(meth.p$chr,meth.p$start, sep = ",")
loadings=data.frame(pca$rotation, index=1:length(loading_name))
rownames(loadings)=loading_name
pc1_loading=loadings[,colnames(loadings) %in% c("PC1", "index")]

pc1_min=head(pc1_loading[order(pc1_loading$PC1),],10)
pc1_max=tail(pc1_loading[order(pc1_loading$PC1),],10)

pc1_loading_rep=rbind(pc1_min, pc1_max)
pc1_loading_rep$site=rownames(pc1_loading_rep)
pc1_loading_rep$index=NULL

library(tidyr)
pc1_loading_rep=pc1_loading_rep%>%
  separate(site, c("chr", "location"), ",")

pc1_loading_rep=pc1_loading_rep[,c(2,3,1)]

write.csv(pc1_loading_rep, "pc1_loading.csv", row.names = F, quote = F)

## PCA on each generation
# extract axes 
axes=as.data.frame(pca$x)
# only use PC2 and PC3
axes=axes[,c(2,3)]

## Parent generation
info.p=info[info$Generation=="P",]
# calculate %var explained
var.p=as.data.frame(pca$x)
var.p=var.p[rownames(var.p) %in% info.p$Label, ]
var.p=apply(var.p, 2, var)
var.p=var.p/sum(var.p)
var.p=var.p[2:3]
# PC2         PC3 
# 0.005090635 0.063732006 

axes.p=axes[rownames(axes) %in% info.p$Label,]
axes.p$Label=rownames(axes.p)

axes.p=merge(axes.p, info.p, by="Label")
axes.p$Group=factor(c("BI", "HL", "BI", "KL", "BI", "KL", "BI", "BI", "KL", "BI", "HL"))
axes.p$Line=axes.p$Group
axes.p$Habitat=factor(axes.p$Habitat)

g1.p=ggplot(axes.p, aes(PC2, PC3))
g2.p=g1.p + geom_point(aes(shape=Line, fill=Habitat, color=Habitat), size=3)+
  scale_shape_manual(values = c(21:22, 24))+
  labs(x="PC2", y="PC3")
g3.p=g2.p + theme_bw()+theme(panel.grid=element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())
g3.p

## F1 generation
info.f1=info[info$Generation=="F1",]
# calculate %var explained
var.f1=as.data.frame(pca$x)
var.f1=var.f1[rownames(var.f1) %in% info.f1$Label, ]
var.f1=apply(var.f1, 2, var)
var.f1=var.f1/sum(var.f1)
var.f1=var.f1[2:3]
# PC2        PC3 
# 0.03087134 0.02641961 

axes.f1=axes[rownames(axes) %in% info.f1$Label,]
axes.f1$Label=rownames(axes.f1)

axes.f1=merge(axes.f1, info.f1, by="Label")

axes.f1$Line=axes.f1$Origin

g1.f1=ggplot(axes.f1, aes(PC2, PC3))
g2.f1=g1.f1 + geom_point(aes(shape=Line, color=Family, fill=Family), size=2)+
  scale_shape_manual(values = c(22,24))+
  labs(x="PC2", y="PC3")
g3.f1=g2.f1 + theme_bw()+theme(panel.grid=element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank())
g3.f1

# F2 generation
info.f2=info[info$Generation=="F2",]
# calculate %var explained
var.f2=as.data.frame(pca$x)
var.f2=var.f2[rownames(var.f2) %in% info.f2$Label, ]
var.f2=apply(var.f2, 2, var)
var.f2=var.f2/sum(var.f2)
var.f2=var.f2[2:3]
# PC2         PC3 
# 0.008608972 0.026590484 

axes.f2=axes[rownames(axes) %in% info.f2$Label,]
axes.f2$Label=rownames(axes.f2)

axes.f2=merge(axes.f2, info.f2, by="Label")
axes.f2$Line=factor(axes.f2$Origin)

g1.f2=ggplot(axes.f2, aes(PC2, PC3))
g2.f2=g1.f2 + geom_point(aes(shape=Line), size=2)+
  scale_shape_manual(values = c(22,24))+
  labs(x="PC2", y="PC3")
g3.f2=g2.f2 + theme_bw()+theme(panel.grid=element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank())
g3.f2

library(cowplot)
pca_plot=plot_grid(g3, g3.p, g3.f1, g3.f2,
                   labels = c("a", "b", "c", "d"),
                   align = "vh")
pca_plot




## Supplemental material figure S1------
library(ggbiplot)
# By generation and line
Line=info$Origin
Line=factor(Line, levels = c("BI", "HL", "KL"), ordered = T)
Generation=info$Generation
Generation=factor(Generation, levels = c("P", "F1", "F2"), ordered = T)
Label=info$Label
Label[c(1:83,88:94)]=NA

s1a=ggbiplot(pca,
              obs.scale = 1, 
              var.scale = 1,
              varname.size = 0,
              var.axes = F)

s1a=s1a + geom_point(aes(color=Generation, fill=Generation, shape=Line), size=3)+
  scale_color_manual(values = c("#999999", "#009E73", "#D55E00"))+
  scale_fill_manual(values = c("#999999", "#009E73", "#D55E00"))+
  scale_shape_manual(values = c(21:22, 24))+
  geom_text(aes(label=Label), hjust=-0.15, size=3)
s1a=s1a + theme_bw()+theme(panel.grid=element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank())
s1a
# By sequencing lane
Lane=info$Lane
Lane=factor(Lane, levels = c("1", "2", "3", "4"), ordered = T)
s2a=ggbiplot(pca,
             obs.scale = 1, 
             var.scale = 1,
             varname.size = 0,
             var.axes = F)

s2a=s2a + geom_point(aes(color=Generation, fill=Generation, shape=Lane), size=3)+
  scale_color_manual(values = c("#999999", "#009E73", "#D55E00"))+
  scale_fill_manual(values = c("#999999", "#009E73", "#D55E00"))+
  scale_shape_manual(values = c(21:24))+
  geom_text(aes(label=Label), hjust=-0.15, size=3)
s2a=s2a + theme_bw()+theme(panel.grid=element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.background = element_blank())
s2a

library(cowplot)
Fig.S1=plot_grid(s1a, s2a,
                 labels = c("a", "b"))
Fig.S1

# track the methylation levels of the 891 DMCs identified in parental generation across F1 and F2 generation----------
CpG.all.meth=as.data.frame(perc.meth.all)
p_dmcs_loc=read.table("dmcs_parental.txt")
colnames(p_dmcs_loc)=c("chr", "start", "end")
p_dmcs_loc_in_all_meth=findOverlaps(as(p_dmcs_loc, "GRanges"), as(unite_norm_10x_CT_GA_nosex_DB, "GRanges"))
p_dmcs_loc_in_all_meth=data.frame(p_dmcs_loc_in_all_meth)

CpG.p.dmcs.in.all.meth=CpG.all.meth[p_dmcs_loc_in_all_meth$subjectHits,]
axes.all$Label==colnames(CpG.p.dmcs.in.all.meth)

library(pheatmap)
sampleinfo=data.frame(Line=axes.all$Origin,
                      Generation=axes.all$Generation)
rownames(sampleinfo)=colnames(CpG.p.dmcs.in.all.meth)
# levels(sampleinfo$Habitat)[2]="Freshwater"
ann_colors=list(Generation=c(P="#00BFC4", F1="#F8766D", F2="#009E73"),
                Line=c(BI="#000000", HL="#CCCCCC", KL="#999999"))
p.heatmap=pheatmap(CpG.p.dmcs.in.all.meth, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   show_colnames = TRUE,
                   cluster_cols = TRUE,
                   border_color = NA,
                   scale = "row",
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   clustering_method = "complete",
                   annotation_col = sampleinfo,
                   annotation_colors = ann_colors,
                   annotation_names_col = F)

# comparing the changes of methylation levels of the 891 DMCs across generations
index_p_BI=as.numeric(rownames(axes.all[axes.all$Generation=="P" & axes.all$Origin=="BI",]))
index_p_KL=as.numeric(rownames(axes.all[axes.all$Generation=="P" & axes.all$Origin=="KL",]))
index_p_HL=as.numeric(rownames(axes.all[axes.all$Generation=="P" & axes.all$Origin=="HL",]))
index_f1_KL=as.numeric(rownames(axes.all[axes.all$Generation=="F1" & axes.all$Origin=="KL",]))
index_f1_HL=as.numeric(rownames(axes.all[axes.all$Generation=="F1" & axes.all$Origin=="HL",]))
index_f2_KL=as.numeric(rownames(axes.all[axes.all$Generation=="F2" & axes.all$Origin=="KL",]))
index_f2_HL=as.numeric(rownames(axes.all[axes.all$Generation=="F2" & axes.all$Origin=="HL",]))

meth_p_BI=CpG.p.dmcs.in.all.meth[,index_p_BI]
meth_p_KL=CpG.p.dmcs.in.all.meth[,index_p_KL]
meth_p_HL=CpG.p.dmcs.in.all.meth[,index_p_HL]
meth_f1_KL=CpG.p.dmcs.in.all.meth[,index_f1_KL]
meth_f1_HL=CpG.p.dmcs.in.all.meth[,index_f1_HL]
meth_f2_KL=CpG.p.dmcs.in.all.meth[,index_f2_KL]
meth_f2_HL=CpG.p.dmcs.in.all.meth[,index_f2_HL]

meth_p_BI$mean=apply(meth_p_BI, 1, mean)
meth_p_KL$mean=apply(meth_p_KL, 1, mean)
meth_p_HL$mean=apply(meth_p_HL, 1, mean)
meth_f1_KL$mean=apply(meth_f1_KL, 1, mean)
meth_f1_HL$mean=apply(meth_f1_HL, 1, mean)
meth_f2_KL$mean=apply(meth_f2_KL, 1, mean)
meth_f2_HL$mean=apply(meth_f2_HL, 1, mean)

meth_f1_KL$between=ifelse(meth_f1_KL$mean>meth_p_BI$mean & meth_f1_KL$mean<meth_p_KL$mean, "yes",
                          ifelse(meth_f1_KL$mean<meth_p_BI$mean & meth_f1_KL$mean>meth_p_KL$mean, "yes", "no"))

nrow(meth_f1_KL[meth_f1_KL$between=="yes",]) # 694
nrow(meth_f1_KL[meth_f1_KL$between=="no",]) # 197

library(DescTools)
expected=c(0.5,1-0.5)
observed=c(694,197)

GTest(x=observed,
      p=expected,
      correct="none") # G = 293.76, X-squared df = 1, p-value < 2.2e-16

# sex bias
meth_f1_KL$diff=abs(meth_f1_KL$mean-meth_p_BI$mean)-abs(meth_f1_KL$mean-meth_p_KL$mean)
meth_f1_KL$bias=ifelse(meth_f1_KL$diff>0 & meth_f1_KL$between=="yes", "father", 
                       ifelse(meth_f1_KL$diff<0 & meth_f1_KL$between=="yes", "mother", "no_bias"))

nrow(meth_f1_KL[meth_f1_KL$bias=="father",]) # 357
nrow(meth_f1_KL[meth_f1_KL$bias=="mother",]) # 337
nrow(meth_f1_KL[meth_f1_KL$bias=="no_bias",]) # 197

library(DescTools)
expected=c(0.5,1-0.5)
observed=c(357,337)

GTest(x=observed,
      p=expected,
      correct="none") # G = 0.57645, X-squared df = 1, p-value = 0.4477
# no bias in sex in KL F1

meth_f1_HL$between=ifelse(meth_f1_HL$mean>meth_p_BI$mean & meth_f1_HL$mean<meth_p_HL$mean, "yes",
                          ifelse(meth_f1_HL$mean<meth_p_BI$mean & meth_f1_HL$mean>meth_p_HL$mean, "yes", "no"))

nrow(meth_f1_HL[meth_f1_HL$between=="yes",]) # 610
nrow(meth_f1_HL[meth_f1_HL$between=="no",]) # 281

library(DescTools)
expected=c(0.5,1-0.5)
observed=c(610,281)

GTest(x=observed,
      p=expected,
      correct="none") # G = 124.41, X-squared df = 1, p-value < 2.2e-16

meth_f1_HL$diff=abs(meth_f1_HL$mean-meth_p_BI$mean)-abs(meth_f1_HL$mean-meth_p_HL$mean)
meth_f1_HL$bias=ifelse(meth_f1_HL$diff>0 & meth_f1_HL$between=="yes", "father", 
                       ifelse(meth_f1_HL$diff<0 & meth_f1_HL$between=="yes", "mother", "no_bias"))

nrow(meth_f1_HL[meth_f1_HL$bias=="father",]) # 280
nrow(meth_f1_HL[meth_f1_HL$bias=="mother",]) # 330
nrow(meth_f1_HL[meth_f1_HL$bias=="no_bias",]) # 281

library(DescTools)
expected=c(0.5,1-0.5)
observed=c(280,330)

GTest(x=observed,
      p=expected,
      correct="none") # G = 4.103, X-squared df = 1, p-value = 0.04281 
# more sites bias toward mother in HL F1

### F2 vs. parental
meth_f2_KL$between=ifelse(meth_f2_KL$mean>meth_p_BI$mean & meth_f2_KL$mean<meth_p_KL$mean, "yes",
                          ifelse(meth_f2_KL$mean<meth_p_BI$mean & meth_f2_KL$mean>meth_p_KL$mean, "yes", "no"))

nrow(meth_f2_KL[meth_f2_KL$between=="yes",]) # 487
nrow(meth_f2_KL[meth_f2_KL$between=="no",]) # 404

meth_f2_HL$between=ifelse(meth_f2_HL$mean>meth_p_BI$mean & meth_f2_HL$mean<meth_p_HL$mean, "yes",
                          ifelse(meth_f2_HL$mean<meth_p_BI$mean & meth_f2_HL$mean>meth_p_HL$mean, "yes", "no"))

nrow(meth_f2_HL[meth_f2_HL$between=="yes",]) # 578
nrow(meth_f2_HL[meth_f2_HL$between=="no",]) # 313

### F2 vs. F1
meth_f2_KL$betweenf1=ifelse(meth_f2_KL$mean>meth_f1_KL$KL1F & meth_f2_KL$mean<meth_f1_KL$KL1M, "yes",
                          ifelse(meth_f2_KL$mean>meth_f1_KL$KL1F & meth_f2_KL$mean<meth_f1_KL$KL1M, "yes", "no"))

nrow(meth_f2_KL[meth_f2_KL$betweenf1=="yes",]) # 219
nrow(meth_f2_KL[meth_f2_KL$betweenf1=="no",]) # 672

meth_f2_HL$betweenf1=ifelse(meth_f2_HL$mean>meth_f1_HL$HL2M3M & meth_f2_HL$mean<meth_f1_HL$HL2F, "yes",
                          ifelse(meth_f2_HL$mean<meth_f1_HL$HL2M3M & meth_f2_HL$mean>meth_f1_HL$HL2F, "yes", "no"))

nrow(meth_f2_HL[meth_f2_HL$betweenf1=="yes",]) # 399
nrow(meth_f2_HL[meth_f2_HL$betweenf1=="no",]) # 492

# sex bias analysis
meth_f2_KL$diff=abs(meth_f2_KL$mean-meth_f1_KL$KL1F)-abs(meth_f2_KL$mean-meth_f1_KL$KL1M)
meth_f2_KL$bias=ifelse(meth_f2_KL$diff>0 & meth_f2_KL$betweenf1=="yes", "father", 
                       ifelse(meth_f2_KL$diff<0 & meth_f2_KL$betweenf1=="yes", "mother", "no_bias"))

nrow(meth_f2_KL[meth_f2_KL$bias=="father",]) # 95
nrow(meth_f2_KL[meth_f2_KL$bias=="mother",]) # 124
nrow(meth_f2_KL[meth_f2_KL$bias=="no_bias",]) # 672

library(DescTools)
expected=c(0.5,1-0.5)
observed=c(95,124)

GTest(x=observed,
      p=expected,
      correct="none") # G = 3.8515, X-squared df = 1, p-value = 0.0497
# more sites bias toward mother in KL F2

meth_f2_HL$diff=abs(meth_f2_HL$mean-meth_f1_HL$HL2F)-abs(meth_f2_HL$mean-meth_f1_HL$HL2M3M)
meth_f2_HL$bias=ifelse(meth_f2_HL$diff>0 & meth_f2_HL$betweenf1=="yes", "father", 
                       ifelse(meth_f2_HL$diff<0 & meth_f2_HL$betweenf1=="yes", "mother", "no_bias"))

nrow(meth_f2_HL[meth_f2_HL$bias=="father",]) # 202
nrow(meth_f2_HL[meth_f2_HL$bias=="mother",]) # 197
nrow(meth_f2_HL[meth_f2_HL$bias=="no_bias",]) # 492

library(DescTools)
expected=c(0.5,1-0.5)
observed=c(202,197)

GTest(x=observed,
      p=expected,
      correct="none") # G = 0.062658, X-squared df = 1, p-value = 0.8023
# no bias in sex in HL F2




# Following codes are not used, because PCAs were performed seperated on each generation-------
# PCA on parental samples--------
info.p=info[info$Generation=="P",]
meth.p=reorganize(unite_norm_10x_CT_GA_nosex_DB,
                  sample.ids = as.character(info.p$Label),
                  treatment = c(rep(0,5),
                                rep(1,6))) # First make a fake treatment label, and adjust in next step

meth.p@sample.ids # Check order of samples

meth.p@treatment=c(1,0,1,0,1,0,1,1,0,1,0) # Marine is 1, Fresh is 0
meth.p.1=as(meth.p, "GRanges")

meth.count.p=regionCounts(meth.p, regions = meth.p.1)
perc.meth.p=percMethylation(meth.count.p)

pca.p=prcomp(t(perc.meth.p), center = T)
summary(pca.p)

# Add population info to pca matrix
df.p=as.data.frame(pca.p$x)
library(dplyr)
info.p=info.p %>%
  slice(match(rownames(df.p), Label))

info.p$Group=factor(c("BI", "HL", "BI", "KL", "BI", "KL", "BI", "BI", "KL", "BI", "HL"))
Strain.P=info.p$Group
levels(Strain.P)=c("BI", "HL", "KL")
Habitat=factor(info.p$Habitat)

g1.p=ggbiplot(pca.p,
              obs.scale = 1, 
              var.scale = 1,
              varname.size = 0,
              var.axes = F)

g2.p=g1.p + geom_point(aes(shape=Strain.P, fill=Habitat, color=Habitat), size=3)+
  scale_shape_manual(values = c(21:22, 24))
g3.p=g2.p + theme_bw()+theme(panel.grid=element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())
g3.p

# PCA on F1 samples------
info.f1=info[info$Generation=="F1",]
meth.f1=reorganize(unite_norm_10x_CT_GA_nosex_DB,
                   sample.ids = as.character(info.f1$Label),
                   treatment = c(rep(0,9),
                                 rep(1,10))) # First make a fake treatment label, and adjust in next step

meth.f1@sample.ids # Check order of samples

meth.f1@treatment=c(rep(0,7), rep(1,12)) # HL is 0, KL is 1
meth.f1.1=as(meth.f1, "GRanges")

meth.count.f1=regionCounts(meth.f1, regions = meth.f1.1)
perc.meth.f1=percMethylation(meth.count.f1)

pca.f1=prcomp(t(perc.meth.f1), center = T)
summary(pca.f1)

# Add population info to pca matrix
df.f1=as.data.frame(pca.f1$x)
library(dplyr)
info.f1=info.f1 %>%
  slice(match(rownames(df.f1), Label))

Population.f1=factor(info.f1$Group)
Strain.F1=info.f1$Origin
Family=info.f1$Family

g1.f1=ggbiplot(pca.f1, 
               obs.scale = 1, 
               var.scale = 1,
               varname.size = 0,
               var.axes = F)

g2.f1=g1.f1 + geom_point(aes(shape=Strain.F1, color=Family, fill=Family), size=2)+
  scale_shape_manual(values = c(22,24))
g3.f1=g2.f1 + theme_bw()+theme(panel.grid=element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank())
g3.f1

# PCA on F2-------
info.f2=info[info$Generation=="F2",]
meth.f2=reorganize(unite_norm_10x_CT_GA_nosex_DB,
                  sample.ids = as.character(info.f2$Label),
                  treatment = c(rep(0,32),
                                rep(1,32))) # First make a fake treatment label, and adjust in next step

meth.f2@sample.ids # Check order of samples

meth.f2@treatment=c(rep(0,28),
                    rep(1,36)) # HL is 0, KL is 1
meth.f2.1=as(meth.f2, "GRanges")

meth.count.f2=regionCounts(meth.f2, regions = meth.f2.1)
perc.meth.f2=percMethylation(meth.count.f2)

pca.f2=prcomp(t(perc.meth.f2), center = T)
summary(pca.f2)

# Add population info to pca matrix
df.f2=as.data.frame(pca.f2$x)
library(dplyr)
info.f2=info.f2 %>%
  slice(match(rownames(df.f2), Label))

Strain.F2=factor(info.f2$Origin)
Lane=factor(info.f2$Lane)

g1.f2=ggbiplot(pca.f2, 
               obs.scale = 1, 
               var.scale = 1,
               varname.size = 0,
               var.axes = F)

g2.f2=g1.f2 + geom_point(aes(color=Lane, shape=Strain.F2, fill=Lane), size=2)+
  scale_shape_manual(values = c(22,24))
g3.f2=g2.f2 + theme_bw()+theme(panel.grid=element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank())
g3.f2
