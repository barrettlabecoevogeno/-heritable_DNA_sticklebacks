load("./F1 vs F2 same strain/meth.RData")
info=read.csv("./stickleback_metadata.csv")
library(methylKit)
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

# The C/T and A/G location is produced before purging high LD sites, using parental samples only.
# load CT SNPs and duplicate column and save as bedfile, result_pat_C_T_maf.bed
CT_maf= read.csv(file="./F1 vs F2 same strain/result_pat_C_T.txt", sep="\t", header=FALSE)
CT_maf = cbind(CT_maf,V3=rep(CT_maf$V2))
write.table(CT_maf,file="./F1 vs F2 same strain/result_pat_C_T_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_CT_bed<- bed_to_granges(file="./F1 vs F2 same strain/result_pat_C_T_mafSTARTEND.bed")

# load GA SNPs and duplicate column and save as bedfile, result_pat_G_A_maf.bed is created in step 5.2
AG_maf= read.csv(file="./F1 vs F2 same strain/result_pat_G_A.txt", sep="\t", header=FALSE)
AG_maf = cbind(AG_maf,V3=rep(AG_maf$V2))
write.table(AG_maf,file="./F1 vs F2 same strain/result_pat_G_A_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_GA_bed<- bed_to_granges(file="./F1 vs F2 same strain/result_pat_G_A_mafSTARTEND.bed")

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

# Calculate Euclidian distance matrix using the tiling windows
# get methylation percentage for all CpG sites
tiles=tileMethylCounts(unite_norm_10x_CT_GA_nosex_DB, win.size = 1000, step.size = 1000) # 7840 windows
meth.regions=regionCounts(unite_norm_10x_CT_GA_nosex_DB, regions = as(tiles, "GRanges"))
perc.meth.regions=percMethylation(meth.regions)

# compile index for each hybrid lines in F1 and F2
info$Label=as.character(info$Label)

index_hl_p=info[info$Generation=="P" & info$Origin=="HL",]$Label
index_kl_p=info[info$Generation=="P" & info$Origin=="KL",]$Label
index_bi_p=info[info$Generation=="P" & info$Origin=="BI",]$Label

index_hl_f1=info[info$Generation=="F1" & info$Origin=="HL",]$Label
index_hl_f2=info[info$Generation=="F2" & info$Origin=="HL",]$Label

index_kl_f1=info[info$Generation=="F1" & info$Origin=="KL" & info$Family=="KL_F1_1",]$Label
index_kl_f2=info[info$Generation=="F2" & info$Origin=="KL",]$Label

# extract methylation matrix for each hybrid line
perc.hl.p=data.frame(perc.meth.regions[,colnames(perc.meth.regions) %in% index_hl_p])
perc.kl.p=data.frame(perc.meth.regions[,colnames(perc.meth.regions) %in% index_kl_p])
perc.bi.p=data.frame(perc.meth.regions[,colnames(perc.meth.regions) %in% index_bi_p])

perc.hl.f1=data.frame(perc.meth.regions[,colnames(perc.meth.regions) %in% index_hl_f1])
perc.hl.f2=data.frame(perc.meth.regions[,colnames(perc.meth.regions) %in% index_hl_f2])

perc.kl.f1=data.frame(perc.meth.regions[,colnames(perc.meth.regions) %in% index_kl_f1])
perc.kl.f2=data.frame(perc.meth.regions[,colnames(perc.meth.regions) %in% index_kl_f2])


# SD calculation
sd.hl.p=data.frame(sd=apply(perc.hl.p, 1, sd)/100)
sd.kl.p=data.frame(sd=apply(perc.kl.p, 1, sd)/100)
sd.bi.p=data.frame(sd=apply(perc.bi.p, 1, sd)/100)

sd.hl.f1=data.frame(sd=apply(perc.hl.f1, 1, sd)/100)
sd.hl.f1$sample="HL_F1"
sd.hl.f2=data.frame(sd=apply(perc.hl.f2, 1, sd)/100)
sd.hl.f2$sample="HL_F2"
sd.kl.f1=data.frame(sd=apply(perc.kl.f1, 1, sd)/100)
sd.kl.f1$sample="KL_F1"
sd.kl.f2=data.frame(sd=apply(perc.kl.f2, 1, sd)/100)
sd.kl.f2$sample="KL_F2"

sd=rbind(sd.hl.f1, sd.hl.f2, sd.kl.f1, sd.kl.f2)

wilcox.test(sd.hl.f1$sd, sd.hl.f2$sd) # W = 29446000, p-value = 5.613e-06, increased variation in f2
wilcox.test(sd.kl.f1$sd, sd.kl.f2$sd) # W = 30814000, p-value = 0.7752, decreased variation in f2

t.test(sd.kl.f2$sd, sd.kl.p$sd)


# mean.plot=ggplot(mean, aes(x=sample, y=mean/100, fill=sample))+
#  geom_boxplot()+
#  theme_bw()+
#  geom_boxplot()+
#  theme(legend.position = "none",
#        panel.grid=element_blank(),
#        axis.text.x=element_text(size = 10, color = "black"),
#        axis.text.y=element_text(size = 10, color = "black"),
#        axis.title.y = element_text(size = 12, color = "black"))+
#  labs(x="", y="Mean methylation levels of individual fish")

# mean.plot
# Plot a violin plot
# library(tidyr)
# hl.f1=gather(perc.hl.f1, sample, meth, HL2F:HL5M)
# hl.f1$id="HL_F1"
# hl.f2=gather(perc.hl.f2, sample, meth, HL1:HL8)
# hl.f2$id="HL_F2"

# kl.f1=gather(perc.kl.f1, sample, meth, KL1F:KL7M)
# kl.f1$id="KL_F1"
# kl.f2=gather(perc.kl.f2, sample, meth, KL1:KL9)
# kl.f2$id="KL_F2"

# violin=rbind(hl.f1, hl.f2, kl.f1, kl.f2)
# violin$id=factor(violin$id, levels = c("HL_F1", "HL_F2", "KL_F1", "KL_F2"), ordered = T)
# library(ggplot2)
# p=ggplot(violin, aes(x=id, y=meth/100))+
#  geom_violin(aes(fill=id))+
#  theme_bw()+
#  theme(legend.position = "none",
#        panel.grid=element_blank(),
#        axis.text.x=element_text(size = 10, color = "black"),
#        axis.text.y=element_text(size = 10, color = "black"),
#        axis.title.y = element_text(size = 12, color = "black"))+
#  labs(x="", y="Methylation levels")
# p

# summary(mean.hl.f1, mean.hl.f2, mean.kl.f1, mean.kl.f2)

# hl.f1.sd=apply(perc.hl.f1,1,sd)
# hl.f2.sd=apply(perc.hl.f2,1,sd)

# kl.f1.sd=apply(perc.kl.f1,1,sd)
# kl.f2.sd=apply(perc.kl.f2,1,sd)

library(ggplot2)
# sd=data.frame(sd=c(hl.f1.sd, hl.f2.sd, kl.f1.sd, kl.f2.sd),
#               id=c(rep("HL_F1", 52940), rep("HL_F2", 52940), rep("KL_F1", 52940), rep("KL_F2", 52940)))
sd.plot=ggplot(sd, aes(x=sample, y=sd))+
  geom_boxplot()+
  scale_y_continuous(limits = c(0, 0.4))+
  scale_fill_brewer()+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size = 12, color = "black"),
        axis.text.y=element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))+
  labs(x="", y="Standard deviation of methylation levels in 1kb windows")
sd.plot


# Cumulative curve
p1=ggplot(violin, aes(x=meth/100, color=id))+
  stat_ecdf(geom = "step")+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size = 10, color = "black"),
        axis.text.y=element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))+
  labs(x="Methylation levels", y="Fraction of data")

p1

ks.test(violin[violin$id=="HL_F1",]$meth, violin[violin$id=="HL_F2",]$meth)
ks.test(violin[violin$id=="KL_F1",]$meth, violin[violin$id=="KL_F2",]$meth)
ks.test(violin[violin$id=="HL_F1",]$meth, violin[violin$id=="KL_F1",]$meth)
ks.test(violin[violin$id=="HL_F2",]$meth, violin[violin$id=="KL_F2",]$meth)

library(cowplot)
p2=plot_grid(p, p0, p1,
             labels = c("a", "b", "c"),
             nrow = 1,
             align = "vh")
p2

# Alternatives to plot cumulative distribution curves, do no use.
plot(ecdf(violin[violin$id=="HL_F1",]$meth),
     col="green")
lines(ecdf(violin[violin$id=="HL_F2",]$meth),
      col="blue")
lines(ecdf(violin[violin$id=="KL_F1",]$meth),
      col="red")
lines(ecdf(violin[violin$id=="KL_F2",]$meth),
      col="orange")

