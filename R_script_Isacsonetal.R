###############################################################################
###                                                                         ###
###             Code for Isacson et al., investigating                      ###
###              spermborne small RNA in IVF couples                        ###
###              Code written by Signe Isacson, 2024                        ###
###                                                                         ###
###############################################################################


# load all necessary packages


library(seqpac)
library(patchwork)
library(UpSetR)
library(tidyverse)
library(ggpubr)
library(readxl)
library(plotROC)


# ---- preprocessing of files ----


counts1 <- make_counts(input="D:/2023/KIPF/sample_batch1")
pheno1 <- data.frame(row.names=colnames(counts1$counts),
                     sampleID=colnames(counts1$counts))
pac1 <- make_PAC(counts1, pheno1)
pac1 <- PAC_filter(pac1, threshold=1, coverage=50, size=c(16,75), stat=TRUE)
gc()

counts2 <- make_counts(input="D:/2023/KIPF/sample_batch2")
pheno2 <- data.frame(row.names=colnames(counts2$counts),
                     sampleID=colnames(counts2$counts))
pac2 <- make_PAC(counts2, pheno2)
pac2 <- PAC_filter(pac2, threshold=1, coverage=50, size=c(16,75), stat=TRUE)
gc()

pheno_merge <- rbind(pheno1, pheno2)
pheno_merge$sampleID <- sub("_L.*","", pheno_merge$sampleID)
rownames(pheno_merge) <- pheno_merge$sampleID

count_merge <- merge(pac1@Counts, pac2@Counts, by=0, all=TRUE)
rownames(count_merge)<-count_merge$Row.names
count_merge<-count_merge[,-1]
colnames(count_merge) <- sub("_L.*","", colnames(count_merge))
count_merge[is.na(count_merge)]<-0

pac <- make_PAC(count_merge, pheno_merge)


# ---- prepare journal information ----


pheno_patient <- as.data.frame(read_excel("~/2024/KIPF/Masterexcel_2024.xlsx", sheet = "New_R"))
pheno_merge2 <- merge(pac@Pheno, pheno_patient, by.x="sampleID", by.y="couplecode_R", all.x = TRUE)
pheno_merge2 <- pheno_merge2[match(pac@Pheno$sampleID, pheno_merge2$sampleID),]
pac@Pheno<-cbind(pac@Pheno, pheno_merge2)
table(pac@Pheno$flowcell)
hist(pac@Pheno$BMI)
save(pac, file="master_pac.Rdata")
load("master_pac.Rdata")


# filter and clean up pac object ----


pac_f <- PAC_filter(pac, threshold=1, coverage=90, size=c(16,75), stat=TRUE)
pac_f <- PAC_norm(pac_f)
#     Remove samples without oocytes
pac_h <- PAC_filter(pac_f, threshold=1, coverage=90, norm="cpm", stat=TRUE,
                    pheno_target = list("Fertilityrate_group", c("yes", "no")))


# ---- map to references ----


map_reanno(pac_h, 
           output_path="D:/OUT",
           ref_paths=list(genome="D:/Genomes/HUMAN/genomes/Homo_sapiens.GRCh38.dna.toplevel.fa",
                          mitochondria="D:/Genomes/HUMAN/mito/mitoncbi.fa"),
           mismatches=3)

reanno_genome <- make_reanno(reanno_path="D:/OUT", 
                             PAC=pac_h, 
                             mis_fasta_check = TRUE)

pac_h <- add_reanno(reanno_genome,
                    mismatches = 3, merge_pac=pac_h)

map_reanno(pac_h, 
           ref_paths=list(ens="D:/Genomes/HUMAN/ncRNA/Homo_sapiens.GRCh38.ncrna.fa",
                          rRNA="D:/Genomes/HUMAN/rRNA/rnacentral_all_rrna.fa",
                          piRNA="D:/Genomes/HUMAN/piRNA/hsa_piRbase_2.0.fa",
                          mito="D:/Genomes/HUMAN/mito/mitoncbi.fa",
                          tRNA="D:/Genomes/HUMAN/tRNA/trnaensembl.fa"),
           output_path="D:/Genomes/out",
           type="internal", 
           mismatches=1,  
           import="biotype", 
           threads=1)

reanno_biotype <- make_reanno(reanno_path="D:/Genomes/out",
                              PAC=pac_h, 
                              mis_fasta_check = TRUE)

bio_search <- list(ens=c("snRNA", "snoRNA", "rRNA", "miRNA", "tRNA", "Mt_tRNA",
                         "scaRNA", "lncRNA", "Y_RNA", "RNY1", "RN7SL", "RNY4", "RNY3", "RNY5",
                         "7SK", "vault", "SRP"),
                   rRNA=c("URS"),
                   piRNA=c("piR"),
                   mito=c("NC_012920"),
                   tRNA=c("Homo_", "Mt_tRNA", "mt_tRNA"))

pac_h <- add_reanno(reanno_biotype, bio_search=bio_search, 
                  type="biotype", bio_perfect=TRUE,
                  merge_pac = pac_h,
                  mismatches = 1)

pac_h<-simplify_reanno(input=pac_h, 
                       hierarchy = list(rRNA="ens_rRNA|rRNA",
                                        miRNA="ens_miRNA",
                                        mito="mito|tRNA_mt_tRNA|tRNA_Mt_tRNA",
                                        tRNA="ens_tRNA|tRNA_Homo_",
                                        ribonucleoproteins=c("ens_SRP|ens_RN7SL|ens_7SK|ens_vault|ens_snoRNA|ens_scaRNA|ens_snRNA|ens_Y_RNA|ens_RNY4|ens_RNY1|ens_RNY3"),
                                        lncRNA="ens_lncRNA",
                                        miscRNA="ens_misc_RNA",
                                        piRNA="piRNA"),
                       mismatches = 1,
                       bio_name = "Biotypes_1mm_n",
                       merge_pac=TRUE)

table(pac_h@Anno$Biotypes_1mm_n)
save(pac_h, file="pac_annotated.Rdata")


# Supplementary Figure 1


pca <- PAC_pca(pac_h)
pca$graphs$PC1_PC2|pca$graphs$PC1_PC3|pca$graphs$PC2_PC3
sb <- PAC_stackbar(pac_h, anno_target = list("Biotypes_1mm_n"))


# ---- perform DEA ----


dsq_sp <- PAC_deseq(pac_h, model = ~ Spermcount + flowcell , pheno_target = list("Spermcount", c("High", "Low")))
dsq_sm <- PAC_deseq(pac_h, model = ~ more_or_equal_to_5_million_motile_sperm + flowcell, pheno_target = list("more_or_equal_to_5_million_motile_sperm", 
                                                                                                             c("yes", "no")))
pac_h@Pheno$fertilityrate<-as.numeric(pac_h@Pheno$fertilityrate)
pac_h@Pheno$Fertilityrate_group2 <- ifelse(pac_h@Pheno$fertilityrate>=70, "High", "Low")
dsq_f <- PAC_deseq(pac_h, model = ~ Fertilityrate_group2 + flowcell, pheno_target = list("Fertilityrate_group2", c("High", "Low")))
dsq_utr <- PAC_deseq(pac_h, model = ~ Mod_ut_num_group + flowcell, pheno_target = list("Mod_ut_num_group", c("High", "Low")))


# ---- Figure 2 ----


## Figure 2 a

dsq_objs <- list(d1=dsq_sp,
                 d3=dsq_f,
                 d4=dsq_utr,
                 d2=dsq_sm)


for(i in names(dsq_objs)){
  res <- dsq_objs[[i]]$result
  lbl <- names(res)[3]
  res$FC<-res[,3]
  res$DE <- ifelse(res$FC>1 & res$padj<0.1, "UP", "")
  res$DE <- ifelse(res$FC<(-1) & res$padj<0.1, "DOWN", res$DE)
  res <- res[,c("FC", "DE", "Biotypes_1mm_n", "padj")]
  res <- na.omit(res)
  res$col <- ifelse(res$DE %in% "UP", res$Biotypes_1mm_n, "")
  res$col <- ifelse(res$DE %in% "DOWN", res$Biotypes_1mm_n, res$col)
  plot<-ggplot(res, aes(FC, -log10(padj), color=col)) + 
    geom_point(size=3) +
    scale_color_manual(values=c("darkgrey", "#d8c2c2", "#8e6e53",
                                "#d16e6e",  "#f79256","#fbd1a2",
                                "#7dcfb6","#00b2ca", "#1d4e89",
                                "#F0F415","#F7F6C5","#F0F0F0",
                                "#F0F0F1", "#F0F0F2","#F0F0F3"))+
    theme_classic()+
    ggtitle(label=(lbl))+
    xlim(-4,4)+
    geom_hline(yintercept = 1, linetype="dashed")+
    geom_vline(xintercept = c(-1,1), linetype="dashed")+
    theme(text=element_text(size=12), legend.position="bottom")
  dsq_objs[[i]]$plot<-plot
  dsq_objs[[i]]$res <- res
}


dsq_objs$d1$plot + dsq_objs$d3$plot + dsq_objs$d4$plot


## Figure 2 b


pac_h <- PAC_summary(pac_h, norm="cpm")

pac_scu <- PAC_filter(pac_h, anno_target = rownames(dsq_objs$d1$res[dsq_objs$d1$res$DE %in% "UP",]))
pie_upsc <- PAC_pie(pac_scu, anno_target = list("Biotypes_1mm_n"), summary = "all", labels = "all", norm="cpm")

pac_scd <- PAC_filter(pac_h, anno_target = rownames(dsq_objs$d1$res[dsq_objs$d1$res$DE %in% "DOWN",]))
pie_downsc <- PAC_pie(pac_scd, anno_target = list("Biotypes_1mm_n"), summary = "all", labels = "all", norm="cpm")

pac_eu <- PAC_filter(pac_h, anno_target = rownames(dsq_objs$d4$res[dsq_objs$d4$res$DE %in% "UP",]))
pie_upeq <- PAC_pie(pac_eu, anno_target = list("Biotypes_1mm_n"), summary = "all", labels = "all", norm="cpm")

pac_ed <- PAC_filter(pac_h, anno_target = rownames(dsq_objs$d4$res[dsq_objs$d4$res$DE %in% "DOWN",]))
pie_downeq <- PAC_pie(pac_ed, anno_target = list("Biotypes_1mm_n"), summary = "all", labels = "all", norm="cpm")

pac_fu <- PAC_filter(pac_h, anno_target = rownames(dsq_objs$d3$res[dsq_objs$d3$res$DE %in% "DOWN",]))
pie_downfert <- PAC_pie(pac_fu, anno_target = list("Biotypes_1mm_n"),  summary = "all", labels = "all", norm="cpm")

# Supplementary Figure 2


all_rows <- unique(c(rownames(dsq_objs$d1$res[dsq_objs$d1$res$DE %in% "UP",]), 
                     rownames(dsq_objs$d1$res[dsq_objs$d1$res$DE %in% "DOWN",]),
                     rownames(dsq_objs$d2$res[dsq_objs$d2$res$DE %in% "UP",]),
                     rownames(dsq_objs$d2$res[dsq_objs$d2$res$DE %in% "DOWN",])))

binary_matrix <- data.frame(
  row_name = all_rows,
  sc_up = as.integer(all_rows %in% rownames(dsq_objs$d1$res[dsq_objs$d1$res$DE %in% "UP",])),
  sc_dwn = as.integer(all_rows %in% rownames(dsq_objs$d1$res[dsq_objs$d1$res$DE %in% "DOWN",])),
  sm_up = as.integer(all_rows %in% rownames(dsq_objs$d2$res[dsq_objs$d2$res$DE %in% "UP",])),
  sm_dwn = as.integer(all_rows %in% rownames(dsq_objs$d2$res[dsq_objs$d2$res$DE %in% "DOWN",])))


upset_data <- binary_matrix[, -1]
upset(upset_data, line.size = 1.4, nsets=20, text.scale = c(1.3, 1.3, 1, 1, 2, 1.5), order.by="freq")
dsq_objs$d1$plot + dsq_objs$d2$plot


# Figure 3


## Figure 3 a


df <- PAC_filter(pac_scu, anno_target = list("Biotypes_1mm_n", "mito"))
map <- PAC_mapper(df, 
                  ref="D:/Genomes/HUMAN/mito/mitoncbi_genelevel.fa",
                  mismatches = 1)

test <- tRNA_class(df, map)
df@Anno[!rownames(df@Anno) %in% rownames(test@Anno),]

tst <- data.frame(test@norm$cpm, seq=test@Anno$tRNA_ref)
tst <- aggregate(tst[,1:70], by=list(tst$seq), FUN=sum)
plot_df <- gather(tst, key="sample", value="cpm", -Group.1)
plot_df$group <- rep(df@Pheno$Spermcount, each=nrow(tst))

plot_df_subset <- plot_df[plot_df$Group.1 %in% c("MT-TH-His","MT-TS1-Ser1", "MT-TQ-Glu"),]

ggplot(plot_df_subset, aes(x=Group.1, y=cpm, fill=group))+
  geom_boxplot(size=1)+
  geom_point(position = position_dodge(width = .75), size=2)+
  theme_classic()+
  scale_fill_manual(values=c("#eba468", "#85aa80"))+
  stat_compare_means()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  facet_wrap(~Group.1, scales = "free", ncol=1)


# Supplementary Figure 3 a


ggplot(plot_df, aes(x=Group.1, y=cpm, fill=group))+
  geom_boxplot(size=1)+
  geom_point(position = position_dodge(width = .75), size=2)+
  theme_classic()+
  scale_fill_manual(values=c("#eba468", "#85aa80"))+
  stat_compare_means()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  facet_wrap(~Group.1, scales = "free", ncol=3)


## Figure 3 b


corrplot <- (tst[tst$Group.1 %in% "MT-TS1-Ser1",2:71])
corrplot[2,] <- (df@Pheno$sperm_conc_milj_per_ml)
rownames(corrplot) <- c("MT_TS1_Ser1", "spermc")
corrplot <- t(corrplot)

ggplot(corrplot, aes(x=MT_TS1_Ser1, y=spermc)) + 
  geom_point(size=5)+
  geom_smooth(method=lm, color="black")+
  theme(axis.line=element_line(size=3), 
        axis.text.x = element_text(size=46),
        axis.text.y = element_text(size=16))+
  theme_classic() +
  stat_cor(method = "pearson", label.x = 3, label.y = 30)

## Figure 3 c


df_roc <- pac_h@Pheno[,c("sperm_conc_milj_per_ml", "Mod_ut_num_group")]
df_roc$spermc_bin <- ifelse(df_roc$sperm_conc_milj_per_ml<=16, 0, 1)

df_roc <- cbind(df_roc, corrplot[,1])

colnames(df_roc) <- c("sperm_conc_milj_per_ml", "Mod_ut_num_group", "spermc_bin", "TS1")
roc1<-ggplot(df_roc, aes(d = spermc_bin, m = TS1)) + geom_roc() + style_roc()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30))
calc_auc(roc1)


## Figure 3 d


df <- PAC_filter(pac_scd, anno_target = list("Biotypes_1mm_n", "ribonucleoproteins"))
df@Anno

map <- PAC_mapper(df, 
                  ref="D:/Genomes/HUMAN/yrna_snrna.fa",
                  mismatches = 1)

test <- tRNA_class(df, map)

df@Anno[!rownames(df@Anno) %in% rownames(test@Anno),]

tst <- data.frame(test@norm$cpm, seq=test@Anno$tRNA_ref)
tst <- aggregate(tst[,1:70], by=list(tst$seq), FUN=sum)
plot_df <- gather(tst, key="sample", value="cpm", -Group.1)
plot_df$group <- rep(df@Pheno$Spermcount, each=nrow(tst))

plot_df_sub <- plot_df[plot_df$Group.1 %in% c("RNY1", "RNY3", "RNY4"),]

ggplot(plot_df_sub, aes(x=Group.1, y=cpm, fill=group))+
  geom_boxplot(size=1)+
  geom_point(position = position_dodge(width = .75), size=2)+
  theme_classic()+
  scale_fill_manual(values=c("#eba468", "#85aa80"))+
  stat_compare_means()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  facet_wrap(~Group.1, scales = "free", ncol=1)


# Supplementary figure 3 b

ggplot(plot_df, aes(x=Group.1, y=cpm, fill=group))+
  geom_boxplot(size=1)+
  geom_point(position = position_dodge(width = .75), size=2)+
  theme_classic()+
  scale_fill_manual(values=c("#eba468", "#85aa80"))+
  stat_compare_means()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  facet_wrap(~Group.1, scales = "free", ncol=3)

## Figure 3 e

corrplot <- (tst[tst$Group.1 %in% "RNY4",2:71])
corrplot[2,] <- (df@Pheno$sperm_conc_milj_per_ml)
rownames(corrplot) <- c("RNY4", "spermc")
corrplot <- t(corrplot)

ggplot(corrplot, aes(x=RNY4, y=spermc)) + 
  geom_point()+
  geom_smooth(method=lm, color="black")+
  labs(title="Ser1 vs sperm conc.")+
  theme_classic() +
  stat_cor(method = "pearson", label.x = 3, label.y = 30)

## Figure 3 f


df_roc <- pac_h@Pheno[,c("sperm_conc_milj_per_ml", "Mod_ut_num_group")]
df_roc$spermc_bin <- ifelse(df_roc$sperm_conc_milj_per_ml<=16, 1, 0)

df_roc <- cbind(df_roc, corrplot[,1])

colnames(df_roc) <- c("sperm_conc_milj_per_ml", "Mod_ut_num_group", "spermc_bin", "RNY4")
roc1<-ggplot(df_roc, aes(d = spermc_bin, m = RNY4)) + geom_roc() + style_roc()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30))
calc_auc(roc1)


# Figure 4


samples_to_use <- rownames(pac_h@Pheno[pac_h@Pheno$Mod_ut_num_group %in% c("High", "Low"),])


# Fig 4 A


df <- PAC_filter(pac_eu, anno_target = list("Biotypes_1mm_n", "miRNA"))
map <- PAC_mapper(df, 
                  ref="D:/Genomes/HUMAN/miRNA/miRbase_2024.fa",
                  mismatches = 1)

test <- tRNA_class(df, map)
tst <- data.frame(test@norm$cpm, seq=test@Anno$type)
tst <- aggregate(tst[,1:70], by=list(tst$seq), FUN=sum)
plot_df <- gather(tst, key="sample", value="cpm", -Group.1)
plot_df$group <- rep(df@Pheno$Mod_ut_num_group, each=nrow(tst))
plot_df$alc <- rep(df@Pheno$Mod_ut_num_group, each=nrow(tst))
plot_df <- plot_df[plot_df$group %in% c("High", "Low"),]

plot_df_subset <- plot_df[plot_df$Group.1 %in% c("let-hsa-let-7g MI0000433", 
                                                 "mir-320b;mir-hsa-mir-320a MI0000542",
                                                 "mir-hsa-mir-30d MI0000255"),]

ggplot(plot_df_subset, aes(x=Group.1, y=cpm, fill=group))+
  geom_boxplot(size=1)+
  geom_point(position = position_dodge(width = .75), size=2)+
  theme_classic()+
  scale_fill_manual(values=c("#eba468", "#85aa80"))+
  stat_compare_means()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  facet_wrap(~Group.1, scales = "free", ncol=1)

## Supplementary Figure 4 a

ggplot(plot_df, aes(x=Group.1, y=cpm, fill=group))+
  geom_boxplot(size=1)+
  geom_point(position = position_dodge(width = .75), size=2)+
  theme_classic()+
  scale_fill_manual(values=c("#eba468", "#85aa80"))+
  stat_compare_means()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  facet_wrap(~Group.1, scales = "free", ncol=5)

## Figure 4 b


corrplot <- plot_df[plot_df$Group.1 %in% "let-hsa-let-7g MI0000433",]
corrplot <- corrplot[corrplot$sample %in% samples_to_use,]
pheno2 <- pac_h@Pheno[pac_h@Pheno$sampleID %in% samples_to_use,]
corrplot2 <- cbind(corrplot, pheno2)
corrplot2$utr <- corrplot2$`Good_quality_embryos_retrieved_ per_fertilized (new modified utilization rate)`
corrplot2 <- corrplot2[,c(1:5,87:90)]
corrplot2$utr <- as.numeric(corrplot2$utr)

ggplot(corrplot2, aes(x=cpm, y=utr)) + 
  geom_point()+
  geom_smooth(method=lm, color="black")+
  labs(title="Ser1 vs sperm conc.")+
  theme_classic() +
  stat_cor(method = "pearson", label.x = 3, label.y = 30)


##Figure 4 c


df_roc <- pac_h@Pheno[,c("sperm_conc_milj_per_ml", "Mod_ut_num_group")]
df_roc$mod_bin <- ifelse(df_roc$Mod_ut_num_group=="High", 1, 2)
df_roc$mod_bin <- ifelse(df_roc$Mod_ut_num_group=="Low", 0, df_roc$mod_bin)
df_roc <- df_roc[!df_roc$mod_bin==2,]
df_roc <- cbind(df_roc, corrplot2$cpm)

colnames(df_roc) <- c("sperm_conc_milj_per_ml", "Mod_ut_num_group", "mod_bin", "RNY4")
roc1<-ggplot(df_roc, aes(d = mod_bin, m = RNY4)) + geom_roc() + style_roc()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30))
calc_auc(roc1)


## Supplementary figure 4 c


corrplot <- plot_df[plot_df$Group.1 %in% "mir-hsa-mir-30d MI0000255",]
corrplot <- corrplot[corrplot$sample %in% samples_to_use,]
pheno2 <- pac_h@Pheno[pac_h@Pheno$sampleID %in% samples_to_use,]
corrplot2 <- cbind(corrplot, pheno2)
corrplot2$utr <- corrplot2$`Good_quality_embryos_retrieved_ per_fertilized (new modified utilization rate)`
corrplot2 <- corrplot2[,c(1:5,87:90)]
corrplot2$utr <- as.numeric(corrplot2$utr)

ggplot(corrplot2, aes(x=cpm, y=utr)) + 
  geom_point()+
  geom_smooth(method=lm, color="black")+
  labs(title="Ser1 vs sperm conc.")+
  theme_classic() +
  stat_cor(method = "pearson", label.x = 3, label.y = 30)

df_roc <- pac_h@Pheno[,c("sperm_conc_milj_per_ml", "Mod_ut_num_group")]
df_roc$mod_bin <- ifelse(df_roc$Mod_ut_num_group=="High", 1, 2)
df_roc$mod_bin <- ifelse(df_roc$Mod_ut_num_group=="Low", 0, df_roc$mod_bin)
df_roc <- df_roc[!df_roc$mod_bin==2,]
df_roc <- cbind(df_roc, corrplot2$cpm)

colnames(df_roc) <- c("sperm_conc_milj_per_ml", "Mod_ut_num_group", "mod_bin", "RNY4")
roc1<-ggplot(df_roc, aes(d = mod_bin, m = RNY4)) + geom_roc() + style_roc()+
  theme(axis.line=element_line(size=1), 
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30))
calc_auc(roc1)


## Figure 4 d


# Supplementary Figure 4 b
#all rsRNA

## Figure 4 e


## Figure 4 f





# Supplementary Figure 5


#read in output from GO term analysis


GOplot<-clipr::read_clip_tbl()

colnames(GOplot)<-c("GO", "ref", "up1", "exp", "over", "foldenr", "rawp", "FDR")
GOplot$foldenr<-as.numeric(GOplot$foldenr)
GOplot$FDR<-as.numeric(GOplot$FDR)
GOplot <- GOplot[order(GOplot$FDR, decreasing = F),]
GOplot <- GOplot[order(GOplot$foldenr, decreasing = T),]
GOplot$GO <- factor(GOplot$GO, levels=GOplot$GO)

GOplot <- GOplot[1:20,]

ggplot(GOplot, aes(x=foldenr, y=GO, color=FDR, size=foldenr))+
  geom_point()+
  theme_classic()+
  scale_y_discrete(limits=rev)
