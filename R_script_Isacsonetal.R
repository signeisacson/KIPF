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


# ---- perform DEA ----

dsq_sp_bmi <- PAC_deseq(pac_h, model = ~ Spermcount + flowcell + BMI_bin , pheno_target = list("Spermcount", c("High", "Low")))
dsq_sp <- PAC_deseq(pac_h, model = ~ Spermcount + flowcell , pheno_target = list("Spermcount", c("High", "Low")))


dsq_sp <- PAC_deseq(pac_h, model = ~ Spermcount + flowcell , pheno_target = list("Spermcount", c("High", "Low")))
dsq_sm <- PAC_deseq(pac_h, model = ~ more_or_equal_to_5_million_motile_sperm + flowcell, pheno_target = list("more_or_equal_to_5_million_motile_sperm", 
                                                                                                             c("yes", "no")))
pac_h@Pheno$fertilityrate<-as.numeric(pac_h@Pheno$fertilityrate)
pac_h@Pheno$Fertilityrate_group2 <- ifelse(pac_h@Pheno$fertilityrate>=70, "High", "Low")
dsq_f <- PAC_deseq(pac_h, model = ~ Fertilityrate_group2 + flowcell, pheno_target = list("Fertilityrate_group2", c("High", "Low")))
dsq_utr <- PAC_deseq(pac_h, model = ~ Mod_ut_num_group + flowcell, pheno_target = list("Mod_ut_num_group", c("High", "Low")))

dsq_utr_allB <- PAC_deseq(pac_h, model = ~ Mod_ut_num_group + flowcell + BMI, pheno_target = list("Mod_ut_num_group", c("High", "Low")))
dsq_utr_bmi <- PAC_deseq(pac_h, model = ~ Mod_ut_num_group + flowcell + BMI_bin, pheno_target = list("Mod_ut_num_group", c("High", "Low")))
dsq_utr <- PAC_deseq(pac_h, model = ~ Mod_ut_num_group + flowcell, pheno_target = list("Mod_ut_num_group", c("High", "Low")))


dsq_utr$plots$volcano | dsq_utr_bmi$plots$volcano
