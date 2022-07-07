## 050121 ##

# Hypothesis: 12-13 Di-HOME-treated macrophages have reduced global chromatin accessibility than vehicle.

library(lme4)
library(data.table)
library(dplyr)
library(ggpubr)
library(limma)
library(edgeR)
library(ggvenn)
library(pheatmap)

## ATAC-seq analysis ====
# 1) Check the fragment lengths for each subject. ====

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/Fragment_Lengths")

plot_frag_length <- function(filename, name){
  
  frag_length <- read.table(filename, header = FALSE, stringsAsFactors = FALSE)
  frag_length_1500max <- frag_length[which(frag_length[,2] == 50):which(frag_length[,2] == 1500),]
  colnames(frag_length_1500max) <- c(paste0(name, "_counts"), "Length")
  return(frag_length_1500max)
  
}

frag_length_ta1_1500 <- plot_frag_length("SampleTA1_fragment_length_count_mapq10.txt", "TA1")
frag_length_ta2_1500 <- plot_frag_length("SampleTA2_fragment_length_count_mapq10.txt", "TA2")
frag_length_ta3_1500 <- plot_frag_length("SampleTA3_fragment_length_count_mapq10.txt", "TA3")
frag_length_ta4_1500 <- plot_frag_length("SampleTA4_fragment_length_count_mapq10.txt", "TA4")
frag_length_tb1_1500 <- plot_frag_length("SampleTB1_fragment_length_count_mapq10.txt", "TB1")
frag_length_tb2_1500 <- plot_frag_length("SampleTB2_fragment_length_count_mapq10.txt", "TB2")
frag_length_tb3_1500 <- plot_frag_length("SampleTB3_fragment_length_count_mapq10.txt", "TB3")
frag_length_tb4_1500 <- plot_frag_length("SampleTB4_fragment_length_count_mapq10.txt", "TB4")
frag_length_tc1_1500 <- plot_frag_length("SampleTC1_fragment_length_count_mapq10.txt", "TC1")
frag_length_tc2_1500 <- plot_frag_length("SampleTC2_fragment_length_count_mapq10.txt", "TC2")
frag_length_tc3_1500 <- plot_frag_length("SampleTC3_fragment_length_count_mapq10.txt", "TC3")
frag_length_tc4_1500 <- plot_frag_length("SampleTC4_fragment_length_count_mapq10.txt", "TC4")

frag_length_list <- list(frag_length_ta1_1500, frag_length_ta2_1500, frag_length_ta3_1500, frag_length_ta4_1500, frag_length_tb1_1500, frag_length_tb2_1500, frag_length_tb3_1500, frag_length_tb4_1500, frag_length_tc1_1500, frag_length_tc2_1500, frag_length_tc3_1500, frag_length_tc4_1500)

frag_length_all <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Length"), frag_length_list)

plot(as.numeric(frag_length_all$TA1_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TA1 (Vehicle)")
plot(as.numeric(frag_length_all$TA2_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TA2 (12,13-diHOME)")
plot(as.numeric(frag_length_all$TA3_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TA3 (LPS)")
plot(as.numeric(frag_length_all$TA4_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TA4 (Vehicle)")
plot(as.numeric(frag_length_all$TB1_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TB1 (12,13-diHOME)")
plot(as.numeric(frag_length_all$TB2_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TB2 (LPS)")
plot(as.numeric(frag_length_all$TB3_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TB3 (Vehicle)")
plot(as.numeric(frag_length_all$TB4_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TB4 (12,13-diHOME)")
plot(as.numeric(frag_length_all$TC1_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TC1 (LPS)")
plot(as.numeric(frag_length_all$TC2_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TC2 (Vehicle)")
plot(as.numeric(frag_length_all$TC3_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TC3 (12,13-diHOME)")
plot(as.numeric(frag_length_all$TC4_counts) ~ as.numeric(frag_length_all$Length), pch = 16, col = "Blue", xlab = "Fragment length (bp)", main = "TC4 (LPS)")

# Yes, there are peaks of read counts in the observed locations, suggestive of good quality ATAC-seq data

# 2) Let's analyze the peak counts for the consensus peaks ====

# Load in the data and remove the non-autosomal consensus peaks

all_consensus_peaks <- read.table("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/All_peaks_countMatrix_vehicle_dihome.txt", header = TRUE) # 135963
non_autosomal <- grepl("chrY", all_consensus_peaks[,1]) | grepl("chrX", all_consensus_peaks[,1])

all_consensus_peaks_auto <- all_consensus_peaks[!non_autosomal,]
all_consensus_peaks_annotations <- all_consensus_peaks_auto[,1:6]
all_consensus_peaks_data <- all_consensus_peaks_auto[,7:ncol(all_consensus_peaks_auto)]
rename_col <- gsub(".sorted.nomt.marked.nodup.rmmulti10.bam", "", colnames(all_consensus_peaks_data))
#rename_col2 <- gsub("Sample", "", rename_col)
colnames(all_consensus_peaks_data) <- rename_col
rownames(all_consensus_peaks_data) <- all_consensus_peaks_annotations[,1]

treatment_info <- as.data.frame(cbind(c("TA1", "TA2", "TA3", "TA4", "TB1", "TB2", "TB3", "TB4", "TC1", "TC2", "TC3", "TC4"), c(rep(c("A", "B", "C", "D"), each = 3)), c(rep(c("vehicle", "dihome", "lps"), 4))))
colnames(treatment_info) <- c("ID", "Donor", "Treatment")
treatment_info_vehicle_dihome <- treatment_info[!treatment_info$Treatment == "lps",]

# Let's print out the cpm-adjusted file that includes the samples treated with LPS

cpm_raw_all <- cpm(all_consensus_peaks_data, log = FALSE)
write.table(cpm_raw_all, "/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/GEO/GEO_Human_ATAC_12.txt", quote = FALSE, sep = "\t", )

# Will only analyze vehicle and 12,13-diHOME samples (excluding the LPS-treated samples)

all_consensus_peaks_data_vehicle_dihome <- all_consensus_peaks_data[,!treatment_info$Treatment == "lps"]

# Convert to CPM and keep peaks that had CPM > 1 in at least 4/12 samples (25%)

cpm_raw <- cpm(all_consensus_peaks_data_vehicle_dihome, log=FALSE)
lcpm_raw <- cpm(all_consensus_peaks_data_vehicle_dihome, log=TRUE)

# There seem to be similarly less differentially accessible regions in macs treated with dihome and lps.

cpm_raw_noted <- cpm_raw>1
sample_peaks <- colSums(cpm_raw_noted)
treatment_info_vehicle_dihome$peaks <- sample_peaks 

treatment_info_vehicle_dihome$Treatment <- factor(treatment_info_vehicle_dihome$Treatment, levels = c("vehicle", "dihome"))

comparisons_broad_nums <- list(c('vehicle', 'dihome'))

broad_nums_treatment <- ggplot(data = treatment_info_vehicle_dihome, aes(x = as.factor(Treatment), y = as.numeric(peaks), fill = as.factor(Treatment))) + 
  geom_boxplot() + xlab("") + ylab("Number of Regions") + ggtitle("") + 
  theme(axis.text.x = element_text(size = 10, vjust = -2.5), 
        axis.title.y=element_text(size = 10), 
        axis.text.y = element_text(size = 8, hjust = 1), 
        axis.title=element_text(size = 10), 
        plot.title = element_text(size = 12, hjust = 0.5), 
        legend.position = "none",
        panel.background = element_rect(fill = NA, color = "black")) + 
  scale_x_discrete(labels=c("vehicle" = "Vehicle", "dihome" = "12,13-diHOME, 50uM")) + scale_fill_manual(values=c("#CCCCCC", "#E69F00")) + geom_point() +
  scale_y_continuous(limits = c(111300, 121650)) + stat_compare_means(comparisons = comparisons_broad_nums, method = "t.test", method.args = list(paired = TRUE, alternative = "two.sided")) 


pdf("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/Figures/boxplot_sites_treatment_050921.pdf", width = 4, height = 4)
broad_nums_treatment
dev.off()

# Let's compare the variances of the two samples (vehicle and 12,13-diHOME treated with KS test)

var.test(peaks ~ Treatment, treatment_info_vehicle_dihome, alternative = "two.sided") # P = 0.002294

# 3) Perform PCA ====

## PCA of 132,340 peaks

sum.PC <- prcomp(t(cpm_raw), scale=FALSE, center=TRUE)
sumsum <- summary(sum.PC)
prop.var<-sumsum$importance
prop.var[1:3,1:5] 

##To see if covariates are correlated with a PC (looking at PC1-5)

Donor <- as.factor(treatment_info_vehicle_dihome$Donor)
Treatment <- as.factor(treatment_info_vehicle_dihome$Treatment)

apply(sum.PC$x, 2, function(col)anova(lm(as.numeric(col) ~ Donor))$'Pr(>F)'[1])
apply(sum.PC$x, 2, function(col)anova(lmer(as.numeric(col) ~ Treatment + (1|Donor)))$'Pr(>F)'[1])

par(mfrow=c(1,2))
plot(sum.PC$x[,1] ~ sum.PC$x[,2], col = Donor, pch = 16, cex = 1.25, xlab = "PC2 (23.2%)", ylab = "PC1 (52.5%)", cex.lab = 1.2, main = "Donor")
legend(-700, 1100, pch=c(16,16,16), col=c("black", "red", "green", "blue"), c("A", "B", "C", "D"), bty="o",  cex=1)

plot(sum.PC$x[,1] ~ sum.PC$x[,2], col = rep(c("#999999", "#E69F00"), 4), pch = 16, cex = 1.25, xlab = "PC2 (23.2%)", ylab = "PC1 (52.5%)", cex.lab = 1.2, main = "Treatment")
legend(-700, 1100, pch=c(16,16,16), col=c("#E69F00", "#999999"), c("12-13-DiHOME", "Untreated"), bty="o",  cex=1)

# 4) How many peaks are lost or gained after 12,13-DiHOME treatment? ====

# how many are present in each treatment?

cpm_raw_vehicle <- cpm_raw[,treatment_info_vehicle_dihome$Treatment == "vehicle"]
cpm_raw_dihome <- cpm_raw[,treatment_info_vehicle_dihome$Treatment == "dihome"]

vehicle_sites <- rownames(cpm_raw_vehicle[rowSums(cpm_raw_vehicle>1) > 3,])
dihome_sites <- rownames(cpm_raw_dihome[rowSums(cpm_raw_dihome>1) > 3,])
sites_venn <- list("Vehicle" = vehicle_sites, "12,13-diHOME, 50uM" = dihome_sites)

ggvenn(
  sites_venn, 
  fill_color = c("#999999", "#E69F00"),
  stroke_size = 0.5, set_name_size = 4
)

gained_sites <- dihome_sites[!dihome_sites %in% vehicle_sites]
lost_sites <- vehicle_sites[!vehicle_sites %in% dihome_sites]

# how many are present in all four dihome treated and absent in all four vehicle?

vehicle_sites_strict <- rownames(cpm_raw_vehicle[rowSums(cpm_raw_vehicle>1) == 4 & rowSums(cpm_raw_dihome>1) == 0,]) # 43
dihome_sites_strict <- rownames(cpm_raw_vehicle[rowSums(cpm_raw_dihome>1) == 4 & rowSums(cpm_raw_vehicle>1) == 0,]) # 6

vehicle_sites_strict_ordered <- names(sort(rowMeans(cpm_raw_vehicle[vehicle_sites_strict,]), decreasing = TRUE))
dihome_sites_strict_ordered <- names(sort(rowMeans(cpm_raw_dihome[dihome_sites_strict,]), decreasing = FALSE))

sites_strict_ordered <- c(vehicle_sites_strict_ordered, dihome_sites_strict_ordered)

all_cpm_raw_ordered <- cbind(cpm_raw_vehicle, cpm_raw_dihome)
all_cpm_raw_ordered_lostorgained <- all_cpm_raw_ordered[sites_strict_ordered,]

treatment_ordered = data.frame("Treatment" = c(rep("Vehicle", 4), rep("12,13-diHOME", 4)))
rownames(treatment_ordered) <- colnames(all_cpm_raw_ordered_lostorgained)
ann_colors = list(Treatment = c("Vehicle"="#999999","12,13-diHOME"="#E69F00"))

pheatmap::pheatmap(all_cpm_raw_ordered_lostorgained, annotation_col = treatment_ordered, cluster_cols = F, cluster_rows = F,annotation_colors=ann_colors)

pdf("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/Figures/43lost_6gained_heatmap_050921.pdf", width = 8, height = 7)
pheatmap::pheatmap(all_cpm_raw_ordered_lostorgained, annotation_col = treatment_ordered, cluster_cols = F, cluster_rows = F,annotation_colors=ann_colors)
dev.off()

# Write out these four files

gained_sites_bed <- cbind(t(as.data.frame(strsplit(gained_sites, "[.]"))), gained_sites)
lost_sites_bed <- cbind(t(as.data.frame(strsplit(lost_sites, "[.]"))), lost_sites)
gained_sites_strict_bed <- cbind(t(as.data.frame(strsplit(dihome_sites_strict, "[.]"))), dihome_sites_strict)
lost_sites_strict_bed <- cbind(t(as.data.frame(strsplit(vehicle_sites_strict, "[.]"))), vehicle_sites_strict)

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs")
write.table(gained_sites_bed, "gained_sites_5447_051021.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(lost_sites_bed, "lost_sites_7569_051021.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(gained_sites_strict_bed, "gained_sites_strict_6_051021.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(lost_sites_strict_bed, "lost_sites_strict_43_051021.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


# 5) Let's perform differential chromatin accessibility analysis (no significant results) ====

# Keep at least present in 3 samples

all_consensus_peaks_data_present <- all_consensus_peaks_data_vehicle_dihome[rowSums(cpm_raw_noted) >= 2,]

## Notes on voom: 
#As noted in the help file for voom, voom can take any of the arguments that are passed onto lmFit, including a blocking variable and correlation.  Using the limma and edgeR packages, an analysis on read counts using voom and a blocking variable could look something like this:
nf <- calcNormFactors(Counts)
design <- model.matrix(~Group)
y <- voom(Counts,design,lib.size=colSums(Counts)*nf)
corfit <- duplicateCorrelation(y,design,block=Subject)
y <- voom(Counts,design,plot=TRUE,lib.size=colSums(Counts)*nf,block=Subject,correlation=corfit$consensus)
fit <- lmFit(y,design,block=Subject,correlation=corfit$consensus)
fit <- eBayes(fit)

# Perform TMM normalization

# To calculate the TMM normalization factors, create a DGElist using the edgeR package: 
dge.nl.nonorm <- DGEList(counts=all_consensus_peaks_data_present)
# Perform TMM normalization using the calcNormFactors function:
dge.nl <- calcNormFactors(dge.nl.nonorm)

# Apply voom transformation
v.nl <- voom(dge.nl ,design=NULL,plot=TRUE)
v.gx<-v.nl$E

# Run duplicateCorrelation on this first voom object: 
design<-model.matrix(~0 + Treatment, groups=treatment_info_vehicle_dihome)
indiv_effect <- duplicateCorrelation(v.gx, design, block = treatment_info_vehicle_dihome$Donor)

# Run voom a second time:
vnl_design2<-voom(dge.nl, design, plot=TRUE, block=treatment_info_vehicle_dihome$Donor, correlation=indiv_effect$consensus.correlation)

# Create a second duplicateCorrelation: 
indiv_effect2 <- duplicateCorrelation(vnl_design2, design, block = treatment_info_vehicle_dihome$Donor)

fit <- lmFit(vnl_design2, design, block = treatment_info_vehicle_dihome$Donor,
             correlation = indiv_effect2$consensus.correlation)
cont_mat <- makeContrasts(dihome = Treatmentdihome - Treatmentvehicle,
                          levels = design)
fit2 <- contrasts.fit(fit, cont_mat)
fit2 <- eBayes(fit2)
results <- decideTests(fit2, p.value = 0.10)
vennDiagram(results, main="FDR 10%")

# Don't really see differentially accessible sites after adjusting for multiple testing

sum(fit2$p.value[,1] < 0.05) # 611

dihome_nom0.05 <- names(fit2$p.value[fit2$p.value[,1] < 0.05,])

dihome_nom0.05_coefficients <- fit2$coefficients[dihome_nom0.05,1] #394/544 have increased accessibility

dihome_nom0.05_coefficients_up <- names(dihome_nom0.05_coefficients)[dihome_nom0.05_coefficients > 0] # 247 had increased accessibility
dihome_nom0.05_coefficients_down <- names(dihome_nom0.05_coefficients)[dihome_nom0.05_coefficients < 0] # 364 had decreased accessibility

dihome_nom0.05_coefficients_up_bed <- cbind(t(as.data.frame(strsplit(dihome_nom0.05_coefficients_up, "[.]"))), dihome_nom0.05_coefficients_up)
dihome_nom0.05_coefficients_down_bed <- cbind(t(as.data.frame(strsplit(dihome_nom0.05_coefficients_down, "[.]"))), dihome_nom0.05_coefficients_down)

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs")
write.table(dihome_nom0.05_coefficients_up_bed, "dihome_nom0.05_coefficients_up247_051021.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(dihome_nom0.05_coefficients_down_bed, "dihome_nom0.05_coefficients_down364_051021.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# 6) Test global differences between the treatment groups (negative) ====

fc_data <- topTable(fit2, number = nrow(fit2))
dihome_fc <- fc_data[,1]

# I don't think there are dramatic differences in chromatin accessibility
# 65190 with greater accessibility and 63772 with lower accessibility

# 7) Let's obtain annotations for the different regions ====

# Let's get Venn Diagrams of the annotations

par(mfrow=c(2,2))

# lost_sites

slices_overlap_lost <- c(86, 0, 27, 108, 1, 122, 4072, 3032, 106, 15, 0, 0)
lbls_overlap_lost <- c("", "", "", "TTS", "", "exon", "intron", "intergenic", "promoter", "", "", "")
pie(slices_overlap_lost, labels = lbls_overlap_lost, main="Lost sites\n(N=7,569)")

# gained sites

slices_overlap_gained <- c(47, 0, 25, 91, 3, 108, 2701, 2289, 160, 22, 0, 0)
lbls_overlap_gained <- c("", "", "", "TTS", "", "exon", "intron", "intergenic", "promoter", "", "", "")
pie(slices_overlap_gained, labels = lbls_overlap_gained, main="Gained sites\n(N=5,447)")

# lost_sites strict

slices_overlap_lost_strict <- c(0, 0, 1, 0, 0, 0, 23, 19, 0, 0, 0, 0)
lbls_overlap_lost_strict <- c("", "", "ncRNA", "", "", "", "intron", "intergenic", "", "", "", "")
pie(slices_overlap_lost_strict, labels = lbls_overlap_lost_strict, main="Lost sites (strict)\n(N=43)")

# gained sites strict

slices_overlap_gained_strict <- c(0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0)
lbls_overlap_gained_strict <- c("", "", "", "", "", "", "intron", "intergenic", "", "", "", "")
pie(slices_overlap_gained_strict, labels = lbls_overlap_gained_strict, main="Gained sites (strict)\n(N=6)")

# downreg_sites

slices_overlap_downreg <- c(4, 0, 1, 5, 0, 2, 193, 158, 1, 0, 0, 0)
lbls_overlap_downreg <- c("", "", "", "", "", "", "intron", "intergenic", "", "", "", "")
pie(slices_overlap_downreg, labels = lbls_overlap_downreg, main="Downregulated sites")

# upreg_sites

slices_overlap_upreg <- c(0, 0, 1, 3, 0, 12, 44, 35, 149, 3, 0, 0)
lbls_overlap_upreg <- c("", "", "", "", "", "exon", "intron", "intergenic", "promoter", "", "", "")
pie(slices_overlap_upreg, labels = lbls_overlap_upreg, main="Upregulated sites")

## plot the transcription factor binding sites

perc <- c(3.32, 10.36, 9.46, 4.65, 2.55, 9.00, 8.17, 3.82)
tfs <- c("ISRE", "PRDM1", "IRF3", "IRF2", "ISRE", "PRDM1", "IRF3", "IRF2")
treat_control <- c(rep("Target", 4), rep("Background", 4))
tfbs_analysis <- as.data.frame(cbind(perc, tfs, treat_control))
tfbs_analysis$tfs <- factor(tfbs_analysis$tfs, levels = c("IRF2", "IRF3", "PRDM1", "ISRE"))

tfbs_plot <- ggplot(data=tfbs_analysis, aes(x=as.factor(tfs), y=as.numeric(perc), fill=as.factor(treat_control))) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() +
  xlab("") + ylab("% of regions containing motif") +
  scale_fill_manual(values=c("#000000", "#660000")) +
  theme(axis.title.y=element_text(size = 10), 
        axis.text.y = element_text(size = 10, hjust = 1, face="bold"), 
        axis.title=element_text(size = 10), 
        plot.title = element_text(size = 12, hjust = 0.5), 
        legend.position = "none",
        panel.background = element_rect(fill = NA, color = "black"))
  







## human RNA-seq analysis ====
# Goal: Perform differential expression of genes in dihome study

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_human")


#BiocManager::install("biomaRt")
library(biomaRt)
library(plyr)
library(limma)
library(edgeR)
library(WGCNA)
library(flashClust)
library(lmerTest)
library(pheatmap)

# 1) Import data

genecount <- read.table("rnaseq_Mf_counts.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
pheno_data <- read.table("pheno_data.txt", stringsAsFactors = FALSE,  header = TRUE)

#### 1) Remove X, Y, and M genes ====

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "uswest")
genes <- rownames(genecount)
genes_clean <- gsub("\\..*","",genes)
G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name", "gene_biotype") ,values = genes_clean,mart = mart)
G_protein <- G_list[G_list$gene_biotype == "protein_coding",]
G_protein_auto <- G_protein[G_protein$chromosome_name == "1" | G_protein$chromosome_name == "2" | G_protein$chromosome_name == "3" | G_protein$chromosome_name == "4" | G_protein$chromosome_name == "5" | G_protein$chromosome_name == "6" | G_protein$chromosome_name == "7" | G_protein$chromosome_name == "8" | G_protein$chromosome_name == "9" | G_protein$chromosome_name == "10" | G_protein$chromosome_name == "11" | G_protein$chromosome_name == "12" | G_protein$chromosome_name == "13" | G_protein$chromosome_name == "14" | G_protein$chromosome_name == "15" | G_protein$chromosome_name == "16" | G_protein$chromosome_name == "17" | G_protein$chromosome_name == "18" | G_protein$chromosome_name == "19" | G_protein$chromosome_name == "20" | G_protein$chromosome_name == "21" | G_protein$chromosome_name == "22",]
G_protein_auto_complete <- G_protein_auto[G_protein_auto$hgnc_symbol != "",] # recent version of ensembl includes genes that have no hgnc name; also need to remove 4 genes that are duplicated
G_protein_auto_complete_final <- G_protein_auto_complete[-which(G_protein_auto_complete$hgnc_symbol %in% c("TBCE", "PINX1", "HERC3", "POLR2J3", "SIGLEC5")),]

write.csv(G_protein_auto_complete_final, "biomart_ids_human_dihome_071621.csv", quote = FALSE, row.names = FALSE)

# 18439 protein-coding and autosomal genes

# Rename ENSG to the human gene names

prot_aut_code <- intersect(G_protein_auto_complete_final$ensembl_gene_id, genes_clean)
prot_aut_code_index <- match(prot_aut_code, genes_clean)
ENSG_replace <- mapvalues(genes_clean, from=G_protein_auto_complete_final$ensembl_gene_id, to=G_protein_auto_complete_final$hgnc_symbol)

genecount_data_prot_aut <- genecount[prot_aut_code_index,]
ENSG_replace_prot_aut <- ENSG_replace[prot_aut_code_index]

rownames(genecount_data_prot_aut) <- ENSG_replace_prot_aut

# 18438 protein-coding and autosomal genes

#### 2) Pre-processing from the RNAseq analysis is as easy as 1-2-3 with limma, Glimma, and edgeR (subsetting lowly expressed genes) ====

## a) Check the density of log-RPM values for raw pre-filtered data

genecount_data_prot_aut[is.na(genecount_data_prot_aut)] <- 0

# How many genes have a count of 0 in all samples?  1812/18438 genes
table(rowSums(genecount_data_prot_aut==0)==24)

## b) Remove lowly expressed genes

cpm_raw <- cpm(genecount_data_prot_aut, log=FALSE)
lcpm_raw <- cpm(genecount_data_prot_aut, log=TRUE)

# Also remove lowly expressed genes (those absent in > 25% of subjects)
exp_genes <- rownames(cpm_raw[rowSums(cpm_raw>1)>=3,]) # 11960/18438 genes present in 20% of our sample; needed a high cutoff of 10 because of the small library sizes
cpm_raw_exp <- cpm_raw[match(exp_genes, rownames(cpm_raw)),]
genecount_data_prot_aut_exp <- genecount_data_prot_aut[match(exp_genes, rownames(cpm_raw)),]
lcpm_exp <- cpm(genecount_data_prot_aut_exp, log=TRUE)

setwd("/Users/caraporsche/Desktop/diHOME RNAseq/")
pdf("dihome_density_plots_061721.pdf", height = 6, width = 8)
par(mfrow=c(1,2))
plot(density(lcpm_raw[,1]), col = "Blue", lwd=2, las=2, main = "", xlab="", ylim=c(0,0.45))
title(main="A. Raw data", xlab="Log-CPM")
abline(v=0, lty=3)
for (i in 2:ncol(lcpm_raw)){
  den <- density(lcpm_raw[,i])
  den <- density(lcpm_raw[,i])
}
plot(density(lcpm_exp[,1]), col = "Purple", lwd=2, las=2, main = "", xlab="", ylim=c(0,0.45))
title(main="B. Filtered data", xlab="Log-CPM")
abline(v=0, lty=3)
for (i in 2:ncol(lcpm_exp)){
  den <- density(lcpm_exp[,i])
  lines(den$x, den$y, lwd=2, col = "Purple")
}
dev.off() # not that big of a difference but 18442-10477 = 7965 genes removed for low expression though

## c) Normalize the data

# To calculate the TMM normalization factors, create a DGElist using the edgeR package: 
dge.nl.nonorm <- DGEList(counts=genecount_data_prot_aut_exp)
# Perform TMM normalization using the calcNormFactors function:
dge.nl <- calcNormFactors(dge.nl.nonorm)

# Boxplots of the unnormalized and normalized values: 
setwd("/Users/caraporsche/Desktop/diHOME RNAseq/")
pdf("dihome_norm_unnorm_boxplots_061721.pdf")
lcpm_nonorm <- cpm(dge.nl.nonorm, log=TRUE)
boxplot(lcpm_nonorm[,1:24], las=2, main="")
title(main="Unnormalized data, 1-24",ylab="Log-CPM")

lcpm_norm <- cpm(dge.nl, log=TRUE)
boxplot(lcpm_norm[,1:24], las=2, main="")
title(main="Normalized data (TMM), 1-24",ylab="Log-CPM")
dev.off()

# TMM normalization and voom

v.nl <- voom(dge.nl,design=NULL,plot=TRUE)
save(v.nl, file="/Users/caraporsche/Desktop/diHOME RNAseq/dihome_human_mac_RNAseq_070121.Rdata")

#### 3) PCA =====

# Let's print out the cpm-adjusted file that includes the samples treated with LPS

load("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_human/dihome_human_mac_RNAseq_070121.Rdata")
v.nl_human_geo <- v.nl$E
colnames(v.nl_human_geo) <- c("DiHOME24hA", "DiHOME24hB", "DiHOME4hA", "DiHOME4hB", "DiHOME8hA", "DiHOME8hB", "DiHOMEp24hA", "DiHOMEp24hB", "DiHOMEp4hA", "DiHOMEp4hB", "DiHOMEp8hA", "DiHOMEp8hB", "DMSO24hA", "DMSO24hB", "DMSO4hA", "DMSO4hB", "DMSO8hA", "DMSO8hB", "Media24hA", "Media24hB", "Media4hA", "Media4hB", "Media8hA", "Media8hB")
write.table(v.nl_human_geo, "/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/GEO/GEO_Human_RNAseq_24.txt", quote = FALSE, sep = "\t", )

# A) Raw

treatment <- as.factor(pheno_data$Treatment)
allergen <- as.factor(pheno_data$Allergen)
timepoint <- as.factor(pheno_data$Timepoint)
patient <- as.factor(pheno_data$Patient)

sum.PC <- prcomp(t(v.nl$E), scale=FALSE, center=TRUE)
sumsum <- summary(sum.PC)
prop.var<-sumsum$importance
prop.var[1:3,1:10] 

pc_list<-c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
m.covars<-cbind.data.frame(treatment, allergen, timepoint, patient)

pval.pca1=matrix(ncol=ncol(m.covars), nrow=10)
colnames(pval.pca1)=colnames(m.covars)
rownames(pval.pca1)=paste(pc_list, " ", "(", round(prop.var[2,1:length(pc_list)], 3), ")", sep = "")

for(i in 1:ncol(m.covars))
{
  for(j in 1:length(pc_list))
  {
    data1= lm(sum.PC$x[,j]~m.covars[,i])
    pval.pca1[j,i]=anova(data1)$'Pr(>F)'[1]
  }
}

pval.pca1[pval.pca1 > 0.05] <- "ns"

setwd("/Users/caraporsche/Desktop/diHOME RNAseq/")
write.csv(pval.pca1, file = "PCA_Raw_061721.csv", quote=FALSE)

#### 4) Identify differentially expressed genes between treatments ====

# treatments: diHOME, DMSO, Media
# allergens: Pnut, NoPnut
# timepoints: 4h, 8h, 24h
# patients: A, B

# comparisons: diHOME vs DMSO, dihome+peanut vs DMSO, and dihome + peanut vs dihome

# A) Linear mixed effects model, controlling for patient and timepoint

load("/Users/caraporsche/Desktop/diHOME RNAseq/dihome_human_mac_RNAseq_070121.Rdata")
voom_norm <- v.nl$E
voom_norm_nomedia <- voom_norm[,1:18]

pheno_data_nomedia <- pheno_data[1:18,]
pheno_data_nomedia$treatment_allergen <- paste(pheno_data_nomedia$Treatment, pheno_data_nomedia$Allergen, sep = "_")

## Need three separate sets of gene expression and phenotype files

#### phenotype data

pheno_data_nomedia_dihomepnut_dmso <- pheno_data_nomedia[!pheno_data_nomedia$treatment_allergen %in% "diHOME_NoPnut",]
pheno_data_nomedia_dihomepnut_dihome <- pheno_data_nomedia[!pheno_data_nomedia$treatment_allergen %in% "DMSO_NoPnut",]
pheno_data_nomedia_dihome_dmso <- pheno_data_nomedia[!pheno_data_nomedia$treatment_allergen %in% "diHOME_Pnut",]

#### ordering treatment variables; dihome+peanut or dihome or peanut should be second for interpreting the fold changes

pheno_data_nomedia_dihomepnut_dmso$treatment_allergen <- factor(pheno_data_nomedia_dihomepnut_dmso$treatment_allergen, levels = c("DMSO_NoPnut", "diHOME_Pnut"))
pheno_data_nomedia_dihomepnut_dihome$treatment_allergen <- factor(pheno_data_nomedia_dihomepnut_dihome$treatment_allergen, levels = c("diHOME_NoPnut", "diHOME_Pnut"))
pheno_data_nomedia_dihome_dmso$treatment_allergen <- factor(pheno_data_nomedia_dihome_dmso$treatment_allergen, levels = c("DMSO_NoPnut", "diHOME_NoPnut"))

#### gene expression data

voom_norm_nomedia_dihomepnut_dmso <- voom_norm_nomedia[,!pheno_data_nomedia$treatment_allergen %in% "diHOME_NoPnut"]
voom_norm_nomedia_dihomepnut_dihome <- voom_norm_nomedia[,!pheno_data_nomedia$treatment_allergen %in% "DMSO_NoPnut"]
voom_norm_nomedia_dihome_dmso <- voom_norm_nomedia[,!pheno_data_nomedia$treatment_allergen %in% "diHOME_Pnut"]

#### log2 fold change

lme_coef_dihomepnut_dmso <- apply(voom_norm_nomedia_dihomepnut_dmso, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + as.factor(Timepoint) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dmso))$coefficients[2,1])
lme_coef_dihomepnut_dihome <- apply(voom_norm_nomedia_dihomepnut_dihome, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + as.factor(Timepoint) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dihome))$coefficients[2,1])
lme_coef_dihome_dmso <- apply(voom_norm_nomedia_dihome_dmso, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + as.factor(Timepoint) + (1|Patient), data = pheno_data_nomedia_dihome_dmso))$coefficients[2,1])

#### pvalue 

lme_pvals_dihomepnut_dmso <- apply(voom_norm_nomedia_dihomepnut_dmso, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + as.factor(Timepoint) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dmso))$coefficients[2,5])
lme_pvals_dihomepnut_dihome <- apply(voom_norm_nomedia_dihomepnut_dihome, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + as.factor(Timepoint) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dihome))$coefficients[2,5])
lme_pvals_dihome_dmso <- apply(voom_norm_nomedia_dihome_dmso, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + as.factor(Timepoint) + (1|Patient), data = pheno_data_nomedia_dihome_dmso))$coefficients[2,5])

#### adjusted p value

lme_pvals_dihomepnut_dmso_adjusted <- p.adjust(lme_pvals_dihomepnut_dmso, method = "fdr") # 1308
lme_pvals_dihomepnut_dihome_adjusted <- p.adjust(lme_pvals_dihomepnut_dihome, method = "fdr") # 164
lme_pvals_dihome_dmso_adjusted <- p.adjust(lme_pvals_dihome_dmso, method = "fdr") # 20

dihomepnut_dmso_pvals_data <- cbind(lme_coef_dihomepnut_dmso, lme_pvals_dihomepnut_dmso, lme_pvals_dihomepnut_dmso_adjusted)
colnames(dihomepnut_dmso_pvals_data) <- c("log2FC", "Unadjusted", "Adjusted")
write.csv(dihomepnut_dmso_pvals_data, "dihomepnut_dmso_pvals_data_071421.csv", quote = FALSE)

dihomepnut_dihome_pvals_data <- cbind(lme_coef_dihomepnut_dihome, lme_pvals_dihomepnut_dihome, lme_pvals_dihomepnut_dihome_adjusted)
colnames(dihomepnut_dihome_pvals_data) <- c("log2FC", "Unadjusted", "Adjusted")
write.csv(dihomepnut_dihome_pvals_data, "dihomepnut_dihome_pvals_data_071421.csv", quote = FALSE)

dihome_dmso_pvals_data <- cbind(lme_coef_dihome_dmso, lme_pvals_dihome_dmso, lme_pvals_dihome_dmso_adjusted)
colnames(dihome_dmso_pvals_data) <- c("log2FC", "Unadjusted", "Adjusted")
write.csv(dihome_dmso_pvals_data, "dihome_dmso_pvals_data_071421.csv", quote = FALSE)

# B) Do the same three comparisons for each timepoint (4, 8, 24)

## 4 hours

pheno_data_nomedia_dihomepnut_dmso_4 <- pheno_data_nomedia_dihomepnut_dmso[pheno_data_nomedia_dihomepnut_dmso$Timepoint %in% "4h",]
pheno_data_nomedia_dihomepnut_dihome_4 <- pheno_data_nomedia_dihomepnut_dihome[pheno_data_nomedia_dihomepnut_dihome$Timepoint %in% "4h",]
pheno_data_nomedia_dihome_dmso_4 <- pheno_data_nomedia_dihome_dmso[pheno_data_nomedia_dihome_dmso$Timepoint %in% "4h",]

voom_norm_nomedia_dihomepnut_dmso_4 <- voom_norm_nomedia_dihomepnut_dmso[,pheno_data_nomedia_dihomepnut_dmso$Timepoint %in% "4h"]
voom_norm_nomedia_dihomepnut_dihome_4 <- voom_norm_nomedia_dihomepnut_dihome[,pheno_data_nomedia_dihomepnut_dihome$Timepoint %in% "4h"]
voom_norm_nomedia_dihome_dmso_4 <- voom_norm_nomedia_dihome_dmso[,pheno_data_nomedia_dihome_dmso$Timepoint %in% "4h"]

lme_coef_dihomepnut_dmso_4 <- apply(voom_norm_nomedia_dihomepnut_dmso_4, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dmso_4))$coefficients[2,1])
lme_coef_dihomepnut_dihome_4 <- apply(voom_norm_nomedia_dihomepnut_dihome_4, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dihome_4))$coefficients[2,1])
lme_coef_dihome_dmso_4 <- apply(voom_norm_nomedia_dihome_dmso_4, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihome_dmso_4))$coefficients[2,1])

lme_pvals_dihomepnut_dmso_4 <- apply(voom_norm_nomedia_dihomepnut_dmso_4, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dmso_4))$coefficients[2,5])
lme_pvals_dihomepnut_dihome_4 <- apply(voom_norm_nomedia_dihomepnut_dihome_4, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dihome_4))$coefficients[2,5])
lme_pvals_dihome_dmso_4 <- apply(voom_norm_nomedia_dihome_dmso_4, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihome_dmso_4))$coefficients[2,5])

adjust_save_pvals <- function(coef, unadjusted_pvals, output_file){
  
  adjusted_p <- p.adjust(unadjusted_pvals, method = "fdr")
  pval_data <- cbind(coef, unadjusted_pvals, adjusted_p)
  colnames(pval_data) <- c("log2FC", "Unadjusted", "Adjusted")
  write.csv(pval_data, output_file, quote = FALSE)
  return(pval_data)
  
}

setwd("/Users/caraporsche/Desktop/diHOME RNAseq")
lme_pvals_dihomepnut_dmso_4_adjusted <- adjust_save_pvals(lme_coef_dihomepnut_dmso_4, lme_pvals_dihomepnut_dmso_4, "dihomepnut_dmso_4_pvals_data_070121.csv")
lme_pvals_dihomepnut_dihome_4_adjusted <- adjust_save_pvals(lme_coef_dihomepnut_dihome_4, lme_pvals_dihomepnut_dihome_4, "dihomepnut_dihome_4_pvals_data_070121.csv")
lme_pvals_dihome_dmso_4_adjusted <- adjust_save_pvals(lme_coef_dihome_dmso_4, lme_pvals_dihome_dmso_4, "dihome_dmso_4_pvals_data_070121.csv")

## 8 hours

pheno_data_nomedia_dihomepnut_dmso_8 <- pheno_data_nomedia_dihomepnut_dmso[pheno_data_nomedia_dihomepnut_dmso$Timepoint %in% "8h",]
pheno_data_nomedia_dihomepnut_dihome_8 <- pheno_data_nomedia_dihomepnut_dihome[pheno_data_nomedia_dihomepnut_dihome$Timepoint %in% "8h",]
pheno_data_nomedia_dihome_dmso_8 <- pheno_data_nomedia_dihome_dmso[pheno_data_nomedia_dihome_dmso$Timepoint %in% "8h",]

voom_norm_nomedia_dihomepnut_dmso_8 <- voom_norm_nomedia_dihomepnut_dmso[,pheno_data_nomedia_dihomepnut_dmso$Timepoint %in% "8h"]
voom_norm_nomedia_dihomepnut_dihome_8 <- voom_norm_nomedia_dihomepnut_dihome[,pheno_data_nomedia_dihomepnut_dihome$Timepoint %in% "8h"]
voom_norm_nomedia_dihome_dmso_8 <- voom_norm_nomedia_dihome_dmso[,pheno_data_nomedia_dihome_dmso$Timepoint %in% "8h"]

lme_coef_dihomepnut_dmso_8 <- apply(voom_norm_nomedia_dihomepnut_dmso_8, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dmso_8))$coefficients[2,1])
lme_coef_dihomepnut_dihome_8 <- apply(voom_norm_nomedia_dihomepnut_dihome_8, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dihome_8))$coefficients[2,1])
lme_coef_dihome_dmso_8 <- apply(voom_norm_nomedia_dihome_dmso_8, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihome_dmso_8))$coefficients[2,1])

lme_pvals_dihomepnut_dmso_8 <- apply(voom_norm_nomedia_dihomepnut_dmso_8, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dmso_8))$coefficients[2,5])
lme_pvals_dihomepnut_dihome_8 <- apply(voom_norm_nomedia_dihomepnut_dihome_8, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dihome_8))$coefficients[2,5])
lme_pvals_dihome_dmso_8 <- apply(voom_norm_nomedia_dihome_dmso_8, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihome_dmso_8))$coefficients[2,5])

lme_pvals_dihomepnut_dmso_8_adjusted <- adjust_save_pvals(lme_coef_dihomepnut_dmso_8, lme_pvals_dihomepnut_dmso_8, "dihomepnut_dmso_8_pvals_data_070121.csv")
lme_pvals_dihomepnut_dihome_8_adjusted <- adjust_save_pvals(lme_coef_dihomepnut_dihome_8, lme_pvals_dihomepnut_dihome_8, "dihomepnut_dihome_8_pvals_data_070121.csv")
lme_pvals_dihome_dmso_8_adjusted <- adjust_save_pvals(lme_coef_dihome_dmso_8, lme_pvals_dihome_dmso_8, "dihome_dmso_8_pvals_data_070121.csv")

## 24 hours

pheno_data_nomedia_dihomepnut_dmso_24 <- pheno_data_nomedia_dihomepnut_dmso[pheno_data_nomedia_dihomepnut_dmso$Timepoint %in% "24h",]
pheno_data_nomedia_dihomepnut_dihome_24 <- pheno_data_nomedia_dihomepnut_dihome[pheno_data_nomedia_dihomepnut_dihome$Timepoint %in% "24h",]
pheno_data_nomedia_dihome_dmso_24 <- pheno_data_nomedia_dihome_dmso[pheno_data_nomedia_dihome_dmso$Timepoint %in% "24h",]

voom_norm_nomedia_dihomepnut_dmso_24 <- voom_norm_nomedia_dihomepnut_dmso[,pheno_data_nomedia_dihomepnut_dmso$Timepoint %in% "24h"]
voom_norm_nomedia_dihomepnut_dihome_24 <- voom_norm_nomedia_dihomepnut_dihome[,pheno_data_nomedia_dihomepnut_dihome$Timepoint %in% "24h"]
voom_norm_nomedia_dihome_dmso_24 <- voom_norm_nomedia_dihome_dmso[,pheno_data_nomedia_dihome_dmso$Timepoint %in% "24h"]

lme_coef_dihomepnut_dmso_24 <- apply(voom_norm_nomedia_dihomepnut_dmso_24, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dmso_24))$coefficients[2,1])
lme_coef_dihomepnut_dihome_24 <- apply(voom_norm_nomedia_dihomepnut_dihome_24, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dihome_24))$coefficients[2,1])
lme_coef_dihome_dmso_24 <- apply(voom_norm_nomedia_dihome_dmso_24, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihome_dmso_24))$coefficients[2,1])

lme_pvals_dihomepnut_dmso_24 <- apply(voom_norm_nomedia_dihomepnut_dmso_24, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dmso_24))$coefficients[2,5])
lme_pvals_dihomepnut_dihome_24 <- apply(voom_norm_nomedia_dihomepnut_dihome_24, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihomepnut_dihome_24))$coefficients[2,5])
lme_pvals_dihome_dmso_24 <- apply(voom_norm_nomedia_dihome_dmso_24, 1, function(x) summary(lmer(as.numeric(x) ~ as.factor(treatment_allergen) + (1|Patient), data = pheno_data_nomedia_dihome_dmso_24))$coefficients[2,5])

lme_pvals_dihomepnut_dmso_24_adjusted <- adjust_save_pvals(lme_coef_dihomepnut_dmso_24, lme_pvals_dihomepnut_dmso_24, "dihomepnut_dmso_24_pvals_data_070121.csv")
lme_pvals_dihomepnut_dihome_24_adjusted <- adjust_save_pvals(lme_coef_dihomepnut_dihome_24, lme_pvals_dihomepnut_dihome_24, "dihomepnut_dihome_24_pvals_data_070121.csv")
lme_pvals_dihome_dmso_24_adjusted <- adjust_save_pvals(lme_coef_dihome_dmso_24, lme_pvals_dihome_dmso_24, "dihome_dmso_24_pvals_data_070121.csv")

# C) Identify DE genes for the overall tests in (A)

dihomepnut_dmso_results <- read.csv("dihomepnut_dmso_pvals_data_070121.csv", header = TRUE, row.names = 1) # 1308
dihomepnut_dihome_results <- read.csv("dihomepnut_dihome_pvals_data_070121.csv", header = TRUE, row.names = 1) # 164; will not use
dihome_dmso_results <- read.csv("dihome_dmso_pvals_data_070121.csv", header = TRUE, row.names = 1) # 20

dihomepnut_dmso_de_genes <- rownames(dihomepnut_dmso_results[dihomepnut_dmso_results[,3] < 0.05,])
dihome_dmso_de_genes <- rownames(dihome_dmso_results[dihome_dmso_results[,2] < 0.05,])
shared_dihomepnut_dihome_genes <- intersect(dihomepnut_dmso_de_genes, dihome_dmso_de_genes)

# 1296 unique DE genes for dihomepnut vs dmso, 8 unique DE genes for dihome vs dmso, and 12 shared

# IL10, TGFB, IL6, IL1B, TNFA, NFKB

# D) Identify DE genes for timepoint-specific tests 

### dihomepnut vs dmso results only

dihomepnut_dmso_4_results <- read.csv("dihomepnut_dmso_4_pvals_data_070121.csv", header = TRUE, row.names = 1) # 1
dihomepnut_dmso_8_results <- read.csv("dihomepnut_dmso_8_pvals_data_070121.csv", header = TRUE, row.names = 1) # 0
dihomepnut_dmso_24_results <- read.csv("dihomepnut_dmso_24_pvals_data_070121.csv", header = TRUE, row.names = 1) # 0

### dihome vs dmso results only

dihome_dmso_4_results <- read.csv("dihome_dmso_4_pvals_data_070121.csv", header = TRUE, row.names = 1) # 3
dihome_dmso_8_results <- read.csv("dihome_dmso_8_pvals_data_070121.csv", header = TRUE, row.names = 1) # 0
dihome_dmso_24_results <- read.csv("dihome_dmso_24_pvals_data_070121.csv", header = TRUE, row.names = 1) # 0

# no differentially expressed genes at most timepoints; lack of statistical power to detect de genes

# E) Heat map (Figure 5) excluding media

# Extract the significantly differentially expressed genes after adjusting for timepoints/patient

dihomepnut_dmso_data <- read.csv("/Users/caraporsche/Desktop/diHOME RNAseq/dihomepnut_dmso_pvals_data_070121.csv")
dihomepnut_dmso_data_sig <- dihomepnut_dmso_data[dihomepnut_dmso_data[,4] < 0.05,]
dihomepnut_dmso_data_sig_order <- dihomepnut_dmso_data_sig[order(dihomepnut_dmso_data_sig[,3]),]
dihomepnut_dmso_data_sig_order_top100 <- dihomepnut_dmso_data_sig_order[1:100,]

voom_norm_nomedia_top100 <- voom_norm_nomedia[rownames(voom_norm_nomedia) %in% dihomepnut_dmso_data_sig_order_top100[,1],]
voom_norm_nomedia_top100_match <- voom_norm_nomedia_top100[match(dihomepnut_dmso_data_sig_order_top100[,1], rownames(voom_norm_nomedia_top100)),]

annotated_treatment <- as.data.frame(cbind(c("24h", "24h", "4h", "4h", "8h", "8h", "24h", "24h", "4h", "4h", "8h", "8h", "24h", "24h", "4h", "4h", "8h", "8h"),
                                           c(rep("12,13-diHOME", 6), rep("12,13-diHOME + Peanut", 6), rep("DMSO", 6))))
rownames(annotated_treatment) <- colnames(voom_norm_nomedia_top100_match)
colnames(annotated_treatment) <- c("Timepoint", "Treatment")
annotated_treatment$Timepoint <- factor(annotated_treatment$Timepoint, levels = c("4h", "8h", "24h"))
annotated_treatment$Treatment <- factor(annotated_treatment$Treatment, levels = c("DMSO", "12,13-diHOME", "12,13-diHOME + Peanut"))

# plot heatmap

heatmap_top100 <- pheatmap::pheatmap(voom_norm_nomedia_top100_match, annotation_col = annotated_treatment, cluster_cols = T, cluster_rows = T)
hc_top100 <- hclust(dist(t(voom_norm_nomedia_top100_match)))
hc_top100_neworder <- c(hc_top100$order[6:length(hc_top100$order)], hc_top100$order[1:5])
hc_top100$order <- hc_top100_neworder
hc_top100$labels <- rep(c("A", "B"), 9)

heatmap_top100_reordered <- pheatmap::pheatmap(voom_norm_nomedia_top100_match, 
                                               annotation_col = annotated_treatment, 
                                               cluster_cols = hc_top100, 
                                               cluster_rows = T,
                                               labels_col = hc_top100$labels)

pdf("human_rnaseq_top100_de_heatmap_082021.pdf", height = 5, width = 7)
heatmap_top100_reordered
dev.off() 

#### 5) Pathway analysis ====

## A) Installation process

library(rWikiPathways)

load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}

cytoscapePing()

installApp('WikiPathways') 
installApp('CyTargetLinker') 
installApp('stringApp') 
installApp('enrichmentMap')

# When using bitr, toType parameter should be "SYMBOL".

## B) Load data and obtain original ENSG ids and convert to ensembl ids using bitr

G_protein_auto_complete_final2 <- G_protein_auto_complete_final[,1:2]

# dihome only data

dihome_dmso_de_data <- read.csv("dihome_dmso_pvals_data_071421.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

dihome_dmso_de_data_final <- cbind(rownames(dihome_dmso_de_data), dihome_dmso_de_data)
colnames(dihome_dmso_de_data_final) <- c("hgnc_symbol", colnames(dihome_dmso_de_data))
dihome_dmso_de_data_combined <- merge(G_protein_auto_complete_final2, dihome_dmso_de_data_final, by = "hgnc_symbol", all = FALSE, sort = FALSE)

dihome_dmso_entrez_ids <- bitr(dihome_dmso_de_data_combined[,2], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
dihome_dmso_entrez_ids_nondups <- dihome_dmso_entrez_ids[-which(duplicated(dihome_dmso_entrez_ids[,1])),]
dihome_dmso_entrez_ids_nondups2 <- dihome_dmso_entrez_ids_nondups[-which(duplicated(dihome_dmso_entrez_ids_nondups[,2])),]
dihome_dmso_de_data_combined_entrez <- merge(dihome_dmso_de_data_combined, dihome_dmso_entrez_ids_nondups2, by.x = "ensembl_gene_id", by.y = "ENSEMBL", all = FALSE, sort = FALSE)

# dihome + peanut data

dihomepnut_dmso_de_data <- read.csv("dihomepnut_dmso_pvals_data_070121.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

dihomepnut_dmso_de_data_final <- cbind(rownames(dihomepnut_dmso_de_data), dihomepnut_dmso_de_data)
colnames(dihomepnut_dmso_de_data_final) <- c("hgnc_symbol", colnames(dihomepnut_dmso_de_data))
dihomepnut_dmso_de_data_combined <- merge(G_protein_auto_complete_final2, dihomepnut_dmso_de_data_final, by = "hgnc_symbol", all = FALSE, sort = FALSE)

dihomepnut_dmso_entrez_ids <- bitr(dihomepnut_dmso_de_data_combined[,2], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
dihomepnut_dmso_entrez_ids_nondups <- dihomepnut_dmso_entrez_ids[-which(duplicated(dihomepnut_dmso_entrez_ids[,1])),]
dihomepnut_dmso_entrez_ids_nondups2 <- dihomepnut_dmso_entrez_ids_nondups[-which(duplicated(dihomepnut_dmso_entrez_ids_nondups[,2])),]
dihomepnut_dmso_de_data_combined_entrez <- merge(dihomepnut_dmso_de_data_combined, dihomepnut_dmso_entrez_ids_nondups2, by.x = "ensembl_gene_id", by.y = "ENSEMBL", all = FALSE, sort = FALSE)

## C) Pathway analysis using GSEA (need to reconvert gene names to ensembl)

wp.hs.gmt <- "wikipathways-202110610-gmt-Homo_sapiens.gmt" # downloaded the human .gmt file manually and placed in directory
wp2gene <- readPathwayGMT(wp.hs.gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
wpid2gene
wpid2name

# dihome only data

dihome_dmso_de_data_combined_entrez$fcsign <- sign(as.numeric(dihome_dmso_de_data_combined_entrez$log2FC))
dihome_dmso_de_data_combined_entrez$logfdr <- -log10(as.numeric(dihome_dmso_de_data_combined_entrez$Unadjusted))
dihome_dmso_de_data_combined_entrez$sig <- dihome_dmso_de_data_combined_entrez$logfdr/dihome_dmso_de_data_combined_entrez$fcsign
dihome_dmso_gsea_sig_expr <- dihome_dmso_de_data_combined_entrez$sig
names(dihome_dmso_gsea_sig_expr) <- dihome_dmso_de_data_combined_entrez$ENTREZID
dihome_dmso_gsea_sig_expr_sorted <- sort(dihome_dmso_gsea_sig_expr, decreasing = TRUE)

dihome_dmso_gwp_sig_expr <- clusterProfiler::GSEA(
  dihome_dmso_gsea_sig_expr_sorted,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

dihome_dmso_gwp_sig_expr2 <- setReadable(dihome_dmso_gwp_sig_expr, org.Hs.eg.db, keyType = "ENTREZID")

dihome_dmso_gwp_sig_expr_df = data.frame(ID=dihome_dmso_gwp_sig_expr2$ID,
                                         Description=dihome_dmso_gwp_sig_expr2$Description,
                                         enrichmentScore=dihome_dmso_gwp_sig_expr2$enrichmentScore,
                                         NES=dihome_dmso_gwp_sig_expr2$NES,
                                         pvalue=dihome_dmso_gwp_sig_expr2$pvalue,
                                         p.adjust=dihome_dmso_gwp_sig_expr2$p.adjust,
                                         rank=dihome_dmso_gwp_sig_expr2$rank,
                                         leading_edge=dihome_dmso_gwp_sig_expr2$leading_edge,
                                         core_enrichment=dihome_dmso_gwp_sig_expr2$core_enrichment
)

dihome_dmso_gwp_sig_expr_df[which(dihome_dmso_gwp_sig_expr_df$NES > 1),]
dihome_dmso_gwp_sig_expr_df[which(dihome_dmso_gwp_sig_expr_df$NES < 1),]

setwd("/Users/caraporsche/Desktop/diHOME RNAseq")
write.csv(dihome_dmso_gwp_sig_expr_df, "dihome_only_gsea_20pathways_082021.csv", row.names = FALSE, quote = FALSE)
write.csv(dihome_dmso_gwp_sig_expr_df[which(dihome_dmso_gwp_sig_expr_df$NES > 1),], "dihome_only_gsea_13uppathways_082021.csv", row.names = FALSE, quote = FALSE)
write.csv(dihome_dmso_gwp_sig_expr_df[which(dihome_dmso_gwp_sig_expr_df$NES < 1),], "dihome_only_gsea_7downpathways_082021.csv", row.names = FALSE, quote = FALSE)

# dihomepnut only data

dihomepnut_dmso_de_data_combined_entrez$fcsign <- sign(as.numeric(dihomepnut_dmso_de_data_combined_entrez$log2FC))
dihomepnut_dmso_de_data_combined_entrez$logfdr <- -log10(as.numeric(dihomepnut_dmso_de_data_combined_entrez$Unadjusted))
dihomepnut_dmso_de_data_combined_entrez$sig <- dihomepnut_dmso_de_data_combined_entrez$logfdr/dihomepnut_dmso_de_data_combined_entrez$fcsign
dihomepnut_dmso_gsea_sig_expr <- dihomepnut_dmso_de_data_combined_entrez$sig
names(dihomepnut_dmso_gsea_sig_expr) <- dihomepnut_dmso_de_data_combined_entrez$ENTREZID
dihomepnut_dmso_gsea_sig_expr_sorted <- sort(dihomepnut_dmso_gsea_sig_expr, decreasing = TRUE)

dihomepnut_dmso_gwp_sig_expr <- clusterProfiler::GSEA(
  dihomepnut_dmso_gsea_sig_expr_sorted,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

dihomepnut_dmso_gwp_sig_expr2 <- setReadable(dihomepnut_dmso_gwp_sig_expr, org.Hs.eg.db, keyType = "ENTREZID")

dihomepnut_dmso_gwp_sig_expr_df = data.frame(ID=dihomepnut_dmso_gwp_sig_expr2$ID,
                                             Description=dihomepnut_dmso_gwp_sig_expr2$Description,
                                             enrichmentScore=dihomepnut_dmso_gwp_sig_expr2$enrichmentScore,
                                             NES=dihomepnut_dmso_gwp_sig_expr2$NES,
                                             pvalue=dihomepnut_dmso_gwp_sig_expr2$pvalue,
                                             p.adjust=dihomepnut_dmso_gwp_sig_expr2$p.adjust,
                                             rank=dihomepnut_dmso_gwp_sig_expr2$rank,
                                             leading_edge=dihomepnut_dmso_gwp_sig_expr2$leading_edge,
                                             core_enrichment=dihomepnut_dmso_gwp_sig_expr2$core_enrichment
)

dihomepnut_dmso_gwp_sig_expr_df[which(dihomepnut_dmso_gwp_sig_expr_df$NES > 1),]
dihomepnut_dmso_gwp_sig_expr_df[which(dihomepnut_dmso_gwp_sig_expr_df$NES < 1),]

setwd("/Users/caraporsche/Desktop/diHOME RNAseq")
write.csv(dihomepnut_dmso_gwp_sig_expr_df, "dihomepnut_gsea_50pathways_082021.csv", row.names = FALSE, quote = FALSE)
write.csv(dihomepnut_dmso_gwp_sig_expr_df[which(dihomepnut_dmso_gwp_sig_expr_df$NES > 1),], "dihomepnut_gsea_39uppathways_082021.csv", row.names = FALSE, quote = FALSE)
write.csv(dihomepnut_dmso_gwp_sig_expr_df[which(dihomepnut_dmso_gwp_sig_expr_df$NES < 1),], "dihomepnut_gsea_11downpathways_082021.csv", row.names = FALSE, quote = FALSE)

## D) Which pathways are unique to dihome vs dihome+pnut?

dihome_pathways <- dihome_dmso_gwp_sig_expr_df[,2] # 20
dihomepnut_pathways <- dihomepnut_dmso_gwp_sig_expr_df[,2] # 50

shared_pathways <- intersect(dihome_pathways, dihomepnut_pathways) # 12
dihome_only_pathways <- dihome_pathways[!dihome_pathways %in% shared_pathways] # 8
dihomepnut_only_pathways <- dihomepnut_pathways[!dihomepnut_pathways %in% shared_pathways] # 38

#### 6) Visualize unique pathways ====

# visualize the networks with Cytoscape (dihome+pnut)

cytoscapePing()

dihomepnut_dmso_de_data_combined_entrez_mock <- dihomepnut_dmso_de_data_combined_entrez
dihomepnut_dmso_de_data_combined_entrez_mock$log2FC[dihomepnut_dmso_de_data_combined_entrez_mock$log2FC > 1] <- 1
dihomepnut_dmso_de_data_combined_entrez_mock$log2FC[dihomepnut_dmso_de_data_combined_entrez_mock$log2FC < -1] <- -1

loadTableData(dihomepnut_dmso_de_data_combined_entrez_mock, data.key.column = "ensembl_gene_id", table.key.column = "Ensembl")

min_dihomepnut = min(dihomepnut_dmso_de_data_combined_entrez_mock["log2FC"],na.rm=TRUE)
max_dihomepnut = max(dihomepnut_dmso_de_data_combined_entrez_mock["log2FC"],na.rm=TRUE)
abs_dihomepnut = max(abs(min_dihomepnut),max_dihomepnut)
data_values = c(-abs_dihomepnut,0,abs_dihomepnut)

display.brewer.all(length(data_values), colorblindFriendly=TRUE, type="div") # div,qual,seq,all
node.colors <- c(rev(brewer.pal(length(data_values), "RdBu")))

setNodeColorMapping("log2FC", data_values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")

lapply(dihomepnut_dmso_gwp_sig_expr_df[,1], function (x) {
  commandsRun(paste0('wikipathways import-as-pathway id=',x))
  loadTableData(dihomepnut_dmso_de_data_combined_entrez_mock, data.key.column = "ensembl_gene_id", table.key.column = "Ensembl")
  setNodeColorMapping("log2FC", data_values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
})

# visualize the networks with Cytoscape (dhome)

dihome_dmso_de_data_combined_entrez_mock <- dihome_dmso_de_data_combined_entrez
dihome_dmso_de_data_combined_entrez_mock$log2FC[dihome_dmso_de_data_combined_entrez_mock$log2FC > 1] <- 1
dihome_dmso_de_data_combined_entrez_mock$log2FC[dihome_dmso_de_data_combined_entrez_mock$log2FC < -1] <- -1

min_dihome = min(dihome_dmso_de_data_combined_entrez_mock["log2FC"],na.rm=TRUE)
max_dihome = max(dihome_dmso_de_data_combined_entrez_mock["log2FC"],na.rm=TRUE)
abs_dihome = max(abs(min_dihome),max_dihome)
data_values_dihome_only = c(-abs_dihome,0,abs_dihome)

display.brewer.all(length(data_values_dihome_only), colorblindFriendly=TRUE, type="div") # div,qual,seq,all
node.colors.dihome_only <- c(rev(brewer.pal(length(data_values_dihome_only), "RdBu")))

lapply(dihome_dmso_gwp_sig_expr_df[,1], function (x) {
  commandsRun(paste0('wikipathways import-as-pathway id=',x))
  loadTableData(dihome_dmso_de_data_combined_entrez_mock, data.key.column = "ensembl_gene_id", table.key.column = "Ensembl")
  setNodeColorMapping("log2FC", data_values_dihome_only, node.colors.dihome_only, default.color = "#FFFFFF", style.name = "WikiPathways")
})

saveSession("cytoscape_dihome_gsea_analysis_071621.cys") # 45 pathways for dihomepnut and 19 for dihome alone

## mouse RNA-seq analysis ====

library(biomaRt)
library(plyr)
library(limma)
library(edgeR)
library(WGCNA)
library(flashClust)
library(stringr)

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse")

treatments <- read.table("targets_mlung.txt", header = TRUE)
data <- read.table("mlung_rnaseq_counts.txt", header = TRUE)

pheno_data <- as.data.frame(t(as.data.frame(str_split(str_remove(treatments[,1], "_LA"), "_"))))
rownames(pheno_data) <- colnames(data)
colnames(pheno_data) <- c("DiHOME", "Allergen")
pheno_data$Phenotype <- paste(pheno_data[,1], pheno_data[,2], sep = "_")

#### 1) Remove X, Y, and M genes and extract mouse genes ====

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "mmusculus_gene_ensembl", 
                   mirror = "useast")
genes <- rownames(data)
G_list <- getBM(filters = "ensembl_gene_id_version", attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version", "external_gene_name", "description", "chromosome_name", "hgnc_id", "hgnc_symbol", "gene_biotype", "entrezgene_id"), values = genes, mart = mart)
G_list_unique <- G_list[!duplicated(G_list[,2]),] # all 55416 enmusg ids detected
G_protein <- G_list_unique[G_list_unique$gene_biotype == "protein_coding",] # 21885/55416 protein-coding
G_protein_auto <- G_protein[G_protein$chromosome_name == "1" | G_protein$chromosome_name == "2" | G_protein$chromosome_name == "3" | G_protein$chromosome_name == "4" | G_protein$chromosome_name == "5" | G_protein$chromosome_name == "6" | G_protein$chromosome_name == "7" | G_protein$chromosome_name == "8" | G_protein$chromosome_name == "9" | G_protein$chromosome_name == "10" | G_protein$chromosome_name == "11" | G_protein$chromosome_name == "12" | G_protein$chromosome_name == "13" | G_protein$chromosome_name == "14" | G_protein$chromosome_name == "15" | G_protein$chromosome_name == "16" | G_protein$chromosome_name == "17" | G_protein$chromosome_name == "18" | G_protein$chromosome_name == "19" | G_protein$chromosome_name == "20" | G_protein$chromosome_name == "21" | G_protein$chromosome_name == "22",] # 20718/21885 on autosomes
G_protein_auto_complete <- G_protein_auto[!duplicated(G_protein_auto$external_gene_name),] # 20683/20704 unique murine gene names

# 20,683 protein-coding and autosomal genes

# Rename ENSG to the mouse gene names

prot_aut_code <- intersect(G_protein_auto_complete$ensembl_gene_id_version, genes)
prot_aut_code_index <- match(prot_aut_code, genes)
ENSG_replace <- mapvalues(genes, from=G_protein_auto_complete$ensembl_gene_id_version, to=G_protein_auto_complete$external_gene_name)

data_renamed <- data[prot_aut_code_index,]
ENSG_replace_prot_aut <- ENSG_replace[prot_aut_code_index]

rownames(data_renamed) <- ENSG_replace_prot_aut

# write the mice biomart ids data frame to a mapping file

write.csv(G_protein_auto_complete, "/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/biomart_ids_mouse_dihome_070821.csv", row.names = FALSE, quote = FALSE)

#### 2) Pre-processing from the RNAseq analysis is as easy as 1-2-3 with limma, Glimma, and edgeR (subsetting lowly expressed genes) ====

## a) Check the density of log-RPM values for raw pre-filtered data

data_renamed[is.na(data_renamed)] <- 0

# How many genes have a count of 0 in all samples?  3165/20683 genes
table(rowSums(data_renamed==0)==12)

## b) Remove lowly expressed genes

cpm_raw <- cpm(data_renamed, log=FALSE)
lcpm_raw <- cpm(data_renamed, log=TRUE)

# Also remove lowly expressed genes (those absent in > 3 subjects)
exp_genes <- rownames(cpm_raw[rowSums(cpm_raw>1)>=3,]) # 14136/20683 genes present in at least 3 subjects
cpm_raw_exp <- cpm_raw[match(exp_genes, rownames(cpm_raw)),]
genecount_data_prot_aut_exp <- data_renamed[match(exp_genes, rownames(cpm_raw)),]
lcpm_exp <- cpm(genecount_data_prot_aut_exp, log=TRUE)

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Preprocessing")
pdf("DiHOME_mouse_062321.pdf", height = 6, width = 8)
par(mfrow=c(1,2))
plot(density(lcpm_raw[,1]), col = "Blue", lwd=2, las=2, main = "", xlab="", ylim=c(0,0.25))
title(main="A. Raw data", xlab="Log-CPM")
abline(v=0, lty=3)
for (i in 2:ncol(lcpm_raw)){
  den <- density(lcpm_raw[,i])
  lines(den$x, den$y, lwd=2, col = "Blue")
}
plot(density(lcpm_exp[,1]), col = "Purple", lwd=2, las=2, main = "", xlab="", ylim=c(0,0.25))
title(main="B. Filtered data", xlab="Log-CPM")
abline(v=0, lty=3)
for (i in 2:ncol(lcpm_exp)){
  den <- density(lcpm_exp[,i])
  lines(den$x, den$y, lwd=2, col = "Purple")
}
dev.off()

## c) Normalize the data

# To calculate the TMM normalization factors, create a DGElist using the edgeR package: 
dge.nl.nonorm <- DGEList(counts=genecount_data_prot_aut_exp)
# Perform TMM normalization using the calcNormFactors function:
dge.nl <- calcNormFactors(dge.nl.nonorm)

# Boxplots of the unnormalized and normalized values: 
setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Preprocessing")
pdf("DiHOME_mouse_preQC_norm_unnorm_boxplots_062321.pdf")
lcpm_nonorm <- cpm(dge.nl.nonorm, log=TRUE)
boxplot(lcpm_nonorm[,1:12], las=2, main="")
title(main="Unnormalized data, 1-12",ylab="Log-CPM")

lcpm_norm <- cpm(dge.nl, log=TRUE)
boxplot(lcpm_norm[,1:12], las=2, main="")
title(main="Normalized data (TMM), 1-12",ylab="Log-CPM")
dev.off()

# TMM normalization and voom

v.nl <- voom(dge.nl,design=NULL,plot=TRUE)
write.csv(v.nl$E, "/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Results/dihome_mouse_eh_cra_rnaseq_voom_norm_073121.csv")

#### 3) PCA =====

# Let's print out the cpm-adjusted file that includes the samples treated with LPS

v.nl <- read.csv("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Results/dihome_mouse_eh_cra_rnaseq_voom_norm_073121.csv", stringsAsFactors = FALSE, row.names = 1)
v.nl_geo <- v.nl
colnames(v.nl_geo) <- gsub("\\.", "-", colnames(v.nl))
write.table(v.nl_geo, "/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/GEO/GEO_Mouse_RNAseq_12.txt", quote = FALSE, sep = "\t", )

# A) Raw

DiHOME <- as.factor(pheno_data[,1])
Allergen <- as.factor(pheno_data[,2])
Phenotype <- as.factor(pheno_data[,3])

sum.PC <- prcomp(t(v.nl$E), scale=FALSE, center=TRUE)
sumsum <- summary(sum.PC)
prop.var<-sumsum$importance
prop.var[1:3,1:10] 

pc_list<-c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
m.covars<-cbind.data.frame(DiHOME, Allergen, Phenotype)

pval.pca1=matrix(ncol=ncol(m.covars), nrow=10)
colnames(pval.pca1)=colnames(m.covars)
rownames(pval.pca1)=paste(pc_list, " ", "(", round(prop.var[2,1:length(pc_list)], 3), ")", sep = "")

for(i in 1:ncol(m.covars))
{
  for(j in 1:length(pc_list))
  {
    data1= lm(sum.PC$x[,j]~m.covars[,i])
    pval.pca1[j,i]=anova(data1)$'Pr(>F)'[1]
  }
}

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Preprocessing")
write.csv(pval.pca1, file = "PCA_Raw_062321.csv", quote=FALSE)

ggplot(pheno_data, aes(x=sum.PC$x[,1], y=sum.PC$x[,2], color=Phenotype)) +
  geom_point()

#### 4) Identify differentially expressed genes between treatments ====

voom_norm <- v.nl$E

pheno_data$pheno <- paste(pheno_data[,1], pheno_data[,2], sep = "_")

# all mice are different; treated with...

design <- model.matrix(~ 0 + as.factor(paste(pheno_data[,1], pheno_data[,2], sep = "_")))
colnames(design) <- c("EH_CRA", "EH", "CRA", "Vehicle")
lmfit <- lmFit(voom_norm, design)
design_fit <- makeContrasts(EH-Vehicle, CRA-Vehicle, EH_CRA-Vehicle, EH_CRA-EH, EH_CRA-CRA, levels = design)
fit_contrasts <- contrasts.fit(lmfit, design_fit)
fit_contrasts_ebayes <- eBayes(fit_contrasts)
log2fc_data <- topTable(fit_contrasts_ebayes, number = nrow(fit_contrasts_ebayes), sort.by="none")

pvals <- fit_contrasts_ebayes$p.value
pvals_adj <- apply(fit_contrasts_ebayes$p.value, 2, function(x) p.adjust(x, "fdr"))
sum(pvals_adj[,1] < 0.05) # EH vs vehicle; 2 genes
sum(pvals_adj[,2] < 0.05) # CRA vs vehicle; 480 genes
sum(pvals_adj[,3] < 0.05) # EH_CRA vs vehicle; 614 genes
sum(pvals_adj[,4] < 0.05) # EH_CRA vs EH; 514 genes
sum(pvals_adj[,5] < 0.05) # EH_CRA vs CRA; 0 genes

# we only care about the first three tests because we're using a venn diagram to identify EH+CRA and CRA-specific genes...

venn_data <- cbind(log2fc_data[,1], pvals[,1], pvals_adj[,1], log2fc_data[,2], pvals[,2], pvals_adj[,2], log2fc_data[,3], pvals[,3], pvals_adj[,3])
colnames(venn_data) <- c("EH_vehicle_log2FC", "EH_vehicle_p", "EH_vehicle_p_adj", "CRA_vehicle_log2FC", "CRA_vehicle_p", "CRA_vehicle_p_adj", "EH_CRA_vehicle_log2FC", "EH_CRA_vehicle_p", "EH_CRA_vehicle_p_adj")

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Results")
write.csv(venn_data, "dihome_mouse_venn_data_3compare_070921.csv", quote = FALSE)

eh_vs_vehicle_genes <- rownames(venn_data)[venn_data[,3] < 0.05]
cra_vs_vehicle_genes <- rownames(venn_data)[venn_data[,6] < 0.05]
eh_cra_vs_vehicle_genes <- rownames(venn_data)[venn_data[,9] < 0.05]

# need to print the data to supplementary tables for significant DE genes for each comparison

eh_vs_vehicle_data <- venn_data[rownames(venn_data) %in% eh_vs_vehicle_genes, 1:3]
cra_vs_vehicle_data <- venn_data[rownames(venn_data) %in% cra_vs_vehicle_genes, 4:6]
eh_cra_vs_vehicle_data <- venn_data[rownames(venn_data) %in% eh_cra_vs_vehicle_genes, 7:9]

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Results")
write.csv(eh_vs_vehicle_data, "dihome_mouse_3eh_vs_vehicle_082421.csv", quote = FALSE)
write.csv(cra_vs_vehicle_data, "dihome_mouse_cra_vs_vehicle_082421.csv", quote = FALSE)
write.csv(eh_cra_vs_vehicle_data, "dihome_mouse_3eh_cra_vs_vehicle_082421.csv", quote = FALSE)

intersect(eh_vs_vehicle_genes, cra_vs_vehicle_genes)
intersect(eh_vs_vehicle_genes, eh_cra_vs_vehicle_genes)
intersect(cra_vs_vehicle_genes, eh_cra_vs_vehicle_genes)

intersect(intersect(eh_vs_vehicle_genes, cra_vs_vehicle_genes), eh_cra_vs_vehicle_genes) # 2 in all three

length(!intersect(intersect(eh_vs_vehicle_genes, cra_vs_vehicle_genes), eh_cra_vs_vehicle_genes) %in% eh_vs_vehicle_genes)
sum(!intersect(intersect(eh_vs_vehicle_genes, cra_vs_vehicle_genes), eh_cra_vs_vehicle_genes) %in% cra_vs_vehicle_genes)
sum(!intersect(intersect(eh_vs_vehicle_genes, cra_vs_vehicle_genes), eh_cra_vs_vehicle_genes) %in% eh_cra_vs_vehicle_genes)

sum(!cra_vs_vehicle_genes %in% intersect(cra_vs_vehicle_genes, eh_cra_vs_vehicle_genes))
sum(!eh_cra_vs_vehicle_genes %in% intersect(cra_vs_vehicle_genes, eh_cra_vs_vehicle_genes))

# extract and print the gene lists

cra_only_genes <- cra_vs_vehicle_genes[!cra_vs_vehicle_genes %in% intersect(cra_vs_vehicle_genes, eh_cra_vs_vehicle_genes)]
eh_cra_only_genes <- eh_cra_vs_vehicle_genes[!eh_cra_vs_vehicle_genes %in% intersect(cra_vs_vehicle_genes, eh_cra_vs_vehicle_genes)]

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Results")
write.csv(cra_only_genes, "CRA_vs_vehicle_only_genes_070721.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(eh_cra_only_genes, "EH_CRA_vs_vehicle_only_genes_070721.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)

#### 5) Perform pathway analysis ====

## A) Installation process

library(rWikiPathways)

load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}

cytoscapePing() #this will tell you if you're able to successfully connect to Cytoscape or not

installApp('WikiPathways') 
installApp('CyTargetLinker') 
installApp('stringApp') 
installApp('enrichmentMap')

## B) Perform enrichment (need to do for both [EH+CRA] and for [CRA])

#### combined data

G_protein_auto_complete_final <- G_protein_auto_complete[,c("ensembl_gene_id", "external_gene_name", "entrezgene_id")]
venn_data_final <- cbind(rownames(venn_data), venn_data)
colnames(venn_data_final) <- c("external_gene_name", colnames(venn_data))

data_combined <- merge(G_protein_auto_complete_final, venn_data_final, by = "external_gene_name", all = FALSE, sort = FALSE)

entrez_ids <- bitr(data_combined[,2],fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Mm.eg.db) # 0.63% of input gene IDs are fail to map
data_combined_entrez_keep <- merge(data_combined, entrez_ids, by.x = "ensembl_gene_id", by.y = "ENSEMBL", all = FALSE, sort = FALSE) # retain 14112 genes

eh_only_full_data <- data_combined_entrez_keep[,c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "EH_vehicle_log2FC", "EH_vehicle_p", "EH_vehicle_p_adj")]
cra_only_full_data <- data_combined_entrez_keep[,c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "CRA_vehicle_log2FC", "CRA_vehicle_p", "CRA_vehicle_p_adj")]
eh_cra_only_full_data <- data_combined_entrez_keep[,c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "EH_CRA_vehicle_log2FC", "EH_CRA_vehicle_p", "EH_CRA_vehicle_p_adj")]

#### extract entrez to gene/pathway names

wp.mm.gmt <- rWikiPathways::downloadPathwayArchive(organism="Mus musculus", format = "gmt")
listOrganisms()

wp2gene <- readPathwayGMT(wp.mm.gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
wpid2gene
wpid2name

#### GSEA analysis (cra)

cra_only_full_data$fcsign <- sign(as.numeric(cra_only_full_data$CRA_vehicle_log2FC))
cra_only_full_data$logfdr <- -log10(as.numeric(cra_only_full_data$CRA_vehicle_p))
cra_only_full_data$sig <- cra_only_full_data$logfdr/cra_only_full_data$fcsign
cra_only_gsea_sig_expr <- cra_only_full_data$sig
names(cra_only_gsea_sig_expr) <- cra_only_full_data$entrezgene_id
cra_only_gsea_sig_expr_sorted <- sort(cra_only_gsea_sig_expr, decreasing = TRUE)

cra_only_gwp_sig_expr <- GSEA(
  cra_only_gsea_sig_expr_sorted,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name
)

cra_only_gwp_sig_expr2 <- setReadable(cra_only_gwp_sig_expr, org.Mm.eg.db, keyType = "ENTREZID")

cra_only_gwp_sig_expr_df <- data.frame(ID=cra_only_gwp_sig_expr2$ID,
                                       Description=cra_only_gwp_sig_expr2$Description,
                                       enrichmentScore=cra_only_gwp_sig_expr2$enrichmentScore,
                                       NES=cra_only_gwp_sig_expr2$NES,
                                       pvalue=cra_only_gwp_sig_expr2$pvalue,
                                       p.adjust=cra_only_gwp_sig_expr2$p.adjust,
                                       rank=cra_only_gwp_sig_expr2$rank,
                                       leading_edge=cra_only_gwp_sig_expr2$leading_edge,
                                       core_enrichment=cra_only_gwp_sig_expr2$core_enrichment
)

cra_only_gwp_sig_expr_df[which(cra_only_gwp_sig_expr_df$NES > 1),] # 39 pathways enriched for upregulated cockroach antigen genes
cra_only_gwp_sig_expr_df[which(cra_only_gwp_sig_expr_df$NES < -1),] # 12 pathways enriched for downregulated cockroach antigen genes

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Results")
write.csv(cra_only_gwp_sig_expr_df, "cra_only_gsea_51pathways_082021.csv", row.names = FALSE, quote = FALSE)

#### GSEA analysis (eh_cra)

eh_cra_only_full_data$fcsign <- sign(as.numeric(eh_cra_only_full_data$EH_CRA_vehicle_log2FC))
eh_cra_only_full_data$logfdr <- -log10(as.numeric(eh_cra_only_full_data$EH_CRA_vehicle_p))
eh_cra_only_full_data$sig <- eh_cra_only_full_data$logfdr/eh_cra_only_full_data$fcsign
eh_cra_only_gsea_sig_expr <- eh_cra_only_full_data$sig
names(eh_cra_only_gsea_sig_expr) <- eh_cra_only_full_data$entrezgene_id
eh_cra_only_gsea_sig_expr_sorted <- sort(eh_cra_only_gsea_sig_expr, decreasing = TRUE)

eh_cra_only_gwp_sig_expr <- GSEA(
  eh_cra_only_gsea_sig_expr_sorted,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name
)

eh_cra_only_gwp_sig_expr2 <- setReadable(eh_cra_only_gwp_sig_expr, org.Mm.eg.db, keyType = "ENTREZID")

eh_cra_only_gwp_sig_expr_df <- data.frame(ID=eh_cra_only_gwp_sig_expr2$ID,
                                       Description=eh_cra_only_gwp_sig_expr2$Description,
                                       enrichmentScore=eh_cra_only_gwp_sig_expr2$enrichmentScore,
                                       NES=eh_cra_only_gwp_sig_expr2$NES,
                                       pvalue=eh_cra_only_gwp_sig_expr2$pvalue,
                                       p.adjust=eh_cra_only_gwp_sig_expr2$p.adjust,
                                       rank=eh_cra_only_gwp_sig_expr2$rank,
                                       leading_edge=eh_cra_only_gwp_sig_expr2$leading_edge,
                                       core_enrichment=eh_cra_only_gwp_sig_expr2$core_enrichment
)

eh_cra_only_gwp_sig_expr_df[which(eh_cra_only_gwp_sig_expr_df$NES > 1),] # 35 pathways enriched for upregulated cockroach antigen genes
eh_cra_only_gwp_sig_expr_df[which(eh_cra_only_gwp_sig_expr_df$NES < -1),] # 5 pathways enriched for downregulated cockroach antigen genes

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Results")
write.csv(eh_cra_only_gwp_sig_expr_df, "eh_cra_gsea_40pathways_082021.csv", row.names = FALSE, quote = FALSE)

# which pathways are unique to cra only and eh_cra only?

cra_only_pathways <- cra_only_gwp_sig_expr_df[,2] # 51 pathways
eh_cra_only_pathways <- eh_cra_only_gwp_sig_expr_df[,2] # 40 pathways
cra_eh_cra_shared_pathways <- intersect(cra_only_pathways, eh_cra_only_pathways) # 35 pathways

cra_only_unique_pathways <- cra_only_pathways[!cra_only_pathways %in% unique(eh_cra_only_pathways)] # 16
eh_cra_only_unique_pathways <- eh_cra_only_pathways[!eh_cra_only_pathways %in% unique(cra_only_pathways)] # 5

# unique pathways

eh_cra_only_gwp_sig_expr_df_unique <- eh_cra_only_gwp_sig_expr_df[eh_cra_only_gwp_sig_expr_df[,2] %in% eh_cra_only_unique_pathways,]
cra_only_gwp_sig_expr_df_unique <- cra_only_gwp_sig_expr_df[cra_only_gwp_sig_expr_df[,2] %in% cra_only_unique_pathways,]

# shared pathways

shared_pathways_eh_cra_cra_only <- eh_cra_only_gwp_sig_expr_df[eh_cra_only_gwp_sig_expr_df[,2] %in% cra_eh_cra_shared_pathways,]

setwd("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse")
write.csv(eh_cra_only_gwp_sig_expr_df_unique, "EH_CRA_Mouse_Unique_Pathways_072621.csv", quote = FALSE, row.names = FALSE)
write.csv(cra_only_gwp_sig_expr_df_unique, "CRA_Mouse_Unique_Pathways_072621.csv", quote = FALSE, row.names = FALSE)
write.csv(shared_pathways_eh_cra_cra_only, "EH_CRA_and_CRA_Mouse_Shared_Pathways_072621.csv", quote = FALSE, row.names = FALSE)

read.csv("EH_CRA_Mouse_Unique_Pathways_072621.csv")

#### 6) Visualize eh_cra only pathways (5 of them) ====

cytoscapePing()

RCy3::commandsRun('wikipathways import-as-pathway id=WP37') 

eh_cra_only_full_data_mock <- eh_cra_only_full_data
eh_cra_only_full_data_mock$EH_CRA_vehicle_log2FC <- as.numeric(eh_cra_only_full_data_mock$EH_CRA_vehicle_log2FC)
eh_cra_only_full_data_mock$EH_CRA_vehicle_log2FC[eh_cra_only_full_data_mock$EH_CRA_vehicle_log2FC > 1] <- 1
eh_cra_only_full_data_mock$EH_CRA_vehicle_log2FC[eh_cra_only_full_data_mock$EH_CRA_vehicle_log2FC < -1] <- -1

loadTableData(eh_cra_only_full_data_mock, data.key.column = "ensembl_gene_id", table.key.column = "Ensembl")

min_eh_cra_only_full_data = min(as.numeric(eh_cra_only_full_data_mock[,"EH_CRA_vehicle_log2FC"]),na.rm=TRUE)
max_eh_cra_only_full_data = max(as.numeric(eh_cra_only_full_data_mock[,"EH_CRA_vehicle_log2FC"]),na.rm=TRUE)
abs_eh_cra_only_full_data = max(abs(min_eh_cra_only_full_data), max_eh_cra_only_full_data)
data.values = c(-abs_eh_cra_only_full_data,0,abs_eh_cra_only_full_data)

display.brewer.all(length(data.values), colorblindFriendly=TRUE, type="div") # div,qual,seq,all
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))

setNodeColorMapping("EH_CRA_vehicle_log2FC", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")

lapply(eh_cra_only_gwp_sig_expr_df_unique$ID, function (x) {
  commandsRun(paste0('wikipathways import-as-pathway id=',x))
  loadTableData(eh_cra_only_full_data_mock, data.key.column = "ensembl_gene_id", table.key.column = "Ensembl")
  setNodeColorMapping("EH_CRA_vehicle_log2FC", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
  
})

#### 7) Need to make a heatmap of the gene expression for the top 100 genes of EH_CRA vs vehicle ====

# Need to load the differentially expressed genes for the test comparing EH+CRA vs vehicle and the gene expression matrix

venn_data_eh_cra_only <- venn_data[venn_data[,9] < 0.05 & venn_data[,3] > 0.05 & venn_data[,6] > 0.05,]
venn_data_eh_cra_ordered <- venn_data_eh_cra_only[order(venn_data_eh_cra_only[,8]),]
top100_eh_cra <- rownames(venn_data_eh_cra_ordered[1:100,])

voom_norm_top100_eh_cra <- voom_norm[rownames(voom_norm) %in% top100_eh_cra,]
colnames(voom_norm_top100_eh_cra) <- substr(colnames(voom_norm_top100_eh_cra), 1, 3)

# plot the heatmap

pheno_table <- as.data.frame(c(rep("3EH + CRA", 3), rep("CRA", 3), rep("3EH", 3), rep("Vehicle", 3)))
rownames(pheno_table) <- colnames(voom_norm_top100_eh_cra)
colnames(pheno_table) <- c("Treatment")
pheno_table[,1] <- factor(pheno_table[,1], levels = c("Vehicle", "3EH", "CRA", "3EH + CRA"))

heatmap_top100 <- pheatmap::pheatmap(voom_norm_top100_eh_cra, annotation_col = pheno_table, cluster_cols = T, cluster_rows = T)
hc_top100 <- hclust(dist(t(voom_norm_top100_eh_cra)))
hc_top100_neworder <- hc_top100$order[c(3:5,1:2,6,10:12,7:9)]
hc_top100$order <- hc_top100_neworder
hc_top100$labels <- paste0("S", hc_top100_neworder)

heatmap_top100_reordered <- pheatmap::pheatmap(voom_norm_top100_eh_cra, 
                                               annotation_col = pheno_table, 
                                               cluster_cols = hc_top100, 
                                               cluster_rows = T,
                                               labels_col = hc_top100$labels)

pdf("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Results/boxplot_top100_ehcra_genes_treatment_080321.pdf", width = 6, height = 7)
heatmap_top100_reordered
dev.off()


#### 8) PCA (can'tn do) of expression differences between groups ====

voom_norm_top100_eh_cra
pheno_table

sum.PC <- prcomp(t(voom_norm_top100_eh_cra), scale=FALSE, center=TRUE)
sumsum <- summary(sum.PC)
prop.var<-sumsum$importance
prop.var[1:3,1:10] 

pc_list<-c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

pval.pca1=matrix(ncol=ncol(pheno_table), nrow=10)
colnames(pval.pca1)=colnames(pheno_table)
rownames(pval.pca1)=paste(pc_list, " ", "(", round(prop.var[2,1:length(pc_list)], 3), ")", sep = "")

for(i in 1:ncol(pheno_table))
{
  for(j in 1:length(pc_list))
  {
    data1= lm(sum.PC$x[,j]~pheno_table[,i])
    pval.pca1[j,i]=anova(data1)$'Pr(>F)'[1]
  }
}

pca_plot_treatment_top100 <- ggplot(pheno_table, aes(x=sum.PC$x[,1], y=sum.PC$x[,2], color=Treatment)) +
  geom_point(size = 3) + theme_bw() + xlab("PC1 (65.3%)") + ylab("PC2 (16.6%)") + theme(legend.position="bottom")

pdf("/Users/kmagnaye/Documents/Lynch_Lab/12-13DiHOME_macs/dihome_mouse/Results/pca_plot_top100_ehcra_genes_treatment_101221.pdf", width = 4.5, height = 4.5)
pca_plot_treatment_top100
dev.off()

# Let's try PCoA using Bray distance

otu_mock_top100 <- otu_table(voom_norm_top100_eh_cra, taxa_are_rows = TRUE)

otu_mock_top100_bray <- phyloseq::distance(otu_mock_top100, method = "bray")
vegan::adonis2(stool_prefmt_beta_bray ~ Responder, data = data.frame(stool_16s_vst_baseline@sam_data)) # P = 0.52, R2 = 0.055

pcoa_stool_postfmt <- ordinate(stool_16s_vst_followup, method="PCoA", distance="unifrac")
plot_ordination(stool_16s_vst_followup, pcoa_stool_postfmt, type="samples", color="Responder", axes = c(1,2)) + scale_color_manual(labels = c("Non-responder", "Responder"), values = c("#FF6600", "#6666FF")) + theme(legend.position = "none")




