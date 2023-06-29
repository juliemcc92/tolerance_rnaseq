library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(fgsea)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(glue)
library(purrr)
library(biomaRt)
library(VennDiagram)
library(RColorBrewer)
library(PCAtools)
library(tayloRswift)

switch(Sys.info()[['user']],
       julie = {fig.file.path <- "C:/Users/julie/Box/tolerance_bulk/figures"
       data.file.path <- "C:/Users/julie/Box/tolerance_bulk/rawdata"},
       stop("I don't recognize your username, type Sys.info() to find out what it is.")
)

#import genes.results files
samples <- c("24_no_stim", "25_bb_restim", "26_pam_restim", "27_bb_6hr", "28_pam_6hr", "29_no_stim", "30_bb_restim", "31_pam_restim", "32_bb_6hr", "33_pam_6hr")
all_results <- lapply(samples, function(x) {
  df <- read.csv(glue("{data.file.path}/{x}.genes.results"), sep = "\t")
  df <- df[, c("gene_id", "expected_count")]
  colnames(df) <- c("ensembl_gene_id", glue("expected_count_{x}"))  
  return(df)
})

count_matrix <- reduce(all_results, left_join, by= "ensembl_gene_id")
write.csv(count_matrix, "C:/Users/julie/Box/tolerance_bulk/rawdata/count_matrix_tolerance.csv")


# Gets gene information fron ensembl
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

# Matches all ensembl gene IDs to mgi symbols
ref_df <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=count_matrix$ensembl_gene_id,mart= mart)

count_matrix <- left_join(count_matrix, ref_df, by = "ensembl_gene_id")

#get rid of duplicate mgi symbol, make mgi symbol row names
count_matrix <- distinct(count_matrix, mgi_symbol, .keep_all= TRUE)
count_matrix <- count_matrix[count_matrix$mgi_symbol != "",]
count_matrix <- count_matrix[!is.na(count_matrix$mgi_symbol),]
rownames(count_matrix) <- count_matrix$mgi_symbol
count_matrix$ensembl_gene_id <- NULL
count_matrix$mgi_symbol <- NULL
count_matrix <- round(count_matrix)

#load metadata
tolerance_metadata <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/tolerance_metadata_2exps.csv", header= TRUE, sep= ",")
rownames(tolerance_metadata)<-tolerance_metadata$X
tolerance_metadata$X <- NULL
all(colnames(count_matrix) == rownames(tolerance_metadata))

#create bb stim count matrix
bb_stim_counts <- count_matrix[, grepl("bb_6hr|no", colnames(count_matrix))]

#load bb stim metadata
bb_stim_metadata <- tolerance_metadata[grepl("bb_6hr|no", rownames(tolerance_metadata)), ]


#DESEQ2 FOR BB STIM AND RESTIM
#create unstim to 6hr deseq data set for bb
bb_stim_dds <- DESeqDataSetFromMatrix(countData = bb_stim_counts, colData = bb_stim_metadata, design = ~ stim)
keep <- rowSums(counts(bb_stim_dds)) >=10
bb_stim_dds <- bb_stim_dds[keep, ]

#Make unstim the reference
bb_stim_dds$stim <- relevel(bb_stim_dds$stim, ref = "none")

#DEseq analysis unstim to 6hr bb
bb_stim_dds <- DESeq(bb_stim_dds)
bb_stim_res <- results(bb_stim_dds)

#make unstim to 6hr deseq results into df
bb_stim_res_data <- as.data.frame(bb_stim_res)

#get rid of blank padj values
bb_stim_res_data <- bb_stim_res_data[!is.na(bb_stim_res_data$padj), ]

#export deseq results for unstim to 6hr
write.csv(bb_stim_res_data, "C:/Users/julie/Box/tolerance_bulk/data/2_exps_bb_stim_resdata.csv")

#DESEQ2 FOR BB ACUTE TO RESTIM
#create bb restim count matrix
bb_restim_counts <- count_matrix[, grepl("bb_6hr|bb_restim", colnames(count_matrix))]

#load bb restim metadata
bb_restim_metadata <- tolerance_metadata[grepl("bb_6hr|bb_restim", rownames(tolerance_metadata)), ]
bb_restim_metadata$time <- factor(bb_restim_metadata$time)

#create stim to restim deseq data set for bb
bb_restim_dds <- DESeqDataSetFromMatrix(countData = bb_restim_counts, colData = bb_restim_metadata, design = ~ time)
keep <- rowSums(counts(bb_restim_dds)) >=10
bb_restim_dds <- bb_restim_dds[keep, ]

#Make stim the reference
bb_restim_dds$time <- relevel(bb_restim_dds$time, ref = 1)

#DEseq analysis stim to restim bb
bb_restim_dds <- DESeq(bb_restim_dds)
bb_restim_res <- results(bb_restim_dds)

#make stim to restim deseq results into df
bb_restim_res_data <- as.data.frame(bb_restim_res)

#get rid of blank padj values
bb_restim_res_data <- bb_restim_res_data[!is.na(bb_restim_res_data$padj), ]

#export deseq results
write.csv(bb_restim_res_data, "C:/Users/julie/Box/tolerance_bulk/data/2_exps_bb_restim_resdata.csv")


#combine BB stim and restim deseq files
#load files
bb_stim_labelled_2exps <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/2_exps_bb_stim_resdata_labelled.csv", sep = ",", stringsAsFactors = FALSE)
bb_restim_labelled_2exps <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/2_exps_bb_restim_resdata_labelled.csv", sep = ",", stringsAsFactors = FALSE)

#merge stim and restim and keep all rows but combine those that have overlapping genes
bb_merge_all_2exps <- full_join(bb_stim_labelled_2exps, bb_restim_labelled_2exps, by = "gene", all = TRUE)
write.csv(bb_merge_all_2exps, "C:/Users/julie/Box/tolerance_bulk/data/2exp_bb_merge_all.csv")

#categories for quadrants- FIX THIS CODE FOR 2 EXPS
bb_tolerizable_2exps <- filter(bb_merge_all_2exps, stim_log2FoldChange > 0.5, restim_log2FoldChange < -0.5, stim_padj <=0.1, restim_padj <=0.1)
write.csv(bb_tolerizable_2exps, "C:/Users/julie/Box/tolerance_bulk/data/bb_tolerizable_2exps.csv")

bb_hyperinduced_2exps <- filter(bb_merge_all_2exps, stim_log2FoldChange > 0.5, restim_log2FoldChange > 0.5, stim_padj <=0.1, restim_padj <=0.1)
write.csv(bb_hyperinduced_2exps, "C:/Users/julie/Box/tolerance_bulk/data/bb_hyperinduced_2exps.csv")

bb_neutral_2exps <- filter(bb_merge_all_2exps, stim_log2FoldChange <0.5 & stim_log2FoldChange >-0.5, restim_log2FoldChange <0.5 & restim_log2FoldChange >-0.5, stim_padj <=0.1, restim_padj<=0.1)
write.csv(bb_neutral_2exps, "C:/Users/julie/Box/tolerance_bulk/data/bb_neutral_2exps.csv")

bb_nontolerizable_2exps <- filter(bb_merge_all_2exps, stim_log2FoldChange >0.5, restim_log2FoldChange <0.5 & restim_log2FoldChange >-0.5, stim_padj <=0.1, restim_padj <=0.1)
write.csv(bb_nontolerizable_2exps, "C:/Users/julie/Box/tolerance_bulk/data/bb_nontolerized_2exps.csv")

bb_down_down_2exps <- filter(bb_merge_all_2exps, stim_log2FoldChange < -0.5, restim_log2FoldChange < -0.5 , stim_padj <=0.1, restim_padj <=0.1)
write.csv(bb_down_down_2exps, "C:/Users/julie/Box/tolerance_bulk/data/bb_down_down_2exps.csv")

bb_down_up_2exps <-filter(bb_merge_all_2exps, stim_log2FoldChange < -0.5, restim_log2FoldChange >0.5 & restim_log2FoldChange >-0.5, stim_padj <=0.1, restim_padj <=0.1)
write.csv(bb_down_up_2exps, "C:/Users/julie/Box/tolerance_bulk/data/bb_down_up_2exps.csv")


#STIM VOLCANO PLOT-NEGATIVE REGULATORS
# Need to set the zeros to something more than zero
bb_stim_res_data[bb_stim_res_data$pvalue == 0,]$pvalue <- sort(unique(bb_stim_res_data$pvalue))[2]
bb_stim_res_data$gene <- rownames(bb_stim_res_data)
# Label top 10 for now
up_label <- bb_stim_res_data %>%
  dplyr::filter(log2FoldChange > 0) %>%
  arrange(pvalue, desc(log2FoldChange)) %>%
  head(15) %>%
  .$gene
down_label <- bb_stim_res_data %>%
  dplyr::filter(log2FoldChange < 0) %>%
  arrange(pvalue, log2FoldChange) %>%
  head(15) %>%
  .$gene
all_label <- c(up_label, down_label)

all_label <- c("Atf3", "Peli3", "Nlrc5", "Cd47", "Irak3", "Lrrk2", "Tank", "Tollip", "Dhx58", "Il10", "Prkcd", "Cryab", "Il15", "Tyk2", "Socs1", "Nfe2l2", "Acod1")


ggplot(bb_stim_res_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(data = bb_stim_res_data %>%
               dplyr::filter(abs(log2FoldChange < 0.5) | padj > 0.05)) +
  geom_point(data = bb_stim_res_data %>%
               dplyr::filter(log2FoldChange > 0.5, padj < 0.05), color = "red") +
  geom_point(data = bb_stim_res_data %>%
               dplyr::filter(log2FoldChange < -0.5, padj < 0.05), color = "blue") +
  geom_label_repel(data= bb_stim_res_data %>%
                    dplyr::filter(gene %in% all_label), 
                  aes(label = gene), min.segment.length = 0.2) +
  xlab("log2FoldChange Stim") +
  ylab("-log10(pvalue) Stim") +
  ggtitle("Described Negative Regulators of Inflammation") +
  theme_classic(base_size = 20)
ggsave(paste(fig.file.path, "stim_vs_unstim_2_exps_volcano_small_gene_list.pdf", sep = "/"), width = 8,height = 6,
       useDingbats = FALSE)

# Can we fit more genes to label
extra_label_genes <- c("Traf3", "Tmem173", "Nlrp6", "Stat3", "Dusp1", "Cyld", "Nlrp12", "Smpdl3b", "Rad23a", "Foxo3", "Tnfaip8l2", "Dok3", "Nlrp6", "Prmt1", "Usp4", "Elf1", "Stap2", "Traf1", "Pkn1", "Nxn", "Smad7", "Smad6")


ggplot(bb_stim_res_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(data = bb_stim_res_data %>%
               dplyr::filter(abs(log2FoldChange < 0.5) | padj > 0.05)) +
  geom_point(data = bb_stim_res_data %>%
               dplyr::filter(log2FoldChange > 0.5, padj < 0.05), color = "red") +
  geom_point(data = bb_stim_res_data %>%
               dplyr::filter(log2FoldChange < -0.5, padj < 0.05), color = "blue") +
  geom_label_repel(data= bb_stim_res_data %>%
                     dplyr::filter(gene %in% all_label), 
                   aes(label = gene), min.segment.length = 0.2) +
  geom_label_repel(data= bb_stim_res_data %>%
                     dplyr::filter(gene %in% extra_label_genes), 
                   aes(label = gene), min.segment.length = 0.2) +
  ggtitle("Stim vs. unstim : 2 exps") +
  theme_classic(base_size = 20)
ggsave(paste(fig.file.path, "stim_vs_unstim_2_exps_volcano_large_gene_list.pdf", sep = "/"), width = 8,height = 6,
       useDingbats = FALSE)


#4 quadrant plot- FIX THIS
#scatter plot to show all bb genes and highlight different categories
bb_all_genes_2exps <- ggplot(bb_merge_all_2exps, aes(x = stim_log2FoldChange, y = restim_log2FoldChange)) +
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj >0.1 | bb_merge_all_2exps$restim_padj >0.1 | 
                                   bb_merge_all_2exps$stim_log2FoldChange > -0.5 & bb_merge_all_2exps$stim_log2FoldChange <0.5 |
                                   bb_merge_all_2exps$restim_log2FoldChange > -0.5 & bb_merge_all_2exps$restim_log2FoldChange <0.5, ], color = "grey")+
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj < 0.1 & bb_merge_all_2exps$restim_padj <0.1 &  
                                   bb_merge_all_2exps$stim_log2FoldChange > 0.5 & bb_merge_all_2exps$restim_log2FoldChange > 0.5, ], color = "green")+
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                   bb_merge_all_2exps$stim_log2FoldChange > 0.5 & bb_merge_all_2exps$restim_log2FoldChange < -0.5, ], color = "red")+                                
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                   bb_merge_all_2exps$stim_log2FoldChange < -0.5 & bb_merge_all_2exps$restim_log2FoldChange > 0.5, ], color = "orange")+                   
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                    bb_merge_all_2exps$stim_log2FoldChange < -0.5 & bb_merge_all_2exps$restim_log2FoldChange < -0.5, ], color = "purple")+ 
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                   bb_merge_all_2exps$stim_log2FoldChange < -0.5 &
                                     bb_merge_all_2exps$restim_log2FoldChange > -0.5 & bb_merge_all_2exps$restim_log2FoldChange <0.5, ], color = "blue")+
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                   bb_merge_all_2exps$stim_log2FoldChange >0.5 &
                                   bb_merge_all_2exps$restim_log2FoldChange > -0.5 & bb_merge_all_2exps$restim_log2FoldChange <0.5, ], color = "blue")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic(base_size = 20)+
  xlab("Log FC stimulation") + ylab("Log FC restimulation")

#geom_text_repel(data = bb_merge_all[bb_merge_all$padj_stim <= 0.1 & bb_merge_all$padj_restim <= 0.1 & 
#bb_merge_all$log2FoldChange_stim >5 & bb_merge_all$log2FoldChange_restim< 0, ], aes(label = gene))+

ggsave(paste(fig.file.path, "bb_quadrants_2exps.pdf", sep = "/"), width = 10, height = 10, useDingbats = FALSE )



### NEW Version,get rid of the blue category ###
ggplot(bb_merge_all_2exps, aes(x = stim_log2FoldChange, y = restim_log2FoldChange)) +
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj >0.1 | bb_merge_all_2exps$restim_padj >0.1 | 
                                         bb_merge_all_2exps$stim_log2FoldChange > -0.5 & bb_merge_all_2exps$stim_log2FoldChange <0.5 |
                                         bb_merge_all_2exps$restim_log2FoldChange > -0.5 & bb_merge_all_2exps$restim_log2FoldChange <0.5, ], color = "grey")+
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj < 0.1 & bb_merge_all_2exps$restim_padj <0.1 &  
                                         bb_merge_all_2exps$stim_log2FoldChange > 0.5 & bb_merge_all_2exps$restim_log2FoldChange > 0.5, ], color = "green")+
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                         bb_merge_all_2exps$stim_log2FoldChange > 0.5 & bb_merge_all_2exps$restim_log2FoldChange < -0.5, ], color = "magenta")+                                
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                         bb_merge_all_2exps$stim_log2FoldChange < -0.5 & bb_merge_all_2exps$restim_log2FoldChange > 0.5, ], color = "orange")+                   
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                         bb_merge_all_2exps$stim_log2FoldChange < -0.5 & bb_merge_all_2exps$restim_log2FoldChange < -0.5, ], color = "purple")+ 
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                         bb_merge_all_2exps$stim_log2FoldChange < -0.5 &
                                         bb_merge_all_2exps$restim_log2FoldChange > -0.5 & bb_merge_all_2exps$restim_log2FoldChange <0.5, ], color = "grey")+
  geom_point(data = bb_merge_all_2exps[bb_merge_all_2exps$stim_padj <0.1 & bb_merge_all_2exps$restim_padj <0.1 & 
                                         bb_merge_all_2exps$stim_log2FoldChange >0.5 &
                                         bb_merge_all_2exps$restim_log2FoldChange > -0.5 & bb_merge_all_2exps$restim_log2FoldChange <0.5, ], color = "grey")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic(base_size = 30)+
  xlab("Log FC stimulation") + ylab("Log FC restimulation")
ggsave(paste(fig.file.path, "bb_quadrants_2exps_magenta.pdf", sep = "/"), width = 10, height = 10, useDingbats = FALSE )


#ggplots for each subpopulation of interest
tol_genes_to_label <- c("Il12b", "Gfi1", "Il27", "Ifit1", "Mx2", "Nox1", "Cxcl10", "Ifit2", "Tnf", "Isg15",
                        "Nod2", "Ccrl2", "Il6", "Lipg", "Ifi205", "Cxcl9", "Nlrp3", "Tlr2", "Il1b", "Il18")
plot_bb_tolerizable <-ggplot(bb_tolerizable_2exps, aes(x = stim_log2FoldChange, y = restim_log2FoldChange)) +
  geom_point(data = bb_tolerizable_2exps, color = "red")+
  geom_label_repel(data = bb_tolerizable_2exps[bb_tolerizable_2exps$stim_log2FoldChange >2 & bb_tolerizable_2exps$restim_log2FoldChange < -2, ], aes(label = gene))+
  theme_classic(base_size = 30)+
  xlab("Log FC stimulation") + ylab("Log FC restimulation")
ggsave(paste(fig.file.path, "bb_tolerizable.pdf", sep = "/"), width = 10, height = 10, useDingbats = FALSE )

#with labels
plot_bb_tolerizable <-ggplot(bb_tolerizable_2exps, aes(x = stim_log2FoldChange, y = restim_log2FoldChange)) +
  geom_point(data = bb_tolerizable_2exps, color = "magenta")+
  geom_label_repel(data = bb_tolerizable_2exps %>%
                    dplyr::filter(gene %in% tol_genes_to_label), aes(label = gene), 
                  size= 6, min.segment.length = 0.2)+
  theme_classic(base_size = 30)+
  xlab("Log FC stimulation") + ylab("Log FC restimulation")
ggsave(paste(fig.file.path, "bb_tolerizable_specific_genes_magenta.pdf", sep = "/"), plot_bb_tolerizable, width = 8, height = 8, useDingbats = FALSE )


plot_bb_hyperinduced <-ggplot(bb_hyperinduced_2exps, aes(x = stim_log2FoldChange, y = restim_log2FoldChange)) +
  geom_point(data = bb_hyperinduced_2exps, color = "green")+
  geom_text_repel(data = bb_hyperinduced_2exps[bb_hyperinduced_2exps$stim_log2FoldChange >2 & bb_hyperinduced_2exps$restim_log2FoldChange > 2, ], aes(label = gene))+
  theme_classic(base_size = 20)+
  xlab("Log FC stimulation") + ylab("Log FC restimulation")
ggsave(paste(fig.file.path, "bb_hyperinduced.pdf", sep = "/"), width = 10, height = 10, useDingbats = FALSE )

# Now one for hyperinduced with specific labels
hyp_ind_genes_to_label <- c("Ly6i", "Nqo1", "Prdx6", "Prdx5", "Prdx1", "Hmox1", "Gclm", "Il10", "Slamf1", "Ereg", "Ifitm1", "Il36g", "Socs1", "Csf3")
plot_bb_hyperinduced <-ggplot(bb_hyperinduced_2exps, aes(x = stim_log2FoldChange, y = restim_log2FoldChange)) +
  geom_point(data = bb_hyperinduced_2exps, fill = "green", color = "black", shape = 21)+
  geom_label_repel(data = bb_hyperinduced_2exps %>%
                    dplyr::filter(gene %in% hyp_ind_genes_to_label), aes(label = gene), 
                  size = 6, min.segment.length = 0.2)+
  theme_classic(base_size = 30)+
  xlab("Log FC stimulation") + ylab("Log FC restimulation")
ggsave(paste(fig.file.path, "bb_hyperinduced_specific_genes_labeled_boxes.pdf", sep = "/"), width = 8, height = 8, useDingbats = FALSE )

#plot hypersuppressed
#just having top genes labelled via geom_point

plot_bb_down_down <-ggplot(bb_down_down_2exps, aes(x = stim_log2FoldChange, y = restim_log2FoldChange)) +
  geom_point(data = bb_down_down_2exps, fill = "purple", color = "black", shape = 21) +
  geom_label_repel(data = bb_down_down_2exps[bb_down_down_2exps$stim_log2FoldChange < -2 & bb_down_down_2exps$restim_log2FoldChange < -2, ], aes(label = gene))+
  theme_classic(base_size = 30)+
  xlab("Log FC stimulation") + ylab("Log FC restimulation")
ggsave(paste(fig.file.path, "hypersuppressed_boxes_borders.pdf", sep = "/"), width = 8, height = 8, useDingbats = FALSE )

#plot secondary induced genes
plot_bb_down_up <-ggplot(bb_down_up_2exps, aes(x = stim_log2FoldChange, y = restim_log2FoldChange)) +
  geom_point(data = bb_down_up_2exps, color = "orange")+
  geom_text_repel(data = bb_down_up_2exps[bb_down_up_2exps$stim_log2FoldChange < -2 & bb_down_up_2exps$restim_log2FoldChange > 2, ], aes(label = gene))+
  theme_classic(base_size = 30)+
  xlab("Log FC stimulation") + ylab("Log FC restimulation")
ggsave(paste(fig.file.path, "bb_down_up.pdf", sep = "/"), width = 10, height = 10, useDingbats = FALSE )

sec_ind_genes_to_label <- c("Hoxa1", "Ccr10", "Glis3", "Eno3", "Etv1", "Trim37",
                            "Irf2bp2", "Pagr1a", "Gm28539", "Gpr35")

ggplot(bb_down_up_2exps, aes(x = stim_log2FoldChange, y = restim_log2FoldChange)) +
  geom_point(data = bb_down_up_2exps, fill = "orange", color = "black", shape =21)+
  geom_label_repel(data = bb_down_up_2exps %>%
                    dplyr::filter(gene %in% sec_ind_genes_to_label), 
                  aes(label = gene), size = 6, min.segment.length = 0.2)+
  theme_classic(base_size = 30)+
  xlab("Log FC stimulation") + ylab("Log FC restimulation")
ggsave(paste(fig.file.path, "secondary_induced_specific_genes_labeled_boxes_borders.pdf", sep = "/"), width = 8, height = 8, useDingbats = FALSE )



#GSEA
#import deseq files for stim and restim
bb_stim_data <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/2_exps_bb_stim_resdata.csv", sep = "/", stringsAsFactors = FALSE, row.names=1)
bb_restim_data <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/2_exps_bb_restim_resdata.csv", sep = "/", stringsAsFactors = FALSE, row.names=1)

mousetohuman <- read.csv(paste(data.file.path, "human_to_mouse_genes.csv", sep = "/"), stringsAsFactors = FALSE)
bb_stim_data$MGI.symbol <- rownames(bb_stim_data)

#restim
mousetohuman <- read.csv(paste(data.file.path, "human_to_mouse_genes.csv", sep = "/"), stringsAsFactors = FALSE)
bb_restim_data$MGI.symbol <- rownames(bb_restim_data)

#all gene sets file msigdb
gene_sets <- gmtPathways(paste(data.file.path, "all_genesets_msigdb_symbols.gmt", sep = "/" ))

#gsea stim
bb_stim_gsea <- bb_stim_data %>%
  dplyr::select(MGI.symbol, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(MGI.symbol) %>%
  summarize(stat=mean(stat)) %>%
  left_join(mousetohuman, by = "MGI.symbol") %>%
  dplyr::select(HGNC.symbol, stat) %>%
  na.omit()
ranks_stim <- deframe(bb_stim_gsea)

tolerance_gene_sets <- gene_sets[grep("HALLMARK|BIOCARTA|REACTOME|BP|CC|MF", names(gene_sets))]

fgseaRes_stim <- fgsea(pathways= tolerance_gene_sets, stats=ranks_stim, nperm=10000)

fgseaResTidy_stim <- fgseaRes_stim %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy_stim$leadingEdge <- sapply(fgseaResTidy_stim$leadingEdge, function(x) paste(x, collapse = " "))
write.csv(fgseaResTidy_stim ,"C:/Users/julie/Box/tolerance_bulk/data/stim_gsea_table.csv", row.names=FALSE)

#plot gsea results stim
ggplot(fgseaResTidy_stim[fgseaResTidy_stim$padj < 0.05, ], aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA Pathways Stim") + 
  theme_minimal()

#gsea restim
bb_restim_gsea <- bb_restim_data %>%
  dplyr::select(MGI.symbol, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(MGI.symbol) %>%
  summarize(stat=mean(stat)) %>%
  left_join(mousetohuman, by = "MGI.symbol") %>%
  dplyr::select(HGNC.symbol, stat) %>%
  na.omit()
ranks_restim <- deframe(bb_restim_gsea)

fgseaRes_restim <- fgsea(pathways= tolerance_gene_sets, stats=ranks_restim, nperm=10000)

fgseaResTidy_restim <- fgseaRes_restim %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy_restim$leadingEdge <- sapply(fgseaResTidy_restim$leadingEdge, function(x) paste(x, collapse = " "))
write.csv(fgseaResTidy_restim ,"C:/Users/julie/Box/tolerance_bulk/data/restim_gsea_table.csv", row.names=FALSE)

#plot gsea results restim
ggplot(fgseaResTidy_restim[fgseaResTidy_restim$padj < 0.05, ], aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA Pathways Restim") + 
  theme_minimal()

#import stim gsea and restim gsea csv files- NEED TO COMBINE CODE FOR THIS INTO THIS FILE
#NOTE- I ADDED _STIM AND _RESTIM AFTER HEADINGS FOR THESE FILES IN EXCEL, EXCEPT FOR "PATHWAY"
gsea_stim_table_labeled <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/stim_gsea_res_labelled_2exps.csv", sep = ",", stringsAsFactors = FALSE)
gsea_restim_table_labeled <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/restim_gsea_res_labelled_2exps.csv", sep = ",", stringsAsFactors = FALSE)

#combine GSEA files into 1 file
gsea_merge_all <- full_join(gsea_stim_table_labeled, gsea_restim_table_labeled, by = "pathway", all = TRUE)
write.csv(gsea_merge_all, "C:/Users/julie/Box/tolerance_bulk/data/gsea_merge_all.csv")

#separate out by category
#add rank column which adds up NES_stim and NES_restim
#move rank column after pathway column
gsea_tolerizable_2exps <- filter(gsea_merge_all, NES_stim > 0, NES_restim < 0, padj_stim <=0.1, padj_restim <=0.1)
gsea_tolerizable_2exps$rank <- abs(gsea_tolerizable_2exps$NES_stim) + abs(gsea_tolerizable_2exps$NES_restim)
gsea_tolerizable_2exps <- gsea_tolerizable_2exps %>% relocate(rank, .after = pathway)
write.csv(gsea_tolerizable_2exps, "C:/Users/julie/Box/tolerance_bulk/data/gsea_tolerizable_2exps.csv")

gsea_hyperinduced_2exps <- filter(gsea_merge_all, NES_stim > 0, NES_restim > 0, padj_stim <=0.1, padj_restim <=0.1)
gsea_hyperinduced_2exps$rank <- abs(gsea_hyperinduced_2exps$NES_stim) + abs(gsea_hyperinduced_2exps$NES_restim)
gsea_hyperinduced_2exps <- gsea_hyperinduced_2exps %>% relocate(rank, .after = pathway)
write.csv(gsea_hyperinduced_2exps, "C:/Users/julie/Box/tolerance_bulk/data/gsea_hyperinduced_2exps.csv")

gsea_down_down_2exps <- filter(gsea_merge_all, NES_stim < 0, NES_restim < 0 , padj_stim <=0.1, padj_restim <=0.1)
gsea_down_down_2exps$rank <- abs(gsea_down_down_2exps$NES_stim) + abs(gsea_down_down_2exps$NES_restim)
gsea_down_down_2exps <- gsea_down_down_2exps %>% relocate(rank, .after = pathway)
write.csv(gsea_down_down_2exps, "C:/Users/julie/Box/tolerance_bulk/data/gsea_down_down_2exps.csv")

gsea_down_up_2exps <-filter(gsea_merge_all, NES_stim < 0, NES_restim > 0 , padj_stim <=0.1, padj_restim <=0.1)
gsea_down_up_2exps$rank <- abs(gsea_down_up_2exps$NES_stim) + abs(gsea_down_up_2exps$NES_restim)
gsea_down_up_2exps <- gsea_down_up_2exps %>% relocate(rank, .after = pathway)
write.csv(gsea_down_up_2exps, "C:/Users/julie/Box/tolerance_bulk/data/gsea_down_up_2exps.csv")

#highlight pathways in gsea merge all that I want to put in plot
highlight <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "GOBP_LIPOPOLYSACCHARIDE_MEDIATED_SIGNALING_PATHWAY", "GOBP_RESPONSE_TO_INTERFERON_GAMMA", "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES",
               "GOBP_ACUTE_PHASE_RESPONSE", "GOBP_DETOXIFICATION", "GOMF_ANTIOXIDANT_ACTIVITY", "BIOCARTA_PROTEASOME_PATHWAY",
               "GNF2_CCNA2", "GNF2_CCNB2", "MANNO_MIDBRAIN_NEUROTYPES_HPROGBP", "MEBARKI_HCC_PROGENITOR_FZD8CRD_UP",
               "MOLENAAR_TARGETS_OF_CCND1_AND_CDK4_DN", "REACTOME_METALLOTHIONEINS_BIND_METALS", "GOBP_FATTY_ACID_CATABOLIC_PROCESS",
               "GOBP_MONOCARBOXYLIC_ACID_CATABOLIC_PROCESS", "GOBP_FATTY_ACID_BETA_OXIDATION", "GOBP_INTRACILIARY_TRANSPORT")
gsea_merge_all$label <- ifelse(gsea_merge_all$pathway %in% highlight, TRUE, FALSE)

#PLOT THE GSEA DATA
#plot GSEA data with specific pathways from Tanja
#Tolerized
tolerized_gsea_pathways<- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "GOBP_DEFENSE_RESPONSE_TO_SYMBIONT", "HALLMARK_INFLAMMATORY_RESPONSE", "GOBP_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY", "GOBP_CELLULAR_RESPONSE_TO_INTERLEUKIN_1", "GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE", "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY", "GOBP_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "GOBP_CYTOPLASMIC_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY")
tolerized_data <- read.csv("C:/Users/julie/Box/tolerance_bulk/data/gsea_tolerizable_2exps.csv", row.names = 1)
stim_data <- tolerized_data %>%
  dplyr::select(pathway, rank, pval_stim, NES_stim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "Acute") %>%
  filter(pathway %in% tolerized_gsea_pathways)

restim_data <- tolerized_data %>%
  dplyr::select(pathway, rank, pval_restim, NES_restim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "Memory") %>%
  filter(pathway %in% tolerized_gsea_pathways)

plot_df <- rbind(stim_data, restim_data)
plot_df$treatment <- factor(plot_df$treatment, levels = c("Acute", "Memory"))

ggplot(plot_df, aes(x = NES, y = pathway, fill = -log10(pval))) +
  geom_point(pch = 21, size= 7) +
  scale_fill_distiller(palette="YlOrRd", trans= "reverse")+
  facet_wrap(~treatment) +
  scale_x_continuous(limits = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic(base_size = 20)+
  theme(panel.background = element_rect(fill = NA, color = "black"))
ggsave("C:/Users/julie/Box/tolerance_bulk/figures/gsea_tolerizable_2exps_specificpathways.pdf", useDingbats = FALSE, width = 15, height = 5)


#TOLERIZED GSEA plot
# Order by rank, take top 10, make dot plot
n_pathways <- 10
tolerized_data <- read.csv("C:/Users/julie/Box/tolerance_bulk/data/gsea_tolerizable_2exps.csv", row.names = 1)

stim_data <- tolerized_data %>%
  dplyr::select(pathway, rank, pval_stim, NES_stim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "stimulation") %>%
  head(n_pathways)

restim_data <- tolerized_data %>%
  dplyr::select(pathway, rank, pval_restim, NES_restim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "restimulation") %>%
  head(n_pathways)

plot_df <- rbind(stim_data, restim_data)
plot_df$treatment <- factor(plot_df$treatment, levels = c("stimulation", "restimulation"))

ggplot(plot_df, aes(x = NES, y = pathway, fill = -log10(pval), size = -log10(pval))) +
  geom_point(pch = 21) +
  scale_fill_distiller(palette="YlOrRd", trans= "reverse")+
  facet_wrap(~treatment) +
  scale_x_continuous(limits = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic(base_size = 20)+
  theme(panel.background = element_rect(fill = NA, color = "black"))
ggsave("C:/Users/julie/Box/tolerance_bulk/figures/gsea_tolerizable_2exps_dotplot.pdf", useDingbats = FALSE, width = 15, height = 5)

#HYPERINDUCED GSEA plot
#Hyperinduced plot with specific pathways
hyperinduced_pathways <- c("REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES", "GOBP_RESPONSE_TO_TOXIC_SUBSTANCE","GOBP_TRIPEPTIDE_TRANSPORT", "GOBP_SUPEROXIDE_METABOLIC_PROCESS", "HALLMARK_MTORC1_SIGNALING", "GOBP_ORGANIC_ACID_TRANSMEMBRANE_TRANSPORT", "GOMF_PEROXIREDOXIN_ACTIVITY", "GOBP_HEPATIC_STELLATE_CELL_ACTIVATION", "GOBP_CARBOXYLIC_ACID_TRANSPORT", "GOBP_RESPONSE_TO_OLEIC_ACID")
hyperinduced_data <- read.csv("C:/Users/julie/Box/tolerance_bulk/data/gsea_hyperinduced_2exps.csv", row.names = 1)

stim_data <- hyperinduced_data %>%
  dplyr::select(pathway, rank, pval_stim, NES_stim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "Acute") %>%
  filter(pathway %in% hyperinduced_pathways)

restim_data <- hyperinduced_data %>%
  dplyr::select(pathway, rank, pval_restim, NES_restim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "Memory") %>%
  filter(pathway %in% hyperinduced_pathways)

plot_df <- rbind(stim_data, restim_data)
plot_df$treatment <- factor(plot_df$treatment, levels = c("Acute", "Memory"))

ggplot(plot_df, aes(x = NES, y = pathway, fill = -log10(pval))) +
  geom_point(pch = 21, size = 7) +
  scale_fill_distiller(palette="YlOrRd", trans= "reverse")+
  facet_wrap(~treatment) +
  scale_x_continuous(limits = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic(base_size = 20)+
  theme(panel.background = element_rect(fill = NA, color = "black"))
ggsave("C:/Users/julie/Box/tolerance_bulk/figures/gsea_hyperinduced_2exps_specificpathways.pdf", useDingbats = FALSE, width = 15, height = 5)


# Order by rank, take top 10, make dot plot
n_pathways <- 10
hyperinduced_data <- read.csv("C:/Users/julie/Box/tolerance_bulk/data/gsea_hyperinduced_2exps.csv", row.names = 1)

stim_data <- hyperinduced_data %>%
  dplyr::select(pathway, rank, pval_stim, NES_stim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "stimulation") %>%
  head(n_pathways)

restim_data <- hyperinduced_data %>%
  dplyr::select(pathway, rank, pval_restim, NES_restim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "restimulation") %>%
  head(n_pathways)

plot_df <- rbind(stim_data, restim_data)
plot_df$treatment <- factor(plot_df$treatment, levels = c("stimulation", "restimulation"))

ggplot(plot_df, aes(x = NES, y = pathway, fill = -log10(pval), size = -log10(pval))) +
  geom_point(pch = 21, size = 7) +
  scale_fill_distiller(palette="YlOrRd", trans= "reverse")+
  facet_wrap(~treatment) +
  scale_x_continuous(limits = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic(base_size = 20)+
  theme(panel.background = element_rect(fill = NA, color = "black"))
ggsave("C:/Users/julie/Box/tolerance_bulk/figures/gsea_hyperinduced_2exps_dotplot.pdf", useDingbats = FALSE, width = 15, height = 5)

#SECONDARY INDUCED GSEA PLOT (DOWN UP)
#HYPERINDUCED GSEA plot
# Order by rank, take top 10, make dot plot
n_pathways <- 10
secondary_induced_data <- read.csv("C:/Users/julie/Box/tolerance_bulk/data/gsea_down_up_2exps.csv", row.names = 1)

stim_data <- secondary_induced_data %>%
  dplyr::select(pathway, rank, pval_stim, NES_stim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "Acute") %>%
  head(n_pathways)

restim_data <- secondary_induced_data %>%
  dplyr::select(pathway, rank, pval_restim, NES_restim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "Memory") %>%
  head(n_pathways)

plot_df <- rbind(stim_data, restim_data)
plot_df$treatment <- factor(plot_df$treatment, levels = c("Acute", "Memory"))

ggplot(plot_df, aes(x = NES, y = pathway, fill = -log10(pval), size = -log10(pval))) +
  geom_point(pch = 21, size = 7) +
  scale_fill_distiller(palette="YlOrRd", trans= "reverse")+
  facet_wrap(~treatment) +
  scale_x_continuous(limits = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic(base_size = 20)+
  theme(panel.background = element_rect(fill = NA, color = "black"))
ggsave("C:/Users/julie/Box/tolerance_bulk/figures/gsea_secondary_induced_2exps_dotplot.pdf", useDingbats = FALSE, width = 15, height = 5)

#hypersuppressed GSEA plot
n_pathways <- 10
hyper_suppressed_data <- read.csv("C:/Users/julie/Box/tolerance_bulk/data/gsea_down_down_2exps.csv", row.names = 1)

stim_data <- hyper_suppressed_data %>%
  dplyr::select(pathway, rank, pval_stim, NES_stim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "Acute") %>%
  head(n_pathways)

restim_data <- hyper_suppressed_data %>%
  dplyr::select(pathway, rank, pval_restim, NES_restim) %>%
  arrange(desc(rank)) %>%
  `colnames<-`(c("pathway", "rank", "pval", "NES")) %>%
  mutate(treatment = "Memory") %>%
  head(n_pathways)

plot_df <- rbind(stim_data, restim_data)
plot_df$treatment <- factor(plot_df$treatment, levels = c("Acute", "Memory"))

ggplot(plot_df, aes(x = NES, y = pathway, fill = -log10(pval))) +
  geom_point(pch = 21, size= 7) +
  scale_fill_distiller(palette="YlOrRd", trans= "reverse")+
  facet_wrap(~treatment) +
  scale_x_continuous(limits = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic(base_size = 20)+
  theme(panel.background = element_rect(fill = NA, color = "black"))
ggsave("C:/Users/julie/Box/tolerance_bulk/figures/gsea_hyper_suppressed_2exps_dotplot.pdf", useDingbats = FALSE, width = 20, height = 5)


#Merge IPA stim and restim files into 1 document
#import files
ipa_stim_2exps <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/IPA_stim_upstream_regulators.csv", sep = ",", stringsAsFactors = FALSE)
ipa_restim_2exps <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/IPA_restim_upstream_analysis.csv", sep = ",", stringsAsFactors = FALSE)


#merge stim and restim and keep all rows but combine those that have overlapping genes
ipa_upstream_merge <- ipa_merge_all_2exps <- full_join(ipa_stim_2exps, ipa_restim_2exps, by = "Upstream_Regulator", all = TRUE)
write.csv(ipa_upstream_merge, "C:/Users/julie/Box/tolerance_bulk/data/ipa_upstream_regulators_merged.csv")


#TABLES FOR SUPPLEMENTARY FIGURES
#deseq results
bb_stim_final_supp <- bb_stim_res_data
bb_stim_final_supp$contrast <- "Acute"

bb_restim_final_supp <- bb_restim_res_data
bb_restim_final_supp$contrast <- "Memory"

bb_stim_restim_supp <- rbind(bb_stim_final_supp, bb_restim_final_supp)
write.csv(bb_stim_restim_supp, "C:/Users/julie/Box/tolerance_bulk/data/combined_deseqresults_table.csv")

#fgsea results
bb_gsea_stim_final_supp <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/stim_gsea_res_2exps.csv", sep = ",", stringsAsFactors = FALSE)
bb_gsea_stim_final_supp$contrast <- "Acute"
bb_gsea_restim_final_supp <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/data/restim_gsea_res_2exps.csv", sep = ",", stringsAsFactors = FALSE)
bb_gsea_restim_final_supp$contrast <- "Memory"

bb_gsea_final_supp <- rbind(bb_gsea_stim_final_supp, bb_gsea_restim_final_supp)
write.csv(bb_gsea_final_supp, "C:/Users/julie/Box/tolerance_bulk/data/combined_gsea_results_table.csv")


#import genes.results files and get TPM values
#import genes.results files

samples <- c("24_no_stim", "25_bb_restim", "27_bb_6hr", "29_no_stim", "30_bb_restim", "32_bb_6hr")
all_TPM <- lapply(samples, function(x) {
  df <- read.csv(glue("{data.file.path}/{x}.genes.results"), sep = "\t")
  df <- df[, c("gene_id", "TPM")]
  colnames(df) <- c("ensembl_gene_id", glue("TPM_{x}"))  
  return(df)
})

tpm_matrix <- reduce(all_TPM, left_join, by= "ensembl_gene_id")

write.csv(tpm_matrix, "C:/Users/julie/Box/tolerance_bulk/data/tpm_values.csv")

#PCAPLOT

# Get rid of all genes that have no expression
expressed_genes <- apply(tpm_matrix, 1, function(x) if (max(x) > 0) TRUE else FALSE)

tpm_matrix <- tpm_matrix[expressed_genes,]
tpm_matrix$ensembl_gene_id <- NULL

# Remove TCR genes
for (regex in c("^TRB", "^TRAV", "^TRAC")) {
  tpm_matrix <- tpm_matrix[!grepl(regex, rownames(tpm_matrix)),]
  
}

metadata <- read.csv(file = "C:/Users/julie/Box/tolerance_bulk/tolerance_pca_metadata_2exps.csv", header= TRUE, sep= ",")
rownames(metadata)<-metadata$X
metadata$X <- NULL


if(all(colnames(tpm_matrix) == rownames(metadata))) {
  var_gene_pca <- pca(tpm_matrix, metadata = metadata, scale = TRUE, removeVar = 0.9)
}

#plot of all samples colored by stim
biplot(var_gene_pca, colby = "treatment", legendPosition = "right", 
       colkey = c("none" = 'purple', "stimulation" = 'orange', "re-stimulation" = 'red'), 
       gridlines.minor = FALSE,
       pointSize = 4,
       lab = NULL,
       axisLabSize = 15)
ggsave("C:/Users/julie/Box/tolerance_bulk/figures/pca_2exps.pdf", useDingbats = FALSE, width = 8, height = 5)



pairsplot(var_gene_pca, colby= "stim",
          colkey = c("none" = 'purple', "bb" = 'orange', "pam" = "red"), 
          lab = rownames(metadata),
          gridlines.minor = FALSE,
          axisLabSize = 15)