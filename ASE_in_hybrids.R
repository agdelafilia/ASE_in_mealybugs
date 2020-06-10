rm(list=ls())
ls()
library(dplyr)
library(tidyr)
library(Rmisc)
library(ggplot2)
library(readr)
library(lattice)
library(grid)
library(gridExtra)

##################################

##### Import and organise raw data

# import raw output from ASEReadCounter

s1.hybrid <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/S1.csv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
s3.hybrid <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/S3.csv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
s6.hybrid <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/S6.csv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
t1.hybrid <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/T1.csv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
t3.hybrid <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/T3.csv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
t6.hybrid <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/T6.csv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)

# import expression data from RSEM

S1_A_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S1_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S1_A_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S1_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S1_B_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S1_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S1_B_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S1_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S3_A_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S3_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S3_A_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S3_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S3_B_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S3_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S3_B_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S3_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S6_A_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S6_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S6_A_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S6_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S6_B_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S6_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S6_B_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S6_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T1_A_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T1_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T1_A_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T1_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T1_B_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T1_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T1_B_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T1_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T3_A_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T3_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T3_A_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T3_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T3_B_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T3_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T3_B_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T3_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T6_A_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T6_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T6_A_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T6_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T6_B_6_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T6_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T6_B_7_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T6_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)

# merge datasets and obtain a list of expressed genes in soma and testis (TPM >= 1)

s_tpm = S1_A_6_genes_results[c(1,6)]
colnames(s_tpm)[colnames(s_tpm)=="TPM"] <- "S1_A_6"
s_tpm$S1_A_7 = S1_A_7_genes_results$TPM
s_tpm$S1_B_6 = S1_B_6_genes_results$TPM
s_tpm$S1_B_7 = S1_B_7_genes_results$TPM
s_tpm$S3_A_6 = S3_A_6_genes_results$TPM
s_tpm$S3_A_7 = S3_A_7_genes_results$TPM
s_tpm$S3_B_6 = S3_B_6_genes_results$TPM
s_tpm$S3_B_7 = S3_B_7_genes_results$TPM
s_tpm$S6_A_6 = S6_A_6_genes_results$TPM
s_tpm$S6_A_7 = S6_A_7_genes_results$TPM
s_tpm$S6_B_6 = S6_B_6_genes_results$TPM
s_tpm$S6_B_7 = S6_B_7_genes_results$TPM

t_tpm = T1_A_6_genes_results[c(1,6)]
colnames(t_tpm)[colnames(t_tpm)=="TPM"] <- "T1_A_6"
t_tpm$T1_A_7 = T1_A_7_genes_results$TPM
t_tpm$T1_B_6 = T1_B_6_genes_results$TPM
t_tpm$T1_B_7 = T1_B_7_genes_results$TPM
t_tpm$T3_A_6 = T3_A_6_genes_results$TPM
t_tpm$T3_A_7 = T3_A_7_genes_results$TPM
t_tpm$T3_B_6 = T3_B_6_genes_results$TPM
t_tpm$T3_B_7 = T3_B_7_genes_results$TPM
t_tpm$T6_A_6 = T6_A_6_genes_results$TPM
t_tpm$T6_A_7 = T6_A_7_genes_results$TPM
t_tpm$T6_B_6 = T6_B_6_genes_results$TPM
t_tpm$T6_B_7 = T6_B_7_genes_results$TPM

s_tpm$s_TPM <- rowMeans(s_tpm[,-1]) # average expression across all samples
t_tpm$t_TPM <- rowMeans(t_tpm[,-1]) # average expression across all samples
colnames(s_tpm)[colnames(s_tpm)=="gene_id"] <- "gene"
colnames(t_tpm)[colnames(t_tpm)=="gene_id"] <- "gene"

s_tpm_expressed <- s_tpm[s_tpm$s_TPM >= 1,]
t_tpm_expressed <- t_tpm[t_tpm$t_TPM >= 1,]
s_tpm_not_expressed <- s_tpm[s_tpm$s_TPM < 1,]
t_tpm_not_expressed <- t_tpm[t_tpm$t_TPM < 1,]

s_tpm_expressed_join=s_tpm_expressed[c(1,14)] # prepare dfs to intersect with ASE data
t_tpm_expressed_join=t_tpm_expressed[c(1,14)] # prepare dfs to intersect with ASE data

s_tpm_not_expressed_join=s_tpm_not_expressed[c(1,14)]
t_tpm_not_expressed_join=t_tpm_not_expressed[c(1,14)]

# examine ASE counts

nrow(s1.hybrid)
nrow(s3.hybrid)
nrow(s6.hybrid)
nrow(t1.hybrid)
nrow(t3.hybrid)
nrow(t6.hybrid)

# add a variant column

s1.hybrid$variant=paste(s1.hybrid$contig,s1.hybrid$position,sep=":")
s3.hybrid$variant=paste(s3.hybrid$contig,s3.hybrid$position,sep=":")
s6.hybrid$variant=paste(s6.hybrid$contig,s6.hybrid$position,sep=":")
t1.hybrid$variant=paste(t1.hybrid$contig,t1.hybrid$position,sep=":")
t3.hybrid$variant=paste(t3.hybrid$contig,t3.hybrid$position,sep=":")
t6.hybrid$variant=paste(t6.hybrid$contig,t6.hybrid$position,sep=":")

# prepare dfs to merge by tissue

s1.hybrid_to_merge = select(s1.hybrid,variant,refCount,altCount,totalCount)
colnames(s1.hybrid_to_merge) = c("variant","refCount.s1","altCount.s1","totalCount.s1")
s3.hybrid_to_merge = select(s3.hybrid,variant,refCount,altCount,totalCount)
colnames(s3.hybrid_to_merge) = c("variant","refCount.s3","altCount.s3","totalCount.s3")
s6.hybrid_to_merge = select(s6.hybrid,variant,refCount,altCount,totalCount)
colnames(s6.hybrid_to_merge) = c("variant","refCount.s6","altCount.s6","totalCount.s6")

t1.hybrid_to_merge = select(t1.hybrid,variant,refCount,altCount,totalCount)
colnames(t1.hybrid_to_merge) = c("variant","refCount.t1","altCount.t1","totalCount.t1")
t3.hybrid_to_merge = select(t3.hybrid,variant,refCount,altCount,totalCount)
colnames(t3.hybrid_to_merge) = c("variant","refCount.t3","altCount.t3","totalCount.t3")
t6.hybrid_to_merge = select(t6.hybrid,variant,refCount,altCount,totalCount)
colnames(t6.hybrid_to_merge) = c("variant","refCount.t6","altCount.t6","totalCount.t6")

# merge all and obtain a full list of variants present in at least 1 sample

all_merged_hybrid.0 = full_join(s1.hybrid_to_merge,s3.hybrid_to_merge,by="variant")
all_merged_hybrid.00 = full_join(all_merged_hybrid.0,s6.hybrid_to_merge,by="variant")
all_merged_hybrid.000 = full_join(all_merged_hybrid.00,t1.hybrid_to_merge,by="variant")
all_merged_hybrid.0000 = full_join(all_merged_hybrid.000,t3.hybrid_to_merge,by="variant")
all_merged_hybrid_raw = full_join(all_merged_hybrid.0000,t6.hybrid_to_merge,by="variant")
all_merged_hybrid_for_status = all_merged_hybrid_raw[c(1)]
all_merged_hybrid_for_status$contig <- sub(":.*", '', all_merged_hybrid_for_status$variant)
nrow(all_merged_hybrid_for_status)

##### Assign genome annotation features to variants

# start by assigning orphan variants (SNPs on contigs without annotated features)

contigs_with_anno <- read_csv("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/contigs_with_anno.bed", 
                              col_names = FALSE) # list of contigs with annotation features 
colnames(contigs_with_anno) <- c("contig")
contigs_with_anno_list <- list(unique(contigs_with_anno)) # a list of all contigs with annotation features

all_merged_hybrid_for_status_orphan_0 = all_merged_hybrid_for_status
all_merged_hybrid_for_status_orphan_0$status <- ifelse(all_merged_hybrid_for_status_orphan_0$contig %in% contigs_with_anno$contig == TRUE,NA,"orphan") # declare a variant orphan if its contig is not in the list
all_merged_hybrid_for_status_orphan = all_merged_hybrid_for_status_orphan_0[complete.cases(all_merged_hybrid_for_status_orphan_0),] # keep orphan variants only
all_merged_hybrid_for_status_orphan$CDS <- NA # assign NA for CDS/mRNA
all_merged_hybrid_for_status_orphan$mRNA <- NA
nrow(all_merged_hybrid_for_status_orphan)
nrow(anti_join(all_merged_hybrid_for_status_orphan,all_merged_hybrid_for_status,by="variant")) # sanity check, should be 0

# to assign rest of categories: retrieve counts intersected with annotation (obtained with BEDTools)

s1.hybrid_ASE_anno_R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/S1_anno_for_R.bed", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
s3.hybrid_ASE_anno_R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/S3_anno_for_R.bed", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
s6.hybrid_ASE_anno_R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/S6_anno_for_R.bed", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
t1.hybrid_ASE_anno_R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/T1_anno_for_R.bed", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
t3.hybrid_ASE_anno_R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/T3_anno_for_R.bed", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
t6.hybrid_ASE_anno_R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/T6_anno_for_R.bed", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)

anno_all_hybrid = rbind(s1.hybrid_ASE_anno_R,s3.hybrid_ASE_anno_R,s6.hybrid_ASE_anno_R,t1.hybrid_ASE_anno_R,t3.hybrid_ASE_anno_R,t6.hybrid_ASE_anno_R) # merge all
nrow(anno_all_hybrid) - nrow(s1.hybrid_ASE_anno_R) - nrow(s3.hybrid_ASE_anno_R) - nrow(s6.hybrid_ASE_anno_R) - nrow(t1.hybrid_ASE_anno_R) - nrow(t3.hybrid_ASE_anno_R) - nrow(t6.hybrid_ASE_anno_R) # sanity check, should be 0
anno_all_unique_hybrid <- distinct(anno_all_hybrid) # remove duplicates
colnames(anno_all_unique_hybrid) = c("contig","position","anno","type","info")
anno_all_unique_hybrid$variant=paste(anno_all_unique_hybrid$contig,anno_all_unique_hybrid$position,sep=":")
ddply(anno_all_unique_hybrid,c("type"),summarise,N = length(variant)) # inspect

# start with exonic sites, so that exonic SNPs that map to other features (e.g. introns) will be given priority

exonic.sites.hybrid <- anno_all_unique_hybrid[anno_all_unique_hybrid$type == "CDS",]
exonic.sites.hybrid.merge_0 <- exonic.sites.hybrid[c(6,5)]
exonic.sites.hybrid.merge = separate(exonic.sites.hybrid.merge_0,info,into = c("CDS","mRNA"),sep = ";") # obtain information from the annotation field
exonic.sites.hybrid.merge$CDS <- sub('ID=', '', exonic.sites.hybrid.merge$CDS) # assign CDS
exonic.sites.hybrid.merge$mRNA <- sub('Parent=', '', exonic.sites.hybrid.merge$mRNA) # assign mRNA

all_merged_hybrid_for_status_exon=ddply(exonic.sites.hybrid.merge,c("variant"), summarize,
                                      CDS=paste(CDS,collapse=","), 
                                      mRNA=paste(mRNA,collapse=",")) # collapse annotations for variants with >1 feature

all_merged_hybrid_for_status_exon$contig <- sub(":.*", '', all_merged_hybrid_for_status_exon$variant)
all_merged_hybrid_for_status_exon$status <- "exonic"
nrow(all_merged_hybrid_for_status_exon)

# now assign all intronic sites

intronic.sites.hybrid <- anno_all_unique_hybrid[anno_all_unique_hybrid$type == "intron",]
intronic.sites.hybrid.merge_0 <- intronic.sites.hybrid[c(6,5)]
intronic.sites.hybrid.merge_0$mRNA <- sub('Parent=', '', intronic.sites.hybrid.merge_0$info)

intronic.sites.hybrid.merge.collapsed = ddply(intronic.sites.hybrid.merge_0,c("variant"), summarize,
                                            mRNA=paste(mRNA,collapse=",")) # collapse annotations for variants with >1 feature

all_merged_hybrid_for_status_intron=anti_join(intronic.sites.hybrid.merge.collapsed,all_merged_hybrid_for_status_exon,by="variant") # exclude intronic variants if they already map to exons
all_merged_hybrid_for_status_intron$contig <- sub(":.*", '', all_merged_hybrid_for_status_intron$variant)
all_merged_hybrid_for_status_intron$status <- "intronic"
all_merged_hybrid_for_status_intron$CDS <- NA
nrow(all_merged_hybrid_for_status_intron)

# assign intergenic sites

all_merged_hybrid_for_status_non_intergenic = rbind(all_merged_hybrid_for_status_exon,all_merged_hybrid_for_status_intron, all_merged_hybrid_for_status_orphan) # merge all variants with exonic/intronic/orphan status
all_merged_hybrid_for_status_intergenic_0 = anti_join(all_merged_hybrid_for_status,all_merged_hybrid_for_status_non_intergenic,by="variant") # assign intergenic status to all remaining variants
all_merged_hybrid_for_status_intergenic = all_merged_hybrid_for_status_intergenic_0[c(1)]
all_merged_hybrid_for_status_intergenic$status <- "intergenic"
all_merged_hybrid_for_status_intergenic$CDS <- NA
all_merged_hybrid_for_status_intergenic$mRNA <- NA
all_merged_hybrid_for_status_intergenic$contig <- sub(":.*", '', all_merged_hybrid_for_status_intergenic$variant)
nrow(all_merged_hybrid_for_status_intergenic) + nrow(all_merged_hybrid_for_status_non_intergenic)

# merge all

all_merged_hybrid_anno_0 = rbind(all_merged_hybrid_for_status_non_intergenic,all_merged_hybrid_for_status_intergenic)
all_merged_hybrid_anno = full_join(all_merged_hybrid_for_status,all_merged_hybrid_anno_0)
annotation_counts_hybrid = ddply(all_merged_hybrid_anno,c("status"),summarise,N = length(variant))

# intersect annotation and read counts

all_merged_hybrid_complete_not_ordered = full_join(all_merged_hybrid_anno,all_merged_hybrid_raw,by="variant","contig") # combine annotation with ASE counts
all_merged_hybrid_complete_not_ordered$order1 = as.integer(gsub(".*_", '', all_merged_hybrid_complete_not_ordered$contig)) # order by contig and position
all_merged_hybrid_complete_not_ordered$order2 = as.integer(gsub(".*:", '', all_merged_hybrid_complete_not_ordered$variant))
all_merged_hybrid_complete = all_merged_hybrid_complete_not_ordered[order(all_merged_hybrid_complete_not_ordered$order1, all_merged_hybrid_complete_not_ordered$order2),]
nrow(all_merged_hybrid_complete)

##### Filter variants

# initial filter: keep variants shared between all tissue replicates

all_merged_hybrid_complete_soma=all_merged_hybrid_complete[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)]
all_merged_hybrid_complete_testis=all_merged_hybrid_complete[c(1,2,3,4,5,15,16,17,18,19,20,21,22,23)]

shared_raw_hybrid_soma = all_merged_hybrid_complete_soma[!is.na(all_merged_hybrid_complete_soma$totalCount.s1) 
                                                     & !is.na(all_merged_hybrid_complete_soma$totalCount.s3) 
                                                     & !is.na(all_merged_hybrid_complete_soma$totalCount.s6),]
shared_raw_hybrid_testis = all_merged_hybrid_complete_testis[!is.na(all_merged_hybrid_complete_testis$totalCount.t1) 
                                                         & !is.na(all_merged_hybrid_complete_testis$totalCount.t3) 
                                                         & !is.na(all_merged_hybrid_complete_testis$totalCount.t6),]

nrow(shared_raw_hybrid_soma)
nrow(shared_raw_hybrid_testis)
shared_raw_hybrid_soma_counts = ddply(shared_raw_hybrid_soma,c("status"),summarise,N = length(variant))
shared_raw_hybrid_testis_counts = ddply(shared_raw_hybrid_testis,c("status"),summarise,N = length(variant))

# first filter: remove SNPs in which valid read depth >90% of total read depth (total count + counts from multiple mapping reads + non ref, non alt bases)

s1.hybrid$prop_other=(s1.hybrid$otherBases+s1.hybrid$lowMAPQDepth)/(s1.hybrid$totalCount+s1.hybrid$otherBases+s1.hybrid$lowMAPQDepth)
s3.hybrid$prop_other=(s3.hybrid$otherBases+s3.hybrid$lowMAPQDepth)/(s3.hybrid$totalCount+s3.hybrid$otherBases+s3.hybrid$lowMAPQDepth)
s6.hybrid$prop_other=(s6.hybrid$otherBases+s6.hybrid$lowMAPQDepth)/(s6.hybrid$totalCount+s6.hybrid$otherBases+s6.hybrid$lowMAPQDepth)
t1.hybrid$prop_other=(t1.hybrid$otherBases+t1.hybrid$lowMAPQDepth)/(t1.hybrid$totalCount+t1.hybrid$otherBases+t1.hybrid$lowMAPQDepth)
t3.hybrid$prop_other=(t3.hybrid$otherBases+t3.hybrid$lowMAPQDepth)/(t3.hybrid$totalCount+t3.hybrid$otherBases+t3.hybrid$lowMAPQDepth)
t6.hybrid$prop_other=(t6.hybrid$otherBases+t6.hybrid$lowMAPQDepth)/(t6.hybrid$totalCount+t6.hybrid$otherBases+t6.hybrid$lowMAPQDepth)

s1.hybrid.other.bases = s1.hybrid[(s1.hybrid$prop_other > 0.10), ]
s3.hybrid.other.bases = s3.hybrid[(s3.hybrid$prop_other > 0.10), ]
s6.hybrid.other.bases = s6.hybrid[(s6.hybrid$prop_other > 0.10), ]
t1.hybrid.other.bases = t1.hybrid[(t1.hybrid$prop_other > 0.10), ]
t3.hybrid.other.bases = t3.hybrid[(t3.hybrid$prop_other > 0.10), ]
t6.hybrid.other.bases = t6.hybrid[(t6.hybrid$prop_other > 0.10), ]

soma.hybrid.other.bases = rbind(s1.hybrid.other.bases,s3.hybrid.other.bases,s6.hybrid.other.bases) # merge variants failing filter in at least one replicate
testis.hybrid.other.bases = rbind(t1.hybrid.other.bases,t3.hybrid.other.bases,t6.hybrid.other.bases)
soma.hybrid.other.bases.id = unique(soma.hybrid.other.bases[c(14)]) # remove duplicate entries
testis.hybrid.other.bases.id = unique(testis.hybrid.other.bases[c(14)])
nrow(soma.hybrid.other.bases.id)
nrow(testis.hybrid.other.bases.id)

first_filter_hybrid_soma = anti_join(shared_raw_hybrid_soma,soma.hybrid.other.bases.id, by=c("variant")) # remove variants failing filter
first_filter_hybrid_testis = anti_join(shared_raw_hybrid_testis,testis.hybrid.other.bases.id, by=c("variant"))
nrow(first_filter_hybrid_soma)
nrow(first_filter_hybrid_testis)
first_filter_hybrid_soma_counts = ddply(first_filter_hybrid_soma,c("status"),summarise,N = length(variant))
first_filter_hybrid_testis_counts = ddply(first_filter_hybrid_testis,c("status"),summarise,N = length(variant))

# final filter: remove sites that are not monomorphic in P. citri or that are different from the reference

pc_asa_hybrid <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/5_ase/results/ASE_output/PC_pure_ASE_sites_check.csv",
                            "\t", escape_double = FALSE, trim_ws = TRUE) # import raw ASE counts from PC transcriptomes
pc_asa_hybrid$variant=paste(pc_asa_hybrid$contig,pc_asa_hybrid$position,sep=":")
pc_asa_hybrid$bias=pc_asa_hybrid$refCount/pc_asa_hybrid$rawDepth # estimate the proportion of citri bases/variant using raw depth rather than total count (to take other bases into account)

snp.plot.pure.PC = ggplot(pc_asa_hybrid, aes(bias)) +
  geom_histogram(breaks=seq(0,1, by=0.01), col="gray1", fill="gray1") +
  labs(title="Pure PC males", y="SNPs", x = "Prop of C bases") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  geom_vline(xintercept = 0.5,linetype="dotted", color = "gray40") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "white",size=0.01)) 
 
pc_hybrid.non.ref.95 = pc_asa_hybrid[(pc_asa_hybrid$bias < 0.95), ] # identify variants failing filter 

final_filter_hybrid_soma = anti_join(first_filter_hybrid_soma,pc_hybrid.non.ref.95, by=c("variant")) # remove variants failing filter
final_filter_hybrid_testis = anti_join(first_filter_hybrid_testis,pc_hybrid.non.ref.95, by=c("variant"))
nrow(final_filter_hybrid_soma)
nrow(final_filter_hybrid_testis)

# collect stats for final set of variants

final_filter_hybrid_soma_counts = ddply(final_filter_hybrid_soma,c("status"),summarise,N = length(variant))
final_filter_hybrid_testis_counts = ddply(final_filter_hybrid_testis,c("status"),summarise,N = length(variant))
final_filter_hybrid_soma$depth_mean = round((final_filter_hybrid_soma$totalCount.s1 + final_filter_hybrid_soma$totalCount.s3 + final_filter_hybrid_soma$totalCount.s6)/3,1)
final_filter_hybrid_testis$depth_mean = round((final_filter_hybrid_testis$totalCount.t1 + final_filter_hybrid_testis$totalCount.t3 + final_filter_hybrid_testis$totalCount.t6)/3,1)

# calculate biases to maternal allele 

final_filter_hybrid_soma$bias.s1 = round(final_filter_hybrid_soma$refCount.s1/final_filter_hybrid_soma$totalCount.s1,3)
final_filter_hybrid_soma$bias.s3 = round(final_filter_hybrid_soma$refCount.s3/final_filter_hybrid_soma$totalCount.s3,3)
final_filter_hybrid_soma$bias.s6 = round(final_filter_hybrid_soma$refCount.s6/final_filter_hybrid_soma$totalCount.s6,3)
final_filter_hybrid_soma$bias_mean = round((final_filter_hybrid_soma$bias.s1 + final_filter_hybrid_soma$bias.s3 + final_filter_hybrid_soma$bias.s6)/3,3)

final_filter_hybrid_testis$bias.t1 = round(final_filter_hybrid_testis$refCount.t1/final_filter_hybrid_testis$totalCount.t1,3)
final_filter_hybrid_testis$bias.t3 = round(final_filter_hybrid_testis$refCount.t3/final_filter_hybrid_testis$totalCount.t3,3)
final_filter_hybrid_testis$bias.t6 = round(final_filter_hybrid_testis$refCount.t6/final_filter_hybrid_testis$totalCount.t6,3)
final_filter_hybrid_testis$bias_mean = round((final_filter_hybrid_testis$bias.t1 + final_filter_hybrid_testis$bias.t3 + final_filter_hybrid_testis$bias.t6)/3,3)

read.depth.hybrid.soma = ddply(final_filter_hybrid_soma, c("status"), summarise,
                             N = length(variant), perc = length(variant)*100/nrow(final_filter_hybrid_soma),
                             mean.depth = mean(depth_mean), sd.depth = sd(depth_mean),median.depth=median(depth_mean),
                             mean.bias = mean(bias_mean), sd.bias = sd(bias_mean),median.bias=median(bias_mean))
read.depth.hybrid.testis = ddply(final_filter_hybrid_testis, c("status"), summarise,
                               N = length(variant), perc = length(variant)*100/nrow(final_filter_hybrid_testis),
                               mean.depth = mean(depth_mean), sd.depth = sd(depth_mean),median.depth=median(depth_mean),
                               mean.bias = mean(bias_mean), sd.bias = sd(bias_mean),median.bias=median(bias_mean))

# plot ASE at SNP level (histogram)

nrow(final_filter_hybrid_soma)
nrow(final_filter_hybrid_testis)
nrow(final_filter_hybrid_soma)
snp.plot.hybrid.soma = ggplot(final_filter_hybrid_soma, aes(bias_mean)) +
  geom_histogram(breaks=seq(0,1, by=0.01), col="orange1", fill="orange1") +
  labs(title="CF soma", y="SNPs", x = expression(paste(p[m]))) +
  scale_y_continuous(breaks=c(0,10000,20000,30000),limits=c(0,30000)) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  geom_vline(xintercept = 0.5,linetype="dotted", color = "gray40") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

snp.plot.hybrid.testis = ggplot(final_filter_hybrid_testis, aes(bias_mean)) +
  geom_histogram(breaks=seq(0,1, by=0.01), col="steelblue4", fill="steelblue4") +
  labs(title="CF testis", y="SNPs", x = expression(paste(p[m]))) +
  scale_y_continuous(breaks=c(0,30000,60000,90000),limits=c(0,90000)) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  geom_vline(xintercept = 0.5,linetype="dotted", color = "gray40") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

fig.1b <-grid.arrange(snp.plot.hybrid.soma,snp.plot.hybrid.testis,nrow=2)

# plot ASE bias at SNP level (barplot) according to bias category

count.hybrid.soma.mat = nrow(final_filter_hybrid_soma[(final_filter_hybrid_soma$bias_mean >= 0.950),])
count.hybrid.soma.mat.bias = nrow(final_filter_hybrid_soma[(final_filter_hybrid_soma$bias_mean > 0.650 & final_filter_hybrid_soma$bias_mean < 0.95),])
count.hybrid.soma.bip = nrow(final_filter_hybrid_soma[(final_filter_hybrid_soma$bias_mean <= 0.650 & final_filter_hybrid_soma$bias_mean >= 0.350),])
count.hybrid.soma.pat.bias = nrow(final_filter_hybrid_soma[(final_filter_hybrid_soma$bias_mean < 0.350 & final_filter_hybrid_soma$bias_mean > 0.050),])
count.hybrid.soma.pat = nrow(final_filter_hybrid_soma[(final_filter_hybrid_soma$bias_mean <= 0.050),])
total.hybrid.soma = nrow(final_filter_hybrid_soma)

final_filter_hybrid_soma.exon = final_filter_hybrid_soma[(final_filter_hybrid_soma$status == "exonic"), ]
count.hybrid.soma.mat.exon = nrow(final_filter_hybrid_soma.exon[(final_filter_hybrid_soma.exon$bias_mean >= 0.950),])
count.hybrid.soma.mat.bias.exon = nrow(final_filter_hybrid_soma.exon[(final_filter_hybrid_soma.exon$bias_mean > 0.650 & final_filter_hybrid_soma.exon$bias_mean < 0.950),])
count.hybrid.soma.bip.exon = nrow(final_filter_hybrid_soma.exon[(final_filter_hybrid_soma.exon$bias_mean <= 0.650 & final_filter_hybrid_soma.exon$bias_mean >= 0.350),])
count.hybrid.soma.pat.bias.exon = nrow(final_filter_hybrid_soma.exon[(final_filter_hybrid_soma.exon$bias_mean < 0.350 & final_filter_hybrid_soma.exon$bias_mean > 0.050),])
count.hybrid.soma.pat.exon = nrow(final_filter_hybrid_soma.exon[(final_filter_hybrid_soma.exon$bias_mean <= 0.050),])
total.hybrid.soma.exon = nrow(final_filter_hybrid_soma.exon)

bias.hybrid.soma <- data.frame("count" = c(count.hybrid.soma.mat,count.hybrid.soma.mat.exon,count.hybrid.soma.mat.bias,count.hybrid.soma.mat.bias.exon,
                                         count.hybrid.soma.bip,count.hybrid.soma.bip.exon,count.hybrid.soma.pat.bias,count.hybrid.soma.pat.bias.exon,
                                         count.hybrid.soma.pat,count.hybrid.soma.pat.exon),
                             "perc" = c(count.hybrid.soma.mat*100/total.hybrid.soma,count.hybrid.soma.mat.exon*100/total.hybrid.soma.exon,
                                        count.hybrid.soma.mat.bias*100/total.hybrid.soma,count.hybrid.soma.mat.bias.exon*100/total.hybrid.soma.exon,
                                        count.hybrid.soma.bip*100/total.hybrid.soma,count.hybrid.soma.bip.exon*100/total.hybrid.soma.exon,
                                        count.hybrid.soma.pat.bias*100/total.hybrid.soma,count.hybrid.soma.pat.bias.exon*100/total.hybrid.soma.exon,
                                        count.hybrid.soma.pat*100/total.hybrid.soma,count.hybrid.soma.pat.exon*100/total.hybrid.soma.exon),
                             "Type"=c("All","Exonic","All","Exonic","All","Exonic","All","Exonic","All","Exonic"),
                             "bias"=c("M","M","MB","MB","B","B",
                                      "PB","PB","P","P"))

count.hybrid.testis.mat = nrow(final_filter_hybrid_testis[(final_filter_hybrid_testis$bias_mean >= 0.950),])
count.hybrid.testis.mat.bias = nrow(final_filter_hybrid_testis[(final_filter_hybrid_testis$bias_mean > 0.650 & final_filter_hybrid_testis$bias_mean < 0.95),])
count.hybrid.testis.bip = nrow(final_filter_hybrid_testis[(final_filter_hybrid_testis$bias_mean <= 0.650 & final_filter_hybrid_testis$bias_mean >= 0.350),])
count.hybrid.testis.pat.bias = nrow(final_filter_hybrid_testis[(final_filter_hybrid_testis$bias_mean < 0.350 & final_filter_hybrid_testis$bias_mean > 0.050),])
count.hybrid.testis.pat = nrow(final_filter_hybrid_testis[(final_filter_hybrid_testis$bias_mean <= 0.050),])
total.hybrid.testis = nrow(final_filter_hybrid_testis)

final_filter_hybrid_testis.exon = final_filter_hybrid_testis[(final_filter_hybrid_testis$status == "exonic"), ]
count.hybrid.testis.mat.exon = nrow(final_filter_hybrid_testis.exon[(final_filter_hybrid_testis.exon$bias_mean >= 0.950),])
count.hybrid.testis.mat.bias.exon = nrow(final_filter_hybrid_testis.exon[(final_filter_hybrid_testis.exon$bias_mean > 0.650 & final_filter_hybrid_testis.exon$bias_mean < 0.950),])
count.hybrid.testis.bip.exon = nrow(final_filter_hybrid_testis.exon[(final_filter_hybrid_testis.exon$bias_mean <= 0.650 & final_filter_hybrid_testis.exon$bias_mean >= 0.350),])
count.hybrid.testis.pat.bias.exon = nrow(final_filter_hybrid_testis.exon[(final_filter_hybrid_testis.exon$bias_mean < 0.350 & final_filter_hybrid_testis.exon$bias_mean > 0.050),])
count.hybrid.testis.pat.exon = nrow(final_filter_hybrid_testis.exon[(final_filter_hybrid_testis.exon$bias_mean <= 0.050),])
total.hybrid.testis.exon = nrow(final_filter_hybrid_testis.exon)

bias.hybrid.testis <- data.frame("count" = c(count.hybrid.testis.mat,count.hybrid.testis.mat.exon,count.hybrid.testis.mat.bias,count.hybrid.testis.mat.bias.exon,
                                           count.hybrid.testis.bip,count.hybrid.testis.bip.exon,count.hybrid.testis.pat.bias,count.hybrid.testis.pat.bias.exon,
                                           count.hybrid.testis.pat,count.hybrid.testis.pat.exon),
                               "perc" = c(count.hybrid.testis.mat*100/total.hybrid.testis,count.hybrid.testis.mat.exon*100/total.hybrid.testis.exon,
                                          count.hybrid.testis.mat.bias*100/total.hybrid.testis,count.hybrid.testis.mat.bias.exon*100/total.hybrid.testis.exon,
                                          count.hybrid.testis.bip*100/total.hybrid.testis,count.hybrid.testis.bip.exon*100/total.hybrid.testis.exon,
                                          count.hybrid.testis.pat.bias*100/total.hybrid.testis,count.hybrid.testis.pat.bias.exon*100/total.hybrid.testis.exon,
                                          count.hybrid.testis.pat*100/total.hybrid.testis,count.hybrid.testis.pat.exon*100/total.hybrid.testis.exon),
                               "Type"=c("All","Exonic","All","Exonic","All","Exonic","All","Exonic","All","Exonic"),
                               "bias"=c("M","M","MB","MB","B","B",
                                        "PB","PB","P","P"))

bias.plot.hybrid.soma = ggplot(bias.hybrid.soma, aes(bias, perc, fill = Type)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(name="Type",labels = c("All", "Exonic")) +
  scale_x_discrete(limit = c("M", "MB", "B","PB","P")) +
  geom_text(aes(label=count),position = position_dodge(0.9),vjust=-1) +
  labs(title="SNP expression patterns in soma", y="\n%", x = "\nCategory of bias") +
  scale_fill_manual(values=c('orange1','moccasin')) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "black",size=0.01))

bias.plot.hybrid.testis = ggplot(bias.hybrid.testis, aes(bias, perc, fill = Type)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(name="Type",labels = c("All", "Exonic")) +
  scale_x_discrete(limit = c("M", "MB", "B","PB","P")) +
  geom_text(aes(label=count),position = position_dodge(0.9),vjust=-1) +
  labs(title="SNP expression patterns in testis", y="\n%", x = "\nCategory of bias") +
  scale_fill_manual(values=c('steelblue4','lightblue')) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "black",size=0.01))

##### Assign genes to SNPs

anno_all_unique_hybrid_exonic = anno_all_unique_hybrid[anno_all_unique_hybrid$type == "CDS",] # select exonic SNPs
anno_all_unique_hybrid_exonic_for_merge = anno_all_unique_hybrid_exonic[c(6,3)]
anno_all_unique_hybrid_exonic_for_merge$gene = sub('\\..*', '', anno_all_unique_hybrid_exonic_for_merge$anno)  # create a gene column
nrow(anno_all_unique_hybrid_exonic_for_merge) 
length(unique(anno_all_unique_hybrid_exonic_for_merge$variant))

anno_hybrid_exonic_for_merge_soma = anno_all_unique_hybrid_exonic_for_merge[anno_all_unique_hybrid_exonic_for_merge$variant %in% final_filter_hybrid_soma.exon$variant,] # keep annotation records for exonic sites passing filter
anno_hybrid_exonic_for_merge_testis = anno_all_unique_hybrid_exonic_for_merge[anno_all_unique_hybrid_exonic_for_merge$variant %in% final_filter_hybrid_testis.exon$variant,]

length(unique(anno_hybrid_exonic_for_merge_soma$variant)) - nrow(final_filter_hybrid_soma[final_filter_hybrid_soma$status == "exonic",]) # sanity check
length(unique(anno_hybrid_exonic_for_merge_testis$variant)) - nrow(final_filter_hybrid_testis[final_filter_hybrid_testis$status == "exonic",]) # sanity check

anno_hybrid_exonic_for_merge_soma[duplicated(anno_hybrid_exonic_for_merge_soma$variant),] # inspect duplicated records
anno_hybrid_exonic_for_merge_testis[duplicated(anno_hybrid_exonic_for_merge_testis$variant),]

anno_hybrid_gene_for_merge_soma=ddply(anno_hybrid_exonic_for_merge_soma,c("variant","gene"), summarize,annotation=paste(anno,collapse=",")) # collapse into genes
anno_hybrid_gene_for_merge_testis=ddply(anno_hybrid_exonic_for_merge_testis,c("variant","gene"), summarize,annotation=paste(anno,collapse=","))

nrow(anno_hybrid_gene_for_merge_soma)
nrow(anno_hybrid_exonic_for_merge_soma)
length(unique(anno_hybrid_exonic_for_merge_soma$variant)) # unique variants
length(unique(anno_hybrid_gene_for_merge_soma$variant))
anno_hybrid_gene_for_merge_soma[duplicated(anno_hybrid_gene_for_merge_soma$variant),] # inspect duplicated records (SNPs that map to different genes)

nrow(anno_hybrid_gene_for_merge_testis)
nrow(anno_hybrid_exonic_for_merge_testis)
length(unique(anno_hybrid_exonic_for_merge_testis$variant))
length(unique(anno_hybrid_gene_for_merge_testis$variant))
anno_hybrid_gene_for_merge_testis[duplicated(anno_hybrid_gene_for_merge_testis$variant),]

# intersect ASE counts and gene information

final_filter_hybrid_soma_exon_for_merge = final_filter_hybrid_soma.exon[c(1,2,6,7,8,9,10,11,12,13,14)]
final_filter_hybrid_testis_exon_for_merge = final_filter_hybrid_testis.exon[c(1,2,6,7,8,9,10,11,12,13,14)]

variants_with_genes_hybrid_soma = full_join(final_filter_hybrid_soma_exon_for_merge,anno_hybrid_gene_for_merge_soma,by="variant")
variants_with_genes_hybrid_testis = full_join(final_filter_hybrid_testis_exon_for_merge,anno_hybrid_gene_for_merge_testis,by="variant")
nrow(variants_with_genes_hybrid_soma)
nrow(variants_with_genes_hybrid_testis)
length(unique(variants_with_genes_hybrid_soma$variant))
length(unique(variants_with_genes_hybrid_testis$variant))

# collect ASE counts at gene level

genes_raw_hybrid_soma = ddply(variants_with_genes_hybrid_soma, c("gene"), summarise,
                            ref.s1 = sum(refCount.s1),
                            alt.s1 = sum(altCount.s1),
                            total.s1 = sum(totalCount.s1),
                            ref.s3 = sum(refCount.s3),
                            alt.s3 = sum(altCount.s3),
                            total.s3 = sum(totalCount.s3),
                            ref.s6 = sum(refCount.s6),
                            alt.s6 = sum(altCount.s6),
                            total.s6 = sum(totalCount.s6),
                            N = length(variant))

genes_raw_hybrid_testis = ddply(variants_with_genes_hybrid_testis, c("gene"), summarise,
                              ref.t1 = sum(refCount.t1),
                              alt.t1 = sum(altCount.t1),
                              total.t1 = sum(totalCount.t1),
                              ref.t3 = sum(refCount.t3),
                              alt.t3 = sum(altCount.t3),
                              total.t3 = sum(totalCount.t3),
                              ref.t6 = sum(refCount.t6),
                              alt.t6 = sum(altCount.t6),
                              total.t6 = sum(totalCount.t6),
                              N = length(variant))

nrow(genes_raw_hybrid_soma)
nrow(genes_raw_hybrid_testis)

##### Filter genes

## SOMA

# first filter: remove genes covered by a single site unless average read depth across replicates > 100

genes_raw_hybrid_soma_one_SNP = genes_raw_hybrid_soma[genes_raw_hybrid_soma$N == 1,] # select genes covered by a single SNP
genes_raw_hybrid_soma_one_SNP$depth_mean = round((genes_raw_hybrid_soma_one_SNP$total.s1 + genes_raw_hybrid_soma_one_SNP$total.s3 + genes_raw_hybrid_soma_one_SNP$total.s6)/3,1) # obtain average read depth
genes_raw_hybrid_soma_one_SNP_less_than_100 = genes_raw_hybrid_soma_one_SNP[(genes_raw_hybrid_soma_one_SNP$depth_mean < 100), ]
nrow(genes_raw_hybrid_soma_one_SNP_less_than_100)
genes_first_filter_hybrid_soma = anti_join(genes_raw_hybrid_soma,genes_raw_hybrid_soma_one_SNP_less_than_100,by=c("gene")) # remove genes failing filter
nrow(genes_first_filter_hybrid_soma)

# second filter: G-test of independence to remove genes with high heterogeneity across replicates

library(DescTools)
library("devtools")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("qvalue")

genes_hybrid_soma_test = genes_first_filter_hybrid_soma
nrow(genes_hybrid_soma_test)
genes_hybrid_soma_test$G = 0
genes_hybrid_soma_test$p.value = 0

for(each_gene in 1:nrow(genes_hybrid_soma_test)){
  stat = GTest(as.table(rbind(c(genes_hybrid_soma_test[each_gene,]$ref.s1,genes_hybrid_soma_test[each_gene,]$ref.s3,genes_hybrid_soma_test[each_gene,]$ref.s6),
                              c(genes_hybrid_soma_test[each_gene,]$alt.s1,genes_hybrid_soma_test[each_gene,]$alt.s3,genes_hybrid_soma_test[each_gene,]$alt.s6))))$statistic
  genes_hybrid_soma_test$G[each_gene] = stat
}
for(each_gene in 1:nrow(genes_hybrid_soma_test)){
  p = GTest(as.table(rbind(c(genes_hybrid_soma_test[each_gene,]$ref.s1,genes_hybrid_soma_test[each_gene,]$ref.s3,genes_hybrid_soma_test[each_gene,]$ref.s6),
                           c(genes_hybrid_soma_test[each_gene,]$alt.s1,genes_hybrid_soma_test[each_gene,]$alt.s3,genes_hybrid_soma_test[each_gene,]$alt.s6))))$p.value
  genes_hybrid_soma_test$p.value[each_gene] = p
}

# compute exact binomial tests for each replicate

genes_hybrid_soma_test$s1.p.value = 0
genes_hybrid_soma_test$s3.p.value = 0
genes_hybrid_soma_test$s6.p.value = 0
tail(genes_hybrid_soma_test)

for(each_gene in 1:nrow(genes_hybrid_soma_test)){
  a = binom.test(genes_hybrid_soma_test[each_gene,]$ref.s1,genes_hybrid_soma_test[each_gene,]$total.s1,0.5,alternative="two.sided")$p.value
  genes_hybrid_soma_test$s1.p.value[each_gene] = a
}
for(each_gene in 1:nrow(genes_hybrid_soma_test)){
  b = binom.test(genes_hybrid_soma_test[each_gene,]$ref.s3,genes_hybrid_soma_test[each_gene,]$total.s3,0.5,alternative="two.sided")$p.value
  genes_hybrid_soma_test$s3.p.value[each_gene] = b
}
for(each_gene in 1:nrow(genes_hybrid_soma_test)){
  c = binom.test(genes_hybrid_soma_test[each_gene,]$ref.s6,genes_hybrid_soma_test[each_gene,]$total.s6,0.5,alternative="two.sided")$p.value
  genes_hybrid_soma_test$s6.p.value[each_gene] = c
}

# determine whether genes pass or fail the heterogeneity (Bonferroni correction)

genes_hybrid_soma_test_pass_G = genes_hybrid_soma_test[(genes_hybrid_soma_test$p.value >= 0.05/nrow(genes_hybrid_soma_test)), ] # not significantly heterogeneous
nrow(genes_hybrid_soma_test_pass_G)
genes_hybrid_soma_test_fail_G = genes_hybrid_soma_test[(genes_hybrid_soma_test$p.value <  0.05/nrow(genes_hybrid_soma_test)), ] # significantly heterogeneous
nrow(genes_hybrid_soma_test_fail_G) + nrow(genes_hybrid_soma_test_pass_G)

# rescue genes that failed G test if all replicates agree in significance of exact binomial test (Bonf. corrected) and bias category

genes_hybrid_soma_test_fail_G$all.p.significant = ifelse(((genes_hybrid_soma_test_fail_G$s1.p.value < 0.05/(3*(nrow(genes_hybrid_soma_test)))) &
                                                            (genes_hybrid_soma_test_fail_G$s3.p.value < 0.05/(3*(nrow(genes_hybrid_soma_test)))) &
                                                            (genes_hybrid_soma_test_fail_G$s6.p.value < 0.05/(3*(nrow(genes_hybrid_soma_test))))) |
                                                         ((genes_hybrid_soma_test_fail_G$s1.p.value > 0.05/(3*(nrow(genes_hybrid_soma_test)))) &
                                                            (genes_hybrid_soma_test_fail_G$s3.p.value > 0.05/(3*(nrow(genes_hybrid_soma_test)))) &
                                                            (genes_hybrid_soma_test_fail_G$s6.p.value > 0.05/(3*(nrow(genes_hybrid_soma_test))))), "AGREE", "DISAGREE")

nrow(genes_hybrid_soma_test_fail_G[(genes_hybrid_soma_test_fail_G$all.p.significant == "DISAGREE"), ])

genes_hybrid_soma_test_fail_G$bias.s1 = genes_hybrid_soma_test_fail_G$ref.s1 / genes_hybrid_soma_test_fail_G$total.s1
genes_hybrid_soma_test_fail_G$bias.s3 = genes_hybrid_soma_test_fail_G$ref.s3 / genes_hybrid_soma_test_fail_G$total.s3
genes_hybrid_soma_test_fail_G$bias.s6 = genes_hybrid_soma_test_fail_G$ref.s6 / genes_hybrid_soma_test_fail_G$total.s6

genes_hybrid_soma_test_fail_G$bias.cat.s1 = "B"
genes_hybrid_soma_test_fail_G$bias.cat.s3 = "B"
genes_hybrid_soma_test_fail_G$bias.cat.s6 = "B"

genes_hybrid_soma_test_fail_G$bias.cat.s1[genes_hybrid_soma_test_fail_G$bias.s1<=0.050] <-"P"
genes_hybrid_soma_test_fail_G$bias.cat.s1[genes_hybrid_soma_test_fail_G$bias.s1>0.050 & genes_hybrid_soma_test_fail_G$bias.s1<0.350]<-"PB"
genes_hybrid_soma_test_fail_G$bias.cat.s1[genes_hybrid_soma_test_fail_G$bias.s1>=0.950]<-"M"
genes_hybrid_soma_test_fail_G$bias.cat.s1[genes_hybrid_soma_test_fail_G$bias.s1>0.650 & genes_hybrid_soma_test_fail_G$bias.s1<0.950]<-"MB"

genes_hybrid_soma_test_fail_G$bias.cat.s3[genes_hybrid_soma_test_fail_G$bias.s3<=0.050] <-"P"
genes_hybrid_soma_test_fail_G$bias.cat.s3[genes_hybrid_soma_test_fail_G$bias.s3>0.050 & genes_hybrid_soma_test_fail_G$bias.s3<0.350]<-"PB"
genes_hybrid_soma_test_fail_G$bias.cat.s3[genes_hybrid_soma_test_fail_G$bias.s3>=0.950]<-"M"
genes_hybrid_soma_test_fail_G$bias.cat.s3[genes_hybrid_soma_test_fail_G$bias.s3>0.650 & genes_hybrid_soma_test_fail_G$bias.s3<0.950]<-"MB"

genes_hybrid_soma_test_fail_G$bias.cat.s6[genes_hybrid_soma_test_fail_G$bias.s6<=0.050] <-"P"
genes_hybrid_soma_test_fail_G$bias.cat.s6[genes_hybrid_soma_test_fail_G$bias.s6>0.050 & genes_hybrid_soma_test_fail_G$bias.s6<0.350]<-"PB"
genes_hybrid_soma_test_fail_G$bias.cat.s6[genes_hybrid_soma_test_fail_G$bias.s6>=0.950]<-"M"
genes_hybrid_soma_test_fail_G$bias.cat.s6[genes_hybrid_soma_test_fail_G$bias.s6>0.650 & genes_hybrid_soma_test_fail_G$bias.s6<0.950]<-"MB"

genes_hybrid_soma_test_fail_G$all.bias.merge = paste(genes_hybrid_soma_test_fail_G$bias.cat.s1, genes_hybrid_soma_test_fail_G$bias.cat.s3, genes_hybrid_soma_test_fail_G$bias.cat.s6 , sep='') 
genes_hybrid_soma_test_fail_G$all.bias.cat = ifelse((genes_hybrid_soma_test_fail_G$bias.cat.s1 == "P" & genes_hybrid_soma_test_fail_G$bias.cat.s3 == "P" & genes_hybrid_soma_test_fail_G$bias.cat.s6 == "P") |
                                                    (genes_hybrid_soma_test_fail_G$bias.cat.s1 == "PB" & genes_hybrid_soma_test_fail_G$bias.cat.s3 == "PB" & genes_hybrid_soma_test_fail_G$bias.cat.s6 == "PB") |
                                                    (genes_hybrid_soma_test_fail_G$bias.cat.s1 == "B" & genes_hybrid_soma_test_fail_G$bias.cat.s3 == "B" & genes_hybrid_soma_test_fail_G$bias.cat.s6 == "B") |
                                                    (genes_hybrid_soma_test_fail_G$bias.cat.s1 == "MB" & genes_hybrid_soma_test_fail_G$bias.cat.s3 == "MB" & genes_hybrid_soma_test_fail_G$bias.cat.s6 == "MB") |
                                                    (genes_hybrid_soma_test_fail_G$bias.cat.s1 == "M" & genes_hybrid_soma_test_fail_G$bias.cat.s3 == "M" & genes_hybrid_soma_test_fail_G$bias.cat.s6 == "M"),"AGREE","DISAGREE")


genes_hybrid_soma_test_fail_G$decision = ifelse((genes_hybrid_soma_test_fail_G$all.p.significant == "AGREE") & (genes_hybrid_soma_test_fail_G$all.bias.cat == "AGREE"), "KEEP", "DROP")
genes_hybrid_soma_test_fail_filter = genes_hybrid_soma_test_fail_G[(genes_hybrid_soma_test_fail_G$decision == "DROP"), ]
genes_second_filter_hybrid_soma = anti_join(genes_first_filter_hybrid_soma, genes_hybrid_soma_test_fail_filter, by=c("gene"))
nrow(genes_hybrid_soma_test_fail_filter)
nrow(genes_second_filter_hybrid_soma)

# final filter: remove genes with TPM < 1

genes_final_filter_hybrid_soma = anti_join(genes_second_filter_hybrid_soma, s_tpm_not_expressed_join, by=c("gene"))
nrow(genes_final_filter_hybrid_soma)

## TESTIS

genes_raw_hybrid_testis_one_SNP = genes_raw_hybrid_testis[genes_raw_hybrid_testis$N == 1,]
genes_raw_hybrid_testis_one_SNP$depth_mean = round((genes_raw_hybrid_testis_one_SNP$total.t1 + genes_raw_hybrid_testis_one_SNP$total.t3 + genes_raw_hybrid_testis_one_SNP$total.t6)/3,1)
genes_raw_hybrid_testis_one_SNP_less_than_100 = genes_raw_hybrid_testis_one_SNP[(genes_raw_hybrid_testis_one_SNP$depth_mean < 100), ]
nrow(genes_raw_hybrid_testis_one_SNP_less_than_100)
genes_first_filter_hybrid_testis = anti_join(genes_raw_hybrid_testis,genes_raw_hybrid_testis_one_SNP_less_than_100,by=c("gene"))
nrow(genes_first_filter_hybrid_testis)

genes_hybrid_testis_test = genes_first_filter_hybrid_testis
head(genes_hybrid_testis_test)
genes_hybrid_testis_test$G = 0
genes_hybrid_testis_test$p.value = 0

for(each_gene in 1:nrow(genes_hybrid_testis_test)){
  stat = GTest(as.table(rbind(c(genes_hybrid_testis_test[each_gene,]$ref.t1,genes_hybrid_testis_test[each_gene,]$ref.t3,genes_hybrid_testis_test[each_gene,]$ref.t6),
                              c(genes_hybrid_testis_test[each_gene,]$alt.t1,genes_hybrid_testis_test[each_gene,]$alt.t3,genes_hybrid_testis_test[each_gene,]$alt.t6))))$statistic
  genes_hybrid_testis_test$G[each_gene] = stat
}
for(each_gene in 1:nrow(genes_hybrid_testis_test)){
  p = GTest(as.table(rbind(c(genes_hybrid_testis_test[each_gene,]$ref.t1,genes_hybrid_testis_test[each_gene,]$ref.t3,genes_hybrid_testis_test[each_gene,]$ref.t6),
                           c(genes_hybrid_testis_test[each_gene,]$alt.t1,genes_hybrid_testis_test[each_gene,]$alt.t3,genes_hybrid_testis_test[each_gene,]$alt.t6))))$p.value
  genes_hybrid_testis_test$p.value[each_gene] = p
}

genes_hybrid_testis_test$t1.p.value = 0
genes_hybrid_testis_test$t3.p.value = 0
genes_hybrid_testis_test$t6.p.value = 0
tail(genes_hybrid_testis_test)

for(each_gene in 1:nrow(genes_hybrid_testis_test)){
  a = binom.test(genes_hybrid_testis_test[each_gene,]$ref.t1,genes_hybrid_testis_test[each_gene,]$total.t1,0.5,alternative="two.sided")$p.value
  genes_hybrid_testis_test$t1.p.value[each_gene] = a
}
for(each_gene in 1:nrow(genes_hybrid_testis_test)){
  b = binom.test(genes_hybrid_testis_test[each_gene,]$ref.t3,genes_hybrid_testis_test[each_gene,]$total.t3,0.5,alternative="two.sided")$p.value
  genes_hybrid_testis_test$t3.p.value[each_gene] = b
}
for(each_gene in 1:nrow(genes_hybrid_testis_test)){
  c = binom.test(genes_hybrid_testis_test[each_gene,]$ref.t6,genes_hybrid_testis_test[each_gene,]$total.t6,0.5,alternative="two.sided")$p.value
  genes_hybrid_testis_test$t6.p.value[each_gene] = c
}

genes_hybrid_testis_test_pass_G = genes_hybrid_testis_test[(genes_hybrid_testis_test$p.value >= 0.05/nrow(genes_hybrid_testis_test)), ]
nrow(genes_hybrid_testis_test_pass_G)
genes_hybrid_testis_test_fail_G = genes_hybrid_testis_test[(genes_hybrid_testis_test$p.value <  0.05/nrow(genes_hybrid_testis_test)), ]
nrow(genes_hybrid_testis_test_fail_G) + nrow(genes_hybrid_testis_test_pass_G)

genes_hybrid_testis_test_fail_G$all.p.significant = ifelse(((genes_hybrid_testis_test_fail_G$t1.p.value < 0.05/(3*(nrow(genes_hybrid_testis_test)))) & (genes_hybrid_testis_test_fail_G$t3.p.value < 0.05/(3*(nrow(genes_hybrid_testis_test)))) & (genes_hybrid_testis_test_fail_G$t6.p.value < 0.05/(3*(nrow(genes_hybrid_testis_test))))) |
                                                           ((genes_hybrid_testis_test_fail_G$t1.p.value > 0.05/(3*(nrow(genes_hybrid_testis_test)))) & (genes_hybrid_testis_test_fail_G$t3.p.value > 0.05/(3*(nrow(genes_hybrid_testis_test)))) & (genes_hybrid_testis_test_fail_G$t6.p.value > 0.05/(3*(nrow(genes_hybrid_testis_test))))), "AGREE", "DISAGREE")

nrow(genes_hybrid_testis_test_fail_G[(genes_hybrid_testis_test_fail_G$all.p.significant == "DISAGREE"), ])

genes_hybrid_testis_test_fail_G$bias.t1 = genes_hybrid_testis_test_fail_G$ref.t1 / genes_hybrid_testis_test_fail_G$total.t1
genes_hybrid_testis_test_fail_G$bias.t3 = genes_hybrid_testis_test_fail_G$ref.t3 / genes_hybrid_testis_test_fail_G$total.t3
genes_hybrid_testis_test_fail_G$bias.t6 = genes_hybrid_testis_test_fail_G$ref.t6 / genes_hybrid_testis_test_fail_G$total.t6

genes_hybrid_testis_test_fail_G$bias.cat.t1 = "B"
genes_hybrid_testis_test_fail_G$bias.cat.t3 = "B"
genes_hybrid_testis_test_fail_G$bias.cat.t6 = "B"

genes_hybrid_testis_test_fail_G$bias.cat.t1[genes_hybrid_testis_test_fail_G$bias.t1<=0.050] <-"P"
genes_hybrid_testis_test_fail_G$bias.cat.t1[genes_hybrid_testis_test_fail_G$bias.t1>0.050 & genes_hybrid_testis_test_fail_G$bias.t1<0.350]<-"PB"
genes_hybrid_testis_test_fail_G$bias.cat.t1[genes_hybrid_testis_test_fail_G$bias.t1>=0.950]<-"M"
genes_hybrid_testis_test_fail_G$bias.cat.t1[genes_hybrid_testis_test_fail_G$bias.t1>0.650 & genes_hybrid_testis_test_fail_G$bias.t1<0.950]<-"MB"

genes_hybrid_testis_test_fail_G$bias.cat.t3[genes_hybrid_testis_test_fail_G$bias.t3<=0.050] <-"P"
genes_hybrid_testis_test_fail_G$bias.cat.t3[genes_hybrid_testis_test_fail_G$bias.t3>0.050 & genes_hybrid_testis_test_fail_G$bias.t3<0.350]<-"PB"
genes_hybrid_testis_test_fail_G$bias.cat.t3[genes_hybrid_testis_test_fail_G$bias.t3>=0.950]<-"M"
genes_hybrid_testis_test_fail_G$bias.cat.t3[genes_hybrid_testis_test_fail_G$bias.t3>0.650 & genes_hybrid_testis_test_fail_G$bias.t3<0.950]<-"MB"

genes_hybrid_testis_test_fail_G$bias.cat.t6[genes_hybrid_testis_test_fail_G$bias.t6<=0.050] <-"P"
genes_hybrid_testis_test_fail_G$bias.cat.t6[genes_hybrid_testis_test_fail_G$bias.t6>0.050 & genes_hybrid_testis_test_fail_G$bias.t6<0.350]<-"PB"
genes_hybrid_testis_test_fail_G$bias.cat.t6[genes_hybrid_testis_test_fail_G$bias.t6>=0.950]<-"M"
genes_hybrid_testis_test_fail_G$bias.cat.t6[genes_hybrid_testis_test_fail_G$bias.t6>0.650 & genes_hybrid_testis_test_fail_G$bias.t6<0.950]<-"MB"

genes_hybrid_testis_test_fail_G$all.bias.merge = paste(genes_hybrid_testis_test_fail_G$bias.cat.t1, genes_hybrid_testis_test_fail_G$bias.cat.t3, genes_hybrid_testis_test_fail_G$bias.cat.t6 , sep='') 
genes_hybrid_testis_test_fail_G$all.bias.cat = ifelse((genes_hybrid_testis_test_fail_G$bias.cat.t1 == "P" & genes_hybrid_testis_test_fail_G$bias.cat.t3 == "P" & genes_hybrid_testis_test_fail_G$bias.cat.t6 == "P") |
                                                      (genes_hybrid_testis_test_fail_G$bias.cat.t1 == "PB" & genes_hybrid_testis_test_fail_G$bias.cat.t3 == "PB" & genes_hybrid_testis_test_fail_G$bias.cat.t6 == "PB") |
                                                      (genes_hybrid_testis_test_fail_G$bias.cat.t1 == "B" & genes_hybrid_testis_test_fail_G$bias.cat.t3 == "B" & genes_hybrid_testis_test_fail_G$bias.cat.t6 == "B") |
                                                      (genes_hybrid_testis_test_fail_G$bias.cat.t1 == "MB" & genes_hybrid_testis_test_fail_G$bias.cat.t3 == "MB" & genes_hybrid_testis_test_fail_G$bias.cat.t6 == "MB") |
                                                      (genes_hybrid_testis_test_fail_G$bias.cat.t1 == "M" & genes_hybrid_testis_test_fail_G$bias.cat.t3 == "M" & genes_hybrid_testis_test_fail_G$bias.cat.t6 == "M"),"AGREE","DISAGREE")


genes_hybrid_testis_test_fail_G$decision = ifelse((genes_hybrid_testis_test_fail_G$all.p.significant == "AGREE") & (genes_hybrid_testis_test_fail_G$all.bias.cat == "AGREE"), "KEEP", "DROP")
genes_hybrid_testis_test_fail_filter = genes_hybrid_testis_test_fail_G[(genes_hybrid_testis_test_fail_G$decision == "DROP"), ]
genes_second_filter_hybrid_testis = anti_join(genes_first_filter_hybrid_testis, genes_hybrid_testis_test_fail_filter, by=c("gene"))

genes_final_filter_hybrid_testis = anti_join(genes_second_filter_hybrid_testis, t_tpm_not_expressed_join, by=c("gene"))
nrow(genes_final_filter_hybrid_testis)

##### Test for imprinted genes and assign category of bias

# SOMA

# pool counts across replicates and calculate gene bias

genes_final_filter_hybrid_soma$ref.all= genes_final_filter_hybrid_soma$ref.s1 + genes_final_filter_hybrid_soma$ref.s3 + genes_final_filter_hybrid_soma$ref.s6
genes_final_filter_hybrid_soma$alt.all= genes_final_filter_hybrid_soma$alt.s1 + genes_final_filter_hybrid_soma$alt.s3 + genes_final_filter_hybrid_soma$alt.s6
genes_final_filter_hybrid_soma$total.all = genes_final_filter_hybrid_soma$ref.all + genes_final_filter_hybrid_soma$alt.all
genes_final_filter_hybrid_soma$bias = genes_final_filter_hybrid_soma$ref.all / genes_final_filter_hybrid_soma$total.all
genes_final_filter_hybrid_soma$exact = 0
genes_final_filter_hybrid_soma$p.value = 0

# repeat exact binomial test with pooled counts

for(each_gene in 1:nrow(genes_final_filter_hybrid_soma)){
  zz = binom.test(genes_final_filter_hybrid_soma[each_gene,]$ref.all,genes_final_filter_hybrid_soma[each_gene,]$total.all,0.5,alternative="two.sided")$statistic
  genes_final_filter_hybrid_soma$exact[each_gene] = zz
}
for(each_gene in 1:nrow(genes_final_filter_hybrid_soma)){
  yy = binom.test(genes_final_filter_hybrid_soma[each_gene,]$ref.all,genes_final_filter_hybrid_soma[each_gene,]$total.all,0.5,alternative="two.sided")$p.value
  genes_final_filter_hybrid_soma$p.value[each_gene] = yy
}

genes_final_filter_hybrid_soma$bias.yn = 0
genes_final_filter_hybrid_soma$bias.yn = ifelse(genes_final_filter_hybrid_soma$p.value>(0.05/nrow(genes_final_filter_hybrid_soma)), "NS", "S")

genes_final_filter_hybrid_soma$bias.cat = "biparental"
genes_final_filter_hybrid_soma$bias.cat[genes_final_filter_hybrid_soma$bias.yn=="S" & genes_final_filter_hybrid_soma$bias<=0.050]<-"paternal.only"
genes_final_filter_hybrid_soma$bias.cat[genes_final_filter_hybrid_soma$bias.yn=="S" & genes_final_filter_hybrid_soma$bias>0.050 & genes_final_filter_hybrid_soma$bias<0.350]<-"paternal.bias"
genes_final_filter_hybrid_soma$bias.cat[genes_final_filter_hybrid_soma$bias.yn=="S" & genes_final_filter_hybrid_soma$bias>=0.950]<-"maternal.only"
genes_final_filter_hybrid_soma$bias.cat[genes_final_filter_hybrid_soma$bias.yn=="S" & genes_final_filter_hybrid_soma$bias>0.650 & genes_final_filter_hybrid_soma$bias<0.950]<-"maternal.bias"

genes_final_filter_hybrid_soma_counts = ddply(genes_final_filter_hybrid_soma, c("bias.cat"), summarise,
                                            N = length(gene),
                                            perc = length(gene)*100/nrow(genes_final_filter_hybrid_soma))


## TESTIS

genes_final_filter_hybrid_testis$ref.all= genes_final_filter_hybrid_testis$ref.t1 + genes_final_filter_hybrid_testis$ref.t3 + genes_final_filter_hybrid_testis$ref.t6
genes_final_filter_hybrid_testis$alt.all= genes_final_filter_hybrid_testis$alt.t1 + genes_final_filter_hybrid_testis$alt.t3 + genes_final_filter_hybrid_testis$alt.t6
genes_final_filter_hybrid_testis$total.all = genes_final_filter_hybrid_testis$ref.all + genes_final_filter_hybrid_testis$alt.all
genes_final_filter_hybrid_testis$bias = genes_final_filter_hybrid_testis$ref.all / genes_final_filter_hybrid_testis$total.all
genes_final_filter_hybrid_testis$exact = 0
genes_final_filter_hybrid_testis$p.value = 0

for(each_gene in 1:nrow(genes_final_filter_hybrid_testis)){
  zz = binom.test(genes_final_filter_hybrid_testis[each_gene,]$ref.all,genes_final_filter_hybrid_testis[each_gene,]$total.all,0.5,alternative="two.sided")$statistic
  genes_final_filter_hybrid_testis$exact[each_gene] = zz
}
for(each_gene in 1:nrow(genes_final_filter_hybrid_testis)){
  yy = binom.test(genes_final_filter_hybrid_testis[each_gene,]$ref.all,genes_final_filter_hybrid_testis[each_gene,]$total.all,0.5,alternative="two.sided")$p.value
  genes_final_filter_hybrid_testis$p.value[each_gene] = yy
}

genes_final_filter_hybrid_testis$bias.yn = 0
genes_final_filter_hybrid_testis$bias.yn = ifelse(genes_final_filter_hybrid_testis$p.value>(0.05/nrow(genes_final_filter_hybrid_testis)), "NS", "S")

genes_final_filter_hybrid_testis$bias.cat = "biparental"
genes_final_filter_hybrid_testis$bias.cat[genes_final_filter_hybrid_testis$bias.yn=="S" & genes_final_filter_hybrid_testis$bias<=0.050]<-"paternal.only"
genes_final_filter_hybrid_testis$bias.cat[genes_final_filter_hybrid_testis$bias.yn=="S" & genes_final_filter_hybrid_testis$bias>0.050 & genes_final_filter_hybrid_testis$bias<0.350]<-"paternal.bias"
genes_final_filter_hybrid_testis$bias.cat[genes_final_filter_hybrid_testis$bias.yn=="S" & genes_final_filter_hybrid_testis$bias>=0.950]<-"maternal.only"
genes_final_filter_hybrid_testis$bias.cat[genes_final_filter_hybrid_testis$bias.yn=="S" & genes_final_filter_hybrid_testis$bias>0.650 & genes_final_filter_hybrid_testis$bias<0.950]<-"maternal.bias"

genes_final_filter_hybrid_testis_counts = ddply(genes_final_filter_hybrid_testis, c("bias.cat"), summarise,
                                              N = length(gene),
                                              perc = length(gene)*100/nrow(genes_final_filter_hybrid_testis))

###### export

#write.csv(genes_final_filter_hybrid_soma, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_final_filter_hybrid_soma.csv")
#write.csv(genes_final_filter_hybrid_testis, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_final_filter_hybrid_testis.csv")

###### export for go

genes_biparental_pat_soma_hybrids <- genes_final_filter_hybrid_soma[genes_final_filter_hybrid_soma$bias.cat == 'biparental' | genes_final_filter_hybrid_soma$bias.cat == 'paternal.bias', ]
genes_biparental_pat_testis_hybrids <- genes_final_filter_hybrid_testis[genes_final_filter_hybrid_testis$bias.cat == 'biparental' | genes_final_filter_hybrid_testis$bias.cat == 'paternal.bias', ]
genes_nonmat_soma_hybrids <- genes_final_filter_hybrid_soma[genes_final_filter_hybrid_soma$bias.cat != 'maternal.only', ]
genes_nonmat_testis_hybrids <- genes_final_filter_hybrid_testis[genes_final_filter_hybrid_testis$bias.cat != 'maternal.only', ]

#write.csv(genes_biparental_pat_soma_hybrids, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_biparental_pat_soma_hybrids.csv")
#write.csv(genes_biparental_pat_testis_hybrids, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_biparental_pat_testis_hybrids.csv")
#write.csv(genes_nonmat_soma_hybrids, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_nonmat_soma_hybrids.csv")
#write.csv(genes_nonmat_testis_hybrids, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_nonmat_testis_hybrids.csv")

###### plot

genes_hybrid_soma_pat_only = genes_final_filter_hybrid_soma_export[(genes_final_filter_hybrid_soma_export$bias.cat == "paternal.only"), ]
genes_hybrid_soma_pat_bias = genes_final_filter_hybrid_soma_export[(genes_final_filter_hybrid_soma_export$bias.cat == "paternal.bias"), ]
genes_hybrid_soma_bip = genes_final_filter_hybrid_soma_export[(genes_final_filter_hybrid_soma_export$bias.cat == "biparental"), ]
genes_hybrid_soma_mat_bias = genes_final_filter_hybrid_soma_export[(genes_final_filter_hybrid_soma_export$bias.cat == "maternal.bias"), ]
genes_hybrid_soma_mat_only = genes_final_filter_hybrid_soma_export[(genes_final_filter_hybrid_soma_export$bias.cat == "maternal.only"), ]

genes_hybrid_testis_pat_only = genes_final_filter_hybrid_testis_export[(genes_final_filter_hybrid_testis_export$bias.cat == "paternal.only"), ]
genes_hybrid_testis_pat_bias = genes_final_filter_hybrid_testis_export[(genes_final_filter_hybrid_testis_export$bias.cat == "paternal.bias"), ]
genes_hybrid_testis_bip = genes_final_filter_hybrid_testis_export[(genes_final_filter_hybrid_testis_export$bias.cat == "biparental"), ]
genes_hybrid_testis_mat_bias = genes_final_filter_hybrid_testis_export[(genes_final_filter_hybrid_testis_export$bias.cat == "maternal.bias"), ]
genes_hybrid_testis_mat_only = genes_final_filter_hybrid_testis_export[(genes_final_filter_hybrid_testis_export$bias.cat == "maternal.only"), ]

genes.hybrid.soma.counts.plot <- ggplot(genes_final_filter_hybrid_soma_counts, aes(bias.cat, perc)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge(),colour="black",fill="orange1") +
  scale_x_discrete(limit = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"),
                   labels = c("M","MB","B","PB","P")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  labs(title="CF soma", y="%", x = "Category of bias") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

genes.hybrid.testis.counts.plot <- ggplot(genes_final_filter_hybrid_testis_counts, aes(bias.cat, perc)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge(),colour="black",fill="steelblue4") +
  scale_x_discrete(limit = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"),
                   labels = c("M","MB","B","PB","P")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  labs(title="CF testis", y="%", x = "Category of bias") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

genes.soma.scatterplot = ggplot(genes_final_filter_hybrid_soma, aes(x=(1+ref.all), y=(1+alt.all),group=bias.cat)) + geom_point(aes(colour=bias.cat,shape=bias.cat),size=1.8) +
  scale_x_log10(limits = c(1,10000000),breaks=c(1,10,100,1000,10000,100000,1000000,10000000),
                labels=c(bquote(0^""),10,bquote(10^2),bquote(10^3),bquote(10^4),bquote(10^5),bquote(10^6),bquote(10^7))) +
  scale_y_log10(limits = c(1,10000000),breaks=c(1,10,100,1000,10000,100000,1000000,10000000),
                labels=c(bquote(0^""),10,bquote(10^2),bquote(10^3),bquote(10^4),bquote(10^5),bquote(10^6),bquote(10^7))) +
  scale_shape_manual(name="Bias",values=c(17,16,16,15,15),
                     breaks=c("maternal.only","maternal.bias","biparental","paternal.bias","paternal.only"),
                     labels=c("M", "MB", "B","PB","P")) + 
  scale_colour_manual(name="Bias",values=c("#B9770E","#F8C471","#FAD7A0","#F39C12","#D35400"),
                      breaks=c("maternal.only","maternal.bias","biparental","paternal.bias","paternal.only"),
                      labels=c("M", "MB", "B","PB","P")) +
  labs(title="CF soma", y="Paternal SNP counts", x = "Maternal SNP counts") +
  geom_abline(slope = 1,linetype="dotted", color = "gray50") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))  +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

genes.testis.scatterplot = ggplot(genes_final_filter_hybrid_testis, aes(x=(1+ref.all), y=(1+alt.all),group=bias.cat)) + geom_point(aes(colour=bias.cat,shape=bias.cat),size=1.8) +
  scale_x_log10(limits = c(1,10000000),breaks=c(1,10,100,1000,10000,100000,1000000,10000000),
                labels=c(bquote(0^""),10,bquote(10^2),bquote(10^3),bquote(10^4),bquote(10^5),bquote(10^6),bquote(10^7))) +
  scale_y_log10(limits = c(1,10000000),breaks=c(1,10,100,1000,10000,100000,1000000,10000000),
                labels=c(bquote(0^""),10,bquote(10^2),bquote(10^3),bquote(10^4),bquote(10^5),bquote(10^6),bquote(10^7))) +
  scale_shape_manual(name="Bias",values=c(17,16,16,15,15), 
                     breaks=c("maternal.only","maternal.bias","biparental","paternal.bias","paternal.only"),
                     labels=c("M", "MB", "B","PB","P")) + 
  scale_colour_manual(name="Bias",values=c("#2980B9","#85C1E9","#AED6F1","#3498DB","#2471A3"),
                      breaks=c("maternal.only","maternal.bias","biparental","paternal.bias","paternal.only"),
                      labels=c("M", "MB", "B","PB","P")) +
  labs(title="CF testis", y="Paternal SNP counts", x = "Maternal SNP counts") +
  geom_abline(slope = 1,linetype="dotted", color = "gray50") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))  +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

fig.1c <- grid.arrange(genes.hybrid.soma.counts.plot,genes.soma.scatterplot,genes.hybrid.testis.counts.plot,genes.testis.scatterplot,ncol=2)

###### examine expression

s.gene.tpm = join(genes_final_filter_hybrid_soma,s_tpm_expressed_join,by="gene")
t.gene.tpm = join(genes_final_filter_hybrid_testis,t_tpm_expressed_join,by="gene")

ggplot(s.gene.tpm, aes(bias,log10(s_TPM))) + geom_point() + geom_smooth(method=loess) 
ggplot(t.gene.tpm, aes(bias,log10(t_TPM))) + geom_point() + geom_smooth(method=loess)

s.gene.tpm.stats = ddply(s.gene.tpm, c("bias.cat"), summarise, mean=mean(s_TPM),sd = sd(s_TPM),sem = sd(s_TPM)/sqrt(length(s_TPM)))
t.gene.tpm.stats = ddply(t.gene.tpm, c("bias.cat"), summarise, mean=mean(t_TPM),sd = sd(t_TPM),sem = sd(t_TPM)/sqrt(length(t_TPM)))

nrow(s.gene.tpm)
nrow(t.gene.tpm)

s.gene.tpm.boxplot = ggplot(s.gene.tpm, aes(bias.cat, log10(s_TPM))) +
  geom_jitter(size=0.5,width = 0.25,alpha=1,aes(colour=bias.cat)) +
  geom_boxplot(alpha=0.75,outlier.shape = NA) +
  scale_x_discrete("Expression bias to maternal genome",limits=c("paternal.only","paternal.bias","biparental","maternal.bias","maternal.only"),
                   labels=c(" P ", " PB ", " B "," MB "," M ")) +
  labs(y="TPM (log10)") +
  scale_fill_manual(name="Category of bias", breaks = c("5.paternal.only","4.paternal.bias","3.biparental","2.maternal.bias","1.maternal.only"),
                    labels = c("P", "PB", "B","MB","M"),
                    values = c("#B9770E","#F8C471","#FAD7A0","#F39C12","#D35400"),guide=FALSE) +
  scale_colour_manual(name="Category of bias", breaks = c("5.paternal.only","4.paternal.bias","3.biparental","2.maternal.bias","1.maternal.only"),
                      labels = c("P", "PB", "B","MB","M"),
                      values = c("#B9770E","#F8C471","#FAD7A0","#F39C12","#D35400"),guide=FALSE) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))  +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

t.gene.tpm.boxplot = ggplot(t.gene.tpm, aes(bias.cat, log10(t_TPM))) +
  geom_jitter(size=0.5,width = 0.25,alpha=1,aes(colour=bias.cat)) +
  geom_boxplot(alpha=0.5,outlier.shape = NA) +
  scale_x_discrete("Expression bias to maternal genome",limits=c("paternal.only","paternal.bias","biparental","maternal.bias","maternal.only"),
                   labels=c(" P ", " PB ", " B "," MB "," M ")) +
  labs(y="TPM (log10)") +
  scale_fill_manual(name="Category of bias", breaks = c("5.paternal.only","4.paternal.bias","3.biparental","2.maternal.bias","1.maternal.only"),
                    labels = c("P", "PB", "B","MB","M"),
                    values = c("#B9770E","#F8C471","#FAD7A0","#F39C12","#D35400"),guide=FALSE) +
  scale_colour_manual(name="Category of bias", breaks = c("5.paternal.only","4.paternal.bias","3.biparental","2.maternal.bias","1.maternal.only"),
                      labels = c("P", "PB", "B","MB","M"),
                      values = c("#2980B9","#85C1E9","#AED6F1","#3498DB","#2471A3"),guide=FALSE) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "black",size=0.01))

grid.arrange(s.gene.tpm.boxplot,t.gene.tpm.boxplot,ncol=2)

s.gene.tpm.barplot = ggplot(s.gene.tpm.stats, aes(bias.cat, mean, fill = bias.cat)) +
  geom_bar(stat="identity", colour="black") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete("Expression bias to maternal genome",limits=c("paternal.only","paternal.bias","biparental","maternal.bias","maternal.only"),
                   labels=c(" P ", " PB ", " B "," MB "," M ")) +
  labs(y="Average TPM") +
  scale_fill_manual(name="Category of bias", breaks = c("paternal.only","paternal.bias","biparental","maternal.bias","maternal.only"),
                    labels = c("P", "PB", "B","MB","M"),
                    values = c("#B9770E","#F8C471","#FAD7A0","#F39C12","#D35400"),guide=FALSE) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "black",size=0.01))

t.gene.tpm.barplot = ggplot(t.gene.tpm.stats, aes(bias.cat, mean, fill = bias.cat)) +
  geom_bar(stat="identity", colour="black") +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete("Category of bias",limits=c("paternal.only","paternal.bias","biparental","maternal.bias","maternal.only"),
                   labels=c(" P ", " PB ", " B "," MB "," M ")) +
  labs(y="Average TPM") +
  scale_fill_manual(name="Category of bias", breaks = c("paternal.only","paternal.bias","biparental","maternal.bias","maternal.only"),
                    labels = c(" P ", " PB ", " B "," MB "," M "),
                    values = c("#2980B9","#85C1E9","#AED6F1","#3498DB","#2471A3"),guide=FALSE) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "black",size=0.01))

s.vs.tmp = ggplot(s.gene.tpm, aes(x=bias, y=log10(s_TPM))) + geom_point(size=0.5) + stat_smooth(colour="orange1") +
  labs(y = "TMP (log10)",
       x="Expression bias to maternal genome") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "black",size=0.01))

t.vs.tmp = ggplot(t.gene.tpm, aes(x=bias, y=log10(t_TPM))) + geom_point(size=0.5) + stat_smooth(colour="steelblue4") +
  labs(y = "TMP (log10)",
       x="Expression bias to maternal genome") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "black",size=0.01))

grid.arrange(s.vs.tmp,t.vs.tmp,s.gene.tpm.boxplot,t.gene.tpm.boxplot,nrow=2)

##### examine

# test consistency across replicates

genes_final_filter_hybrid_soma_consistency <- genes_final_filter_hybrid_soma[c(1,2,4,5,7,8,10)]
genes_final_filter_hybrid_testis_consistency <- genes_final_filter_hybrid_testis[c(1,2,4,5,7,8,10)]

genes_final_filter_hybrid_soma_consistency$bias1 <- genes_final_filter_hybrid_soma_consistency$ref.s1 / genes_final_filter_hybrid_soma_consistency$total.s1
genes_final_filter_hybrid_soma_consistency$bias3 <- genes_final_filter_hybrid_soma_consistency$ref.s3 / genes_final_filter_hybrid_soma_consistency$total.s3
genes_final_filter_hybrid_soma_consistency$bias6 <- genes_final_filter_hybrid_soma_consistency$ref.s6 / genes_final_filter_hybrid_soma_consistency$total.s6

genes_final_filter_hybrid_soma_consistency$diff13 <- abs(genes_final_filter_hybrid_soma_consistency$bias1 - genes_final_filter_hybrid_soma_consistency$bias3)
genes_final_filter_hybrid_soma_consistency$diff16 <- abs(genes_final_filter_hybrid_soma_consistency$bias1 - genes_final_filter_hybrid_soma_consistency$bias6)
genes_final_filter_hybrid_soma_consistency$diff36 <- abs(genes_final_filter_hybrid_soma_consistency$bias3 - genes_final_filter_hybrid_soma_consistency$bias6)

genes_final_filter_hybrid_soma_consistency$gene.diff <- round((genes_final_filter_hybrid_soma_consistency$diff13 +
                                                                 genes_final_filter_hybrid_soma_consistency$diff16 +
                                                                 genes_final_filter_hybrid_soma_consistency$diff36)/3,3)

mean(genes_final_filter_hybrid_soma_consistency$gene.diff)

genes_final_filter_hybrid_testis_consistency$bias1 <- genes_final_filter_hybrid_testis_consistency$ref.t1 / genes_final_filter_hybrid_testis_consistency$total.t1
genes_final_filter_hybrid_testis_consistency$bias3 <- genes_final_filter_hybrid_testis_consistency$ref.t3 / genes_final_filter_hybrid_testis_consistency$total.t3
genes_final_filter_hybrid_testis_consistency$bias6 <- genes_final_filter_hybrid_testis_consistency$ref.t6 / genes_final_filter_hybrid_testis_consistency$total.t6

genes_final_filter_hybrid_testis_consistency$diff13 <- abs(genes_final_filter_hybrid_testis_consistency$bias1 - genes_final_filter_hybrid_testis_consistency$bias3)
genes_final_filter_hybrid_testis_consistency$diff16 <- abs(genes_final_filter_hybrid_testis_consistency$bias1 - genes_final_filter_hybrid_testis_consistency$bias6)
genes_final_filter_hybrid_testis_consistency$diff36 <- abs(genes_final_filter_hybrid_testis_consistency$bias3 - genes_final_filter_hybrid_testis_consistency$bias6)

genes_final_filter_hybrid_testis_consistency$gene.diff <- round((genes_final_filter_hybrid_testis_consistency$diff13 +
                                                                 genes_final_filter_hybrid_testis_consistency$diff16 +
                                                                 genes_final_filter_hybrid_testis_consistency$diff36)/3,3)

mean(genes_final_filter_hybrid_testis_consistency$gene.diff)

# SNP density

mean(genes_final_filter_hybrid_soma$N)
mean(genes_final_filter_hybrid_testis$N)

find.mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
}

find.mode(genes_final_filter_hybrid_soma$N)
find.mode(genes_final_filter_hybrid_testis$N)

# correlation between tissues

genes_final_filter_hybrid_soma_corr <- genes_final_filter_hybrid_soma[c(1,15,19)]
genes_final_filter_hybrid_testis_corr <- genes_final_filter_hybrid_testis[c(1,15,19)]
nrow(genes_final_filter_hybrid_testis_corr)
names(genes_final_filter_hybrid_soma_corr) <- c("gene", "s.bias","s.cat")
names(genes_final_filter_hybrid_testis_corr) <- c("gene", "t.bias","t.cat")

tissue_correlation_hybrid <- merge(genes_final_filter_hybrid_soma_corr, genes_final_filter_hybrid_testis_corr, by=c("gene"))
tissue_correlation_hybrid$agree <- ifelse(tissue_correlation_hybrid$s.cat == tissue_correlation_hybrid$t.cat,"y","n")
count(tissue_correlation_hybrid$agree)
cor(tissue_correlation_hybrid$s.bias, tissue_correlation_hybrid$t.bias, method=c("spearman"))

fig.1d = ggplot(tissue_correlation_hybrid, aes(x=s.bias, y=t.bias)) + geom_point(size=0.5,aes(colour=agree)) + 
  labs(title="",x = expression(paste(p[m]," in soma",paste())), y = expression(paste(p[m]," in testis",paste()))) +
  scale_colour_manual(values=c("darkslategray","darkturquoise")) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position = "none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "black",size=0))

tissue_correlation_hybrid_count = ddply(tissue_correlation_hybrid,c("s.cat","t.cat"),summarise,N = length(gene)) # build overlap matrix

tissue_correlation_hybrid_count_s = ddply(tissue_correlation_hybrid,c("s.cat"),summarise,N = length(gene)) # build overlap matrix
tissue_correlation_hybrid_count_t = ddply(tissue_correlation_hybrid,c("t.cat"),summarise,N = length(gene)) # build overlap matrix
new_order <- c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only")
tissue_correlation_hybrid_count_s[match(new_order,tissue_correlation_hybrid_count_s$s.cat),]
tissue_correlation_hybrid_count_t[match(new_order,tissue_correlation_hybrid_count_t$t.cat),]

fig.1e <- ggplot(tissue_correlation_hybrid_count, aes(s.cat, t.cat)) + geom_point(aes(size = 5*N),alpha=0.2,
                                                                        colour=c("darkturquoise","darkslategray","darkslategray",
                                                                               "darkslategray","darkslategray","darkturquoise",
                                                                               "darkslategray","darkslategray","darkturquoise",
                                                                               "darkslategray","darkslategray","darkturquoise","darkturquoise")) +
  geom_text(aes(label = N)) + scale_size(range = c(5, 20)) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  scale_x_discrete(limits = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"), labels = c("M\n(770)", "MB\n(3034)","B\n(195)","PB\n(5)","P\n(1)")) +
  scale_y_discrete(limits = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"), labels = c("M\n(3172)", "MB\n(816)","B\n(13)","PB\n(3)","P\n(1)")) + 
  labs(title="",x = "Soma", y = "Testis") +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0.01),panel.grid.major = element_line(colour = "black",size=0.01))
