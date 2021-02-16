rm(list=ls())
ls()
library(plry)
library(dplyr)
library(tidyr)
library(Rmisc)
library(ggplot2)
library(readr)
library(lattice)
library(grid)
library(gridExtra)
library(reshape2)

### import files

# reciprocal best blast hits

dmela.blast <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/all_methods/blast/DMEL.vs.PCIT.1e-05.25.rbbh.txt",
                          "\t", escape_double = FALSE, trim_ws = TRUE)
aphid.blast <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/all_methods/blast/APISUM.vs.PCITbis.1e-05.25.rbbh.txt",
                         "\t", escape_double = FALSE, trim_ws = TRUE)

# D. mel annotation

dmela.isoforms.to.genes <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/all_methods/blast/DMELA.isoforms.to.genes.txt",
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
dmel.annotation <- read_table2("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/all_methods/blast/DMEL_annotation.tsv", 
                               col_names = FALSE)
# A. pisum annotation

apisum.annotation <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/all_methods/blast/APISUM_annotation.txt", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
# P. citri

pcitri.annotation <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/all_methods/blast/Planococcus_citri_Pcitri.v0.proteins.fa.grouped.by.transcript.tsv", 
    "\t", escape_double = FALSE, trim_ws = TRUE)
pcitri.uniprot.raw <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/all_methods/blast/Planococcus_citri_Pcitri.v0.proteins.fa.blastp.uniprot_sprot.1e-10.tsv", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)

# Genes with ASE

genes.hybrid.soma <- read_csv("~/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_soma.csv")
genes.hybrid.testis <- read_csv("~/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_testis.csv")
genes.intra <- read_csv("~/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes.intra.combined.csv")

# obtain genes of interest

genes.hybrid.soma.of.interest <- genes.hybrid.soma[genes.hybrid.soma$bias.cat == "biparental" |
                                                     genes.hybrid.soma$bias.cat == "paternal.bias" |
                                                     genes.hybrid.soma$bias.cat == "paternal.only",][c("gene","bias","bias.cat")]
genes.hybrid.soma.of.interest <- genes.hybrid.soma.of.interest[order(genes.hybrid.soma.of.interest$bias),]

genes.hybrid.testis.of.interest <- genes.hybrid.testis[genes.hybrid.testis$bias.cat == "biparental" |
                                                     genes.hybrid.testis$bias.cat == "paternal.bias" |
                                                     genes.hybrid.testis$bias.cat == "paternal.only",][c("gene","bias","bias.cat")]
genes.hybrid.testis.of.interest <- genes.hybrid.testis.of.interest[order(genes.hybrid.testis.of.interest$bias),]

genes.intra.of.interest <- genes.intra[genes.intra$bias.cat.intra == "biparental" |
                                                         genes.intra$bias.cat.intra == "paternal.bias" |
                                                         genes.intra$bias.cat.intra == "paternal.only",][c("gene","bias.cat.WxC","bias.cat.WxB","bias.cat.CxB","bias.intra","bias.cat.intra")]
genes.intra.of.interest <- genes.intra.of.interest[order(genes.intra.of.interest$bias.intra),]

# do they have annotation features?

colnames(pcitri.annotation)[1] <- "gene"
pcitri.annotation$gene = sub('\\..*', '', pcitri.annotation$gene)

genes.hybrid.soma.of.interest$annotated <- ifelse(genes.hybrid.soma.of.interest$gene %in% pcitri.annotation$gene == TRUE,"Yes","No")
genes.hybrid.testis.of.interest$annotated <- ifelse(genes.hybrid.testis.of.interest$gene %in% pcitri.annotation$gene == TRUE,"Yes","No")
genes.intra.of.interest$annotated <- ifelse(genes.intra.of.interest$gene %in% pcitri.annotation$gene == TRUE,"Yes","No")
count(genes.intra.of.interest$annotated)

# process D. mel files

colnames(dmela.blast)[1] <- "dmela.isoform"
colnames(dmela.blast)[2] <- "gene"
colnames(dmela.blast)[3] <- "length.diff.to.dmel"
dmela.blast$gene = sub('\\..*', '', dmela.blast$gene)
dmela.blast<-dmela.blast[c(1,2,3)]

long.format.dmela.isoforms.to.genes.0 <- melt(dmela.isoforms.to.genes, c("dmela.mRNA"))
long.format.dmela.isoforms.to.genes <- long.format.dmela.isoforms.to.genes.0[c(1,3)]
colnames(long.format.dmela.isoforms.to.genes)[2] <- "dmela.isoform"
dmela.gene.prot <- long.format.dmela.isoforms.to.genes[complete.cases(long.format.dmela.isoforms.to.genes),]
dmela.blast <- join(dmela.blast, dmela.gene.prot, by="dmela.isoform")
dmela.blast <- dmela.blast[c("gene","dmela.mRNA")]

colnames(dmel.annotation)[1] <- "dmela.gene"
colnames(dmel.annotation)[2] <- "dmela.mRNA"
dmel.annotation <- dmel.annotation[c(1,2)]
dmela.blast.all <- join(dmela.blast,dmel.annotation,by="dmela.mRNA")
dmela.blast.all <- dmela.blast.all[c("gene","dmela.gene")]

# assign ortholog in D. mel

genes.hybrid.soma.of.interest.dmel <- left_join(genes.hybrid.soma.of.interest,dmela.blast.all,by="gene")
genes.hybrid.testis.of.interest.dmel <- left_join(genes.hybrid.testis.of.interest,dmela.blast.all,by="gene")
genes.intra.of.interest.dmel <- left_join(genes.intra.of.interest,dmela.blast.all,by="gene")

# process A. pisum files

colnames(aphid.blast)[1] <- "apisum.gene"
colnames(aphid.blast)[2] <- "gene"
colnames(aphid.blast)[3] <- "length.diff.to.apisum"
aphid.blast$gene = sub('\\..*', '', aphid.blast$gene)
apisum.blast<-aphid.blast[c(1,2,3)]

colnames(apisum.annotation)[1] <- "apisum.gene"
colnames(apisum.annotation)[3] <- "apisum.anno"
apisum.annotation <- apisum.annotation[c(1,3)]
apisum.annotation.unique <- apisum.annotation[complete.cases(apisum.annotation),]

apisum.blast.all <- join(apisum.blast, apisum.annotation.unique,by="apisum.gene")
apisum.blast.all <- apisum.blast.all[c("gene","apisum.gene","apisum.anno")]
apisum.blast.all.collapsed <- ddply(apisum.blast.all,c("gene"),summarize,
                                     apisum.gene = paste(apisum.gene ,collapse=","), 
                                     apisum.anno = paste(apisum.anno,collapse=","))

# assign ortholog in A. pisum

genes.hybrid.soma.of.interest.both <- left_join(genes.hybrid.soma.of.interest.dmel, apisum.blast.all.collapsed ,by="gene")
genes.hybrid.testis.of.interest.both <- left_join(genes.hybrid.testis.of.interest.dmel, apisum.blast.all.collapsed ,by="gene")
genes.intra.of.interest.both <- left_join(genes.intra.of.interest.dmel, apisum.blast.all.collapsed ,by="gene")

# export

#write.csv(genes.hybrid.soma.of.interest.both, file = "/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_annotated/genes.hybrid.soma.of.interest.both.csv")
#write.csv(genes.hybrid.testis.of.interest.both, file = "/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_annotated/genes.hybrid.testis.of.interest.both.csv")
#write.csv(genes.intra.of.interest.both, file = "/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_annotated/genes.intra.of.interest.both.csv")

# what does Swiss Uniprot say?

colnames(pcitri.uniprot.raw)[1] <- "gene"
colnames(pcitri.uniprot.raw)[11] <- "E"
colnames(pcitri.uniprot.raw)[15] <- "swiss.prot.annotation"

# that's a long list! use duplicated() to retain first (best) match

pcitri.uniprot.best.hit <- pcitri.uniprot.raw[!duplicated(pcitri.uniprot.raw$gene),]
pcitri.uniprot.best.hit$gene <- sub('\\..*', '', pcitri.uniprot.best.hit$gene)

pcitri.uniprot <- pcitri.uniprot.best.hit[c(1,11,15)]
pcitri.uniprot.0 <- separate(pcitri.uniprot, swiss.prot.annotation, sep=" OS=", c("anno", "taxon"))
pcitri.uniprot.00 <- separate(pcitri.uniprot.0, taxon,sep=" PE=", c("taxon", "else"))[c(1,2,3,4)]
pcitri.uniprot.000 <- separate(pcitri.uniprot.00, taxon,sep=" GN=", c("taxon", "matched_gene"))[c(1,2,3,4,5)]
pcitri.uniprot.0000 <- pcitri.uniprot.000[c(1,2,4,5,3)]
pcitri.uniprot.00000 <- separate(pcitri.uniprot.0000, anno,sep="\\s", c("discard", "uniprot.anno"),extra = "merge") # to separate by first space

pcitri.uniprot.000000 <- pcitri.uniprot.00000[!duplicated(pcitri.uniprot.00000$gene),]
pcitri.uniprot.merge <- pcitri.uniprot.000000[c(1,6,4,3,2)]
count(pcitri.uniprot.merge$taxon) # do best matches make sense? check taxon counts

# merge with ASE genes
nrow(genes.hybrid.soma.of.interest.complete.anno)
genes.hybrid.soma.of.interest.complete.anno <- left_join(genes.hybrid.soma.of.interest.both, pcitri.uniprot.merge, by="gene")
genes.hybrid.testis.of.interest.complete.anno <- left_join(genes.hybrid.testis.of.interest.both, pcitri.uniprot.merge, by="gene")
genes.intra.of.interest.complete.anno <- left_join(genes.intra.of.interest.both, pcitri.uniprot.merge, by="gene")

# extract bias only to merge
genes.hybrid.soma.of.interest.bias <- genes.hybrid.soma.of.interest.complete.anno[c("gene","bias.cat")]
colnames(genes.hybrid.soma.of.interest.bias)[2] <- "bias.cat.soma"
genes.hybrid.testis.of.interest.bias <- genes.hybrid.testis.of.interest.complete.anno[c("gene","bias.cat")]
colnames(genes.hybrid.testis.of.interest.bias)[2] <- "bias.cat.testis"
genes.intra.of.interest.bias <- genes.intra.of.interest.complete.anno[c("gene","bias.cat.intra")]

# merge
genes.intra.of.interest.complete.anno.0 <- left_join(genes.intra.of.interest.complete.anno, genes.hybrid.soma.of.interest.bias,by="gene")
genes.intra.complete.anno.export <- left_join(genes.intra.of.interest.complete.anno.0, genes.hybrid.testis.of.interest.bias,by="gene")

new.anno.soma <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_annotated/genes.hybrid.soma.full.annotation.csv",
                          ",", escape_double = FALSE, trim_ws = TRUE)
new.anno.soma <- new.anno.soma[c("gene","uniprot.anno")]

genes.intra.complete.anno.export.change.anno <- left_join(genes.intra.complete.anno.export, new.anno.soma, by=c("gene"))
genes.intra.complete.anno.export.change.anno$uniprot.anno <- ifelse(is.na(genes.intra.complete.anno.export.change.anno$uniprot.anno.y),
                                                                    genes.intra.complete.anno.export.change.anno$uniprot.anno.x,genes.intra.complete.anno.export.change.anno$uniprot.anno.y)

genes.intra.complete.anno.export <- genes.intra.complete.anno.export

#write.csv(genes.hybrid.soma.of.interest.complete.anno, file = "/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_annotated/genes.hybrid.soma.complete.anno.export.csv")
#write.csv(genes.hybrid.testis.of.interest.complete.anno, file = "/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_annotated/genes.hybrid.testis.complete.anno.export.csv")
write.csv(genes.intra.complete.anno.export, file = "/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/genes_annotated/genes.intra.complete.anno.export.csv")
