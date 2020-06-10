rm(list=ls())
ls()
library("dplyr")
library("tidyr")
library("Rmisc")
library("ggplot2")
library("readr")
library("lattice")
library("grid")
library("gridExtra")

###################################

##### Import and organise raw data

# import raw output from ASEReadCounter

WC01.cw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC01.cw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WC01.wc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC01.wc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WC04.cw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC04.cw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WC04.wc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC04.wc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WC06.cw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC06.cw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WC06.wc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC06.wc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WC19.cw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC19.cw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WC19.wc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC19.wc.csv","\t", escape_double = FALSE, trim_ws = TRUE)

CW01.cw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW01.cw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CW01.wc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW01.wc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CW04.cw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW04.cw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CW04.wc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW04.wc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CW10.cw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW10.cw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CW10.wc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW10.wc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CW11.cw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW11.cw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CW11.wc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW11.wc.csv","\t", escape_double = FALSE, trim_ws = TRUE)

WB02.bw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB02.bw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WB02.wb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB02.wb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WB03.bw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB03.bw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WB03.wb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB03.wb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WB08.bw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB08.bw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
WB08.wb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB08.wb.csv","\t", escape_double = FALSE, trim_ws = TRUE)

BW02.bw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW02.bw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BW02.wb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW02.wb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BW03.bw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW03.bw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BW03.wb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW03.wb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BW05.bw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW05.bw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BW05.wb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW05.wb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BW12.bw.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW12.bw.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BW12.wb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW12.wb.csv","\t", escape_double = FALSE, trim_ws = TRUE)

CB01.bc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB01.bc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CB01.cb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB01.cb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CB04.bc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB04.bc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CB04.cb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB04.cb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CB12.bc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB12.bc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CB12.cb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB12.cb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CB20.bc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB20.bc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
CB20.cb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB20.cb.csv","\t", escape_double = FALSE, trim_ws = TRUE)

BC03.bc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC03.bc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BC03.cb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC03.cb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BC04.bc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC04.bc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BC04.cb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC04.cb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BC13.bc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC13.bc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BC13.cb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC13.cb.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BC17.bc.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC17.bc.csv","\t", escape_double = FALSE, trim_ws = TRUE)
BC17.cb.raw <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC17.cb.csv","\t", escape_double = FALSE, trim_ws = TRUE)

# import expression data from RSEM

WC01_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WC1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WC04_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WC4.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WC06_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WC6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WC19_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WC19.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CW01_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CW1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CW04_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CW4.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CW10_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CW10.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CW11_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CW11.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WB02_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WB2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WB03_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WB3.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WB08_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WB8.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BW02_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BW2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BW03_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BW3.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BW05_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BW5.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BW12_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BW12.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CB01_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CB1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CB04_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CB4.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CB12_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CB12.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CB20_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CB20.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BC03_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BC3.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BC04_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BC4.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BC13_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BC13.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BC17_genes_results <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BC17.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)

# explore raw output

# W and C

WC01.cw.raw.hist = ggplot(WC01.cw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="WC01.cw",y=nrow(WC01.cw.raw)) 
WC01.wc.raw.hist = ggplot(WC01.wc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="WC01.wc",y=nrow(WC01.wc.raw)) 
WC04.cw.raw.hist = ggplot(WC04.cw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="WC04.cw",y=nrow(WC04.cw.raw)) 
WC04.wc.raw.hist = ggplot(WC04.wc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="WC04.wc",y=nrow(WC04.wc.raw)) 
WC06.cw.raw.hist = ggplot(WC06.cw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="WC06.cw",y=nrow(WC06.cw.raw)) 
WC06.wc.raw.hist = ggplot(WC06.wc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="WC06.wc",y=nrow(WC06.wc.raw)) 
WC19.cw.raw.hist = ggplot(WC19.cw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="WC19.cw",y=nrow(WC19.cw.raw)) 
WC19.wc.raw.hist = ggplot(WC19.wc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="WC19.wc",y=nrow(WC19.wc.raw))

CW01.cw.raw.hist = ggplot(CW01.cw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="CW01.cw",y=nrow(CW01.cw.raw)) 
CW01.wc.raw.hist = ggplot(CW01.wc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="CW01.wc",y=nrow(CW01.wc.raw)) 
CW04.cw.raw.hist = ggplot(CW04.cw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="CW04.cw",y=nrow(CW04.cw.raw)) 
CW04.wc.raw.hist = ggplot(CW04.wc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="CW04.wc",y=nrow(CW04.wc.raw)) 
CW10.cw.raw.hist = ggplot(CW10.cw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="CW10.cw",y=nrow(CW10.cw.raw)) 
CW10.wc.raw.hist = ggplot(CW10.wc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="CW10.wc",y=nrow(CW10.wc.raw)) 
CW11.cw.raw.hist = ggplot(CW11.cw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="CW11.cw",y=nrow(CW11.cw.raw)) 
CW11.wc.raw.hist = ggplot(CW11.wc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="CW11.wc",y=nrow(CW11.wc.raw))

#grid.arrange(WC01.cw.raw.hist,WC01.wc.raw.hist,WC04.cw.raw.hist,WC04.wc.raw.hist,WC06.cw.raw.hist,WC06.wc.raw.hist,WC19.cw.raw.hist,WC19.wc.raw.hist,
    #         CW01.cw.raw.hist,CW01.wc.raw.hist,CW04.cw.raw.hist,CW04.wc.raw.hist,CW10.cw.raw.hist,CW10.wc.raw.hist,CW11.cw.raw.hist,CW11.wc.raw.hist,nrow=2)

# W and B

WB02.bw.raw.hist = ggplot(WB02.bw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="WB02.bw",y=nrow(WB02.bw.raw)) 
WB02.wb.raw.hist = ggplot(WB02.wb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="WB02.wb",y=nrow(WB02.wb.raw)) 
WB03.bw.raw.hist = ggplot(WB03.bw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="WB03.bw",y=nrow(WB03.bw.raw)) 
WB03.wb.raw.hist = ggplot(WB03.wb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="WB03.wb",y=nrow(WB03.wb.raw)) 
WB08.bw.raw.hist = ggplot(WB08.bw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="WB08.bw",y=nrow(WB08.bw.raw)) 
WB08.wb.raw.hist = ggplot(WB08.wb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="WB08.wb",y=nrow(WB08.wb.raw))

BW02.bw.raw.hist = ggplot(BW02.bw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="BW02.bw",y=nrow(BW02.bw.raw)) 
BW02.wb.raw.hist = ggplot(BW02.wb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="BW02.wb",y=nrow(BW02.wb.raw)) 
BW03.bw.raw.hist = ggplot(BW03.bw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="BW03.bw",y=nrow(BW03.bw.raw)) 
BW03.wb.raw.hist = ggplot(BW03.wb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="BW03.wb",y=nrow(BW03.wb.raw)) 
BW05.bw.raw.hist = ggplot(BW05.bw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="BW05.bw",y=nrow(BW05.bw.raw)) 
BW05.wb.raw.hist = ggplot(BW05.wb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="BW05.wb",y=nrow(BW05.wb.raw)) 
BW12.bw.raw.hist = ggplot(BW12.bw.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="BW12.bw",y=nrow(BW12.bw.raw)) 
BW12.wb.raw.hist = ggplot(BW12.wb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="BW12.wb",y=nrow(BW12.wb.raw))

blank = grid.rect(gp=gpar(col="white"))

#grid.arrange(WB02.bw.raw.hist,WB02.wb.raw.hist,WB03.bw.raw.hist,WB03.wb.raw.hist,WB08.bw.raw.hist,WB08.wb.raw.hist,blank,blank,
    #         BW02.bw.raw.hist,BW02.wb.raw.hist,BW03.bw.raw.hist,BW03.wb.raw.hist,BW05.bw.raw.hist,BW05.wb.raw.hist,BW12.bw.raw.hist,BW12.wb.raw.hist,nrow=2)

# C and B

CB01.bc.raw.hist = ggplot(CB01.bc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="CB01.bc",y=nrow(CB01.bc.raw)) 
CB01.cb.raw.hist = ggplot(CB01.cb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="CB01.cb",y=nrow(CB01.cb.raw)) 
CB04.bc.raw.hist = ggplot(CB04.bc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="CB04.bc",y=nrow(CB04.bc.raw)) 
CB04.cb.raw.hist = ggplot(CB04.cb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="CB04.cb",y=nrow(CB04.cb.raw)) 
CB12.bc.raw.hist = ggplot(CB12.bc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="CB12.bc",y=nrow(CB12.bc.raw)) 
CB12.cb.raw.hist = ggplot(CB12.cb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="CB12.cb",y=nrow(CB12.cb.raw)) 
CB20.bc.raw.hist = ggplot(CB20.bc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="CB20.bc",y=nrow(CB20.bc.raw)) 
CB20.cb.raw.hist = ggplot(CB20.cb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="CB20.cb",y=nrow(CB20.cb.raw))

BC03.bc.raw.hist = ggplot(BC03.bc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="BC03.bc",y=nrow(BC03.bc.raw)) 
BC03.cb.raw.hist = ggplot(BC03.cb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="BC03.cb",y=nrow(BC03.cb.raw)) 
BC04.bc.raw.hist = ggplot(BC04.bc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="BC04.bc",y=nrow(BC04.bc.raw)) 
BC04.cb.raw.hist = ggplot(BC04.cb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="BC04.cb",y=nrow(BC04.cb.raw)) 
BC13.bc.raw.hist = ggplot(BC13.bc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="BC13.bc",y=nrow(BC13.bc.raw)) 
BC13.cb.raw.hist = ggplot(BC13.cb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="BC13.cb",y=nrow(BC13.cb.raw)) 
BC17.bc.raw.hist = ggplot(BC17.bc.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="blue") + labs(title="BC17.bc",y=nrow(BC17.bc.raw)) 
BC17.cb.raw.hist = ggplot(BC17.cb.raw, aes(refCount/totalCount)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="red1") + labs(title="BC17.cb",y=nrow(BC17.cb.raw))

#grid.arrange(CB01.bc.raw.hist,CB01.cb.raw.hist,CB04.bc.raw.hist,CB04.cb.raw.hist,CB12.bc.raw.hist,CB12.cb.raw.hist,CB20.bc.raw.hist,CB20.cb.raw.hist,
        #     BC03.bc.raw.hist,BC03.cb.raw.hist,BC04.bc.raw.hist,BC04.cb.raw.hist,BC13.bc.raw.hist,BC13.cb.raw.hist,BC17.bc.raw.hist,BC17.cb.raw.hist,nrow=2)

#add an origin (which vcf file did the variant originate from?) and variant column

WC01.cw.edit<-WC01.cw.raw
WC01.wc.edit<-WC01.wc.raw
WC04.cw.edit<-WC04.cw.raw
WC04.wc.edit<-WC04.wc.raw
WC06.cw.edit<-WC06.cw.raw
WC06.wc.edit<-WC06.wc.raw
WC19.cw.edit<-WC19.cw.raw
WC19.wc.edit<-WC19.wc.raw
CW01.cw.edit<-CW01.cw.raw
CW01.wc.edit<-CW01.wc.raw
CW04.cw.edit<-CW04.cw.raw
CW04.wc.edit<-CW04.wc.raw
CW10.cw.edit<-CW10.cw.raw
CW10.wc.edit<-CW10.wc.raw
CW11.cw.edit<-CW11.cw.raw
CW11.wc.edit<-CW11.wc.raw
WB02.bw.edit<-WB02.bw.raw
WB02.wb.edit<-WB02.wb.raw
WB03.bw.edit<-WB03.bw.raw
WB03.wb.edit<-WB03.wb.raw
WB08.bw.edit<-WB08.bw.raw
WB08.wb.edit<-WB08.wb.raw
BW02.bw.edit<-BW02.bw.raw
BW02.wb.edit<-BW02.wb.raw
BW03.bw.edit<-BW03.bw.raw
BW03.wb.edit<-BW03.wb.raw
BW05.bw.edit<-BW05.bw.raw
BW05.wb.edit<-BW05.wb.raw
BW12.bw.edit<-BW12.bw.raw
BW12.wb.edit<-BW12.wb.raw
CB01.bc.edit<-CB01.bc.raw
CB01.cb.edit<-CB01.cb.raw
CB04.bc.edit<-CB04.bc.raw
CB04.cb.edit<-CB04.cb.raw
CB12.bc.edit<-CB12.bc.raw
CB12.cb.edit<-CB12.cb.raw
CB20.bc.edit<-CB20.bc.raw
CB20.cb.edit<-CB20.cb.raw
BC03.bc.edit<-BC03.bc.raw
BC03.cb.edit<-BC03.cb.raw
BC04.bc.edit<-BC04.bc.raw
BC04.cb.edit<-BC04.cb.raw
BC13.bc.edit<-BC13.bc.raw
BC13.cb.edit<-BC13.cb.raw
BC17.bc.edit<-BC17.bc.raw
BC17.cb.edit<-BC17.cb.raw

WC01.cw.edit$origin = "cw"
WC01.wc.edit$origin = "wc"
WC04.cw.edit$origin = "cw"
WC04.wc.edit$origin = "wc"
WC06.cw.edit$origin = "cw"
WC06.wc.edit$origin = "wc"
WC19.cw.edit$origin = "cw"
WC19.wc.edit$origin = "wc"
CW01.cw.edit$origin = "cw"
CW01.wc.edit$origin = "wc"
CW04.cw.edit$origin = "cw"
CW04.wc.edit$origin = "wc"
CW10.cw.edit$origin = "cw"
CW10.wc.edit$origin = "wc"
CW11.cw.edit$origin = "cw"
CW11.wc.edit$origin = "wc"
WB02.bw.edit$origin = "bw"
WB02.wb.edit$origin = "wb"
WB03.bw.edit$origin = "bw"
WB03.wb.edit$origin = "wb"
WB08.bw.edit$origin = "bw"
WB08.wb.edit$origin = "wb"
BW02.bw.edit$origin = "bw"
BW02.wb.edit$origin = "wb"
BW03.bw.edit$origin = "bw"
BW03.wb.edit$origin = "wb"
BW05.bw.edit$origin = "bw"
BW05.wb.edit$origin = "wb"
BW12.bw.edit$origin = "bw"
BW12.wb.edit$origin = "wb"
CB01.bc.edit$origin = "bc"
CB01.cb.edit$origin = "cb"
CB04.bc.edit$origin = "bc"
CB04.cb.edit$origin = "cb"
CB12.bc.edit$origin = "bc"
CB12.cb.edit$origin = "cb"
CB20.bc.edit$origin = "bc"
CB20.cb.edit$origin = "cb"
BC03.bc.edit$origin = "bc"
BC03.cb.edit$origin = "cb"
BC04.bc.edit$origin = "bc"
BC04.cb.edit$origin = "cb"
BC13.bc.edit$origin = "bc"
BC13.cb.edit$origin = "cb"
BC17.bc.edit$origin = "bc"
BC17.cb.edit$origin = "cb"

WC01.cw.edit$variant=paste(WC01.cw.edit$contig,WC01.cw.edit$position,sep=":")
WC01.wc.edit$variant=paste(WC01.wc.edit$contig,WC01.wc.edit$position,sep=":")
WC04.cw.edit$variant=paste(WC04.cw.edit$contig,WC04.cw.edit$position,sep=":")
WC04.wc.edit$variant=paste(WC04.wc.edit$contig,WC04.wc.edit$position,sep=":")
WC06.cw.edit$variant=paste(WC06.cw.edit$contig,WC06.cw.edit$position,sep=":")
WC06.wc.edit$variant=paste(WC06.wc.edit$contig,WC06.wc.edit$position,sep=":")
WC19.cw.edit$variant=paste(WC19.cw.edit$contig,WC19.cw.edit$position,sep=":")
WC19.wc.edit$variant=paste(WC19.wc.edit$contig,WC19.wc.edit$position,sep=":")
CW01.cw.edit$variant=paste(CW01.cw.edit$contig,CW01.cw.edit$position,sep=":")
CW01.wc.edit$variant=paste(CW01.wc.edit$contig,CW01.wc.edit$position,sep=":")
CW04.cw.edit$variant=paste(CW04.cw.edit$contig,CW04.cw.edit$position,sep=":")
CW04.wc.edit$variant=paste(CW04.wc.edit$contig,CW04.wc.edit$position,sep=":")
CW10.cw.edit$variant=paste(CW10.cw.edit$contig,CW10.cw.edit$position,sep=":")
CW10.wc.edit$variant=paste(CW10.wc.edit$contig,CW10.wc.edit$position,sep=":")
CW11.cw.edit$variant=paste(CW11.cw.edit$contig,CW11.cw.edit$position,sep=":")
CW11.wc.edit$variant=paste(CW11.wc.edit$contig,CW11.wc.edit$position,sep=":")
WB02.bw.edit$variant=paste(WB02.bw.edit$contig,WB02.bw.edit$position,sep=":")
WB02.wb.edit$variant=paste(WB02.wb.edit$contig,WB02.wb.edit$position,sep=":")
WB03.bw.edit$variant=paste(WB03.bw.edit$contig,WB03.bw.edit$position,sep=":")
WB03.wb.edit$variant=paste(WB03.wb.edit$contig,WB03.wb.edit$position,sep=":")
WB08.bw.edit$variant=paste(WB08.bw.edit$contig,WB08.bw.edit$position,sep=":")
WB08.wb.edit$variant=paste(WB08.wb.edit$contig,WB08.wb.edit$position,sep=":")
BW02.bw.edit$variant=paste(BW02.bw.edit$contig,BW02.bw.edit$position,sep=":")
BW02.wb.edit$variant=paste(BW02.wb.edit$contig,BW02.wb.edit$position,sep=":")
BW03.bw.edit$variant=paste(BW03.bw.edit$contig,BW03.bw.edit$position,sep=":")
BW03.wb.edit$variant=paste(BW03.wb.edit$contig,BW03.wb.edit$position,sep=":")
BW05.bw.edit$variant=paste(BW05.bw.edit$contig,BW05.bw.edit$position,sep=":")
BW05.wb.edit$variant=paste(BW05.wb.edit$contig,BW05.wb.edit$position,sep=":")
BW12.bw.edit$variant=paste(BW12.bw.edit$contig,BW12.bw.edit$position,sep=":")
BW12.wb.edit$variant=paste(BW12.wb.edit$contig,BW12.wb.edit$position,sep=":")
CB01.bc.edit$variant=paste(CB01.bc.edit$contig,CB01.bc.edit$position,sep=":")
CB01.cb.edit$variant=paste(CB01.cb.edit$contig,CB01.cb.edit$position,sep=":")
CB04.bc.edit$variant=paste(CB04.bc.edit$contig,CB04.bc.edit$position,sep=":")
CB04.cb.edit$variant=paste(CB04.cb.edit$contig,CB04.cb.edit$position,sep=":")
CB12.bc.edit$variant=paste(CB12.bc.edit$contig,CB12.bc.edit$position,sep=":")
CB12.cb.edit$variant=paste(CB12.cb.edit$contig,CB12.cb.edit$position,sep=":")
CB20.bc.edit$variant=paste(CB20.bc.edit$contig,CB20.bc.edit$position,sep=":")
CB20.cb.edit$variant=paste(CB20.cb.edit$contig,CB20.cb.edit$position,sep=":")
BC03.bc.edit$variant=paste(BC03.bc.edit$contig,BC03.bc.edit$position,sep=":")
BC03.cb.edit$variant=paste(BC03.cb.edit$contig,BC03.cb.edit$position,sep=":")
BC04.bc.edit$variant=paste(BC04.bc.edit$contig,BC04.bc.edit$position,sep=":")
BC04.cb.edit$variant=paste(BC04.cb.edit$contig,BC04.cb.edit$position,sep=":")
BC13.bc.edit$variant=paste(BC13.bc.edit$contig,BC13.bc.edit$position,sep=":")
BC13.cb.edit$variant=paste(BC13.cb.edit$contig,BC13.cb.edit$position,sep=":")
BC17.bc.edit$variant=paste(BC17.bc.edit$contig,BC17.bc.edit$position,sep=":")
BC17.cb.edit$variant=paste(BC17.cb.edit$contig,BC17.cb.edit$position,sep=":")

# assign maternal and paternal counts for each df depending on genotype/origin
  
names(WC01.cw.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WC01.wc.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WC04.cw.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WC04.wc.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WC06.cw.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WC06.wc.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WC19.cw.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WC19.wc.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CW01.cw.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CW01.wc.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CW04.cw.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CW04.wc.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CW10.cw.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CW10.wc.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CW11.cw.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CW11.wc.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WB02.bw.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WB02.wb.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WB03.bw.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WB03.wb.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WB08.bw.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(WB08.wb.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BW02.bw.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BW02.wb.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BW03.bw.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BW03.wb.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BW05.bw.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BW05.wb.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BW12.bw.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BW12.wb.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CB01.bc.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CB01.cb.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CB04.bc.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CB04.cb.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CB12.bc.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CB12.cb.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CB20.bc.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(CB20.cb.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BC03.bc.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BC03.cb.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BC04.bc.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BC04.cb.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BC13.bc.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BC13.cb.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BC17.bc.edit) <- c("contig","position","variantID","refAllele","altAllele","pat.count","mat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")
names(BC17.cb.edit) <- c("contig","position","variantID","refAllele","altAllele","mat.count","pat.count","total.count","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","origin","variant")

# join both df for each samples and reorder

WC01 <- rbind(WC01.cw.edit, WC01.wc.edit)
WC04 <- rbind(WC04.cw.edit, WC04.wc.edit)
WC06 <- rbind(WC06.cw.edit, WC06.wc.edit)
WC19 <- rbind(WC19.cw.edit, WC19.wc.edit)
CW01 <- rbind(CW01.cw.edit, CW01.wc.edit)
CW04 <- rbind(CW04.cw.edit, CW04.wc.edit)
CW10 <- rbind(CW10.cw.edit, CW10.wc.edit)
CW11 <- rbind(CW11.cw.edit, CW11.wc.edit)
WB02 <- rbind(WB02.bw.edit, WB02.wb.edit)
WB03 <- rbind(WB03.bw.edit, WB03.wb.edit)
WB08 <- rbind(WB08.bw.edit, WB08.wb.edit)
BW02 <- rbind(BW02.bw.edit, BW02.wb.edit)
BW03 <- rbind(BW03.bw.edit, BW03.wb.edit)
BW05 <- rbind(BW05.bw.edit, BW05.wb.edit)
BW12 <- rbind(BW12.bw.edit, BW12.wb.edit)
CB01 <- rbind(CB01.bc.edit, CB01.cb.edit)
CB04 <- rbind(CB04.bc.edit, CB04.cb.edit)
CB12 <- rbind(CB12.bc.edit, CB12.cb.edit)
CB20 <- rbind(CB20.bc.edit, CB20.cb.edit)
BC03 <- rbind(BC03.bc.edit, BC03.cb.edit)
BC04 <- rbind(BC04.bc.edit, BC04.cb.edit)
BC13 <- rbind(BC13.bc.edit, BC13.cb.edit)
BC17 <- rbind(BC17.bc.edit, BC17.cb.edit)

WC01$order1 <- as.integer(gsub(".*_", '',WC01$contig))
WC01$order2 <- as.integer(gsub(".*:", '',WC01$variant))
WC04$order1 <- as.integer(gsub(".*_", '',WC04$contig))
WC04$order2 <- as.integer(gsub(".*:", '',WC04$variant))
WC06$order1 <- as.integer(gsub(".*_", '',WC06$contig))
WC06$order2 <- as.integer(gsub(".*:", '',WC06$variant))
WC19$order1 <- as.integer(gsub(".*_", '',WC19$contig))
WC19$order2 <- as.integer(gsub(".*:", '',WC19$variant))
CW01$order1 <- as.integer(gsub(".*_", '',CW01$contig))
CW01$order2 <- as.integer(gsub(".*:", '',CW01$variant))
CW04$order1 <- as.integer(gsub(".*_", '',CW04$contig))
CW04$order2 <- as.integer(gsub(".*:", '',CW04$variant))
CW10$order1 <- as.integer(gsub(".*_", '',CW10$contig))
CW10$order2 <- as.integer(gsub(".*:", '',CW10$variant))
CW11$order1 <- as.integer(gsub(".*_", '',CW11$contig))
CW11$order2 <- as.integer(gsub(".*:", '',CW11$variant))
WB02$order1 <- as.integer(gsub(".*_", '',WB02$contig))
WB02$order2 <- as.integer(gsub(".*:", '',WB02$variant))
WB03$order1 <- as.integer(gsub(".*_", '',WB03$contig))
WB03$order2 <- as.integer(gsub(".*:", '',WB03$variant))
WB08$order1 <- as.integer(gsub(".*_", '',WB08$contig))
WB08$order2 <- as.integer(gsub(".*:", '',WB08$variant))
BW02$order1 <- as.integer(gsub(".*_", '',BW02$contig))
BW02$order2 <- as.integer(gsub(".*:", '',BW02$variant))
BW03$order1 <- as.integer(gsub(".*_", '',BW03$contig))
BW03$order2 <- as.integer(gsub(".*:", '',BW03$variant))
BW05$order1 <- as.integer(gsub(".*_", '',BW05$contig))
BW05$order2 <- as.integer(gsub(".*:", '',BW05$variant))
BW12$order1 <- as.integer(gsub(".*_", '',BW12$contig))
BW12$order2 <- as.integer(gsub(".*:", '',BW12$variant))
CB01$order1 <- as.integer(gsub(".*_", '',CB01$contig))
CB01$order2 <- as.integer(gsub(".*:", '',CB01$variant))
CB04$order1 <- as.integer(gsub(".*_", '',CB04$contig))
CB04$order2 <- as.integer(gsub(".*:", '',CB04$variant))
CB12$order1 <- as.integer(gsub(".*_", '',CB12$contig))
CB12$order2 <- as.integer(gsub(".*:", '',CB12$variant))
CB20$order1 <- as.integer(gsub(".*_", '',CB20$contig))
CB20$order2 <- as.integer(gsub(".*:", '',CB20$variant))
BC03$order1 <- as.integer(gsub(".*_", '',BC03$contig))
BC03$order2 <- as.integer(gsub(".*:", '',BC03$variant))
BC04$order1 <- as.integer(gsub(".*_", '',BC04$contig))
BC04$order2 <- as.integer(gsub(".*:", '',BC04$variant))
BC13$order1 <- as.integer(gsub(".*_", '',BC13$contig))
BC13$order2 <- as.integer(gsub(".*:", '',BC13$variant))
BC17$order1 <- as.integer(gsub(".*_", '',BC17$contig))
BC17$order2 <- as.integer(gsub(".*:", '',BC17$variant))

WC01 <- WC01[order(WC01$order1,WC01$order2),]
WC04 <- WC04[order(WC04$order1,WC04$order2),]
WC06 <- WC06[order(WC06$order1,WC06$order2),]
WC19 <- WC19[order(WC19$order1,WC19$order2),]
CW01 <- CW01[order(CW01$order1,CW01$order2),]
CW04 <- CW04[order(CW04$order1,CW04$order2),]
CW10 <- CW10[order(CW10$order1,CW10$order2),]
CW11 <- CW11[order(CW11$order1,CW11$order2),]
WB02 <- WB02[order(WB02$order1,WB02$order2),]
WB03 <- WB03[order(WB03$order1,WB03$order2),]
WB08 <- WB08[order(WB08$order1,WB08$order2),]
BW02 <- BW02[order(BW02$order1,BW02$order2),]
BW03 <- BW03[order(BW03$order1,BW03$order2),]
BW05 <- BW05[order(BW05$order1,BW05$order2),]
BW12 <- BW12[order(BW12$order1,BW12$order2),]
CB01 <- CB01[order(CB01$order1,CB01$order2),]
CB04 <- CB04[order(CB04$order1,CB04$order2),]
CB12 <- CB12[order(CB12$order1,CB12$order2),]
CB20 <- CB20[order(CB20$order1,CB20$order2),]
BC03 <- BC03[order(BC03$order1,BC03$order2),]
BC04 <- BC04[order(BC04$order1,BC04$order2),]
BC13 <- BC13[order(BC13$order1,BC13$order2),]
BC17 <- BC17[order(BC17$order1,BC17$order2),]

col_order <- c("contig","position","variant","refAllele","altAllele","mat.count","pat.count", "total.count","lowMAPQDepth", "lowBaseQDepth","rawDepth","otherBases","improperPairs","origin")

WC01.ind <- WC01[, col_order]
WC04.ind <- WC04[, col_order]
WC06.ind <- WC06[, col_order]
WC19.ind <- WC19[, col_order]
CW01.ind <- CW01[, col_order]
CW04.ind <- CW04[, col_order]
CW10.ind <- CW10[, col_order]
CW11.ind <- CW11[, col_order]
WB02.ind <- WB02[, col_order]
WB03.ind <- WB03[, col_order]
WB08.ind <- WB08[, col_order]
BW02.ind <- BW02[, col_order]
BW03.ind <- BW03[, col_order]
BW05.ind <- BW05[, col_order]
BW12.ind <- BW12[, col_order]
CB01.ind <- CB01[, col_order]
CB04.ind <- CB04[, col_order]
CB12.ind <- CB12[, col_order]
CB20.ind <- CB20[, col_order]
BC03.ind <- BC03[, col_order]
BC04.ind <- BC04[, col_order]
BC13.ind <- BC13[, col_order]
BC17.ind <- BC17[, col_order]

WC01.ind.hist = ggplot(WC01.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum4") + labs(title="WC01",y=nrow(WC01.ind))
WC04.ind.hist = ggplot(WC04.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum4") + labs(title="WC04",y=nrow(WC04.ind))
WC06.ind.hist = ggplot(WC06.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum4") + labs(title="WC06",y=nrow(WC06.ind))
WC19.ind.hist = ggplot(WC19.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum4") + labs(title="WC19",y=nrow(WC19.ind))
CW01.ind.hist = ggplot(CW01.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum2") + labs(title="CW01",y=nrow(CW01.ind))
CW04.ind.hist = ggplot(CW04.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum2") + labs(title="CW04",y=nrow(CW04.ind))
CW10.ind.hist = ggplot(CW10.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum2") + labs(title="CW10",y=nrow(CW10.ind))
CW11.ind.hist = ggplot(CW11.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum2") + labs(title="CW11",y=nrow(CW11.ind))
WB02.ind.hist = ggplot(WB02.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="cyan4") + labs(title="WB02",y=nrow(WB02.ind))
WB03.ind.hist = ggplot(WB03.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="cyan4") + labs(title="WB03",y=nrow(WB03.ind))
WB08.ind.hist = ggplot(WB08.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="cyan4") + labs(title="WB08",y=nrow(WB08.ind))
BW02.ind.hist = ggplot(BW02.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="cyan3") + labs(title="BW02",y=nrow(BW02.ind))
BW03.ind.hist = ggplot(BW03.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="cyan3") + labs(title="BW03",y=nrow(BW03.ind))
BW05.ind.hist = ggplot(BW05.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="cyan3") + labs(title="BW05",y=nrow(BW05.ind))
BW12.ind.hist = ggplot(BW12.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="cyan3") + labs(title="BW12",y=nrow(BW12.ind))
CB01.ind.hist = ggplot(CB01.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold3") + labs(title="CB01",y=nrow(CB01.ind))
CB04.ind.hist = ggplot(CB04.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold3") + labs(title="CB04",y=nrow(CB04.ind))
CB12.ind.hist = ggplot(CB12.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold3") + labs(title="CB12",y=nrow(CB12.ind))
CB20.ind.hist = ggplot(CB20.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold3") + labs(title="CB20",y=nrow(CB20.ind))
BC03.ind.hist = ggplot(BC03.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold4") + labs(title="BC03",y=nrow(BC03.ind))
BC04.ind.hist = ggplot(BC04.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold4") + labs(title="BC04",y=nrow(BC04.ind))
BC13.ind.hist = ggplot(BC13.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold4") + labs(title="BC13",y=nrow(BC13.ind))
BC17.ind.hist = ggplot(BC17.ind, aes(mat.count/total.count)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold4") + labs(title="BC17",y=nrow(BC17.ind))

grid.arrange(WC01.ind.hist,WC04.ind.hist,WC06.ind.hist,WC19.ind.hist,CW01.ind.hist,CW04.ind.hist,CW10.ind.hist,CW11.ind.hist,
            WB02.ind.hist,WB03.ind.hist,WB08.ind.hist,blank,BW02.ind.hist,BW03.ind.hist,BW05.ind.hist,BW12.ind.hist,
            CB01.ind.hist,CB04.ind.hist,CB12.ind.hist,CB20.ind.hist,BC03.ind.hist,BC04.ind.hist,BC13.ind.hist,BC17.ind.hist,nrow=6)

# prepare dfs to merge by genotype

WC01.to.merge = select(WC01.ind,variant,mat.count,pat.count,total.count)
WC04.to.merge = select(WC04.ind,variant,mat.count,pat.count,total.count)
WC06.to.merge = select(WC06.ind,variant,mat.count,pat.count,total.count)
WC19.to.merge = select(WC19.ind,variant,mat.count,pat.count,total.count)
CW01.to.merge = select(CW01.ind,variant,mat.count,pat.count,total.count)
CW04.to.merge = select(CW04.ind,variant,mat.count,pat.count,total.count)
CW10.to.merge = select(CW10.ind,variant,mat.count,pat.count,total.count)
CW11.to.merge = select(CW11.ind,variant,mat.count,pat.count,total.count)
WB02.to.merge = select(WB02.ind,variant,mat.count,pat.count,total.count)
WB03.to.merge = select(WB03.ind,variant,mat.count,pat.count,total.count)
WB08.to.merge = select(WB08.ind,variant,mat.count,pat.count,total.count)
BW02.to.merge = select(BW02.ind,variant,mat.count,pat.count,total.count)
BW03.to.merge = select(BW03.ind,variant,mat.count,pat.count,total.count)
BW05.to.merge = select(BW05.ind,variant,mat.count,pat.count,total.count)
BW12.to.merge = select(BW12.ind,variant,mat.count,pat.count,total.count)
CB01.to.merge = select(CB01.ind,variant,mat.count,pat.count,total.count)
CB12.to.merge = select(CB12.ind,variant,mat.count,pat.count,total.count)
CB20.to.merge = select(CB20.ind,variant,mat.count,pat.count,total.count)
BC03.to.merge = select(BC03.ind,variant,mat.count,pat.count,total.count)
BC04.to.merge = select(BC04.ind,variant,mat.count,pat.count,total.count)
BC13.to.merge = select(BC13.ind,variant,mat.count,pat.count,total.count)
BC17.to.merge = select(BC17.ind,variant,mat.count,pat.count,total.count)

colnames(WC01.to.merge) = c("variant","mat.count.WC01","pat.count.WC01","total.count.WC01")
colnames(WC04.to.merge) = c("variant","mat.count.WC04","pat.count.WC04","total.count.WC04")
colnames(WC06.to.merge) = c("variant","mat.count.WC06","pat.count.WC06","total.count.WC06")
colnames(WC19.to.merge) = c("variant","mat.count.WC19","pat.count.WC19","total.count.WC19")
colnames(CW01.to.merge) = c("variant","mat.count.CW01","pat.count.CW01","total.count.CW01")
colnames(CW04.to.merge) = c("variant","mat.count.CW04","pat.count.CW04","total.count.CW04")
colnames(CW10.to.merge) = c("variant","mat.count.CW10","pat.count.CW10","total.count.CW10")
colnames(CW11.to.merge) = c("variant","mat.count.CW11","pat.count.CW11","total.count.CW11")
colnames(WB02.to.merge) = c("variant","mat.count.WB02","pat.count.WB02","total.count.WB02")
colnames(WB03.to.merge) = c("variant","mat.count.WB03","pat.count.WB03","total.count.WB03")
colnames(WB08.to.merge) = c("variant","mat.count.WB08","pat.count.WB08","total.count.WB08")
colnames(BW02.to.merge) = c("variant","mat.count.BW02","pat.count.BW02","total.count.BW02")
colnames(BW03.to.merge) = c("variant","mat.count.BW03","pat.count.BW03","total.count.BW03")
colnames(BW05.to.merge) = c("variant","mat.count.BW05","pat.count.BW05","total.count.BW05")
colnames(BW12.to.merge) = c("variant","mat.count.BW12","pat.count.BW12","total.count.BW12")
colnames(CB01.to.merge) = c("variant","mat.count.CB01","pat.count.CB01","total.count.CB01")
colnames(CB12.to.merge) = c("variant","mat.count.CB12","pat.count.CB12","total.count.CB12")
colnames(CB20.to.merge) = c("variant","mat.count.CB20","pat.count.CB20","total.count.CB20")
colnames(BC03.to.merge) = c("variant","mat.count.BC03","pat.count.BC03","total.count.BC03")
colnames(BC04.to.merge) = c("variant","mat.count.BC04","pat.count.BC04","total.count.BC04")
colnames(BC13.to.merge) = c("variant","mat.count.BC13","pat.count.BC13","total.count.BC13")
colnames(BC17.to.merge) = c("variant","mat.count.BC17","pat.count.BC17","total.count.BC17")

# merge by genotype and obtain a list of variants present in at least one replicate

wc.0 = full_join(WC01.to.merge,WC04.to.merge,by="variant")
wc.00 = full_join(wc.0,WC06.to.merge,by="variant")
wc.all = full_join(wc.00,WC19.to.merge,by="variant")
cw.0 = full_join(CW01.to.merge,CW04.to.merge,by="variant")
cw.00 = full_join(cw.0,CW10.to.merge,by="variant")
cw.all = full_join(cw.00,CW11.to.merge,by="variant")

WC.all = full_join(wc.all,cw.all,by="variant")

wb.0 = full_join(WB02.to.merge,WB03.to.merge,by="variant")
wb.all = full_join(wb.0,WB08.to.merge,by="variant")
bw.0 = full_join(BW02.to.merge,BW03.to.merge,by="variant")
bw.00 = full_join(bw.0,BW05.to.merge,by="variant")
bw.all = full_join(bw.00,BW12.to.merge,by="variant")

WB.all = full_join(wb.all,bw.all,by="variant")

cb.0 = full_join(CB01.to.merge,CB12.to.merge,by="variant")
cb.all = full_join(cb.0,CB20.to.merge,by="variant")
bc.0 = full_join(BC03.to.merge,BC04.to.merge,by="variant")
bc.00 = full_join(bc.0,BC13.to.merge,by="variant")
bc.all = full_join(bc.00,BC17.to.merge,by="variant")

CB.all = full_join(cb.all,bc.all,by="variant")

nrow(WC.all)
nrow(WB.all)
nrow(CB.all)

##### Assign genome annotation features to variants

# start by assigning orphan variants (SNPs on contigs without annotated features)

contigs_with_anno <- read_csv("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/6_annotation/genome_annotation/contigs_with_anno.bed", 
                              col_names = FALSE) 
colnames(contigs_with_anno) <- c("contig")
contigs_with_anno_list <- list(unique(contigs_with_anno))

WC.all.status = WC.all[c(1)]
WB.all.status = WB.all[c(1)]
CB.all.status = CB.all[c(1)]

WC.all.status.orphan.0 = WC.all.status
WB.all.status.orphan.0 = WB.all.status
CB.all.status.orphan.0 = CB.all.status

WC.all.status.orphan.0$contig = sub(":.*", '', WC.all.status.orphan.0$variant)
WB.all.status.orphan.0$contig = sub(":.*", '', WB.all.status.orphan.0$variant)
CB.all.status.orphan.0$contig = sub(":.*", '', CB.all.status.orphan.0$variant)

WC.all.status.orphan.0$status = ifelse(WC.all.status.orphan.0$contig %in% contigs_with_anno$contig == TRUE,NA,"orphan")
WB.all.status.orphan.0$status = ifelse(WB.all.status.orphan.0$contig %in% contigs_with_anno$contig == TRUE,NA,"orphan")
CB.all.status.orphan.0$status = ifelse(CB.all.status.orphan.0$contig %in% contigs_with_anno$contig == TRUE,NA,"orphan")

WC.all.status.orphan.0 = WC.all.status.orphan.0[complete.cases(WC.all.status.orphan.0),]
WB.all.status.orphan.0 = WB.all.status.orphan.0[complete.cases(WB.all.status.orphan.0),]
CB.all.status.orphan.0 = CB.all.status.orphan.0[complete.cases(CB.all.status.orphan.0),]

WC.all.status.orphan <- WC.all.status.orphan.0[c(1,3)]
WB.all.status.orphan <- WB.all.status.orphan.0[c(1,3)]
CB.all.status.orphan <- CB.all.status.orphan.0[c(1,3)]

WC.all.status.orphan$CDS <- NA
WB.all.status.orphan$CDS <- NA
CB.all.status.orphan$CDS <- NA

WC.all.status.orphan$mRNA <- NA
WB.all.status.orphan$mRNA <- NA
CB.all.status.orphan$mRNA <- NA

# to assign rest of categories: retrieve counts intersected with annotation

BC03.bc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC03.bc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC03.cb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC03.cb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC04.bc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC04.bc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC04.cb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC04.cb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC13.bc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC13.bc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC13.cb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC13.cb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC17.bc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC17.bc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC17.cb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BC17.cb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BW02.bw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW02.bw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BW02.wb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW02.wb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BW03.bw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW03.bw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BW03.wb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW03.wb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BW05.bw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW05.bw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BW05.wb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW05.wb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BW12.bw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW12.bw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BW12.wb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/BW12.wb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CB01.bc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB01.bc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CB01.cb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB01.cb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CB04.bc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB04.bc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CB04.cb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB04.cb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CB12.bc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB12.bc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CB12.cb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB12.cb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CB20.bc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB20.bc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CB20.cb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CB20.cb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CW01.cw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW01.cw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CW01.wc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW01.wc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CW04.cw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW04.cw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CW04.wc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW04.wc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CW10.cw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW10.cw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CW10.wc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW10.wc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CW11.cw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW11.cw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
CW11.wc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/CW11.wc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WB02.bw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB02.bw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WB02.wb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB02.wb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WB03.bw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB03.bw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WB03.wb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB03.wb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WB08.bw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB08.bw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WB08.wb.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WB08.wb.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WC01.cw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC01.cw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WC01.wc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC01.wc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WC04.cw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC04.cw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WC04.wc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC04.wc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WC06.cw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC06.cw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WC06.wc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC06.wc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WC19.cw.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC19.cw.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
WC19.wc.anno.R <- read_delim("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/5_ase/results/ASE_output/WC19.wc.anno.intra.R.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

WC.anno.all = rbind(WC01.cw.anno.R,WC01.wc.anno.R,WC04.cw.anno.R,WC04.wc.anno.R,WC06.cw.anno.R,WC06.wc.anno.R,WC19.cw.anno.R,WC19.wc.anno.R,
                    CW01.cw.anno.R,CW01.wc.anno.R,CW04.cw.anno.R,CW04.wc.anno.R,CW10.cw.anno.R,CW10.wc.anno.R,CW11.cw.anno.R,CW11.wc.anno.R)
WB.anno.all = rbind(WB02.bw.anno.R,WB02.wb.anno.R,WB03.bw.anno.R,WB03.wb.anno.R,WB08.bw.anno.R,WB08.wb.anno.R,
                    BW02.bw.anno.R,BW02.wb.anno.R,BW03.bw.anno.R,BW03.wb.anno.R,BW05.bw.anno.R,BW05.wb.anno.R,BW12.bw.anno.R,BW12.wb.anno.R)
CB.anno.all = rbind(CB01.bc.anno.R,CB01.cb.anno.R,CB12.bc.anno.R,CB12.cb.anno.R,CB20.bc.anno.R,CB20.cb.anno.R,
                    BC03.bc.anno.R,BC03.cb.anno.R,BC04.bc.anno.R,BC04.cb.anno.R,BC13.bc.anno.R,BC13.cb.anno.R,BC17.bc.anno.R,BC17.cb.anno.R)

WC.anno.all.unique <- distinct(WC.anno.all)
WB.anno.all.unique <- distinct(WB.anno.all)
CB.anno.all.unique <- distinct(CB.anno.all)

colnames(WC.anno.all.unique) = c("contig","position","anno","type","info")
colnames(WB.anno.all.unique) = c("contig","position","anno","type","info")
colnames(CB.anno.all.unique) = c("contig","position","anno","type","info")

WC.anno.all.unique$variant = paste(WC.anno.all.unique$contig,WC.anno.all.unique$position,sep=":")
WB.anno.all.unique$variant = paste(WB.anno.all.unique$contig,WB.anno.all.unique$position,sep=":")
CB.anno.all.unique$variant = paste(CB.anno.all.unique$contig,CB.anno.all.unique$position,sep=":")

nrow(WC.anno.all.unique)
nrow(WB.anno.all.unique)
nrow(CB.anno.all.unique)

# start with exonic sites, so that exonic SNPs that map to other features (e.g. introns) will be given priority

exonic.sites.WC <- WC.anno.all.unique[WC.anno.all.unique$type == "CDS",]
exonic.sites.WB <- WB.anno.all.unique[WB.anno.all.unique$type == "CDS",]
exonic.sites.CB <- CB.anno.all.unique[CB.anno.all.unique$type == "CDS",]

exonic.sites.WC.merge <- exonic.sites.WC[c(6,5)]
exonic.sites.WB.merge <- exonic.sites.WB[c(6,5)]
exonic.sites.CB.merge <- exonic.sites.CB[c(6,5)]

exonic.sites.WC.merge <- separate(exonic.sites.WC.merge,info,into = c("CDS","mRNA"),sep = ";") # obtain information from the annotation field
exonic.sites.WB.merge <- separate(exonic.sites.WB.merge,info,into = c("CDS","mRNA"),sep = ";")
exonic.sites.CB.merge <- separate(exonic.sites.CB.merge,info,into = c("CDS","mRNA"),sep = ";")

exonic.sites.WC.merge$CDS <- sub('ID=', '', exonic.sites.WC.merge$CDS)
exonic.sites.WB.merge$CDS <- sub('ID=', '', exonic.sites.WB.merge$CDS)
exonic.sites.CB.merge$CDS <- sub('ID=', '', exonic.sites.CB.merge$CDS)
exonic.sites.WC.merge$mRNA <- sub('Parent=', '', exonic.sites.WC.merge$mRNA)
exonic.sites.WB.merge$mRNA <- sub('Parent=', '', exonic.sites.WB.merge$mRNA)
exonic.sites.CB.merge$mRNA <- sub('Parent=', '', exonic.sites.CB.merge$mRNA)

WC.all.status.exon=ddply(exonic.sites.WC.merge,c("variant"), summarize,
                                      CDS=paste(CDS,collapse=","), 
                                      mRNA=paste(mRNA,collapse=",")) # collapse annotations for variants with >1 feature

WB.all.status.exon=ddply(exonic.sites.WB.merge,c("variant"), summarize,
                                      CDS=paste(CDS,collapse=","), 
                                      mRNA=paste(mRNA,collapse=","))

CB.all.status.exon=ddply(exonic.sites.CB.merge,c("variant"), summarize,
                                      CDS=paste(CDS,collapse=","), 
                                      mRNA=paste(mRNA,collapse=","))

WC.all.status.exon$status <- "exonic"
WB.all.status.exon$status <- "exonic"
CB.all.status.exon$status <- "exonic"

# now assign all intronic sites

intronic.sites.WC <- WC.anno.all.unique[WC.anno.all.unique$type == "intron",]
intronic.sites.WB <- WB.anno.all.unique[WB.anno.all.unique$type == "intron",]
intronic.sites.CB <- CB.anno.all.unique[CB.anno.all.unique$type == "intron",]

intronic.sites.WC.merge <- intronic.sites.WC[c(6,5)]
intronic.sites.WB.merge <- intronic.sites.WB[c(6,5)]
intronic.sites.CB.merge <- intronic.sites.CB[c(6,5)]

intronic.sites.WC.merge$mRNA <- sub('Parent=', '', intronic.sites.WC.merge$info)
intronic.sites.WB.merge$mRNA <- sub('Parent=', '', intronic.sites.WB.merge$info)
intronic.sites.CB.merge$mRNA <- sub('Parent=', '', intronic.sites.CB.merge$info)

intronic.sites.WC.merge.collapsed = ddply(intronic.sites.WC.merge,c("variant"), summarize,mRNA=paste(mRNA,collapse=",")) # collapse annotations for variants with >1 feature
intronic.sites.WB.merge.collapsed = ddply(intronic.sites.WB.merge,c("variant"), summarize,mRNA=paste(mRNA,collapse=","))
intronic.sites.CB.merge.collapsed = ddply(intronic.sites.CB.merge,c("variant"), summarize,mRNA=paste(mRNA,collapse=","))

WC.all.status.intron = anti_join(intronic.sites.WC.merge.collapsed,WC.all.status.exon,by="variant") # exclude intronic variants if they already map to exons
WB.all.status.intron = anti_join(intronic.sites.WB.merge.collapsed,WB.all.status.exon,by="variant")
CB.all.status.intron = anti_join(intronic.sites.CB.merge.collapsed,CB.all.status.exon,by="variant")

WC.all.status.intron$status <- "intronic"
WB.all.status.intron$status <- "intronic"
CB.all.status.intron$status <- "intronic"

WC.all.status.intron$CDS <- NA
WB.all.status.intron$CDS <- NA
CB.all.status.intron$CDS <- NA

# assign intergenic sites

WC.all.status.non.intergenic = rbind(WC.all.status.exon,WC.all.status.intron,WC.all.status.orphan) # merge all variants with exonic/intronic/orphan status
WB.all.status.non.intergenic = rbind(WB.all.status.exon,WB.all.status.intron,WB.all.status.orphan)
CB.all.status.non.intergenic = rbind(CB.all.status.exon,CB.all.status.intron,CB.all.status.orphan)

WC.all.status.intergenic.0 = anti_join(WC.all.status,WC.all.status.non.intergenic,by="variant") # assign intergenic status to all remaining variants
WB.all.status.intergenic.0 = anti_join(WB.all.status,WB.all.status.non.intergenic,by="variant")
CB.all.status.intergenic.0 = anti_join(CB.all.status,CB.all.status.non.intergenic,by="variant")

WC.all.status.intergenic <- WC.all.status.intergenic.0[c(1)]
WB.all.status.intergenic <- WB.all.status.intergenic.0[c(1)]
CB.all.status.intergenic <- CB.all.status.intergenic.0[c(1)]

WC.all.status.intergenic$CDS <- NA
WB.all.status.intergenic$CDS <- NA
CB.all.status.intergenic$CDS <- NA

WC.all.status.intergenic$mRNA <- NA
WB.all.status.intergenic$mRNA <- NA
CB.all.status.intergenic$mRNA <- NA

WC.all.status.intergenic$status <- "intergenic"
WB.all.status.intergenic$status <- "intergenic"
CB.all.status.intergenic$status <- "intergenic"

nrow(WC.all.status.exon)+nrow(WC.all.status.intron)+nrow(WC.all.status.intergenic)+nrow(WC.all.status.orphan) - nrow(WC.all) # sanity check
nrow(WB.all.status.exon)+nrow(WB.all.status.intron)+nrow(WB.all.status.intergenic)+nrow(WB.all.status.orphan) - nrow(WB.all)
nrow(CB.all.status.exon)+nrow(CB.all.status.intron)+nrow(CB.all.status.intergenic)+nrow(CB.all.status.orphan) - nrow(CB.all)

# merge all

WC.all.anno = rbind(WC.all.status.non.intergenic,WC.all.status.intergenic)
WB.all.anno = rbind(WB.all.status.non.intergenic,WB.all.status.intergenic)
CB.all.anno = rbind(CB.all.status.non.intergenic,CB.all.status.intergenic)

WC.annotation.counts = ddply(WC.all.anno,c("status"),summarise,N = length(variant)) # get counts
WB.annotation.counts = ddply(WB.all.anno,c("status"),summarise,N = length(variant)) # get counts
CB.annotation.counts = ddply(CB.all.anno,c("status"),summarise,N = length(variant)) # get counts

# intersect annotation and read counts

WC.all.complete.not.ordered = full_join(WC.all.anno,WC.all,by="variant")
WB.all.complete.not.ordered = full_join(WB.all.anno,WB.all,by="variant")
CB.all.complete.not.ordered = full_join(CB.all.anno,CB.all,by="variant")

nrow(WC.all.complete.not.ordered)
nrow(WB.all.complete.not.ordered)
nrow(CB.all.complete.not.ordered)

WC.all.complete.not.ordered$order1 = as.integer(gsub(".*_", '', (sub(":.*", '', WC.all.complete.not.ordered$variant))))
WB.all.complete.not.ordered$order1 = as.integer(gsub(".*_", '', (sub(":.*", '', WB.all.complete.not.ordered$variant))))
CB.all.complete.not.ordered$order1 = as.integer(gsub(".*_", '', (sub(":.*", '', CB.all.complete.not.ordered$variant))))

WC.all.complete.not.ordered$order2 = as.integer(gsub(".*:", '', WC.all.complete.not.ordered$variant))
WB.all.complete.not.ordered$order2 = as.integer(gsub(".*:", '', WB.all.complete.not.ordered$variant))
CB.all.complete.not.ordered$order2 = as.integer(gsub(".*:", '', CB.all.complete.not.ordered$variant))

WC.all.complete = WC.all.complete.not.ordered[order(WC.all.complete.not.ordered$order1, WC.all.complete.not.ordered$order2),]
WB.all.complete = WB.all.complete.not.ordered[order(WB.all.complete.not.ordered$order1, WB.all.complete.not.ordered$order2),]
CB.all.complete = CB.all.complete.not.ordered[order(CB.all.complete.not.ordered$order1, CB.all.complete.not.ordered$order2),]

nrow(WC.all.complete)
nrow(WB.all.complete)
nrow(CB.all.complete)

##### Filter variants

# initial filter: keep variants shared in at least 3 of 3/4 reciprocal genotype replicates

WC.all.complete$in.WC01 = ifelse(is.na(WC.all.complete$total.count.WC01),0,1)
WC.all.complete$in.WC04 = ifelse(is.na(WC.all.complete$total.count.WC04),0,1)
WC.all.complete$in.WC06 = ifelse(is.na(WC.all.complete$total.count.WC06),0,1)
WC.all.complete$in.WC19 = ifelse(is.na(WC.all.complete$total.count.WC19),0,1)
WC.all.complete$in.CW01 = ifelse(is.na(WC.all.complete$total.count.CW01),0,1)
WC.all.complete$in.CW04 = ifelse(is.na(WC.all.complete$total.count.CW04),0,1)
WC.all.complete$in.CW10 = ifelse(is.na(WC.all.complete$total.count.CW10),0,1)
WC.all.complete$in.CW11 = ifelse(is.na(WC.all.complete$total.count.CW11),0,1)

WC.all.complete$num.rep.WC = WC.all.complete$in.WC01+WC.all.complete$in.WC04+WC.all.complete$in.WC06+WC.all.complete$in.WC19
WC.all.complete$num.rep.CW = WC.all.complete$in.CW01+WC.all.complete$in.CW04+WC.all.complete$in.CW10+WC.all.complete$in.CW11
WC.all.complete$shared = ifelse((WC.all.complete$num.rep.WC)>=3 & (WC.all.complete$num.rep.CW)>=3,"Y","N")
WC.shared.complete <- WC.all.complete[WC.all.complete$shared == "Y",]

WB.all.complete$in.WB02 = ifelse(is.na(WB.all.complete$total.count.WB02),0,1)
WB.all.complete$in.WB03 = ifelse(is.na(WB.all.complete$total.count.WB03),0,1)
WB.all.complete$in.WB08 = ifelse(is.na(WB.all.complete$total.count.WB08),0,1)
WB.all.complete$in.BW02 = ifelse(is.na(WB.all.complete$total.count.BW02),0,1)
WB.all.complete$in.BW03 = ifelse(is.na(WB.all.complete$total.count.BW03),0,1)
WB.all.complete$in.BW05 = ifelse(is.na(WB.all.complete$total.count.BW05),0,1)
WB.all.complete$in.BW12 = ifelse(is.na(WB.all.complete$total.count.BW12),0,1)

WB.all.complete$num.rep.WB = WB.all.complete$in.WB02+WB.all.complete$in.WB03+WB.all.complete$in.WB08
WB.all.complete$num.rep.BW = WB.all.complete$in.BW02+WB.all.complete$in.BW03+WB.all.complete$in.BW05+WB.all.complete$in.BW12
WB.all.complete$shared = ifelse((WB.all.complete$num.rep.WB)>=3 & (WB.all.complete$num.rep.BW)>=3,"Y","N")
WB.shared.complete <- WB.all.complete[WB.all.complete$shared == "Y",]

CB.all.complete$in.CB01 = ifelse(is.na(CB.all.complete$total.count.CB01),0,1)
CB.all.complete$in.CB12 = ifelse(is.na(CB.all.complete$total.count.CB12),0,1)
CB.all.complete$in.CB20 = ifelse(is.na(CB.all.complete$total.count.CB20),0,1)
CB.all.complete$in.BC03 = ifelse(is.na(CB.all.complete$total.count.BC03),0,1)
CB.all.complete$in.BC04 = ifelse(is.na(CB.all.complete$total.count.BC04),0,1)
CB.all.complete$in.BC13 = ifelse(is.na(CB.all.complete$total.count.BC13),0,1)
CB.all.complete$in.BC17 = ifelse(is.na(CB.all.complete$total.count.BC17),0,1)

CB.all.complete$num.rep.CB = CB.all.complete$in.CB01+CB.all.complete$in.CB12+CB.all.complete$in.CB20
CB.all.complete$num.rep.BC = CB.all.complete$in.BC03+CB.all.complete$in.BC04+CB.all.complete$in.BC13+CB.all.complete$in.BC17
CB.all.complete$shared = ifelse((CB.all.complete$num.rep.CB)>=3 & (CB.all.complete$num.rep.BC)>=3,"Y","N")
CB.shared.complete <- CB.all.complete[CB.all.complete$shared == "Y",]

nrow(WC.shared.complete)
nrow(WB.shared.complete)
nrow(CB.shared.complete)

WC.shared.annotation.counts = ddply(WC.shared.complete,c("status"),summarise,N = length(variant),perc=N*100/nrow(WC.shared.complete)) # get counts
WB.shared.annotation.counts = ddply(WB.shared.complete,c("status"),summarise,N = length(variant),perc=N*100/nrow(WB.shared.complete)) # get counts
CB.shared.annotation.counts = ddply(CB.shared.complete,c("status"),summarise,N = length(variant),perc=N*100/nrow(CB.shared.complete)) # get counts

# filter: remove SNPs in which valid read depth >90% of total read depth

WC01.ind$prop_other=(WC01.ind$otherBases+WC01.ind$lowMAPQDepth)/(WC01.ind$total.count+WC01.ind$otherBases+WC01.ind$lowMAPQDepth)
WC04.ind$prop_other=(WC04.ind$otherBases+WC04.ind$lowMAPQDepth)/(WC04.ind$total.count+WC04.ind$otherBases+WC04.ind$lowMAPQDepth)
WC06.ind$prop_other=(WC06.ind$otherBases+WC06.ind$lowMAPQDepth)/(WC06.ind$total.count+WC06.ind$otherBases+WC06.ind$lowMAPQDepth)
WC19.ind$prop_other=(WC19.ind$otherBases+WC19.ind$lowMAPQDepth)/(WC19.ind$total.count+WC19.ind$otherBases+WC19.ind$lowMAPQDepth)
CW01.ind$prop_other=(CW01.ind$otherBases+CW01.ind$lowMAPQDepth)/(CW01.ind$total.count+CW01.ind$otherBases+CW01.ind$lowMAPQDepth)
CW04.ind$prop_other=(CW04.ind$otherBases+CW04.ind$lowMAPQDepth)/(CW04.ind$total.count+CW04.ind$otherBases+CW04.ind$lowMAPQDepth)
CW10.ind$prop_other=(CW10.ind$otherBases+CW10.ind$lowMAPQDepth)/(CW10.ind$total.count+CW10.ind$otherBases+CW10.ind$lowMAPQDepth)
CW11.ind$prop_other=(CW11.ind$otherBases+CW11.ind$lowMAPQDepth)/(CW11.ind$total.count+CW11.ind$otherBases+CW11.ind$lowMAPQDepth)
WB02.ind$prop_other=(WB02.ind$otherBases+WB02.ind$lowMAPQDepth)/(WB02.ind$total.count+WB02.ind$otherBases+WB02.ind$lowMAPQDepth)
WB03.ind$prop_other=(WB03.ind$otherBases+WB03.ind$lowMAPQDepth)/(WB03.ind$total.count+WB03.ind$otherBases+WB03.ind$lowMAPQDepth)
WB08.ind$prop_other=(WB08.ind$otherBases+WB08.ind$lowMAPQDepth)/(WB08.ind$total.count+WB08.ind$otherBases+WB08.ind$lowMAPQDepth)
BW02.ind$prop_other=(BW02.ind$otherBases+BW02.ind$lowMAPQDepth)/(BW02.ind$total.count+BW02.ind$otherBases+BW02.ind$lowMAPQDepth)
BW03.ind$prop_other=(BW03.ind$otherBases+BW03.ind$lowMAPQDepth)/(BW03.ind$total.count+BW03.ind$otherBases+BW03.ind$lowMAPQDepth)
BW05.ind$prop_other=(BW05.ind$otherBases+BW05.ind$lowMAPQDepth)/(BW05.ind$total.count+BW05.ind$otherBases+BW05.ind$lowMAPQDepth)
BW12.ind$prop_other=(BW12.ind$otherBases+BW12.ind$lowMAPQDepth)/(BW12.ind$total.count+BW12.ind$otherBases+BW12.ind$lowMAPQDepth)
CB01.ind$prop_other=(CB01.ind$otherBases+CB01.ind$lowMAPQDepth)/(CB01.ind$total.count+CB01.ind$otherBases+CB01.ind$lowMAPQDepth)
CB12.ind$prop_other=(CB12.ind$otherBases+CB12.ind$lowMAPQDepth)/(CB12.ind$total.count+CB12.ind$otherBases+CB12.ind$lowMAPQDepth)
CB20.ind$prop_other=(CB20.ind$otherBases+CB20.ind$lowMAPQDepth)/(CB20.ind$total.count+CB20.ind$otherBases+CB20.ind$lowMAPQDepth)
BC03.ind$prop_other=(BC03.ind$otherBases+BC03.ind$lowMAPQDepth)/(BC03.ind$total.count+BC03.ind$otherBases+BC03.ind$lowMAPQDepth)
BC04.ind$prop_other=(BC04.ind$otherBases+BC04.ind$lowMAPQDepth)/(BC04.ind$total.count+BC04.ind$otherBases+BC04.ind$lowMAPQDepth)
BC13.ind$prop_other=(BC13.ind$otherBases+BC13.ind$lowMAPQDepth)/(BC13.ind$total.count+BC13.ind$otherBases+BC13.ind$lowMAPQDepth)
BC17.ind$prop_other=(BC17.ind$otherBases+BC17.ind$lowMAPQDepth)/(BC17.ind$total.count+BC17.ind$otherBases+BC17.ind$lowMAPQDepth)

WC01.other.bases <- WC01.ind[(WC01.ind$prop_other > 0.10),]
WC04.other.bases <- WC04.ind[(WC04.ind$prop_other > 0.10),]
WC06.other.bases <- WC06.ind[(WC06.ind$prop_other > 0.10),]
WC19.other.bases <- WC19.ind[(WC19.ind$prop_other > 0.10),]
CW01.other.bases <- CW01.ind[(CW01.ind$prop_other > 0.10),]
CW04.other.bases <- CW04.ind[(CW04.ind$prop_other > 0.10),]
CW10.other.bases <- CW10.ind[(CW10.ind$prop_other > 0.10),]
CW11.other.bases <- CW11.ind[(CW11.ind$prop_other > 0.10),]
WB02.other.bases <- WB02.ind[(WB02.ind$prop_other > 0.10),]
WB03.other.bases <- WB03.ind[(WB03.ind$prop_other > 0.10),]
WB08.other.bases <- WB08.ind[(WB08.ind$prop_other > 0.10),]
BW02.other.bases <- BW02.ind[(BW02.ind$prop_other > 0.10),]
BW03.other.bases <- BW03.ind[(BW03.ind$prop_other > 0.10),]
BW05.other.bases <- BW05.ind[(BW05.ind$prop_other > 0.10),]
BW12.other.bases <- BW12.ind[(BW12.ind$prop_other > 0.10),]
CB01.other.bases <- CB01.ind[(CB01.ind$prop_other > 0.10),]
CB12.other.bases <- CB12.ind[(CB12.ind$prop_other > 0.10),]
CB20.other.bases <- CB20.ind[(CB20.ind$prop_other > 0.10),]
BC03.other.bases <- BC03.ind[(BC03.ind$prop_other > 0.10),]
BC04.other.bases <- BC04.ind[(BC04.ind$prop_other > 0.10),]
BC13.other.bases <- BC13.ind[(BC13.ind$prop_other > 0.10),]
BC17.other.bases <- BC17.ind[(BC17.ind$prop_other > 0.10),]

WC.other.bases = rbind(WC01.other.bases, WC04.other.bases, WC06.other.bases, WC19.other.bases, CW01.other.bases, CW04.other.bases, CW10.other.bases, CW11.other.bases)
WB.other.bases = rbind(WB02.other.bases, WB03.other.bases, WB08.other.bases, BW02.other.bases, BW03.other.bases, BW05.other.bases, BW12.other.bases)
CB.other.bases = rbind(CB01.other.bases, CB12.other.bases, CB20.other.bases, BC03.other.bases, BC04.other.bases, BC13.other.bases, BC17.other.bases)

WC.other.bases.id = unique(WC.other.bases[c(3)])
WB.other.bases.id = unique(WB.other.bases[c(3)])
CB.other.bases.id = unique(CB.other.bases[c(3)])

WC.shared.filtered = anti_join(WC.shared.complete,WC.other.bases,by=c("variant"))
WB.shared.filtered = anti_join(WB.shared.complete,WB.other.bases,by=c("variant"))
CB.shared.filtered = anti_join(CB.shared.complete,CB.other.bases,by=c("variant"))

WC.shared.filtered.annotation.counts = ddply(WC.shared.filtered,c("status"),summarise,N = length(variant),perc=N*100/nrow(WC.shared.filtered)) # get counts
WB.shared.filtered.annotation.counts = ddply(WB.shared.filtered,c("status"),summarise,N = length(variant),perc=N*100/nrow(WB.shared.filtered)) # get counts
CB.shared.filtered.annotation.counts = ddply(CB.shared.filtered,c("status"),summarise,N = length(variant),perc=N*100/nrow(CB.shared.filtered)) # get counts

nrow(WC.shared.filtered)
nrow(WB.shared.filtered)
nrow(CB.shared.filtered)

##################################

##### collect stats

# define a function which maps NA to 0 like this:
  
na2zero <- function(x) ifelse(is.na(x), 0, x)

WC.shared.filtered$bias.WC01 = round(WC.shared.filtered$mat.count.WC01/WC.shared.filtered$total.count.WC01,3)
WC.shared.filtered$bias.WC04 = round(WC.shared.filtered$mat.count.WC04/WC.shared.filtered$total.count.WC04,3)
WC.shared.filtered$bias.WC06 = round(WC.shared.filtered$mat.count.WC06/WC.shared.filtered$total.count.WC06,3)
WC.shared.filtered$bias.WC19 = round(WC.shared.filtered$mat.count.WC19/WC.shared.filtered$total.count.WC19,3)
WC.shared.filtered$bias.CW01 = round(WC.shared.filtered$mat.count.CW01/WC.shared.filtered$total.count.CW01,3)
WC.shared.filtered$bias.CW04 = round(WC.shared.filtered$mat.count.CW04/WC.shared.filtered$total.count.CW04,3)
WC.shared.filtered$bias.CW10 = round(WC.shared.filtered$mat.count.CW10/WC.shared.filtered$total.count.CW10,3)
WC.shared.filtered$bias.CW11 = round(WC.shared.filtered$mat.count.CW11/WC.shared.filtered$total.count.CW11,3)

# obtain average bias per reciprocal genotype

WC.shared.filtered$bias.WC = round(((na2zero(WC.shared.filtered$bias.WC01)+na2zero(WC.shared.filtered$bias.WC04)+na2zero(WC.shared.filtered$bias.WC06)+na2zero(WC.shared.filtered$bias.WC19))/na2zero(WC.shared.filtered$num.rep.WC)),3)
WC.shared.filtered$bias.CW = round(((na2zero(WC.shared.filtered$bias.CW01)+na2zero(WC.shared.filtered$bias.CW04)+na2zero(WC.shared.filtered$bias.CW10)+na2zero(WC.shared.filtered$bias.CW11))/na2zero(WC.shared.filtered$num.rep.CW)),3)

WB.shared.filtered$bias.WB02 = round(WB.shared.filtered$mat.count.WB02/WB.shared.filtered$total.count.WB02,3)
WB.shared.filtered$bias.WB03 = round(WB.shared.filtered$mat.count.WB03/WB.shared.filtered$total.count.WB03,3)
WB.shared.filtered$bias.WB08 = round(WB.shared.filtered$mat.count.WB08/WB.shared.filtered$total.count.WB08,3)
WB.shared.filtered$bias.BW02 = round(WB.shared.filtered$mat.count.BW02/WB.shared.filtered$total.count.BW02,3)
WB.shared.filtered$bias.BW03 = round(WB.shared.filtered$mat.count.BW03/WB.shared.filtered$total.count.BW03,3)
WB.shared.filtered$bias.BW05 = round(WB.shared.filtered$mat.count.BW05/WB.shared.filtered$total.count.BW05,3)
WB.shared.filtered$bias.BW12 = round(WB.shared.filtered$mat.count.BW12/WB.shared.filtered$total.count.BW12,3)

WB.shared.filtered$bias.WB = round(((na2zero(WB.shared.filtered$bias.WB02)+na2zero(WB.shared.filtered$bias.WB03)+na2zero(WB.shared.filtered$bias.WB08))/na2zero(WB.shared.filtered$num.rep.WB)),3)
WB.shared.filtered$bias.BW = round(((na2zero(WB.shared.filtered$bias.BW02)+na2zero(WB.shared.filtered$bias.BW03)+na2zero(WB.shared.filtered$bias.BW05)+na2zero(WB.shared.filtered$bias.BW12))/na2zero(WB.shared.filtered$num.rep.BW)),3)

CB.shared.filtered$bias.CB01 = round(CB.shared.filtered$mat.count.CB01/CB.shared.filtered$total.count.CB01,3)
CB.shared.filtered$bias.CB12 = round(CB.shared.filtered$mat.count.CB12/CB.shared.filtered$total.count.CB12,3)
CB.shared.filtered$bias.CB20 = round(CB.shared.filtered$mat.count.CB20/CB.shared.filtered$total.count.CB20,3)
CB.shared.filtered$bias.BC03 = round(CB.shared.filtered$mat.count.BC03/CB.shared.filtered$total.count.BC03,3)
CB.shared.filtered$bias.BC04 = round(CB.shared.filtered$mat.count.BC04/CB.shared.filtered$total.count.BC04,3)
CB.shared.filtered$bias.BC13 = round(CB.shared.filtered$mat.count.BC13/CB.shared.filtered$total.count.BC13,3)
CB.shared.filtered$bias.BC17 = round(CB.shared.filtered$mat.count.BC17/CB.shared.filtered$total.count.BC17,3)

CB.shared.filtered$bias.CB = round(((na2zero(CB.shared.filtered$bias.CB01)+na2zero(CB.shared.filtered$bias.CB12)+na2zero(CB.shared.filtered$bias.CB20))/na2zero(CB.shared.filtered$num.rep.CB)),3)
CB.shared.filtered$bias.BC = round(((na2zero(CB.shared.filtered$bias.BC03)+na2zero(CB.shared.filtered$bias.BC04)+na2zero(CB.shared.filtered$bias.BC13)+na2zero(CB.shared.filtered$bias.BC17))/na2zero(CB.shared.filtered$num.rep.BC)),3)

mean(WC.shared.filtered$bias.WC)
mean(WC.shared.filtered$bias.CW)
mean(WB.shared.filtered$bias.WB)
mean(WB.shared.filtered$bias.BW)
mean(CB.shared.filtered$bias.CB)
mean(CB.shared.filtered$bias.BC)

sd(WC.shared.filtered$bias.WC)
sd(WC.shared.filtered$bias.CW)
sd(WB.shared.filtered$bias.WB)
sd(WB.shared.filtered$bias.BW)
sd(CB.shared.filtered$bias.CB)
sd(CB.shared.filtered$bias.BC)

WC.bias.summary = ddply(WC.shared.filtered, c("status"), summarise,
                             N = length(variant), perc = length(variant)*100/nrow(WC.shared.filtered),
                             bias.WC = mean(bias.WC),sd.WC=sd(bias.WC),
                             bias.CW = mean(bias.CW),sd.CW=sd(bias.CW))

WB.bias.summary = ddply(WB.shared.filtered, c("status"), summarise,
                        N = length(variant), perc = length(variant)*100/nrow(WB.shared.filtered),
                        bias.WB = mean(bias.WB),sd.WB=sd(bias.WB),
                        bias.BW = mean(bias.BW),sd.BW=sd(bias.BW))

CB.bias.summary = ddply(CB.shared.filtered, c("status"), summarise,
                        N = length(variant), perc = length(variant)*100/nrow(CB.shared.filtered),
                        bias.CB = mean(bias.CB),sd.CB=sd(bias.CB),
                        bias.BC = mean(bias.BC),sd.BC=sd(bias.BC))

WC.exonic.only = WC.shared.filtered[WC.shared.filtered$status == "exonic",]
WB.exonic.only = WB.shared.filtered[WB.shared.filtered$status == "exonic",]
CB.exonic.only = CB.shared.filtered[CB.shared.filtered$status == "exonic",]

nrow(WC.exonic.only) / nrow(WC.shared.filtered)
nrow(WB.exonic.only) / nrow(WB.shared.filtered)
nrow(CB.exonic.only) / nrow(CB.shared.filtered)

nrow(WC.shared.filtered)
nrow(WB.shared.filtered)
nrow(CB.shared.filtered)

# fig. 2b: plot ASE at SNP level (histogram)

WC.bias.SNP.hist = ggplot(WC.shared.filtered, aes(bias.WC)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum4") + labs(title="WC",x="", y="SNPs") + ylim(0, 5000) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  geom_vline(xintercept = 0.5,linetype="dotted", color = "gray40") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

CW.bias.SNP.hist = ggplot(WC.shared.filtered, aes(bias.CW)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="plum2") + labs(title="CW",x="", y="SNPs") + ylim(0, 5000) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  geom_vline(xintercept = 0.5,linetype="dotted", color = "gray40") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

WB.bias.SNP.hist = ggplot(WB.shared.filtered, aes(bias.WB)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="cyan4") + labs(title="WB",x="", y="SNPs") + ylim(0, 5000) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  geom_vline(xintercept = 0.5,linetype="dotted", color = "gray40") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

BW.bias.SNP.hist = ggplot(WB.shared.filtered, aes(bias.BW)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="cyan3") + labs(title="BW",x="", y="SNPs") + ylim(0, 5000) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  geom_vline(xintercept = 0.5,linetype="dotted", color = "gray40") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

CB.bias.SNP.hist = ggplot(CB.shared.filtered, aes(bias.CB)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold4") + labs(title="CB",x="", y="SNPs") + ylim(0, 200) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  geom_vline(xintercept = 0.5,linetype="dotted", color = "gray40") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

BC.bias.SNP.hist = ggplot(CB.shared.filtered, aes(bias.BC)) + geom_histogram(breaks=seq(0,1, by=0.02),fill="gold3") + labs(title="BC",x="", y="SNPs") + ylim(0, 200) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  geom_vline(xintercept = 0.5,linetype="dotted", color = "gray40") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

#ggplot(WC.exonic.only, aes(bias.WC,bias.CW)) + geom_point()
#ggplot(WB.exonic.only, aes(bias.WB,bias.BW)) + geom_point()
#ggplot(CB.exonic.only, aes(bias.CB,bias.BC)) + geom_point()

# plot ASE bias at SNP level (barplot)

nrow(WC.exonic.only)
nrow(WB.exonic.only)
nrow(CB.exonic.only)

count.WC.mat = nrow(WC.exonic.only[(WC.exonic.only$bias.WC >= 0.950),])
count.WC.mat.bias = nrow(WC.exonic.only[(WC.exonic.only$bias.WC > 0.650 & WC.exonic.only$bias.WC < 0.950),])
count.WC.bip = nrow(WC.exonic.only[(WC.exonic.only$bias.WC <= 0.650 & WC.exonic.only$bias.WC >= 0.350),])
count.WC.pat.bias = nrow(WC.exonic.only[(WC.exonic.only$bias.WC < 0.350 & WC.exonic.only$bias.WC > 0.050),])
count.WC.pat = nrow(WC.exonic.only[(WC.exonic.only$bias.WC <= 0.050),])
count.WC.total = nrow(WC.exonic.only)

counts.SNP.WC <- data.frame("count" = c(count.WC.mat,count.WC.mat.bias,count.WC.bip,count.WC.pat.bias,count.WC.pat),
                             "perc" = c(count.WC.mat*100/count.WC.total,count.WC.mat.bias*100/count.WC.total,
                                        count.WC.bip*100/count.WC.total,count.WC.pat.bias*100/count.WC.total,count.WC.pat*100/count.WC.total),
                             "bias"=c("M","MB","B","PB","P"))

count.CW.mat = nrow(WC.exonic.only[(WC.exonic.only$bias.CW >= 0.950),])
count.CW.mat.bias = nrow(WC.exonic.only[(WC.exonic.only$bias.CW > 0.650 & WC.exonic.only$bias.CW < 0.950),])
count.CW.bip = nrow(WC.exonic.only[(WC.exonic.only$bias.CW <= 0.650 & WC.exonic.only$bias.CW >= 0.350),])
count.CW.pat.bias = nrow(WC.exonic.only[(WC.exonic.only$bias.CW < 0.350 & WC.exonic.only$bias.CW > 0.050),])
count.CW.pat = nrow(WC.exonic.only[(WC.exonic.only$bias.CW <= 0.050),])
count.CW.total = nrow(WC.exonic.only)

counts.SNP.CW <- data.frame("count" = c(count.CW.mat,count.CW.mat.bias,count.CW.bip,count.CW.pat.bias,count.CW.pat),
                            "perc" = c(count.CW.mat*100/count.CW.total,count.CW.mat.bias*100/count.CW.total,
                                       count.CW.bip*100/count.CW.total,count.CW.pat.bias*100/count.CW.total,count.CW.pat*100/count.CW.total),
                            "bias"=c("M","MB","B","PB","P"))

count.WB.mat = nrow(WB.exonic.only[(WB.exonic.only$bias.WB >= 0.950),])
count.WB.mat.bias = nrow(WB.exonic.only[(WB.exonic.only$bias.WB > 0.650 & WB.exonic.only$bias.WB < 0.950),])
count.WB.bip = nrow(WB.exonic.only[(WB.exonic.only$bias.WB <= 0.650 & WB.exonic.only$bias.WB >= 0.350),])
count.WB.pat.bias = nrow(WB.exonic.only[(WB.exonic.only$bias.WB < 0.350 & WB.exonic.only$bias.WB > 0.050),])
count.WB.pat = nrow(WB.exonic.only[(WB.exonic.only$bias.WB <= 0.050),])
count.WB.total = nrow(WB.exonic.only)

counts.SNP.WB <- data.frame("count" = c(count.WB.mat,count.WB.mat.bias,count.WB.bip,count.WB.pat.bias,count.WB.pat),
                            "perc" = c(count.WB.mat*100/count.WB.total,count.WB.mat.bias*100/count.WB.total,
                                       count.WB.bip*100/count.WB.total,count.WB.pat.bias*100/count.WB.total,count.WB.pat*100/count.WB.total),
                            "bias"=c("M","MB","B","PB","P"))

count.BW.mat = nrow(WB.exonic.only[(WB.exonic.only$bias.BW >= 0.950),])
count.BW.mat.bias = nrow(WB.exonic.only[(WB.exonic.only$bias.BW > 0.650 & WB.exonic.only$bias.BW < 0.95),])
count.BW.bip = nrow(WB.exonic.only[(WB.exonic.only$bias.BW <= 0.650 & WB.exonic.only$bias.BW >= 0.350),])
count.BW.pat.bias = nrow(WB.exonic.only[(WB.exonic.only$bias.BW < 0.350 & WB.exonic.only$bias.BW > 0.050),])
count.BW.pat = nrow(WB.exonic.only[(WB.exonic.only$bias.BW <= 0.050),])
count.BW.total = nrow(WB.exonic.only)

counts.SNP.BW <- data.frame("count" = c(count.BW.mat,count.BW.mat.bias,count.BW.bip,count.BW.pat.bias,count.BW.pat),
                            "perc" = c(count.BW.mat*100/count.BW.total,count.BW.mat.bias*100/count.BW.total,
                                       count.BW.bip*100/count.BW.total,count.BW.pat.bias*100/count.BW.total,count.BW.pat*100/count.BW.total),
                            "bias"=c("M","MB","B","PB","P"))

count.CB.mat = nrow(CB.exonic.only[(CB.exonic.only$bias.CB >= 0.950),])
count.CB.mat.bias = nrow(CB.exonic.only[(CB.exonic.only$bias.CB > 0.650 & CB.exonic.only$bias.CB < 0.95),])
count.CB.bip = nrow(CB.exonic.only[(CB.exonic.only$bias.CB <= 0.650 & CB.exonic.only$bias.CB >= 0.350),])
count.CB.pat.bias = nrow(CB.exonic.only[(CB.exonic.only$bias.CB < 0.350 & CB.exonic.only$bias.CB > 0.050),])
count.CB.pat = nrow(CB.exonic.only[(CB.exonic.only$bias.CB <= 0.050),])
count.CB.total = nrow(CB.exonic.only)

counts.SNP.CB <- data.frame("count" = c(count.CB.mat,count.CB.mat.bias,count.CB.bip,count.CB.pat.bias,count.CB.pat),
                            "perc" = c(count.CB.mat*100/count.CB.total,count.CB.mat.bias*100/count.CB.total,
                                       count.CB.bip*100/count.CB.total,count.CB.pat.bias*100/count.CB.total,count.CB.pat*100/count.CB.total),
                            "bias"=c("M","MB","B","PB","P"))

count.BC.mat = nrow(CB.exonic.only[(CB.exonic.only$bias.BC >= 0.950),])
count.BC.mat.bias = nrow(CB.exonic.only[(CB.exonic.only$bias.BC > 0.650 & CB.exonic.only$bias.BC < 0.95),])
count.BC.bip = nrow(CB.exonic.only[(CB.exonic.only$bias.BC <= 0.650 & CB.exonic.only$bias.BC >= 0.350),])
count.BC.pat.bias = nrow(CB.exonic.only[(CB.exonic.only$bias.BC < 0.350 & CB.exonic.only$bias.BC > 0.050),])
count.BC.pat = nrow(CB.exonic.only[(CB.exonic.only$bias.BC <= 0.050),])
count.BC.total = nrow(CB.exonic.only)

counts.SNP.BC <- data.frame("count" = c(count.BC.mat,count.BC.mat.bias,count.BC.bip,count.BC.pat.bias,count.BC.pat),
                            "perc" = c(count.BC.mat*100/count.BC.total,count.BC.mat.bias*100/count.BC.total,
                                       count.BC.bip*100/count.BC.total,count.BC.pat.bias*100/count.BC.total,count.BC.pat*100/count.BC.total),
                            "bias"=c("M","MB","B","PB","P"))

# get gene counts

WC.anno.all.unique.exonic <- WC.anno.all.unique[WC.anno.all.unique$type == "CDS", ] # should be equal to exonic.sites.WC
WB.anno.all.unique.exonic <- WB.anno.all.unique[WB.anno.all.unique$type == "CDS", ] # should be equal to exonic.sites.WB
CB.anno.all.unique.exonic <- CB.anno.all.unique[CB.anno.all.unique$type == "CDS", ] # should be equal to exonic.sites.CB

WC.anno.all.unique.exonic.for.merge = WC.anno.all.unique.exonic[c(6,3)]
WB.anno.all.unique.exonic.for.merge = WB.anno.all.unique.exonic[c(6,3)]
CB.anno.all.unique.exonic.for.merge = CB.anno.all.unique.exonic[c(6,3)]

WC.anno.all.unique.exonic.for.merge$gene = sub('\\..*', '', WC.anno.all.unique.exonic.for.merge$anno)
WB.anno.all.unique.exonic.for.merge$gene = sub('\\..*', '', WB.anno.all.unique.exonic.for.merge$anno)
CB.anno.all.unique.exonic.for.merge$gene = sub('\\..*', '', CB.anno.all.unique.exonic.for.merge$anno)

WC.exonic.for.merge <- WC.anno.all.unique.exonic.for.merge[WC.anno.all.unique.exonic.for.merge$variant %in% WC.exonic.only$variant,]
WB.exonic.for.merge <- WB.anno.all.unique.exonic.for.merge[WB.anno.all.unique.exonic.for.merge$variant %in% WB.exonic.only$variant,]
CB.exonic.for.merge <- CB.anno.all.unique.exonic.for.merge[CB.anno.all.unique.exonic.for.merge$variant %in% CB.exonic.only$variant,]

length(unique(WC.exonic.for.merge$variant))
length(unique(WB.exonic.for.merge$variant))
length(unique(CB.exonic.for.merge$variant))

WC.gene.for.merge=ddply(WC.exonic.for.merge,c("variant","gene"), summarize, annotation=paste(anno,collapse=","))
WB.gene.for.merge=ddply(WB.exonic.for.merge,c("variant","gene"), summarize, annotation=paste(anno,collapse=","))
CB.gene.for.merge=ddply(CB.exonic.for.merge,c("variant","gene"), summarize, annotation=paste(anno,collapse=","))

nrow(WC.gene.for.merge)
nrow(WC.exonic.for.merge)
length(unique(WC.exonic.for.merge$variant))
WC.gene.for.merge[duplicated(WC.gene.for.merge$variant),] # inspect duplicated records (if any)

nrow(WB.gene.for.merge)
nrow(WB.exonic.for.merge)
length(unique(WB.exonic.for.merge$variant))
WB.gene.for.merge[duplicated(WB.gene.for.merge$variant),] # inspect duplicated records (if any)

nrow(CB.gene.for.merge)
nrow(CB.exonic.for.merge)
length(unique(CB.exonic.for.merge$variant))
CB.gene.for.merge[duplicated(CB.gene.for.merge$variant),] # inspect duplicated records (if any)

WC.exonic.only.to.merge <- WC.exonic.only[c(1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,39,40)]
WB.exonic.only.to.merge <- WB.exonic.only[c(1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,35,36)]
CB.exonic.only.to.merge <- CB.exonic.only[c(1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,35,36)]

WC.variants.with.genes = full_join(WC.exonic.only.to.merge,WC.gene.for.merge,by="variant")
WB.variants.with.genes = full_join(WB.exonic.only.to.merge,WB.gene.for.merge,by="variant")
CB.variants.with.genes = full_join(CB.exonic.only.to.merge,CB.gene.for.merge,by="variant")

genes.raw.WC = ddply(WC.variants.with.genes, c("gene"), summarise,
                     mat.WC01 = sum(na2zero(mat.count.WC01)), pat.WC01 = sum(na2zero(pat.count.WC01)), total.WC01 = sum(na2zero(total.count.WC01)),
                     mat.WC04 = sum(na2zero(mat.count.WC04)), pat.WC04 = sum(na2zero(pat.count.WC04)), total.WC04 = sum(na2zero(total.count.WC04)),
                     mat.WC06 = sum(na2zero(mat.count.WC06)), pat.WC06 = sum(na2zero(pat.count.WC06)), total.WC06 = sum(na2zero(total.count.WC06)),
                     mat.WC19 = sum(na2zero(mat.count.WC19)), pat.WC19 = sum(na2zero(pat.count.WC19)), total.WC19 = sum(na2zero(total.count.WC19)),
                     mat.CW01 = sum(na2zero(mat.count.CW01)), pat.CW01 = sum(na2zero(pat.count.CW01)), total.CW01 = sum(na2zero(total.count.CW01)),
                     mat.CW04 = sum(na2zero(mat.count.CW04)), pat.CW04 = sum(na2zero(pat.count.CW04)), total.CW04 = sum(na2zero(total.count.CW04)),
                     mat.CW10 = sum(na2zero(mat.count.CW10)), pat.CW10 = sum(na2zero(pat.count.CW10)), total.CW10 = sum(na2zero(total.count.CW10)),
                     mat.CW11 = sum(na2zero(mat.count.CW11)), pat.CW11 = sum(na2zero(pat.count.CW11)), total.CW11 = sum(na2zero(total.count.CW11)),
                     N = length(variant)) # all NAs will be relabelled to 0

genes.raw.WB = ddply(WB.variants.with.genes, c("gene"), summarise,
                     mat.WB02 = sum(na2zero(mat.count.WB02)), pat.WB02 = sum(na2zero(pat.count.WB02)), total.WB02 = sum(na2zero(total.count.WB02)),
                     mat.WB03 = sum(na2zero(mat.count.WB03)), pat.WB03 = sum(na2zero(pat.count.WB03)), total.WB03 = sum(na2zero(total.count.WB03)),
                     mat.WB08 = sum(na2zero(mat.count.WB08)), pat.WB08 = sum(na2zero(pat.count.WB08)), total.WB08 = sum(na2zero(total.count.WB08)),
                     mat.BW02 = sum(na2zero(mat.count.BW02)), pat.BW02 = sum(na2zero(pat.count.BW02)), total.BW02 = sum(na2zero(total.count.BW02)),
                     mat.BW03 = sum(na2zero(mat.count.BW03)), pat.BW03 = sum(na2zero(pat.count.BW03)), total.BW03 = sum(na2zero(total.count.BW03)),
                     mat.BW05 = sum(na2zero(mat.count.BW05)), pat.BW05 = sum(na2zero(pat.count.BW05)), total.BW05 = sum(na2zero(total.count.BW05)),
                     mat.BW12 = sum(na2zero(mat.count.BW12)), pat.BW12 = sum(na2zero(pat.count.BW12)), total.BW12 = sum(na2zero(total.count.BW12)),
                     N = length(variant)) # all NAs will be relabelled to 0

genes.raw.CB = ddply(CB.variants.with.genes, c("gene"), summarise,
                     mat.CB01 = sum(na2zero(mat.count.CB01)), pat.CB01 = sum(na2zero(pat.count.CB01)), total.CB01 = sum(na2zero(total.count.CB01)),
                     mat.CB12 = sum(na2zero(mat.count.CB12)), pat.CB12 = sum(na2zero(pat.count.CB12)), total.CB12 = sum(na2zero(total.count.CB12)),
                     mat.CB20 = sum(na2zero(mat.count.CB20)), pat.CB20 = sum(na2zero(pat.count.CB20)), total.CB20 = sum(na2zero(total.count.CB20)),
                     mat.BC03 = sum(na2zero(mat.count.BC03)), pat.BC03 = sum(na2zero(pat.count.BC03)), total.BC03 = sum(na2zero(total.count.BC03)),
                     mat.BC04 = sum(na2zero(mat.count.BC04)), pat.BC04 = sum(na2zero(pat.count.BC04)), total.BC04 = sum(na2zero(total.count.BC04)),
                     mat.BC13 = sum(na2zero(mat.count.BC13)), pat.BC13 = sum(na2zero(pat.count.BC13)), total.BC13 = sum(na2zero(total.count.BC13)),
                     mat.BC17 = sum(na2zero(mat.count.BC17)), pat.BC17 = sum(na2zero(pat.count.BC17)), total.BC17 = sum(na2zero(total.count.BC17)),
                     N = length(variant)) # all NAs will be relabelled to 0

nrow(genes.raw.WC)
nrow(genes.raw.WB)
nrow(genes.raw.CB)

genes.raw.WC$num.rep.WC = ifelse(genes.raw.WC$total.WC01 == 0,0,1) + ifelse(genes.raw.WC$total.WC04 == 0,0,1) + ifelse(genes.raw.WC$total.WC06 == 0,0,1) + ifelse(genes.raw.WC$total.WC19 == 0,0,1)
genes.raw.WC$num.rep.CW = ifelse(genes.raw.WC$total.CW01 == 0,0,1) + ifelse(genes.raw.WC$total.CW04 == 0,0,1) + ifelse(genes.raw.WC$total.CW10 == 0,0,1) + ifelse(genes.raw.WC$total.CW11 == 0,0,1)
genes.raw.WB$num.rep.WB = ifelse(genes.raw.WB$total.WB02 == 0,0,1) + ifelse(genes.raw.WB$total.WB03 == 0,0,1) + ifelse(genes.raw.WB$total.WB08 == 0,0,1)
genes.raw.WB$num.rep.BW = ifelse(genes.raw.WB$total.BW02 == 0,0,1) + ifelse(genes.raw.WB$total.BW03 == 0,0,1) + ifelse(genes.raw.WB$total.BW05 == 0,0,1) + ifelse(genes.raw.WB$total.BW12 == 0,0,1)
genes.raw.CB$num.rep.CB = ifelse(genes.raw.CB$total.CB01 == 0,0,1) + ifelse(genes.raw.CB$total.CB12 == 0,0,1) + ifelse(genes.raw.CB$total.CB20 == 0,0,1)
genes.raw.CB$num.rep.BC = ifelse(genes.raw.CB$total.BC03 == 0,0,1) + ifelse(genes.raw.CB$total.BC04 == 0,0,1) + ifelse(genes.raw.CB$total.BC13 == 0,0,1) + ifelse(genes.raw.CB$total.BC17 == 0,0,1)

##### Filter genes

# first filter: remove genes covered by a single site unless average read depth across replicates > 100

genes.raw.WC.one.SNP <- genes.raw.WC[genes.raw.WC$N == 1,]
genes.raw.WB.one.SNP <- genes.raw.WB[genes.raw.WB$N == 1,]
genes.raw.CB.one.SNP <- genes.raw.CB[genes.raw.CB$N == 1,]

genes.raw.WC.one.SNP$WC_depth_mean = round((genes.raw.WC.one.SNP$total.WC01 + genes.raw.WC.one.SNP$total.WC04 + genes.raw.WC.one.SNP$total.WC06 + genes.raw.WC.one.SNP$total.WC19) / genes.raw.WC.one.SNP$num.rep.WC,1)
genes.raw.WC.one.SNP$CW_depth_mean = round((genes.raw.WC.one.SNP$total.CW01 + genes.raw.WC.one.SNP$total.CW04 + genes.raw.WC.one.SNP$total.CW10 + genes.raw.WC.one.SNP$total.CW11) /  genes.raw.WC.one.SNP$num.rep.CW,1)
genes.raw.WB.one.SNP$WB_depth_mean = round((genes.raw.WB.one.SNP$total.WB02 + genes.raw.WB.one.SNP$total.WB03 + genes.raw.WB.one.SNP$total.WB08) / genes.raw.WB.one.SNP$num.rep.WB,1)
genes.raw.WB.one.SNP$BW_depth_mean = round((genes.raw.WB.one.SNP$total.BW02 + genes.raw.WB.one.SNP$total.BW03 + genes.raw.WB.one.SNP$total.BW05 + genes.raw.WB.one.SNP$total.BW12) / genes.raw.WB.one.SNP$num.rep.BW,1)
genes.raw.CB.one.SNP$CB_depth_mean = round((genes.raw.CB.one.SNP$total.CB01 + genes.raw.CB.one.SNP$total.CB12 + genes.raw.CB.one.SNP$total.CB20) / genes.raw.CB.one.SNP$num.rep.CB,1)
genes.raw.CB.one.SNP$BC_depth_mean = round((genes.raw.CB.one.SNP$total.BC03 + genes.raw.CB.one.SNP$total.BC04 + genes.raw.CB.one.SNP$total.BC13 + genes.raw.CB.one.SNP$total.BC17) / genes.raw.CB.one.SNP$num.rep.BC,1)

genes.raw.WC.one.SNP.less.than.100 <- genes.raw.WC.one.SNP[(genes.raw.WC.one.SNP$WC_depth_mean < 100) | (genes.raw.WC.one.SNP$CW_depth_mean < 100), ] # filter out if at least one reciprocal genotype has an average read depth below 100
genes.raw.WB.one.SNP.less.than.100 <- genes.raw.WB.one.SNP[(genes.raw.WB.one.SNP$WB_depth_mean < 100) | (genes.raw.WB.one.SNP$BW_depth_mean < 100), ]
genes.raw.CB.one.SNP.less.than.100 <- genes.raw.CB.one.SNP[(genes.raw.CB.one.SNP$CB_depth_mean < 100) | (genes.raw.CB.one.SNP$BC_depth_mean < 100), ]

nrow(genes.raw.WC.one.SNP.less.than.100)
nrow(genes.raw.WB.one.SNP.less.than.100)
nrow(genes.raw.CB.one.SNP.less.than.100)

genes.first.filter.WC <- anti_join(genes.raw.WC,genes.raw.WC.one.SNP.less.than.100,by=c("gene"))
genes.first.filter.WB <- anti_join(genes.raw.WB,genes.raw.WB.one.SNP.less.than.100,by=c("gene"))
genes.first.filter.CB <- anti_join(genes.raw.CB,genes.raw.CB.one.SNP.less.than.100,by=c("gene"))

# second filter: G-test of independence to remove genes with high heterogeneity across replicates

library("DescTools")
library("devtools")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("qvalue")

genes.first.filter.WC.test <- genes.first.filter.WC[c("gene","mat.WC01","pat.WC01","total.WC01","mat.WC04","pat.WC04","total.WC04","mat.WC06","pat.WC06","total.WC06","mat.WC19","pat.WC19","total.WC19","num.rep.WC")]
genes.first.filter.CW.test <- genes.first.filter.WC[c("gene","mat.CW01","pat.CW01","total.CW01","mat.CW04","pat.CW04","total.CW04","mat.CW10","pat.CW10","total.CW10","mat.CW11","pat.CW11","total.CW11","num.rep.CW")]
genes.first.filter.WB.test <- genes.first.filter.WB[c("gene","mat.WB02","pat.WB02","total.WB02","mat.WB03","pat.WB03","total.WB03","mat.WB08","pat.WB08","total.WB08","num.rep.WB")]
genes.first.filter.BW.test <- genes.first.filter.WB[c("gene","mat.BW02","pat.BW02","total.BW02","mat.BW03","pat.BW03","total.BW03","mat.BW05","pat.BW05","total.BW05","mat.BW12","pat.BW12","total.BW12","num.rep.BW")]
genes.first.filter.CB.test <- genes.first.filter.CB[c("gene","mat.CB01","pat.CB01","total.CB01","mat.CB12","pat.CB12","total.CB12","mat.CB20","pat.CB20","total.CB20","num.rep.CB")]
genes.first.filter.BC.test <- genes.first.filter.CB[c("gene","mat.BC03","pat.BC03","total.BC03","mat.BC04","pat.BC04","total.BC04","mat.BC13","pat.BC13","total.BC13","mat.BC17","pat.BC17","total.BC17","num.rep.BC")]

genes.first.filter.WC.test$G.WC = 0
genes.first.filter.CW.test$G.CW = 0
genes.first.filter.WB.test$G.WB = 0
genes.first.filter.BW.test$G.BW = 0
genes.first.filter.CB.test$G.CB = 0
genes.first.filter.BC.test$G.BC = 0

genes.first.filter.WC.test$p.value.WC = 0
genes.first.filter.CW.test$p.value.CW = 0
genes.first.filter.WB.test$p.value.WB = 0
genes.first.filter.BW.test$p.value.BW = 0
genes.first.filter.CB.test$p.value.CB = 0
genes.first.filter.BC.test$p.value.BC = 0

### WYE x CP

# WC

# split by number of replicates

genes.first.filter.WC.test.no.miss <- genes.first.filter.WC.test[genes.first.filter.WC.test$num.rep.WC == 4, ]
genes.first.filter.WC.test.no.WC01 <- genes.first.filter.WC.test[genes.first.filter.WC.test$total.WC01 == 0, ]
genes.first.filter.WC.test.no.WC04 <- genes.first.filter.WC.test[genes.first.filter.WC.test$total.WC04 == 0, ]
genes.first.filter.WC.test.no.WC06 <- genes.first.filter.WC.test[genes.first.filter.WC.test$total.WC06 == 0, ]
genes.first.filter.WC.test.no.WC19 <- genes.first.filter.WC.test[genes.first.filter.WC.test$total.WC19 == 0, ]

# no missing replicates

# G test

for(each_gene in 1:nrow(genes.first.filter.WC.test.no.miss)){
  stat = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC19),
                              c(genes.first.filter.WC.test.no.miss[each_gene,]$pat.WC01,genes.first.filter.WC.test.no.miss[each_gene,]$pat.WC04,genes.first.filter.WC.test.no.miss[each_gene,]$pat.WC06,genes.first.filter.WC.test.no.miss[each_gene,]$pat.WC19))))$statistic
  genes.first.filter.WC.test.no.miss$G.WC[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.WC.test.no.miss)){
  p = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC19),
                              c(genes.first.filter.WC.test.no.miss[each_gene,]$pat.WC01,genes.first.filter.WC.test.no.miss[each_gene,]$pat.WC04,genes.first.filter.WC.test.no.miss[each_gene,]$pat.WC06,genes.first.filter.WC.test.no.miss[each_gene,]$pat.WC19))))$p.value
 genes.first.filter.WC.test.no.miss$p.value.WC[each_gene] = p
}

# individual exact binomial test

genes.first.filter.WC.test.no.miss$WC01.p.value = 0
genes.first.filter.WC.test.no.miss$WC04.p.value = 0
genes.first.filter.WC.test.no.miss$WC06.p.value = 0
genes.first.filter.WC.test.no.miss$WC19.p.value = 0

for(each_gene in 1:nrow(genes.first.filter.WC.test.no.miss)){
  a = binom.test(genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.miss[each_gene,]$total.WC01,0.5,alternative="two.sided")$p.value
  genes.first.filter.WC.test.no.miss$WC01.p.value[each_gene] = a
}
for(each_gene in 1:nrow(genes.first.filter.WC.test.no.miss)){
  b = binom.test(genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.miss[each_gene,]$total.WC04,0.5,alternative="two.sided")$p.value
  genes.first.filter.WC.test.no.miss$WC04.p.value[each_gene] = b
}
for(each_gene in 1:nrow(genes.first.filter.WC.test.no.miss)){
  c = binom.test(genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.miss[each_gene,]$total.WC06,0.5,alternative="two.sided")$p.value
  genes.first.filter.WC.test.no.miss$WC06.p.value[each_gene] = c
}
for(each_gene in 1:nrow(genes.first.filter.WC.test.no.miss)){
  d = binom.test(genes.first.filter.WC.test.no.miss[each_gene,]$mat.WC19,genes.first.filter.WC.test.no.miss[each_gene,]$total.WC19,0.5,alternative="two.sided")$p.value
  genes.first.filter.WC.test.no.miss$WC19.p.value[each_gene] = d
}

# individual category of bias

genes.first.filter.WC.test.no.miss$bias.WC01 = genes.first.filter.WC.test.no.miss$mat.WC01 / genes.first.filter.WC.test.no.miss$total.WC01
genes.first.filter.WC.test.no.miss$bias.WC04 = genes.first.filter.WC.test.no.miss$mat.WC04 / genes.first.filter.WC.test.no.miss$total.WC04
genes.first.filter.WC.test.no.miss$bias.WC06 = genes.first.filter.WC.test.no.miss$mat.WC06 / genes.first.filter.WC.test.no.miss$total.WC06
genes.first.filter.WC.test.no.miss$bias.WC19 = genes.first.filter.WC.test.no.miss$mat.WC19 / genes.first.filter.WC.test.no.miss$total.WC19

genes.first.filter.WC.test.no.miss$bias.cat.WC01 = "B"
genes.first.filter.WC.test.no.miss$bias.cat.WC04 = "B"
genes.first.filter.WC.test.no.miss$bias.cat.WC06 = "B"
genes.first.filter.WC.test.no.miss$bias.cat.WC19 = "B"

genes.first.filter.WC.test.no.miss$bias.cat.WC01[genes.first.filter.WC.test.no.miss$bias.WC01<=0.050] <-"P"
genes.first.filter.WC.test.no.miss$bias.cat.WC01[genes.first.filter.WC.test.no.miss$bias.WC01>0.050 & genes.first.filter.WC.test.no.miss$bias.WC01<0.350]<-"PB"
genes.first.filter.WC.test.no.miss$bias.cat.WC01[genes.first.filter.WC.test.no.miss$bias.WC01>=0.950]<-"M"
genes.first.filter.WC.test.no.miss$bias.cat.WC01[genes.first.filter.WC.test.no.miss$bias.WC01>0.650 & genes.first.filter.WC.test.no.miss$bias.WC01<0.950]<-"MB"

genes.first.filter.WC.test.no.miss$bias.cat.WC04[genes.first.filter.WC.test.no.miss$bias.WC04<=0.050] <-"P"
genes.first.filter.WC.test.no.miss$bias.cat.WC04[genes.first.filter.WC.test.no.miss$bias.WC04>0.050 & genes.first.filter.WC.test.no.miss$bias.WC04<0.350]<-"PB"
genes.first.filter.WC.test.no.miss$bias.cat.WC04[genes.first.filter.WC.test.no.miss$bias.WC04>=0.950]<-"M"
genes.first.filter.WC.test.no.miss$bias.cat.WC04[genes.first.filter.WC.test.no.miss$bias.WC04>0.650 & genes.first.filter.WC.test.no.miss$bias.WC04<0.950]<-"MB"

genes.first.filter.WC.test.no.miss$bias.cat.WC06[genes.first.filter.WC.test.no.miss$bias.WC06<=0.050] <-"P"
genes.first.filter.WC.test.no.miss$bias.cat.WC06[genes.first.filter.WC.test.no.miss$bias.WC06>0.050 & genes.first.filter.WC.test.no.miss$bias.WC06<0.350]<-"PB"
genes.first.filter.WC.test.no.miss$bias.cat.WC06[genes.first.filter.WC.test.no.miss$bias.WC06>=0.950]<-"M"
genes.first.filter.WC.test.no.miss$bias.cat.WC06[genes.first.filter.WC.test.no.miss$bias.WC06>0.650 & genes.first.filter.WC.test.no.miss$bias.WC06<0.950]<-"MB"

genes.first.filter.WC.test.no.miss$bias.cat.WC19[genes.first.filter.WC.test.no.miss$bias.WC19<=0.050] <-"P"
genes.first.filter.WC.test.no.miss$bias.cat.WC19[genes.first.filter.WC.test.no.miss$bias.WC19>0.050 & genes.first.filter.WC.test.no.miss$bias.WC19<0.350]<-"PB"
genes.first.filter.WC.test.no.miss$bias.cat.WC19[genes.first.filter.WC.test.no.miss$bias.WC19>=0.950]<-"M"
genes.first.filter.WC.test.no.miss$bias.cat.WC19[genes.first.filter.WC.test.no.miss$bias.WC19>0.650 & genes.first.filter.WC.test.no.miss$bias.WC19<0.950]<-"MB"

# decide

genes.first.filter.WC.test.no.miss$pass.G = ifelse(genes.first.filter.WC.test.no.miss$p.value.WC < 0.05/nrow(genes.first.filter.WC.test),"FAIL","PASS")
genes.first.filter.WC.test.no.miss$all.p.significant = ifelse((genes.first.filter.WC.test.no.miss$WC01.p.value < 0.05/nrow(genes.first.filter.WC.test) &
                                                                  genes.first.filter.WC.test.no.miss$WC04.p.value < 0.05/nrow(genes.first.filter.WC.test) &
                                                                     genes.first.filter.WC.test.no.miss$WC06.p.value < 0.05/nrow(genes.first.filter.WC.test) &
                                                                        genes.first.filter.WC.test.no.miss$WC19.p.value < 0.05/nrow(genes.first.filter.WC.test)) |
                                                                        (genes.first.filter.WC.test.no.miss$WC01.p.value >= 0.05/nrow(genes.first.filter.WC.test) &
                                                                            genes.first.filter.WC.test.no.miss$WC04.p.value >= 0.05/nrow(genes.first.filter.WC.test) &
                                                                               genes.first.filter.WC.test.no.miss$WC06.p.value >= 0.05/nrow(genes.first.filter.WC.test) &
                                                                                  genes.first.filter.WC.test.no.miss$WC19.p.value >= 0.05/nrow(genes.first.filter.WC.test)),
                                                                                "AGREE","DISAGREE")

genes.first.filter.WC.test.no.miss$all.bias.cat = ifelse((genes.first.filter.WC.test.no.miss$bias.cat.WC01 == "P" & genes.first.filter.WC.test.no.miss$bias.cat.WC04 == "P" & genes.first.filter.WC.test.no.miss$bias.cat.WC06 == "P" & genes.first.filter.WC.test.no.miss$bias.cat.WC19 == "P") |
                                                    (genes.first.filter.WC.test.no.miss$bias.cat.WC01 == "PB" & genes.first.filter.WC.test.no.miss$bias.cat.WC04 == "PB" & genes.first.filter.WC.test.no.miss$bias.cat.WC06 == "PB" & genes.first.filter.WC.test.no.miss$bias.cat.WC19 == "PB") |
                                                    (genes.first.filter.WC.test.no.miss$bias.cat.WC01 == "B" & genes.first.filter.WC.test.no.miss$bias.cat.WC04 == "B" & genes.first.filter.WC.test.no.miss$bias.cat.WC06 == "B" & genes.first.filter.WC.test.no.miss$bias.cat.WC19 == "B") |
                                                    (genes.first.filter.WC.test.no.miss$bias.cat.WC01 == "MB" & genes.first.filter.WC.test.no.miss$bias.cat.WC04 == "MB" & genes.first.filter.WC.test.no.miss$bias.cat.WC06 == "MB" & genes.first.filter.WC.test.no.miss$bias.cat.WC19 == "MB") |
                                                    (genes.first.filter.WC.test.no.miss$bias.cat.WC01 == "M" & genes.first.filter.WC.test.no.miss$bias.cat.WC04 == "M" & genes.first.filter.WC.test.no.miss$bias.cat.WC06 == "M" & genes.first.filter.WC.test.no.miss$bias.cat.WC19 == "M"),"AGREE","DISAGREE")

genes.first.filter.WC.test.no.miss.fail.filter = genes.first.filter.WC.test.no.miss[genes.first.filter.WC.test.no.miss$pass.G == "FAIL" & (genes.first.filter.WC.test.no.miss$all.p.significant == "DISAGREE" | genes.first.filter.WC.test.no.miss$all.bias.cat == "DISAGREE"),]
genes.first.filter.WC.test.no.miss.pass.filter = anti_join(genes.first.filter.WC.test.no.miss,genes.first.filter.WC.test.no.miss.fail.filter,by="gene")

# no WC01

for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC01)){
  stat = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.WC01[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.WC01[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.WC01[each_gene,]$mat.WC19),
                              c(genes.first.filter.WC.test.no.WC01[each_gene,]$pat.WC04,genes.first.filter.WC.test.no.WC01[each_gene,]$pat.WC06,genes.first.filter.WC.test.no.WC01[each_gene,]$pat.WC19))))$statistic
  genes.first.filter.WC.test.no.WC01$G.WC[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC01)){
  p = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.WC01[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.WC01[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.WC01[each_gene,]$mat.WC19),
                           c(genes.first.filter.WC.test.no.WC01[each_gene,]$pat.WC04,genes.first.filter.WC.test.no.WC01[each_gene,]$pat.WC06,genes.first.filter.WC.test.no.WC01[each_gene,]$pat.WC19))))$p.value
  genes.first.filter.WC.test.no.WC01$p.value.WC[each_gene] = p
}

genes.first.filter.WC.test.no.WC01$pass.G = ifelse(genes.first.filter.WC.test.no.WC01$p.value.WC < 0.05/nrow(genes.first.filter.WC.test),"FAIL","PASS")

#genes.first.filter.WC.test.no.WC01$WC01.p.value = NA
#genes.first.filter.WC.test.no.WC01$WC04.p.value = 0
#genes.first.filter.WC.test.no.WC01$WC06.p.value = 0
#genes.first.filter.WC.test.no.WC01$WC19.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC01)){
#  b = binom.test(genes.first.filter.WC.test.no.WC01[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.WC01[each_gene,]$total.WC04,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC01$WC04.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC01)){
#  c = binom.test(genes.first.filter.WC.test.no.WC01[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.WC01[each_gene,]$total.WC06,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC01$WC06.p.value[each_gene] = c
#}
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC01)){
#  d = binom.test(genes.first.filter.WC.test.no.WC01[each_gene,]$mat.WC19,genes.first.filter.WC.test.no.WC01[each_gene,]$total.WC19,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC01$WC19.p.value[each_gene] = d
#}

# no WC04

for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC04)){
  stat = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.WC04[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.WC04[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.WC04[each_gene,]$mat.WC19),
                              c(genes.first.filter.WC.test.no.WC04[each_gene,]$pat.WC01,genes.first.filter.WC.test.no.WC04[each_gene,]$pat.WC06,genes.first.filter.WC.test.no.WC04[each_gene,]$pat.WC19))))$statistic
  genes.first.filter.WC.test.no.WC04$G.WC[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC04)){
  p = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.WC04[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.WC04[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.WC04[each_gene,]$mat.WC19),
                           c(genes.first.filter.WC.test.no.WC04[each_gene,]$pat.WC01,genes.first.filter.WC.test.no.WC04[each_gene,]$pat.WC06,genes.first.filter.WC.test.no.WC04[each_gene,]$pat.WC19))))$p.value
  genes.first.filter.WC.test.no.WC04$p.value.WC[each_gene] = p
}

genes.first.filter.WC.test.no.WC04$pass.G = ifelse(genes.first.filter.WC.test.no.WC04$p.value.WC < 0.05/nrow(genes.first.filter.WC.test),"FAIL","PASS")

#genes.first.filter.WC.test.no.WC04$WC01.p.value = 0
#genes.first.filter.WC.test.no.WC04$WC04.p.value = NA
#genes.first.filter.WC.test.no.WC04$WC06.p.value = 0
#genes.first.filter.WC.test.no.WC04$WC19.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC04)){
#  a = binom.test(genes.first.filter.WC.test.no.WC04[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.WC04[each_gene,]$total.WC01,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC04$WC01.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC04)){
#  c = binom.test(genes.first.filter.WC.test.no.WC04[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.WC04[each_gene,]$total.WC06,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC04$WC06.p.value[each_gene] = c
#}
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC04)){
#  d = binom.test(genes.first.filter.WC.test.no.WC04[each_gene,]$mat.WC19,genes.first.filter.WC.test.no.WC04[each_gene,]$total.WC19,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC04$WC19.p.value[each_gene] = d
#}

# no WC06

for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC06)){
  stat = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.WC06[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.WC06[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.WC06[each_gene,]$mat.WC19),
                              c(genes.first.filter.WC.test.no.WC06[each_gene,]$pat.WC01,genes.first.filter.WC.test.no.WC06[each_gene,]$pat.WC04,genes.first.filter.WC.test.no.WC06[each_gene,]$pat.WC19))))$statistic
  genes.first.filter.WC.test.no.WC06$G.WC[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC06)){
  p = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.WC06[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.WC06[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.WC06[each_gene,]$mat.WC19),
                           c(genes.first.filter.WC.test.no.WC06[each_gene,]$pat.WC01,genes.first.filter.WC.test.no.WC06[each_gene,]$pat.WC04,genes.first.filter.WC.test.no.WC06[each_gene,]$pat.WC19))))$p.value
  genes.first.filter.WC.test.no.WC06$p.value.WC[each_gene] = p
}

genes.first.filter.WC.test.no.WC06$pass.G = ifelse(genes.first.filter.WC.test.no.WC06$p.value.WC < 0.05/nrow(genes.first.filter.WC.test),"FAIL","PASS")

#genes.first.filter.WC.test.no.WC06$WC01.p.value = 0
#genes.first.filter.WC.test.no.WC06$WC04.p.value = 0
#genes.first.filter.WC.test.no.WC06$WC06.p.value = NA
#genes.first.filter.WC.test.no.WC06$WC19.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC06)){
#  a = binom.test(genes.first.filter.WC.test.no.WC06[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.WC06[each_gene,]$total.WC01,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC06$WC01.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC06)){
#  b = binom.test(genes.first.filter.WC.test.no.WC06[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.WC06[each_gene,]$total.WC04,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC06$WC04.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC06)){
#  d = binom.test(genes.first.filter.WC.test.no.WC06[each_gene,]$mat.WC19,genes.first.filter.WC.test.no.WC06[each_gene,]$total.WC19,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC06$WC19.p.value[each_gene] = d
#}

# no WC19

for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC19)){
  stat = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.WC19[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.WC19[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.WC19[each_gene,]$mat.WC06),
                              c(genes.first.filter.WC.test.no.WC19[each_gene,]$pat.WC01,genes.first.filter.WC.test.no.WC19[each_gene,]$pat.WC04,genes.first.filter.WC.test.no.WC19[each_gene,]$pat.WC06))))$statistic
  genes.first.filter.WC.test.no.WC19$G.WC[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC19)){
  p = GTest(as.table(rbind(c(genes.first.filter.WC.test.no.WC19[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.WC19[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.WC19[each_gene,]$mat.WC06),
                           c(genes.first.filter.WC.test.no.WC19[each_gene,]$pat.WC01,genes.first.filter.WC.test.no.WC19[each_gene,]$pat.WC04,genes.first.filter.WC.test.no.WC19[each_gene,]$pat.WC06))))$p.value
  genes.first.filter.WC.test.no.WC19$p.value.WC[each_gene] = p
}

genes.first.filter.WC.test.no.WC19$pass.G = ifelse(genes.first.filter.WC.test.no.WC19$p.value.WC < 0.05/nrow(genes.first.filter.WC.test),"FAIL","PASS")

#genes.first.filter.WC.test.no.WC19$WC01.p.value = 0
#genes.first.filter.WC.test.no.WC19$WC04.p.value = 0
#genes.first.filter.WC.test.no.WC19$WC06.p.value = 0
#genes.first.filter.WC.test.no.WC19$WC19.p.value = NA
#
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC19)){
#  a = binom.test(genes.first.filter.WC.test.no.WC19[each_gene,]$mat.WC01,genes.first.filter.WC.test.no.WC19[each_gene,]$total.WC01,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC19$WC01.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC19)){
#  b = binom.test(genes.first.filter.WC.test.no.WC19[each_gene,]$mat.WC04,genes.first.filter.WC.test.no.WC19[each_gene,]$total.WC04,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC19$WC04.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.WC.test.no.WC19)){
#  c = binom.test(genes.first.filter.WC.test.no.WC19[each_gene,]$mat.WC06,genes.first.filter.WC.test.no.WC19[each_gene,]$total.WC06,0.5,alternative="two.sided")$p.value
#  genes.first.filter.WC.test.no.WC19$WC06.p.value[each_gene] = c
#}

### CW

# CW

# split by number of replicates

genes.first.filter.CW.test.no.miss <- genes.first.filter.CW.test[genes.first.filter.CW.test$num.rep.CW == 4, ]
genes.first.filter.CW.test.no.CW01 <- genes.first.filter.CW.test[genes.first.filter.CW.test$total.CW01 == 0, ]
genes.first.filter.CW.test.no.CW04 <- genes.first.filter.CW.test[genes.first.filter.CW.test$total.CW04 == 0, ]
genes.first.filter.CW.test.no.CW10 <- genes.first.filter.CW.test[genes.first.filter.CW.test$total.CW10 == 0, ]
genes.first.filter.CW.test.no.CW11 <- genes.first.filter.CW.test[genes.first.filter.CW.test$total.CW11 == 0, ]

# no missing replicates

# G test

for(each_gene in 1:nrow(genes.first.filter.CW.test.no.miss)){
  stat = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW11),
                              c(genes.first.filter.CW.test.no.miss[each_gene,]$pat.CW01,genes.first.filter.CW.test.no.miss[each_gene,]$pat.CW04,genes.first.filter.CW.test.no.miss[each_gene,]$pat.CW10,genes.first.filter.CW.test.no.miss[each_gene,]$pat.CW11))))$statistic
  genes.first.filter.CW.test.no.miss$G.CW[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.CW.test.no.miss)){
  p = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW11),
                           c(genes.first.filter.CW.test.no.miss[each_gene,]$pat.CW01,genes.first.filter.CW.test.no.miss[each_gene,]$pat.CW04,genes.first.filter.CW.test.no.miss[each_gene,]$pat.CW10,genes.first.filter.CW.test.no.miss[each_gene,]$pat.CW11))))$p.value
  genes.first.filter.CW.test.no.miss$p.value.CW[each_gene] = p
}

# individual exact binomial test

genes.first.filter.CW.test.no.miss$CW01.p.value = 0
genes.first.filter.CW.test.no.miss$CW04.p.value = 0
genes.first.filter.CW.test.no.miss$CW10.p.value = 0
genes.first.filter.CW.test.no.miss$CW11.p.value = 0

for(each_gene in 1:nrow(genes.first.filter.CW.test.no.miss)){
  a = binom.test(genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.miss[each_gene,]$total.CW01,0.5,alternative="two.sided")$p.value
  genes.first.filter.CW.test.no.miss$CW01.p.value[each_gene] = a
}
for(each_gene in 1:nrow(genes.first.filter.CW.test.no.miss)){
  b = binom.test(genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.miss[each_gene,]$total.CW04,0.5,alternative="two.sided")$p.value
  genes.first.filter.CW.test.no.miss$CW04.p.value[each_gene] = b
}
for(each_gene in 1:nrow(genes.first.filter.CW.test.no.miss)){
  c = binom.test(genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.miss[each_gene,]$total.CW10,0.5,alternative="two.sided")$p.value
  genes.first.filter.CW.test.no.miss$CW10.p.value[each_gene] = c
}
for(each_gene in 1:nrow(genes.first.filter.CW.test.no.miss)){
  d = binom.test(genes.first.filter.CW.test.no.miss[each_gene,]$mat.CW11,genes.first.filter.CW.test.no.miss[each_gene,]$total.CW11,0.5,alternative="two.sided")$p.value
  genes.first.filter.CW.test.no.miss$CW11.p.value[each_gene] = d
}

# individual category of bias

genes.first.filter.CW.test.no.miss$bias.CW01 = genes.first.filter.CW.test.no.miss$mat.CW01 / genes.first.filter.CW.test.no.miss$total.CW01
genes.first.filter.CW.test.no.miss$bias.CW04 = genes.first.filter.CW.test.no.miss$mat.CW04 / genes.first.filter.CW.test.no.miss$total.CW04
genes.first.filter.CW.test.no.miss$bias.CW10 = genes.first.filter.CW.test.no.miss$mat.CW10 / genes.first.filter.CW.test.no.miss$total.CW10
genes.first.filter.CW.test.no.miss$bias.CW11 = genes.first.filter.CW.test.no.miss$mat.CW11 / genes.first.filter.CW.test.no.miss$total.CW11

genes.first.filter.CW.test.no.miss$bias.cat.CW01 = "B"
genes.first.filter.CW.test.no.miss$bias.cat.CW04 = "B"
genes.first.filter.CW.test.no.miss$bias.cat.CW10 = "B"
genes.first.filter.CW.test.no.miss$bias.cat.CW11 = "B"

genes.first.filter.CW.test.no.miss$bias.cat.CW01[genes.first.filter.CW.test.no.miss$bias.CW01<=0.050] <-"P"
genes.first.filter.CW.test.no.miss$bias.cat.CW01[genes.first.filter.CW.test.no.miss$bias.CW01>0.050 & genes.first.filter.CW.test.no.miss$bias.CW01<0.350]<-"PB"
genes.first.filter.CW.test.no.miss$bias.cat.CW01[genes.first.filter.CW.test.no.miss$bias.CW01>=0.950]<-"M"
genes.first.filter.CW.test.no.miss$bias.cat.CW01[genes.first.filter.CW.test.no.miss$bias.CW01>0.650 & genes.first.filter.CW.test.no.miss$bias.CW01<0.950]<-"MB"

genes.first.filter.CW.test.no.miss$bias.cat.CW04[genes.first.filter.CW.test.no.miss$bias.CW04<=0.050] <-"P"
genes.first.filter.CW.test.no.miss$bias.cat.CW04[genes.first.filter.CW.test.no.miss$bias.CW04>0.050 & genes.first.filter.CW.test.no.miss$bias.CW04<0.350]<-"PB"
genes.first.filter.CW.test.no.miss$bias.cat.CW04[genes.first.filter.CW.test.no.miss$bias.CW04>=0.950]<-"M"
genes.first.filter.CW.test.no.miss$bias.cat.CW04[genes.first.filter.CW.test.no.miss$bias.CW04>0.650 & genes.first.filter.CW.test.no.miss$bias.CW04<0.950]<-"MB"

genes.first.filter.CW.test.no.miss$bias.cat.CW10[genes.first.filter.CW.test.no.miss$bias.CW10<=0.050] <-"P"
genes.first.filter.CW.test.no.miss$bias.cat.CW10[genes.first.filter.CW.test.no.miss$bias.CW10>0.050 & genes.first.filter.CW.test.no.miss$bias.CW10<0.350]<-"PB"
genes.first.filter.CW.test.no.miss$bias.cat.CW10[genes.first.filter.CW.test.no.miss$bias.CW10>=0.950]<-"M"
genes.first.filter.CW.test.no.miss$bias.cat.CW10[genes.first.filter.CW.test.no.miss$bias.CW10>0.650 & genes.first.filter.CW.test.no.miss$bias.CW10<0.950]<-"MB"

genes.first.filter.CW.test.no.miss$bias.cat.CW11[genes.first.filter.CW.test.no.miss$bias.CW11<=0.050] <-"P"
genes.first.filter.CW.test.no.miss$bias.cat.CW11[genes.first.filter.CW.test.no.miss$bias.CW11>0.050 & genes.first.filter.CW.test.no.miss$bias.CW11<0.350]<-"PB"
genes.first.filter.CW.test.no.miss$bias.cat.CW11[genes.first.filter.CW.test.no.miss$bias.CW11>=0.950]<-"M"
genes.first.filter.CW.test.no.miss$bias.cat.CW11[genes.first.filter.CW.test.no.miss$bias.CW11>0.650 & genes.first.filter.CW.test.no.miss$bias.CW11<0.950]<-"MB"

# decide

genes.first.filter.CW.test.no.miss$pass.G = ifelse(genes.first.filter.CW.test.no.miss$p.value.CW < 0.05/nrow(genes.first.filter.CW.test),"FAIL","PASS")
genes.first.filter.CW.test.no.miss$all.p.significant = ifelse((genes.first.filter.CW.test.no.miss$CW01.p.value < 0.05/nrow(genes.first.filter.CW.test) &
                                                                 genes.first.filter.CW.test.no.miss$CW04.p.value < 0.05/nrow(genes.first.filter.CW.test) &
                                                                 genes.first.filter.CW.test.no.miss$CW10.p.value < 0.05/nrow(genes.first.filter.CW.test) &
                                                                 genes.first.filter.CW.test.no.miss$CW11.p.value < 0.05/nrow(genes.first.filter.CW.test)) |
                                                                (genes.first.filter.CW.test.no.miss$CW01.p.value >= 0.05/nrow(genes.first.filter.CW.test) &
                                                                   genes.first.filter.CW.test.no.miss$CW04.p.value >= 0.05/nrow(genes.first.filter.CW.test) &
                                                                   genes.first.filter.CW.test.no.miss$CW10.p.value >= 0.05/nrow(genes.first.filter.CW.test) &
                                                                   genes.first.filter.CW.test.no.miss$CW11.p.value >= 0.05/nrow(genes.first.filter.CW.test)),
                                                              "AGREE","DISAGREE")

genes.first.filter.CW.test.no.miss$all.bias.cat = ifelse((genes.first.filter.CW.test.no.miss$bias.cat.CW01 == "P" & genes.first.filter.CW.test.no.miss$bias.cat.CW04 == "P" & genes.first.filter.CW.test.no.miss$bias.cat.CW10 == "P" & genes.first.filter.CW.test.no.miss$bias.cat.CW11 == "P") |
                                                           (genes.first.filter.CW.test.no.miss$bias.cat.CW01 == "PB" & genes.first.filter.CW.test.no.miss$bias.cat.CW04 == "PB" & genes.first.filter.CW.test.no.miss$bias.cat.CW10 == "PB" & genes.first.filter.CW.test.no.miss$bias.cat.CW11 == "PB") |
                                                           (genes.first.filter.CW.test.no.miss$bias.cat.CW01 == "B" & genes.first.filter.CW.test.no.miss$bias.cat.CW04 == "B" & genes.first.filter.CW.test.no.miss$bias.cat.CW10 == "B" & genes.first.filter.CW.test.no.miss$bias.cat.CW11 == "B") |
                                                           (genes.first.filter.CW.test.no.miss$bias.cat.CW01 == "MB" & genes.first.filter.CW.test.no.miss$bias.cat.CW04 == "MB" & genes.first.filter.CW.test.no.miss$bias.cat.CW10 == "MB" & genes.first.filter.CW.test.no.miss$bias.cat.CW11 == "MB") |
                                                           (genes.first.filter.CW.test.no.miss$bias.cat.CW01 == "M" & genes.first.filter.CW.test.no.miss$bias.cat.CW04 == "M" & genes.first.filter.CW.test.no.miss$bias.cat.CW10 == "M" & genes.first.filter.CW.test.no.miss$bias.cat.CW11 == "M"),"AGREE","DISAGREE")

genes.first.filter.CW.test.no.miss.fail.filter = genes.first.filter.CW.test.no.miss[genes.first.filter.CW.test.no.miss$pass.G == "FAIL" & (genes.first.filter.CW.test.no.miss$all.p.significant == "DISAGREE" | genes.first.filter.CW.test.no.miss$all.bias.cat == "DISAGREE"),]
genes.first.filter.CW.test.no.miss.pass.filter = anti_join(genes.first.filter.CW.test.no.miss,genes.first.filter.CW.test.no.miss.fail.filter,by="gene")

# no CW01

for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW01)){
  stat = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.CW01[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.CW01[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.CW01[each_gene,]$mat.CW11),
                              c(genes.first.filter.CW.test.no.CW01[each_gene,]$pat.CW04,genes.first.filter.CW.test.no.CW01[each_gene,]$pat.CW10,genes.first.filter.CW.test.no.CW01[each_gene,]$pat.CW11))))$statistic
  genes.first.filter.CW.test.no.CW01$G.CW[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW01)){
  p = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.CW01[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.CW01[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.CW01[each_gene,]$mat.CW11),
                           c(genes.first.filter.CW.test.no.CW01[each_gene,]$pat.CW04,genes.first.filter.CW.test.no.CW01[each_gene,]$pat.CW10,genes.first.filter.CW.test.no.CW01[each_gene,]$pat.CW11))))$p.value
  genes.first.filter.CW.test.no.CW01$p.value.CW[each_gene] = p
}

genes.first.filter.CW.test.no.CW01$pass.G = ifelse(genes.first.filter.CW.test.no.CW01$p.value.CW < 0.05/nrow(genes.first.filter.CW.test),"FAIL","PASS")

#genes.first.filter.CW.test.no.CW01$CW01.p.value = NA
#genes.first.filter.CW.test.no.CW01$CW04.p.value = 0
#genes.first.filter.CW.test.no.CW01$CW10.p.value = 0
#genes.first.filter.CW.test.no.CW01$CW11.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW01)){
#  b = binom.test(genes.first.filter.CW.test.no.CW01[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.CW01[each_gene,]$total.CW04,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW01$CW04.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW01)){
#  c = binom.test(genes.first.filter.CW.test.no.CW01[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.CW01[each_gene,]$total.CW10,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW01$CW10.p.value[each_gene] = c
#}
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW01)){
#  d = binom.test(genes.first.filter.CW.test.no.CW01[each_gene,]$mat.CW11,genes.first.filter.CW.test.no.CW01[each_gene,]$total.CW11,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW01$CW11.p.value[each_gene] = d
#}

# no CW04

for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW04)){
  stat = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.CW04[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.CW04[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.CW04[each_gene,]$mat.CW11),
                              c(genes.first.filter.CW.test.no.CW04[each_gene,]$pat.CW01,genes.first.filter.CW.test.no.CW04[each_gene,]$pat.CW10,genes.first.filter.CW.test.no.CW04[each_gene,]$pat.CW11))))$statistic
  genes.first.filter.CW.test.no.CW04$G.CW[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW04)){
  p = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.CW04[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.CW04[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.CW04[each_gene,]$mat.CW11),
                           c(genes.first.filter.CW.test.no.CW04[each_gene,]$pat.CW01,genes.first.filter.CW.test.no.CW04[each_gene,]$pat.CW10,genes.first.filter.CW.test.no.CW04[each_gene,]$pat.CW11))))$p.value
  genes.first.filter.CW.test.no.CW04$p.value.CW[each_gene] = p
}

genes.first.filter.CW.test.no.CW04$pass.G = ifelse(genes.first.filter.CW.test.no.CW04$p.value.CW < 0.05/nrow(genes.first.filter.CW.test),"FAIL","PASS")

#genes.first.filter.CW.test.no.CW04$CW01.p.value = 0
#genes.first.filter.CW.test.no.CW04$CW04.p.value = NA
#genes.first.filter.CW.test.no.CW04$CW10.p.value = 0
#genes.first.filter.CW.test.no.CW04$CW11.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW04)){
#  a = binom.test(genes.first.filter.CW.test.no.CW04[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.CW04[each_gene,]$total.CW01,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW04$CW01.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW04)){
#  c = binom.test(genes.first.filter.CW.test.no.CW04[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.CW04[each_gene,]$total.CW10,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW04$CW10.p.value[each_gene] = c
#}
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW04)){
#  d = binom.test(genes.first.filter.CW.test.no.CW04[each_gene,]$mat.CW11,genes.first.filter.CW.test.no.CW04[each_gene,]$total.CW11,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW04$CW11.p.value[each_gene] = d
#}

# no CW10

for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW10)){
  stat = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.CW10[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.CW10[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.CW10[each_gene,]$mat.CW11),
                              c(genes.first.filter.CW.test.no.CW10[each_gene,]$pat.CW01,genes.first.filter.CW.test.no.CW10[each_gene,]$pat.CW04,genes.first.filter.CW.test.no.CW10[each_gene,]$pat.CW11))))$statistic
  genes.first.filter.CW.test.no.CW10$G.CW[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW10)){
  p = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.CW10[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.CW10[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.CW10[each_gene,]$mat.CW11),
                           c(genes.first.filter.CW.test.no.CW10[each_gene,]$pat.CW01,genes.first.filter.CW.test.no.CW10[each_gene,]$pat.CW04,genes.first.filter.CW.test.no.CW10[each_gene,]$pat.CW11))))$p.value
  genes.first.filter.CW.test.no.CW10$p.value.CW[each_gene] = p
}

genes.first.filter.CW.test.no.CW10$pass.G = ifelse(genes.first.filter.CW.test.no.CW10$p.value.CW < 0.05/nrow(genes.first.filter.CW.test),"FAIL","PASS")

#genes.first.filter.CW.test.no.CW10$CW01.p.value = 0
#genes.first.filter.CW.test.no.CW10$CW04.p.value = 0
#genes.first.filter.CW.test.no.CW10$CW10.p.value = NA
#genes.first.filter.CW.test.no.CW10$CW11.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW10)){
#  a = binom.test(genes.first.filter.CW.test.no.CW10[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.CW10[each_gene,]$total.CW01,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW10$CW01.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW10)){
#  b = binom.test(genes.first.filter.CW.test.no.CW10[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.CW10[each_gene,]$total.CW04,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW10$CW04.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW10)){
#  d = binom.test(genes.first.filter.CW.test.no.CW10[each_gene,]$mat.CW11,genes.first.filter.CW.test.no.CW10[each_gene,]$total.CW11,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW10$CW11.p.value[each_gene] = d
#}

# no CW11

for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW11)){
  stat = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.CW11[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.CW11[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.CW11[each_gene,]$mat.CW10),
                              c(genes.first.filter.CW.test.no.CW11[each_gene,]$pat.CW01,genes.first.filter.CW.test.no.CW11[each_gene,]$pat.CW04,genes.first.filter.CW.test.no.CW11[each_gene,]$pat.CW10))))$statistic
  genes.first.filter.CW.test.no.CW11$G.CW[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW11)){
  p = GTest(as.table(rbind(c(genes.first.filter.CW.test.no.CW11[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.CW11[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.CW11[each_gene,]$mat.CW10),
                           c(genes.first.filter.CW.test.no.CW11[each_gene,]$pat.CW01,genes.first.filter.CW.test.no.CW11[each_gene,]$pat.CW04,genes.first.filter.CW.test.no.CW11[each_gene,]$pat.CW10))))$p.value
  genes.first.filter.CW.test.no.CW11$p.value.CW[each_gene] = p
}

genes.first.filter.CW.test.no.CW11$pass.G = ifelse(genes.first.filter.CW.test.no.CW11$p.value.CW < 0.05/nrow(genes.first.filter.CW.test),"FAIL","PASS")

#genes.first.filter.CW.test.no.CW11$CW01.p.value = 0
#genes.first.filter.CW.test.no.CW11$CW04.p.value = 0
#genes.first.filter.CW.test.no.CW11$CW10.p.value = 0
#genes.first.filter.CW.test.no.CW11$CW11.p.value = NA
#
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW11)){
#  a = binom.test(genes.first.filter.CW.test.no.CW11[each_gene,]$mat.CW01,genes.first.filter.CW.test.no.CW11[each_gene,]$total.CW01,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW11$CW01.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW11)){
#  b = binom.test(genes.first.filter.CW.test.no.CW11[each_gene,]$mat.CW04,genes.first.filter.CW.test.no.CW11[each_gene,]$total.CW04,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW11$CW04.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.CW.test.no.CW11)){
#  c = binom.test(genes.first.filter.CW.test.no.CW11[each_gene,]$mat.CW10,genes.first.filter.CW.test.no.CW11[each_gene,]$total.CW10,0.5,alternative="two.sided")$p.value
#  genes.first.filter.CW.test.no.CW11$CW10.p.value[each_gene] = c
#}

nrow(genes.first.filter.WC.test.no.miss.fail.filter)
nrow(genes.first.filter.CW.test.no.miss.fail.filter)
genes.second.filter.WC.fail = rbind(genes.first.filter.WC.test.no.miss.fail.filter[c("gene")],genes.first.filter.CW.test.no.miss.fail.filter[c("gene")])
length(unique(genes.second.filter.WC.fail$gene))
nrow(genes.second.filter.WC.fail)
genes.second.filter.WC = anti_join(genes.first.filter.WC,genes.second.filter.WC.fail,by="gene")

### WYE x BGOX

## WB

# split by number of replicates

genes.first.filter.WB.test.no.miss <- genes.first.filter.WB.test[genes.first.filter.WB.test$num.rep.WB == 3, ] # no need to split, only three replicates

# G test

for(each_gene in 1:nrow(genes.first.filter.WB.test.no.miss)){
  stat = GTest(as.table(rbind(c(genes.first.filter.WB.test.no.miss[each_gene,]$mat.WB02,genes.first.filter.WB.test.no.miss[each_gene,]$mat.WB03,genes.first.filter.WB.test.no.miss[each_gene,]$mat.WB08),
                              c(genes.first.filter.WB.test.no.miss[each_gene,]$pat.WB02,genes.first.filter.WB.test.no.miss[each_gene,]$pat.WB03,genes.first.filter.WB.test.no.miss[each_gene,]$pat.WB08))))$statistic
  genes.first.filter.WB.test.no.miss$G.WB[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.WB.test.no.miss)){
  p = GTest(as.table(rbind(c(genes.first.filter.WB.test.no.miss[each_gene,]$mat.WB02,genes.first.filter.WB.test.no.miss[each_gene,]$mat.WB03,genes.first.filter.WB.test.no.miss[each_gene,]$mat.WB08),
                           c(genes.first.filter.WB.test.no.miss[each_gene,]$pat.WB02,genes.first.filter.WB.test.no.miss[each_gene,]$pat.WB03,genes.first.filter.WB.test.no.miss[each_gene,]$pat.WB08))))$p.value
  genes.first.filter.WB.test.no.miss$p.value.WB[each_gene] = p
}

# individual exact binomial test

genes.first.filter.WB.test.no.miss$WB02.p.value = 0
genes.first.filter.WB.test.no.miss$WB03.p.value = 0
genes.first.filter.WB.test.no.miss$WB08.p.value = 0

for(each_gene in 1:nrow(genes.first.filter.WB.test.no.miss)){
  a = binom.test(genes.first.filter.WB.test.no.miss[each_gene,]$mat.WB02,genes.first.filter.WB.test.no.miss[each_gene,]$total.WB02,0.5,alternative="two.sided")$p.value
  genes.first.filter.WB.test.no.miss$WB02.p.value[each_gene] = a
}
for(each_gene in 1:nrow(genes.first.filter.WB.test.no.miss)){
  b = binom.test(genes.first.filter.WB.test.no.miss[each_gene,]$mat.WB03,genes.first.filter.WB.test.no.miss[each_gene,]$total.WB03,0.5,alternative="two.sided")$p.value
  genes.first.filter.WB.test.no.miss$WB03.p.value[each_gene] = b
}
for(each_gene in 1:nrow(genes.first.filter.WB.test.no.miss)){
  c = binom.test(genes.first.filter.WB.test.no.miss[each_gene,]$mat.WB08,genes.first.filter.WB.test.no.miss[each_gene,]$total.WB08,0.5,alternative="two.sided")$p.value
  genes.first.filter.WB.test.no.miss$WB08.p.value[each_gene] = c
}

# individual category of bias

genes.first.filter.WB.test.no.miss$bias.WB02 = genes.first.filter.WB.test.no.miss$mat.WB02 / genes.first.filter.WB.test.no.miss$total.WB02
genes.first.filter.WB.test.no.miss$bias.WB03 = genes.first.filter.WB.test.no.miss$mat.WB03 / genes.first.filter.WB.test.no.miss$total.WB03
genes.first.filter.WB.test.no.miss$bias.WB08 = genes.first.filter.WB.test.no.miss$mat.WB08 / genes.first.filter.WB.test.no.miss$total.WB08

genes.first.filter.WB.test.no.miss$bias.cat.WB02 = "B"
genes.first.filter.WB.test.no.miss$bias.cat.WB03 = "B"
genes.first.filter.WB.test.no.miss$bias.cat.WB08 = "B"

genes.first.filter.WB.test.no.miss$bias.cat.WB02[genes.first.filter.WB.test.no.miss$bias.WB02<=0.050] <-"P"
genes.first.filter.WB.test.no.miss$bias.cat.WB02[genes.first.filter.WB.test.no.miss$bias.WB02>0.050 & genes.first.filter.WB.test.no.miss$bias.WB02<0.350]<-"PB"
genes.first.filter.WB.test.no.miss$bias.cat.WB02[genes.first.filter.WB.test.no.miss$bias.WB02>=0.950]<-"M"
genes.first.filter.WB.test.no.miss$bias.cat.WB02[genes.first.filter.WB.test.no.miss$bias.WB02>0.650 & genes.first.filter.WB.test.no.miss$bias.WB02<0.950]<-"MB"

genes.first.filter.WB.test.no.miss$bias.cat.WB03[genes.first.filter.WB.test.no.miss$bias.WB03<=0.050] <-"P"
genes.first.filter.WB.test.no.miss$bias.cat.WB03[genes.first.filter.WB.test.no.miss$bias.WB03>0.050 & genes.first.filter.WB.test.no.miss$bias.WB03<0.350]<-"PB"
genes.first.filter.WB.test.no.miss$bias.cat.WB03[genes.first.filter.WB.test.no.miss$bias.WB03>=0.950]<-"M"
genes.first.filter.WB.test.no.miss$bias.cat.WB03[genes.first.filter.WB.test.no.miss$bias.WB03>0.650 & genes.first.filter.WB.test.no.miss$bias.WB03<0.950]<-"MB"

genes.first.filter.WB.test.no.miss$bias.cat.WB08[genes.first.filter.WB.test.no.miss$bias.WB08<=0.050] <-"P"
genes.first.filter.WB.test.no.miss$bias.cat.WB08[genes.first.filter.WB.test.no.miss$bias.WB08>0.050 & genes.first.filter.WB.test.no.miss$bias.WB08<0.350]<-"PB"
genes.first.filter.WB.test.no.miss$bias.cat.WB08[genes.first.filter.WB.test.no.miss$bias.WB08>=0.950]<-"M"
genes.first.filter.WB.test.no.miss$bias.cat.WB08[genes.first.filter.WB.test.no.miss$bias.WB08>0.650 & genes.first.filter.WB.test.no.miss$bias.WB08<0.950]<-"MB"

# decide

genes.first.filter.WB.test.no.miss$pass.G = ifelse(genes.first.filter.WB.test.no.miss$p.value.WB < 0.05/nrow(genes.first.filter.WB.test),"FAIL","PASS")
genes.first.filter.WB.test.no.miss$all.p.significant = ifelse((genes.first.filter.WB.test.no.miss$WB02.p.value < 0.05/nrow(genes.first.filter.WB.test) &
                                                                 genes.first.filter.WB.test.no.miss$WB03.p.value < 0.05/nrow(genes.first.filter.WB.test) &
                                                                 genes.first.filter.WB.test.no.miss$WB08.p.value < 0.05/nrow(genes.first.filter.WB.test)) |
                                                                (genes.first.filter.WB.test.no.miss$WB02.p.value >= 0.05/nrow(genes.first.filter.WB.test) &
                                                                   genes.first.filter.WB.test.no.miss$WB03.p.value >= 0.05/nrow(genes.first.filter.WB.test) &
                                                                   genes.first.filter.WB.test.no.miss$WB08.p.value >= 0.05/nrow(genes.first.filter.WB.test)),"AGREE","DISAGREE")

genes.first.filter.WB.test.no.miss$all.bias.cat = ifelse((genes.first.filter.WB.test.no.miss$bias.cat.WB02 == "P" & genes.first.filter.WB.test.no.miss$bias.cat.WB03 == "P" & genes.first.filter.WB.test.no.miss$bias.cat.WB08 == "P") |
                                                           (genes.first.filter.WB.test.no.miss$bias.cat.WB02 == "PB" & genes.first.filter.WB.test.no.miss$bias.cat.WB03 == "PB" & genes.first.filter.WB.test.no.miss$bias.cat.WB08 == "PB") |
                                                           (genes.first.filter.WB.test.no.miss$bias.cat.WB02 == "B" & genes.first.filter.WB.test.no.miss$bias.cat.WB03 == "B" & genes.first.filter.WB.test.no.miss$bias.cat.WB08 == "B") |
                                                           (genes.first.filter.WB.test.no.miss$bias.cat.WB02 == "MB" & genes.first.filter.WB.test.no.miss$bias.cat.WB03 == "MB" & genes.first.filter.WB.test.no.miss$bias.cat.WB08 == "MB") |
                                                           (genes.first.filter.WB.test.no.miss$bias.cat.WB02 == "M" & genes.first.filter.WB.test.no.miss$bias.cat.WB03 == "M" & genes.first.filter.WB.test.no.miss$bias.cat.WB08 == "M"),"AGREE","DISAGREE")

genes.first.filter.WB.test.no.miss.fail.filter = genes.first.filter.WB.test.no.miss[genes.first.filter.WB.test.no.miss$pass.G == "FAIL" & (genes.first.filter.WB.test.no.miss$all.p.significant == "DISAGREE" | genes.first.filter.WB.test.no.miss$all.bias.cat == "DISAGREE"),]
genes.first.filter.WB.test.no.miss.pass.filter = anti_join(genes.first.filter.WB.test.no.miss,genes.first.filter.WB.test.no.miss.fail.filter,by="gene")

### BW

# BW

# split by number of replicates

genes.first.filter.BW.test.no.miss <- genes.first.filter.BW.test[genes.first.filter.BW.test$num.rep.BW == 4, ]
genes.first.filter.BW.test.no.BW02 <- genes.first.filter.BW.test[genes.first.filter.BW.test$total.BW02 == 0, ]
genes.first.filter.BW.test.no.BW03 <- genes.first.filter.BW.test[genes.first.filter.BW.test$total.BW03 == 0, ]
genes.first.filter.BW.test.no.BW05 <- genes.first.filter.BW.test[genes.first.filter.BW.test$total.BW05 == 0, ]
genes.first.filter.BW.test.no.BW12 <- genes.first.filter.BW.test[genes.first.filter.BW.test$total.BW12 == 0, ]

nrow(genes.first.filter.BW.test.no.miss)
nrow(genes.first.filter.BW.test.no.BW02)
nrow(genes.first.filter.BW.test.no.BW03)
nrow(genes.first.filter.BW.test.no.BW05)
nrow(genes.first.filter.BW.test.no.BW12)

# no missing replicates

# G test

for(each_gene in 1:nrow(genes.first.filter.BW.test.no.miss)){
  stat = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW12),
                              c(genes.first.filter.BW.test.no.miss[each_gene,]$pat.BW02,genes.first.filter.BW.test.no.miss[each_gene,]$pat.BW03,genes.first.filter.BW.test.no.miss[each_gene,]$pat.BW05,genes.first.filter.BW.test.no.miss[each_gene,]$pat.BW12))))$statistic
  genes.first.filter.BW.test.no.miss$G.BW[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.BW.test.no.miss)){
  p = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW12),
                           c(genes.first.filter.BW.test.no.miss[each_gene,]$pat.BW02,genes.first.filter.BW.test.no.miss[each_gene,]$pat.BW03,genes.first.filter.BW.test.no.miss[each_gene,]$pat.BW05,genes.first.filter.BW.test.no.miss[each_gene,]$pat.BW12))))$p.value
  genes.first.filter.BW.test.no.miss$p.value.BW[each_gene] = p
}

# individual exact binomial test

genes.first.filter.BW.test.no.miss$BW02.p.value = 0
genes.first.filter.BW.test.no.miss$BW03.p.value = 0
genes.first.filter.BW.test.no.miss$BW05.p.value = 0
genes.first.filter.BW.test.no.miss$BW12.p.value = 0

for(each_gene in 1:nrow(genes.first.filter.BW.test.no.miss)){
  a = binom.test(genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.miss[each_gene,]$total.BW02,0.5,alternative="two.sided")$p.value
  genes.first.filter.BW.test.no.miss$BW02.p.value[each_gene] = a
}
for(each_gene in 1:nrow(genes.first.filter.BW.test.no.miss)){
  b = binom.test(genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.miss[each_gene,]$total.BW03,0.5,alternative="two.sided")$p.value
  genes.first.filter.BW.test.no.miss$BW03.p.value[each_gene] = b
}
for(each_gene in 1:nrow(genes.first.filter.BW.test.no.miss)){
  c = binom.test(genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.miss[each_gene,]$total.BW05,0.5,alternative="two.sided")$p.value
  genes.first.filter.BW.test.no.miss$BW05.p.value[each_gene] = c
}
for(each_gene in 1:nrow(genes.first.filter.BW.test.no.miss)){
  d = binom.test(genes.first.filter.BW.test.no.miss[each_gene,]$mat.BW12,genes.first.filter.BW.test.no.miss[each_gene,]$total.BW12,0.5,alternative="two.sided")$p.value
  genes.first.filter.BW.test.no.miss$BW12.p.value[each_gene] = d
}

# individual category of bias

genes.first.filter.BW.test.no.miss$bias.BW02 = genes.first.filter.BW.test.no.miss$mat.BW02 / genes.first.filter.BW.test.no.miss$total.BW02
genes.first.filter.BW.test.no.miss$bias.BW03 = genes.first.filter.BW.test.no.miss$mat.BW03 / genes.first.filter.BW.test.no.miss$total.BW03
genes.first.filter.BW.test.no.miss$bias.BW05 = genes.first.filter.BW.test.no.miss$mat.BW05 / genes.first.filter.BW.test.no.miss$total.BW05
genes.first.filter.BW.test.no.miss$bias.BW12 = genes.first.filter.BW.test.no.miss$mat.BW12 / genes.first.filter.BW.test.no.miss$total.BW12

genes.first.filter.BW.test.no.miss$bias.cat.BW02 = "B"
genes.first.filter.BW.test.no.miss$bias.cat.BW03 = "B"
genes.first.filter.BW.test.no.miss$bias.cat.BW05 = "B"
genes.first.filter.BW.test.no.miss$bias.cat.BW12 = "B"

genes.first.filter.BW.test.no.miss$bias.cat.BW02[genes.first.filter.BW.test.no.miss$bias.BW02<=0.050] <-"P"
genes.first.filter.BW.test.no.miss$bias.cat.BW02[genes.first.filter.BW.test.no.miss$bias.BW02>0.050 & genes.first.filter.BW.test.no.miss$bias.BW02<0.350]<-"PB"
genes.first.filter.BW.test.no.miss$bias.cat.BW02[genes.first.filter.BW.test.no.miss$bias.BW02>=0.950]<-"M"
genes.first.filter.BW.test.no.miss$bias.cat.BW02[genes.first.filter.BW.test.no.miss$bias.BW02>0.650 & genes.first.filter.BW.test.no.miss$bias.BW02<0.950]<-"MB"

genes.first.filter.BW.test.no.miss$bias.cat.BW03[genes.first.filter.BW.test.no.miss$bias.BW03<=0.050] <-"P"
genes.first.filter.BW.test.no.miss$bias.cat.BW03[genes.first.filter.BW.test.no.miss$bias.BW03>0.050 & genes.first.filter.BW.test.no.miss$bias.BW03<0.350]<-"PB"
genes.first.filter.BW.test.no.miss$bias.cat.BW03[genes.first.filter.BW.test.no.miss$bias.BW03>=0.950]<-"M"
genes.first.filter.BW.test.no.miss$bias.cat.BW03[genes.first.filter.BW.test.no.miss$bias.BW03>0.650 & genes.first.filter.BW.test.no.miss$bias.BW03<0.950]<-"MB"

genes.first.filter.BW.test.no.miss$bias.cat.BW05[genes.first.filter.BW.test.no.miss$bias.BW05<=0.050] <-"P"
genes.first.filter.BW.test.no.miss$bias.cat.BW05[genes.first.filter.BW.test.no.miss$bias.BW05>0.050 & genes.first.filter.BW.test.no.miss$bias.BW05<0.350]<-"PB"
genes.first.filter.BW.test.no.miss$bias.cat.BW05[genes.first.filter.BW.test.no.miss$bias.BW05>=0.950]<-"M"
genes.first.filter.BW.test.no.miss$bias.cat.BW05[genes.first.filter.BW.test.no.miss$bias.BW05>0.650 & genes.first.filter.BW.test.no.miss$bias.BW05<0.950]<-"MB"

genes.first.filter.BW.test.no.miss$bias.cat.BW12[genes.first.filter.BW.test.no.miss$bias.BW12<=0.050] <-"P"
genes.first.filter.BW.test.no.miss$bias.cat.BW12[genes.first.filter.BW.test.no.miss$bias.BW12>0.050 & genes.first.filter.BW.test.no.miss$bias.BW12<0.350]<-"PB"
genes.first.filter.BW.test.no.miss$bias.cat.BW12[genes.first.filter.BW.test.no.miss$bias.BW12>=0.950]<-"M"
genes.first.filter.BW.test.no.miss$bias.cat.BW12[genes.first.filter.BW.test.no.miss$bias.BW12>0.650 & genes.first.filter.BW.test.no.miss$bias.BW12<0.950]<-"MB"

# decide

genes.first.filter.BW.test.no.miss$pass.G = ifelse(genes.first.filter.BW.test.no.miss$p.value.BW < 0.05/nrow(genes.first.filter.BW.test),"FAIL","PASS")
genes.first.filter.BW.test.no.miss$all.p.significant = ifelse((genes.first.filter.BW.test.no.miss$BW02.p.value < 0.05/nrow(genes.first.filter.BW.test) &
                                                                 genes.first.filter.BW.test.no.miss$BW03.p.value < 0.05/nrow(genes.first.filter.BW.test) &
                                                                 genes.first.filter.BW.test.no.miss$BW05.p.value < 0.05/nrow(genes.first.filter.BW.test) &
                                                                 genes.first.filter.BW.test.no.miss$BW12.p.value < 0.05/nrow(genes.first.filter.BW.test)) |
                                                                (genes.first.filter.BW.test.no.miss$BW02.p.value >= 0.05/nrow(genes.first.filter.BW.test) &
                                                                   genes.first.filter.BW.test.no.miss$BW03.p.value >= 0.05/nrow(genes.first.filter.BW.test) &
                                                                   genes.first.filter.BW.test.no.miss$BW05.p.value >= 0.05/nrow(genes.first.filter.BW.test) &
                                                                   genes.first.filter.BW.test.no.miss$BW12.p.value >= 0.05/nrow(genes.first.filter.BW.test)),
                                                              "AGREE","DISAGREE")

genes.first.filter.BW.test.no.miss$all.bias.cat = ifelse((genes.first.filter.BW.test.no.miss$bias.cat.BW02 == "P" & genes.first.filter.BW.test.no.miss$bias.cat.BW03 == "P" & genes.first.filter.BW.test.no.miss$bias.cat.BW05 == "P" & genes.first.filter.BW.test.no.miss$bias.cat.BW12 == "P") |
                                                           (genes.first.filter.BW.test.no.miss$bias.cat.BW02 == "PB" & genes.first.filter.BW.test.no.miss$bias.cat.BW03 == "PB" & genes.first.filter.BW.test.no.miss$bias.cat.BW05 == "PB" & genes.first.filter.BW.test.no.miss$bias.cat.BW12 == "PB") |
                                                           (genes.first.filter.BW.test.no.miss$bias.cat.BW02 == "B" & genes.first.filter.BW.test.no.miss$bias.cat.BW03 == "B" & genes.first.filter.BW.test.no.miss$bias.cat.BW05 == "B" & genes.first.filter.BW.test.no.miss$bias.cat.BW12 == "B") |
                                                           (genes.first.filter.BW.test.no.miss$bias.cat.BW02 == "MB" & genes.first.filter.BW.test.no.miss$bias.cat.BW03 == "MB" & genes.first.filter.BW.test.no.miss$bias.cat.BW05 == "MB" & genes.first.filter.BW.test.no.miss$bias.cat.BW12 == "MB") |
                                                           (genes.first.filter.BW.test.no.miss$bias.cat.BW02 == "M" & genes.first.filter.BW.test.no.miss$bias.cat.BW03 == "M" & genes.first.filter.BW.test.no.miss$bias.cat.BW05 == "M" & genes.first.filter.BW.test.no.miss$bias.cat.BW12 == "M"),"AGREE","DISAGREE")

genes.first.filter.BW.test.no.miss.fail.filter = genes.first.filter.BW.test.no.miss[genes.first.filter.BW.test.no.miss$pass.G == "FAIL" & (genes.first.filter.BW.test.no.miss$all.p.significant == "DISAGREE" | genes.first.filter.BW.test.no.miss$all.bias.cat == "DISAGREE"),]
genes.first.filter.BW.test.no.miss.pass.filter = anti_join(genes.first.filter.BW.test.no.miss,genes.first.filter.BW.test.no.miss.fail.filter,by="gene")

# no BW02

for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW02)){
  stat = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.BW02[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.BW02[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.BW02[each_gene,]$mat.BW12),
                              c(genes.first.filter.BW.test.no.BW02[each_gene,]$pat.BW03,genes.first.filter.BW.test.no.BW02[each_gene,]$pat.BW05,genes.first.filter.BW.test.no.BW02[each_gene,]$pat.BW12))))$statistic
  genes.first.filter.BW.test.no.BW02$G.BW[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW02)){
  p = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.BW02[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.BW02[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.BW02[each_gene,]$mat.BW12),
                           c(genes.first.filter.BW.test.no.BW02[each_gene,]$pat.BW03,genes.first.filter.BW.test.no.BW02[each_gene,]$pat.BW05,genes.first.filter.BW.test.no.BW02[each_gene,]$pat.BW12))))$p.value
  genes.first.filter.BW.test.no.BW02$p.value.BW[each_gene] = p
}

genes.first.filter.BW.test.no.BW02$pass.G = ifelse(genes.first.filter.BW.test.no.BW02$p.value.BW < 0.05/nrow(genes.first.filter.BW.test),"FAIL","PASS")

#genes.first.filter.BW.test.no.BW02$BW02.p.value = NA
#genes.first.filter.BW.test.no.BW02$BW03.p.value = 0
#genes.first.filter.BW.test.no.BW02$BW05.p.value = 0
#genes.first.filter.BW.test.no.BW02$BW12.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW02)){
#  b = binom.test(genes.first.filter.BW.test.no.BW02[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.BW02[each_gene,]$total.BW03,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW02$BW03.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW02)){
#  c = binom.test(genes.first.filter.BW.test.no.BW02[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.BW02[each_gene,]$total.BW05,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW02$BW05.p.value[each_gene] = c
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW02)){
#  d = binom.test(genes.first.filter.BW.test.no.BW02[each_gene,]$mat.BW12,genes.first.filter.BW.test.no.BW02[each_gene,]$total.BW12,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW02$BW12.p.value[each_gene] = d
#}

# no BW03

for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW03)){
  stat = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.BW03[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.BW03[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.BW03[each_gene,]$mat.BW12),
                              c(genes.first.filter.BW.test.no.BW03[each_gene,]$pat.BW02,genes.first.filter.BW.test.no.BW03[each_gene,]$pat.BW05,genes.first.filter.BW.test.no.BW03[each_gene,]$pat.BW12))))$statistic
  genes.first.filter.BW.test.no.BW03$G.BW[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW03)){
  p = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.BW03[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.BW03[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.BW03[each_gene,]$mat.BW12),
                           c(genes.first.filter.BW.test.no.BW03[each_gene,]$pat.BW02,genes.first.filter.BW.test.no.BW03[each_gene,]$pat.BW05,genes.first.filter.BW.test.no.BW03[each_gene,]$pat.BW12))))$p.value
  genes.first.filter.BW.test.no.BW03$p.value.BW[each_gene] = p
}

genes.first.filter.BW.test.no.BW03$pass.G = ifelse(genes.first.filter.BW.test.no.BW03$p.value.BW < 0.05/nrow(genes.first.filter.BW.test),"FAIL","PASS")

#genes.first.filter.BW.test.no.BW03$BW02.p.value = 0
#genes.first.filter.BW.test.no.BW03$BW03.p.value = NA
#genes.first.filter.BW.test.no.BW03$BW05.p.value = 0
#genes.first.filter.BW.test.no.BW03$BW12.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW03)){
#  a = binom.test(genes.first.filter.BW.test.no.BW03[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.BW03[each_gene,]$total.BW02,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW03$BW02.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW03)){
#  c = binom.test(genes.first.filter.BW.test.no.BW03[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.BW03[each_gene,]$total.BW05,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW03$BW05.p.value[each_gene] = c
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW03)){
#  d = binom.test(genes.first.filter.BW.test.no.BW03[each_gene,]$mat.BW12,genes.first.filter.BW.test.no.BW03[each_gene,]$total.BW12,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW03$BW12.p.value[each_gene] = d
#}

# no BW05

#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW05)){
#  stat = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.BW05[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.BW05[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.BW05[each_gene,]$mat.BW12),
#                              c(genes.first.filter.BW.test.no.BW05[each_gene,]$pat.BW02,genes.first.filter.BW.test.no.BW05[each_gene,]$pat.BW03,genes.first.filter.BW.test.no.BW05[each_gene,]$pat.BW12))))$statistic
#  genes.first.filter.BW.test.no.BW05$G.BW[each_gene] = stat
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW05)){
#  p = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.BW05[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.BW05[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.BW05[each_gene,]$mat.BW12),
#                           c(genes.first.filter.BW.test.no.BW05[each_gene,]$pat.BW02,genes.first.filter.BW.test.no.BW05[each_gene,]$pat.BW03,genes.first.filter.BW.test.no.BW05[each_gene,]$pat.BW12))))$p.value
#  genes.first.filter.BW.test.no.BW05$p.value.BW[each_gene] = p
#}
#
#genes.first.filter.BW.test.no.BW05$pass.G = ifelse(genes.first.filter.BW.test.no.BW05$p.value.BW < 0.05/nrow(genes.first.filter.BW.test),"FAIL","PASS")
#
#genes.first.filter.BW.test.no.BW05$BW02.p.value = 0
#genes.first.filter.BW.test.no.BW05$BW03.p.value = 0
#genes.first.filter.BW.test.no.BW05$BW05.p.value = NA
#genes.first.filter.BW.test.no.BW05$BW12.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW05)){
#  a = binom.test(genes.first.filter.BW.test.no.BW05[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.BW05[each_gene,]$total.BW02,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW05$BW02.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW05)){
#  b = binom.test(genes.first.filter.BW.test.no.BW05[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.BW05[each_gene,]$total.BW03,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW05$BW03.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW05)){
#  d = binom.test(genes.first.filter.BW.test.no.BW05[each_gene,]$mat.BW12,genes.first.filter.BW.test.no.BW05[each_gene,]$total.BW12,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW05$BW12.p.value[each_gene] = d
#}

# no BW12

#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW12)){
#  stat = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.BW12[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.BW12[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.BW12[each_gene,]$mat.BW05),
#                              c(genes.first.filter.BW.test.no.BW12[each_gene,]$pat.BW02,genes.first.filter.BW.test.no.BW12[each_gene,]$pat.BW03,genes.first.filter.BW.test.no.BW12[each_gene,]$pat.BW05))))$statistic
#  genes.first.filter.BW.test.no.BW12$G.BW[each_gene] = stat
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW12)){
#  p = GTest(as.table(rbind(c(genes.first.filter.BW.test.no.BW12[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.BW12[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.BW12[each_gene,]$mat.BW05),
#                           c(genes.first.filter.BW.test.no.BW12[each_gene,]$pat.BW02,genes.first.filter.BW.test.no.BW12[each_gene,]$pat.BW03,genes.first.filter.BW.test.no.BW12[each_gene,]$pat.BW05))))$p.value
#  genes.first.filter.BW.test.no.BW12$p.value.BW[each_gene] = p
#}
#
#genes.first.filter.BW.test.no.BW12$pass.G = ifelse(genes.first.filter.BW.test.no.BW12$p.value.BW < 0.05/nrow(genes.first.filter.BW.test),"FAIL","PASS")

#genes.first.filter.BW.test.no.BW12$BW02.p.value = 0
#genes.first.filter.BW.test.no.BW12$BW03.p.value = 0
#genes.first.filter.BW.test.no.BW12$BW05.p.value = 0
#genes.first.filter.BW.test.no.BW12$BW12.p.value = NA
#
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW12)){
#  a = binom.test(genes.first.filter.BW.test.no.BW12[each_gene,]$mat.BW02,genes.first.filter.BW.test.no.BW12[each_gene,]$total.BW02,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW12$BW02.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW12)){
#  b = binom.test(genes.first.filter.BW.test.no.BW12[each_gene,]$mat.BW03,genes.first.filter.BW.test.no.BW12[each_gene,]$total.BW03,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW12$BW03.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.BW.test.no.BW12)){
#  c = binom.test(genes.first.filter.BW.test.no.BW12[each_gene,]$mat.BW05,genes.first.filter.BW.test.no.BW12[each_gene,]$total.BW05,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BW.test.no.BW12$BW05.p.value[each_gene] = c
#}

nrow(genes.first.filter.WB.test.no.miss.fail.filter)
nrow(genes.first.filter.BW.test.no.miss.fail.filter)
genes.second.filter.WB.fail = rbind(genes.first.filter.WB.test.no.miss.fail.filter[c("gene")],genes.first.filter.BW.test.no.miss.fail.filter[c("gene")])
length(unique(genes.second.filter.WB.fail$gene))
nrow(genes.first.filter.WB) - length(unique(genes.second.filter.WB.fail$gene))
genes.second.filter.WB = anti_join(genes.first.filter.WB,genes.second.filter.WB.fail,by="gene")

### CP x BGOX

## CB

# split by number of replicates

genes.first.filter.CB.test.no.miss <- genes.first.filter.CB.test[genes.first.filter.CB.test$num.rep.CB == 3, ] # no need to split, only three replicates

# G test

for(each_gene in 1:nrow(genes.first.filter.CB.test.no.miss)){
  stat = GTest(as.table(rbind(c(genes.first.filter.CB.test.no.miss[each_gene,]$mat.CB01,genes.first.filter.CB.test.no.miss[each_gene,]$mat.CB12,genes.first.filter.CB.test.no.miss[each_gene,]$mat.CB20),
                              c(genes.first.filter.CB.test.no.miss[each_gene,]$pat.CB01,genes.first.filter.CB.test.no.miss[each_gene,]$pat.CB12,genes.first.filter.CB.test.no.miss[each_gene,]$pat.CB20))))$statistic
  genes.first.filter.CB.test.no.miss$G.CB[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.CB.test.no.miss)){
  p = GTest(as.table(rbind(c(genes.first.filter.CB.test.no.miss[each_gene,]$mat.CB01,genes.first.filter.CB.test.no.miss[each_gene,]$mat.CB12,genes.first.filter.CB.test.no.miss[each_gene,]$mat.CB20),
                           c(genes.first.filter.CB.test.no.miss[each_gene,]$pat.CB01,genes.first.filter.CB.test.no.miss[each_gene,]$pat.CB12,genes.first.filter.CB.test.no.miss[each_gene,]$pat.CB20))))$p.value
  genes.first.filter.CB.test.no.miss$p.value.CB[each_gene] = p
}


# individual exact binomial test

genes.first.filter.CB.test.no.miss$CB01.p.value = 0
genes.first.filter.CB.test.no.miss$CB12.p.value = 0
genes.first.filter.CB.test.no.miss$CB20.p.value = 0

for(each_gene in 1:nrow(genes.first.filter.CB.test.no.miss)){
  a = binom.test(genes.first.filter.CB.test.no.miss[each_gene,]$mat.CB01,genes.first.filter.CB.test.no.miss[each_gene,]$total.CB01,0.5,alternative="two.sided")$p.value
  genes.first.filter.CB.test.no.miss$CB01.p.value[each_gene] = a
}
for(each_gene in 1:nrow(genes.first.filter.CB.test.no.miss)){
  b = binom.test(genes.first.filter.CB.test.no.miss[each_gene,]$mat.CB12,genes.first.filter.CB.test.no.miss[each_gene,]$total.CB12,0.5,alternative="two.sided")$p.value
  genes.first.filter.CB.test.no.miss$CB12.p.value[each_gene] = b
}
for(each_gene in 1:nrow(genes.first.filter.CB.test.no.miss)){
  c = binom.test(genes.first.filter.CB.test.no.miss[each_gene,]$mat.CB20,genes.first.filter.CB.test.no.miss[each_gene,]$total.CB20,0.5,alternative="two.sided")$p.value
  genes.first.filter.CB.test.no.miss$CB20.p.value[each_gene] = c
}

# individual category of bias

genes.first.filter.CB.test.no.miss$bias.CB01 = genes.first.filter.CB.test.no.miss$mat.CB01 / genes.first.filter.CB.test.no.miss$total.CB01
genes.first.filter.CB.test.no.miss$bias.CB12 = genes.first.filter.CB.test.no.miss$mat.CB12 / genes.first.filter.CB.test.no.miss$total.CB12
genes.first.filter.CB.test.no.miss$bias.CB20 = genes.first.filter.CB.test.no.miss$mat.CB20 / genes.first.filter.CB.test.no.miss$total.CB20

genes.first.filter.CB.test.no.miss$bias.cat.CB01 = "B"
genes.first.filter.CB.test.no.miss$bias.cat.CB12 = "B"
genes.first.filter.CB.test.no.miss$bias.cat.CB20 = "B"

genes.first.filter.CB.test.no.miss$bias.cat.CB01[genes.first.filter.CB.test.no.miss$bias.CB01<=0.050] <-"P"
genes.first.filter.CB.test.no.miss$bias.cat.CB01[genes.first.filter.CB.test.no.miss$bias.CB01>0.050 & genes.first.filter.CB.test.no.miss$bias.CB01<0.350]<-"PB"
genes.first.filter.CB.test.no.miss$bias.cat.CB01[genes.first.filter.CB.test.no.miss$bias.CB01>=0.950]<-"M"
genes.first.filter.CB.test.no.miss$bias.cat.CB01[genes.first.filter.CB.test.no.miss$bias.CB01>0.650 & genes.first.filter.CB.test.no.miss$bias.CB01<0.950]<-"MB"

genes.first.filter.CB.test.no.miss$bias.cat.CB12[genes.first.filter.CB.test.no.miss$bias.CB12<=0.050] <-"P"
genes.first.filter.CB.test.no.miss$bias.cat.CB12[genes.first.filter.CB.test.no.miss$bias.CB12>0.050 & genes.first.filter.CB.test.no.miss$bias.CB12<0.350]<-"PB"
genes.first.filter.CB.test.no.miss$bias.cat.CB12[genes.first.filter.CB.test.no.miss$bias.CB12>=0.950]<-"M"
genes.first.filter.CB.test.no.miss$bias.cat.CB12[genes.first.filter.CB.test.no.miss$bias.CB12>0.650 & genes.first.filter.CB.test.no.miss$bias.CB12<0.950]<-"MB"

genes.first.filter.CB.test.no.miss$bias.cat.CB20[genes.first.filter.CB.test.no.miss$bias.CB20<=0.050] <-"P"
genes.first.filter.CB.test.no.miss$bias.cat.CB20[genes.first.filter.CB.test.no.miss$bias.CB20>0.050 & genes.first.filter.CB.test.no.miss$bias.CB20<0.350]<-"PB"
genes.first.filter.CB.test.no.miss$bias.cat.CB20[genes.first.filter.CB.test.no.miss$bias.CB20>=0.950]<-"M"
genes.first.filter.CB.test.no.miss$bias.cat.CB20[genes.first.filter.CB.test.no.miss$bias.CB20>0.650 & genes.first.filter.CB.test.no.miss$bias.CB20<0.950]<-"MB"

# decide

genes.first.filter.CB.test.no.miss$pass.G = ifelse(genes.first.filter.CB.test.no.miss$p.value.CB < 0.05/nrow(genes.first.filter.CB.test),"FAIL","PASS")
genes.first.filter.CB.test.no.miss$all.p.significant = ifelse((genes.first.filter.CB.test.no.miss$CB01.p.value < 0.05/nrow(genes.first.filter.CB.test) &
                                                                 genes.first.filter.CB.test.no.miss$CB12.p.value < 0.05/nrow(genes.first.filter.CB.test) &
                                                                 genes.first.filter.CB.test.no.miss$CB20.p.value < 0.05/nrow(genes.first.filter.CB.test)) |
                                                                (genes.first.filter.CB.test.no.miss$CB01.p.value >= 0.05/nrow(genes.first.filter.CB.test) &
                                                                   genes.first.filter.CB.test.no.miss$CB12.p.value >= 0.05/nrow(genes.first.filter.CB.test) &
                                                                   genes.first.filter.CB.test.no.miss$CB20.p.value >= 0.05/nrow(genes.first.filter.CB.test)),"AGREE","DISAGREE")

genes.first.filter.CB.test.no.miss$all.bias.cat = ifelse((genes.first.filter.CB.test.no.miss$bias.cat.CB01 == "P" & genes.first.filter.CB.test.no.miss$bias.cat.CB12 == "P" & genes.first.filter.CB.test.no.miss$bias.cat.CB20 == "P") |
                                                           (genes.first.filter.CB.test.no.miss$bias.cat.CB01 == "PB" & genes.first.filter.CB.test.no.miss$bias.cat.CB12 == "PB" & genes.first.filter.CB.test.no.miss$bias.cat.CB20 == "PB") |
                                                           (genes.first.filter.CB.test.no.miss$bias.cat.CB01 == "B" & genes.first.filter.CB.test.no.miss$bias.cat.CB12 == "B" & genes.first.filter.CB.test.no.miss$bias.cat.CB20 == "B") |
                                                           (genes.first.filter.CB.test.no.miss$bias.cat.CB01 == "MB" & genes.first.filter.CB.test.no.miss$bias.cat.CB12 == "MB" & genes.first.filter.CB.test.no.miss$bias.cat.CB20 == "MB") |
                                                           (genes.first.filter.CB.test.no.miss$bias.cat.CB01 == "M" & genes.first.filter.CB.test.no.miss$bias.cat.CB12 == "M" & genes.first.filter.CB.test.no.miss$bias.cat.CB20 == "M"),"AGREE","DISAGREE")

genes.first.filter.CB.test.no.miss.fail.filter = genes.first.filter.CB.test.no.miss[genes.first.filter.CB.test.no.miss$pass.G == "FAIL" & (genes.first.filter.CB.test.no.miss$all.p.significant == "DISAGREE" | genes.first.filter.CB.test.no.miss$all.bias.cat == "DISAGREE"),]

nrow(genes.first.filter.CB.test.no.miss.fail.filter)
nrow(genes.first.filter.CB.test.no.miss.pass.filter)
genes.first.filter.CB.test.no.miss.pass.filter = anti_join(genes.first.filter.CB.test.no.miss,genes.first.filter.CB.test.no.miss.fail.filter,by="gene")

### BC

# BC

# split by number of replicates

genes.first.filter.BC.test.no.miss <- genes.first.filter.BC.test[genes.first.filter.BC.test$num.rep.BC == 4, ]
genes.first.filter.BC.test.no.BC03 <- genes.first.filter.BC.test[genes.first.filter.BC.test$total.BC03 == 0, ]
genes.first.filter.BC.test.no.BC04 <- genes.first.filter.BC.test[genes.first.filter.BC.test$total.BC04 == 0, ]
genes.first.filter.BC.test.no.BC13 <- genes.first.filter.BC.test[genes.first.filter.BC.test$total.BC13 == 0, ]
genes.first.filter.BC.test.no.BC17 <- genes.first.filter.BC.test[genes.first.filter.BC.test$total.BC17 == 0, ]

nrow(genes.first.filter.BC.test.no.miss)
nrow(genes.first.filter.BC.test.no.BC03)
nrow(genes.first.filter.BC.test.no.BC04)
nrow(genes.first.filter.BC.test.no.BC13)
nrow(genes.first.filter.BC.test.no.BC17)

# no missing replicates

# G test

for(each_gene in 1:nrow(genes.first.filter.BC.test.no.miss)){
  stat = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC17),
                              c(genes.first.filter.BC.test.no.miss[each_gene,]$pat.BC03,genes.first.filter.BC.test.no.miss[each_gene,]$pat.BC04,genes.first.filter.BC.test.no.miss[each_gene,]$pat.BC13,genes.first.filter.BC.test.no.miss[each_gene,]$pat.BC17))))$statistic
  genes.first.filter.BC.test.no.miss$G.BC[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.BC.test.no.miss)){
  p = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC17),
                           c(genes.first.filter.BC.test.no.miss[each_gene,]$pat.BC03,genes.first.filter.BC.test.no.miss[each_gene,]$pat.BC04,genes.first.filter.BC.test.no.miss[each_gene,]$pat.BC13,genes.first.filter.BC.test.no.miss[each_gene,]$pat.BC17))))$p.value
  genes.first.filter.BC.test.no.miss$p.value.BC[each_gene] = p
}

# individual exact binomial test

genes.first.filter.BC.test.no.miss$BC03.p.value = 0
genes.first.filter.BC.test.no.miss$BC04.p.value = 0
genes.first.filter.BC.test.no.miss$BC13.p.value = 0
genes.first.filter.BC.test.no.miss$BC17.p.value = 0

for(each_gene in 1:nrow(genes.first.filter.BC.test.no.miss)){
  a = binom.test(genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.miss[each_gene,]$total.BC03,0.5,alternative="two.sided")$p.value
  genes.first.filter.BC.test.no.miss$BC03.p.value[each_gene] = a
}
for(each_gene in 1:nrow(genes.first.filter.BC.test.no.miss)){
  b = binom.test(genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.miss[each_gene,]$total.BC04,0.5,alternative="two.sided")$p.value
  genes.first.filter.BC.test.no.miss$BC04.p.value[each_gene] = b
}
for(each_gene in 1:nrow(genes.first.filter.BC.test.no.miss)){
  c = binom.test(genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.miss[each_gene,]$total.BC13,0.5,alternative="two.sided")$p.value
  genes.first.filter.BC.test.no.miss$BC13.p.value[each_gene] = c
}
for(each_gene in 1:nrow(genes.first.filter.BC.test.no.miss)){
  d = binom.test(genes.first.filter.BC.test.no.miss[each_gene,]$mat.BC17,genes.first.filter.BC.test.no.miss[each_gene,]$total.BC17,0.5,alternative="two.sided")$p.value
  genes.first.filter.BC.test.no.miss$BC17.p.value[each_gene] = d
}

# individual category of bias

genes.first.filter.BC.test.no.miss$bias.BC03 = genes.first.filter.BC.test.no.miss$mat.BC03 / genes.first.filter.BC.test.no.miss$total.BC03
genes.first.filter.BC.test.no.miss$bias.BC04 = genes.first.filter.BC.test.no.miss$mat.BC04 / genes.first.filter.BC.test.no.miss$total.BC04
genes.first.filter.BC.test.no.miss$bias.BC13 = genes.first.filter.BC.test.no.miss$mat.BC13 / genes.first.filter.BC.test.no.miss$total.BC13
genes.first.filter.BC.test.no.miss$bias.BC17 = genes.first.filter.BC.test.no.miss$mat.BC17 / genes.first.filter.BC.test.no.miss$total.BC17

genes.first.filter.BC.test.no.miss$bias.cat.BC03 = "B"
genes.first.filter.BC.test.no.miss$bias.cat.BC04 = "B"
genes.first.filter.BC.test.no.miss$bias.cat.BC13 = "B"
genes.first.filter.BC.test.no.miss$bias.cat.BC17 = "B"

genes.first.filter.BC.test.no.miss$bias.cat.BC03[genes.first.filter.BC.test.no.miss$bias.BC03<=0.050] <-"P"
genes.first.filter.BC.test.no.miss$bias.cat.BC03[genes.first.filter.BC.test.no.miss$bias.BC03>0.050 & genes.first.filter.BC.test.no.miss$bias.BC03<0.350]<-"PB"
genes.first.filter.BC.test.no.miss$bias.cat.BC03[genes.first.filter.BC.test.no.miss$bias.BC03>=0.950]<-"M"
genes.first.filter.BC.test.no.miss$bias.cat.BC03[genes.first.filter.BC.test.no.miss$bias.BC03>0.650 & genes.first.filter.BC.test.no.miss$bias.BC03<0.950]<-"MB"

genes.first.filter.BC.test.no.miss$bias.cat.BC04[genes.first.filter.BC.test.no.miss$bias.BC04<=0.050] <-"P"
genes.first.filter.BC.test.no.miss$bias.cat.BC04[genes.first.filter.BC.test.no.miss$bias.BC04>0.050 & genes.first.filter.BC.test.no.miss$bias.BC04<0.350]<-"PB"
genes.first.filter.BC.test.no.miss$bias.cat.BC04[genes.first.filter.BC.test.no.miss$bias.BC04>=0.950]<-"M"
genes.first.filter.BC.test.no.miss$bias.cat.BC04[genes.first.filter.BC.test.no.miss$bias.BC04>0.650 & genes.first.filter.BC.test.no.miss$bias.BC04<0.950]<-"MB"

genes.first.filter.BC.test.no.miss$bias.cat.BC13[genes.first.filter.BC.test.no.miss$bias.BC13<=0.050] <-"P"
genes.first.filter.BC.test.no.miss$bias.cat.BC13[genes.first.filter.BC.test.no.miss$bias.BC13>0.050 & genes.first.filter.BC.test.no.miss$bias.BC13<0.350]<-"PB"
genes.first.filter.BC.test.no.miss$bias.cat.BC13[genes.first.filter.BC.test.no.miss$bias.BC13>=0.950]<-"M"
genes.first.filter.BC.test.no.miss$bias.cat.BC13[genes.first.filter.BC.test.no.miss$bias.BC13>0.650 & genes.first.filter.BC.test.no.miss$bias.BC13<0.950]<-"MB"

genes.first.filter.BC.test.no.miss$bias.cat.BC17[genes.first.filter.BC.test.no.miss$bias.BC17<=0.050] <-"P"
genes.first.filter.BC.test.no.miss$bias.cat.BC17[genes.first.filter.BC.test.no.miss$bias.BC17>0.050 & genes.first.filter.BC.test.no.miss$bias.BC17<0.350]<-"PB"
genes.first.filter.BC.test.no.miss$bias.cat.BC17[genes.first.filter.BC.test.no.miss$bias.BC17>=0.950]<-"M"
genes.first.filter.BC.test.no.miss$bias.cat.BC17[genes.first.filter.BC.test.no.miss$bias.BC17>0.650 & genes.first.filter.BC.test.no.miss$bias.BC17<0.950]<-"MB"

# decide

genes.first.filter.BC.test.no.miss$pass.G = ifelse(genes.first.filter.BC.test.no.miss$p.value.BC < 0.05/nrow(genes.first.filter.BC.test),"FAIL","PASS")
genes.first.filter.BC.test.no.miss$all.p.significant = ifelse((genes.first.filter.BC.test.no.miss$BC03.p.value < 0.05/nrow(genes.first.filter.BC.test) &
                                                                 genes.first.filter.BC.test.no.miss$BC04.p.value < 0.05/nrow(genes.first.filter.BC.test) &
                                                                 genes.first.filter.BC.test.no.miss$BC13.p.value < 0.05/nrow(genes.first.filter.BC.test) &
                                                                 genes.first.filter.BC.test.no.miss$BC17.p.value < 0.05/nrow(genes.first.filter.BC.test)) |
                                                                (genes.first.filter.BC.test.no.miss$BC03.p.value >= 0.05/nrow(genes.first.filter.BC.test) &
                                                                   genes.first.filter.BC.test.no.miss$BC04.p.value >= 0.05/nrow(genes.first.filter.BC.test) &
                                                                   genes.first.filter.BC.test.no.miss$BC13.p.value >= 0.05/nrow(genes.first.filter.BC.test) &
                                                                   genes.first.filter.BC.test.no.miss$BC17.p.value >= 0.05/nrow(genes.first.filter.BC.test)),
                                                              "AGREE","DISAGREE")

genes.first.filter.BC.test.no.miss$all.bias.cat = ifelse((genes.first.filter.BC.test.no.miss$bias.cat.BC03 == "P" & genes.first.filter.BC.test.no.miss$bias.cat.BC04 == "P" & genes.first.filter.BC.test.no.miss$bias.cat.BC13 == "P" & genes.first.filter.BC.test.no.miss$bias.cat.BC17 == "P") |
                                                           (genes.first.filter.BC.test.no.miss$bias.cat.BC03 == "PB" & genes.first.filter.BC.test.no.miss$bias.cat.BC04 == "PB" & genes.first.filter.BC.test.no.miss$bias.cat.BC13 == "PB" & genes.first.filter.BC.test.no.miss$bias.cat.BC17 == "PB") |
                                                           (genes.first.filter.BC.test.no.miss$bias.cat.BC03 == "B" & genes.first.filter.BC.test.no.miss$bias.cat.BC04 == "B" & genes.first.filter.BC.test.no.miss$bias.cat.BC13 == "B" & genes.first.filter.BC.test.no.miss$bias.cat.BC17 == "B") |
                                                           (genes.first.filter.BC.test.no.miss$bias.cat.BC03 == "MB" & genes.first.filter.BC.test.no.miss$bias.cat.BC04 == "MB" & genes.first.filter.BC.test.no.miss$bias.cat.BC13 == "MB" & genes.first.filter.BC.test.no.miss$bias.cat.BC17 == "MB") |
                                                           (genes.first.filter.BC.test.no.miss$bias.cat.BC03 == "M" & genes.first.filter.BC.test.no.miss$bias.cat.BC04 == "M" & genes.first.filter.BC.test.no.miss$bias.cat.BC13 == "M" & genes.first.filter.BC.test.no.miss$bias.cat.BC17 == "M"),"AGREE","DISAGREE")

genes.first.filter.BC.test.no.miss.fail.filter = genes.first.filter.BC.test.no.miss[genes.first.filter.BC.test.no.miss$pass.G == "FAIL" & (genes.first.filter.BC.test.no.miss$all.p.significant == "DISAGREE" | genes.first.filter.BC.test.no.miss$all.bias.cat == "DISAGREE"),]
genes.first.filter.BC.test.no.miss.pass.filter = anti_join(genes.first.filter.BC.test.no.miss,genes.first.filter.BC.test.no.miss.fail.filter,by="gene")

# no BC03

for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC03)){
  stat = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.BC03[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.BC03[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.BC03[each_gene,]$mat.BC17),
                              c(genes.first.filter.BC.test.no.BC03[each_gene,]$pat.BC04,genes.first.filter.BC.test.no.BC03[each_gene,]$pat.BC13,genes.first.filter.BC.test.no.BC03[each_gene,]$pat.BC17))))$statistic
  genes.first.filter.BC.test.no.BC03$G.BC[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC03)){
  p = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.BC03[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.BC03[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.BC03[each_gene,]$mat.BC17),
                           c(genes.first.filter.BC.test.no.BC03[each_gene,]$pat.BC04,genes.first.filter.BC.test.no.BC03[each_gene,]$pat.BC13,genes.first.filter.BC.test.no.BC03[each_gene,]$pat.BC17))))$p.value
  genes.first.filter.BC.test.no.BC03$p.value.BC[each_gene] = p
}

genes.first.filter.BC.test.no.BC03$pass.G = ifelse(genes.first.filter.BC.test.no.BC03$p.value.BC < 0.05/nrow(genes.first.filter.BC.test),"FAIL","PASS")

#genes.first.filter.BC.test.no.BC03$BC03.p.value = NA
#genes.first.filter.BC.test.no.BC03$BC04.p.value = 0
#genes.first.filter.BC.test.no.BC03$BC13.p.value = 0
#genes.first.filter.BC.test.no.BC03$BC17.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC03)){
#  b = binom.test(genes.first.filter.BC.test.no.BC03[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.BC03[each_gene,]$total.BC04,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC03$BC04.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC03)){
#  c = binom.test(genes.first.filter.BC.test.no.BC03[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.BC03[each_gene,]$total.BC13,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC03$BC13.p.value[each_gene] = c
#}
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC03)){
#  d = binom.test(genes.first.filter.BC.test.no.BC03[each_gene,]$mat.BC17,genes.first.filter.BC.test.no.BC03[each_gene,]$total.BC17,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC03$BC17.p.value[each_gene] = d
#}

# no BC04

for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC04)){
  stat = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.BC04[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.BC04[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.BC04[each_gene,]$mat.BC17),
                              c(genes.first.filter.BC.test.no.BC04[each_gene,]$pat.BC03,genes.first.filter.BC.test.no.BC04[each_gene,]$pat.BC13,genes.first.filter.BC.test.no.BC04[each_gene,]$pat.BC17))))$statistic
  genes.first.filter.BC.test.no.BC04$G.BC[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC04)){
  p = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.BC04[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.BC04[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.BC04[each_gene,]$mat.BC17),
                           c(genes.first.filter.BC.test.no.BC04[each_gene,]$pat.BC03,genes.first.filter.BC.test.no.BC04[each_gene,]$pat.BC13,genes.first.filter.BC.test.no.BC04[each_gene,]$pat.BC17))))$p.value
  genes.first.filter.BC.test.no.BC04$p.value.BC[each_gene] = p
}

genes.first.filter.BC.test.no.BC04$pass.G = ifelse(genes.first.filter.BC.test.no.BC04$p.value.BC < 0.05/nrow(genes.first.filter.BC.test),"FAIL","PASS")

#genes.first.filter.BC.test.no.BC04$BC03.p.value = 0
#genes.first.filter.BC.test.no.BC04$BC04.p.value = NA
#genes.first.filter.BC.test.no.BC04$BC13.p.value = 0
#genes.first.filter.BC.test.no.BC04$BC17.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC04)){
#  a = binom.test(genes.first.filter.BC.test.no.BC04[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.BC04[each_gene,]$total.BC03,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC04$BC03.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC04)){
#  c = binom.test(genes.first.filter.BC.test.no.BC04[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.BC04[each_gene,]$total.BC13,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC04$BC13.p.value[each_gene] = c
#}
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC04)){
#  d = binom.test(genes.first.filter.BC.test.no.BC04[each_gene,]$mat.BC17,genes.first.filter.BC.test.no.BC04[each_gene,]$total.BC17,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC04$BC17.p.value[each_gene] = d
#}

# no BC13

#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC13)){
#  stat = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.BC13[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.BC13[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.BC13[each_gene,]$mat.BC17),
#                              c(genes.first.filter.BC.test.no.BC13[each_gene,]$pat.BC03,genes.first.filter.BC.test.no.BC13[each_gene,]$pat.BC04,genes.first.filter.BC.test.no.BC13[each_gene,]$pat.BC17))))$statistic
#  genes.first.filter.BC.test.no.BC13$G.BC[each_gene] = stat
#}
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC13)){
#  p = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.BC13[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.BC13[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.BC13[each_gene,]$mat.BC17),
#                           c(genes.first.filter.BC.test.no.BC13[each_gene,]$pat.BC03,genes.first.filter.BC.test.no.BC13[each_gene,]$pat.BC04,genes.first.filter.BC.test.no.BC13[each_gene,]$pat.BC17))))$p.value
#  genes.first.filter.BC.test.no.BC13$p.value.BC[each_gene] = p
#}
#
#genes.first.filter.BC.test.no.BC13$pass.G = ifelse(genes.first.filter.BC.test.no.BC13$p.value.BC < 0.05/nrow(genes.first.filter.BC.test),"FAIL","PASS")
#
#genes.first.filter.BC.test.no.BC13$BC03.p.value = 0
#genes.first.filter.BC.test.no.BC13$BC04.p.value = 0
#genes.first.filter.BC.test.no.BC13$BC13.p.value = NA
#genes.first.filter.BC.test.no.BC13$BC17.p.value = 0
#
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC13)){
#  a = binom.test(genes.first.filter.BC.test.no.BC13[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.BC13[each_gene,]$total.BC03,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC13$BC03.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC13)){
#  b = binom.test(genes.first.filter.BC.test.no.BC13[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.BC13[each_gene,]$total.BC04,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC13$BC04.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC13)){
#  d = binom.test(genes.first.filter.BC.test.no.BC13[each_gene,]$mat.BC17,genes.first.filter.BC.test.no.BC13[each_gene,]$total.BC17,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC13$BC17.p.value[each_gene] = d
#}

# no BC17

for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC17)){
  stat = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.BC17[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.BC17[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.BC17[each_gene,]$mat.BC13),
                              c(genes.first.filter.BC.test.no.BC17[each_gene,]$pat.BC03,genes.first.filter.BC.test.no.BC17[each_gene,]$pat.BC04,genes.first.filter.BC.test.no.BC17[each_gene,]$pat.BC13))))$statistic
  genes.first.filter.BC.test.no.BC17$G.BC[each_gene] = stat
}
for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC17)){
  p = GTest(as.table(rbind(c(genes.first.filter.BC.test.no.BC17[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.BC17[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.BC17[each_gene,]$mat.BC13),
                           c(genes.first.filter.BC.test.no.BC17[each_gene,]$pat.BC03,genes.first.filter.BC.test.no.BC17[each_gene,]$pat.BC04,genes.first.filter.BC.test.no.BC17[each_gene,]$pat.BC13))))$p.value
  genes.first.filter.BC.test.no.BC17$p.value.BC[each_gene] = p
}

genes.first.filter.BC.test.no.BC17$pass.G = ifelse(genes.first.filter.BC.test.no.BC17$p.value.BC < 0.05/nrow(genes.first.filter.BC.test),"FAIL","PASS")

#genes.first.filter.BC.test.no.BC17$BC03.p.value = 0
#genes.first.filter.BC.test.no.BC17$BC04.p.value = 0
#genes.first.filter.BC.test.no.BC17$BC13.p.value = 0
#genes.first.filter.BC.test.no.BC17$BC17.p.value = NA
#
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC17)){
#  a = binom.test(genes.first.filter.BC.test.no.BC17[each_gene,]$mat.BC03,genes.first.filter.BC.test.no.BC17[each_gene,]$total.BC03,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC17$BC03.p.value[each_gene] = a
#}
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC17)){
#  b = binom.test(genes.first.filter.BC.test.no.BC17[each_gene,]$mat.BC04,genes.first.filter.BC.test.no.BC17[each_gene,]$total.BC04,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC17$BC04.p.value[each_gene] = b
#}
#for(each_gene in 1:nrow(genes.first.filter.BC.test.no.BC17)){
#  c = binom.test(genes.first.filter.BC.test.no.BC17[each_gene,]$mat.BC13,genes.first.filter.BC.test.no.BC17[each_gene,]$total.BC13,0.5,alternative="two.sided")$p.value
#  genes.first.filter.BC.test.no.BC17$BC13.p.value[each_gene] = c
#}

nrow(genes.first.filter.CB.test.no.miss.fail.filter)
nrow(genes.first.filter.BC.test.no.miss.fail.filter)
genes.second.filter.CB.fail = rbind(genes.first.filter.CB.test.no.miss.fail.filter[c("gene")],genes.first.filter.BC.test.no.miss.fail.filter[c("gene")])
length(unique(genes.second.filter.CB.fail$gene))
nrow(genes.second.filter.CB)
genes.second.filter.CB = anti_join(genes.first.filter.CB,genes.second.filter.CB.fail,by="gene")

# Final filter: remove genes with TPM < 1

intra_tpm = WC01_genes_results[c(1,6)]
colnames(intra_tpm)[colnames(intra_tpm)=="TPM"] <- "WC01"
intra_tpm$WC04 = WC04_genes_results$TPM
intra_tpm$WC06 = WC06_genes_results$TPM
intra_tpm$WC19 = WC19_genes_results$TPM
intra_tpm$CW01 = CW01_genes_results$TPM
intra_tpm$CW04 = CW04_genes_results$TPM
intra_tpm$CW10 = CW10_genes_results$TPM
intra_tpm$CW11 = CW11_genes_results$TPM
intra_tpm$WB02 = WB02_genes_results$TPM
intra_tpm$WB03 = WB03_genes_results$TPM
intra_tpm$WB08 = WB08_genes_results$TPM
intra_tpm$BW02 = BW02_genes_results$TPM
intra_tpm$BW03 = BW03_genes_results$TPM
intra_tpm$BW05 = BW05_genes_results$TPM
intra_tpm$BW12 = BW12_genes_results$TPM
intra_tpm$CB01 = CB01_genes_results$TPM
intra_tpm$CB12 = CB12_genes_results$TPM
intra_tpm$CB20 = CB20_genes_results$TPM
intra_tpm$BC03 = BC03_genes_results$TPM
intra_tpm$BC04 = BC04_genes_results$TPM
intra_tpm$BC13 = BC13_genes_results$TPM
intra_tpm$BC17 = BC17_genes_results$TPM

intra_tpm$mean.tpm.WC <- rowMeans(intra_tpm[,2:9])
intra_tpm$mean.tpm.WB <- rowMeans(intra_tpm[,10:16])
intra_tpm$mean.tpm.CB <- rowMeans(intra_tpm[,17:23])

colnames(intra_tpm)[colnames(intra_tpm)=="gene_id"] <- "gene"

WC_tpm_expressed <- intra_tpm[intra_tpm$mean.tpm.WC >= 1,]
WB_tpm_expressed <- intra_tpm[intra_tpm$mean.tpm.WB >= 1,]
CB_tpm_expressed <- intra_tpm[intra_tpm$mean.tpm.CB >= 1,]

WC_tpm_not_expressed <- intra_tpm[intra_tpm$mean.tpm.WC < 1,]
WB_tpm_not_expressed <- intra_tpm[intra_tpm$mean.tpm.WB < 1,]
CB_tpm_not_expressed <- intra_tpm[intra_tpm$mean.tpm.CB < 1,]

genes.final.filter.WC = anti_join(genes.second.filter.WC,WC_tpm_not_expressed,by="gene")
genes.final.filter.WB = anti_join(genes.second.filter.WB,WB_tpm_not_expressed,by="gene")
genes.final.filter.CB = anti_join(genes.second.filter.CB,CB_tpm_not_expressed,by="gene")

nrow(genes.final.filter.WC)
nrow(genes.final.filter.WB)
nrow(genes.final.filter.CB)

########## calculate

# bias

genes.final.filter.WC$mat.WC = genes.final.filter.WC$mat.WC01 + genes.final.filter.WC$mat.WC04 + genes.final.filter.WC$mat.WC06 + genes.final.filter.WC$mat.WC19
genes.final.filter.WC$mat.CW = genes.final.filter.WC$mat.CW01 + genes.final.filter.WC$mat.CW04 + genes.final.filter.WC$mat.CW10 + genes.final.filter.WC$mat.CW11
genes.final.filter.WB$mat.WB = genes.final.filter.WB$mat.WB02 + genes.final.filter.WB$mat.WB03 + genes.final.filter.WB$mat.WB08
genes.final.filter.WB$mat.BW = genes.final.filter.WB$mat.BW02 + genes.final.filter.WB$mat.BW03 + genes.final.filter.WB$mat.BW05 + genes.final.filter.WB$mat.BW12
genes.final.filter.CB$mat.CB = genes.final.filter.CB$mat.CB01 + genes.final.filter.CB$mat.CB12 + genes.final.filter.CB$mat.CB20
genes.final.filter.CB$mat.BC = genes.final.filter.CB$mat.BC03 + genes.final.filter.CB$mat.BC04 + genes.final.filter.CB$mat.BC13 + genes.final.filter.CB$mat.BC17

genes.final.filter.WC$total.WC = genes.final.filter.WC$total.WC01 + genes.final.filter.WC$total.WC04 + genes.final.filter.WC$total.WC06 + genes.final.filter.WC$total.WC19
genes.final.filter.WC$total.CW = genes.final.filter.WC$total.CW01 + genes.final.filter.WC$total.CW04 + genes.final.filter.WC$total.CW10 + genes.final.filter.WC$total.CW11
genes.final.filter.WB$total.WB = genes.final.filter.WB$total.WB02 + genes.final.filter.WB$total.WB03 + genes.final.filter.WB$total.WB08
genes.final.filter.WB$total.BW = genes.final.filter.WB$total.BW02 + genes.final.filter.WB$total.BW03 + genes.final.filter.WB$total.BW05 + genes.final.filter.WB$total.BW12
genes.final.filter.CB$total.CB = genes.final.filter.CB$total.CB01 + genes.final.filter.CB$total.CB12 + genes.final.filter.CB$total.CB20
genes.final.filter.CB$total.BC = genes.final.filter.CB$total.BC03 + genes.final.filter.CB$total.BC04 + genes.final.filter.CB$total.BC13 + genes.final.filter.CB$total.BC17

genes.final.filter.WC$bias.WC = genes.final.filter.WC$mat.WC / genes.final.filter.WC$total.WC
genes.final.filter.WC$bias.CW = genes.final.filter.WC$mat.CW / genes.final.filter.WC$total.CW
genes.final.filter.WB$bias.WB = genes.final.filter.WB$mat.WB / genes.final.filter.WB$total.WB
genes.final.filter.WB$bias.BW = genes.final.filter.WB$mat.BW / genes.final.filter.WB$total.BW
genes.final.filter.CB$bias.CB = genes.final.filter.CB$mat.CB / genes.final.filter.CB$total.CB
genes.final.filter.CB$bias.BC = genes.final.filter.CB$mat.BC / genes.final.filter.CB$total.BC

genes.final.filter.WC$p.value.WC = 0
genes.final.filter.WC$p.value.CW = 0
genes.final.filter.WB$p.value.WB = 0
genes.final.filter.WB$p.value.BW = 0
genes.final.filter.CB$p.value.CB = 0
genes.final.filter.CB$p.value.BC = 0

for(each_gene in 1:nrow(genes.final.filter.WC)){
  yy = binom.test(genes.final.filter.WC[each_gene,]$mat.WC,genes.final.filter.WC[each_gene,]$total.WC,0.5,alternative="two.sided")$p.value
  genes.final.filter.WC$p.value.WC[each_gene] = yy
}
for(each_gene in 1:nrow(genes.final.filter.WC)){
  yy = binom.test(genes.final.filter.WC[each_gene,]$mat.CW,genes.final.filter.WC[each_gene,]$total.CW,0.5,alternative="two.sided")$p.value
  genes.final.filter.WC$p.value.CW[each_gene] = yy
}

for(each_gene in 1:nrow(genes.final.filter.WB)){
  yy = binom.test(genes.final.filter.WB[each_gene,]$mat.WB,genes.final.filter.WB[each_gene,]$total.WB,0.5,alternative="two.sided")$p.value
  genes.final.filter.WB$p.value.WB[each_gene] = yy
}
for(each_gene in 1:nrow(genes.final.filter.WB)){
  yy = binom.test(genes.final.filter.WB[each_gene,]$mat.BW,genes.final.filter.WB[each_gene,]$total.BW,0.5,alternative="two.sided")$p.value
  genes.final.filter.WB$p.value.BW[each_gene] = yy
}

for(each_gene in 1:nrow(genes.final.filter.CB)){
  yy = binom.test(genes.final.filter.CB[each_gene,]$mat.CB,genes.final.filter.CB[each_gene,]$total.CB,0.5,alternative="two.sided")$p.value
  genes.final.filter.CB$p.value.CB[each_gene] = yy
}
for(each_gene in 1:nrow(genes.final.filter.CB)){
  yy = binom.test(genes.final.filter.CB[each_gene,]$mat.BC,genes.final.filter.CB[each_gene,]$total.BC,0.5,alternative="two.sided")$p.value
  genes.final.filter.CB$p.value.BC[each_gene] = yy
}

genes.final.filter.WC$bias.yn.WC = 0
genes.final.filter.WC$bias.yn.CW = 0
genes.final.filter.WB$bias.yn.WB = 0
genes.final.filter.WB$bias.yn.BW = 0
genes.final.filter.CB$bias.yn.CB = 0
genes.final.filter.CB$bias.yn.BC = 0

genes.final.filter.WC$bias.yn.WC = ifelse(genes.final.filter.WC$p.value.WC>(0.05/nrow(genes.final.filter.WC)), "NS", "S")
genes.final.filter.WC$bias.yn.CW = ifelse(genes.final.filter.WC$p.value.CW>(0.05/nrow(genes.final.filter.WC)), "NS", "S")
genes.final.filter.WB$bias.yn.WB = ifelse(genes.final.filter.WB$p.value.WB>(0.05/nrow(genes.final.filter.WB)), "NS", "S")
genes.final.filter.WB$bias.yn.BW = ifelse(genes.final.filter.WB$p.value.BW>(0.05/nrow(genes.final.filter.WB)), "NS", "S")
genes.final.filter.CB$bias.yn.CB = ifelse(genes.final.filter.CB$p.value.CB>(0.05/nrow(genes.final.filter.CB)), "NS", "S")
genes.final.filter.CB$bias.yn.BC = ifelse(genes.final.filter.CB$p.value.BC>(0.05/nrow(genes.final.filter.CB)), "NS", "S")

genes.final.filter.WC$bias.cat.WC = "biparental"
genes.final.filter.WC$bias.cat.WC[genes.final.filter.WC$bias.yn.WC=="S" & genes.final.filter.WC$bias.WC<=0.050]<-"paternal.only"
genes.final.filter.WC$bias.cat.WC[genes.final.filter.WC$bias.yn.WC=="S" & genes.final.filter.WC$bias.WC>0.050 & genes.final.filter.WC$bias.WC<0.350]<-"paternal.bias"
genes.final.filter.WC$bias.cat.WC[genes.final.filter.WC$bias.yn.WC=="S" & genes.final.filter.WC$bias.WC>=0.950]<-"maternal.only"
genes.final.filter.WC$bias.cat.WC[genes.final.filter.WC$bias.yn.WC=="S" & genes.final.filter.WC$bias.WC>0.650 & genes.final.filter.WC$bias.WC<0.950]<-"maternal.bias"

genes.final.filter.WC$bias.cat.CW = "biparental"
genes.final.filter.WC$bias.cat.CW[genes.final.filter.WC$bias.yn.CW=="S" & genes.final.filter.WC$bias.CW<=0.050]<-"paternal.only"
genes.final.filter.WC$bias.cat.CW[genes.final.filter.WC$bias.yn.CW=="S" & genes.final.filter.WC$bias.CW>0.050 & genes.final.filter.WC$bias.CW<0.350]<-"paternal.bias"
genes.final.filter.WC$bias.cat.CW[genes.final.filter.WC$bias.yn.CW=="S" & genes.final.filter.WC$bias.CW>=0.950]<-"maternal.only"
genes.final.filter.WC$bias.cat.CW[genes.final.filter.WC$bias.yn.CW=="S" & genes.final.filter.WC$bias.CW>0.650 & genes.final.filter.WC$bias.CW<0.950]<-"maternal.bias"

genes.final.filter.WB$bias.cat.WB = "biparental"
genes.final.filter.WB$bias.cat.WB[genes.final.filter.WB$bias.yn.WB=="S" & genes.final.filter.WB$bias.WB<=0.050]<-"paternal.only"
genes.final.filter.WB$bias.cat.WB[genes.final.filter.WB$bias.yn.WB=="S" & genes.final.filter.WB$bias.WB>0.050 & genes.final.filter.WB$bias.WB<0.350]<-"paternal.bias"
genes.final.filter.WB$bias.cat.WB[genes.final.filter.WB$bias.yn.WB=="S" & genes.final.filter.WB$bias.WB>=0.950]<-"maternal.only"
genes.final.filter.WB$bias.cat.WB[genes.final.filter.WB$bias.yn.WB=="S" & genes.final.filter.WB$bias.WB>0.650 & genes.final.filter.WB$bias.WB<0.950]<-"maternal.bias"

genes.final.filter.WB$bias.cat.BW = "biparental"
genes.final.filter.WB$bias.cat.BW[genes.final.filter.WB$bias.yn.BW=="S" & genes.final.filter.WB$bias.BW<=0.050]<-"paternal.only"
genes.final.filter.WB$bias.cat.BW[genes.final.filter.WB$bias.yn.BW=="S" & genes.final.filter.WB$bias.BW>0.050 & genes.final.filter.WB$bias.BW<0.350]<-"paternal.bias"
genes.final.filter.WB$bias.cat.BW[genes.final.filter.WB$bias.yn.BW=="S" & genes.final.filter.WB$bias.BW>=0.950]<-"maternal.only"
genes.final.filter.WB$bias.cat.BW[genes.final.filter.WB$bias.yn.BW=="S" & genes.final.filter.WB$bias.BW>0.650 & genes.final.filter.WB$bias.BW<0.950]<-"maternal.bias"

genes.final.filter.CB$bias.cat.CB = "biparental"
genes.final.filter.CB$bias.cat.CB[genes.final.filter.CB$bias.yn.CB=="S" & genes.final.filter.CB$bias.CB<=0.050]<-"paternal.only"
genes.final.filter.CB$bias.cat.CB[genes.final.filter.CB$bias.yn.CB=="S" & genes.final.filter.CB$bias.CB>0.050 & genes.final.filter.CB$bias.CB<0.350]<-"paternal.bias"
genes.final.filter.CB$bias.cat.CB[genes.final.filter.CB$bias.yn.CB=="S" & genes.final.filter.CB$bias.CB>=0.950]<-"maternal.only"
genes.final.filter.CB$bias.cat.CB[genes.final.filter.CB$bias.yn.CB=="S" & genes.final.filter.CB$bias.CB>0.650 & genes.final.filter.CB$bias.CB<0.950]<-"maternal.bias"

genes.final.filter.CB$bias.cat.BC = "biparental"
genes.final.filter.CB$bias.cat.BC[genes.final.filter.CB$bias.yn.BC=="S" & genes.final.filter.CB$bias.BC<=0.050]<-"paternal.only"
genes.final.filter.CB$bias.cat.BC[genes.final.filter.CB$bias.yn.BC=="S" & genes.final.filter.CB$bias.BC>0.050 & genes.final.filter.CB$bias.BC<0.350]<-"paternal.bias"
genes.final.filter.CB$bias.cat.BC[genes.final.filter.CB$bias.yn.BC=="S" & genes.final.filter.CB$bias.BC>=0.950]<-"maternal.only"
genes.final.filter.CB$bias.cat.BC[genes.final.filter.CB$bias.yn.BC=="S" & genes.final.filter.CB$bias.BC>0.650 & genes.final.filter.CB$bias.BC<0.950]<-"maternal.bias"

bias.WC.counts <- ddply(genes.final.filter.WC, c("bias.cat.WC"), summarise,N = length(gene), perc = length(gene)*100/nrow(genes.final.filter.WC))
bias.CW.counts <- ddply(genes.final.filter.WC, c("bias.cat.CW"), summarise,N = length(gene), perc = length(gene)*100/nrow(genes.final.filter.WC))
bias.WB.counts <- ddply(genes.final.filter.WB, c("bias.cat.WB"), summarise,N = length(gene), perc = length(gene)*100/nrow(genes.final.filter.WB))
bias.BW.counts <- ddply(genes.final.filter.WB, c("bias.cat.BW"), summarise,N = length(gene), perc = length(gene)*100/nrow(genes.final.filter.WB))
bias.CB.counts <- ddply(genes.final.filter.CB, c("bias.cat.CB"), summarise,N = length(gene), perc = length(gene)*100/nrow(genes.final.filter.CB))
bias.BC.counts <- ddply(genes.final.filter.CB, c("bias.cat.BC"), summarise,N = length(gene), perc = length(gene)*100/nrow(genes.final.filter.CB))

bias.WC.counts
bias.CW.counts
bias.WB.counts
bias.BW.counts
bias.CB.counts
bias.BC.counts

##### Incorporate reciprocal information

genes.final.filter.WC$bias <- (genes.final.filter.WC$bias.WC + genes.final.filter.WC$bias.CW)/2
genes.final.filter.WB$bias <- (genes.final.filter.WB$bias.WB + genes.final.filter.WB$bias.BW)/2
genes.final.filter.CB$bias <- (genes.final.filter.CB$bias.CB + genes.final.filter.CB$bias.BC)/2

genes.final.filter.WC$rec.diff <- abs(genes.final.filter.WC$bias.WC - genes.final.filter.WC$bias.CW)
genes.final.filter.WB$rec.diff <- abs(genes.final.filter.WB$bias.WB - genes.final.filter.WB$bias.BW)
genes.final.filter.CB$rec.diff <- abs(genes.final.filter.CB$bias.CB - genes.final.filter.CB$bias.BC)

mean(genes.final.filter.WC$rec.diff)
mean(genes.final.filter.WB$rec.diff)
mean(genes.final.filter.CB$rec.diff)

cor(genes.final.filter.WC$bias.WC, genes.final.filter.WC$bias.CW,method="spearman")
cor(genes.final.filter.WB$bias.WB, genes.final.filter.WB$bias.BW,method="spearman")
cor(genes.final.filter.CB$bias.CB, genes.final.filter.CB$bias.BC,method="spearman")

count(ifelse(genes.final.filter.WC$bias.cat.WC == genes.final.filter.WC$bias.cat.CW,"SAME","DIFFERENT"))
count(ifelse(genes.final.filter.WB$bias.cat.WB == genes.final.filter.WB$bias.cat.BW,"SAME","DIFFERENT"))
count(ifelse(genes.final.filter.CB$bias.cat.CB == genes.final.filter.CB$bias.cat.BC,"SAME","DIFFERENT"))

genes.final.filter.WC$concordance <- ifelse(((genes.final.filter.WC$bias.cat.WC == genes.final.filter.WC$bias.cat.CW) | (genes.final.filter.WC$rec.diff < mean(genes.final.filter.WC$rec.diff)*2)),"AGREE","DISAGREE")
genes.final.filter.WB$concordance <- ifelse(((genes.final.filter.WB$bias.cat.WB == genes.final.filter.WB$bias.cat.BW) | (genes.final.filter.WB$rec.diff < mean(genes.final.filter.WB$rec.diff)*2)),"AGREE","DISAGREE")
genes.final.filter.CB$concordance <- ifelse(((genes.final.filter.CB$bias.cat.CB == genes.final.filter.CB$bias.cat.BC) | (genes.final.filter.CB$rec.diff < mean(genes.final.filter.CB$rec.diff)*2)),"AGREE","DISAGREE")

genes.final.filter.WC$bias.cat = "DISAGREE"
genes.final.filter.WC$bias.cat[genes.final.filter.WC$concordance=="AGREE" & genes.final.filter.WC$bias<=0.050]<-"paternal.only"
genes.final.filter.WC$bias.cat[genes.final.filter.WC$concordance=="AGREE" & genes.final.filter.WC$bias>0.050 & genes.final.filter.WC$bias<0.350]<-"paternal.bias"
genes.final.filter.WC$bias.cat[genes.final.filter.WC$concordance=="AGREE" & genes.final.filter.WC$bias>=0.950]<-"maternal.only"
genes.final.filter.WC$bias.cat[genes.final.filter.WC$concordance=="AGREE" & genes.final.filter.WC$bias>0.650 & genes.final.filter.WC$bias<0.950]<-"maternal.bias"
genes.final.filter.WC$bias.cat[genes.final.filter.WC$concordance=="AGREE" & (genes.final.filter.WC$bias.cat.WC == "biparental" | genes.final.filter.WC$bias.cat.CW == "biparental")]<-"biparental"

genes.final.filter.WB$bias.cat = "DISAGREE"
genes.final.filter.WB$bias.cat[genes.final.filter.WB$concordance=="AGREE" & genes.final.filter.WB$bias<=0.050]<-"paternal.only"
genes.final.filter.WB$bias.cat[genes.final.filter.WB$concordance=="AGREE" & genes.final.filter.WB$bias>0.050 & genes.final.filter.WB$bias<0.350]<-"paternal.bias"
genes.final.filter.WB$bias.cat[genes.final.filter.WB$concordance=="AGREE" & genes.final.filter.WB$bias>=0.950]<-"maternal.only"
genes.final.filter.WB$bias.cat[genes.final.filter.WB$concordance=="AGREE" & genes.final.filter.WB$bias>0.650 & genes.final.filter.WB$bias<0.950]<-"maternal.bias"
genes.final.filter.WB$bias.cat[genes.final.filter.WB$concordance=="AGREE" & (genes.final.filter.WB$bias.cat.WB == "biparental" | genes.final.filter.WB$bias.cat.BW == "biparental")]<-"biparental"

genes.final.filter.CB$bias.cat = "DISAGREE"
genes.final.filter.CB$bias.cat[genes.final.filter.CB$concordance=="AGREE" & genes.final.filter.CB$bias<=0.050]<-"paternal.only"
genes.final.filter.CB$bias.cat[genes.final.filter.CB$concordance=="AGREE" & genes.final.filter.CB$bias>0.050 & genes.final.filter.CB$bias<0.350]<-"paternal.bias"
genes.final.filter.CB$bias.cat[genes.final.filter.CB$concordance=="AGREE" & genes.final.filter.CB$bias>=0.950]<-"maternal.only"
genes.final.filter.CB$bias.cat[genes.final.filter.CB$concordance=="AGREE" & genes.final.filter.CB$bias>0.650 & genes.final.filter.CB$bias<0.950]<-"maternal.bias"
genes.final.filter.CB$bias.cat[genes.final.filter.CB$concordance=="AGREE" & (genes.final.filter.CB$bias.cat.CB == "biparental" | genes.final.filter.CB$bias.cat.BC == "biparental")]<-"biparental"

count(genes.final.filter.WC$bias.cat)
count(genes.final.filter.WB$bias.cat)
count(genes.final.filter.CB$bias.cat)
#write.csv(genes.final.filter.WC, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_intra_WC.csv")
#write.csv(genes.final.filter.WB, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_intra_WB.csv")
#write.csv(genes.final.filter.CB, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_intra_CB.csv")

# plot

bias.WC.counts.plot <- ggplot(bias.WC.counts, aes(bias.cat.WC, perc)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge(),fill="plum4") +
  scale_x_discrete(limit = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"),
                   labels = c("M","MB","B","PB","P")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  labs(title="WC", y="%", x = "Category of bias") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme()

bias.CW.counts.plot <- ggplot(bias.CW.counts, aes(bias.cat.CW, perc)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge(),fill="plum2") +
  scale_x_discrete(limit = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"),
                   labels = c("M","MB","B","PB","P")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  labs(title="CW", y="%", x = "Category of bias") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme()

bias.WB.counts.plot <- ggplot(bias.WB.counts, aes(bias.cat.WB, perc)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge(),fill="cyan4") +
  scale_x_discrete(limit = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"),
                   labels = c("M","MB","B","PB","P")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  labs(title="WB", y="%", x = "Category of bias") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme()

bias.BW.counts.plot <- ggplot(bias.BW.counts, aes(bias.cat.BW, perc)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge(),fill="cyan3") +
  scale_x_discrete(limit = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"),
                   labels = c("M","MB","B","PB","P")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  labs(title="BW", y="%", x = "Category of bias") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme()

bias.CB.counts.plot <- ggplot(bias.CB.counts, aes(bias.cat.CB, perc)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge(),fill="gold4") +
  scale_x_discrete(limit = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"),
                   labels = c("M","MB","B","PB","P")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  labs(title="CB", y="%", x = "Category of bias") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme()

bias.BC.counts.plot <- ggplot(bias.BC.counts, aes(bias.cat.BC, perc)) + ylim(c(0,90)) +
  geom_bar(stat="identity", position=position_dodge(),fill="gold3") +
  scale_x_discrete(limit = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"),
                   labels = c("M","MB","B","PB","P")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  labs(title="BC", y="%", x = "Category of bias") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme()

#grid.arrange(bias.WC.counts.plot,bias.CW.counts.plot,bias.WB.counts.plot,bias.BW.counts.plot,bias.CB.counts.plot,bias.BC.counts.plot,ncol=2)

WxC.plot = ggplot(genes.final.filter.WC, aes(bias.WC,bias.CW)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x="WC bias", y="CW bias") +
  geom_rect(data=NULL,aes(xmin=0.65,xmax=1.01,ymin=0.65,ymax=1.01),
            fill="#FFE9ED",alpha=0.9,color="white") +
  geom_rect(data=NULL,aes(xmin=-0.01,xmax=0.35,ymin=-0.01,ymax=0.35),
            fill="#F3FDFF",color="white",alpha=0.9) +
  geom_vline(xintercept = 0.35,linetype="dotted", color = "gray40") +
  geom_vline(xintercept = 0.65,linetype="dotted", color = "gray40") +
  geom_hline(yintercept = 0.65,linetype="dotted", color = "gray40") +
  geom_hline(yintercept = 0.35,linetype="dotted", color = "gray40") + 
  geom_point(aes(colour=bias.cat),size=0.5) + 
  labs(title="WxC",x="WC bias", y="CW bias") +
  scale_colour_manual(breaks = c("biparental","maternal.bias","maternal.only","DISAGREE"),
                      labels = c("B", "MB", "M",""),values=c("magenta","lightcyan4","magenta3","magenta4")) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

WxB.plot = ggplot(genes.final.filter.WB, aes(bias.WB,bias.BW)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x="WB bias", y="BW bias") +
  geom_rect(data=NULL,aes(xmin=0.65,xmax=1.01,ymin=0.65,ymax=1.01),
            fill="#FFE9ED",alpha=0.9,color="white") +
  geom_rect(data=NULL,aes(xmin=-0.01,xmax=0.35,ymin=-0.01,ymax=0.35),
            fill="#F3FDFF",color="white",alpha=0.9) +
  geom_vline(xintercept = 0.35,linetype="dotted", color = "gray40") +
  geom_vline(xintercept = 0.65,linetype="dotted", color = "gray40") +
  geom_hline(yintercept = 0.65,linetype="dotted", color = "gray40") +
  geom_hline(yintercept = 0.35,linetype="dotted", color = "gray40") + 
  geom_point(aes(colour=bias.cat),size=0.5) + 
  labs(title="WxB",x="WC bias", y="CW bias") +
  scale_colour_manual(breaks = c("biparental","maternal.bias","maternal.only","DISAGREE"),
                      labels = c("B", "MB", "M",""),values=c("aquamarine3","lightcyan4","aquamarine4","darkslategray")) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

CxB.plot = ggplot(genes.final.filter.CB, aes(bias.CB,bias.BC)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x="CB bias", y="BC bias") +
  geom_rect(data=NULL,aes(xmin=0.65,xmax=1.01,ymin=0.65,ymax=1.01),
            fill="#FFE9ED",alpha=0.9,color="white") +
  geom_rect(data=NULL,aes(xmin=-0.01,xmax=0.35,ymin=-0.01,ymax=0.35),
            fill="#F3FDFF",color="white",alpha=0.9) +
  geom_vline(xintercept = 0.35,linetype="dotted", color = "gray40") +
  geom_vline(xintercept = 0.65,linetype="dotted", color = "gray40") +
  geom_hline(yintercept = 0.65,linetype="dotted", color = "gray40") +
  geom_hline(yintercept = 0.35,linetype="dotted", color = "gray40") + 
  geom_point(aes(colour=bias.cat),size=1) + 
  labs(title="WxB",x="WC bias", y="CW bias") +
  scale_colour_manual(breaks = c("biparental","maternal.bias","maternal.only","DISAGREE"),
                      labels = c("B", "MB", "M",""),values=c("darkgoldenrod2","lightcyan4","darkgoldenrod3","tan4")) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

genes.final.filter.WC.rec.counts = ddply(genes.final.filter.WC,c("bias.cat"),summarise,N = length(gene), perc = length(gene)*100/nrow(genes.final.filter.WC))
genes.final.filter.WB.rec.counts = ddply(genes.final.filter.WB,c("bias.cat"),summarise,N = length(gene), perc = length(gene)*100/nrow(genes.final.filter.WB))
genes.final.filter.CB.rec.counts = ddply(genes.final.filter.CB,c("bias.cat"),summarise,N = length(gene), perc = length(gene)*100/nrow(genes.final.filter.CB))

genes.final.filter.WC.rec.counts$bias.cat.order <- revalue(genes.final.filter.WC.rec.counts$bias.cat, c("biparental"="3B",
                                                                                                 "maternal.bias"="2MB",
                                                                                                 "maternal.only"="1M",
                                                                                                 "DISAGREE"="4D"))
genes.final.filter.WB.rec.counts$bias.cat.order <- revalue(genes.final.filter.WB.rec.counts$bias.cat, c("biparental"="3B",
                                                                                                        "maternal.bias"="2MB",
                                                                                                        "maternal.only"="1M",
                                                                                                        "DISAGREE"="4D"))
genes.final.filter.CB.rec.counts$bias.cat.order <- revalue(genes.final.filter.CB.rec.counts$bias.cat, c("biparental"="3B",
                                                                                                        "maternal.bias"="2MB",
                                                                                                        "maternal.only"="1M",
                                                                                                        "DISAGREE"="4D"))

genes.final.filter.WC.rec.counts <- genes.final.filter.WC.rec.counts[order(genes.final.filter.WC.rec.counts$bias.cat.order),]
genes.final.filter.WB.rec.counts <- genes.final.filter.WB.rec.counts[order(genes.final.filter.WB.rec.counts$bias.cat.order),]
genes.final.filter.CB.rec.counts <- genes.final.filter.CB.rec.counts[order(genes.final.filter.CB.rec.counts$bias.cat.order),]


WxC.plot2 = ggplot(genes.final.filter.WC.rec.counts,aes(x=bias.cat.order,y=perc,fill=bias.cat.order)) +
  geom_bar(stat="identity",colour="black") +
  scale_x_discrete(breaks = c("1M","2MB","3B","4D"), labels = c("M", "MB", "B","No POE bias")) +
  geom_text(aes(label=N), position = position_dodge(0.9),vjust = -1)  +
  scale_fill_manual(breaks = c("1M","2MB","3B","4D"), labels = c("M", "MB", "B",""),
                    values=c("magenta4","magenta3","magenta","lightcyan4")) +
  scale_y_continuous(limits=c(0,100)) +
  labs(title="WxC",y="%", x = "ASE category") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

WxB.plot2 = ggplot(genes.final.filter.WB.rec.counts,aes(x=bias.cat.order,y=perc,fill=bias.cat.order)) +
  geom_bar(stat="identity",colour="black") +
  scale_x_discrete(breaks = c("1M","2MB","3B","4D"), labels = c("M", "MB", "B","No POE bias")) +
  geom_text(aes(label=N), position = position_dodge(0.9),vjust = -1)  +
  scale_fill_manual(breaks = c("1M","2MB","3B","4D"), labels = c("M", "MB", "B",""),
                    values=c("darkslategray","aquamarine4","aquamarine3","lightcyan4")) +
  scale_y_continuous(limits=c(0,100)) +
  labs(title="WxB",y="%", x = "ASE category") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

CxB.plot2 = ggplot(genes.final.filter.CB.rec.counts,aes(x=bias.cat.order,y=perc,fill=bias.cat.order)) +
  geom_bar(stat="identity",colour="black") +
  scale_x_discrete(breaks = c("1M","2MB","3B","4D"), labels = c("M", "MB", "B","No POE bias")) +
  geom_text(aes(label=N), position = position_dodge(0.9),vjust = -1)  +
  scale_fill_manual(breaks = c("1M","2MB","3B","4D"), labels = c("M", "MB", "B",""),
                    values=c("tan4","darkgoldenrod3","darkgoldenrod2","lightcyan4")) +
  scale_y_continuous(limits=c(0,100)) +
  labs(title="CxB",y="%", x = "ASE category") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

aa <- grid.arrange(WC.bias.SNP.hist,CW.bias.SNP.hist,ncol=1)
bb <- grid.arrange(WB.bias.SNP.hist,BW.bias.SNP.hist,ncol=1)
cc <- grid.arrange(CB.bias.SNP.hist,BC.bias.SNP.hist,ncol=1)

# supplinfo

#png("/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/manuscript/supplinfo.intra.genotypes.counts.jpg",
     #width = 3800, height = 3400, units = 'px', res = 300)
#grid.arrange(bias.WC.counts.plot,bias.CW.counts.plot,bias.WB.counts.plot,bias.BW.counts.plot,bias.CB.counts.plot,bias.BC.counts.plot,ncol=2)
#dev.off()

# hybrids vs intraspecific: a comparison

genes.final.filter.WC
genes.final.filter.WB
genes.final.filter.CB

# import hybrid datasets

genes_final_filter_hybrid_soma <- read_csv("Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_soma.csv")
genes_final_filter_hybrid_testis <- read_csv("Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_testis.csv")

genes.hybrid.soma <- genes_final_filter_hybrid_soma[c("gene","bias","bias.cat")]
genes.hybrid.testis <- genes_final_filter_hybrid_testis[c("gene","bias","bias.cat")]

genes.WxC <- genes.final.filter.WC[c("gene","bias","bias.cat")]
colnames(genes.WxC) = c("gene","bias.WxC","bias.cat.WxC")
genes.WxB <- genes.final.filter.WB[c("gene","bias","bias.cat")]
colnames(genes.WxB) = c("gene","bias.WxB","bias.cat.WxB")
genes.CxB <- genes.final.filter.CB[c("gene","bias","bias.cat")]
colnames(genes.CxB) = c("gene","bias.CxB","bias.cat.CxB")

genes.intra.0 <- full_join(genes.WxC,genes.WxB,by="gene")
genes.intra <- full_join(genes.intra.0,genes.CxB,by="gene")

# average pm across genotypes

genes.intra$num.gen = ifelse(is.na(genes.intra$bias.WxC),0,1) + ifelse(is.na(genes.intra$bias.WxB),0,1) + ifelse(is.na(genes.intra$bias.CxB),0,1)
genes.intra$bias.intra = round(((na2zero(genes.intra$bias.WxC)+na2zero(genes.intra$bias.WxB)+na2zero(genes.intra$bias.CxB))/genes.intra$num.gen),3)

genes.intra.biparental <- subset(genes.intra, bias.cat.WxC == "biparental" | bias.cat.WxB == "biparental" | bias.cat.CxB == "biparental")
genes.intra.biparental$bias.cat.intra <- "biparental"
genes.intra.disagree <- subset(genes.intra, (bias.cat.WxC == "DISAGREE" | is.na(bias.cat.WxC)) &
                                 (bias.cat.WxB == "DISAGREE" | is.na(bias.cat.WxB)) &
                                 (bias.cat.CxB == "DISAGREE" | is.na(bias.cat.CxB)))
genes.intra.disagree$bias.cat.intra <- "DISAGREE"

genes.intra.maternal <- anti_join(genes.intra,rbind(genes.intra.biparental,genes.intra.disagree),by=c("gene"))
genes.intra.maternal$bias.cat.intra <- ifelse(genes.intra.maternal$bias.intra >= 0.95,"maternal.only","maternal.bias")

genes.intra.combined <- rbind(genes.intra.biparental,genes.intra.maternal,genes.intra.disagree)
#write.csv(genes.intra.combined, file = "/Users/delafilia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes.intra.combined.csv")

# cross soma and testis-only genes with intraspecific genes

genes.hybrid.soma.only <- anti_join(genes.hybrid.soma,genes.hybrid.testis,by="gene")
genes.hybrid.testis.only <- anti_join(genes.hybrid.testis,genes.hybrid.soma,by="gene")

genes.hybrid.soma.intra <- inner_join(genes.hybrid.soma.only,genes.intra.combined,by="gene")
genes.hybrid.testis.intra <- inner_join(genes.hybrid.testis.only,genes.intra.combined,by="gene")

nrow(genes.hybrid.soma.intra)
nrow(genes.hybrid.testis.intra)

# check correlation
cor(genes.hybrid.soma.intra$bias, genes.hybrid.soma.intra$bias.intra, method=c("spearman"))
cor(genes.hybrid.testis.intra$bias, genes.hybrid.testis.intra$bias.intra, method=c("spearman"))

genes.hybrid.soma.intra$p.diff <- genes.hybrid.soma.intra$bias - genes.hybrid.soma.intra$bias.intra
genes.hybrid.testis.intra$p.diff <- genes.hybrid.testis.intra$bias - genes.hybrid.testis.intra$bias.intra

sd(genes.hybrid.testis.intra$p.diff)*100

hist.diff.soma.intra <- ggplot(data=genes.hybrid.soma.intra, aes(p.diff)) +
  geom_rect(data=NULL,aes(xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf),
            fill="snow1",alpha=0.9,color="white") +
  geom_rect(data=NULL,aes(xmin=0,xmax=-Inf,ymin=-Inf,ymax=Inf),
            fill="snow2",color="white",alpha=0.9) + 
  geom_vline(xintercept = 0,color = "black") +
  geom_histogram(breaks=seq(-1,1, by=0.05), fill=c('orange1')) +
  labs(title="",y="Gene counts", x = expression(paste("Hybrid soma ",p[m]," - intraspecific ",p[m]))) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

hist.diff.testis.intra <- ggplot(data=genes.hybrid.testis.intra, aes(p.diff)) +
  geom_rect(data=NULL,aes(xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf),
            fill="snow1",alpha=0.9,color="white") +
  geom_rect(data=NULL,aes(xmin=0,xmax=-Inf,ymin=-Inf,ymax=Inf),
            fill="snow2",color="white",alpha=0.9) +
  geom_vline(xintercept = 0,color = "black") +
  geom_histogram(breaks=seq(-1,1, by=0.05), fill=c('steelblue4')) + ylim(0,100) +
  labs(title="",y="Gene counts", x = expression(paste("Hybrid testis ",p[m]," - intraspecific ",p[m]))) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

ggplot(genes.hybrid.soma.intra,aes(bias.intra,bias)) + geom_point() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
ggplot(genes.hybrid.testis.intra,aes(bias.intra,bias)) + geom_point() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
