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

##### import all expression counts (pure_WYE, hybrid, intraspecific)

PC_M_1_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/PC_M_1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
PC_M_2_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/PC_M_2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
PC_M_3_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/PC_M_3.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)

S1_A_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S1_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S1_A_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S1_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S1_B_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S1_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S1_B_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S1_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S3_A_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S3_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S3_A_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S3_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S3_B_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S3_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S3_B_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S3_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S6_A_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S6_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S6_A_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S6_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S6_B_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S6_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
S6_B_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/S6_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T1_A_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T1_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T1_A_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T1_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T1_B_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T1_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T1_B_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T1_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T3_A_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T3_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T3_A_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T3_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T3_B_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T3_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T3_B_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T3_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T6_A_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T6_A_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T6_A_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T6_A_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T6_B_6_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T6_B_6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
T6_B_7_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/T6_B_7.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)

WC01_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WC1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WC04_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WC4.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WC06_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WC6.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WC19_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WC19.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CW01_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CW1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CW04_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CW4.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CW10_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CW10.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CW11_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CW11.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WB02_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WB2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WB03_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WB3.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
WB08_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/WB8.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BW02_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BW2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BW03_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BW3.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BW05_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BW5.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BW12_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BW12.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CB01_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CB1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CB04_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CB4.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CB12_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CB12.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
CB20_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/CB20.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BC03_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BC3.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BC04_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BC4.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BC13_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BC13.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BC17_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/intraspecific/3_rsem/rsem_output/BC17.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)

##### get averages

# pure WYE

wye_tpm = PC_M_1_genes_results[c(1,6)]
colnames(wye_tpm)[colnames(wye_tpm)=="TPM"] <- "PC_M_1"
wye_tpm$PC_M_2 = PC_M_2_genes_results$TPM
wye_tpm$PC_M_3 = PC_M_3_genes_results$TPM
colnames(wye_tpm)[colnames(wye_tpm)=="gene_id"] <- "gene"
wye_tpm$wye_TPM <- rowMeans(wye_tpm[,-1])

# hybrids

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

s_tpm$s_TPM <- rowMeans(s_tpm[,-1])
t_tpm$t_TPM <- rowMeans(t_tpm[,-1])
colnames(s_tpm)[colnames(s_tpm)=="gene_id"] <- "gene"
colnames(t_tpm)[colnames(t_tpm)=="gene_id"] <- "gene"

# intra

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

intra_tpm$mean.tpm.WC <- rowMeans(intra_tpm[,2:5])
intra_tpm$mean.tpm.CW <- rowMeans(intra_tpm[,6:9])
intra_tpm$mean.tpm.WB <- rowMeans(intra_tpm[,10:13])
intra_tpm$mean.tpm.BW <- rowMeans(intra_tpm[,14:16])
intra_tpm$mean.tpm.CB <- rowMeans(intra_tpm[,17:19])
intra_tpm$mean.tpm.BC <- rowMeans(intra_tpm[,20:23])
intra_tpm$intra_TPM <- rowMeans(intra_tpm[,2:23])

colnames(intra_tpm)[colnames(intra_tpm)=="gene_id"] <- "gene"

##### merge

nrow(wye_tpm)
nrow(s_tpm)
nrow(t_tpm)
nrow(intra_tpm)

wye_tpm_merge <- wye_tpm[c("gene","wye_TPM")]
s_tpm_merge <- s_tpm[c("gene","s_TPM")]
t_tpm_merge <- t_tpm[c("gene","t_TPM")]

expression_avg0 <- merge(wye_tpm_merge,s_tpm_merge,by="gene")
expression_avg00 <- merge(expression_avg0,t_tpm_merge,by="gene")
expression_avg <- merge(expression_avg00,intra_tpm,by="gene")

##### examine correlation between TPM

cor(expression_avg$wye_TPM, expression_avg$s_TPM, method=c("spearman"))
cor(expression_avg$wye_TPM, expression_avg$t_TPM, method=c("spearman"))
cor(expression_avg$wye_TPM, expression_avg$intra_TPM, method=c("spearman"))

##### import ASE results

hybrid.genes.soma <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_soma.csv", ",", escape_double = FALSE, trim_ws = TRUE)
hybrid.genes.testis <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_testis.csv", ",", escape_double = FALSE, trim_ws = TRUE)
intra.genes <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes.intra.combined.csv", ",", escape_double = FALSE, trim_ws = TRUE)

hybrid.genes.soma.TPM <- merge(hybrid.genes.soma,expression_avg,by="gene")
hybrid.genes.testis.TPM <- merge(hybrid.genes.testis,expression_avg,by="gene")
intra.genes.TPM <- left_join(intra.genes,expression_avg,by="gene")
nrow(intra.genes.TPM)

# build dataframe for intraspecific model

intra.genes.TPM <- intra.genes.TPM[c("gene","bias.WxC","bias.cat.WxC","bias.WxB","bias.cat.WxB","bias.CxB","bias.cat.CxB","bias.intra","bias.cat.intra","wye_TPM","intra_TPM","s_TPM","t_TPM")]
hist(log10(intra.genes.TPM$wye_TPM))

model1 <- lm(log10(wye_TPM) ~ bias.intra, data = intra.genes.TPM)
summary(model1)
par(mfrow=c(2,2))
plot(model1)
AIC(model1)

ggplot(intra.genes.TPM, aes(x = bias.intra, y = log10(intra.genes.TPM$wye_TPM))) + 
  geom_point() +
  stat_smooth(method = "lm")

intra.genes.TPM$bias.intra2 <- intra.genes.TPM$bias.intra^2
model2 <- lm(log10(wye_TPM) ~ bias.intra + bias.intra2, data = intra.genes.TPM)
par(mfrow=c(2,2))
anova(model2)
summary(model2)
anova(model2,model1,test="Chisq")

nrow(intra.genes.TPM)

supplinfo.tmp.vs.bias <- ggplot(intra.genes.TPM, aes(x = bias.intra, y = log10(intra.genes.TPM$wye_TPM))) + 
  geom_point(size=0.5,width = 0.2) +
  labs(title="", y="TPM (log10)", x = expression(paste(p[m]))) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2) , col = "red") + 
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))
intra.genes.TPM.no.POE <- intra.genes.TPM[intra.genes.TPM$bias.cat.intra != "DISAGREE",]
nrow(intra.genes.TPM.no.POE)
nrow(intra.genes.TPM)

supplinfo.tmp.vs.bias.cat <- ggplot(intra.genes.TPM.no.POE, aes(x=bias.cat.intra, y=log10(wye_TPM)))  +
  geom_jitter(size=0.5,width = 0.25,alpha=1,aes(colour=bias.cat.intra)) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  scale_x_discrete("ASE category",limits=c("biparental","maternal.bias","maternal.only"),
                   labels=c(" B "," MB "," M ")) +
  scale_fill_manual(name="ASE category", breaks = c("biparental","maternal.bias","maternal.only"),
                    labels = c("B", "MB", "M"),
                    values = c("gray60","gray70","gray80"),guide=FALSE) +
  scale_colour_manual(name="ASE category", breaks = c("biparental","maternal.bias","maternal.only"),
                      labels = c("B", "MB", "M"),
                      values = c("gray60","gray70","gray80"),guide=FALSE)  +
  labs(title="",x="ASE expression", y ="log10 (TPM)") + 
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

grid.arrange(supplinfo.tmp.vs.bias,supplinfo.tmp.vs.bias.cat,ncol=1)

#png("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/manuscript/supplinfo.tmp.vs.bias.jpg",
    #width = 1800, height = 1800, units = 'px', res = 300)
#supplinfo.tmp.vs.bias
#dev.off()

png("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/manuscript/supplinfo.tmp.vs.bias.cat.revision.jpg",
width = 1400, height = 2800, units = 'px', res = 300)
grid.arrange(supplinfo.tmp.vs.bias,supplinfo.tmp.vs.bias.cat,ncol=1)
dev.off()
