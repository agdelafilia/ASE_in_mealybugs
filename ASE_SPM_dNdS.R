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
library(PMCMR)

#here we will take existing Dn/Ds estimates for genes, parse them with the sex-bias and parent-of-origin datasets to explore the interplay between paternal-genome silencing, sex-specificity, and evolutionary rates
setwd("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/all_methods/andrew")

# citsex is the sex specific expression data, 3 reps of males, 3 of females, as you'll see below I'll use these data to calculate SPM the specificity metric
citsex <- read.csv("P_citri_FPKMs_trimmed.csv",header=T,stringsAsFactors=F)

# use new matrix Andrew sent me; the new matrix does not contain duplicated genes (Andrew is sendind dn/ds for the longest trasncripts)
pcit2 <- read.csv("p_citri_matrix_v2.5.csv", stringsAsFactors=F)

# let's generate info about sex specificity
# Let's check tho: correlations are consistently high so this shouldn't be a problem
cor(citsex[-c(1)],method=c("spearman"))

#normalize, add the reps and divide by the number to get average male or female expression
f.adj<-rowSums(citsex[,2:4])/3
m.adj<-rowSums(citsex[,5:7])/3

#now let's do SPM. it's a metric from 0 to 1 with 0 = no expression in the focal sex, 1 = 100% in the focal sex
#from female perspective, female specific (1) to male (0)
# the SPM formula is focal^2 / (sum of squares of all sets)
f.spm<-(f.adj^2)/((m.adj^2)+(f.adj^2))
ann.spm<-as.data.frame(cbind(citsex[,1],as.numeric(f.spm)))
colnames(ann.spm)<-c("Gene","SPM.Female")
#now we have a per-gene list of sex-specificity, stored as ann.spm

#plot sex-specificity for suppl. info
hist(f.spm, main = substitute(paste("Distribution of sex-specific expression in ",italic("P. citri"))),
     xlab="Male (0) to female-specific (1)", ylab="")
                              
#png("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/manuscript/supplinfo.sex.sp.jpg",
#    width = 1800, height = 1800, units = 'px', res = 300)
#hist(f.spm, main = substitute(paste("Distribution of sex-specific expression in ",italic("P. citri"))),
#     xlab="Male (0) to female-specific (1)", ylab="")
#dev.off()
#note Male SPM is just the compliment of F.SPM, and is redundant with only 2 sexes

#we'll merge this specifity data with the Dn/Ds data
cit3<-merge(pcit2,ann.spm,by="Gene")

#because R is weird with data-types, we need to ensure that the SPM numbers are being correctly interpreted as double precision numbers
cit3$SPM.Female<-as.double(as.character(cit3$SPM.Female))

#now we want to combine this with Andres's final pass of parent-of-origin expression 
finalcit<-read.csv("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes.intra.combined.csv", stringsAsFactors=F)
?predict.glm
#we need the "gene" column to be capitialized to exactly match for merging
colnames(finalcit)[1]<-"Gene"

# remove DISAGREE genes
cit4 <- finalcit[finalcit$bias.cat.intra != "DISAGREE",]
count(cit4$bias.cat.intra)
cit5<-merge(cit3,cit4,by="Gene")
nrow(cit5)

# import TPM counts for pure citri males (we have FPKM above but for consistency, but let's use TPM for gene expression levels)
PC_M_1_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/PC_M_1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
PC_M_2_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/PC_M_2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
PC_M_3_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/PC_M_3.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
wye_tpm = PC_M_1_genes_results[c(1,6)]
colnames(wye_tpm)[colnames(wye_tpm)=="TPM"] <- "PC_M_1"
wye_tpm$PC_M_2 = PC_M_2_genes_results$TPM
wye_tpm$PC_M_3 = PC_M_3_genes_results$TPM
colnames(wye_tpm)[colnames(wye_tpm)=="gene_id"] <- "gene"
wye_tpm$wye_TPM <- rowMeans(wye_tpm[,-1])
wye_tpm<-wye_tpm[c(1,5)]
colnames(wye_tpm)[1] <- "Gene"
colnames(wye_tpm)[2] <- "TPM"
cit5 <- merge(cit5,wye_tpm)

# import total read counts from intraspecific males to use as weights for the binomial glm
WC <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_intra_WC.csv", ",", escape_double = FALSE, trim_ws = TRUE)
WB <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_intra_WB.csv", ",", escape_double = FALSE, trim_ws = TRUE)
CB <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_intra_CB.csv", ",", escape_double = FALSE, trim_ws = TRUE)

WC <- WC[c("gene","total.WC","total.CW")]
WB <- WB[c("gene","total.WB","total.BW")]
CB <- CB[c("gene","total.CB","total.BC")]

WC$weight.WC <- WC$total.WC + WC$total.CW
WB$weight.WB <- WB$total.WB + WB$total.BW
CB$weight.CB <- CB$total.CB + CB$total.BC

intra <- cit5[c("Gene","bias.cat.WxC","bias.cat.WxB","bias.cat.CxB","num.gen")]
colnames(intra) <- c("gene","WxC","WxB","CxB","num.gen")
intra <- left_join(intra,WC[c(1,4)],by="gene")
intra <- left_join(intra,WB[c(1,4)],by="gene")
intra <- left_join(intra,CB[c(1,4)],by="gene")

intra$weight.WC <- ifelse(is.na(intra$weight.WC),0,intra$weight.WC)
intra$weight.WB <- ifelse(is.na(intra$weight.WB),0,intra$weight.WB)
intra$weight.CB <- ifelse(is.na(intra$weight.CB),0,intra$weight.CB)

intra$weight <- intra$weight.WC + intra$weight.WB + intra$weight.CB
intra.weights <- intra[c("gene","weight")]
colnames(intra.weights)[1] <- "Gene"
cit5 <- merge(cit5,intra.weights)

# a glm and the quadratic term to ask about the relationship between sex-specificity and parent-of-origin-bias
plot(cit5$SPM.Female,cit5$bias.intra)
cit5$bias.intra2 <- cit5$bias.intra^2
te.binomial2 <- glm(SPM.Female ~ bias.intra + bias.intra2,family="quasibinomial",weights=weight,data=cit5)

summary(te.binomial2)
anova(te.binomial2)
anova(te.binomial2,test="F")

#library("rsq")
#rsq(te.binomial2)
glm.plot <- ggplot(cit5, aes(x = bias.intra, y = SPM.Female)) + 
  geom_point(aes(colour=SPM.Female),size=0.5) +
  scale_colour_gradient2(low="deepskyblue2", mid = "darkgrey", high="firebrick2",midpoint=0.5) +
  labs(title="", y="SPM", x = expression(paste(p[m]))) +
  stat_smooth(method = "glm", formula = y ~ x + I(x^2), col = "firebrick4") +
  coord_cartesian(ylim=c(0,1)) + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))
nrow(cit5)

# is the hump-shape due to background expression of female-specific genes in males?

# import TPM counts for pure citri males (we have FPKM above but for consistency let's use TPM for gene expression levels)
PC_M_1_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/PC_M_1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
PC_M_2_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/PC_M_2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
PC_M_3_genes_results <- read_delim("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/hybrid/3_rsem/rsem_output/PC_M_3.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
wye_tpm = PC_M_1_genes_results[c(1,6)]
colnames(wye_tpm)[colnames(wye_tpm)=="TPM"] <- "PC_M_1"
wye_tpm$PC_M_2 = PC_M_2_genes_results$TPM
wye_tpm$PC_M_3 = PC_M_3_genes_results$TPM
colnames(wye_tpm)[colnames(wye_tpm)=="gene_id"] <- "gene"
wye_tpm$wye_TPM <- rowMeans(wye_tpm[,-1])
wye_tpm<-wye_tpm[c(1,5)]
colnames(wye_tpm)[1] <- "Gene"
colnames(wye_tpm)[2] <- "TPM"

cit5 <- merge(cit5,wye_tpm)
spm.tpm <- ggplot(cit5, aes(y = log10(TPM), x = SPM.Female)) + 
  geom_point(aes(colour=SPM.Female),size=0.5) +
  scale_colour_gradient2(low="deepskyblue2", mid = "darkgrey", high="firebrick2",midpoint=0.5) +
  labs(title="", x="SPM", y = "TPM (log10)") +
  geom_smooth(col = "firebrick4", se = F) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

# let's check if the relationship holds with only moderately expressed genes and above

cit5.cutoff<-cit5[cit5$TPM >= 10,]
nrow(cit5.cutoff)
te.binomial2.10 <- glm(SPM.Female ~ bias.intra + bias.intra2,family="quasibinomial",weights=weight,data=cit5.cutoff)
summary(te.binomial2.10)
anova(te.binomial2.10)
anova(te.binomial2.10,test="F")

cutoff10 <- ggplot(cit5.cutoff, aes(x = bias.intra, y = SPM.Female)) + 
  geom_point(aes(colour=SPM.Female),size=0.5) +
  scale_colour_gradient2(low="deepskyblue2", mid = "darkgrey", high="firebrick2",midpoint=0.5) +
  labs(title="TPM >= 10 (1,721 genes)", y="SPM", x = expression(paste(p[m]))) +
  stat_smooth(method = "glm", formula = y ~ x + I(x^2), col = "firebrick4") +
  coord_cartesian(ylim=c(0,1)) + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

cit5.cutoff<-cit5[cit5$TPM >= 100,]
nrow(cit5.cutoff)
te.binomial2.100 <- glm(SPM.Female ~ bias.intra + bias.intra2,family="quasibinomial",weights=weight,data=cit5.cutoff)
summary(te.binomial2.100)
anova(te.binomial2.100)
anova(te.binomial2.100,test="F")

cutoff100 <- ggplot(cit5.cutoff, aes(x = bias.intra, y = SPM.Female)) + 
  geom_point(aes(colour=SPM.Female),size=0.5) +
  scale_colour_gradient2(low="deepskyblue2", mid = "darkgrey", high="firebrick2",midpoint=0.5) +
  labs(title="TPM >= 100 (267 genes)", y="SPM", x = expression(paste(p[m]))) +
  stat_smooth(method = "glm", formula = y ~ x + I(x^2), col = "firebrick4") +
  coord_cartesian(ylim=c(0,1)) + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

library(patchwork)
fig.s9 <- spm.tpm / (cutoff10 | cutoff100)

#Now let's ask how Dn/Ds changes between categories
cit5<-cit5[which(cit5$bias.cat.intra!="DISAGREE"),]
cit5$bias.cat.intra<-as.factor(cit5$bias.cat.intra)

#we'll use a Kruskal-Wallis test (non-parametric ANOVA equivalent) as we have obvious categories (parent-bias) and a continuous response that is not likely normally distributed (Dn/Ds).
kruskal.test(dN.dS~bias.cat.intra, data=cit5)

#to do posthoc-testing, we'll use a Nemenyi test from the PMCMR package
library("PMCMR")
posthoc.kruskal.nemenyi.test(dN.dS ~ bias.cat.intra, data=cit5)

# ANDRES: Now let's look at how Dn/Ds changes between categories in hybrids

scit<-read.csv("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_soma.csv", stringsAsFactors=F)
tcit<-read.csv("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_testis.csv", stringsAsFactors=F)

#we need the "gene" column to be capitialized to exactly match for merging
colnames(scit)[1]<-"Gene"
colnames(tcit)[1]<-"Gene"

# merge with cit3
scit.dnds0<-merge(cit3,scit,by="Gene")
tcit.dnds0<-merge(cit3,tcit,by="Gene")

#keep only B, MB, M

scit.dnds <- scit.dnds0[(scit.dnds0$bias.cat == "biparental" | scit.dnds0$bias.cat == "maternal.bias" | scit.dnds0$bias.cat == "maternal.only"),]
tcit.dnds <- tcit.dnds0[(tcit.dnds0$bias.cat == "biparental" | tcit.dnds0$bias.cat == "maternal.bias" | tcit.dnds0$bias.cat == "maternal.only"),]

#Now let's ask how Dn/Ds changes between categories
scit.dnds$bias.cat<-as.factor(scit.dnds$bias.cat)
tcit.dnds$bias.cat<-as.factor(tcit.dnds$bias.cat)

kruskal.test(dN.dS~bias.cat, data=scit.dnds)
posthoc.kruskal.nemenyi.test(dN.dS ~ bias.cat, data=scit.dnds)

kruskal.test(dN.dS~bias.cat, data=tcit.dnds)
posthoc.kruskal.nemenyi.test(dN.dS ~ bias.cat, data=tcit.dnds)

# plot log10(dn/ds)

dnds.vs.bias.intra <- ggplot(data=cit5, aes(x=bias.cat.intra, y=log10(dN.dS+0.001)))  +
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
  geom_text(x=1.5,y=0.36,   label="4e-4" ,size=3) +
  geom_text(x=2,y=0.56,   label="2e-3"  ,size=3) +
  geom_text(x=2.5,y=0.76, label="1" ,size=3) +
  geom_segment(aes(x=1,xend=2,y=0.30,yend=0.30),size=0.1) +
  geom_segment(aes(x=1,xend=3,y=0.50,yend=0.50),size=0.1) +
  geom_segment(aes(x=2,xend=3,y=0.70,yend=0.70),size=0.1) +
  labs(title="",x="ASE category", y ="log10 (dN/dS)") + 
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

dnds.vs.bias.soma <- ggplot(data=scit.dnds, aes(x=bias.cat, y=log10(dN.dS+0.001)))  +
  geom_jitter(size=0.5,width = 0.25,alpha=1,aes(colour=bias.cat)) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  scale_x_discrete("ASE category",limits=c("biparental","maternal.bias","maternal.only"),
                   labels=c(" B "," MB "," M ")) +
  scale_fill_manual(name="ASE category", breaks = c("biparental","maternal.bias","maternal.only"),
                    labels = c("B", "MB", "M"),guide=FALSE) +
  scale_colour_manual(name="ASE category", breaks = c("biparental","maternal.bias","maternal.only"),
                      labels = c("B", "MB", "M"),
                      values = c("#B9770E","#F8C471","#FAD7A0"),guide=FALSE)  +
  geom_text(x=1.5,y=0.36,   label="2e-5" ,size=3) +
  geom_text(x=2,y=0.56,   label="0.17"  ,size=3) +
  geom_text(x=2.5,y=0.76, label="1e-6" ,size=3) +
  geom_segment(aes(x=1,xend=2,y=0.30,yend=0.30),size=0.1) +
  geom_segment(aes(x=1,xend=3,y=0.50,yend=0.50),size=0.1) +
  geom_segment(aes(x=2,xend=3,y=0.70,yend=0.70),size=0.1) +
  labs(title="",x="ASE category", y ="log10 (dN/dS)") + 
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

dnds.vs.bias.testis <- ggplot(data=tcit.dnds, aes(x=bias.cat, y=log10(dN.dS+0.001)))  +
  geom_jitter(size=0.5,width = 0.25,alpha=1,aes(colour=bias.cat)) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  scale_x_discrete("ASE category",limits=c("biparental","maternal.bias","maternal.only"),
                   labels=c(" B "," MB "," M ")) +
  scale_fill_manual(name="ASE category", breaks = c("biparental","maternal.bias","maternal.only"),
                    labels = c("B", "MB", "M"),guide=FALSE) +
  scale_colour_manual(name="ASE category", breaks = c("biparental","maternal.bias","maternal.only"),
                      labels = c("B", "MB", "M"),
                      values = c("#2980B9","#85C1E9","#AED6F1"),guide=FALSE) +
  geom_text(x=1.5,y=0.36,   label="0.29" ,size=3) +
  geom_text(x=2,y=0.56,   label="0.25"  ,size=3) +
  geom_text(x=2.5,y=0.76, label="0.92" ,size=3) +
  geom_segment(aes(x=1,xend=2,y=0.30,yend=0.30),size=0.1) +
  geom_segment(aes(x=1,xend=3,y=0.50,yend=0.50),size=0.1) +
  geom_segment(aes(x=2,xend=3,y=0.70,yend=0.70),size=0.1) +
  labs(title="",x="ASE expression", y ="log10 (dN/dS)") + 
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + #theme_classic()
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + theme(legend.position="none") +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

#boxplot(dN.dS ~ bias.cat, data=scit.dnds,outline=F,notch=T,ylab="Dn/Ds",las=1)

#tiff("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/manuscript/fig3updated.tiff",
#    width = 2500, height = 2500, units = 'px', res = 300)
#grid.arrange(glm.plot,dnds.vs.bias.intra,dnds.vs.bias.soma,dnds.vs.bias.testis,ncol=2)
#dev.off()

# how do genes cluster in syntheny groups with solenopsis?

# import Andrew's df with gene-to-syntheny group info

syntheny <- read.csv("p_citri_matrix_v4.5.csv",header=T,stringsAsFactors=F)
# how many do we have?
count(syntheny$Chr)

# group non-chr scaffolds
#syntheny$Chr <- ifelse(is.na(syntheny$Chr), "unplaced",syntheny$Chr)
#count(syntheny$Chr)
#
colnames(syntheny)[1]<-"gene"

# import hybrids and intra
hybrid.soma <- read_csv("~/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_soma.csv")
hybrid.testis <- read_csv("~/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes_final_filter_hybrid_testis.csv")
intra <- read_csv("~/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes.intra.combined.csv")

s.chr <- left_join(hybrid.soma[c(1,19)],syntheny[c(1,9)],by="gene")
t.chr <- left_join(hybrid.testis[c(1,19)],syntheny[c(1,9)],by="gene")
i.chr <- left_join(intra[c(1,10)],syntheny[c(1,9)],by="gene")

s.chr$bias.cat.group <- ifelse(s.chr$bias.cat == "paternal.bias" | s.chr$bias.cat == "paternal.only" | s.chr$bias.cat == "biparental", "biparental.extended", s.chr$bias.cat)
t.chr$bias.cat.group <- ifelse(t.chr$bias.cat == "paternal.bias" | t.chr$bias.cat == "paternal.only" | t.chr$bias.cat == "biparental", "biparental.extended", t.chr$bias.cat)


s.chr.counts = ddply(s.chr, c("Chr","bias.cat.group"), summarise,
                                              N = length(gene))

t.chr.counts = ddply(t.chr, c("Chr","bias.cat.group"), summarise,
                     N = length(gene))

i.chr.counts = ddply(i.chr, c("Chr","bias.cat.intra"), summarise,
                     N = length(gene))

s.chr.counts$bias.cat.group <- factor(s.chr.counts$bias.cat.group,levels = c("maternal.only","maternal.bias","biparental.extended"))
t.chr.counts$bias.cat.group <- factor(t.chr.counts$bias.cat.group,levels = c("maternal.only","maternal.bias","biparental.extended"))
i.chr.counts$bias.cat.intra <- factor(i.chr.counts$bias.cat.intra,levels = c("maternal.only","maternal.bias","biparental","DISAGREE"))

s.chr.plot <- ggplot(s.chr.counts, aes(fill=bias.cat, y=N, x=Chr)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="CF soma",x=expression(paste("Synteny group with ", italic("P. solenopsis"))), y ="Gene count") +
  scale_fill_manual(name="ASE category",values=c("gray50", "gray60", "gray70", "gray80", 
                            "gray90"), breaks=c("maternal.only","maternal.bias","biparental","paternal.bias","paternal.only"),
                    labels=c("M", "MB", "B","PB","P")) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

t.chr.plot <- ggplot(t.chr.counts, aes(fill=bias.cat, y=N, x=Chr)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="CF reproductive tract",x=expression(paste("Synteny group with ", italic("P. solenopsis"))), y ="Gene count") +
  scale_fill_manual(name="ASE category",values=c("gray50", "gray60", "gray70", "gray80", 
                                                 "gray90"), breaks=c("maternal.only","maternal.bias","biparental","paternal.bias","paternal.only"),
                    labels=c("M", "MB", "B","PB","P")) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

i.chr.plot <- ggplot(i.chr.counts, aes(fill=bias.cat.intra, y=N, x=Chr)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="Intraspecific males",x=expression(paste("Synteny group with ", italic("P. solenopsis"))), y ="Gene count") +
  scale_fill_manual(name="ASE category",values=c("gray50", "gray60", "gray70", "gray95"), breaks=c("maternal.only","maternal.bias","biparental","DISAGREE"),
                    labels=c("M", "MB", "B","No POE bias")) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme(panel.grid.minor = element_line(colour = "black",size=0),panel.grid.major = element_line(colour = "white",size=0))

#jpeg("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/manuscript/ASE_solenopsis.jpg",
#    width = 2500, height = 2500, units = 'px', res = 300)
grid.arrange(s.chr.plot,t.chr.plot,i.chr.plot,ncol=1)
#dev.off()
