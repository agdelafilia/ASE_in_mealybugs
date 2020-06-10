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

# here we will take existing Dn/Ds estimates for genes, parse them with the sex-bias and parent-of-origin datasets to explore the interplay between paternal-genome silencing, sex-specificity, and evolutionary rates
setwd("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/all_methods/andrew")

# citsex is the sex specific expression data, 3 reps of males, 3 of females
citsex <- read.csv("P_citri_FPKMs_trimmed.csv",header=T,stringsAsFactors=F)

# import dn/ds matrix
pcit2 <- read.csv("p_citri_matrix_v2.5.csv", stringsAsFactors=F)

#normalize, add the reps and divide by the number to get average male or female expression
f.adj<-rowSums(citsex[,2:4])/3
m.adj<-rowSums(citsex[,5:7])/3

# SPM. it's a metric from 0 to 1 with 0 = no expression in the focal sex, 1 = 100% in the focal sex
# from female perspective, female specific (1) to male (0)
# the SPM formula is focal^2 / (sum of squares of all sets)
f.spm<-(f.adj^2)/((m.adj^2)+(f.adj^2))
ann.spm<-as.data.frame(cbind(citsex[,1],as.numeric(f.spm)))
colnames(ann.spm)<-c("Gene","SPM.Female")
#now we have a per-gene list of sex-specificity, stored as ann.spm

#plot sex-specificity for suppl. info
hist(f.spm, main = substitute(paste("Distribution of sex-specific expression in ",italic("P. citri"))),
     xlab="Male (0) to female-specific (1)", ylab="")
                              

# we'll merge this specifity data with the Dn/Ds data
cit3<-merge(pcit2,ann.spm,by="Gene")
cit3$SPM.Female<-as.double(as.character(cit3$SPM.Female))

# combine this with parent-of-origin expression in intraspecific males
finalcit<-read.csv("/Users/agarcia/Documents/genomics/RNA_seq_projects/mealybugs/results/freeze_nov19/genes.intra.combined.csv", stringsAsFactors=F)
colnames(finalcit)[1]<-"Gene"

# remove DISAGREE genes (no true POE pattern)
cit4 <- finalcit[finalcit$bias.cat.intra != "DISAGREE",]
count(cit4$bias.cat.intra)
cit5<-merge(cit3,cit4,by="Gene")
nrow(cit5)

# a glm and the quadratic term to ask about the relationship between sex-specificity and parent-of-origin-bias
plot(cit5$SPM.Female,cit5$bias.intra)
cit5$bias.intra2 <- cit5$bias.intra^2
te.binomial <- glm(SPM.Female ~ bias.intra,family="binomial",data=cit5)
te.binomial2 <- glm(SPM.Female ~ bias.intra + bias.intra2,family="binomial",data=cit5)

summary(te.binomial2)
anova(te.binomial,te.binomial2,test="LRT")

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

# now let's ask how Dn/Ds changes between categories
cit5<-cit5[which(cit5$bias.cat.intra!="DISAGREE"),]
cit5$bias.cat.intra<-as.factor(cit5$bias.cat.intra)

# we'll use a Kruskal-Wallis test (non-parametric ANOVA equivalent) as we have obvious categories (parent-bias) and a continuous response that is not likely normally distributed (Dn/Ds).
kruskal.test(dN.dS~bias.cat.intra, data=cit5)

#to do posthoc-testing, we'll use a Nemenyi test from the PMCMR package
library("PMCMR")
posthoc.kruskal.nemenyi.test(dN.dS ~ bias.cat.intra, data=cit5)

# Ds changes between categories in hybrids

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
