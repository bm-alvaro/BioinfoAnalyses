# Alvaro Benitez, 2023
# This R script is performing a Differential Expression analysis in RNAseq data. Parameter optimization for specific projects is mandatory. 

library(limma)
library(Glimma)
library(edgeR)
library(here)

# Basic transformations
readcount <- read.csv(here("readcount.csv"))
row.names(readcount) <- readcount$gene_id
readcount$gene_id <- NULL
A <- readcount[,c(1:6,19:24)]
d0 <- DGEList(A)

class(d0)
dim(d0)

# One could trim names with this code
smplnames <- substring(colnames(d0), 1, as.numeric(nchar(colnames(d0))-7))
smplnames

colnames(d0) <- smplnames

d0$samples
d0$samples$group

# Factor setup
timepoint <- as.factor(c("dpi0","dpi10","dpi0","dpi10","dpi0","dpi10","dpi0","dpi10","dpi0","dpi10","dpi0","dpi10"))
timepoint
levels(timepoint)

genotype <- as.factor(c("MUC","MUC","MUC","MUC","MUC","MUC","PERS","PERS","PERS","PERS","PERS","PERS"))
genotype

timegt <- paste(timepoint, genotype, sep = "_")
timegt

d0$samples$timepoint <- timepoint
d0$samples$genotype <- genotype
d0$samples$timegt <- timegt

d0
d0$samples

#A CPM value of 1 for a gene equates to having 15 counts in the sample with the lowest sequencing depth 
#(T7_1, library size approx. 15 million) or 20 counts in the sample with the greatest sequencing depth (T14_2, library size approx. 20 million). 
#The log-CPM values will be used for exploratory plots. When log=TRUE, the cpm function adds an offset to the CPM values before converting to the log2-scale. 
#By default, the offset is 2/L where 2 is the prior count and L is the average library size in millions, 
#so the log-CPM values are related to the CPM values by log2(CPM + 2/L). 
#This calculation ensures that any two read counts with identical CPM values will also have identical log-CPM values. 
#The prior count avoids taking the logarithm of zero, and also reduces spurious variability for genes with very low counts 
#by shrinking all the inter-sample log-fold-changes towards zero, something that is helpful for exploratory plotting. 
#For this dataset, the average library size is about 17.23 million, so L approx. 17.23 and the minimum log-CPM value for each sample becomes log2(2/17.2328) = -3.107085. 
#In other words, a count of zero for this data maps to a log-CPM value of -3.11 after adding the prior count or offset

cpm_d0 <- cpm(d0)
lcpm_d0 <- cpm(d0, log=TRUE)

L <- mean(d0$samples$lib.size) * 1e-6
M <- median(d0$samples$lib.size) * 1e-6

summary(lcpm_d0)

table(rowSums(d0$counts==0)==15)

# filterByExpr keeps genes that have a CPM of (10/M) or more in at least *six* samples (number of samples is determined by the smallest group)
keep.exprs <- filterByExpr(d0, group=timepoint)
d01 <- d0[keep.exprs,, keep.lib.sizes=FALSE]
dim(d01)

############ Density of log-CPM values for raw pre-filtered data and post-filtered data ############ 
cContrast <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4", # purple
  "#FF7F00", # orange
  "black", "gold1", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "maroon", "deeppink1", "blue1",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4"
)
pie(rep(1, 18), col = cMaria)

lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(d0)
col <- cContrast
par(mfrow=c(1,2))
plot(density(lcpm_d0[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm_d0[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", smplnames, text.col=col, bty="n")
lcpm_d01 <- cpm(d01, log=TRUE)
plot(density(lcpm_d01[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm_d01[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", smplnames, text.col=col, bty="n")

############ Calculate normalization factors ############ 

d01 <- calcNormFactors(d01, method = "TMM")
d01$samples$norm.factors
d01$samples

glMDSPlot(lcpm_d01, labels=paste(genotype, timepoint, timegt, sep="_"), 
          groups=d01$samples[,c(4,5,6)], launch=TRUE)

############ Design matrix and contrasts ############ 

design <- model.matrix(~0+timegt)
is.fullrank(design)
colnames(design) <- gsub("timegt", "", colnames(design))
design

contr.matrix <- makeContrasts(
  dpi10_M_Mvsdpi10_MUC = dpi10_M_M - dpi10_MUC,
  dpi10_M_Pvsdpi10_MUC = dpi10_M_P - dpi10_MUC,
  levels = colnames(design))
contr.matrix

contrast_matrix_to_save <- data.frame(contr.matrix)
write.table(contrast_matrix_to_save, file = here("contrast_matrix_B.txt"), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

############ Removing heteroscedascity from count data ############ 

par(mfrow=c(1,2))
v <- voom(d01, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))

tfit <- treat(vfit, lfc=log2(1.5))
dt <- decideTests(tfit)
summary(dt)

library(dplyr)
cpepo_gene_descr_ncbi <- read.delim("~/resources/cpepo_gene_descr_ncbi.txt", header=FALSE)
cpepo_gene_descr_ncbi$V1 <- as.character(cpepo_gene_descr_ncbi$V1)

dpi10_MUCvsdpi0_MUC <- topTreat(tfit, coef=1, n=Inf)
dpi10_MUCvsdpi0_MUC <- left_join(dpi10_MUCvsdpi0_MUC, cpepo_gene_descr_ncbi, by = c("id"="V1"))
dpi10_MUCvsdpi0_MUC_p <- subset(dpi10_MUCvsdpi0_MUC, dpi10_MUCvsdpi0_MUC$adj.P.Val < 0.05)
write.table(dpi10_MUCvsdpi0_MUC_p, file = here("dpi10_MUCvsdpi0_MUC_p.csv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(dpi10_MUCvsdpi0_MUC, file = here("dpi10_MUCvsdpi0_MUC_without_filter.csv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=lcpm_d01, groups=treattime, launch=TRUE)