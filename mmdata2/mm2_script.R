## Pre-processsing RNA-seq data

# WD 설정
getwd()
setwd("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2")
setwd("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2/mmdata2")

# data import
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# Read the sample information into R (분석할 샘플 12가지의 정보가 담긴 시트)
sampleinfo <- read.delim("mmdata2/SampleInfo.txt")
View(sampleinfo)
sampleinfo

# Reading in the count data
# Read the data into R
seqdata <- read.delim("mmdata2/GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
head(seqdata)
View(seqdata)
dim(seqdata)
# 각 행은 Entrez gene ID(유전자별로 등록된 고유의 식별번호, Entrez는 NCBI가 운영하는 유전 서열 데이터베이스)의 정보를 제공.

# Data format
# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
View(countdata)
head(countdata)
# samples' name
colnames(countdata)
# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
substr("ThisIsAString", start=1, stop=5)
colnames(countdata) <- substr(colnames(countdata), 1, 7)
View(countdata)
# e.g. MCL1.DG_BC2CTUACXX_ACTTGA_L002_R1 -> MCL1.DG
table(colnames(countdata)==sampleinfo$SampleName)

# Filter some lowly expressed genes
# Obtain CPMs (counts-per-million)
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)
#omitted 
col1sum <- sum(countdata[,1])/1000000
countdata[1,1]/col1sum
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)
# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)
# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

# challenge 1
plot(myCPM[,2],countdata[,2], xlab="CPM", ylab="Raw Count", main=colnames(myCPM)[2], ylim=c(0,50), xlim=c(0,3))
abline(v=0.5, h=10, col="blue")
