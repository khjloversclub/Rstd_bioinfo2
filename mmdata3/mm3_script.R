### Differential Expression of RNA-seq data(RNA-seq 데이터의 차등적 발현)

# WD 설정
getwd()
setwd("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2")
list.files("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2/mmdata3")

## Differential expression with edgeR
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
load("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2/mmdata2/preprocessing.Rdata")

# 차등 발현 분석을 위한 간단한 전처리 과정
## Read the counts from the downloaded data
seqdata <- read.delim("mmdata2/GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
countdata
colnames(countdata) <- substr(colnames(countdata), 1, 7)
countdata

View(countdata)
## Calculate the Counts Per Million measure
myCPM <- cpm(countdata)
## Identify genes with at least 0.5 cpm in at least 2 samples
thresh <- myCPM > 0.5
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
## Convert to an edgeR object
dgeObj <- DGEList(counts.keep)
## Perform TMM normalisation
dgeObj <- calcNormFactors(dgeObj)
## Obtain corrected sample information
sampleinfo <- read.delim("mmdata2/SampleInfo_Corrected.txt")
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group


### Create the design matrix(변수 CellType, Status로 상호 작용이 없는 모델과 있는 모델을 적용해 행렬 생성)
# Create the two variables
group <- as.character(group)
type <- sapply(strsplit(group, ".", fixed=T), function(x) x[1])
status <- sapply(strsplit(group, ".", fixed=T), function(x) x[2])
# Specify a design matrix with an intercept term
design <- model.matrix(~ type + status)
design
# 데이터 탐색
plotMDS(dgeObj, labels=group, cex=0.75, xlim=c(-4, 5))
#분산 추정(공통 분산은 모든 유전자에 대한 평균; dataset 전체의 BCV* 추정)
dgeObj <- estimateCommonDisp(dgeObj)
#유전자별 분산 추정치를 추정하여 평균 계수 크기에 따른 가능한 추세를 허용함.
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)
plotBCV(dgeObj)

## Testing for differential expression
# Fit the linear model
fit <- glmFit(dgeObj, design)
names(fit)
head(coef(fit))
#Conduct likelihood ratio tests for luminal vs basal and show the top genes
# likelihood-ratio test(공분산 검증)
lrt.BvsL <- glmLRT(fit, coef=2)
topTags(lrt.BvsL)

# Contrast (To find the differentially expressed genes btw Pre & Vir)
PvsV <- makeContrasts(statuspregnant-statusvirgin, levels=design) #(Intercept)가 Intercept로 변경된다는 메시지 출력, 크게 중요하지 않음
lrt.pVsV <- glmLRT(fit, contrast=PvsV)
topTags(lrt.pVsV)

save(lrt.BvsL,dgeObj,group,file="mmdata3/DE.Rdata")
