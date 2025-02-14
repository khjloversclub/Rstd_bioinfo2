## Pre-processsing RNA-seq data

# WD 설정
getwd()
setwd("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2")
list.files("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2/mmdata2")

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

# Filter some lowly expressed genes(발현량 낮은 유전자 제거.)
# Obtain CPMs (counts-per-million)
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)
#omitted 
col1sum <- sum(countdata[,1])/1000000
countdata[1,1]/col1sum
# Which values in myCPM are greater than 0.5?
#참고 논문에서는 cpm값이 1 미만일 경우 미발현 유전자로 간주하여 거름.
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




## convert the counts to DGEList object(count data를 저장하는 기능)
dgeObj <- DGEList(counts.keep) # 객체 만들기.
# have a look at dgeObj
dgeObj
# See what slots are stored in dgeObj
names(dgeObj)
# Look at the structure of the list
str(dgeObj)
# Library size information is stored in the samples slot
dgeObj$samples



# QC
## confirm; library sizes & distribution plots
dgeObj$samples$lib.size
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)
#### rep_len(width, NR)에서 다음과 같은 에러가 발생했습니다: NULL을 길이가 0인 객체에 복제할 수 없습니다
########################################뭐가문제일까/ 0214.18:56 지금은 또 괜찮아졌네.

#추가정보: 경고메시지(들):
#  mean.default(width)에서:
#  인자가 수치형 또는 논리형이 아니므로 NA를 반환합니다
# Add a title to the plot
title("Barplot of library sizes")

# Get log2 counts per million; log2 scale에서 read count 분포 확인
logcounts <- cpm(dgeObj,log=TRUE)
# Check  density distributions of raw log-intensities of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
# 파란 수평선으로부터 멀리 떨어진 샘플은 특별히 더 조사할 필요 있음.
## MCL1.LA, MCL1.LB, MCL1.LE, MCL1.LF 4가지가 가장 눈에 띔.





# multi-dimensional scaling plot(MDSplot): 주성분 분석을 시각화하여 데이터에서 가장 큰 변인을 파악.
plotMDS(dgeObj)
# 그룹화 정보에 따라 샘플 분류(색상, 점, label 변경 등)
# We specify the option to let us plot two plots side-by-side
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$CellType) ####### NULL #######

######## troubleshooting with chatgpt #########
head(sampleinfo)  # 데이터 미리보기
str(sampleinfo)   # 데이터 구조 확인
colnames(sampleinfo)  # 컬럼 이름 확인

#CellType의 factor 여부; character or NULL 출력 > factor로 변환
class(sampleinfo$CellType) # character 출력
sampleinfo$CellType <- as.factor(sampleinfo$CellType)
levels(sampleinfo$CellType)  # 다시 확인
################# factor로 변환 성공 #################

## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)
col.cell
sampleinfo$CellType

# Redo the MDS with cell type colouring #################################
plotMDS(dgeObj,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
# Add a title
title("Cell type")

# Similarly for status ##트슛 타임##
levels(sampleinfo$Status)
col.status <- c("blue","red","dark green")[sampleinfo$Status]
col.status
plotMDS(dgeObj,col=col.status)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")
######################범례가0으로나옴#######################
#Status를 factor로 변경
class(sampleinfo$Status) #character 출력
sampleinfo$Status <- as.factor(sampleinfo$Status)
levels(sampleinfo$Status)  # 다시 확인
# 위의 코드 재실행/ 문제 해결 0214.19:42

# Q:올바른 위치에 없는 것으로 보이는 두 샘플을 식별하세요.
## 보이지 않는 두 샘플이 조사대상일 것.