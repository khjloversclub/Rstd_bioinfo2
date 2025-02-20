### Annotation and Visualisation of RNA-seq results
# RNA-seq 결과 주석 및 시각화

# WD 설정
getwd()
setwd("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2")
list.files("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2/mmdata4")

#이전 과정인 DE analysis 관련 객체가 있는지 확인.
suppressPackageStartupMessages(library(edgeR))
load("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2/mmdata3/DE.Rdata") #또는 해당 파일 클릭하여 열기.

# Entrez Gene ID 주석 확인(현재 유일하게 확인 가능한 주석)
results <- as.data.frame(topTags(lrt.BvsL,n = Inf))
results
dim(results)
#plotSmear로 MA plot과 유사한 방식으로 DE 유전자 강조.
library(edgeR)
#summary(de <- decideTestsDGE(lrt.BvsL))

####오류 발생####
#h(simpleError(msg, call))에서 다음과 같은 에러가 발생했습니다: 
#함수 'summary'를 위한 메소드 선택시 인수 'object'를 평가하는데 오류가 발생했습니다: 함수 "decideTestsDGE"를 찾을 수 없습니다

class(lrt.BvsL) #객체 확인; "DGELRT" 출력
summary(de <- decideTests.DGELRT(lrt.BvsL)) #DGELRT로 변경
detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.BvsL, de.tags=detags)

### Adding the annotation to the edgeR result
#org.Mm.eg.db 패키지로 주석 처리하기기
#### source("http://www.bioconductor.org/biocLite.R")    #Page not found
# BiocLite 자체가 없어진 것 같음. 대체할 것을 찾아보기.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db", force = TRUE)
# org.Mm.eg.db만 일단 다운로드

##### BiocLite가 없으므로 무시
#biocLite("org.Mm.eg.db")
# For Human
#biocLite("org.Hs.eg.db")
#####

library(org.Mm.eg.db) #어떤 정보를 원하는지 결정
columns(org.Mm.eg.db) #무엇을 추출할 수 있는지 확인
keytypes(org.Mm.eg.db) #원하는 정보 추출을 위해 key(의 집합으)로 DB 필터링
keys(org.Mm.eg.db, keytype="ENTREZID")[1:10]  #ENTREZID에서 사용할 key 유형 확인

## Build up the query step-by-step(사용할 key의 유효성 검사)
my.keys <- c("50916", "110308","12293")
my.keys %in% keys(org.Mm.eg.db, keytype="ENTREZID")
all(my.keys %in% keys(org.Mm.eg.db, keytype="ENTREZID"))

# 단계별 쿼리 작성 ############################
## to be filled-in interactively during the class.
ann <- select(org.Mm.eg.db, keys = rownames(results), columns = c("ENTREZID", "SYMBOL", "GENENAME"))
# Have a look at the annotation
ann
# symbol은 각 유전자마다 문자형 ID를 부여했다고 이해하면 될 것 같음.
table(ann$ENTREZID==rownames(results))
# 주석 정보를 데이터프레임에 binding
results.annotated <- cbind(results, ann)
results.annotated
#CSV 파일로 저장
write.csv(results.annotated,file="B.PregVsLacResults.csv",row.names=FALSE)

#challenge1: Re-visit the `plotSmear` plot from above and use the `text` function to add labels for the names of the top 200 most DE genes
plotSmear(lrt.BvsL, de.tags=detags)
N <- 200
text(results.annotated$logCPM[1:N],results.annotated$logFC[1:N],labels = results.annotated$SYMBOL[1:N],col="blue")
# text 함수를 사용하여 상위 200개 유전자 레이블 추가하여 플롯 재생성.

# Alternative; Volcano plot(두 그룹 사이에서 DEG를 효과적으로 시각화하는 그래프)
signif <- -log10(results.annotated$FDR)
plot(results.annotated$logFC, signif, pch=16)
points(results.annotated[detags, 'logFC'], -log10(results.annotated[detags, 'FDR']), pch=16, col='red')

#개별 샘플의 발현 수준 살펴보기 by stripchart.
library(RColorBrewer)
par(mfrow=c(1,3))
normCounts <- dgeObj$counts
# Let's look at the first gene in the topTable, Krt5, which has a rowname 50916
stripchart(normCounts["110308",]~group)
# This plot is ugly, let's make it better
stripchart(normCounts["110308",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,col=1:6,method="jitter")
# Let's use nicer colours
nice.col <- brewer.pal(6,name="Dark2")
stripchart(normCounts["110308",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="Krt5")

################glXYPlot; 상호작용 패널 설정##################
library(Glimma)
group2 <- group
levels(group2) <- c("basal.lactate","basal.preg","basal.virgin","lum.lactate", "lum.preg", "lum.virgin")

glXYPlot(x=results$logFC, y=-log10(results$FDR),
         xlab="logFC", ylab="B", main="B.PregVsLacResults",
         counts=normCounts, groups=group2, status=de,
         anno=ann, id.column="ENTREZID", folder="volcano")
######### Second argument should contain the first ######### 포기



## Retrieving Genomic Locations(게놈 위치 검색)
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", force = TRUE) #설치 완료.
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
tx <- TxDb.Mmusculus.UCSC.mm10.knownGene
columns(tx)

#Challenge 1: find the exon coordinates of genes, 50916, 110308, 12293,  using TxDb~ package
keysss<-c('50916', '110308', '12293')
select(tx, keys = keysss, keytype = 'GENEID', columns = c("EXONCHROM","EXONSTART","EXONEND"))

## Overview of GenomicRanges(GR)
library(GenomicRanges)
simple.range <-GRanges("1", IRanges(start=1000,end=2000)) #불특정한 염색체의 시작, 끝의 위치를 지정.
simple.range
chrs <- c("chr13", "chr15","chr5")
start <- c(73000000, 101000000, 15000000)
end <- c(74000000, 102000000, 16000000)
my.ranges <- GRanges(rep(chrs,3), IRanges(start=rep(start,each=3), end=rep(end,each=3)))
# 두 개의 GR 사이의 겹치는 영역 빠르게 식별(단, 객체별 명명 규칙에 주의)
keysss<-c('50916', '110308', '12293')
genePos <- select(tx, keys = keysss, keytype = 'GENEID', columns = c("EXONCHROM","EXONSTART","EXONEND"))
geneRanges <- GRanges(genePos$EXONCHROM, IRanges(genePos$EXONSTART,genePos$EXONEND), GENEID=genePos$GENEID)
geneRanges
findOverlaps(my.ranges,geneRanges)
seqlevelsStyle(geneRanges)
seqlevelsStyle(simple.range)

# 특정 유전자의 구조에 접근하기
exo <- exonsBy(tx,'gene')
range(exo[['50916']])
exo[["50916"]]

##########메모리 정리 한 번씩 해주기###########
rm(tx)            # 기존 객체 제거
gc()              # 메모리 정리
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
tx <- TxDb.Mmusculus.UCSC.mm10.knownGene # 객체 다시 불러오기
##########메모리 정리 한 번씩 해주기###########

## Exporting tracks(트랙 내보내기)

