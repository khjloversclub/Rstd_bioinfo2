# 계속된 커밋 오류로 튜토리얼 1장부터 다시 실행, 대용량 파일은 선택적으로 걸러낼 것

# WD 설정
getwd()
setwd("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2")
setwd("/home/kimhyunjin/Rstd_khj/GitHub_khj/Rstd_bioinfo2/mmdata1")

# 데이터 가져오기
library(Rsubread)
fastq.files <- list.files(path = "./mmdata1", pattern = ".fastq.gz$", full.names = TRUE)
fastq.files

# 인덱스 구축(생략), 바로 reference genome의 염색체1 read alignment)
align(index="mmdata1/chr1_mm10",readfile1=fastq.files) 
# 위 명령어가 웹페이지 상에 누락돼있음, with solutions에서만 확인 가능
args(align)
bam.files <- list.files(path = "./mmdata1", pattern = ".BAM$", full.names = TRUE)
bam.files
props <- propmapped(files=bam.files)
props

## Challenge 1
bam.files.multi <- gsub(".fastq.gz", ".multi.bam", as.character(fastq.files))
align(index="mmdata1/chr1_mm10", readfile1=fastq.files, unique = FALSE, nBestLocations = 6, output_format="BAM", output_file=bam.files.multi)
dir("./mmdata1", pattern="*.multi.bam")
props.multimap <- propmapped(files=bam.files.multi)
props.multimap

## QC
# Extract quality scores
qs <- qualityScores(filename="mmdata1/SRR1552450.fastq.gz",nreads=100)
# Check dimension of qs
dim(qs)
# Check first few elements of qs with head
head(qs)
boxplot(qs)

##Challenge 2
qs.51 <- qualityScores(filename="mmdata1/SRR1552451.fastq.gz",nreads=50)
dim(qs.51)
head(qs.51)
boxplot(qs.51)

# Counting
fc <- featureCounts(bam.files, annot.inbuilt="mm10")
names(fc)
## Take a look at the featurecounts stats
fc$stat
## Take a look at the dimensions to see the number of genes
dim(fc$counts)
## Take a look at the first 6 lines
head(fc$counts)
#show annotation
head(fc$annotation)

# Challenge 3
fc.exon <- featureCounts(bam.files, annot.inbuilt="mm10", useMetaFeatures = FALSE)
dim(fc.exon.multimap$counts)
head(fc.exon.multimap$counts)
fc.exon.multimap <- featureCounts(bam.files.multi, annot.inbuilt="mm10",useMetaFeatures = FALSE,countMultiMappingReads = TRUE)
