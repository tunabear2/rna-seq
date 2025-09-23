# rna-seq
sample 정보
샘플종 : Human , 조직종류 : Chang(HeLa)

Bulk RNA-seq 진행
TdT vs TdT Transwell, prom1_TdT vs prom1_TdT Transwell

샘플당 3개
> 실습 시 3개 sample로 진행하기.

merge나 concatenate 할지 말지 고려하기.

pipeline
1. align reference 준비
Genome FASTA: GRCh38.p14 , GENCODE human v49 (GRCh38.p14) GTF 사용
2. 논문에 들어가야 하니 FastQC 진행 후 원시 데이터 품질 정보 수집
3. Align 전에 FASTQ 데이터 필터링 . fastp or trimmomatic 중에 하나 선택하기. (Hard filtering 시에는 trimmomatic 아니면 fastp)
4. 필터링 이후 fastq 데이터 FastQC 진행.
5. read alignment 
STAR or HISAT2 사용 . STAR는 Suffix Array + Hash 기반 . HISAT2는 Burrows-Wheeler Transform (BWT)기반
작업 환경의 RAM 메모리가 여유로우면 STAR(30~40Gb) 아닐 때는 HISAT2 (10Gb 내외)
6. FeatureCounts 진행 
gene_count.csv 공통 파일 생성

여기까지가 공통 전처리 정
