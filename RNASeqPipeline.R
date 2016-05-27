####################
# DSCS6020 Final project
# Author: Angeline Madrid
# Filename: RNASeqPipeline.R
# Revision date: August 20, 2015
####################
# Purpose: This is a pipeline for loading, preprocessing, and
# analyzing raw RNA-Seq data produced from a Next Generation
# Sequencing experiment performed on a sample of RNA. All steps
# in the pipeline use the Bioconductor software tool for R
# (http://www.bioconductor.org/)
####################
# Outline of pipeline steps:
# 1. Obtain data: download SRA file from GEO website (http://www.ncbi.nlm.nih.gov/gds/),
#     load SRA file, and convert to FASTQ format
# 2. Check quality of data
# 3. Remove adapter sequences
# 4. Filter/trim data based on quality (Phred scores)
# 5. Align data to reference genome/transcriptome
# 6. Annotate transcripts and count mapped read for transcript abundance
####################
####################
# use Bioconductor
source("http://bioconductor.org/biocLite.R")

####################
# 1. Obtain data: download SRA file from GEO website (http://www.ncbi.nlm.nih.gov/gds/)
#     load SRA file, and convert to FASTQ format
####################
biocLite("SRAdb")
library(SRAdb)
sqlfile <- 'SRAmetadb.sqlite'
# it takes a while to load the database
if(!file.exists('SRAmetadb.sqlite')) {
  sra_con <- dbConnect(SQLite(),sqlfile)
  sqlfile<-getSRAdbFile()
}
# to test existance of database
#dbListTables(sra_con)
#dbListFields(sra_con, "study")

# list addresses of FASTQ files associated with SRX015776 (corn experiment #)
# and get ftp address
rs<-listSRAfile( c("SRX015776"), sra_con, fileType="sra")
# download sra data file from ftp site, downloaded files are saved to current
# working directory
sra_file<-getSRAfile( in_acc=c("SRR034098"), sra_con=sra_con, destDir=getwd(),
                      fileType = 'sra',srcType='ftp')

# convert from sra to FASTQ format
# this uses the SRA toolkit http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
# the new fASTQ file is saved to the current working directory
system("/Applications/sratoolkit.2.5.2-mac64/bin/fastq-dump --split-3 SRR034098.sra")

####################
# 2. Check quality of data
####################
biocLite("qrqc")
library(qrqc)
fq_file<-"SRR034098.fastq"
# creates and returns a FASTQSummary object
s.fastq<-readSeqFile(fq_file, type="fastq",max.length=1000)
plot_qual<-plotQuals(s.fastq)
plot_freq_bases<-plotBases(s.fastq)

####################
####################
# END OF OBTAIN DATA
####################
####################
# 3. Remove adapter sequences
####################
biocLite("girafe")
biocLite("ShortRead")
library(girafe)
library(ShortRead)
# load the fastq file
# creates and returns a ShortReadQ object
# readFastq() is from ShortRead pkg
reads_raw<-readFastq(fq_file)
# trim the adapters
# the 3' end adapter is found at 
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM433622
adapter<-"TCGTATGCCGTCTTCTGCTTG"
# trimAdapter() is from girafe pkg
reads_na<-trimAdapter(reads_raw, adapter)

# filter the 5' end
# GTTCAGAGTTCTACAGTCCGACGATC
adapter_2<-"GTTCAGAGTTCTACAGTCCGACGATC"
# trimLRPatterns() from ShortRead pkg
reads_na_2<-trimLRPatterns(Rpattern=adapter_2,subject=reads_na)

outfile_2<-"SRR034098_na_2.fastq"
# writeFastq() is from ShortRead pkg
writeFastq(reads_na_2,outfile_2)
####################
# 4. Filter/trim data based on quality (Phred scores)
####################
biocLite("seqTools")
library(seqTools)
# discarded sequences with quality < 20 and length < 15
num_reads_3_trimQ<-trimFastq("SRR034098_na_2.fastq",
                             outfile="SRR034098_3_trimq.fastq",
                             qualDiscard=20,minSeqLen=15)

####################
####################
# END OF PRE-PROCESSING
####################
####################
# 5. Align data to reference genome/transcriptome
####################
biocLite("Rsubread")
library(Rsubread)
# transcript reference file
# reference file downloaded from http://www.maizegdb.org/assembly
ref_file_trans<-"Zea_mays.AGPv3.22.dna.genome.fa"
# index files are written to the current working directory
# buildindex() is from Rsubread pkg
buildindex(basename="ref_index",reference=ref_file_trans)
# BAM file is written to working directory
# alignment align() is from Rsubread pkg
align(index="ref_index",readfile1="SRR034098_3_trimq.fastq",output_file="align_3_ResTrans.BAM")

####################
# 6. Annotate transcripts and count mapped read for transcript abundance
####################
# annotation reference downloaded from http://www.maizegdb.org/assembly
ann_ref<-"Zea_mays.AGPv3.22.gff3"
fc_3_corn<-featureCounts("align_3_ResTrans.BAM",annot.ext=ann_ref,
                         isGTFAnnotationFile=TRUE,reportReads=TRUE)

####################
####################
# END FILE
####################
####################
