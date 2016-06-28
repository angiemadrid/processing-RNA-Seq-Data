# processing-RNA-Seq-Data
Project for DSCS6020 (NEU): Bioinformatics pipeline using Bioconductor for processing RNA-Seq data 

INTRODUCTION

This project followed a typical bioinformatics pipeline starting at the
beginning by reading and preprocessing raw sequence data, then continuing
with aligning the data with a reference transcriptome, and performing
annotation and simple analysis. This report will discuss data collection, a
description of the data file formats, and a description of the R packages used
in the implementation of this pipeline. It will conclude with details of each
step in the pipeline and a summary of outcomes.
Transcriptome analysis is the process for finding and characterizing genes,
and measuring gene expression (i.e., mRNA composition). One way of
conducting such an analysis is by high-throughput RNA sequencing
experiments (RNA-Seq). The general process for a RNA-Seq experiment is
discussed in the next section. In general, the raw sequencing data are the
output of short sequence reads of amino acids that are contained within a
sample of RNA. The original sample of RNA is fragmented and fed into a
sequencing instrument designed to scan and detect amino acids located
within the fragments of RNA[1]. The size of the output data ranges from
100MB to 10GB and greater. The online database GEO (Gene Expression
Omnibus), found at http://www.ncbi.nlm.nih.gov/geo/, is a public repository of
raw datasets supplied by scientists and researchers performing studies using
sequencing instruments to analyze RNA samples. There are currently
1,552,201 samples in the GEO repository. A raw, unprocessed dataset was
downloaded from GEO for this project (as discussed below).
The tools used to create a pipeline are all found in R packages provided by
Bioconductor (https://www.bioconductor.org/), which is a software tool used by many computational
biologists and bioinformaticians working with biological data. 

WHERE THE DATA COMES FROM

~What is the data~

The type of data used in this project is created when a biological sample of
RNA is processed through a next generation sequencing (NGS) instrument.
RNA samples are made up of sequences composed of four nucleobases: 
adenosine (A), guanosine (G), uracil (U), and cytosine (C). Before the sample
is introduced into the sequencer, the RNA is converted to cDNA
(complementary DNA) through a process called reverse transcription. To
synthesize cDNA from RNA, the RNA sequence is fragmented and stripped of
non-coding RNA. The fragmented sequences are then reversed and all the
uracil molecules are converted into thymine (T). These newly formed,
shortened sequences of DNA build a sequencing library that will be scanned
and read by the sequencing instrument. The actual sequencing procedure
carried out in NGS instruments varies according to manufacturer. However,
in general, once the library sequences are constructed, adapter sequences are
ligated to the ends of the library sequences, which help with the amplification
of the library and to start the sequencing (bio-chemical) reaction. Each
molecule is then fluorescently tagged, and the sequences are imaged as they
pass through the sequencer detector. Within the software of the instrument,
the images are analyzed, bases are detected, and base calls are assigned. The
resulting output is millions of short-sequence reads, which are 30-400 base
pairs (bp) long and correspond to individual cDNA fragments.[2,3]

~How data is collected~

Due to their small size, short reads often do not span the entire cDNA
sequence resulting in uneven coverage of individual transcripts. In addition,
the millions of reads produced from a sequenced sample only provide a small
percent of any given base existing in the library resulting in a sampling
variance. There is a term, depth of coverage, which refers to the percent of
transcripts measured, and in order to increase coverage and reduce sampling
variation, deeper sequencing is required. Wang et al. (2011)[4] suggest that a
minimum of 10 million reads per sample is needed in order to achieve above median
coverage.

~What are the data file formats~

There are a variety of file formats for unprocessed (raw) sequence data used,
depending on which instrument is used for the sequencing. In general, these
files are ASCII text files and contain sequence data (a string of G’s, C’s, A’s
and T’s) and other metadata about the sequence. Researchers conducting
RNA-Seq experiments provide their unprocessed data to NCBI, where they
are housed in the Sequence Read Archive (SRA). SRA files may be accessed
through GEO and are converted and saved in the FASTQ file format, which
specifically, contains the read identifiers (assigned by the sequencer
software), nucleotide sequences (i.e., the reads) and quality scores.

~Phred quality scores~

The Phred quality scores provided in FASTQ files are a measure of the
probability that a base call of a given nucleotide is not correct. For example, a
low Phred score of 10 for a given base means that the probability that this
base call is incorrect is 1 in 10, or 90% accuracy. A typical threshold is a
Phred of 20, meaning the probability of an incorrect base call is 1 in 100, or
99% accuracy. If the sequencer software cannot determine a Phred score, no
score is given and an “N” is assigned. Quality scores are assigned by the
sequencer software using Phred algorithms, which calculate scores according
to signal intensity measurements during sequencing. The encoding of Phred
scores in FASTQ files are the ASCII values of the characters 0 to 93 with an
offset of 33, because ASCII values of 0-32 encode for characters that are nonprintable.
The lowest quality score is indicated by “!”, and the highest by “~”.[5]

BIOCONDUCTOR

Bioconductor (https://www.bioconductor.org/) is an open-source R project for
handling and analyzing high-throughput sequencing data and was used for
the implementation portion of this project. The specific Bioconductor
packages used to complete each step in this project’s pipeline are discussed
below. The majority of Bioconductor’s packages, and all the packages used
here create S4 objects for storage and retrieval of the data, instead of the
more informal S3 class object, a data frame.

~Packages~

SRAdb: This package establishes an SQL database using the data from
GEO’s website and interacts directly with the website for accessing the
datasets. Once the database is loaded in R, the SRA files are accessed and
downloaded using the getSRAfile() function and specifying the desired run
number as an argument, “SRR034098” for this project.

SRA toolkit: NCBI provides a toolkit for converting SRA files to FASTQ files.
The toolkit must be installed on the local computer and is accessed using the
fastq-dump() function through a system call in R, e.g.,:
system("fastq-dump --split-3 SRR034098.sra")
The newly converted FASTQ file is saved in the local working directory.

qrqc: This package handles quality control of sequence data and provides set
of analyses for a summary of the overall quality of the data. A
FASTQSummary-class object is created when a FASTQ file is read using the
readSeqFile() function.

ShortRead: The function readFastq(), from the ShortRead package, reads a
FASTQ file and produces as ShortReadQ object with slots for sequences, IDs,
and quality scores. This package also enables writeFastq(), which was used in
this report for writing output during the preprocessing stage.

giraffe: This package trims adapters through its trimAdapter() method, which
takes a ShortReadQ object as an argument.

seqTools: This package provides a method, trimFastq(), that will trim
sequences in a FASTQ file according to quality score and minimum sequence
length, and write to an output file.

Rsubread: This package performs alignment using the seed-and-vote
algorithm. This algorithm is based on hash tables. Tophat, a widely used
aligner, also uses a hash-table-based algorithm called seed-and-extend. Liao
et al., have shown that seed-and-vote algorithm is more efficient and accurate
than the seed-and-extend method. This package achieves alignment by first
building an index of the reference genome using the buildindex() function,
which takes as an argument a FASTA file containing the reference sequence.
The align() method maps the preprocessed sequences in the FASTQ file to the
index of the reference and creates BAM files. A BAM file is a binary version
of the SAM (Sequence Alignment/Map) file format, a format for storing
sequence alignments. Since BAM files are compressed versions of SAM files,
they are more efficiently used in software pipelines.
This package also provides the function featureCounts(), which handles
annotation and counting mapped reads.

DESCRIPTION OF PIPELINE

~Obtain data~

The dataset used in this project comes from a sequencing study (Nobuta,
2008)[6] that examined RNA features of size distribution, tissue-specific
regulation and sequence conservation in leaves, ears and tassels of maize
(corn). Specifically, this project used the sequencing output from corn tassels.
The size of the SRA file is 85.2Mb and was downloaded to the local hard drive
using Bioconductor’s SRAdb package, which accessed GEO’s online SRA
database. SRA files contain the sequencing output directly from the
sequencing instrument and are in a raw format. 

The data in the SRA need to be converted into FASTQ format before any preprocessing
or analyzing of the data can occur. This is done using the fastqdump
method from NCBI’s SRA toolkit[7].

Each read in the FASTQ file contains a line with a “@” character and a
sequence identifier. The next line is the sequence of bases (G’s, C’s, A’s and
T’s). The third line is also the sequence identifier preceded by a “+” character.
The last line is the sequence of quality scores associated with each base in the
second line.
A quality check of the raw data is needed to ensure that there are not any
problems or biases resulting from the sequencer itself or the sequencing
library.

Using the plotting functions for this object the quality frequencies by position
and per base sequence content can be visualized. The qrqc package also
provides plotting methods.

~Preprocess data~

The next step in the pipeline is the preprocessing of the data. First, the
adapters need to be trimmed from the sequences. Recall that adapter
sequences are ligated to the ends of the library sequences to help with
sequencing (bio-chemical) reaction. It is important to remove any portion of
the adapter sequence from the reads; otherwise the alignment will be
affected. According to the study page for Nobuta, 2008,
(http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM433622) there were
two adapter sequences that need to be trimmed. They are
“TCGTATGCCGTCTTCTGCTTG” for the 3’ adapter and
“GTTCAGAGTTCTACAGTCCGACGATC” for the 5’ adapter.
The FASTQ file containing the raw reads was read into R using the
readFastq() function from Bioconductor’s ShortRead package.

The first adapter was trimmed from the sequences using trimAdapter()
function from the giraffe package. For the second adapter sequence, the
trimLRPatterns() function was used to trim the sequences. Next, sequences
were trimmed by quality using the trimFastq() function from the seqTools
package and written to an output FASTQ file.
Sequences were then discarded according to a quality score threshold of 20
and a length shorter than 15. The Nobuta (2008) study discarded sequences
shorter than 15. Their preprocessing of the data resulted in 2,463,445
sequences remaining. In order to obtain a similar number of sequences, a
quality score of 20 was used for this project. The preprocessing for this project
resulted in 2,371,031 sequences.

~Align data~

The next step after preprocessing is the alignment of the reads to a reference
genome. Mapping short reads to reference in order to align them to genomic
sequence aids in finding transcription structure by creating a map of
transcription structure and expression level of each gene.[3] The reference
genome was downloaded manually from the Maize Genetics and Genomics
Database (http://www.maizegdb.org/assembly) and saved to the working
directory. This file was then passed in to the buildindex() function, which
created a reference index. The method align() aligned the preprocessed FASTQ file to this reference
index producing a BAM file that is ready to be annotated and analyzed.

~Annotate and count mapped reads~

After reads have been aligned to a reference, the next step is to determine the
identities of the genes within the sample. This is done through the
annotation process, which links mapped transcripts to a reference annotation
file that contains genomic information. Transcript abundance may then be
calculated by counting the reads that are linked. The abundance is indicative
of which genes get translated.

The Rsubread package from Bioconductor provides the featureCounts()
function that loads the reference annotation file and assigns the aligned
reads to the genomic features. The reference annotation file (in GTT format)
used in this project was also downloaded from Maize Genetics and Genomics
Database.

In the end there were 12,615 sequences successfully annotated. 

REFERENCES

1 Morgan, Martin (2013) Sequences, Genomes, and Genes in R / Bioconductor
https://www.ebi.ac.uk/training/sites/ebi.ac.uk.training/files/materials/2013/13
1021_HTS/genesandgenomes.pdf

2 https://www.ebi.ac.uk/training/online/course/ebi-next-generationsequencing-
practical-course/what-next-generation-dna-sequencing/illumina-

3 RNA-seqlopedia http://rnaseq.uoregon.edu/

4 Wang, Y, N Ghaffari, CD Johnson, UM Braga-Neto, H Wang, R Chen, H
Zhou. (2011). Evaluation of the coverage and depth of transcriptome by RNASeq
in chickens. BMC Bioinformatics 12 Suppl 10:S5

5 Phred - Quality Base Calling, http://www.phrap.com/phred/

6 Nobuta K et al., "Distinct size distribution of endogeneous siRNAs in maize:
Evidence from deep sequencing in the mop1-1 mutant.", Proc Natl Acad Sci U
S A, 2008 Sep 24;105(39):14958-63

7 SRA Toolkit
http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=so
ftware&s=software

