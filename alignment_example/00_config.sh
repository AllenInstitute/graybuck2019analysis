#!/bin/bash
#11/08/2018   ver 0.1
# Batch 41 MiSeq samples

#Config file helps to set up all variables that might change when running a sequencing analysis. 

#WorkingDirectory src/$executionDate

mainDir="scatac_analysis"
executionDate="2018-11-08" #yyyy-mm-dd
baseFolder=//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis
cd $baseFolder
workingDir=src/$executionDate

# Make sure that the base path is correct before starting the execution
PROJECT_FOLDER="$baseFolder/data/$executionDate" 
if [ ! -d $PROJECT_FOLDER ]; then
	echo "WARNING: The output directory [${PROJECT_FOLDER}] doesn't exist."
    echo "OUTPUT will be saved at: [${PROJECT_FOLDER}]"
    mkdir -p $PROJECT_FOLDER
    mkdir -p $PROJECT_FOLDER/{raw/{bam,fastq},clean}
fi

echo $PROJECT_FOLDER
export PROJECT_FOLDER=$PROJECT_FOLDER

## CONFIGURATIONS

## ALL COMMON OUTPUT FILES

alignmentBAM=$PROJECT_FOLDER/clean/aligned.bam
trimmedBAM=$PROJECT_FOLDER/clean/trimmed.bam
mappedBAM=$PROJECT_FOLDER/clean/mapped.bam

## ALL COMMON DIRECTORIES

# Fastq Directory
FD=$PROJECT_FOLDER/raw/fastq
FASTQ_DIR=$PROJECT_FOLDER/raw/fastq

# BAM Directory
BD=$PROJECT_FOLDER/raw/bam
BAM_DIR=$PROJECT_FOLDER/raw/bam
export BAM_DIR=$BD

# additional BAM file suffix (use .rmd.srt if the files are .rmd.srt.bam, for example)
SUFFIX=.rmd.srt
SUFFIX_BIGWIGS=rmd.srt.bam

# HotSpot output directory (use absolute/complete directory location)
#HD=$PROJECT_FOLDER/clean/hotspot_gt$FS
#HD=$PROJECT_FOLDER/clean/hotspot_${DS}M
HD=$PROJECT_FOLDER/clean/hotspot

# Stats Directory
SD=$PROJECT_FOLDER/clean/stats
STATS_DIR=$PROJECT_FOLDER/clean/stats/
export CONFIG_STATS_DIRECTORY=$SD
export CONFIG_STATS_RES_DIR=$baseFolder/results/$executionDate

# Preseq Directory
PD=$PROJECT_FOLDER/clean/preseq
PRESEQ_DIR=$PROJECT_FOLDER/clean/preseq/
export CONFIG_PRESEQ_DIRECTORY=$PD
export CONFIG_PRESEQ_EXCEL=$PRESEQ_EXCEL
export CONFIG_PRESEQ_RES_DIR=$baseFolder/results/$executionDate

# Insert Sizes Directory
ID=$PROJECT_FOLDER/clean/inserts
INSERTS_DIR=$PROJECT_FOLDER/clean/inserts/
export CONFIG_INSERTS_DIRECTORY=$ID


## ALL COMMON APPLICATIONS


# HotSpot distr location
HSD=//allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/hotspot-distr


# Bowtie Location
BT=//allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/bowtie-1.1.0/bowtie

#Bowtie settings
BS="-p 8 -m 1 -S -X 2000 --chunkmbs 512 --un unaligned.fastq"
BOWTIE_SETTINGS="-p 8 -m 1 -S -X 2000 --chunkmbs 512 --un unaligned.fastq"

##Standardize tools folder, to avoid conficts in this area

# Trim Galore location
TG=//allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/trim_galore
TRIMGALORE_BIN=//allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/trim_galore

# Samtools location
ST=samtools
SAMTOOLS_BIN=samtools

# Preseq location
PS=//allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/preseq-0.1.0.Linux_x86_64/preseq
PRESEQ_BIN=//allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/preseq-0.1.0.Linux_x86_64/preseq

# bamToBed location (from bedtools)
BB=bamToBed
BAMTOBED_BIN=bamToBed

# Collect Insert Size Metrics location (from Picard Tools)
CI="java -jar //allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/picard.jar CollectInsertSizeMetrics"
PICARD_CISM_BIN="java -jar //allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/picard.jar CollectInsertSizeMetrics"

# Check genome and read length compatibility at
# http://www.uwencode.org/proj/hotspot/

#Maybe have a copy inside this folder with all refs?
# Genome Index Location
GL=/data/rnaseqanalysis/RNAseq/indexes/mm10/genome

# Genome
GENOME=mm10

# Genome Read Length
RL=50
