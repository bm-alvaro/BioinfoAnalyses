#!/usr/bin/env bash

# Alvaro Benitez Mateo, 2022

SCRIPT_NAME=$(basename $0)
RED=$(tput setaf 1)
NORMAL=$(tput sgr0)
PICARD=/mnt/home/soft/picard/programs/x86_64/2.7.1

echo "########################################################################################"
echo "############################# Mapping and variant calling ##############################"
echo "########################################################################################"
echo "                              Alvaro Benitez Mateo, 2022"
echo "
"

Help()
{
    echo "Syntax $SCRIPT_NAME [-h help] [-i input_directory] [-o output_directory] [-t threads] [-p ploidy]"
    echo ""
    echo "This script is going to perform the first step out of two for variant calling on input data"
    echo ""
    echo "Options: "
    echo "-i    Path to data directory"
    echo "-o    Path to output directory"
    echo "-t    Number of threads to use"
    echo "-p    Sample ploidy. In case of analysing a pooled sample, use this parameter as an indicator of the number of pooled invididuals (p = #ind * ploidy)"
    echo ""
}

while getopts ":i:o:t:p: :h" option;

    do

    case $option in 
    h)
        Help
        exit;;
    esac

    case $option in 

    i) INPUT=$OPTARG
    ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
    Help >&2
    exit 1
    ;;

    o) OUTPUT=$OPTARG
    ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
    Help >&2
    exit 1
    ;;

    t) THREADS=$OPTARG
    ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
    Help >&2
    exit 1
    ;;

    p) PLOIDY=$OPTARG
    ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
    Help >&2
    exit 1
    ;;

    \?) 
    printf "${RED}Error: Invalid option -%s\n${NORMAL}" "$OPTARG" >&2
    Help >&2
    exit 1
    ;;

    esac

done

echo "The pathing to be use is:
$INPUT as the directory containing the raw data
$OUTPUT as the directory to stream the data output

$THREADS threads to use
$PLOIDY : sample ploidy
"

# Check input directory

if [ ! -d $INPUT ]; then
    printf "${RED}Error: Input directory does not exist\n${NORMAL}"
    exit 0
fi

# Check input files 
fq_files=$(ls $INPUT/*.fq.gz | wc -l);

if [ $fq_files != 0 ]; then
    printf "Raw data is available in -%s\nNext step is to create the output subdirectories" "$INPUT"

else
    printf "${RED}Error: Input directory does not contains fq files\n${NORMAL}"
    exit 0
fi

# Check resources directory and files

RESOURCES=/mnt/home/users/bio_293_ual/alvarobm/Alvaro/resources/RESOURCES

if [ ! -d $RESOURCES ]; then
    printf "${RED}Error: Resources directory does not exist\n${NORMAL}"
    exit 0
fi

resources_files=$(ls $RESOURCES/*.fa | wc -l);

if [ $resources_files != 0 ]; then
    printf "Fasta file available in -%s\n" "$RESOURCES"

else 
    printf "${RED}Error: Input directory does not contains reference genome files\n${NORMAL}"
    exit 0
fi

check_index_ref_bwa=$(ls $RESOURCES/cpepo_genome_v4.1.ann | wc -l);

if [ $check_index_ref_bwa != 0 ]; then
    printf "BWA index available in -%s\n" "$RESOURCES"

else
    printf "${RED}Error: Input directory does not contains indexed reference genome files by BWA, trying to index now\n${NORMAL}"
    bwa index -p $RESOURCES/cpepo_genome_v4.1 $RESOURCES/cpepo_genome_v4.1.fa
fi

check_index_ref_samtools=$(ls $RESOURCES/cpepo_genome_v4.1.fa.fai | wc -l);

if [ $check_index_ref_samtools != 0 ]; then
    printf "SAMTOOLS index available in -%s\n" "$RESOURCES"

else
    printf "${RED}Error: Input directory does not contains indexed reference genome files by SAMTOOLS, trying to index now${NORMAL}"
    samtools faidx $RESOURCES/cpepo_genome_v4.1.fa > $RESOURCES/cpepo_genome_v4.1.fa.fai
fi

check_index_ref_dict=$(ls $RESOURCES/cpepo_genome_v4.1.dict | wc -l);

if [ $check_index_ref_dict != 0 ]; then
    printf "Dict index available in -%s\n" "$RESOURCES"

else 
    printf "${RED}Error: Input directory does not contains Dict file, trying to index now\n${NORMAL}"
    java -jar $PICARD/picard.jar CreateSequenceDictionary \
      R=$RESOURCES/cpepo_genome_v4.1.fa \
      O=$RESOURCES/cpepo_genome_v4.1.dict
fi

# Create output subdirectories if they don't exist yet

mkdir -p $OUTPUT/fastqc/logs
mkdir -p $OUTPUT/samtools/logs
mkdir -p $OUTPUT/gatk-individual/logs
mkdir -p $OUTPUT/gatk-joint/logs

# Defining directory variables

OUT_FASTQC=$OUTPUT/fastqc
LOGS_FASTQC=$OUT_FASTQC/logs
OUT_SAMTOOLS=$OUTPUT/samtools
LOGS_SAMTOOLS=$OUT_SAMTOOLS/logs
OUT_GATK_IND=$OUTPUT/gatk-individual
LOGS_GATK_IND=$OUT_GATK_IND/logs
OUT_GATK_JOINT=$OUTPUT/gatk-joint
LOGS_GATK_JOINT=$OUT_GATK_JOINT/logs

echo "Alignment with bwa"

echo "BWA (mapping step) ---- STARTING"
date
echo ""

for R1 in $INPUT/*1.fq.gz

    do

        FILE_NAME=$(basename $R1 | sed "s/_1.fq.gz//" - )
        FILE_NAME_R1=$(basename $R1 | sed "s/.fq.gz//" - )
        R2=$(echo $R1 | sed "s/1.fq.gz/2.fq.gz/" - )
        FILE_NAME_R2=$(basename $R2 | sed "s/.fq.gz//" - )

        echo "##################################################"

        date
        echo "Left pair: $FILE_NAME_R1"
        echo "Right pair: $FILE_NAME_R2"

        echo "##################################################"

        SAMPLE=$FILE_NAME
        PLATFORM="ILLUMINA"
        LIBRARY="Lib-1"
        TPU="none"
        RGID=${SAMPLE}_${LIBRARY}
        RG="@RG\tID:${RGID}\tSM:${SAMPLE}\tPL:${PLATFORM}\tLB:${LIBRARY}\tPU:${TPU}"

        time bwa mem \
            -R $RG \
            $RESOURCES/cpepo_genome_v4.1 \
            $R1 $R2 \
            -t $THREADS | samtools sort -@ $THREADS -m 2G -O bam -T $OUT_SAMTOOLS/$FILE_NAME.tmp -o $OUT_SAMTOOLS/$FILE_NAME.sorted.bam - 2>> $LOGS_SAMTOOLS/bwa_mem_samtools_sort.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "BWA (mapping step) ---- FINISHED"
date
echo ""

echo "PICARD (identifying duplicates step) ---- STARTING"
date
echo ""

for BAM in $OUT_SAMTOOLS/*.sorted.bam

    do
        
        FILE_NAME=$(basename $BAM | sed "s/.sorted.bam//" - )

        echo "##################################################"

        date
        echo "Sample: $FILE_NAME"

        echo "##################################################"

        time java -jar $PICARD/picard.jar MarkDuplicates \
            I=$BAM \
            O=$OUT_SAMTOOLS/$FILE_NAME.marked_duplicates.bam \
            M=$OUT_SAMTOOLS/$FILE_NAME.marked_dup_metrics.txt 2>> $LOGS_SAMTOOLS/MarkDuplicates.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "PICARD (identifying duplicates step) ---- FINISHED"
date
echo ""

echo "SAMtools (indexing bam-files step) ---- STARTING"
date
echo ""

for BAM in $OUT_SAMTOOLS/*.marked_duplicates.bam

    do

        FILE_NAME=$(basename $BAM | sed "s/.marked_duplicates.bam//" - )

        echo "##################################################"

        date
        echo "Sample: $FILE_NAME"

        echo "##################################################"

        time samtools index -b $BAM 2>>$LOGS_SAMTOOLS/samtools_index.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "SAMtools (indexing bam-files step) ---- FINISHED"
date
echo ""

echo "GATK (HaplotypeCaller step) ---- STARTING"
date
echo ""

module unload picard
module load gatk

for BAM in $OUT_SAMTOOLS/*.marked_duplicates.bam

    do

    FILE_NAME=$(basename $BAM | sed "s/.marked_duplicates.bam//" - )

    echo "##################################################"

    date
    echo "Sample: $FILE_NAME"

    echo "##################################################"

    gatk --java-options "-Xmx144g" HaplotypeCaller \
        -R $RESOURCES/cpepo_genome_v4.1.fa \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -I $BAM \
        -O $OUT_GATK_IND/$FILE_NAME.g.vcf.gz \
        --sample-ploidy $PLOIDY

    echo "Processing of the file $FILE_NAME is finished :-) !!!"
    echo "##################################################"

done

echo "GATK (HaplotypeCaller step) ---- FINISHED"
date
echo ""