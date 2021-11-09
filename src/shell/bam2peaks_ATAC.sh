#!/bin/bash
#$ -S /bin/bash
#$ -cwd

set -e

#======================================================================================
# bam2peak.sh
# 
#   Description
#       This shell script crocess the bam file nad perform MACS2 based peak-calling
#       out put files are stored in the "peakcalling" directory
# 
#   Usage
#       bash bam2peaks.sh -f [input.bam] -r [s or p] -o [out_prefix] -p [cut off p-value]
#       
#   Argument
#       -h                             Display help
#       -f              bam file       Input bam file
#       -r              charactor      read type (s (single-end) or p (paired-end))
#       -o              charactor      prefix for output files
#       -p              numeric        cut off p-value for  macs2 peak-calling
# 
#   output files
#       ${prefix}.noMT.bam                  bam file            mitochondorial genome removed
#       ${prefix}.proecssed.bam             bam file            mitochondorial genome and PCR dupicate removed (if the read type is paired-end) sorted
#       ${prefix}.narrowPeak                bed file            narrow peaks
#       ${prefix}.filt.narrowPeak           bed file            black list removed narrow peaks

#======================================================================================

###------------------------------------------------------------------------------------
### healp
###------------------------------------------------------------------------------------
function usage {
cat <<EOM
Usage: $(basename "$0") -f [input.bam] -r [s or p] -o [out_prefix] -p [cut off p-value].
    -h          Display help
    -f VALUE    file path
    -r VARUE    read type (s (single-end) or p (paired-end))
    -o VALUE    output prefix
    -p VALUE    cut-off p-value for peak calling
EOM

exit 2
}
###------------------------------------------------------------------------------------
### Definition of process for each argument
###------------------------------------------------------------------------------------
while getopts ":f:r:o:p:h" optKey; do
    case "$optKey" in
    f)
        file=${OPTARG}
        ;;
    r)
        read=${OPTARG}
        ;;
    o)
        prefix=${OPTARG}
        ;;
    p)
        pval=${OPTARG}
        ;;
    '-h'|'--help'|* )
        usage
        ;;
    esac
done

echo "INFO :`date`" 1>&2
echo "COMMAND: bam2peak.sh -f $file -r $read -o $prefix -p $pval" 1>&2
echo -e "\n===== `date`: peakcalling directory generate =====" 1>&2
mkdir -p peakcalling

if [ $read = "p" ] ;then
    ###------------------------------------------------------------------------------------
    ### For Paired-end bam processing
    
            # using Harvard ATAC-seq module (python script removeChrom.py)
            # https://github.com/harvardinformatics/ATAC-seq
    
    echo -e "\n===== `date`: paired-end bam file processing =====" 1>&2
    ###------------------------------------------------------------------------------------
    samtools view -h $file | removeChrom.py - - chrM | samtools view -h \
        |samtools collate -@ 4 -O - \
        |samtools fixmate -@ 4 -m - - \
        |samtools sort -@ 4 -O BAM \
        |samtools markdup -@ 4 -r -s - - \
        |samtools sort -@ 4 -o ./peakcalling/${prefix}.processed.bam
        samtools index ./peakcalling/${prefix}.processed.bam
    
elif [ $read = "s" ] ;then
    ###------------------------------------------------------------------------------------
    ### For single-end bam processing
    
            # using Harvard ATAC-seq module (python script removeChrom.py)
            # https://github.com/harvardinformatics/ATAC-seq
    
    echo -e "\n===== `date`: single-end bam file processing ... =====" 1>&2
    ###------------------------------------------------------------------------------------
    samtools view -h $file | removeChrom.py - - chrM | samtools view -h \
        |samtools collate -@ 4 -O - \
        |samtools sort -@ 4 -o ./peakcalling/${prefix}.processed.bam
        samtools index ./peakcalling/${prefix}.processed.bam
fi

###------------------------------------------------------------------------------------
### peak calling

echo -e "\n===== `date`: Peak calling =====" 1>&2
###------------------------------------------------------------------------------------
macs2 callpeak -t ./peakcalling/${prefix}.processed.bam \
    -f BAM \
    -g "hs" \
    -p $pval \
    --shift -37 \
    --extsize 73 \
    -B --SPMR \
    -n ./peakcalling/${prefix}

###------------------------------------------------------------------------------------
### # Filter blacklist regions and unplaced contigs

echo -e "\n===== `date`: Remove blacklist regions =====" 1>&2
###------------------------------------------------------------------------------------

bedtools intersect -v -a ./peakcalling/${prefix}_peaks.narrowPeak -b ./data/hg19_blacklist.bed \
    | grep -P 'chr[\dXY]+[ \t]' > ./peakcalling/${prefix}_peaks.filt.narrowPeak

echo -e "\n===== `date`: Done! =====" 1>&2
