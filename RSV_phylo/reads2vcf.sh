#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=minimap
#SBATCH --time=01:00:00
#SBATCH --partition=prod
#SBATCH --account=kunf0091
#SBATCH --output=assembly.%A_%a.out
#SBATCH --error=assembly.%A_%a.err

#Modify the array size accordingly, you can use the command "wc -l file_list.txt"
#SBATCH --array=1-263
 
#Load the modules
module load regenie_env

# run like:
# sbatch reads2vcf.sh <file_list> <ref_name>
# e.g. sbatch reads2vcf.sh filelist.txt EPI_ISL_1653999.fasta
#Get imput and output file from the list
f=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

#IFS='.' read -r -a f <<< "$i";
outdir=/dpc/kunf0055/RSVB/
wdir=/dpc/kunf0091/RSV/

# Run the command
minimap2 -a $2 $wdir/trimmed_reads/$f.fastq > $f.sam
 
samtools view -bS $f.sam | samtools sort -o $f.s.bam
samtools view -h -F 4 -b $f.s.bam >  $f.s.f.bam
samtools index $f.s.f.bam

#NanoCaller --bam $f.s.f.bam --ref $2 --output ${outdir}/vcf/$f --prefix $f 
#bgzip ${outdir}/nanocaller.vcf
#bcftools index ${outdir}/nanocaller.vcf.gz

# 2. Generate the consensus FASTA
#bcftools consensus -f reference.fasta ${outdir}/nanocaller.vcf.gz > ${outdir}/consensus.fasta
