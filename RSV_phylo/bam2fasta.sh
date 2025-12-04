#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=consensus2
#SBATCH --time=01:00:00
#SBATCH --partition=prod
#SBATCH --account=kunf0091
#SBATCH --output=consensus2.%A_%a.out
#SBATCH --error=consensus2.%A_%a.err

#Modify the array size accordingly, you can use the command "wc -l file_list.txt"
#SBATCH --array=1-263
 
#Load the modules
module load miniconda
conda activate nanocaller_env

# run like:
# sbatch reads2vcf.sh <file_list> <ref_name>
# e.g. sbatch bam2fasta.sh filelist.txt EPI_ISL_1653999.fasta
#Get imput and output file from the list
f=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

#IFS='.' read -r -a f <<< "$i";
outdir=/dpc/kunf0055/RSVB/
wdir=/dpc/kunf0055/RSVB/sam2/

# Run the command
#samtools index $f.s.f.bam

/dpc/kunf0055/RSVB/NanoCaller/NanoCaller --bam ${wdir}/$f.s.f.bam --ref $2 --output ${outdir}/vcf/$f --prefix $f 

# 2. Generate the consensus FASTA
bcftools consensus -f $2 --missing N ${outdir}/vcf/$f/$f.vcf.gz | sed "1s/.*/>$f/" > ${outdir}/consensus2/$f.fasta

#bgzip ${outdir}/nanocaller.vcf
#bcftools index ${outdir}/nanocaller.vcf.gz

# 2. Generate the consensus FASTA
#bcftools consensus -f reference.fasta ${outdir}/nanocaller.vcf.gz > ${outdir}/consensus.fasta
