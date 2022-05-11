#!/bin/bash

## run this script on a qrsh session, not on the head node
## this will submit cluster qsub jobs and monitor them
## example command to submit:
## snakemake -j 9 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time}"
## You can also qsub snakecharmer.sh to get logs from this

#$ -cwd
#$ -V
#$ -j yes
#$ -o logs/
#$ -N "snake"
#$ -l data
#$ -l h_rt=10:00:00


mkdir -p logs

## activate conda environment
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda activate CLIP_mapping


echo "############################### START PIPELINE #############################"
echo $(date)

snakemake --unlock && ##necessary because a lot of killed snakes

snakemake --jobs 12 --cluster-config cluster_config.json --cluster "qsub -cwd -V -j yes -o {cluster.err} -m {cluster.m} -M {cluster.account} -pe smp {cluster.n} -l m_mem_free={cluster.m_mem_free} -l h_rt={cluster.time} -N {cluster.name} -l data" --snakefile Snakefile --latency-wait 600 --rerun-incomplete --configfile config_parpipe.yaml

exit
