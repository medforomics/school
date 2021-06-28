#!/bin/bash
#SBATCH --share
#SBATCH --partition=short
#SBATCH --job-name=goal_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=8:00:00
#SBATCH --mem=65536


help_and_exit() {
  local retval=${1:-1}
  cat <<EOF
${0##*/} [-s|--sample] [-p|--panel] [-c|--config] fastq1 fastq2
EOF
  exit "$retval"
}

configFile="snakemake_goalConfig.json"

while (( $# )); do
  case $1 in
    -s|--sample)  sample=$2; shift ;;
    -p|--panel) panelFolder=$2; shift ;;
    -c|--config) configFile=$2; shift ;;
    -*)        printf 'Unknown option: %q\n\n' "$1"
               help_and_exit 1 ;;
    *)         args+=( "$1" ) ;;
  esac
  shift
done
set -- "${args[@]}"

fq1=$1
fq2=$2
fq1_file="$(basename -- $fq1)"
fq2_file="$(basename -- $fq2)"

# Path to reference files
refs="/data/user/$USER/GOAL/genome/GOAL_GRCh38"

# Directory to write results
working_directory="/data/user/$USER/GOAL/OUTPUT/GRCh38/$sample"


if [ ! -d "$working_directory" ]; then
    mkdir -p $working_directory
fi
if [ ! -f "$working_directory/$fq1_file" ]; then
    cp $fq1 $working_directory
fi
if [ ! -f "$working_directory/$fq2_file" ]; then
    cp $fq2 $working_directory
fi

module load Singularity
module load snakemake

# if not using slurm, then SLURM_CPUS_ON_NODE value must be set manually
# to the number of cores you want the job to use
export SINGULARITYENV_SLURM_CPUS_ON_NODE=$SLURM_CPUS_ON_NODE
snakemake -r --cores all --use-singularity \
    --singularity-prefix /data/user/$USER/GOAL/.snakemake \
    --singularity-args "-B $panelFolder:/panel,$refs:/refs --cleanenv -c " \
    --directory $working_directory \
    --latency-wait 30 \
    --configfile $configFile \
    --config sample=$sample fq1=$fq1_file fq2=$fq2_file
