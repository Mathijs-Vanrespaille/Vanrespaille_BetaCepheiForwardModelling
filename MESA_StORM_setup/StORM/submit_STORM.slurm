#!/bin/bash
#SBATCH --job-name=run_STORM_test
#SBATCH --account=ivsusers
#SBATCH --time=1440
#SBATCH --output=std/stdout_test.log
#SBATCH --error=std/stderr_test.log
#SBATCH --cpus-per-task=1
#SBATCH --hint=multithread
#SBATCH --ntasks=1
#SBATCH --mem=1400
#SBATCH --partition=normal
#SBATCH --qos=normal

export Zini=$1
export Mini=$2
export logD=$3
export aov=$4
export fov=$5
export frot_file=$6
export MESA_output_dir=$7
export STORM_setup_dir=$8
export STORM_binary=$9

source ~/.bashrc

# Run StORM on all output steps of the MESA model with the given parameters using a given rotation rates.
bash run_STORM.sh $Zini $Mini $logD $aov $fov $frot_file $MESA_output_dir $STORM_setup_dir $STORM_binary
