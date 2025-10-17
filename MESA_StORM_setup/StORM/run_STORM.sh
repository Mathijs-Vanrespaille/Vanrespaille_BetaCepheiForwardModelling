#!/bin/bash
export Zini=$1 # Initial metallicity
export Mini=$2 # Initial mass
export logD=$3 # Envelope mixing strength
export aov=$4 # Core step overshoot parameter, implicitely assumed to always be zero below as it is not included in the naming scheme
export fov=$5 # Exponential core overshoot parameter
export frot_file=$6 # Full path to the file listing all relative rotation rates at which to compute StORM output
export MESA_output_dir=$7 # The full path to the top level directory of the MESA model grid
export STORM_setup_dir=$8 # The full path to where the Python scripts are stored. Probably the same directory this file is in
export STORM_binary=$9 # The full path to the (static) StORM binary

# Echo the input parameters for testing and debugging purposes
#echo Zini $Zini
#echo Mini $Mini
#echo logD $logD
#echo aov $aov
#echo fov $fov
#echo frot_file $frot_file
#echo MESA_output_dir $MESA_output_dir
#echo STORM_setup_dir $STORM_setup_dir
#echo STORM_binary $STORM_binary

# Go into the MESA storm directory
export STORM_work_dir=${MESA_output_dir}/Z${Zini}_M${Mini}_logD${logD}_fov${fov}/storm
mkdir -p $STORM_work_dir
cd $STORM_work_dir

# Prepare where to store the output and error messages
mkdir -p job_logs/ job_errs/

# Copy over the necessary files
cp $STORM_binary .
cp ${STORM_setup_dir}/commands_template.txt .
cp ${STORM_setup_dir}/create_commands_per_model.py .
cp $frot_file .
mkdir commands

# Make the inlists
python3 create_commands_per_model.py $frot_file ${MESA_output_dir}/Z${Zini}_M${Mini}_logD${logD}_fov${fov} # No / at the end of this directory!!!

# Execution:
for command_file in commands/*.txt
do
# Set up files to save the terminal output
  command_filename=$(basename ${command_file})
  command_filename=${command_filename::-4}
  log_file=job_logs/$command_filename.log
  err_file=job_errs/$command_filename.err
# Now run StORM, saving the output
  ./storm < $command_file 1>>$log_file 2>>$err_file
done

# Removing the executable as it is quite large
rm storm
# Removing the various support files and scripts
rm $frot_file create_commands_per_model.py commands_template.txt
rm -r commands
