source ~/.bashrc

# Set paths to where the MESA output and the folder with the templates and python scripts to complete them
export STORM_setup_dir=/STER/mathijsv/updatingSiemensGrid/STORM_setup
export MESA_output_dir=/STER/mathijsv/updatingSiemensGrid/MESA_output_fullGrid_fixed

# Set path to the csv summarising on which MESA input parameter runs to perform the oscillation analysis
export STORM_parameters=${STORM_setup_dir}/STORM_parameters.csv
# Set path to the storm executable
export STORM_binary=${STORM_setup_dir}/storm
# Set path to the file containing the dynamical rotation rates
export frot_file=${STORM_setup_dir}/frots.txt

exec < $STORM_parameters
read header
#echo "header is : $header"
while IFS="," read -r Zini Mini logD aov fov
do
  #echo "Zini : $Zini"
  #echo "Mini : $Mini"
  #echo "logD : $logD"
  #echo "aov : $aov"
  #echo "fov : $fov"
  # Run StORM on all output steps of the MESA model with the given parameters using a given rotation rates.
  echo $Zini $Mini $logD $fov
  sbatch submit_STORM.slurm $Zini $Mini $logD $aov $fov $frot_file $MESA_output_dir $STORM_setup_dir $STORM_binary
done
