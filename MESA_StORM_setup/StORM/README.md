This directory contains the scripts and executable to compute the StORM grid after the grid of MESA stellar models has been completed. One slurm job is launched for each evolutionary track in the grid. Each of these jobs computes the StORM output at all available ages and at all requested rotation rates relative to the critical rotation rate. Such a job broadly consists of four steps:

1. Create a storm work directory under the evolutionary track in the grid of stellar models and copy the relevant files over 
2. Convert the `.GYRE` input models to `.GSM` models and prepare the commands inlists for StORM
3. Perform all the StORM computations
4. Compute all the storm output
5. Remove files that are no longer needed

----------------------------

Herein, it is assumed that the stellar models are organised as follows:
- one top level directory
    - a subdirectory for each evolutionary model named as *Z{Z}_M{M}_logD{logD}_fov{fov}*, containing at least
        - history file named as *Z{Z}_M{M}_logD{logD}_fov{fov}.hist*
        - gyre
            - input_models : directory containing a `.GYRE` input model created by MESA for each output step, named as *Z{Z}_M{M}_logD{logD}_fov{fov}_Xc{Xc}_mn{mn}.GYRE*
        - storm\*
            - input_models\* : directory containing an input model for StORM, either created by MESA or with our own conversion script, for each output step, named as *Z{Z}_M{M}_logD{logD}_fov{fov}_Xc{Xc}_mn{mn}.HDF5*
            - output\* : directory containing the StORM output files, named as *storm_Xc{Xc}_mn{mn}_frot{frot}.HDF*

\* Directory and the files therein will be created with the scripts in this directory 

----------------------------

This directory contains the following files: 
 - *commands_template.txt* : A template of the StORM commands input, to be completed with the input model, output file and rotation frequency. 
 - *create_commands_per_model.py* : A script meant to create a commands input file for StORM from the commands template for each output step in the evolutionary track
 - *frots.txt* : A list of the rotation rates relative to the Keplerian critical rotation rate at which to compute StORM output
 - *run_STORM.sh* : Shell script that executes the five steps above on a provided set of model parameters using the scripts at the given paths. 
 - *STORM_parameters.csv* : An overview of the grid parameters at which the. One slurm job is created for each line in this file. 
 - *storm* : The (static) storm executable
 - *submit_all.sh* : File to execute to read all parameters from parameter file set within and submits all jobs
 - *submit_STORM.slurm* : Describes the slurm settings used for each job

 The create_commands_per_model.py script depends on the following packages and was used with the following versions. It is expected to not be strongly dependent on the particular versions of these packages however. 
  - `python (3.8.8)`
  - `numpy (1.23.5)`
  - `h5py (2.10.0)`
