## How to run the MESA work directory

The MESA work directory provided in this repository is run in an unusual manner. The output directory, initial metallicity, initial mass, core overshooting and envelope mixing are not set in the inlist, but rather from the terminal. To set up and run your own model, take the following steps:

1. Copy the rate tables in this repository to a location of your choosing, preferably easy to reach for your MESA executable. 
2. Open *inlist_project_Vanrespaille2025* and set the five `*_cache_dir` paths.
3. Ensure the nuclear network file *Brinkman_pp_cno_alphabackbone.net* is saved in an accessible location
4. Ensure the output paths exist. By default, there are given by `{your_output_dir}/profiles/` and `{your_output_dir}/gyre/input_files`. 
5. Decide whether you want to run a pre-main-sequence model and whether you wish to run the model to the ZAMS or to the TAMS. If you wish to...
  - Run from ZAMS to TAMS: Ensure that the file `preMS/Z{your_Z}_M{your_M}_ZAMS.mod` exists. Open *inlist_project_Vanrespaille2025* and set `create_pre_main_sequence_model = .false.` and `stop_near_zams = .false.`
  - Run from birth to ZAMS: Open *inlist_project_Vanrespaille2025* and set `create_pre_main_sequence_model = .true.` and `stop_near_zams = .true.`. This will create the file `preMS/Z{your_Z}_M{your_M}_ZAMS.mod`
  - Run from birth to TAMS: Open *inlist_project_Vanrespaille2025* and set `create_pre_main_sequence_model = .true.` and `stop_near_zams = .true.`.
  
  Note that this directory comes with a number of ZAMS models included in the preMS subdirectory.

6. Do `./clean` and `./mk`
7. Now run the model using: `./rn <<< $(echo \"{your_output_dir}\", {your_Z}, {your_M}, {your_logD}, {your_aov}, {your_fov})`. 
  
Herein the variables you must set are
 -  `{your_output_dir}` the path to where the output should be saved. Using a full path is recommended. Ensure this path exists. 
 - `{your_Z}` the initial metallicity mass fraction, up to three decimal points. 
 - `{your_M}` the initial mass in solar mass, up to two decimal points. 
 - `{your_logD}` the logarithmic mixing strength at the base of the envelope in cm^2/s, up to one decimal point. 
 - `{your_aov}` the step core overshoot parameter in local scale heights, up to three decimal points. 
 - `{your_fov}` the exponential core overshoot parameter in local scale heights, up to three decimal points. 

`{your_aov}` and `{your_fov}` should not both be non-zero. At least one must always be zero. If `{your_aov}` is zero and `{your_fov}` is non-zero, the model will use exponential core overshoot. If `{your_aov}` is non-zero and `{your_fov}` is zero, the model will use step core overshoot. If both `{your_aov}` and `{your_fov}` are zero, no core boundary mixing is performed. The overshooting around shell convection zones is always done by exponential overshoot with `overshoot_f=0.010`.
