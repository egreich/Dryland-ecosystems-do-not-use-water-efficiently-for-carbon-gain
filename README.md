# DEPART-SAM Model

This repository contains models and code to run the DEPART-SAM model. An explanation of the file structure is below:

1. scripts
  - 00_data_clean.R
    - This script reads in all necessary data and formats the information consistently across sites
    - This includes filling small data gaps
  - 01_model_prep.R
    - preps the data to be in the correct format for the Bayesian model
  - 02_script_SAM.R
    - runs the DEPART-SAM model
  - SAM_function.R
    - read in by 02_script_SAM.R
  - check_convergence_function.R
    - read in by 02_script_SAM.R
  - coda_functions.R
    - functions to help summarize JAGS coda output
  - functions.R
    - helper functions, read in by various scripts

2. models
  - Model_v3.R
    - combined DEPART-SAM model
  - Model_split_v3.R
    - combined DEPART-SAM model, but splits the model to account for a large data gap at US-Vcp

3. input_data
  - read into 00_data_clean.R

4. clean_data
  - read into 01_model_prep.R and 02_script_SAM.R

5. output_dfs
  - output data frames. Contains posterior means and uncertainty for variables of interest.
