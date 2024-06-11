# SimulateSAS
 Final Code and Scripts for the simulate SAS project
 
## Dependencies
You will need the following R libraries for these scripts to run:
- tidyverse 
- data.table
- dplyr 
- abc
- ggExtra (for final visualization)


Additionally, you will need to make sure the:
"ML_functions.R" script is in your working directory. 

## To Run the Pipeline
To simulate a single iteration of SAS -

1. Be sure to run "run_likelihood_cmd.R" on your observed haplotype counts. The output of this file *MUST* be named "ml_output_observed". The script should do this for you, but make sure it's in your working directory for downstream steps. 

2. Run "simSAS_cmd.R". You will additionally need "UKB_mafs_filtered.txt" and 
"UKB_r2_filtered.txt" in your working directory. This script will only simulate SAS for one combination of s and F. I looped over this script 50k times in a shell script to create the simulations for this project so they could easily be batched in parallel to a remote comuting cluster. 

3. Cat together the outputs of simSAS_cmd.R to produce a table of s, F, and sum of squared error values for each simulation. I pulled this locally to bring into R and finish the ABC protocol. 

4. Run ABC on the combined simSAS output table, using R library 'abc'. More details in "run_ABC.R"