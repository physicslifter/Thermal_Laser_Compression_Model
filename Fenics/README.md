# Files for Fenics

## Functions_Notebook.py
Interactive notebook (to be worked w/ in VS code) that explored FEM functions in Fenics

## current_fits_animation.py
Script for visualizing the current fit from the optimization, and updates as the optimization goes on

## data_dict.json
Json file for storing the "real" data from the experiment(s)

## default_values_dict.json
json file containing default values for the FEM setup, including: final time, number of steps, initial temp, peak temp, a, b, k_1 (MgO), rho, c, MgO_length

## eval_to_gen.py
Contains functions for reading optimization_output.csv and rewriting the data to be in terms of the best fit for a generation, instead of evaluations, which optimization_output.csv is in

## func1.py
Base function for running the FEM

## gen.csv
Optimization data in terms of generation

## global_variables.json
JSON file for global variables to be used throughout the optimization process

## geometry_dict.json
json specifying the geometry of the finite element experiment

## img
subfolder for img documenting issues while building the model

## model_results.json
stores the output for the most recent fit

## optimization_output.csv
Output from the optimization

## optimize.py
Script for running the optimization

## opt_func.py
Wrapped function for setting up and running the optimization

## plot.py
Functions for plotting a single FEM run

##  plot_animation.py
Script for animating the progress of the optimization as it runs

## progress.txt
txt file for documenting where I am at in the build process

## prove_working.py
python script to prove that the functions are working as intended. (in notebook format, to be run inside of vs code)

## read_data.py
Function for reading in the data from excel files, cleaning it up, and saving it in json format

## run_model.py
Function for running the FEM model. References Func1.py

## wrapped_run.py
Functions for getting goodness of fit assessment (LEAST SQUARES IN OUR CASE))

## write_dicts.py
Script for writing/cleaning data from excel files into JSON format
