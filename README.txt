R_code
by Jussi Mäkinen
started 10-29-2020

Content:
R_code for running SDMs for few hummingbird species
SDMs are run by dependent r scripts that are for specific analysis steps

Scripts:
SDM
- set species, environment, offsets and inference method
- call Run_analysis for building SDM and visualizing results
- return model characteristics and validation results
- print results with text-file containing model information

Run_analysis
- compile data by calling take_data scripts
- run inference by calling inference method -script
- visualize results by calling inference method specific visualization-script
- return validation results and prediction maps