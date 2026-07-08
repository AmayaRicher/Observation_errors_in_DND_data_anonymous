#Title: Rare but pervasive: false-positive errors bias estimates of fine-scale occupancy probabilities in plant population surveys

#Files
This folder contains:
-a file to run analysis and generate figures and tables from the main text (except simulation analysis): *DND_data_Analysis_and_Figures.qmd*
-a file to run analysis and generate figures and tables from the main text (only simulation analysis): *DND_data_Simulations.R*
-a file containing Supplementary material: *DND_data_SuppMat.qmd*, along with the PDF version *DND_data_SuppMat.pdf*
-a *data* folder containing the raw data
-a *pictures* folder containing pictures for Supplementary material
-a *saved_tables* folder containing the saved outputs from running the *DND_data_Analysis_and_Figures.qmd* script
-a *simulation_outputs* folder containing the outputs from the simulation analysis in the *DND_data_Simulation.R* script.
Note that DND stands short for detection/non-detection.

# Running the code
First, run the *DND_data_Analysis_and_Figures.qmd* script. This script will load the data from the *data* folder, save formatted data in the *saved_tables* folder, fit the models and print the main figures.
Then, run the *DND_data_Simulation.R* script, which contains the code to simulate fine-scale DND data, run occupancy models and print main figures. Because running the code for simulations and model fitting takes a while, outputs of simulations are already saved and stored in *simulation_outputs*, if one does not desire to rerun simulation analysis.
