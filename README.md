Solvation Data Simulation and Optimization Script
This MATLAB script is designed for performing solvation free energy calculations, enabling users to conduct simulations or optimize parameters for solvation free energy data. 
Below is an overview of its functionality, structure, and usage instructions.
---------------------------------------
Author
Francisco Carlos Paes
January 2025
Equipe Thermodynamique et Energie (ThermE)
Laboratoire Réactions et Génie des Procédés (LRGP)
UMR 7274 CNRS - Université de Lorraine
---------------------------------------
Features

Solvation Data Loading:
  - Solvation data for stable molecules, free radicals, and transition states.
  - Solvation free energies of activation for H-abstraction reactions.
  
Input Parameters:
  Critical constants, acentric factors, Twu-91 alpha function parameters.
  COSMO-RS sigma-profiles.
  Functional groups and their contribution parameters.
  Peng-Robinson Equation of State (PR EoS) and COSMO-RS universal constants.

Calculations:

Perform either simulations or optimizations followed by simulations.
Generate results as parity plots and deviation distribution plots.
Outputs:

Results stored in solvation-data/solv-energy-res.dat.
Visual plots for analysis.
Requirements
MATLAB installed with necessary toolboxes.
Ensure all required files are available:
source/ directory containing relevant functions.
solvation-data/ESOLV_DATA.xlsx for solvation data.
properties/Compounds_data.xlsm for compound properties.
sigma-profiles/p_sigma_res.xlsx for sigma profiles.
Usage Instructions
Setup:

Place all necessary files in the same directory structure as described.
Edit the Excel file INPUTS.xlsx to define settings for your calculations.
Run the Script:

Execute the script in MATLAB by running main_script.m.
The script will prompt for calculation options:
Simulation only: Directly simulate solvation energies using existing parameters.
Optimization + simulation: Optimize group contributions before simulation.
None: Exit the program.
Output:

Results will be saved in solvation-data/solv-energy-res.dat.
Plots will be generated for molecules, free radicals, transition states, or activation free energies, depending on the input.
Customization:

Modify the input files to include new solutes or parameters.
Adjust script sections for specific calculation methods or visualizations.
Key Functions
READ_OPTIONS: Reads and loads user-defined options.
PARAMETRIZATION_COSMO: Loads COSMO-RS and PR EoS parameters.
LOAD_GCM: Loads group contribution models.
SIMULATION and SIMULATION_DDG: Performs solvation energy calculations.
OPTIMIZATION_GC: Optimizes group contribution parameters.
Plotting Functions:
PLOT_MOLECULES, PLOT_RADICALS, PLOT_TS, PLOT_ACTIVATION.
Notes
Ensure MATLAB has access to all required directories and files before running the script.
To optimize specific molecule types, follow prompts during the Optimization + simulation step.
Check solvation-data/solv-energy-res.dat for detailed results.
Contact
For questions or issues, please contact Francisco Carlos Paes at the Thermodynamique et Energie (ThermE) team, LRGP.

Enjoy simulating!
