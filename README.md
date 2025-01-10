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

(1) Solvation Data Loading:
  - Solvation data for stable molecules, free radicals, and transition states.
  - Solvation free energies of activation for H-abstraction reactions.
  - 
(2) Input Parameters:
  - Critical constants,
  - acentric factors,
  - Twu-91 alpha function parameters.
  - COSMO-RS sigma-profiles.
  - Functional groups and their contribution parameters.
  - Peng-Robinson Equation of State (PR EoS) and COSMO-RS universal constants.

(3) Calculations:
  - Perform simulations of solvation free energy data
  - Optimizations (contribution constants of function groups) followed by simulations

(4) Outputs:
  - Results stored in solvation-data/solv-energy-res.dat.
  - Visual plots and deviation metrics for analysis.

---------------------------------------
Requirements
MATLAB installed with necessary toolboxes.
Ensure all required files are available:
- source/ directory containing all functions.
- solvation-data/ESOLV_DATA.xlsx for solvation data.
- properties/Compounds_data.xlsm for pure compound inputs.
- sigma-profiles for to access the .dat files of sigma profiles.

---------------------------------------
Usage Instructions

(1) Setup: 
Edit the Excel file INPUTS.xlsx to define settings for your calculations.
- Type of calculation = similation only or optimization + simulation
- Data to simulate = Molecules, Free Radicals, Transition-States, or Solvation Free Energy of Activation
- Level of theory = only BP-TZVPD-FINE is available in this implementation
- Sigma-profiles (solvents) = chose between Quantum-based sigma-profiles or sigma-profiles predicted by group contribution
- Sigma-profiles (solutes)  = ibid
- method for EoS parameters (solvents) = chose between experimental data or group contribution predictions to calculate the Peng-Robinson EoS parameters
- method for EoS parameters (solutes)  = ibid
- Alpha-function = Alpha-function to be used with the Peng-Robinson EoS (Twu-91 or SOave 1972)
You can also modify the universal constants of the mixing rule of the Peng-Robinson equation, as well as the universal constants of the COSMO-RS model
To do this, select either 'Parametrization_Mix_Rule' or 'Parametrization_COSMO' worksheets

(2) Place all necessary files in the same directory structure as described.

Run the Script:

(3) Execute the script in MATLAB by running MAIN_ESOLV.m.
The script will prompt a wait bar to follow the progress of the chosen calculation

(4) Results will be saved in solvation-data/solv-energy-res.dat.
Plots will be generated for molecules, free radicals, transition states, or activation free energies, depending on the input.

---------------------------------------
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
