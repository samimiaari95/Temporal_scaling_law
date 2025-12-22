# Temporal Scaling Law Analysis

This project analyzes temporal scaling relationships in soil water infiltration and exfiltration (drainage) processes using ParFlow simulations.

## Overview

The workflow investigates how soil hydraulic properties, infiltration rates, and soil depth affect the time to reach steady-state conditions during infiltration and subsequent drainage. The analysis derives power-law scaling relationships between non-dimensional parameters.

## Running the Workflow

### Run all steps:
```bash
python main.py --all
```

### Run individual steps:
```bash
python main.py --step 1  # Setup only
python main.py --step 3  # Data processing only
python main.py --step 4  # Plotting only
```

**Note:** Step 2 (ParFlow simulations) must be run separately on the HPC system using the generated batch script.

---

## Key Scientific Outputs

1. **Infiltration scaling law:** λ/SST_i ∝ (α·d·q/K_s)^b
2. **Drainage scaling law:** λ/(SST_d·K_s) ∝ (α·d)^b
3. **Time ratio relationships:** SST_d/SST_i depends on λ and prior infiltration rate
4. **Depth dependencies:** Both infiltration and drainage times scale with soil depth
5. **Infiltration rate effects:** Power-law exponent varies logarithmically with infiltration rate

---

## Main Workflow (main.py)

The `main.py` script orchestrates a 4-step workflow that can be run sequentially or individually:

### Step 1: Simulation Setup
**Function:** `simulation_setup.parflow_namelist_inex()`, `simulation_setup.simulation_bash()`, `simulation_setup.parflow_namelist_inex_fixed()`

**What it does:**
- Generates ParFlow TCL namelist files for infiltration and exfiltration simulations
- Creates test cases with varying hydraulic parameters (van Genuchten parameters, infiltration rates, soil depths)
- Produces two sets of test cases:
  - **Random cases (0-999):** Random combinations of soil types and parameters for broad coverage
  - **Fixed cases (0-99):** Systematic combinations for controlled analysis
- Generates a bash script to run all ParFlow simulations

**Inputs:**
- Template TCL files: `inputs/infiltration.tcl`, `inputs/exfiltration.tcl`
- Van Genuchten parameters: `inputs/WaterFlowParameters.csv` (12 soil types)
- Parameter ranges:
  - Hydraulic conductivity (k): from parameter database
  - Infiltration rate (q): 0.0001-1.0 × k (random) or fixed values (0.0001, 0.00015, 0.001)
  - van Genuchten α (alfa): from parameter database (1/m)
  - van Genuchten n: from parameter database (-)
  - Residual water content (θ_r): from parameter database (-)
  - Saturated water content (θ_s): from parameter database (-)
  - Soil depth (d): 1-10 m (random) or 1,2,3 m (fixed)

**Outputs:**
- Test case directories: `/p/scratch/cslts/miaari1/testcases/test_case{0-999}`
- Test case directories (fixed): `/p/scratch/cslts/miaari1/infexfcases/test_case{0-99}`
- Each directory contains:
  - `infiltration.tcl`: ParFlow namelist for infiltration phase
  - `exfiltration.tcl`: ParFlow namelist for drainage phase
  - `settings.txt`: Parameter values for the case
- Bash script: `parflow_simulations.sh` (batch job script to run all simulations)

---

### Step 2: Run Simulations
**What it does:**
- **Note:** This step is performed externally on HPC system
- Executes ParFlow simulations using the generated TCL files
- Runs paired infiltration-exfiltration simulations:
  1. **Infiltration phase:** Applies constant infiltration rate until steady-state
  2. **Exfiltration phase:** Uses infiltration end state as initial condition, allows drainage until steady-state
- KINSOL solver logs convergence and steady-state detection

**Inputs:**
- `infiltration.tcl` and `exfiltration.tcl` from each test case directory
- ParFlow executable on HPC cluster
- Batch script: `parflow_simulations.sh`

**Outputs (per test case):**
- Pressure field files: `*.out.press.{timestep}.pfb` (ParFlow binary format)
- KINSOL log files:
  - `infiltration.out.kinsol.log`: Convergence log for infiltration
  - `exfiltration.out.kinsol.log`: Convergence log for drainage
- Other ParFlow output files (saturation, permeability fields)

**How to run:**
```bash
sbatch parflow_simulations.sh
```

---

### Step 3: Data Processing
**Function:** `data_processing.steadystate_kinsol_infiltrationdrainage()`

**What it does:**
- Parses KINSOL log files to extract steady-state times
- Identifies when infiltration and drainage reached steady-state conditions
- Reads simulation parameters from `settings.txt`
- Extracts top layer pressure from final drainage state
- Compiles all results into a single CSV database
- Filters out cases that didn't reach steady-state

**Inputs:**
- KINSOL logs: `{case_dir}/infiltration.out.kinsol.log`, `{case_dir}/exfiltration.out.kinsol.log`
- Settings: `{case_dir}/settings.txt`
- Pressure files: `{case_dir}/exfiltration.out.press.00000.pfb`

**Processing logic:**
- Steady-state detection: KINSOL achieves convergence without reinitialization between consecutive time steps
- For each test case:
  1. Parse infiltration log → extract `inf_time` (hours to steady-state)
  2. If infiltration succeeded, parse exfiltration log → extract `exf_time`
  3. Read van Genuchten parameters and domain properties
  4. Extract top layer pressure at end of drainage
  5. Store all data if both phases reached steady-state

**Outputs:**
- CSV file: `outputs/inf_exf_specialcases_toplayer.csv`
- Columns:
  - `case_index`: Test case number
  - `q`: Infiltration rate (m/hr)
  - `k`: Hydraulic conductivity (m/hr)
  - `d`: Soil depth (m)
  - `alfa`: van Genuchten α (1/m)
  - `n`: van Genuchten n (-)
  - `theta_r`: Residual water content (-)
  - `theta_s`: Saturated water content (-)
  - `inf_time`: Infiltration steady-state time (hr)
  - `exf_time`: Exfiltration steady-state time (hr)
  - `toplayer_pressure`: Top layer pressure head at end of drainage (m)

**Console output:**
- Number of successfully processed cases
- Cases where infiltration or exfiltration didn't reach steady-state

---

### Step 4: Analysis and Visualization
**Functions:** Multiple plotting functions from `plotting.py`

**What it does:**
Generates comprehensive visualizations and power-law analyses of the temporal scaling relationships.

#### 4a. Pressure Profile Evolution
**Function:** `plotting.plot_pressure_profile()`
- **Input:** Example pressure field files from `pressureprofileexample/` directory
- **Output:** `pressureprofileexample/pressure_profile.png`
- **Shows:** How pressure profile evolves over time during infiltration (every 10 time steps from 0-210)

#### 4b. Lambda Dependence Analysis
**Function:** `plotting.infexf_lambda_dependence()`
- **Input:** `inputs/drainage_inf_testcases_toplayer.csv`
- **Output:** `outputs/exfinf_lambda.png`
- **Shows:** Ratio of drainage to infiltration time vs. characteristic length scale λ = d - |ψ_top|
- **Purpose:** Demonstrates that drainage takes longer than infiltration for the same soil

#### 4c. Individual Time Dependencies
**Functions:** `plotting.dexf_dependence()`, `plotting.dinf_dependence()`
- **Input:** `inputs/drainage_inf_testcases_toplayer.csv`
- **Outputs:** 
  - `outputs/tdr_vs_d.png` (drainage time vs. depth)
  - `outputs/tinf_vs_d.png` (infiltration time vs. depth)
- **Shows:** How steady-state times scale with soil depth

#### 4d. Infiltration Scaling Law
**Function:** `plotting.inf_solution_plots()`
- **Input:** `inputs/infiltration_pressure_index.csv`
- **Output:** `outputs/inf_adq-k_vs_v_symbolstest.png`
- **Analysis:** 
  - X-axis: α·d·q/K_s (non-dimensional infiltration parameter)
  - Y-axis: λ/SST_i (non-dimensional infiltration velocity)
  - Fits power law: y = a·x^b across all soil types
  - Different colors/markers for different hydraulic conductivities
- **Reports:** Power-law coefficients (a, b) and R² value

#### 4e. Drainage Scaling Law
**Function:** `plotting.exf_solution_plots()`
- **Input:** `inputs/drainage_inf_testcases_toplayer.csv`
- **Output:** `outputs/exf_ad_vs_v-k_symbolstest.png`
- **Analysis:**
  - X-axis: α·d (non-dimensional depth)
  - Y-axis: λ/(SST_d·K_s) (non-dimensional drainage velocity)
  - Fits power law across all soil types
  - Filters cases where k >= q
- **Reports:** Power-law coefficients and R²

#### 4f. Infiltration Rate Effects
**Function:** `plotting.qfit()`
- **Input:** `inputs/exf_intercept.csv`
- **Output:** `outputs/exf_intercept_qfit.png`
- **Analysis:** 
  - Fits separate power laws for different infiltration rates
  - Shows how prior infiltration rate affects subsequent drainage behavior
  - Annotates each fitted curve with equation

#### 4g. Infiltration Rate vs. Slope
**Function:** `plotting.q_vs_slope()`
- **Input:** `inputs/q_vs_fittinglinesslope2.csv`
- **Output:** `outputs/q_vs_slope2.png`
- **Analysis:**
  - Plots power-law slope vs. infiltration rate
  - Fits logarithmic model: slope = a·ln(q) + b
  - Demonstrates systematic relationship between infiltration rate and scaling behavior

---

## Supporting Modules

### file_io.py
Functions for reading/writing ParFlow files, settings files, and KINSOL logs:
- `read_settings_file()`: Parse settings.txt
- `get_toplayer_pressure()`: Extract surface pressure from .pfb files
- `parse_kinsol_log()`: Extract steady-state time from KINSOL logs
- `write_settings_file()`: Write parameter files

### simulation_setup.py
Functions for creating ParFlow configurations:
- `update_namelist_parameters()`: Modify TCL templates with new parameters
- `create_test_case()`: Generate complete test case directory
- `parflow_namelist_inex()`: Generate random parameter combinations
- `parflow_namelist_inex_fixed()`: Generate systematic parameter grid

### data_processing.py
Functions for extracting and compiling results:
- `steadystate_kinsol_infiltrationdrainage()`: Main data extraction function
- `steadystate_kinsol_drain_inf()`: Helper for drainage analysis

### analysis.py
Statistical analysis functions:
- `fitting_func()`: Power-law fitting using log-linearization

### plotting.py
Visualization functions (8 different plot types as described in Step 4)

### utils.py
Utility functions:
- `make_dir()`: Recursive directory creation
- `powerlaw_func()`: Power-law function y = a·x^b
- `linear_law()`: Linear function for log-space fitting
- `read_nc()`: Read NetCDF files
- `delete_files()`: Cleanup utilities


Disclaimer: this code was generated with support from AI