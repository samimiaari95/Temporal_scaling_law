from file_io import DIRPATH, SCRATCHPATH, write_settings_file
import os
from utils import make_dir
import random
import pandas as pd

# ==============================================================================
# CATEGORY 2: SIMULATION SETUP
# Functions for creating ParFlow namelist files and simulation configurations
# ==============================================================================

def update_namelist_parameters(namelist, k, q, alfa, n, theta_r, theta_s, d):
    """
    Update ParFlow namelist with new parameter values.
    
    This function modifies specific lines in the ParFlow TCL namelist to update
    hydraulic parameters for van Genuchten model.
    
    Parameters:
    -----------
    namelist : list
        List of strings representing the namelist file lines
    k : float
        Hydraulic conductivity (m/hr)
    q : float
        Infiltration/exfiltration rate (m/hr)
    alfa : float
        van Genuchten alpha parameter (1/m)
    n : float
        van Genuchten n parameter (-)
    theta_r : float
        Residual water content (-)
    theta_s : float
        Saturated water content (-)
    d : float
        Soil depth (m)
    
    Returns:
    --------
    list
        Updated namelist with new parameter values
    """
    updated_namelist = namelist.copy()
    
    # Line indices for different parameters (hardcoded based on template structure)
    k_line = 64
    q_line = 214
    d_line = 54
    nz_line = 31
    alfa_line1 = 142
    alfa_line2 = 151
    n_line1 = 143
    n_line2 = 152
    porosity_line = 129
    sres_line = 153
    
    # Update hydraulic conductivity
    updated_namelist[k_line] = namelist[k_line].replace("0.01", f"{k}")
    
    # Update infiltration/exfiltration rate
    updated_namelist[q_line] = namelist[q_line].replace("0.001", f"{q}")
    
    # Update van Genuchten parameters
    updated_namelist[alfa_line1] = namelist[alfa_line1].replace("1.0", f"{alfa}")
    updated_namelist[alfa_line2] = namelist[alfa_line2].replace("1.", f"{alfa}")
    updated_namelist[n_line1] = namelist[n_line1].replace("2.0", f"{n}")
    updated_namelist[n_line2] = namelist[n_line2].replace("2.0", f"{n}")
    
    # Update soil hydraulic properties
    updated_namelist[porosity_line] = namelist[porosity_line].replace("0.5", f"{theta_s}")
    updated_namelist[sres_line] = namelist[sres_line].replace("0.2", f"{theta_r/theta_s}")
    
    # Update domain geometry
    updated_namelist[d_line] = namelist[d_line].replace("4.0", f"{d}")
    updated_namelist[nz_line] = namelist[nz_line].replace("40", f"{int(d/0.1)}")
    
    return updated_namelist


def create_test_case(case_dir, namelist_template, settings_dict, tcl_filename="infiltration.tcl"):
    """
    Create a complete test case directory with TCL file and settings.
    
    Parameters:
    -----------
    case_dir : str
        Path to the test case directory to create
    namelist_template : list
        Template namelist to modify
    settings_dict : dict
        Dictionary with simulation parameters
    tcl_filename : str
        Name of the output TCL file (default: "infiltration.tcl")
    """
    # Create case directory
    make_dir(case_dir)
    
    # Update namelist with parameters
    updated_namelist = update_namelist_parameters(
        namelist_template,
        settings_dict['k'],
        settings_dict['q'],
        settings_dict['alfa'],
        settings_dict['n'],
        settings_dict['theta_r'],
        settings_dict['theta_s'],
        settings_dict['d']
    )
    
    # Write settings file
    write_settings_file(case_dir, settings_dict)
    
    # Write TCL file
    tcl_path = os.path.join(case_dir, tcl_filename)
    with open(tcl_path, 'w') as f:
        f.write('\n'.join(updated_namelist))


def simulation_bash():
    """
    Generate bash script for running multiple ParFlow simulations.
    
    Creates a batch script that runs infiltration and exfiltration
    simulations for 100 test cases.
    """
    filepath = os.path.join(DIRPATH, "parflow_simulations.sh")
    
    # Read existing script header
    with open(filepath, "r") as f:
        sbatch = f.readlines()
    
    # Keep only the header (first 19 lines)
    sbatch = sbatch[:19]
    
    # Add commands for each test case
    run_command = "tclsh infiltration.tcl"
    run_command_2 = "tclsh exfiltration.tcl"
    
    for i in range(100):
        case_path = f"cd /p/scratch/cslts/miaari1/infexfcases/test_case{i}"
        sbatch.append(case_path)
        sbatch.append(run_command)
        sbatch.append(run_command_2)
    
    # Write updated script
    with open(filepath, 'w') as f:
        f.write('\n'.join(sbatch))

def parflow_namelist_inex():
    """
    Generate ParFlow namelist files for infiltration and exfiltration simulations.
    
    Creates test cases with both infiltration and drainage (exfiltration) phases.
    Used for studying temporal scaling relationships in both wetting and drying.
    
    Generates test cases 0-999 with low infiltration rates (q/k = 0.0001-1).
    """
    dir_path = os.path.join(SCRATCHPATH, "testcases")
    inf_filepath = os.path.join(DIRPATH, "inputs", "infiltration.tcl")
    exf_filepath = os.path.join(DIRPATH, "inputs", "exfiltration.tcl")
    vg_params_path = os.path.join(DIRPATH, "inputs", "WaterFlowParameters.csv")
    
    # Read namelist templates
    with open(inf_filepath, "r") as f:
        inf_namelist = f.readlines()
    
    with open(exf_filepath, "r") as f:
        exf_namelist = f.readlines()
    
    # Load parameters
    vg_params = pd.read_csv(vg_params_path)
    vg_indexes = list(range(12))
    d_range = list(range(1, 11))
    
    # Generate test cases
    for i in range(1000):
        case_path = os.path.join(dir_path, f"test_case{i}")
        
        # Sample parameters
        vg_index = random.choice(vg_indexes)
        alfa = vg_params["Alpha"].iloc[vg_index]
        n = vg_params["n"].iloc[vg_index]
        theta_r = vg_params["Sr"].iloc[vg_index]
        theta_s = vg_params["Ss"].iloc[vg_index]
        k = vg_params["Ks"].iloc[vg_index]
        
        # Low infiltration rate for exfiltration studies
        q_k_ratio = random.uniform(0.0001, 1)
        q = q_k_ratio * k
        d = random.choice(d_range)
        
        # Create settings
        settings_dict = {
            'k': k, 'q': q, 'alfa': alfa, 'n': n,
            'theta_r': theta_r, 'theta_s': theta_s, 'd': d
        }
        
        # Create infiltration case
        create_test_case(case_path, inf_namelist, settings_dict, "infiltration.tcl")
        
        # Create exfiltration case (no q in namelist, uses initial condition from infiltration)
        exf_namelist_updated = update_namelist_parameters(
            exf_namelist, k, 0, alfa, n, theta_r, theta_s, d
        )
        exf_tcl_path = os.path.join(case_path, "exfiltration.tcl")
        with open(exf_tcl_path, 'w') as f:
            f.write('\n'.join(exf_namelist_updated))


def parflow_namelist_inex_fixed():
    """
    Generate ParFlow namelist files with fixed parameter combinations.
    
    Creates a deterministic set of test cases by systematically varying:
    - Infiltration rate (q)
    - van Genuchten parameters (different soil types)
    - Soil depth (d)
    
    This provides a structured parameter space for analysis.
    """
    dir_path = os.path.join(SCRATCHPATH, "infexfcases")
    inf_filepath = os.path.join(DIRPATH, "inputs", "infiltration.tcl")
    exf_filepath = os.path.join(DIRPATH, "inputs", "exfiltration.tcl")
    vg_params_path = os.path.join(DIRPATH, "inputs", "WaterFlowParameters.csv")
    
    # Read templates
    with open(inf_filepath, "r") as f:
        inf_namelist = f.readlines()
    
    with open(exf_filepath, "r") as f:
        exf_namelist = f.readlines()
    
    # Load parameters
    vg_params = pd.read_csv(vg_params_path)
    vg_indexes = list(range(12))
    
    # Define parameter combinations
    q_vg_list = {
        "0.001": vg_indexes,
        "0.0001": vg_indexes,
        "0.00015": vg_indexes
    }
    dlist = [1, 2, 3]
    
    # Generate test cases systematically
    case_index = 0
    for q, vg_indices in q_vg_list.items():
        for vg_idx in vg_indices:
            for d in dlist:
                # Extract parameters for this soil type
                alfa = vg_params["Alpha"].iloc[vg_idx]
                n = vg_params["n"].iloc[vg_idx]
                theta_r = vg_params["Sr"].iloc[vg_idx]
                theta_s = vg_params["Ss"].iloc[vg_idx]
                k = vg_params["Ks"].iloc[vg_idx]
                
                case_path = os.path.join(dir_path, f"test_case{case_index}")
                
                # Create settings
                settings_dict = {
                    'k': k, 'q': float(q), 'alfa': alfa, 'n': n,
                    'theta_r': theta_r, 'theta_s': theta_s, 'd': d
                }
                
                # Create both infiltration and exfiltration cases
                create_test_case(case_path, inf_namelist, settings_dict, "infiltration.tcl")
                
                exf_namelist_updated = update_namelist_parameters(
                    exf_namelist, k, 0, alfa, n, theta_r, theta_s, d
                )
                exf_tcl_path = os.path.join(case_path, "exfiltration.tcl")
                with open(exf_tcl_path, 'w') as f:
                    f.write('\n'.join(exf_namelist_updated))
                
                case_index += 1
