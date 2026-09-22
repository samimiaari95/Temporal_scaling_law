import shutil

from file_io import DIRPATH, SCRATCHPATH, read_settings_file, get_toplayer_pressure, parse_kinsol_log
import os
import numpy as np
import pandas as pd

# ==============================================================================
# CATEGORY 3: DATA PROCESSING AND EXTRACTION
# Functions for extracting steady-state times and processing simulation results
# ==============================================================================

def steadystate_kinsol_drain_inf(case_index):
    """
    Extract infiltration steady-state time for drainage analysis.
    
    Helper function for drainage simulations that need to know when
    infiltration reached steady state.
    
    Parameters:
    -----------
    case_index : int
        Index of the test case
    
    Returns:
    --------
    float or None
        Infiltration steady-state time, or None if not reached
    """
    # dir_path = os.path.join(SCRATCHPATH, "infexfcases")
    dir_path = os.path.join(SCRATCHPATH, "tolerance_sensitivity", "residualtol10-5")
    case_dir = os.path.join(dir_path, f"test_case{case_index}")
    kinsolfile = os.path.join(case_dir, "infiltration.out.kinsol.log")
    
    ss_time = parse_kinsol_log(kinsolfile)
    
    if ss_time is None:
        print(f"Infiltration not reached: case {case_index}")
    
    return ss_time


def steadystate_kinsol_infiltrationdrainage():
    """
    Extract steady-state times from infiltration-exfiltration simulations.
    
    Processes paired infiltration/exfiltration simulations to extract:
    - Infiltration steady-state time
    - Exfiltration (drainage) steady-state time
    - Simulation parameters
    - Top layer pressure after drainage
    
    This data is used to analyze temporal scaling in drainage processes.
    """
    # dir_path = os.path.join(SCRATCHPATH, "infexfcases")
    dir_path = os.path.join(SCRATCHPATH, "tolerance_sensitivity", "residualtol10-5")
    
    # Initialize data collection
    data = {
        "case_index": [], "q": [], "k": [], "d": [],
        "alfa": [], "n": [], "theta_r": [], "theta_s": [],
        "inf_time": [], "exf_time": [], "toplayer_pressure": []
    }
    
    # Process each test case
    num_cases = len(os.listdir(dir_path))
    count = 0
    for case_index in range(num_cases):
        print(f"Processing case {case_index}")
        
        case_dir = os.path.join(dir_path, f"test_case{case_index}")
        
        # First check if infiltration reached steady state
        inf_time = steadystate_kinsol_drain_inf(case_index)
        
        if inf_time is not None:
            # Process exfiltration
            exf_kinsolfile = os.path.join(case_dir, "exfiltration.out.kinsol.log")
            exf_time = parse_kinsol_log(exf_kinsolfile)
            
            if exf_time is not None:
                # Read settings
                settings = read_settings_file(case_dir)
                
                # Get top layer pressure after drainage
                press_file = os.path.join(case_dir, "exfiltration.out.press.00000.pfb")
                toplayer_pressure = get_toplayer_pressure(press_file)
                
                # Store results
                data["case_index"].append(case_index)
                data["inf_time"].append(inf_time)
                data["exf_time"].append(exf_time)
                data["toplayer_pressure"].append(toplayer_pressure)
                
                for key in ["k", "q", "alfa", "n", "theta_r", "theta_s", "d"]:
                    data[key].append(settings[key])
            else:
                print(f"Exfiltration not reached: case {case_index}")
        if inf_time is None or exf_time is None:
            count += 1
            # write to txt file the case_index of cases that did not reach steady state for exfiltration
            with open("unreached_exfiltration_cases.txt", "a") as f:
                f.write(f"{case_index}\n")
                f.close()
                
    print(count)
    # Save to CSV
    df = pd.DataFrame(data)
    print(f"Processed {len(data['exf_time'])} cases successfully")
    print(df)
    
    output_path = os.path.join(DIRPATH, "inputs", "tolerance_05_inf_exf_times.csv")
    df.to_csv(output_path, index=False)

def fixfailedcases():
    dir_path = os.path.join(SCRATCHPATH, "infexfcases")
    with open(os.path.join(DIRPATH, "unreached_exfiltration_cases.txt"), "r") as f:
        case_indices = [int(line.strip()) for line in f.readlines()]
    #drop duplicates in case_indices
    case_indices = list(set(case_indices))
    case_indices.sort()
    print(f"Fixing {len(case_indices)} cases that did not reach steady state for exfiltration...")

    for case_index in case_indices:
        print(f"Fixing case {case_index}...")
        case_dir = os.path.join(dir_path, f"test_case{case_index}")
        settings_path = os.path.join(case_dir, "settings.txt")
        with open(settings_path, "r") as f:
            settings_lines = f.readlines()

        if len(settings_lines) < 3:
            print(f"Skipping case {case_index}: malformed settings.txt")
            continue

        d_value = float(settings_lines[-1].strip())
        new_d_value = np.random.choice([1, 2, 3, 4])  # Randomly select a new d value

        print(settings_lines)
        settings_lines[-1] = f"{new_d_value}\n"
        print(d_value, new_d_value)
        print(settings_lines)

        with open(settings_path, "w") as f:
            f.writelines(settings_lines)

        # Update infiltration.tcl: replace the pressure line with negative q (preserve indentation)
        tcl_path = os.path.join(case_dir, "infiltration.tcl")
        with open(tcl_path, "r") as f:
            tcl_lines = f.readlines()

        d_line1_prefix = "pfset ComputationalGrid.NZ                      "
        d_line2_prefix = "pfset Geom.domain.Upper.Z                         "
        updated = False
        for line_index, line in enumerate(tcl_lines):
            if line.startswith(d_line1_prefix):
                tcl_lines[line_index] = f"{d_line1_prefix}{int(new_d_value*10)}\n"
                updated = True
            elif line.startswith(d_line2_prefix):
                tcl_lines[line_index] = f"{d_line2_prefix}{int(new_d_value)}\n"
                updated = True
                break

        if not updated:
            print(f"Warning: d line not found in case {case_index}")
            continue

        with open(tcl_path, "w") as f:
            f.writelines(tcl_lines)

        # TODO do the same for exfiltration
        exf_path = os.path.join(case_dir, "exfiltration.tcl")
        with open(exf_path, "r") as f:
            exf_lines = f.readlines()
        updated = False
        for line_index, line in enumerate(exf_lines):
            if line.startswith(d_line1_prefix):
                exf_lines[line_index] = f"{d_line1_prefix}{int(new_d_value*10)}\n"
                updated = True
            elif line.startswith(d_line2_prefix):
                exf_lines[line_index] = f"{d_line2_prefix}{int(new_d_value)}\n"
                updated = True
                break

        if not updated:
            print(f"Warning: d line not found in case {case_index}")
            continue

        with open(exf_path, "w") as f:
            f.writelines(exf_lines)


        # create a new bash script to run all selected cases with longer time for drainage

def select1000cases_andrerunlongertimesfordrainage():
    dir_path = os.path.join(SCRATCHPATH, "infexfcases")
    num_cases = len(os.listdir(dir_path))
    selected_cases = np.random.choice(num_cases, size=1000, replace=False)
    
    # read all settings of selected cases and make a list for each parameter and then plot the distribution of each parameter to see if there are any patterns in the cases that did not reach steady state for drainage. This can help inform which cases to re-run with longer times.
    allsettings = {key: [] for key in ["k", "q", "alfa", "n", "theta_r", "theta_s", "d"]}

    count = 0
    for case_index in selected_cases:
        case_dir = os.path.join(dir_path, f"test_case{case_index}")
        settings = read_settings_file(case_dir)
        for key in allsettings:
            allsettings[key].append(settings[key])

        exf_kinsolfile = os.path.join(case_dir, "exfiltration.out.kinsol.log")
        exf_time = parse_kinsol_log(exf_kinsolfile)
        if exf_time is None:
            print(f"Re-running case {case_index} with longer time for drainage...")
            count += 1
            # Here you would implement the logic to modify the ParFlow input files
            # to allow for a longer simulation time and then re-run the simulation.
            # modify the namelist line "MaxIter" in the exfiltration namelist to a larger value (e.g., 20000) and then run the simulation using the bash script generated in the setup step.
            
            # create a new bash script to run all selected cases with longer time for drainage
        
        # save allsettings to a csv for later analysis
    print(count)
    allsettings_df = pd.DataFrame(allsettings)
    allsettings_output_path = os.path.join(DIRPATH, "outputs", "selected_1000cases_settings.csv")
    # allsettings_df.to_csv(allsettings_output_path, index=False)
