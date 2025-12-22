from file_io import DIRPATH, SCRATCHPATH, read_settings_file, get_toplayer_pressure, parse_kinsol_log
import os
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
    dir_path = os.path.join(SCRATCHPATH, "infexfcases")
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
    dir_path = os.path.join(SCRATCHPATH, "infexfcases")
    
    # Initialize data collection
    data = {
        "case_index": [], "q": [], "k": [], "d": [],
        "alfa": [], "n": [], "theta_r": [], "theta_s": [],
        "inf_time": [], "exf_time": [], "toplayer_pressure": []
    }
    
    # Process each test case
    num_cases = len(os.listdir(dir_path))
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
    
    # Save to CSV
    df = pd.DataFrame(data)
    print(f"Processed {len(data['exf_time'])} cases successfully")
    print(df.head())
    
    output_path = os.path.join(DIRPATH, "outputs", "inf_exf_specialcases_toplayer.csv")
    df.to_csv(output_path, index=False)
