from matplotlib import pyplot as plt
import os
import re
import SLOTH.sloth.IO
# ==============================================================================
# CATEGORY 1: FILE I/O AND DATA READING
# Functions for reading simulation outputs, settings files, and KINSOL logs
# ==============================================================================


# Configure matplotlib defaults
plt.rcParams.update({'font.size': 22})

# Path constants
DIRPATH = os.path.dirname(os.path.realpath(__file__))
SCRATCHPATH = "/p/scratch/cslts/miaari1"



def read_settings_file(case_dir):
    """
    Read simulation settings from a settings.txt file.
    
    Parameters:
    -----------
    case_dir : str
        Path to the case directory containing settings.txt
    
    Returns:
    --------
    dict
        Dictionary containing simulation parameters:
        {k, q, alfa, n, theta_r, theta_s, d}
    """
    settings_path = os.path.join(case_dir, "settings.txt")
    
    with open(settings_path, "r") as f:
        settings = f.readlines()
    
    # Parse settings (skip header line)
    return {
        "k": float(settings[1].strip()),
        "q": float(settings[2].strip()),
        "alfa": float(settings[3].strip()),
        "n": float(settings[4].strip()),
        "theta_r": float(settings[5].strip()),
        "theta_s": float(settings[6].strip()),
        "d": float(settings[7].strip())
    }


def get_toplayer_pressure(pfb_path):
    """
    Extract the top layer pressure from a ParFlow binary file.
    
    Parameters:
    -----------
    pfb_path : str
        Path to the .pfb file
    
    Returns:
    --------
    float
        Surface pressure value at the top layer
    """
    data = SLOTH.sloth.IO.read_pfb(pfb_path)
    surface_pressure = data[-1, 0, 0]
    return surface_pressure


def parse_kinsol_log(kinsolfile_path):
    """
    Parse KINSOL log file to extract steady-state convergence time.
    
    The function looks for the pattern where KINSOL achieves steady state,
    indicated by consecutive time steps without reinitialization.
    
    Parameters:
    -----------
    kinsolfile_path : str
        Path to the KINSOL log file
    
    Returns:
    --------
    float or None
        Time (in hours) when steady state was reached, or None if not reached
    """
    if not os.path.exists(kinsolfile_path):
        return None
    
    with open(kinsolfile_path, "r") as f:
        lines = f.readlines()
    
    # Extract lines containing KINSOL time step information
    time_kinsol = [
        line for line in lines 
        if "KINSOL starting step for time" in line or "KINSolInit nni=    0  fnorm=" in line
    ]
    
    # Find the first occurrence where steady state is reached
    # (no reinitialization between consecutive time steps)
    for i in range(1, len(time_kinsol)):
        current_is_init = "KINSolInit nni=    0  fnorm=" in time_kinsol[i]
        previous_is_init = "KINSolInit nni=    0  fnorm=" in time_kinsol[i-1]
        
        if not current_is_init and not previous_is_init:
            # Steady state reached - extract time from previous entry
            time_matches = re.findall(r"\d+\.\d+", time_kinsol[i-1])
            if time_matches:
                return float(time_matches[0])
    
    return None


def write_settings_file(case_dir, settings_dict):
    """
    Write simulation settings to a settings.txt file.
    
    Parameters:
    -----------
    case_dir : str
        Path to the case directory
    settings_dict : dict
        Dictionary containing simulation parameters
    """
    settings_lines = [
        "k,q,alfa,n,theta_r,theta_s,d",
        f"{settings_dict['k']}",
        f"{settings_dict['q']}",
        f"{settings_dict['alfa']}",
        f"{settings_dict['n']}",
        f"{settings_dict['theta_r']}",
        f"{settings_dict['theta_s']}",
        f"{settings_dict['d']}"
    ]
    
    settings_path = os.path.join(case_dir, "settings.txt")
    with open(settings_path, 'w') as f:
        f.write('\n'.join(settings_lines))
