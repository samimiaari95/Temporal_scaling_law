from file_io import DIRPATH
from analysis import fitting_func
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import SLOTH.sloth.IO
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# ==============================================================================
# CATEGORY 4: PLOTTING AND VISUALIZATION
# Functions for creating plots and figures
# ==============================================================================

def plot_ss_profiles():
    """
    Plot steady-state pressure profiles for different q/k ratios.
    
    Reads pressure profile .pfb files from ss_profilesplot directory
    and creates overlaid profiles showing the effect of infiltration rate
    relative to hydraulic conductivity.
    """
    plt.rcParams.update({'font.size': 22})
    folderpath = os.path.join(DIRPATH, 'ss_profilesplot')
    
    # Create depth array (4m depth, 0.1m resolution)
    z = list(range(41, 1, -1))
    z = [x * 0.1 for x in z]
    
    # Get all pressure profile files
    files = [os.path.join(folderpath, x) for x in os.listdir(folderpath)]
    
    # Plot each profile
    for file in files:
        # Parse filename to extract parameters
        name = os.path.basename(file).replace(".pfb", "")
        namelist = name.split("_")
        q = float(namelist[0].replace("-", "."))
        k = float(namelist[1].replace("-", "."))
        t = namelist[2]
        qk_ratio = q / k
        
        # Read and plot pressure profile
        data = SLOTH.sloth.IO.read_pfb(file)
        plt.plot(data[:, 0, 0], z, color="black")
        plt.annotate(f"q/k={qk_ratio:.3f}", 
                    xy=(data[-1, 0, 0], z[-1]), 
                    color="black")
    
    # Format plot
    plt.gca().invert_yaxis()
    plt.xlim([-4, 0.2])
    plt.xlabel("Pressure (m)")
    plt.ylabel("Soil depth (m)")
    
    # Save figure
    output_path = os.path.join(DIRPATH, 'ss_profilesplot', 'ss_profiles.png')
    plt.savefig(output_path)
    plt.close()


def plot_pressure_profile():
    """
    Plot temporal evolution of pressure profile during infiltration.
    
    Creates a figure showing how the pressure profile develops over time
    from initial condition to steady state. Useful for visualizing the
    infiltration front propagation.
    """
    plt.rcParams.update({'font.size': 22})
    fig, ax = plt.subplots(figsize=(16, 9))
    
    # Path to example simulation
    name = os.path.join(DIRPATH, "pressureprofileexample", "infiltration")
    
    # Create depth array
    pressures = {}
    pressures["z"] = [-1 * z / 10 for z in range(0, 40, 1)]
    
    # Plot profiles at different times
    for t in range(0, 210, 10):
        pfb_file = name + '.out.press.' + ('{:05d}'.format(t)) + '.pfb'
        print(pfb_file)
        
        data = SLOTH.sloth.IO.read_pfb(pfb_file)
        plt.plot(data[:, 0, 0], list(reversed(pressures['z'])), color="black")
        
        pressures[f"time={t}"] = data[:, 0, 0]
    
    # Format plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Pressure head (m)")
    ax.set_ylabel("Depth (m)")
    
    # Save figure
    output_path = os.path.join(DIRPATH, "pressure_profile.png")
    fig.savefig(output_path)
    plt.close()


def infexf_lambda_dependence():
    """
    Plot drainage/infiltration time ratio vs. characteristic length scale.
    
    Analyzes the relationship between the drainage-to-infiltration time ratio
    and the characteristic length scale (lambda = d - |psi_top|).
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Column names
    d_column = "d"
    inf_t_column = "inf_time"
    exf_t_column = "exf_time"
    toplayer_pressure_column = "toplayer_pressure"
    
    # Create figure
    plt.figure(figsize=(16, 9))
    
    # Plot each depth group with different color
    for d_val, group in df.groupby(d_column):
        y = group[exf_t_column] / group[inf_t_column]
        x = group[d_column] - 0.05 - abs(group[toplayer_pressure_column])
        plt.scatter(x, y, label=f"d = {d_val}")
    
    # Format plot
    plt.grid(True)
    # plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("λ (m)")
    plt.ylabel(r"$t_{d} / t_{i}  (-)$")
    plt.legend(title="d values")
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_path, "exfinf_lambda.png"), dpi=300)
    plt.close()


def infexf_dependence():
    """
    Plot drainage and infiltration times vs. characteristic length scale.
    
    Creates a dual y-axis plot showing both drainage and infiltration
    steady-state times as functions of the characteristic length scale.
    Different colors represent different soil depths.
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Column names
    d_column = "d"
    inf_t_column = "inf_time"
    exf_t_column = "exf_time"
    toplayer_pressure_column = "toplayer_pressure"
    
    # Create figure with twin y-axes
    fig, ax1 = plt.subplots(figsize=(16, 9))
    ax2 = ax1.twinx()
    
    # Create color map for different depths
    unique_d = sorted(df[d_column].unique())
    colors = plt.cm.viridis(np.linspace(0, 1, len(unique_d)))
    
    # For combined legend
    handles = []
    labels = []
    
    # Plot each depth group
    for color, d_val in zip(colors, unique_d):
        group = df[df[d_column] == d_val]
        x = group[d_column] - 0.05 - abs(group[toplayer_pressure_column])
        yd = group[exf_t_column]
        yi = group[inf_t_column]
        
        # Plot drainage time (primary y-axis)
        s1 = ax1.scatter(x, yd, color=color, alpha=0.6, 
                        label=fr"$SST_d$ (d={d_val} m)")
        
        # Plot infiltration time (secondary y-axis)
        s2 = ax2.scatter(x, yi, color=color, alpha=0.4, marker='x',
                        label=fr"$SST_i$ (d={d_val} m)")
        
        handles.append((s1, s2))
        labels.append(fr"d={d_val}")
    
    # Format plot
    ax1.set_xlabel("λ (m)")
    ax1.set_ylabel(r"$SST_{d}$ (hr)")
    ax1.set_yscale('log')
    ax1.grid(True, linestyle='--', alpha=0.3)
    
    ax2.set_ylabel(r"$SST_{i}$ (hr)")
    ax2.set_yscale('log')
    
    # Combined legend
    ax1.legend(handles, labels, loc='lower right', fontsize=8,
              handler_map={tuple: HandlerTuple(ndivide=None)})
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "exf_vs_inf.png"), dpi=300)
    plt.close()


def dexf_dependence():
    """
    Plot drainage steady-state time vs. soil depth.
    
    Shows how the time to reach steady state during drainage depends
    on the total soil depth.
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Extract data
    x = df["d"]
    y = df["exf_time"]
    
    # Create plot
    plt.figure()
    plt.scatter(x, y, color='black')
    plt.grid(True)
    plt.yscale("log")
    plt.xlabel("d (m)")
    plt.ylabel(r"$t_{d}$ (hr)")
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_path, r"tdr_vs_d.png"), dpi=300)
    plt.close()


def dinf_dependence():
    """
    Plot infiltration steady-state time vs. soil depth.
    
    Shows how the time to reach steady state during infiltration depends
    on the total soil depth.
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Extract data
    x = df["d"]
    y = df["inf_time"]
    
    # Create plot
    plt.figure()
    plt.scatter(x, y, color='black')
    plt.grid(True)
    plt.yscale("log")
    plt.xlabel("d (m)")
    plt.ylabel(r"$t_{i}$ (hr)")
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_path, r"tinf_vs_d.png"), dpi=300)
    plt.close()


def exf_solution_plots(solution=None):
    """
    Create comprehensive plot showing drainage time scaling relationships.
    
    Plots the non-dimensional drainage velocity (lambda / (SST_d * Ks))
    against the non-dimensional depth (alpha * d) for different soil types.
    Includes power-law fit to the data.
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Get unique soil types (by hydraulic conductivity)
    soil_types = sorted(df["k"].unique())
    
    # Create color map and marker styles
    soil_types_colors = plt.cm.jet(np.linspace(0, 1, len(soil_types)))
    colors_dic = {soil_types[i]: soil_types_colors[i] for i in range(len(soil_types))}
    
    list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5,2), "x"]
    markers = {soil_types[i]: list_markers[i] for i in range(len(soil_types))}
    
    print(f"Number of soil types: {len(soil_types)}")
    print(f"Soil types (Ks values): {soil_types}")
    
    # Create figure
    ax = plt
    ax.figure(figsize=(22, 9))
    
    # Collect all data for fitting
    x_all = []
    y_all = []
    
    # Plot each soil type
    for soil in soil_types:
        df_soil = df[df["k"] == soil]
        
        for i in range(len(df_soil)):
            q = df_soil["q"].iloc[i]
            k = df_soil["k"].iloc[i]
            d = df_soil["d"].iloc[i]
            alfa = df_soil["alfa"].iloc[i]
            exf_t = df_soil["exf_time"].iloc[i]
            toplayer_pressure = df_soil["toplayer_pressure"].iloc[i]
            
            # Filter: only cases where k >= q and exf_t is valid
            if k >= q and exf_t != 0:
                # Compute non-dimensional variables
                if solution == 1:
                    x, y, xlabel, ylabel = exfsolution_1(
                        q, k, d, df_soil["n"].iloc[i], alfa,
                        df_soil["theta_r"].iloc[i], df_soil["theta_s"].iloc[i],
                        df_soil["inf_time"].iloc[i], exf_t, toplayer_pressure
                    )
                elif solution == 2:
                    x, y, xlabel, ylabel = exfsolution_2(
                        q, k, d, df_soil["n"].iloc[i], alfa,
                        df_soil["theta_r"].iloc[i], df_soil["theta_s"].iloc[i],
                        df_soil["inf_time"].iloc[i], exf_t, toplayer_pressure
                    )
                else:
                    raise ValueError("Drainage solution must be identified as 1 or 2")
                
                x_all.append(x)
                y_all.append(y)
                
                ax.scatter(x, y, c=colors_dic[soil], cmap='jet', s=80, marker=markers[k])
    
    # Create legend for soil types
    for k in markers.keys():
        ax.scatter([], [], c=colors_dic[k], s=80, marker=markers[k], 
                  label=f"{np.around(k, 4)}")
    
    ax.legend(title="Ks (m/hr)", loc='upper right', bbox_to_anchor=(1.01, 0.8))
    
    # Perform power-law fit
    x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(x_all, y_all)
    
    plt.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=5)
    
    print(f"Power-law fit: y = {round(a_fit, 2)} * x^{round(b_fit, 2)}")
    print(f"R² = {round(r2, 2)}")
    
    # Format plot
    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    ax.xlabel(f"{xlabel}", fontsize=32)
    ax.ylabel(f"{ylabel}", fontsize=32)
    
    # Save figure
    figname = "exf_ad_vs_v-k_symbols.png" if solution == 2 else "exf_ad_vs_tka-thetasr_symbols.png"
    ax.savefig(os.path.join(output_path, figname))
    plt.close()

def exfsolution_1(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    """
    Alternative formulation for drainage scaling (not currently used).
    
    Computes:
    - x: alpha * d
    - y: (SST_d * Ks * alpha) / (theta_s - theta_r)
    
    This is an alternative non-dimensionalization that was explored but
    not used in the final analysis.
    """
    x = alfa * d
    y = (exf_t * k * alfa) / (theta_s - theta_r)
    xlabel = r'$\alpha\cdot d (-)$'
    ylabel = r'$\frac{t_{d}\cdot K_s\cdot \alpha}{\theta_s-\theta_r} (-)$'
    return x, y, xlabel, ylabel

def exfsolution_2(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    """
    Calculate non-dimensional variables for drainage time scaling.
    
    Computes:
    - x: alpha * d (non-dimensional depth)
    - y: lambda / (SST_d * Ks) (non-dimensional drainage velocity)
    
    Parameters:
    -----------
    q : float
        Infiltration rate (m/hr)
    k : float
        Hydraulic conductivity (m/hr)
    d : float
        Soil depth (m)
    n : float
        van Genuchten n parameter
    alfa : float
        van Genuchten alpha parameter (1/m)
    theta_r, theta_s : float
        Residual and saturated water contents
    inf_t, exf_t : float
        Infiltration and exfiltration steady-state times (hr)
    toplayer_pressure : float
        Pressure at top layer (m)
    
    Returns:
    --------
    tuple
        (x, y, xlabel, ylabel) for plotting
    """
    x = alfa * d
    y = ((d - 0.05) - abs(toplayer_pressure)) / (exf_t * k)
    xlabel = r'$\alpha\cdot d (-)$'
    ylabel = r'$v_{d} / K_{s} (-)$'
    return x, y, xlabel, ylabel


def inf_solution_plots():
    """
    Create comprehensive plot showing infiltration time scaling relationships.
    
    Plots the non-dimensional infiltration velocity (lambda / SST_i)
    against the non-dimensional infiltration parameter (alpha * d * q / Ks).
    Includes power-law fit to the data.
    """
    inf_cases_path = os.path.join(DIRPATH, "inputs", "infiltration_pressure_index.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(inf_cases_path)
    
    # Get unique soil types
    soil_types = sorted(df["k"].unique())
    
    # Create color map and markers
    soil_types_colors = plt.cm.jet(np.linspace(0, 1, len(soil_types)))
    colors_dic = {soil_types[i]: soil_types_colors[i] for i in range(len(soil_types))}
    
    list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5,2), "x"]
    markers = {soil_types[i]: list_markers[i] for i in range(len(soil_types))}
    
    print(f"Number of soil types: {len(soil_types)}")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 9))
    
    print("Plotting points...")
    x_all = []
    y_all = []
    
    # Plot each soil type
    for soil in soil_types:
        df_soil = df[df["k"] == soil]
        
        for i in range(len(df_soil)):
            q = df_soil["q"].iloc[i]
            k = df_soil["k"].iloc[i]
            d = df_soil["d"].iloc[i]
            alfa = df_soil["alfa"].iloc[i]
            inf_t = df_soil["time"].iloc[i]
            toplayer_pressure = df_soil["toplayer_pressure"].iloc[i]
            
            # Filter: only cases where k >= q
            if k >= q:
                # Compute non-dimensional variables
                x, y, xlabel, ylabel = infsolution(
                    q, k, d, df_soil["n"].iloc[i], alfa,
                    df_soil["theta_r"].iloc[i], df_soil["theta_s"].iloc[i],
                    inf_t, toplayer_pressure
                )
                
                x_all.append(x)
                y_all.append(y)
                
                ax.scatter(x, y, c=colors_dic[soil], cmap='jet', s=80, marker=markers[k])
    
    # Create legend
    for k in markers.keys():
        ax.scatter([], [], c=colors_dic[k], s=80, marker=markers[k],
                  label=f"{np.around(k, 4)}")
    
    ax.legend(title="Ks (m/hr)", loc='center right', bbox_to_anchor=(1.01, 0.41))
    
    # Perform power-law fit
    x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(x_all, y_all)
    
    plt.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=5)
    
    print(f"Power-law fit: y = {round(a_fit, 2)} * x^{round(b_fit, 2)}")
    print(f"R² = {round(r2, 2)}")
    print(f"Pearson r² = {round(r2, 2)}")
    
    # Format plot
    ax.grid(True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"{xlabel}", fontsize=32)
    ax.set_ylabel(f"{ylabel}", fontsize=32)
    plt.tight_layout()
    
    # Save figure
    fig.savefig(os.path.join(output_path, "inf_adq-k_vs_v_symbols.png"))
    plt.close()


def infsolution(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    """
    Calculate non-dimensional variables for infiltration time scaling.
    
    Computes:
    - x: alpha * d * q / Ks (non-dimensional infiltration parameter)
    - y: lambda / SST_i (non-dimensional infiltration velocity)
    
    Parameters:
    -----------
    q : float
        Infiltration rate (m/hr)
    k : float
        Hydraulic conductivity (m/hr)
    d : float
        Soil depth (m)
    n : float
        van Genuchten n parameter
    alfa : float
        van Genuchten alpha parameter (1/m)
    theta_r, theta_s : float
        Residual and saturated water contents
    inf_t : float
        Infiltration steady-state time (hr)
    toplayer_pressure : float
        Pressure at top layer (m)
    
    Returns:
    --------
    tuple
        (x, y, xlabel, ylabel) for plotting
    """
    x = alfa * d * q / k
    y = ((d - 0.05) - abs(toplayer_pressure)) / inf_t
    xlabel = r'$\frac{\alpha\cdot d\cdot q}{K_s} (-)$'
    ylabel = r'$v_{i} (m/hr)$'
    return x, y, xlabel, ylabel


def qfit():
    """
    Create detailed plot showing effect of infiltration rate on drainage scaling.
    
    For selected soil types, plots the drainage scaling relationship for
    different infiltration rates and fits power laws to each rate separately.
    Demonstrates how infiltration rate affects subsequent drainage behavior.
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "exf_intercept.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Select specific soil types to analyze
    soil_types = sorted(df["k"].unique())
    markers = {0.0045: "d"}  # Diamond marker for selected soil
    
    # Create figure
    ax = plt
    ax.figure(figsize=(16, 9))
    
    # Compute non-dimensional variables
    df["x"] = df["d"] * df["alfa"]
    df["y"] = ((df["d"] - 0.05) - abs(df["toplayer_pressure"])) / (df["exf_time"] * df["k"])
    
    # Annotation dictionary for fitted curves
    ann_ab = {
        "0.00010.002": [r'$0.225x^{1.356}$', (0.73, 0.234386417)],
        "0.000150.002": [r'$0.234x^{1.385}$', (0.73, 0.338135512)],
        "0.0010.0045": [r'$0.304x^{1.717}$', (2, 0.09250808557112719)],
        "0.00010.0045": [r'$0.125x^{1.306}$', (2, 0.05052146488857635)],
        "0.000150.0045": [r'$0.156x^{1.409}$', (2, 0.05871341875690252)],
        "0.00050.0045": [r'$0.253x^{1.634}$', (2, 0.08165277976910569)]
    }
    
    ind = 0
    colors = ["r", "g", "b", "darkorange"]
    
    # Plot each soil type
    for soil in soil_types:
        if soil == 0.002:
            continue
        
        df_soil = df[df["k"] == soil]
        ax.scatter(df_soil["x"], df_soil["y"], c="k", s=80, marker=markers[soil])
        
        # Fit and plot each infiltration rate separately
        for q in df_soil["q"].unique():
            df_q = df_soil[df_soil["q"] == q]
            
            # Fit power law
            x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(
                df_q["x"].tolist(), df_q["y"].tolist()
            )
            
            # Plot fitted line
            ax.plot(x_fit, y_fit, color=colors[ind], linewidth=1)
            
            print(f"q={q}, Ks={soil}: power law a={a_fit:.3f}, b={b_fit:.3f}, R²={r2:.3f}")
            
            # Add annotations
            key = f"{q}{soil}"
            ax.annotate(f"y={ann_ab[key][0]}", xy=ann_ab[key][1],
                       color=colors[ind], fontsize=15)
            ax.annotate(f"q={q}", 
                       xy=(df_q["x"].tolist()[1] + 0.1, df_q["y"].tolist()[1] - 0.002),
                       color=colors[ind], fontsize=15)
            
            ind += 1
    
    # Create legend
    for k in markers.keys():
        ax.scatter([], [], c="k", s=80, marker=markers[k], label=f"{np.around(k, 4)}")
    
    ax.legend(title="Ks (m/hr)", loc='upper right')
    
    # Format plot
    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    xlabel = r'$\alpha\cdot d (-)$'
    ylabel = r'$v_{d} / K_{s} (-)$'
    ax.xlabel(f"{xlabel}")
    ax.ylabel(f"{ylabel}")
    
    # Save figure
    ax.savefig(os.path.join(output_path, "exf_intercept_qfit.png"))
    plt.close()


def q_vs_slope():
    """
    Analyze relationship between infiltration rate and power-law slope.
    
    Plots the slope of the power-law fits (from different infiltration rates)
    against the infiltration rate itself. Fits a logarithmic model to this
    relationship.
    """
    input_path = os.path.join(DIRPATH, "inputs", "q_vs_fittinglinesslope2.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    
    # Load data
    df = pd.read_csv(input_path)
    q = df['q (m/hr)'].values
    slope = df['fitting line slope'].values
    
    print(f"q values: {q}")
    print(f"q data type: {type(q)}")
    
    # Define logarithmic model
    def log_model(x, a, b):
        """Logarithmic model: y = a * ln(x) + b"""
        return a * np.log(x) + b
    
    # Perform curve fitting
    popt, pcov = curve_fit(log_model, q, slope)
    
    # Get fitted values
    slope_fitted = log_model(q, *popt)
    
    # Calculate R-squared
    r_squared = r2_score(slope, slope_fitted)
    
    # Create plot
    plt.figure(figsize=(16, 9))
    plt.scatter(q, slope, color='k', label='Data Points')
    plt.plot(q, slope_fitted, color='k',
            label=f'y = {popt[0]:.2f}ln(Δq) {popt[1]:.2f}\n$R^2$ = {r_squared:.2f}')
    plt.xscale('log')
    plt.xlabel('Δq (m/hr)')
    plt.ylabel('Power laws slope')
    plt.legend()
    plt.grid(True)
    
    # Save figure
    plt.savefig(os.path.join(output_path, "q_vs_slope2.png"))
    plt.close()