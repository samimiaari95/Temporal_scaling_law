from file_io import DIRPATH, SCRATCHPATH
from analysis import fitting_func
import os
import cmocean
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.legend_handler import HandlerTuple
import SLOTH.sloth.IO
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# ==============================================================================
# CATEGORY 4: PLOTTING AND VISUALIZATION
# Functions for creating plots and figures
# ==============================================================================


def plot_pressure_profile():
    """
    Plot temporal evolution of pressure profile during infiltration.
    
    Creates a figure showing how the pressure profile develops over time
    from initial condition to steady state. Useful for visualizing the
    infiltration front propagation.
    """
    # plt.rcParams.update({'font.size': 16})
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
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(DIRPATH, "pressure_profile.png")
    fig.savefig(output_path, dpi=300)
    plt.close()

def plot_failedprofile():
    name = os.path.join(SCRATCHPATH, "infexfcases", "test_case351", "infiltration")
    depth = [-1 * z / 10 for z in range(0, 70, 1)]
    for t in range(0, 150):
        pfb_file = name + '.out.press.' + ('{:05d}'.format(t)) + '.pfb'
        data = SLOTH.sloth.IO.read_pfb(pfb_file)
        plt.plot(data[:, 0, 0], list(reversed(depth)), color="black")
    print(data)
    plt.savefig(os.path.join(DIRPATH, "failed_profile.png"))
    plt.close()
        

def infexf_lambda_dependence_backup():
    """
    Plot drainage/infiltration time ratio vs. characteristic length scale.
    
    Analyzes the relationship between the drainage-to-infiltration time ratio
    and the characteristic length scale (lambda = d - |psi_top|).
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Column names
    d_column = "d"
    inf_t_column = "inf_time"
    exf_t_column = "exf_time"
    toplayer_pressure_column = "toplayer_pressure"
    
    # Create figure
    plt.figure(figsize=(16, 9))
    colors = cmocean.cm.balance(np.linspace(0, 1, len(df.groupby(d_column))))

    # Plot each depth group with different color
    for d_val, group in df.groupby(d_column):
        y = group[exf_t_column] / group[inf_t_column]
        x = group[d_column] - 0.05 - abs(group[toplayer_pressure_column])
        plt.scatter(x, y, label=f"d = {d_val}", color=colors[list(df.groupby(d_column).groups.keys()).index(d_val)])
    
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

def infexf_lambda_dependence():
    """
    Plot drainage/infiltration time ratio vs. characteristic length scale.
    
    Analyzes the relationship between the drainage-to-infiltration time ratio
    and the characteristic length scale (lambda = d - |psi_top|).
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Column names
    d_column = "d"
    inf_t_column = "inf_time"
    exf_t_column = "exf_time"
    toplayer_pressure_column = "toplayer_pressure"
    
    # Create figure
    fig, ax = plt.subplots(figsize=(7.09, 4))

    depth_markers = {
        1: "o",
        2: "^",
        3: "s",
        4: "P",
        5: "*",
        6: "X",
        7: "D",
        8: "v",
        9: "<",
        10: ">",
    }

    ks_values = df["k"].to_numpy()
    positive_ks = ks_values[ks_values > 0]
    norm = LogNorm(vmin=np.nanmin(positive_ks), vmax=np.nanmax(positive_ks))
    cmap = cmocean.cm.balance
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    # Plot each depth group with a distinct marker and Ks-based color scale
    for d_val in sorted(df[d_column].unique()):
        group = df[df[d_column] == d_val]
        y = group[exf_t_column] / group[inf_t_column]
        x = group[d_column] - 0.05 - abs(group[toplayer_pressure_column])
        marker = depth_markers.get(int(d_val), "o")
        ax.scatter(
            x,
            y,
            label=f"d = {d_val}",
            c=group["k"],
            cmap=cmap,
            norm=norm,
            marker=marker,
            edgecolors="none",
            s=10,
        )
    
    # Format plot
    ax.grid(True)
    # plt.xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("λ (m)")
    ax.set_ylabel(r"$t_{d} / t_{i}  (-)$")
    depth_handles = [
        ax.scatter([], [], color="black", marker=marker, s=10, label=f"d = {depth}")
        for depth, marker in depth_markers.items()
    ]
    ax.legend(handles=depth_handles, ncols=5, title="d values", fontsize='small', labelspacing=0.3, handletextpad=0.5, framealpha=0.6)
    fig.colorbar(sm, ax=ax, label=r"$K_s$ (m/hr)", fraction=0.02, pad=0.03)
    fig.tight_layout()
    
    # Save figure
    fig.savefig(os.path.join(output_path, "exfinf_lambda.pdf"), dpi=500)
    fig.savefig(os.path.join(output_path, "exfinf_lambda.png"), dpi=500)
    plt.close(fig)

def depth_dependence_combined():
    """
    Plot steady-state time vs. soil depth for:
    (a) Infiltration
    (b) Drainage
    
    Produces a single 180 mm wide publication-ready figure.
    """
    # ---- Global font settings ----
    # plt.rcParams.update({
    #     "font.size": 16,
    #     "axes.labelsize": 16,
    #     "xtick.labelsize": 14,
    #     "ytick.labelsize": 14,
    #     "legend.fontsize": 14
    # })
    exf_cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    
    df = pd.read_csv(exf_cases_path)
    
    # Extract data
    x = df["d"]
    y_drain = df["exf_time"]
    y_inf = df["inf_time"]
    
    # --- Figure size ---
    width_mm = 180
    width_in = width_mm / 25.4
    height_in = width_in * 0.5   # good aspect ratio for 2 panels
    
    fig, axes = plt.subplots(1, 2, figsize=(width_in, height_in), sharey=False)
    
    # ---- (a) Infiltration ----
    axes[0].scatter(x, y_inf, color='black', s=10)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("d (m)")
    axes[0].set_ylabel(r"$t_{i}$ (hr)")
    axes[0].grid(True, lw=0.5)
    axes[0].text(0.85, 0.12, "(a)", transform=axes[0].transAxes,
                 fontsize=12, fontweight="bold", va="top")

    # ---- (b) Drainage ----
    axes[1].scatter(x, y_drain, color='black', s=10)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("d (m)")
    axes[1].set_ylabel(r"$t_{d}$ (hr)")
    axes[1].grid(True, lw=0.5)
    axes[1].text(0.85, 0.12, "(b)", transform=axes[1].transAxes,
                 fontsize=12, fontweight="bold", va="top")
    
    
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_path, "td_ti_vs_d_combined.pdf"),
                dpi=500,
                bbox_inches="tight")
    plt.savefig(os.path.join(output_path, "td_ti_vs_d_combined.png"),
                dpi=500,
                bbox_inches="tight")
    
    plt.close()


def compare_inf_exf_times_between_files():
    """
    Compare infiltration and exfiltration times from two input tables.

    Produces a two-panel figure with infiltration time on the left and
    exfiltration time on the right, overlaying the values from
    tolerance_05_inf_exf_times.csv and inf_exf_times_config.csv.
    The x-axis is water table depth from the CSV `d` column, sorted from shallow to deep.
    """
    tolerance_cases_path = os.path.join(DIRPATH, "inputs", "tolerance_05_inf_exf_times.csv")
    config_cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")

    tolerance_df = pd.read_csv(tolerance_cases_path)
    config_df = pd.read_csv(config_cases_path)

    if "case_index" in tolerance_df.columns and "case_index" in config_df.columns:
        tolerance_df = tolerance_df.set_index("case_index")
        config_df = config_df.set_index("case_index")

    tolerance_df = tolerance_df.assign(water_table_depth=tolerance_df["d"])
    config_df = config_df.assign(water_table_depth=config_df["d"])

    tolerance_inf_df = tolerance_df.sort_values("water_table_depth")
    config_inf_df = config_df.sort_values("water_table_depth")
    tolerance_exf_df = tolerance_df.sort_values("water_table_depth")
    config_exf_df = config_df.sort_values("water_table_depth")
    tolerance_inf_df = tolerance_inf_df[tolerance_inf_df["inf_time"] > 1]
    config_inf_df = config_inf_df[config_inf_df["inf_time"] > 1]
    tolerance_exf_df = tolerance_exf_df[tolerance_exf_df["exf_time"] > 1]
    config_exf_df = config_exf_df[config_exf_df["exf_time"] > 1]

    x_inf_config = config_inf_df["water_table_depth"].to_numpy()
    x_exf_config = config_exf_df["water_table_depth"].to_numpy()
    x_inf_tolerance = tolerance_inf_df["water_table_depth"].to_numpy()
    x_exf_tolerance = tolerance_exf_df["water_table_depth"].to_numpy()
    print("plotting it")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)

    axes[0].scatter(x_inf_tolerance, tolerance_inf_df["inf_time"], label=r"ResidualTol=$10^{-5}$", color="tab:blue", marker="o", s=18, alpha=1)
    axes[0].scatter(x_inf_config, config_inf_df["inf_time"], label=r"ResidualTol=$10^{-7}$", color="tab:orange", marker="o", s=18, alpha=0.5)

    axes[1].scatter(x_exf_tolerance, tolerance_exf_df["exf_time"], label=r"ResidualTol=$10^{-5}$", color="tab:blue", marker="o", s=18, alpha=1)
    axes[1].scatter(x_exf_config, config_exf_df["exf_time"], label=r"ResidualTol=$10^{-7}$", color="tab:orange", marker="o", s=18, alpha=0.5)

    axes[0].set_title("Infiltration time")
    axes[1].set_title("Drainage time")

    for ax, ylabel in zip(axes, [r"$t_{i}$ (hr)", r"$t_{d}$ (hr)"]):
        ax.set_yscale("log")
        ax.set_xlabel("water table depth (m)")
        ax.set_ylabel(ylabel)
        ax.grid(True, lw=0.5)
        ax.legend(framealpha=0.7)

    fig.tight_layout()
    # fig.savefig(os.path.join(output_path, "compare_inf_exf_times_between_files.pdf"), dpi=500, bbox_inches="tight")
    fig.savefig(os.path.join(output_path, "compare_inf_exf_times_between_files.png"), dpi=500, bbox_inches="tight")
    plt.close(fig)


def exf_solution_plots(solution=None):
    """
    Create comprehensive plot showing drainage time scaling relationships.
    
    Plots the non-dimensional drainage velocity (lambda / (SST_d * Ks))
    against the non-dimensional depth (alpha * d) for different soil types.
    Includes power-law fit to the data.
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Get unique soil types (by hydraulic conductivity)
    soil_types = sorted(df["k"].unique())
    
    # Create color map and marker styles
    soil_types_colors = cmocean.cm.balance(np.linspace(0, 1, len(soil_types)))
    colors_dic = {soil_types[i]: soil_types_colors[i] for i in range(len(soil_types))}
    
    list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5,2), "x"]
    markers = {soil_types[i]: list_markers[i] for i in range(len(soil_types))}
    
    print(f"Number of soil types: {len(soil_types)}")
    print(f"Soil types (Ks values): {soil_types}")
    
    # Create figure
    ax = plt
    ax.figure(figsize=(7.09, 3.65))  # 180 mm wide figure with good aspect ratio
    
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
                elif solution == 3:
                    x, y, xlabel, ylabel = exfsolution_3(
                        q, k, d, df_soil["n"].iloc[i], alfa,
                        df_soil["theta_r"].iloc[i], df_soil["theta_s"].iloc[i],
                        df_soil["inf_time"].iloc[i], exf_t, toplayer_pressure
                    )
                else:
                    raise ValueError("Drainage solution must be identified as 1, 2, or 3")

                x_all.append(x)
                y_all.append(y)
                
                ax.scatter(x, y, c=colors_dic[soil], cmap='jet', s=20, marker=markers[k])
    
    # Create legend for soil types
    for k in markers.keys():
        ax.scatter([], [], c=colors_dic[k], s=80, marker=markers[k], 
                  label=f"{np.around(k, 4)}")
    
    ax.legend(title="Ks (m/hr)", loc='upper right', bbox_to_anchor=(1.01, 1.01), fontsize=12)
    
    # Perform power-law fit
    x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(x_all, y_all)
    
    #plt.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=5)
    
    print(f"Power-law fit: y = {round(a_fit, 2)} * x^{round(b_fit, 2)}")
    print(f"R² = {round(r2, 2)}")
    
    # Format plot
    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    ax.xlabel(f"{xlabel}")#, fontsize=32)
    ax.ylabel(f"{ylabel}")#, fontsize=32)
    plt.tight_layout()
    
    # Save figure
    figname = "exf_ad_vs_v-k_symbols.png" if solution == 2 else "exf_ad_vs_tka-thetasr_symbols.png"
    if solution == 3:
        figname = "exf_ad_vs_1td-symbols.pdf"
    ax.savefig(os.path.join(output_path, figname), dpi=400)
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

def exfsolution_3(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    """
    Calculate non-dimensional variables for drainage time scaling.
    
    Computes:
    - x: alpha * d (non-dimensional depth)
    - y: 1 / SST_d
    
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
    y = 1 / exf_t
    xlabel = r'$\alpha\cdot d (-)$'
    ylabel = r'$1 / t_{d} (-)$'
    return x, y, xlabel, ylabel


def inf_solution_plots():
    """
    Create comprehensive plot showing infiltration time scaling relationships.
    
    Plots the non-dimensional infiltration velocity (lambda / SST_i)
    against the non-dimensional infiltration parameter (alpha * d * q / Ks).
    Includes power-law fit to the data.
    """
    inf_cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(inf_cases_path)
    
    # Get unique soil types
    soil_types = sorted(df["k"].unique())
    
    # Create color map and markers
    soil_types_colors = cmocean.cm.balance(np.linspace(0, 1, len(soil_types)))
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
            inf_t = df_soil["inf_time"].iloc[i]
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
    fig.savefig(os.path.join(output_path, "inf_adq-k_vs_v_symbols.png"), dpi=300)
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
    ax.figure(figsize=(7.09, 3.65))  # 180 mm wide figure with good aspect ratio
    
    # Compute non-dimensional variables
    df["x"] = df["d"] * df["alfa"]
    df["y"] = ((df["d"] - 0.05) - abs(df["toplayer_pressure"])) / (df["exf_time"] * df["k"])
    
    # Annotation dictionary for fitted curves
    ann_ab = {
        "0.00010.002": [r'$0.225x^{1.356}$', (0.73, 0.234386417)],
        "0.000150.002": [r'$0.234x^{1.385}$', (0.73, 0.338135512)],
        "0.0010.0045": [r'$0.304x^{1.717}$', (2.05, 0.088)],
        "0.00010.0045": [r'$0.125x^{1.306}$', (2.05, 0.048)],
        "0.000150.0045": [r'$0.156x^{1.409}$', (2.05, 0.058)],
        "0.00050.0045": [r'$0.253x^{1.634}$', (2.05, 0.074)]
    }
    
    ind = 0
    # FIX: Replaced beige/low-contrast colormap with 4 highly distinct colors
    colors = ["#003f5c", "#58508d", "#bc5090", "#45987f"]
    
    # Plot each soil type
    for soil in soil_types:
        if soil == 0.002:
            continue
        
        df_soil = df[df["k"] == soil]
        ax.scatter(df_soil["x"], df_soil["y"], c="k", s=10, marker=markers[soil])
        
        # Fit and plot each infiltration rate separately
        for q in df_soil["q"].unique():
            df_q = df_soil[df_soil["q"] == q]
            
            # Fit power law
            x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(
                df_q["x"].tolist(), df_q["y"].tolist()
            )
            
            # Plot fitted line
            ax.plot(x_fit, y_fit, color=colors[ind], linewidth=2)
            
            print(f"q={q}, Ks={soil}: power law a={a_fit:.3f}, b={b_fit:.3f}, R²={r2:.3f}")
            
            # FIX: Merged the q label into the equation block at the roomy left-side (x=2)
            key = f"{q}{soil}"
            combined_text = f"y={ann_ab[key][0]}  $q={q}$"
            ax.annotate(combined_text, xy=ann_ab[key][1],
                       color=colors[ind], fontsize=11, va='center')
            
            ind += 1
    
    # Create legend
    for k in markers.keys():
        ax.scatter([], [], c="k", s=10, marker=markers[k], label=f"{np.around(k, 4)}")
    
    ax.legend(title="$K_s$ (m/hr)", loc='upper right')
    
    # Format plot
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.xscale("log")
    ax.yscale("log")
    xlabel = r"$d' (-)$"
    ylabel = r'$\frac{v_{d}}{K_{s}} (-)$'
    ax.xlabel(f"{xlabel}")
    ax.ylabel(f"{ylabel}")
    plt.tight_layout(pad=0.2)
    
    # Save figure
    ax.savefig(os.path.join(output_path, "exf_intercept_qfit.pdf"), dpi=500)
    ax.savefig(os.path.join(output_path, "exf_intercept_qfit.png"), dpi=500)
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
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_path, "q_vs_slope2.png"), dpi=300)
    plt.close()