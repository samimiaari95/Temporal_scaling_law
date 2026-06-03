import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cmocean

from analysis import fitting_func
from file_io import DIRPATH


plt.rcParams.update({"font.size": 20})

CSV_PATH = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
OUTPUT_DIR = os.path.join(DIRPATH, "outputs", "all_combinations")


def load_data(csv_path=CSV_PATH, time_col="inf_time"):
    """Load the CSV and build the dimensionless variables used for the plots."""
    df = pd.read_csv(csv_path)

    required_columns = {"q", "k", "d", "alfa", "n", "theta_s", "theta_r", time_col}
    missing = required_columns.difference(df.columns)
    if missing:
        raise ValueError(f"CSV is missing required columns: {sorted(missing)}")

    if time_col == "inf_time":
        time_label = r"t'_i"
    elif time_col == "exf_time":
        time_label = r"t'_d"
    else:
        time_label = time_col

    data = pd.DataFrame(
        {
            "q": df["q"],
            "d": df["d"],
            "q_prime": df["q"] / df["k"],
            "d_prime": df["alfa"] * df["d"],
            "n_theta_s": df["n"] * (df["theta_s"] - df["theta_r"]),
            "time_value": df[time_col],
            "k": df["k"],
            "time_label": time_label,
        }
    )

    return data.dropna()


def _latex_label(name):
    labels = {
        "q_prime": r"q'",
        "d_prime": r"d'",
        "n_theta_s": r"n'",
    }
    return labels[name]


def _slug(name):
    return (
        name.replace("q_prime", "qprime")
        .replace("d_prime", "dprime")
        .replace("n_theta_s", "nthetas")
        .replace("/", "_over_")
        .replace("*", "_times_")
        .replace(" ", "")
        .replace("(", "")
        .replace(")", "")
    )


def _panel_letter(index):
    return f"({chr(ord('a') + index)})"


def build_combinations(data, include_reciprocal=False):
    """Create the requested x-axis combinations from q', d', and ntheta_s."""
    variables = {
        "q_prime": data["q_prime"].to_numpy(),
        "d_prime": data["d_prime"].to_numpy(),
        "n_theta_s": data["n_theta_s"].to_numpy(),
    }

    combos = []

    # One variable divided by another.
    # ordered_pairs = [
    #     ("q_prime", "d_prime"),
    #     ("q_prime", "n_theta_s"),
    #     ("d_prime", "q_prime"),
    #     ("d_prime", "n_theta_s"),
    #     ("n_theta_s", "q_prime"),
    #     ("n_theta_s", "d_prime"),
    # ]
    ordered_pairs = [
        ("q_prime", "d_prime"),
        ("q_prime", "n_theta_s"),
        ("d_prime", "n_theta_s")
    ]
    for numerator, denominator in ordered_pairs:
        combos.append(
            {
                "name": f"{numerator}_over_{denominator}",
                "xlabel": rf"$\frac{{{_latex_label(numerator)}}}{{{_latex_label(denominator)}}}$",
                "x": variables[numerator] / variables[denominator],
            }
        )

    # One variable divided by the product of the remaining two.
    for numerator in variables:
        if "d" in numerator:
            continue  # skip d'
        remaining = [key for key in variables if key != numerator]
        combos.append(
            {
                "name": f"{numerator}_over_{remaining[0]}_{remaining[1]}",
                "xlabel": rf"$\frac{{{_latex_label(numerator)}}}{{{_latex_label(remaining[0])}\,{_latex_label(remaining[1])}}}$",
                "x": variables[numerator] / (variables[remaining[0]] * variables[remaining[1]]),
            }
        )

    # Multiplication of all variables.
    combos.append(
        {
            "name": "product_all",
            "xlabel": r"$q'\, d'\, n'$",
            "x": variables["q_prime"] * variables["d_prime"] * variables["n_theta_s"],
        }
    )

    # if include_reciprocal:
    #     combos.append(
    #         {
    #             "name": "reciprocal_product_all",
    #             "xlabel": r"$1/(q'\, d'\, n')$",
    #             "x": 1.0 / (variables["q_prime"] * variables["d_prime"] * variables["n_theta_s"]),
    #         }
    #     )

    return combos


def load_velocity_data(csv_path=CSV_PATH, time_col="exf_time"):
    """Load data for the additional velocity-based group: v_d/K_s vs combinations of d' and n'."""
    df = pd.read_csv(csv_path)

    required_columns = {"q", "k", "d", "alfa", "n", "theta_s", "theta_r", "toplayer_pressure", time_col}
    missing = required_columns.difference(df.columns)
    if missing:
        raise ValueError(f"CSV is missing required columns: {sorted(missing)}")

    if time_col == "exf_time":
        time_label = r"t'_d"
    else:
        time_label = time_col

    hydraulic_head = df["d"] - np.abs(df["toplayer_pressure"])
    v_d = hydraulic_head / df[time_col]

    data = pd.DataFrame(
        {
            "d_prime": df["alfa"] * df["d"],
            "n_theta_s": df["n"] * (df["theta_s"] - df["theta_r"]),
            "y": v_d / df["k"],
            "k": df["k"],
            "time_label": time_label,
        }
    )

    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    return data


def build_velocity_combinations(data, include_reciprocal=False):
    """Create combinations for variables d' and n'."""
    variables = {
        "d_prime": data["d_prime"].to_numpy(),
        "n_theta_s": data["n_theta_s"].to_numpy(),
    }

    combos = []

    ordered_pairs = [
        ("d_prime", "n_theta_s")
    ]
    for numerator, denominator in ordered_pairs:
        combos.append(
            {
                "name": f"{numerator}_over_{denominator}",
                "xlabel": rf"$\frac{{{_latex_label(numerator)}}}{{{_latex_label(denominator)}}}$",
                "x": variables[numerator] / variables[denominator],
            }
        )

    combos.append(
        {
            "name": "d_prime_only",
            "xlabel": r"$d'$",
            "x": variables["d_prime"],
        }
    )
    combos.append(
        {
            "name": "product_all",
            "xlabel": r"$d'\, n'$",
            "x": variables["d_prime"] * variables["n_theta_s"],
        }
    )

    return combos


def _plot_single(ax, x_values, y_values, k_values, xlabel, y_label, colors_dic, markers, panel_label=None):
    mask = (x_values > 0) & (y_values > 0)
    x_plot = x_values[mask]
    y_plot = y_values[mask]
    k_plot = k_values[mask]

    if x_plot.size < 3:
        ax.axis("off")
        return None

    x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(x_plot, y_plot)

    for soil in sorted(np.unique(k_plot)):
        soil_mask = k_plot == soil
        ax.scatter(
            x_plot[soil_mask],
            y_plot[soil_mask],
            c=[colors_dic[soil]],
            s=80,
            marker=markers[soil],
            alpha=0.7,
        )

    ax.plot(x_fit, y_fit, color="k", linewidth=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(y_label)
    ax.set_title(rf"$y = {a_fit:.2f}x^{{{b_fit:.2f}}},\ R^2 = {r2:.2f}$")
    if panel_label is not None:
        ax.text(
            0.03,
            0.94,
            panel_label,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=18,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8),
        )
    ax.grid(True, which="both", alpha=0.3)
    return {
        "a_fit": a_fit,
        "b_fit": b_fit,
        "r2": r2,
    }


def plot_all_combinations(
    csv_path=CSV_PATH,
    output_dir=OUTPUT_DIR,
    save_individual=True,
    time_col="inf_time",
    y_mode="ks",
):
    data = load_data(csv_path, time_col=time_col)
    time_label = data["time_label"].iloc[0]
    if y_mode == "ks":
        y_values = data["k"].to_numpy() * data["time_value"].to_numpy() / data["d"].to_numpy()
        y_label = rf"$K_s {time_label}$"
        output_tag = "ks"
        include_reciprocal = False
    elif y_mode == "q":
        y_values = data["q"].to_numpy() * data["time_value"].to_numpy() / data["d"].to_numpy()
        y_label = rf"$q {time_label}$"
        output_tag = "q"
        include_reciprocal = True
    else:
        raise ValueError("y_mode must be 'ks' or 'q'")
    soil_types = sorted(data["k"].unique())
    soil_type_colors = cmocean.cm.balance(np.linspace(0, 1, len(soil_types)))
    colors_dic = {soil_types[i]: soil_type_colors[i] for i in range(len(soil_types))}

    list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5, 2), "x"]
    markers = {soil_types[i]: list_markers[i] for i in range(len(soil_types))}

    combos = build_combinations(data, include_reciprocal=include_reciprocal)

    time_output_dir = os.path.join(output_dir, f"{time_col}_{output_tag}")
    os.makedirs(time_output_dir, exist_ok=True)

    n_plots = len(combos)
    ncols = 3
    nrows = int(np.ceil(n_plots / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(18, 5.5 * nrows))
    axes = np.atleast_1d(axes).ravel()

    results = []
    for idx, combo in enumerate(combos):
        ax = axes[idx]
        panel_label = _panel_letter(idx)
        result = _plot_single(
            ax,
            combo["x"],
            y_values,
            data["k"].to_numpy(),
            combo["xlabel"],
            y_label,
            colors_dic,
            markers,
            panel_label=panel_label,
        )
        results.append((combo["name"], result))

        if save_individual and result is not None:
            fig_ind, ax_ind = plt.subplots(figsize=(8, 6))
            _plot_single(
                ax_ind,
                combo["x"],
                y_values,
                data["k"].to_numpy(),
                combo["xlabel"],
                y_label,
                colors_dic,
                markers,
                panel_label=panel_label,
            )
            fig_ind.tight_layout()
            fig_ind.savefig(os.path.join(time_output_dir, f"{_slug(combo['name'])}.png"), dpi=300)
            plt.close(fig_ind)

    for idx in range(n_plots, len(axes)):
        axes[idx].axis("off")

    legend_handles = []
    for k in markers.keys():
        handle = axes[0].scatter([], [], c=[colors_dic[k]], s=80, marker=markers[k], label=f"{np.around(k, 4)}")
        legend_handles.append(handle)

    fig.legend(
        handles=legend_handles,
        title=r"$K_s$ (m/hr)",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.10),
        ncol=4,
    )
    # fig.suptitle(rf"All combinations of $q'$, $d'$, and $n'$ against {y_label}", y=0.995)
    fig.tight_layout(rect=[0, 0.06, 1, 0.98])
    fig.savefig(os.path.join(time_output_dir, f"all_combinations_grid_{time_col}_{output_tag}.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    return results


def plot_velocity_group(csv_path=CSV_PATH, output_dir=OUTPUT_DIR, save_individual=True, time_col="exf_time"):
    """Additional group: y = v_d/K_s vs combinations of d' and n'."""
    data = load_velocity_data(csv_path=csv_path, time_col=time_col)
    y_values = data["y"].to_numpy()
    y_label = r"$v_d/K_s$"

    soil_types = sorted(data["k"].unique())
    soil_type_colors = cmocean.cm.balance(np.linspace(0, 1, len(soil_types)))
    colors_dic = {soil_types[i]: soil_type_colors[i] for i in range(len(soil_types))}

    list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5, 2), "x"]
    markers = {soil_types[i]: list_markers[i] for i in range(len(soil_types))}

    combos = build_velocity_combinations(data, include_reciprocal=True)

    velocity_output_dir = os.path.join(output_dir, f"{time_col}_vd_ks")
    os.makedirs(velocity_output_dir, exist_ok=True)

    n_plots = len(combos)
    ncols = 3
    nrows = int(np.ceil(n_plots / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(18, 5.5 * nrows))
    axes = np.atleast_1d(axes).ravel()

    results = []
    for idx, combo in enumerate(combos):
        ax = axes[idx]
        panel_label = _panel_letter(idx)
        result = _plot_single(
            ax,
            combo["x"],
            y_values,
            data["k"].to_numpy(),
            combo["xlabel"],
            y_label,
            colors_dic,
            markers,
            panel_label=panel_label,
        )
        results.append((combo["name"], result))

        if save_individual and result is not None:
            fig_ind, ax_ind = plt.subplots(figsize=(8, 6))
            _plot_single(
                ax_ind,
                combo["x"],
                y_values,
                data["k"].to_numpy(),
                combo["xlabel"],
                y_label,
                colors_dic,
                markers,
                panel_label=panel_label,
            )
            fig_ind.tight_layout()
            fig_ind.savefig(os.path.join(velocity_output_dir, f"{_slug(combo['name'])}.png"), dpi=300)
            plt.close(fig_ind)

    for idx in range(n_plots, len(axes)):
        axes[idx].axis("off")

    legend_handles = []
    for k in markers.keys():
        handle = axes[0].scatter([], [], c=[colors_dic[k]], s=80, marker=markers[k], label=f"{np.around(k, 4)}")
        legend_handles.append(handle)

    fig.legend(
        handles=legend_handles,
        title=r"$K_s$ (m/hr)",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.2),
        ncol=4,
    )
    # fig.suptitle(r"All combinations of $d'$ and $n'$ against $v_d/K_s$", y=0.995)
    fig.tight_layout(rect=[0, 0.06, 1, 0.98])
    fig.savefig(os.path.join(velocity_output_dir, f"all_combinations_grid_{time_col}_vd_ks.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    return results


def main():
    for time_col in ("inf_time", "exf_time"):
        for y_mode in ("ks", "q"):
            results = plot_all_combinations(time_col=time_col, y_mode=y_mode)
            for name, result in results:
                if result is None:
                    print(f"{time_col}::{y_mode}::{name}: insufficient data")
                else:
                    print(
                        f"{time_col}::{y_mode}::{name}: y = {result['a_fit']:.3e} x^{result['b_fit']:.3f}, R^2 = {result['r2']:.3f}"
                    )

    velocity_results = plot_velocity_group(time_col="exf_time")
    for name, result in velocity_results:
        if result is None:
            print(f"exf_time::vd_ks::{name}: insufficient data")
        else:
            print(
                f"exf_time::vd_ks::{name}: y = {result['a_fit']:.3e} x^{result['b_fit']:.3f}, R^2 = {result['r2']:.3f}"
            )


if __name__ == "__main__":
    main()