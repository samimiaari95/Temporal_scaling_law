from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import math
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp, anderson
import os
import matplotlib
import matplotlib.pyplot as plt
from utils import powerlaw_func, linear_law

plt.rcParams.update({'font.size': 22})


def fitting_func(xlist, ylist):
    # linearize
    y_lin = np.log(ylist)
    x_lin = np.log(xlist)
    # Fit the function
    params, covariance = curve_fit(linear_law, x_lin, y_lin)
    a_fit, b_fit = params

    # fitting accuracy
    y_fit = [linear_law(x, a_fit, b_fit) for x in x_lin]
    R_square = r2_score(y_lin, y_fit)
    r = np.corrcoef(y_lin, y_fit)
    r2 = r[0][1]**2

    # back transform to power law
    a_fit = np.exp(a_fit)

    x_fit = [min(xlist), max(xlist)]
    y_fit = [powerlaw_func(x, a_fit, b_fit) for x in x_fit]
    return x_fit, y_fit, a_fit, b_fit, r2

def plot8(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = alfa*d
    y = ((d-0.05)-abs(toplayer_pressure))/(exf_t*k)

    xlabel = r'$\alpha\cdot d (-)$'
    ylabel = r'$\frac{d-|\psi_{toplayer}|}{t\cdot K_s} (-)$'
    return x, y, xlabel, ylabel

def qfit():
    exf_cases_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "inputs", "exf_intercept.csv")
    output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "outputs")

    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    inf_t_column = "inf_time"
    exf_t_column = "exf_time"
    toplayer_pressure_column = "toplayer_pressure"

    df = pd.read_csv(exf_cases_path)

    soil_types = list(df[k_column].unique())
    soil_types.sort()

    markers = {0.002:"P", 0.0045:"d"}
    markers = {0.0045:"d"}

    ax = plt
    ax.figure(figsize=(16,9))

    df["x"] = df[d_column]*df[alfa_column]
    df["y"] = ((df[d_column]-0.05)-abs(df[toplayer_pressure_column]))/(df[exf_t_column]*df[k_column])
    # df["y"] = ((df[d_column]-0.05)-abs(df[toplayer_pressure_column]))#df[exf_t_column]
    alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"]
    ann_ab = {"0.00010.002": [r'$0.225x^{1.356}$', (0.73, 0.234386417)],
            "0.000150.002": [r'$0.234x^{1.385}$', (0.73, 0.338135512)],
            "0.0010.0045": [r'$0.304x^{1.717}$', (2, 0.09250808557112719)],
            "0.00010.0045": [r'$0.125x^{1.306}$', (2, 0.05052146488857635)],
            "0.000150.0045": [r'$0.156x^{1.409}$', (2, 0.05871341875690252)],
            "0.00050.0045": [r'$0.253x^{1.634}$', (2, 0.08165277976910569)]
            }
    ind=0
    colors = ["r", "g", "b", "darkorange"]
    for soil in soil_types:
        if soil == 0.002:
            continue
        x_axis = []
        y_axis = []
        df_soil = df[df[k_column]==soil]
        ax.scatter(df_soil["x"], df_soil["y"], c="k",s=80, marker=markers[soil])

        for q in df_soil[q_column].unique():
            df_q = df_soil[df_soil[q_column]==q]
            x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(df_q["x"], df_q["y"])
            print(df_q["x"])
            print(type(df_q["x"]))
            #ax.plot(x_fit, y_fit, color="k", linewidth=1)
            ax.plot(x_fit, y_fit, color=colors[ind] ,linewidth=1)
            print(f"{q}{soil} this is {alphabet[ind]} and minx={min(x_fit)} and maxy={max(y_fit)}")
            #ax.annotate(f"R²={round(r2,2)}\nf(x)={round(a_fit, 2)}x^({round(b_fit, 2)})", xy=(min(x_fit), max(y_fit)), color="black")
            anx = min(x_fit)-0.1 if soil==0.0045 else min(x_fit)-0.03
            anx = min(x_fit) if soil==0.0045 else min(x_fit)-0.03
            any = max(y_fit)-0.05 if soil==0.002 and q==0.0001 else max(y_fit)
            #ax.annotate(ann_ab[f"{q}{soil}"][0], xy=ann_ab[f"{q}{soil}"][1], color="black", fontsize=15)
            ax.annotate(f"y={ann_ab[f'{q}{soil}'][0]}", xy=ann_ab[f"{q}{soil}"][1],color=colors[ind], fontsize=15)
            ax.annotate(f"q={q}", xy=(df_q["x"].to_list()[1]+0.1, df_q["y"].to_list()[1]-0.002),color=colors[ind], fontsize=15)
            ind += 1

    for k in markers.keys():
        ax.scatter([],[], c="k",s=80, marker=markers[k], label=f"{np.around(k, 4)}")

    ax.legend(title="Ks (m/hr)", loc='upper right')#, bbox_to_anchor=(1.135, 0.5))


    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    xlabel = r'$\alpha\cdot d (-)$'
    # xlabel = r'$d(m)$'
    ylabel = r'$\frac{λ}{t_{dr}\cdot K_s} (-)$'
    # ylabel = r'$d-|\psi_{toplayer}| (m)$'
    ax.xlabel(f"{xlabel}")
    ax.ylabel(f"{ylabel}")
    ax.savefig(os.path.join(output_path, "exf_intercept_qfit.png"))


def q_vs_slope():
    input_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "inputs", "q_vs_fittinglinesslope2.csv")
    output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "outputs")


    # Extracting the data
    df = pd.read_csv(input_path)
    q = df['q (m/hr)'].values
    slope = df['fitting line slope'].values
    print(q)
    print(type(q))
    # Define a logarithmic model for fitting
    def log_model(x, a, b):
        return a * np.log(x) + b

    # Perform curve fitting
    popt, pcov = curve_fit(log_model, q, slope)

    # Get the fitted values
    slope_fitted = log_model(q, *popt)

    # Calculate R-squared for the fit
    r_squared = r2_score(slope, slope_fitted)

    # Plot the data points and the fitting curve
    plt.figure(figsize=(16,9))
    plt.scatter(q, slope, color='k', label='Data Points')
    plt.plot(q, slope_fitted, color='k', label=f'y = {popt[0]:.2f}ln(Δq) {popt[1]:.2f}\n$R^2$ = {r_squared:.2f}')
    plt.xscale('log')  # since we're dealing with logarithmic fitting
    plt.xlabel('Δq (m/hr)')
    plt.ylabel('Power laws slope')
    #plt.title('Logarithmic Fit of Fitting Line Slope vs. q')
    plt.legend()
    plt.grid(True)

    # Show plot
    plt.savefig(os.path.join(output_path, "q_vs_slope2.png"))

qfit()