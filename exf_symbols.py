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


def plot2(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = alfa*d
    y = (exf_t*k*alfa)/(theta_s-theta_r)
    xlabel = r'$\alpha\cdot d (-)$'
    ylabel = r'$\frac{t_{dr}\cdot K_s\cdot \alpha}{\theta_s-\theta_r} (-)$'
    return x, y, xlabel, ylabel

def plot8(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = alfa*d
    y = ((d-0.05)-abs(toplayer_pressure))/(exf_t*k)
    xlabel = r'$\alpha\cdot d (-)$'
    ylabel = r'$\frac{λ}{t_{dr}\cdot K_s} (-)$'
    return x, y, xlabel, ylabel


exf_cases_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "inputs", "drainage_inf_testcases_toplayer.csv")
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

list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5,2), "x"]
print(len(list_markers))
print(len(soil_types))
print(soil_types)
markers = {soil_types[i]:list_markers[i] for i in range(len(soil_types))}
print(markers)

ax = plt
ax.figure(figsize=(22,9))

x_all = []
y_all = []
colors = []
for soil in soil_types:
    x_axis = []
    y_axis = []
    df_soil = df[df[k_column]==soil]
    for i in range(len(df_soil)):
        q = df_soil[q_column].iloc[i]
        k = df_soil[k_column].iloc[i]
        d = df_soil[d_column].iloc[i]
        n = df_soil[n_column].iloc[i]
        alfa = df_soil[alfa_column].iloc[i]
        theta_r = df_soil[theta_r_column].iloc[i]
        theta_s = df_soil[theta_s_column].iloc[i]
        inf_t = df_soil[inf_t_column].iloc[i]
        exf_t = df_soil[exf_t_column].iloc[i]
        toplayer_pressure = df_soil[toplayer_pressure_column].iloc[i]
        if k>=q and exf_t!=0:# and k!=0.0025 and k!=0.0026 and k!=0.0012 and k!=0.002:
            x, y, xlabel, ylabel = plot2(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure)
            y_axis.append(y)
            x_axis.append(x)
            
            ax.scatter(x, y, c="k",s=80, marker=markers[k])
            # Add text labels for each point
            #labels = f"q={q:.2e}"
            #plt.text(x, y, labels, fontsize=9, ha='right', va='bottom')  # Adjust alignment as needed


    x_all.extend(x_axis)
    y_all.extend(y_axis)
    colors.extend([soil]*len(x_axis))
for k in markers.keys():
    ax.scatter([],[], c="k",s=80, marker=markers[k], label=f"{np.around(k, 4)}")

ax.legend(title="Ks (m/hr)", loc='upper right', bbox_to_anchor=(1.01, 0.8))
# ax.legend(title="Ks (m/hr)", loc='upper right', bbox_to_anchor=(1.008, 1.018))

# linearize
y_lin = np.log(y_all)
x_lin = np.log(x_all)
# Fit the function
params, covariance = curve_fit(linear_law, x_lin, y_lin)
a_fit, b_fit = params

# fitting accuracy
y_fit = [linear_law(x, a_fit, b_fit) for x in x_lin]
R_square = r2_score(y_lin, y_fit)
print(f"R2 = {R_square}")
r = np.corrcoef(y_lin, y_fit)
r2 = r[0][1]**2

# back transform to power law
a_fit = np.exp(a_fit)
print(f"Fitted a: {a_fit} and b:{b_fit}")

x_fit = [min(x_all), max(x_all)]
y_fit = [powerlaw_func(x, a_fit, b_fit) for x in x_fit]
plt.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=5)
#plt.annotate(f"R²={round(r2,2)}\nf(x)={round(a_fit, 2)}x^({round(b_fit, 2)})", xy=(min(x_all), max(y_all)/10), color="black")
print("annotations:")
print(f"r2:{round(r2,2)}")
print(f"f(x)={round(a_fit, 2)}x^({round(b_fit, 2)}")
ax.grid(True)
ax.xscale("log")
ax.yscale("log")
ax.xlabel(f"{xlabel}", fontsize=32)
ax.ylabel(f"{ylabel}", fontsize=32)
ax.savefig(os.path.join(output_path, f"exf_ad_vs_tka-thetasr_symbolstesting.png"))
