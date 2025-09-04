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


# infiltration plot
def plot11(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = alfa*d*q/k
    y = ((d-0.05)-abs(toplayer_pressure))/(inf_t)
    xlabel = r'$\frac{\alpha\cdot d\cdot q}{K_s} (-)$'
    ylabel = r'$\frac{λ}{SST_{i}} (m/hr)$'
    return x, y, xlabel, ylabel

inf_cases_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "inputs", "infiltration_pressure_index.csv")
output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "outputs")

q_column = "q"
k_column = "k"
d_column = "d"
n_column = "n"
alfa_column = "alfa"
theta_r_column = "theta_r"
theta_s_column = "theta_s"
inf_t_column = "time"
toplayer_pressure_column = "toplayer_pressure"

df = pd.read_csv(inf_cases_path)

soil_types = list(df[k_column].unique())
soil_types.sort()
soil_types_colors = plt.cm.jet(np.linspace(0,1,len(soil_types)))
colors_dic = {soil_types[i]:soil_types_colors[i] for i in range(len(soil_types))}
list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5,2), "x"]

print(len(soil_types))
markers = {soil_types[i]:list_markers[i] for i in range(len(soil_types))}

fig, ax = plt.subplots(figsize=(16, 9))
#ax.figure(figsize=(24,11))
print("plotting points...")
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
        toplayer_pressure = df_soil[toplayer_pressure_column].iloc[i]
        if k>=q:# and k!=0.0025 and k!=0.0026 and k!=0.0012 and k!=0.002:# and (q/k)>=0.001:
            x, y, xlabel, ylabel = plot11(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)
            y_axis.append(y)
            x_axis.append(x)

            ax.scatter(x, y, c=colors_dic[soil], cmap='jet',s=80, marker=markers[k])#, label=k)

    x_all.extend(x_axis)
    y_all.extend(y_axis)
    colors.extend([soil]*len(x_axis))
for k in markers.keys():
    # change logarithmic colors for every marker
    ax.scatter([],[], c=colors_dic[k],s=80, marker=markers[k], label=f"{np.around(k, 4)}")

#ax.legend(title="Ks (m/hr)", loc='center right', bbox_to_anchor=(1.135, 0.5))
ax.legend(title="Ks (m/hr)", loc='center right', bbox_to_anchor=(1.01, 0.41))
#ax.tight_layout()

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
print(f"pearson r2: {r2}")

# back transform to power law
a_fit = np.exp(a_fit)
print(f"Fitted a: {a_fit} and b:{b_fit}")

# check kolmogorov-smirnov
#anderson_darling(y_all)
#KS_fit(y_all, y_fit)

x_fit = [min(x_all), max(x_all)]
y_fit = [powerlaw_func(x, a_fit, b_fit) for x in x_fit]
plt.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=5)
#plt.annotate(f"R²={round(r2,2)}\nf(x)={round(a_fit, 2)}x^({round(b_fit, 2)})", xy=(min(x_all), max(y_all)/10), color="black")
print("annotations:")
print(f"r2:{round(r2,2)}")
print(f"f(x)={round(a_fit, 2)}x^({round(b_fit, 2)}")

#ax.scatter(x_all, y_all,c=colors, cmap='jet', s=30, norm=matplotlib.colors.LogNorm(vmin=kmin, vmax=kmax))
ax.grid(True)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(f"{xlabel}", fontsize=32)
ax.set_ylabel(f"{ylabel}", fontsize=32)
#plt.gcf().subplots_adjust(bottom=0.4)
plt.tight_layout()

fig.savefig(os.path.join(output_path, "inf_adq-k_vs_v_symbolstest.png"))
