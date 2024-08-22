from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn import preprocessing
import math
import pandas as pd
import numpy as np
from scipy.stats import powerlaw, kstest, ks_2samp, anderson
import os
import matplotlib
import matplotlib.pyplot as plt
from utils import powerlaw_func, linear_law

plt.rcParams.update({'font.size': 22})

def anderson_darling(data):
    result = anderson(data, dist='power law')
    print(f"Anderson-Darling Test Statistic: {result.statistic}")
    print(f"Critical Values: {result.critical_values}")
    print(f"Significance Levels: {result.significance_level}")
    print(f"Is the data drawn from the specified distribution? {result.statistic < result.critical_values[2]}")

def KS_fit(data, fitted):
    # Perform Kolmogorov-Smirnov test for goodness of fit
    #ks_statistic, ks_p_value = kstest(data, powerlaw_func, args=(alpha,b))
    ks_statistic, ks_p_value = ks_2samp(data, fitted)
    print(f"ks_statistic = {ks_statistic}")
    print(f"ks_p_value = {ks_p_value}")
    return

######################### INFILTRATION INVENTORY ###########################
def plot1(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = q
    y = t/d
    xlabel = "q (m/hr)"
    ylabel = "t/d (hr/m)"
    return x, y, xlabel, ylabel

#def plot2(q, k, d, n, alfa, theta_r, theta_s, t):
#    x = q/k
#    y = t*k/d
#    xlabel = "q/k (-)"
#    ylabel = "t*k/d (-)"
#    return x, y, xlabel, ylabel

#def plot3(q, k, d, n, alfa, theta_r, theta_s, t):
#    x = q/k
#    y = t*q/d
#    xlabel = "q/k (-)"
#    ylabel = "t*q/d (-)"
#    return x, y, xlabel, ylabel

def plot2(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = q/(k*d*alfa)
    y = t*k*alfa/((theta_s-theta_r))
    xlabel = "q/(k*d*α) (-)"
    ylabel = "t*k*α/(θs - θr) (-)"
    return x, y, xlabel, ylabel

def plot3(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = (q/(k))
    y = (t*k/(n*d))**(1-1/n)
    xlabel = "q/k (-)"
    ylabel = "(t*k/(n*d))^(1-1/n) (-)"
    return x, y, xlabel, ylabel

def plot4(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = alfa*d
    y = t*q*alfa
    xlabel = "d*α (-)"
    ylabel = "t*q*α (-)"
    return x, y, xlabel, ylabel

#def plot4star(q, k, d, n, alfa, theta_r, theta_s, t):      # better fit and physical meaning of t*q*alfa
#    x = alfa*d
#    y = t*k*alfa
#    xlabel = "d*α (-)"
#    ylabel = "t*k*α (-)"
#    return x, y, xlabel, ylabel

def plot5(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = (q)/(alfa*d)
    y = t
    xlabel = "q/(d*alpha) (m/hr)"
    ylabel = "t"
    return x, y, xlabel, ylabel

def plot6(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = q/(k*d*alfa)
    y = t/max_time
    xlabel = "q/(k*d*alpha) (-)"
    ylabel = "t/tmax (-)"
    return x, y, xlabel, ylabel

def plot7(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = (q)/(k*alfa*d)
    y = t*alfa*k/(theta_s-theta_r)
    xlabel = "(q)/(k*alfa*d)"
    ylabel = "t*alfa*k/(theta_s-theta_r)"
    return x, y, xlabel, ylabel

def plot8(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = k*(((1)/(q))**(1-1/n)) 
    y = t/max_time
    xlabel = "q/(k*d*alpha) (-)"
    ylabel = "t/tmax (-)"
    return x, y, xlabel, ylabel

def plot9(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = q/k
    y = ((d-0.05)-abs(toplayer_pressure))/t # 0.05 to account for the centroil of the layer depth, while for pressure it's already accounted for in parflow
                                            # layer thickness is 0.1m in all samples
    xlabel = "q/k (-)"
    ylabel = "(d-abs(ptop))/t (m/hr)"
    return x, y, xlabel, ylabel

def plot10(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = q/k
    y = ((d-0.05)-abs(toplayer_pressure))/(t*k)
    xlabel = "q/k (-)"
    ylabel = "(d-abs(ptop))/(t*k) (-)"
    return x, y, xlabel, ylabel

def plot11(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = alfa*d*q/k
    y = ((d-0.05)-abs(toplayer_pressure))/(t)
    xlabel = "alfa*d*q/Ks (-)"
    ylabel = "(d-abs(ptop))/(t) (m/hr)"
    return x, y, xlabel, ylabel

def plot12(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure):
    x = alfa*d#*q
    y = ((d-0.05)-abs(toplayer_pressure))/(t*k)
    xlabel = "alfa*d*q/Ks (-)"
    ylabel = "(d-abs(ptop))/(t*k) (-)"
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
t_column = "time"
index_column = "case_index"
toplayer_pressure_column = "toplayer_pressure"

df = pd.read_csv(inf_cases_path)

soil_types = list(df[k_column].unique())
soil_types.sort()
max_time = max(df[t_column])
#print(max_time)
min_time = min(df[t_column])
#print(min_time)
kmin = 10
kmax = 0

x_all = []
y_all = []
colors = []
for soil in soil_types:
    x_axis = []
    y_axis = []
    df_soil = df[df[k_column]==soil]
    index = 0

    for i in range(len(df_soil)):
        q = df_soil[q_column].iloc[i]
        k = df_soil[k_column].iloc[i]
        d = df_soil[d_column].iloc[i]
        n = df_soil[n_column].iloc[i]
        alfa = df_soil[alfa_column].iloc[i]
        theta_r = df_soil[theta_r_column].iloc[i]
        theta_s = df_soil[theta_s_column].iloc[i]
        t = df_soil[t_column].iloc[i]
        toplayer_pressure = df_soil[toplayer_pressure_column].iloc[i]
        if k>=q and t!=0 and (q/k)>=0.001 and k!=0.0025 and k!=0.0026 and k!=0.002 and k!=0.0012:# and k==0.0007:
            x, y, xlabel, ylabel = plot11(q, k, d, n, alfa, theta_r, theta_s, t, toplayer_pressure)
            y_axis.append(y)
            x_axis.append(x)
            if k>kmax:
                kmax = k
            if k<kmin:
                kmin = k
            index += 1
    x_all.extend(x_axis)
    y_all.extend(y_axis)
    colors.extend([soil]*len(x_axis))

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

ax = plt
ax.figure(figsize=(16,9))
ax.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=5)
ax.annotate(f"R²={round(R_square,2)}\nf(x)={round(a_fit, 2)}x^({round(b_fit, 2)})", xy=(max(x_all)/100, min(y_all)), color="black")
ax.scatter(x_all, y_all,c=colors, cmap='jet', s=30, norm=matplotlib.colors.LogNorm(vmin=0.0001, vmax=1)) #vmin=kmin, vmax=kmax
ax.grid(True)
ax.xscale("log")
ax.yscale("log")
ax.xlabel(f"{xlabel}") # α
ax.ylabel(f"{ylabel}")
ax.colorbar().ax.set_ylabel('Ks (m/hr)')
ax.savefig(os.path.join(output_path, "inf_adq-Ks_vs_dp-t_fitted.png"))

#ax.show()

    
    