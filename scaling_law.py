from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_absolute_error
import math
import pandas as pd
import numpy as np
import netCDF4 as nc
import SLOTH.sloth.toolBox
import SLOTH.sloth.IO
import SLOTH.sloth.PlotLib
import SLOTH.sloth.mapper
import os
import re
import shutil
from datetime import date, timedelta
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from utils import powerlaw_func, linear_law
import copy


plt.rcParams.update({'font.size': 22})

def export_parflow_csv():
    #name='/p/project/cslts/miaari1/python_scripts/parflow/comparewithanalytical/report/2ndinf_logkfit/exfiltration'
    name='/p/project/cslts/miaari1/python_scripts/parflow/comparewithanalytical/report/optimize_2nd/exfiltration'

    pressures = {}
    pressures["z"] = [z/100 for z in range(5,100, 10)]
    timesteps = [0, 1, 3, 5, 10, 15, 20, 30, 50, 75, 100]

    for t in timesteps:
        print(name + '.out.satur.'+('{:05d}'.format(t))+'.pfb')
        data = SLOTH.sloth.IO.read_pfb(name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb')
        pressures[f"time={t}"] = data[:,0,0]

    print(data[:,5,5].shape)
    print(len(pressures["z"]))
    df = pd.DataFrame(pressures)
    #df.to_csv('/p/project/cslts/miaari1/python_scripts/outputs/report/2ndinf_logkfit.csv', index=False)
    df.to_csv('/p/project/cslts/miaari1/python_scripts/outputs/report/optimize_2nd.csv', index=False)

def find_roots_from_plot(L):
    # Define the function
    def f(x):
        y = np.tan(x*L) + 2*x
        return y

    # Create an array of x-values
    x_values = np.linspace(0, 20, 100000000)
    print("calculating y values")
    # Compute the corresponding y-values
    y_values = f(x_values)

    zeros = [0] * len(x_values)
    line1 = LineString(np.column_stack((x_values, y_values)))
    line2 = LineString(np.column_stack((x_values, zeros)))
    print("looking for intersections")
    intersection = line1.intersection(line2)
    print("applying newton-raphson")
    roots = [p.x for p in intersection]
    roots = [x for x in roots if abs(f(x))< 1e-5]
    return roots
    
def residue(lambdaa, L, t, z):
    residue_value = (np.sin(lambdaa*z))*np.sin(lambdaa*L)*(np.exp(-t*(lambdaa*lambdaa)))/(1+(L/2)+(2*L*lambdaa*lambdaa))
    return residue_value

def sum_residue(L, t, z, roots):
    sum_values = []    
    for root in roots:
        resid = residue(root, L, t, z)
        sum_values.append(resid)
    residues = np.sum(sum_values)
    return residues
    
def plot_t0_analytical():
    # define variables in m and hr
    L_star = 1 #m
    Ks = 0.01 #m/hr
    c = 10 #1/m 
    Ss = 0.4
    Sr = 0.06
    qa_star = 0.00 #m/hr
    qb_star = 0.001 #m/hr
    pressure_0 = 0 #m

    # define dimensionless parameters
    L = c*L_star
    qa = qa_star/Ks
    qb = qb_star/Ks

    # get the roots for this L
    roots = find_roots_from_plot(L)

    # define layers
    z_values = [x/100 for x in range(5, L_star*100-4, 1)]
    data = {}
    data = {"z": z_values}

    # calculate for timesteps
    timesteps = [10000] # to make sure to be at steady state at t=1000hr
    for t_star in timesteps:
        data[f"time={t_star}"] = []
        t = (c*Ks*t_star)/(Ss-Sr)
        for z_star in data["z"]:            
            # for z_star=0 bottom layer, z_star=L_star top layer
            z = c*z_star
            
            residues = sum_residue(L, t, z, roots)
            
            K = qb - (qb - np.exp(c*pressure_0))*np.exp(-z) - 4*(qb - qa)*(np.exp((L-z)/2)*np.exp(-t/4))*residues

            K_star = K*Ks
            pressure = (np.log(K_star/Ks))/c
            data[f"time={t_star}"].append(pressure)
            
    df = pd.DataFrame(data)
    y_axis = df["z"].to_list()
    y_axis = [L_star-y for y in y_axis]
    timesteps = [10000]
    for t in timesteps:
        x_axis = df[f"time={t}"].to_list()
        plt.plot(x_axis, y_axis, "b")
        plt.annotate(f"t=0", xy=(x_axis[-1]-0.01, y_axis[-1]), color="b")

def analytical_pressure():
    # define variables in m and hr
    L_star = 1 #m
    Ks = 0.01 #m/hr
    c = 10 #/m
    Ss = 1
    Sr = 0.2
    qa_star = 0.001 #m/hr
    qb_star = 0.009 #m/hr
    pressure_0 = 0 #m

    # define dimensionless parameters
    L = c*L_star
    qa = qa_star/Ks
    qb = qb_star/Ks

    # get the roots for this L
    roots = find_roots_from_plot(L)

    # define layers
    z_values = [x/100 for x in range(5, L_star*100-4, 1)]
    data = {}
    data = {"z": z_values}

    # calculate for timesteps
    timesteps = [1, 3, 75, 100, 150, 200]
    for t_star in timesteps:
        data[f"time={t_star}"] = []
        t = (c*Ks*t_star)/(Ss-Sr)
        for z_star in data["z"]:            
            # for z_star=0 bottom layer, z_star=L_star top layer
            z = c*z_star
            
            residues = sum_residue(L, t, z, roots)
            
            K = qb - (qb - np.exp(c*pressure_0))*np.exp(-z) - 4*(qb - qa)*(np.exp((L-z)/2)*np.exp(-t/4))*residues

            K_star = K*Ks
            pressure = (np.log(K_star/Ks))/c
            data[f"time={t_star}"].append(pressure)
            
    df = pd.DataFrame(data)
    df.to_csv("/p/project/cslts/miaari1/python_scripts/outputs/report/2ndinf_SandY.csv", index=False)
    return

def plot_analytical_numerical():    
    # read from numerical solution
    #numerical_path = "/p/project/cslts/miaari1/python_scripts/outputs/report/2ndinf_logkfit.csv"
    numerical_path = "/p/project/cslts/miaari1/python_scripts/outputs/report/optimize_2nd.csv"
    numerical_solution = pd.read_csv(numerical_path)
    y_axis = numerical_solution["z"].to_list()
    y_axis = [1-x for x in y_axis]
    timesteps = [0, 1, 3, 20]
    line_type = {"0":"-", "1":"--", "3":":", "20":"-.", "150":"-."}
    for t in timesteps:
        x_axis = numerical_solution[f"time={t}"].to_list()
        plt.plot(x_axis, y_axis, f"r{line_type[f'{t}']}")
        if t!=0:
            plt.annotate(f"{t}", xy=(x_axis[-1], y_axis[-1]), color="r")
        elif t==0:
            plt.annotate(f"t={t}", xy=(x_axis[-1], y_axis[-1]), color="r")

    # read from analytical solution
    analytical_path = "/p/project/cslts/miaari1/python_scripts/outputs/report/2ndinf_SandY.csv"
    analytical_solution = pd.read_csv(analytical_path)
    y_axis = analytical_solution["z"].to_list()
    y_axis = [1-x for x in y_axis]
    timesteps = [1, 3, 150]
    for t in timesteps:
        x_axis = analytical_solution[f"time={t}"].to_list()
        plt.plot(x_axis, y_axis, f"b{line_type[f'{t}']}")
        if t==150:
            plt.annotate(f"{t}", xy=(x_axis[-1], y_axis[-10]), color="b")    
        else:
            plt.annotate(f"{t}", xy=(x_axis[-1], y_axis[-1]), color="b")

    plt.gca().invert_yaxis()
    plt.xlim([-0.3, 0])
    plt.ylim([1, 0])
    plt.xlabel("Pressure (m)")
    plt.ylabel("Soil depth (m)")
    line_up, = plt.plot([], 'b', label='Analytical model')
    line_down, = plt.plot([], 'r', label='Numerical model')
    plt.legend(handles=[line_up, line_down])
    plt.show()


def fit_scalinglaw():
    filepath = '/p/project/cslts/miaari1/python_scripts/outputs/testcases_appended.csv' # with watertable depth, n, alfa, theta, Ks from rosetta
    
    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    t_column = "time"
    df = pd.read_csv(filepath)
    x_axis = []
    y_axis = []

    for i in range(len(df)):
        q = df[q_column].iloc[i]
        k = df[k_column].iloc[i]
        d = df[d_column].iloc[i]
        n = df[n_column].iloc[i]
        alfa = df[alfa_column].iloc[i]
        theta_r = df[theta_r_column].iloc[i]
        theta_s = df[theta_s_column].iloc[i]
        t = df[t_column].iloc[i]
        if k>=q and t!=0 and (q/k)>=0.001:
            
            x = q/k
            y = (d/(t*k))
            #*((q/k)**(8.6))
            
            y_axis.append(y)
            x_axis.append(x)
    
    # Fit the function
    params, covariance = curve_fit(powerlaw_func, x_axis, y_axis)
    a_fit, b_fit = params
    print(f"Fitted a: {a_fit} and b:{b_fit}")

    # Generate data points for the fitted curve
    y_fit = [powerlaw_func(x, a_fit, b_fit) for x in x_axis]

    # calculate fitting accuracy
    R_square = r2_score(y_axis, y_fit)
    print(f"R2 = {R_square}")
    MSE = np.square(np.subtract(y_axis,y_fit)).mean()
    print(f"MSE = {MSE*10**-9} x10⁹")
    RMSE = math.sqrt(MSE)
    print(f"RMSE = {RMSE}")
    MAE = mean_absolute_error(y_axis, y_fit)
    print(f"MAE = {MAE}")
    
    x_fit = [min(x_axis), max(x_axis)]
    y_fit = [powerlaw_func(x, a_fit, b_fit) for x in x_fit]

    plt.figure(figsize=(16,9))
    plt.grid(True)
    plt.scatter(x_axis, y_axis, s=30, facecolors='k', edgecolors='k', label="data points")
    plt.plot(x_fit, y_fit, color="red", label="fitted line")
    #plt.annotate(f"R²={round(R_square, 3)}\nf(x)={round(a_fit, 3)}x^({round(b_fit, 3)})", xy=(10, 0.1), color="black")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("q/k=Kr (-)")
    plt.ylabel("(tq/d)*(q/k)^8.6 (-)") # α
    plt.legend()
    plt.show()

def soil_type():
    filepath = '/p/project/cslts/miaari1/python_scripts/outputs/testcases_appended.csv' # with watertable depth, n, alfa, theta, Ks from rosetta
    
    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    t_column = "time"
    df = pd.read_csv(filepath)
        
    soil_types = list(df[k_column].unique())
    soil_types.sort()

    ax = plt
    ax.figure(figsize=(16,9))

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
            t = df_soil[t_column].iloc[i]
            if k>=q and t!=0 and (q/k)>=0.001:
                x = t
                y = k
                y_axis.append(y)
                x_axis.append(x)
        print(soil)
        ax.scatter(x_axis, y_axis, s=30, label=f"Ks={round(soil, 4)}")

    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    ax.xlabel("t (hr)")
    ax.ylabel("k (m/hr)") # α
    ax.legend()
    ax.show()

def soil_type_colormap():
    filepath = '/p/project/cslts/miaari1/python_scripts/outputs/testcases_appended.csv' # with watertable depth, n, alfa, theta, Ks from rosetta
    
    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    t_column = "time"
    df = pd.read_csv(filepath)
        
    soil_types = list(df[k_column].unique())
    soil_types.sort()

    ax = plt
    ax.figure(figsize=(16,9))

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
            if k>=q and t!=0 and (q/k)>=0.001:
                x = alfa*d
                y = t*q*alfa
                #x = q/k
                #y = k*t*alfa/(theta_s-theta_r)
                y_axis.append(y)
                x_axis.append(x)
                index += 1
        print(soil)
        x_all.extend(x_axis)
        y_all.extend(y_axis)
        colors.extend([soil]*len(x_axis))

    ax.scatter(x_all, y_all,c=colors, cmap='jet', s=30, norm=matplotlib.colors.LogNorm())
    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    ax.xlabel("x")
    ax.ylabel("y") # α
    ax.colorbar().ax.set_ylabel('Ks (m/hr)')
    ax.show()

def exf_soil_type_colormap():
    filepath = '/p/project/cslts/miaari1/python_scripts/outputs/testcases_drainage.csv' # with watertable depth, n, alfa, theta, Ks from rosetta
    
    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    t_column = "exf_time"
    df = pd.read_csv(filepath)
        
    soil_types = list(df[k_column].unique())
    soil_types.sort()

    ax = plt
    ax.figure(figsize=(16,9))

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
            if k>=q and t!=0 and (q/k)>=0.001:
                x = alfa*d
                y = t*k*alfa
                #x = q/k
                #y = k*t*alfa/(theta_s-theta_r)
                y_axis.append(y)
                x_axis.append(x)
                index += 1
        print(soil)
        x_all.extend(x_axis)
        y_all.extend(y_axis)
        colors.extend([soil]*len(x_axis))

    ax.scatter(x_all, y_all,c=colors, cmap='jet', s=30, norm=matplotlib.colors.LogNorm())
    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    ax.xlabel("x")
    ax.ylabel("y") # α
    ax.colorbar().ax.set_ylabel('Ks (m/hr)')
    ax.show()

def constantq_constantd():
    filepath1 = '/p/project/cslts/miaari1/python_scripts/outputs/testcases_constantq_constantd.csv' # without watertable depth
    filepath2 = '/p/project/cslts/miaari1/python_scripts/outputs/testcases_constantq_constantd_v3.csv' # without watertable depth

    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    t_column = "time"
    # first file
    df1 = pd.read_csv(filepath1)
    df1["q/k"] = df1[q_column]/df1[k_column]
    df1 = df1[df1["q/k"]<=1]
    df1.reset_index(drop=True, inplace=True)

    # second file
    df2 = pd.read_csv(filepath2)
    df2["q/k"] = df2[q_column]/df2[k_column]
    df2 = df2[df2["q/k"]<=1]
    df2.reset_index(drop=True, inplace=True)

    k1 = df1[k_column].to_list()
    t1 = df1[t_column].to_list()
    k2 = df2[k_column].to_list()
    t2 = df2[t_column].to_list()


    plt.figure(figsize=(16,9))
    plt.grid(True)
    plt.scatter(k1, t1, s=30, facecolors='b', edgecolors='b', label="q=0.0002 m/hr")
    plt.scatter(k2, t2, s=30, facecolors='r', edgecolors='r', label="q=0.001 m/hr")
    plt.xscale("log")
    #plt.yscale("log")
    plt.xlabel("Ks (m/hr)")
    plt.ylabel("t (hr)")
    plt.title("d=4m & constant q")
    plt.legend()
    plt.show()

def sensitivity_alfa_n_k(sens_alfa=None, sens_n=None, sens_k=None, sens_alfan=None):
    def plot5(q, k, d, n, alfa, theta_r, theta_s, t):
        x = (q)/(k*alfa*d)
        y = t
        xlabel = "q/(k*d*alpha)"
        ylabel = "t"
        return x, y, xlabel, ylabel
    inf_cases_path = '/p/project/cslts/miaari1/python_scripts/outputs/infiltration_testcases.csv'
    if sens_alfa:
        title = "Sensitivity analysis for alpha"
        var_alfan_path = '/p/project/cslts/miaari1/python_scripts/parflow/sensitivity_analysis/infiltration_sensitivity_alfa.csv'
    elif sens_n:
        title = "Sensitivity analysis for n"
        var_alfan_path = '/p/project/cslts/miaari1/python_scripts/parflow/sensitivity_analysis/infiltration_sensitivity_n.csv'
    elif sens_k:
        title = "Sensitivity analysis for saturated hydraulic conductivity Ks"
        var_alfan_path = '/p/project/cslts/miaari1/python_scripts/parflow/sensitivity_analysis/infiltration_sensitivity_k.csv'    
    elif sens_alfan:
        title = "Sensitivity analysis for alpha and n"
        var_alfan_path = '/p/project/cslts/miaari1/python_scripts/parflow/sensitivity_analysis/infiltration_sensitivity_alfa_n.csv'    

    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    t_column = "time"
    df = pd.read_csv(inf_cases_path)

    soil_types = list(df[k_column].unique())
    soil_types.sort()

    ax = plt
    ax.figure(figsize=(16,9))
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
            if k>=q and t!=0 and (q/k)>=0.001:
                x, y, xlabel, ylabel = plot5(q, k, d, n, alfa, theta_r, theta_s, t)
                y_axis.append(y)
                x_axis.append(x)
                index += 1

        x_all.extend(x_axis)
        y_all.extend(y_axis)
        colors.extend([soil]*len(x_axis))

    
    ax.scatter(x_all, y_all,c=colors, cmap='jet', s=30, norm=matplotlib.colors.LogNorm())

    var_alfan_df = pd.read_csv(var_alfan_path)
    val_x = []
    val_y = []
    max_err_x = []
    max_err_y = []
    min_err_x = []
    min_err_y = []

    for i in range(len(var_alfan_df)):
        q = var_alfan_df.iloc[i,1]
        k = var_alfan_df.iloc[i,2]
        d = var_alfan_df.iloc[i,3]
        n = var_alfan_df.iloc[i,5]
        alfa = var_alfan_df.iloc[i,4]
        theta_r = var_alfan_df.iloc[i,6]
        theta_s = var_alfan_df.iloc[i,7]
        t = var_alfan_df.iloc[i,8]
        x, y, xlabel, ylabel = plot5(q, k, d, n, alfa, theta_r, theta_s, t)
        val_x.append(x)
        val_y.append(y)

        q = var_alfan_df.iloc[i,9]
        k = var_alfan_df.iloc[i,10]
        d = var_alfan_df.iloc[i,11]
        n = var_alfan_df.iloc[i,13]
        alfa = var_alfan_df.iloc[i,12]
        theta_r = var_alfan_df.iloc[i,14]
        theta_s = var_alfan_df.iloc[i,15]
        t = var_alfan_df.iloc[i,16]
        x, y, xlabel, ylabel = plot5(q, k, d, n, alfa, theta_r, theta_s, t)
        max_err_x.append(x)
        max_err_y.append(y)
        if var_alfan_df.iloc[i,16] != np.float("nan"):
            q = var_alfan_df.iloc[i,17]
            k = var_alfan_df.iloc[i,18]
            d = var_alfan_df.iloc[i,19]
            n = var_alfan_df.iloc[i,21]
            alfa = var_alfan_df.iloc[i,20]
            theta_r = var_alfan_df.iloc[i,22]
            theta_s = var_alfan_df.iloc[i,23]
            t = var_alfan_df.iloc[i,24]
            x, y, xlabel, ylabel = plot5(q, k, d, n, alfa, theta_r, theta_s, t)
            min_err_x.append(x)
            min_err_y.append(y)

    for i in range(len(var_alfan_df)):
        ax.plot([min_err_x[i], val_x[i], max_err_x[i]], [min_err_y[i], val_y[i], max_err_y[i]], marker='+', c="k", linewidth=2.0)#, label=f'Varied Input: {alfa} and {n}')

    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    ax.xlabel(f"{xlabel}") # α
    ax.ylabel(f"{ylabel}")
    ax.colorbar().ax.set_ylabel('Ks (m/hr)')
    ax.title(title)
    ax.show()



def plot_wet_dry_law(**kwargs):
    def get_data(plot_nb):
        exf_cases_path = '/p/project/cslts/miaari1/python_scripts/outputs/drainage_inf_testcases.csv'
        q_column = "q"
        k_column = "k"
        d_column = "d"
        n_column = "n"
        alfa_column = "alfa"
        theta_r_column = "theta_r"
        theta_s_column = "theta_s"
        inf_t_column = "inf_time"
        exf_t_column = "exf_time"
        df = pd.read_csv(exf_cases_path)

        soil_types = list(df[k_column].unique())
        soil_types.sort()

        length = 0
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
                inf_t = df_soil[inf_t_column].iloc[i]
                exf_t = df_soil[exf_t_column].iloc[i]
                if k>=q and exf_t!=0 and (q/k)>=0.001:
                    if plot_nb==0:
                        x = q*alfa*d/k
                        y = inf_t*q/d
                        xlabel = "q*α*d/k (-)"
                        ylabel = "t_wet*q/d (-)"
                    elif plot_nb==1:
                        x = q*alfa*d/k
                        y = exf_t*q/d
                        xlabel = "q*α*d/k (-)"
                        ylabel = "t_dry*q/d (-)"
                    length += 1
                    y_axis.append(y)
                    x_axis.append(x)
                    index += 1

            x_all.extend(x_axis)
            y_all.extend(y_axis)
            colors.extend([soil]*len(x_axis))
        return x_all, y_all, xlabel, ylabel, colors

    fig = plt
    fig.figure(figsize=(16,9))
    gs  = fig.GridSpec(nrows=1,ncols=2)
    gs.update(wspace=0.2, hspace=0.4)
    
    ax = fig.subplot(gs[0, 0])
    x_all, y_all, xlabel, ylabel, colors = get_data(0)
    ax.scatter(x_all, y_all,c=colors, cmap='jet', s=30, norm=matplotlib.colors.LogNorm())
    ax.grid(True)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(f"{xlabel}")
    ax.set_ylabel(f"{ylabel}")

    ax2 = fig.subplot(gs[0, 1])
    x_all, y_all, xlabel, ylabel, colors = get_data(1)
    sp = ax2.scatter(x_all, y_all,c=colors, cmap='jet', s=30, norm=matplotlib.colors.LogNorm())
    ax2.grid(True)
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_xlabel(f"{xlabel}")
    ax2.set_ylabel(f"{ylabel}")


    #sp = ax2.scatter(x_all, y_all,c=colors, cmap='jet', s=30, norm=matplotlib.colors.LogNorm())
    cbar_ax = fig.axes([0.91, 0.15, 0.03, 0.7])
    fig.colorbar(sp, cax=cbar_ax).ax.set_ylabel('Ks (m/hr)')
    fig.show()

def fit_SWCC():
    def van_genuchten_maulem_k(h, alfa, n):
        m = 1-1/n
        k = ((1-((alfa*h)**(n-1))*((1+(alfa*h)**(n))**(-m)))**(2))/((1+(alfa*h)**(n))**(m/2))
        return k#np.log(k)
    def gardner_k(h, c):
        k = np.exp(-c*h)
        return k#np.log(k)
    
    c = 10 #/m
    
    h_values = [h/1000 for h in range(1, 1000)]
    y_values = [gardner_k(x, c) for x in h_values]

    params, covariance = curve_fit(van_genuchten_maulem_k, h_values, y_values, bounds=[[0.01, 0.5],[15, 10]])
    alfa, n = params
    
    print(f"Fitted alfa: {alfa}, n: {n}")

    y_VG = [van_genuchten_maulem_k(h, alfa, n) for h in h_values]
    y_values = [gardner_k(x, c) for x in h_values]
    
    # Plot the original data and the fitted curve
    plt.plot(h_values, y_VG, "r", label='Van Genuchten-Mualem model')
    plt.plot(h_values, y_values, "blue", label='Gardner model')
    plt.ylabel('Relative hydraulic conductivity Kr (-)')
    plt.xlabel('Pressure (m)')
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.grid(True)
    plt.show()


def SWCC_theta():
    def van_genuchten_maulem(h, theta_s, theta_r, alfa, n):
        m = 1-1/n
        theta = theta_r + (theta_s - theta_r)/((1+(alfa*h)**n)**(m))
        return theta
    def gardner(h, theta_s, theta_r, c):
        theta = theta_r + (theta_s - theta_r)*(np.exp(-c*h))
        return theta

    # van genuchten variables obtained from the fitting in previous step Kr vs pressure
    alfa = 3.698059609737854 #/m
    n = 1.8990947179855633 #(-)
    
    #alfa = 3.917989460214275 #/m
    #n = 2.622947993712746 #(-)
    # gardner parameter
    c = 10 #/m
    
    # soil parameters
    theta_s = 0.4
    theta_r = 0.06

    h_values = [h/1000 for h in range(1, 1000)]

    y_G = [gardner(x, theta_s, theta_r, c) for x in h_values]
    y_VGM = [van_genuchten_maulem(x, theta_s, theta_r, alfa, n) for x in h_values]
    
    # Plot the original data and the fitted curve
    plt.plot(h_values, y_VGM, "r", label='Van Genuchten-Mualem model')
    plt.plot(h_values, y_G, "blue", label='Gardner model')
    plt.ylabel('Soil water content θ (m³/m³)')
    plt.xlabel('Pressure (m)')
    plt.xscale("log")
    plt.legend()
    plt.grid(True)
    plt.show()

#export_parflow_csv()
#plot_t0_analytical()
#plot_analytical_numerical()
SWCC_theta()