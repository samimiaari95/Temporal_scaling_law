import math
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_absolute_error


def make_dir(dir_path):
    """
    Creates folders for every directory in the specified path, pass if the directory already exists
    Arguments:
    -   dir_path: path defined as a string
    Returns:
    -   Nothing
    Author: Sami
    """
    # iterate through the defined path
    while not os.path.isdir(dir_path):
        # check if the parent directory exists
        if not os.path.isdir(os.path.dirname(dir_path)):
            make_dir(os.path.dirname(dir_path))
        else:
            os.mkdir(dir_path)

def fit_powerlaw(x_axis, y_axis):
    def powerlaw(h, a, b):
        y = a*(h**b)
        #y = (a)*np.exp(-b*h)
        return y
    
    params, covariance = curve_fit(powerlaw, x_axis, y_axis)
    # Extract the fitted parameters
    a_fit, b_fit = params

    # Print the results
    print(f"Fitted a: {a_fit} and b:{b_fit}")

    y_fit = [powerlaw(x, a_fit, b_fit) for x in x_axis]

    R_square = r2_score(y_axis, y_fit)
    print(f"R2 = {R_square}")
    #return R_square
    MSE = np.square(np.subtract(y_axis,y_fit)).mean()
    print(f"MSE = {MSE*10**-9} x10⁹")
    RMSE = math.sqrt(MSE)
    print(f"RMSE = {RMSE}")
    MAE = mean_absolute_error(y_axis, y_fit)
    print(f"MAE = {MAE}")
    
    x_fit = [min(x_axis), max(x_axis)]
    y_fit = [powerlaw(x, a_fit, b_fit) for x in x_fit]

    return x_fit, y_fit

def plotlog_show(x_axis, y_axis, color="k", linestyle="-", label=None, xlabel="x", ylabel="y", scatter=True):
    plt.figure(figsize=(16,9))
    plt.grid(True)

    if label and scatter:
        plt.scatter(x_axis, y_axis, color=f"{color}{linestyle}", label=label)
    elif not label and scatter:
        plt.scatter(x_axis, y_axis, color=f"{color}{linestyle}")
    elif label and not scatter:
        plt.plot(x_axis, y_axis, color=f"{color}{linestyle}", label=label)
    elif not label and not scatter:
        plt.plot(x_axis, y_axis, color=f"{color}{linestyle}")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if label:
        plt.legend()
    plt.show()

def powerlaw_func(h, a, b):
        y = a*(h**b)
        return y

def linear_law(x, a, b) :
        return a + x * b

def read_nc(filepath, var):
    ncfile = nc.Dataset(filepath)
    return ncfile[var][:]

def open_nc(filepath):
    ncfile = nc.Dataset(filepath)
    variables = ncfile.variables
    for var in variables:
        print(ncfile[var])

def get_prudenceMask(lat2D, lon2D, prudName):
    """ return a prudance mask

    Return a boolean mask-array (True = masked, False = not masked) based on
    a passed set of longitude and latitude values and the name of the prudence
    region.
    The shape of the mask-array is set equal to the shape of input lat2D.
    Source: http://prudence.dmi.dk/public/publications/PSICC/Christensen&Christensen.pdf p.38

    Input values:
    -------------
    lat2D:    ndarray
        2D latitude information for each pixel
    lon2D:    ndarray
        2D longitude information for each pixel
    prudName: str
        Short name of prudence region

    Return value:
    -------------
    prudMask: ndarray
        Ndarray of dtype boolean of the same shape as lat2D.
        True = masked; False = not masked
    """
    if (prudName=='BI'):
        prudMask = np.where((lat2D < 50.0) | (lat2D > 59.0)  | (lon2D < -10.0) | (lon2D >  2.0), False, True)
    elif (prudName=='IP'):
        prudMask = np.where((lat2D < 36.0) | (lat2D > 44.0)  | (lon2D < -10.0) | (lon2D >  3.0), False, True)
    elif (prudName=='FR'):
        prudMask = np.where((lat2D < 44.0) | (lat2D > 50.0)  | (lon2D < -5.0) | (lon2D >  5.0), False, True)
    elif (prudName=='ME'):
        prudMask = np.where((lat2D < 48.0) | (lat2D > 55.0)  | (lon2D < 2.0) | (lon2D >  16.0), False, True)
    elif (prudName=='SC'):
        prudMask = np.where((lat2D < 55.0) | (lat2D > 70.0)  | (lon2D < 5.0) | (lon2D >  30.0), False, True)
    elif (prudName=='AL'):
        prudMask = np.where((lat2D < 44.0) | (lat2D > 48.0)  | (lon2D < 5.0) | (lon2D >  15.0), False, True)
    elif (prudName=='MD'):
        prudMask = np.where((lat2D < 36.0) | (lat2D > 44.0)  | (lon2D < 3.0) | (lon2D >  25.0), False, True)
    elif (prudName=='EA'):
        prudMask = np.where((lat2D < 44.0) | (lat2D > 55.0)  | (lon2D < 16.0) | (lon2D >  30.0), False, True)
    else:
        print(f'prudance region {prudName} not found --> EXIT')
    return prudMask

def get_S4W_basin(lat2D, lon2D, region):
    """ return a prudance mask

    Return a boolean mask-array (True = masked, False = not masked) based on
    a passed set of longitude and latitude values and the name of the prudence
    region.
    The shape of the mask-array is set equal to the shape of input lat2D.
    Source: http://prudence.dmi.dk/public/publications/PSICC/Christensen&Christensen.pdf p.38

    Input values:
    -------------
    lat2D:    ndarray
        2D latitude information for each pixel
    lon2D:    ndarray
        2D longitude information for each pixel
    prudName: str
        Short name of prudence region

    Return value:
    -------------
    prudMask: ndarray
        Ndarray of dtype boolean of the same shape as lat2D.
        True = masked; False = not masked
    """
    if (region=='SEINE'):
        prudMask = np.where((lat2D < 47.0) | (lat2D > 50.0)  | (lon2D < -2.0) | (lon2D >  3.0), False, True)
    elif (region=='IP'):
        prudMask = np.where((lat2D < 36.0) | (lat2D > 44.0)  | (lon2D < -10.0) | (lon2D >  3.0), False, True)
    elif (region=='FR'):
        prudMask = np.where((lat2D < 44.0) | (lat2D > 50.0)  | (lon2D < -5.0) | (lon2D >  5.0), False, True)
    elif (region=='ME'):
        prudMask = np.where((lat2D < 48.0) | (lat2D > 55.0)  | (lon2D < 2.0) | (lon2D >  16.0), False, True)
    elif (region=='SC'):
        prudMask = np.where((lat2D < 55.0) | (lat2D > 70.0)  | (lon2D < 5.0) | (lon2D >  30.0), False, True)
    elif (region=='AL'):
        prudMask = np.where((lat2D < 44.0) | (lat2D > 48.0)  | (lon2D < 5.0) | (lon2D >  15.0), False, True)
    elif (region=='MD'):
        prudMask = np.where((lat2D < 36.0) | (lat2D > 44.0)  | (lon2D < 3.0) | (lon2D >  25.0), False, True)
    elif (region=='EA'):
        prudMask = np.where((lat2D < 44.0) | (lat2D > 55.0)  | (lon2D < 16.0) | (lon2D >  30.0), False, True)
    else:
        print(f'prudance region {region} not found --> EXIT')
    return prudMask

def delete_files(dirpath, key):
    for file in os.listdir(dirpath):
        if key in file:
            os.remove(os.path.join(dirpath, file))