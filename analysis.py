from utils import powerlaw_func, linear_law
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
# ==============================================================================
# CATEGORY 5: ANALYSIS AND CURVE FITTING
# Functions for statistical analysis and power-law fitting
# ==============================================================================

def fitting_func(xlist, ylist):
    """
    Fit a power-law function to data using log-linearization.
    
    Fits the model: y = a * x^b
    by transforming to linear space: ln(y) = ln(a) + b * ln(x)
    
    Parameters:
    -----------
    xlist : list or array
        Independent variable data
    ylist : list or array
        Dependent variable data
    
    Returns:
    --------
    tuple
        (x_fit, y_fit, a_fit, b_fit, r2)
        - x_fit: x values for plotting fitted line
        - y_fit: fitted y values
        - a_fit: power-law coefficient
        - b_fit: power-law exponent
        - r2: coefficient of determination (R²)
    """
    # Linearize by taking logarithms
    y_lin = np.log(ylist)
    x_lin = np.log(xlist)
    
    # Fit linear model in log space
    params, covariance = curve_fit(linear_law, x_lin, y_lin)
    a_fit, b_fit = params
    
    # Calculate fitting accuracy
    y_fit_lin = [linear_law(x, a_fit, b_fit) for x in x_lin]
    R_square = r2_score(y_lin, y_fit_lin)
    
    # Calculate Pearson correlation coefficient
    r = np.corrcoef(y_lin, y_fit_lin)
    r2 = r[0][1]**2
    
    # Back-transform to power law
    a_fit = np.exp(a_fit)
    
    # Generate fitted curve for plotting
    x_fit = [min(xlist), max(xlist)]
    y_fit = [powerlaw_func(x, a_fit, b_fit) for x in x_fit]
    
    return x_fit, y_fit, a_fit, b_fit, r2
