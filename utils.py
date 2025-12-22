import os
import netCDF4 as nc


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

def delete_files(dirpath, key):
    for file in os.listdir(dirpath):
        if key in file:
            os.remove(os.path.join(dirpath, file))