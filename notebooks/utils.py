import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed") # due a RuntimeWarning with numpy.dtype
import pandas as pd
import numpy as np
import math

def mdia_to_xyz_minCurve(deviation):
    """
    A function to convert a well deviation (path of csv file) given in MD[m], incl[deg], azi[deg]
    into an xyz array with x[m], y[m], z[m] using the minimum curvature method
    according to: [drillingformulas.com](http://bit.ly/2MNp7U0)
    """
    # read data
    data = pd.read_csv(deviation, sep=',', header='infer')
    # clean data
    data.drop(columns=['Unnamed: 0'], inplace=True)
    #data['Dogleg [deg/30m]'].replace(np.nan, 0, inplace=True)
    # add columns needed for calculations
    data['Dogleg_rad [rad/30m]'] = np.radians(data['Dogleg [deg/30m]'])
    data['RatioFactor'] = (2 / data['Dogleg_rad [rad/30m]']) * np.tan(data['Dogleg_rad [rad/30m]'] / 2)
    # calculate intervals
    delta_MD = np.array(data['MD[m]'][1:]) - np.array(data['MD[m]'][:-1])
    # get uppers/lowers
    RF_lower = np.array(data['RatioFactor'][1:])
    incl_upper = np.array(data['Inc[deg]'][:-1])
    incl_lower = np.array(data['Inc[deg]'][1:])
    azi_upper = np.array(data['Azi[deg]'][:-1])
    azi_lower = np.array(data['Azi[deg]'][1:])
    # calculate xyz
    east_x = delta_MD / 2 * (np.sin(np.radians(incl_upper)) * np.sin(np.radians(azi_upper)) 
                             + np.sin(np.radians(incl_lower)) * np.sin(np.radians(azi_lower))) * RF_lower
    
    north_y = delta_MD / 2 * (np.sin(np.radians(incl_upper)) * np.cos(np.radians(azi_upper)) 
                             + np.sin(np.radians(incl_lower)) * np.cos(np.radians(azi_lower))) * RF_lower
    
    TVD_z = delta_MD / 2 * (np.cos(np.radians(incl_upper)) + np.cos(np.radians(incl_lower))) * RF_lower
    
    return east_x, north_y, TVD_z
