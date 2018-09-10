import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed") # due a RuntimeWarning with numpy.dtype
import pandas as pd
import numpy as np

def mdia_to_xyz_minCurve(deviation):
    """
    Returns a tuple whose elements are 1D arrays corresponding to the sequences x, y, and tvd_z values
    from a well path deviation file (csv) with columns MD[m], incl[m], azi[deg].
    Uses the minimum curvature method, [drillingformulas.com](http://bit.ly/2MNp7U0)

    Args:
        str. The path to the deviation file (csv format)
    Returns:
        tuple. 3-elements of x, y, and tvd_z points.
    """
    # read data
    data = pd.read_csv(deviation, sep=',', header='infer')
    # clean data
    data.drop(columns=['Unnamed: 0'], inplace=True)
    # add columns needed for calculations
    data['Dogleg_rad [rad/30m]'] = np.radians(data['Dogleg [deg/30m]'])
    data['RatioFactor'] = (2 / data['Dogleg_rad [rad/30m]']) * np.tan(data['Dogleg_rad [rad/30m]'] / 2)
    data['RatioFactor'].replace(np.nan, 1, inplace=True)
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


def top_incAziXYZ(tops_df, survey_edt):
    """
    Takes a tops dataframe and returns same dataframe with addtional columns
    for Inc and Azi based on:
    inc = inc_dev_lower - ((MD_dev_lower - MD_Top) * (delta_incl) / (delta_MD_dev))
    azi = azi_dev_lower - ((MD_dev_lower - MD_Top) * (delta_azi) / (delta_MD_dev))
    """
    for top in tops_df.itertuples():
        current_Top_MD = top[1] 
        for depth in survey_edt.itertuples():
            current_dev_MD = depth[1]
            if current_Top_MD < current_dev_MD:
                # get row above and below top
                enclosing_rows = survey_edt[depth[0] - 1:depth[0] + 1]
                # grab MD, incl, azi values from rows
                MD_dev_low = enclosing_rows[1:2]['MD[m]']
                MD_dev_up = enclosing_rows[:1]['MD[m]']
                Inc_dev_low = enclosing_rows[1:2]['Inc[deg]']
                Inc_dev_up = enclosing_rows[:1]['Inc[deg]']
                Azi_dev_low = enclosing_rows[1:2]['Azi[deg]']
                Azi_dev_up = enclosing_rows[:1]['Azi[deg]']
                # grab x,y,z from rows
                x_dev_low = enclosing_rows[1:2]['x[m]']
                x_dev_up = enclosing_rows[:1]['x[m]']
                y_dev_low = enclosing_rows[1:2]['y[m]']
                y_dev_up = enclosing_rows[:1]['y[m]']
                z_dev_low = enclosing_rows[1:2]['z[m]']
                z_dev_up = enclosing_rows[:1]['z[m]']
                # assign calculated value to correct column of correct row
                # use .loc to avoid warning about assigning to a copy
                tops_df.loc[top[:1][0],'Inc[deg]'] = Inc_dev_low.values[0] - ((MD_dev_low.values[0] - current_Top_MD) 
                                                                * (Inc_dev_low.values[0] - Inc_dev_up.values[0]) 
                                                                / (MD_dev_low.values[0] - MD_dev_up.values[0]))
                tops_df.loc[top[:1][0],'Azi[deg]'] = Azi_dev_low.values[0] - ((MD_dev_low.values[0] - current_Top_MD) 
                                                                * (Azi_dev_low.values[0] - Azi_dev_up.values[0]) 
                                                                / (MD_dev_low.values[0] - MD_dev_up.values[0]))
                
                tops_df.loc[top[:1][0],'x[m]'] = x_dev_low.values[0] - ((MD_dev_low.values[0] - current_Top_MD) 
                                                                * (x_dev_low.values[0] - x_dev_up.values[0]) 
                                                                / (MD_dev_low.values[0] - MD_dev_up.values[0]))
                tops_df.loc[top[:1][0],'y[m]'] = y_dev_low.values[0] - ((MD_dev_low.values[0] - current_Top_MD) 
                                                                * (y_dev_low.values[0] - y_dev_up.values[0]) 
                                                                / (MD_dev_low.values[0] - MD_dev_up.values[0]))
                tops_df.loc[top[:1][0],'z[m]'] = z_dev_low.values[0] - ((MD_dev_low.values[0] - current_Top_MD) 
                                                                * (z_dev_low.values[0] - z_dev_up.values[0]) 
                                                                / (MD_dev_low.values[0] - MD_dev_up.values[0]))
                break
    return tops_df