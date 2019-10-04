import os, sys
import numpy as np

def flag_data_fast(flags_to_use, flag_names, flag_values, flags):
    '''
     Builds a flag mask from specific flags. 
    '''
    flag_bits = np.uint64()
    flag_bits = np.dtype(type(flags[0][0])).type(flag_bits)
    for flag in flags_to_use:
        try:
            flag_bits = flag_bits | flag_values[flag_names.index(flag)]
        except:
            print(flag + " not present")
    
    return (flags & flag_bits) > 0

def BW06_analysis(var, flags_to_use, flag_names, flag_vals, flags, dummy=False):
    '''
     Apply Bailey and Werdell, 2006 match-up criteria
     https://www.sciencedirect.com/science/article/abs/pii/S0034425706000472
    '''
    # int dictionary for return values
    BW_vals = {}
    BW_vals['nvals'] = np.nan
    BW_vals['var_mean'] = np.nan
    BW_vals['var_median'] = np.nan
    BW_vals['var_std'] = np.nan
    BW_vals['filtered_var_mean'] = np.nan
    BW_vals['filtered_var_median'] = np.nan
    BW_vals['filtered_var_std'] = np.nan
    BW_vals['CV'] = np.nan
    
    if dummy:
        return BW_vals
    
    # logical array for good points
    logic = np.ones(np.shape(var))

    # get the flags
    flag_mask = flag_data_fast(flags_to_use, flag_names, flag_vals, flags)
    flag_mask = flag_mask.astype(float)
    flag_mask[flag_mask == 0.0] = np.nan

    # apply the flags
    var[np.isfinite(flag_mask)] = np.nan
    logic[np.isfinite(flag_mask)] = 0.0
    
    # calculate the unfiltered mean, median, std
    BW_vals['var_mean'] = np.nanmean(var)
    BW_vals['var_median'] = np.nanmedian(var)
    BW_vals['var_std'] = np.nanstd(var)

    # check if enough points to continue
    if np.nansum(np.isfinite(var)) < (np.nansum(np.isfinite(logic)) / 2) or \
       np.nansum(np.isfinite(var)) < 5:
        print('Not enough points to calculate filtered values')
        return BW_vals
        
    # calculate filtered mean, media, std
    logic[var > (var+(np.std(var)*1.5))] = 0.0
    logic[var < (var-(np.std(var)*1.5))] = 0.0
    var[logic == 0.0] = np.nan

    BW_vals['nvals'] = np.nansum(np.isfinite(var))
    BW_vals['filtered_var_mean'] = np.nanmean(var)
    BW_vals['filtered_var_median'] = np.nanmedian(var)
    BW_vals['filtered_var_std'] = np.nanstd(var)
    
    # calculate coefficient of variation
    BW_vals['CV'] = BW_vals['filtered_var_std'] / BW_vals['filtered_var_mean']

    return BW_vals