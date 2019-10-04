import os, sys
import numpy as np
import xmltodict
import xarray as xr

# specific tools (which can be found here ../../Hub_tools/)
import image_tools as img
import validation_tools as vt

Default_flags = ['CLOUD', 'CLOUD_AMBIGUOUS', 'CLOUD_MARGIN',\
                 'INVALID', 'COSMETIC', 'SATURATED', 'SUSPECT',\
                 'HISOLZEN', 'HIGHGLINT', 'SNOW_ICE', 'AC_FAIL',\
                 'WHITECAPS']

def extract_OLCI_matchup_data(input_dir, match_lat, match_lon,\
                      insitu_product, pad_size=1,\
                      radiometry_type='Oa%s_reflectance', nchannels=21, BW06=True,
                      flags_to_use=Default_flags):
    '''
     Extracts OLCI data for match-up analysis
    '''
    # get geo-coordinates
    geo_file = os.path.join(input_dir, 'geo_coordinates.nc')
    nc_fid = xr.open_dataset(geo_file)
    LON  = nc_fid.longitude.data
    LAT  = nc_fid.latitude.data
    nc_fid.close()

    # get band centres and names
    xml_file = os.path.join(input_dir, 'xfdumanifest.xml')
    with open(xml_file) as fd:
        doc = xmltodict.parse(fd.read())
        bands = doc['xfdu:XFDU']['metadataSection']['metadataObject'][4]\
                    ['metadataWrap']['xmlData']['olci:olciProductInformation']\
                    ['olci:bandDescriptions']['sentinel3:band']
        
    band_names = []; band_wavs = []
    for band in bands:
        band_names.append(band['@name'])
        band_wavs.append(float(band['sentinel3:centralWavelength']))
    
    # get closest point to our validation point
    dist = img.spheric_dist(LAT, match_lat, LON, match_lon)
    i1, j1 = np.where((dist == np.nanmin(dist)))
    i1 = i1[0] ; j1 = j1[0]

    # get flags
    if 'reflectance' in radiometry_type:
        flag_file = os.path.join(input_dir, 'wqsf.nc')
        flag_var = 'WQSF'
    else:
        flag_file = os.path.join(input_dir, 'qualityFlags.nc')
        flag_var = 'quality_flags'
        
    flag_fid = xr.open_dataset(flag_file)
    flags = flag_fid.get(flag_var).data[i1-pad_size:i1+pad_size+1,j1-pad_size:j1+pad_size+1]
    flag_names = flag_fid.get(flag_var).flag_meanings.split(' ')
    flag_vals = flag_fid.get(flag_var).flag_masks
    flag_fid.close()

    # cut required variables from closest point and surrounding box
    extracted_vars = []
    if 'chl' in insitu_product:
        products = ['chl_nn.nc','chl_oc4me.nc']
        for product in products:
            var_file = os.path.join(input_dir, product)
            if os.path.exists(var_file):
                var_fid = xr.open_dataset(var_file)
                ext_var = var_fid.get(product.split(',')[0].upper())\
                          .data[i1-pad_size:i1+pad_size+1,j1-pad_size:j1+pad_size+1]
                var_fid.close()
                if BW06:
                    proc_var = vt.BW06_analysis(ext_var, flags_to_use,\
                                                flag_names, flag_vals, flags)
                else:
                    proc_var = np.nanmean(ext_var)
            else:
                if BW06:
                    proc_var = vt.BW06_analysis([], [], [], [], [], dummy=True)
                else:
                    proc_var = np.nan
            extracted_vars.append(proc_var)
            
    elif 'rrs' in insitu_product:    
        extracted_vars = [];
        for rad_channel_number in range(1, nchannels+1):
            rad_channel = radiometry_type % (str(rad_channel_number).zfill(2))
            var_file = os.path.join(input_dir, rad_channel + '.nc')
            if os.path.exists(var_file):
                var_fid = xr.open_dataset(var_file)
                ext_var = var_fid.get(rad_channel)\
                          .data[i1-pad_size:i1+pad_size+1,j1-pad_size:j1+pad_size+1]  
                var_fid.close()
                if BW06:
                    proc_var = vt.BW06_analysis(ext_var, flags_to_use,\
                                                flag_names, flag_vals, flags)
                else:
                    proc_var = np.nanmean(ext_var)
            else:
                if BW06:
                    proc_var = vt.BW06_analysis([], [], [], [], [], dummy=True)
                else:
                    proc_var = np.nan
            extracted_vars.append(proc_var) 
    else:
        print('Incorrect product')
    
    return band_wavs, extracted_vars 

def extract_CCI_matchup_data(insitu_product, insitu_record, rad_cols):
    '''
     Extracts CCI data for match-up analysis
    '''
    if 'chl' in insitu_product:
        pass
    elif 'rrs' in insitu_product:
        insitu_vals = []
        # read insitu data    
        for rr in rad_cols:
            try:
                insitu_vals.append(float(insitu_record[rr]))
            except:
                insitu_vals.append(np.nan)
    return insitu_vals