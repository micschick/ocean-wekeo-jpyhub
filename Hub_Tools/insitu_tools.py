import os
import numpy as np

def read_CCI_insitu_file(CCI_directory, insitu_product, insitu_resolution):
    '''
     READ CCI insitu validation data
    '''
    if insitu_product == 'chla':
        input_file = os.path.join(os.getcwd(), CCI_directory, 'insitudb_'\
                     + insitu_product + '.tab')
        max_cols = 7
    else:
        if insitu_resolution == 'fullrange':
            input_file = os.path.join(os.getcwd(), CCI_directory, 'insitudb_'\
                         + insitu_product + '.tab')
            if 'rrs' in insitu_product:
                max_cols = 615
            else:
                max_cols = 644
        else:
            input_file = os.path.join(os.getcwd(), CCI_directory, 'insitudb_'\
                         + insitu_product + '_satbands' + str(insitu_resolution) +'.tab')
            if 'rrs' in insitu_product:
                max_cols = 82
            else:
                max_cols = 275
                
    head_start = 0
    head_end = 0
    with open(input_file, 'r') as infile:
        count = -1
        for line in infile:
            count = count + 1
            if line.startswith('/*'):
                head_start = count
            if line.startswith('*/'):
                head_end = count

    insitu_file_read = np.genfromtxt(input_file, delimiter='\t', skip_header=head_end+1, usecols=range(max_cols), dtype='str')
    insitu_headers = insitu_file_read[0, :]
    insitu_data = insitu_file_read[1:, :]
    
    return insitu_headers, insitu_data

def get_CCI_coordinate_cols(insitu_headers):
    '''
     Get the indices for lat, lon and time columns
    '''
    lat_col = 0
    lon_col = 0
    time_col = 0

    count = -1
    for header in insitu_headers:
        count = count + 1
        if 'Lat' in header:
            lat_col = count
        if 'Lon' in header:
            lon_col = count
        if 'Time' in header:
            time_col = count

    return lat_col, lon_col, time_col

def get_CCI_radiometry_cols(insitu_headers, insitu_resolution):
    '''
     Get the indices for the OLCI band centres and radiometry
    '''
    rad_cols = []
    lambda_cols = []

    count = -1
    for header in insitu_headers:
        count = count + 1
        if 'OLC' in header:
            if 'RRS' in header:
                rad_cols.append(count)
            elif 'Lambda' in header:
                lambda_cols.append(count)

    return rad_cols, lambda_cols