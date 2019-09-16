from ipywidgets import Layout, Box, FloatSlider
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import xarray as xr

def plot_eq(vmin, vmax, **kwargs):

    nchannels = int(len(kwargs)/3)    
    RED_list = []
    GREEN_list = []
    BLUE_list = []
    
    for ii in range(nchannels):
        RED_list.append(kwargs['red'+str(ii)])
        GREEN_list.append(kwargs['green'+str(ii)])
        BLUE_list.append(kwargs['blue'+str(ii)])
            
    REDS = np.asarray(RED_list)
    GREENS = np.asarray(GREEN_list)
    BLUES = np.asarray(BLUE_list)
   
    x = np.linspace(vmin, vmax, nchannels)

    plt.figure(figsize=(15,15))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[1,4])
    
    RG_sensitivity = REDS + GREENS
    RG_sensitivity = (RG_sensitivity - min(RG_sensitivity)) / (max(RG_sensitivity) - min(RG_sensitivity))
        
    # plot 1
    ax1 = plt.subplot(gs[0, 0])
    pr, = plt.plot(x, REDS, 'r')
    pg, = plt.plot(x, GREENS, 'g')
    pb, = plt.plot(x, BLUES, 'b')
    pn, = plt.plot(x, (BLUES+GREENS+REDS)/3, 'k--')
    pbr, = plt.plot(x, RG_sensitivity, '0.5', linestyle='--')
    plt.ylim(0, 1.1)
    plt.xlabel('Channel')
    plt.ylabel('intensity')
    leg1 =plt.legend([pr, pg, pb, pn, pbr],['red','green','blue','norm','r/g sens'])
    plt.xticks(x)
    plt.xlim([vmin, vmax])
    plt.ylim([0, 1.25])
    for ii in range(nchannels):
    	plt.plot([x[ii],x[ii]],[0, 1.25],'k--', linewidth=0.5)
    	plt.text(x[ii],1.275,str(ii),color='0.5')
    
    # plot 2
    ax2 = plt.subplot(gs[1,0])
    SCALING = np.linspace(0,1,len(REDS))
    for ii in range(len(REDS)):
        red_val = (SCALING[ii], REDS[ii], REDS[ii])
        green_val = (SCALING[ii], GREENS[ii], GREENS[ii])
        blue_val = (SCALING[ii], BLUES[ii], BLUES[ii])
        if ii == 0:
            red_tuple = (red_val,)
            green_tuple = (green_val,)
            blue_tuple = (blue_val,)
        else:
        	red_tuple = red_tuple + (red_val,)
        	green_tuple = green_tuple + (green_val,)
        	blue_tuple = blue_tuple + (blue_val,)
                
    cdict1 = {'red': red_tuple,
              'green': green_tuple,
              'blue': blue_tuple}
    thismap = LinearSegmentedColormap('manual', cdict1)

    # import test SST file: downloaded from CMEMS
    SST_file = "global-analysis-forecast-phy-001-024_1552636895524.nc"
    ds1 = xr.open_dataset(SST_file)
    LAT = ds1.latitude
    LON = ds1.longitude
    SST_K = np.squeeze(ds1.thetao.values) + 273.15
    ds1.close()
    SST_K = SST_K[:,1500:2500]
    LON = LON[1500:2500]
    plt.pcolormesh(LON, LAT, SST_K, cmap=thismap, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.show()
    return ax1, ax2, REDS, GREENS, BLUES