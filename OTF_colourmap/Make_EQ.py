from ipywidgets import Layout, Box, FloatSlider
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import xarray as xr
import os

def process_data(VAR, LON, LAT, bbox, varname, reduce, log_var=False):
    '''
     Just some pre-processing. Probably requires extensive testing.
    '''
    if len(np.shape(VAR)) > 2:
        VAR = np.squeeze(VAR[0,:,:])

    if 'sst' in varname and np.nanmin(VAR) > 100:
        VAR = VAR - 273.15

    ii = np.where((LON >= bbox[0]) & (LON <= bbox[1]))[0]
    jj = np.where((LAT >= bbox[2]) & (LAT <= bbox[3]))[0]

    LON = LON[ii[0]:ii[-1]]
    LAT = LAT[jj[0]:jj[-1]]
    VAR = VAR[jj[0]:jj[-1],ii[0]:ii[-1]]

    LON = LON[::reduce]
    LAT = LAT[::reduce]
    VAR = VAR[::reduce,::reduce]
    
    if log_var:
        vlimits = [np.log10(np.nanmin(VAR)-0.01), np.log10(np.nanmax(VAR)+0.01), 255]
        VAR = np.log10(VAR)
    else:
        vlimits = [int(np.nanmin(VAR))-1, int(np.nanmax(VAR))+1, 255]

    return VAR, LON, LAT, vlimits

def build_plot_command(vlimits, red_widgets, green_widgets, blue_widgets):
    
    runCMD = 'iplot = interactive(eq.plot_eq'\
                        + ', varname=fixed(varname)' \
                        + ', var=fixed(VAR), lon=fixed(LON), lat=fixed(LAT)' \
                        + ', vmin=fixed(' + str(vlimits[0]) + ')' \
                        + ', vmax=fixed(' + str(vlimits[1]) + ')' \
                        + ', log_var=fixed(log_var), ' \
                        + ', '.join(red_widgets) + ', ' \
                        + ', '.join(green_widgets) + ', ' \
                        + ', '.join(blue_widgets) + ')'
    return runCMD
    
def make_widgets(channel_red, channel_green, channel_blue):
    '''
     Build the interactive plotter widgets
    '''
    red_widgets = [] ; green_widgets = [] ; blue_widgets = []
    for ii in range(len(channel_red)):
        red_widgets.append('red' + str(ii) + '=' + 'widgets.FloatSlider(value=' \
                           + str(channel_red[ii]) + ', min=0.0, max=1.0, step=0.025)')
        green_widgets.append('green' + str(ii) + '=' + 'widgets.FloatSlider(value=' \
                             + str(channel_green[ii]) + ', min=0.0, max=1.0, step=0.025)')
        blue_widgets.append('blue' + str(ii) + '=' + 'widgets.FloatSlider(value=' \
                            + str(channel_blue[ii]) + ', min=0.0, max=1.0, step=0.025)')
    return red_widgets, green_widgets, blue_widgets

def plot_eq(varname, var, lon, lat, vmin, vmax, log_var, **kwargs):

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

    plt.figure(figsize=(10,10))
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
    plt.xlabel('Channel (top: grey) & Mapped value (bottom: black)')
    plt.ylabel('Channel intensity')
    leg1 = plt.legend([pr, pg, pb, pn, pbr],['red','green','blue','norm','r/g sens'])
    
    if log_var:
        ticks = 10**(x)
        newticks = []
        for tick in ticks:
            newticks.append(float(int(tick*100)/100.))
        plt.xticks(x, newticks)
    else:
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

    # plot variable
    plt.pcolormesh(lon, lat, var, cmap=thismap, vmin=vmin, vmax=vmax)
    
    if log_var:
        cbar = plt.colorbar(ticks=x)
        cbar.ax.set_yticklabels(newticks)
    else:
        cbar = plt.colorbar()
      
    cbar.set_label(varname, fontweight='bold', labelpad=10)

    plt.show()
    
    return ax1, ax2, REDS, GREENS, BLUES