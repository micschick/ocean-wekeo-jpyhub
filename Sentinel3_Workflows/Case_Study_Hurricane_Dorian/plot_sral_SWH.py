#!/Users/blo/Anaconda3/python
'''
 Plots SRAL, SLSTR and OLCI data
'''
import os
import sys
import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import gridspec
import matplotlib.patheffects as path_effects
import matplotlib.ticker as mticker
import glob
from skimage import exposure
import grid_tools as gt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def check_versions():
    print('Versions')
    
def run(swh_file, raster_dir, raster_type, fname, label):

    fsz = 32 

    if raster_type == 'OLCI_L2':
        nc_fid = nc.Dataset(os.path.join(raster_dir,'geo_coordinates.nc'),'r')
        raster_lon = nc_fid.variables['longitude'][:]
        raster_lat = nc_fid.variables['latitude'][:]
        nc_fid.close()

        nc_fid = nc.Dataset(os.path.join(raster_dir,'chl_nn.nc'),'r')
        raster_field = nc_fid.variables['CHL_NN'][:]
        nc_fid.close()
    elif raster_type == 'OLCI_L1':
        nc_fid = nc.Dataset(os.path.join(raster_dir,'geo_coordinates.nc'),'r')
        raster_lon = nc_fid.variables['longitude'][:]
        raster_lat = nc_fid.variables['latitude'][:]
        nc_fid.close()

        ch_file = raster_dir+'/Oa01'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa01_reflectance = ch_fid.variables['Oa01_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa02'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa02_reflectance = ch_fid.variables['Oa02_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa03'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa03_reflectance = ch_fid.variables['Oa03_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa04'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa04_reflectance = ch_fid.variables['Oa04_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa05'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa05_reflectance = ch_fid.variables['Oa05_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa06'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa06_reflectance = ch_fid.variables['Oa06_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa07'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa07_reflectance = ch_fid.variables['Oa07_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa08'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa08_reflectance = ch_fid.variables['Oa08_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa09'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa09_reflectance = ch_fid.variables['Oa09_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa10'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa10_reflectance = ch_fid.variables['Oa10_radiance'][:]
        ch_fid.close()

        ch_file = raster_dir+'/Oa11'+'_radiance.nc'
        ch_fid = nc.Dataset(ch_file,'r')
        Oa11_reflectance = ch_fid.variables['Oa11_radiance'][:]
        ch_fid.close()

        red = np.log10(
            0.05 + 0.01 * Oa01_reflectance + 0.09 * Oa02_reflectance + 0.35 * Oa03_reflectance + 0.04 * Oa04_reflectance + 0.01 * Oa05_reflectance + 0.59 * Oa06_reflectance + 0.85 * Oa07_reflectance + 0.12 * Oa08_reflectance + 0.07 * Oa09_reflectance + 0.04 * Oa10_reflectance)
        green = np.log10(
            0.05 + 0.26 * Oa03_reflectance + 0.21 * Oa04_reflectance + 0.50 * Oa05_reflectance + Oa06_reflectance + 0.38 * Oa07_reflectance + 0.04 * Oa08_reflectance + 0.03 * Oa09_reflectance + 0.02 * Oa10_reflectance)
        blue = np.log10(
            0.05 + 0.07 * Oa01_reflectance + 0.28 * Oa02_reflectance + 1.77 * Oa03_reflectance + 0.47 * Oa04_reflectance + 0.16 * Oa05_reflectance)

    else:
        nc_fid = nc.Dataset(os.path.join(raster_dir,'geodetic_in.nc'),'r')
        raster_lon = nc_fid.variables['longitude_in'][:,100:-100]
        raster_lat = nc_fid.variables['latitude_in'][:,100:-100]
        nc_fid.close()

        #nc_fid = nc.Dataset(os.path.join(raster_dir,'S9_BT_in.nc'),'r')
        #red = nc_fid.variables['S9_BT_in'][:]

        #nc_fid = nc.Dataset(os.path.join(raster_dir,'S8_BT_in.nc'),'r')
        #green = nc_fid.variables['S8_BT_in'][:]

        nc_fid = nc.Dataset(os.path.join(raster_dir,'S7_BT_in.nc'),'r')
        raster_field = nc_fid.variables['S7_BT_in'][:,100:-100]

        raster_field = gt.gap_fill(raster_field)

        nc_fid.close()

    if raster_type == 'OLCI_L1':

        unhitch = True
        hist = False
        trunc = True
        Contrast_val = 1.0 # contrast proxy
        Contrast =[Contrast_val, Contrast_val, Contrast_val]
        Brightness = 1.0 # inverse brightness proxy

        if trunc:
            min_red = np.percentile(red,5)
            max_red = np.percentile(red,95)
            min_green = np.percentile(green,5)
            max_green = np.percentile(green,95)
            min_blue = np.percentile(blue,5)
            max_blue = np.percentile(blue,95)

            red[red<min_red] = min_red
            green[green<min_green] = min_green
            blue[blue<min_blue] = min_blue

            red[red>max_red] = max_red
            green[green>max_green] = max_green
            blue[blue>max_blue] = max_blue

        if unhitch:
            # normalise
            red = (red - np.nanmin(red)) / (np.nanmax(red) - np.nanmin(red))
            green = (green - np.nanmin(green)) / (np.nanmax(green) - np.nanmin(green))
            blue = (blue - np.nanmin(blue)) / (np.nanmax(blue) - np.nanmin(blue))

            # non-linearity: contrast
            red = red ** float(Contrast[0])
            green = green ** float(Contrast[1])
            blue = blue ** float(Contrast[2])
        else:
            minval = np.nanmin([np.nanmin(red), np.nanmin(green), np.nanmin(blue)])
            red = red - minval
            green = green - minval
            blue = blue - minval

            # re-scale1: 1 end
            maxval = np.nanmax([np.nanmax(red), np.nanmax(green), np.nanmax(blue)])
            red = abs(red / maxval)
            green = abs(green / maxval)
            blue = abs(blue / maxval)

            # non-linearity
            red = red ** float(Contrast[0])
            green = green ** float(Contrast[1])
            blue = blue ** float(Contrast[2])

            # re-scale1 : 0 end
            minval = np.nanmin([np.nanmin(red), np.nanmin(green), np.nanmin(blue)])
            red = red - minval
            green = green - minval
            blue = blue - minval

            # re-scale1: 1 end
            maxval = np.nanmax([np.nanmax(red), np.nanmax(green), np.nanmax(blue)])
            red = abs(red / maxval)
            green = abs(green / maxval)
            blue = abs(blue / maxval)

        height = np.shape(red)[0]
        width = np.shape(red)[1]
        DATA = np.zeros((height, width, 3), dtype=np.float32)

        DATA[..., 0] = red
        DATA[..., 1] = green
        DATA[..., 2] = blue

        if hist:
            DATA = exposure.equalize_adapthist(DATA, nbins=64)
            
        mesh_rgb = DATA[:, :-1, :]
        colorTuple = mesh_rgb.reshape((mesh_rgb.shape[0] * mesh_rgb.shape[1]), 3)
        colorTuple = np.insert(colorTuple, 3, 1.0, axis=1)

    # init plot
    land_resolution = '10m'
    land_poly = cfeature.NaturalEarthFeature('physical', 'land', land_resolution,
                                        edgecolor='k',
                                        facecolor=cfeature.COLORS['land'])
    lonmin = -83
    lonmax = -69.95
    latmin = 21
    latmax = 31

    # run through to get max and min for normalisation
    SWH_max = 7.218
    SWH_min = 0.0

    # run through to get max and min for normalisation
    SSHA_max = 0.594
    SSHA_min = -0.156

    nc_fid = nc.Dataset(swh_file,'r')
    LON = nc_fid.variables['lon_01'][:]
    LAT = nc_fid.variables['lat_01'][:]
    SWH = nc_fid.variables['swh_ocean_01_ku'][:]
    SSHA = nc_fid.variables['ssha_01_ku'][:]
    FLAGS = nc_fid.variables['swh_ocean_qual_01_ku'][:]
    RMS = nc_fid.variables['swh_ocean_rms_01_ku'][:]
    RAIN = nc_fid.variables['rain_flag_01_ku'][:]
    SURF_TYPE = nc_fid.variables['surf_type_01'][:]

    RANGE = nc_fid.variables['range_ocean_01_ku'][:]
    ALT = nc_fid.variables['alt_01'][:]
    nc_fid.close()


    SSH = ALT - RANGE

    LON[LON>180] = LON[LON>180]-360

    ii = np.where((LON>lonmin) & (LON<lonmax) & (LAT>latmin) & (LAT<latmax))

    LON = LON[ii]
    LAT = LAT[ii]
    SWH = SWH[ii]
    SSHA = SSHA[ii]
    FLAGS = FLAGS[ii]
    RMS = RMS[ii]
    RAIN = RAIN[ii]
    SURF_TYPE = SURF_TYPE[ii]

    # flagging   
    SWH[RAIN != 0] = np.nan
    SWH[SURF_TYPE != 0] = np.nan
    SWH[FLAGS != 0] = np.nan

    SSHA[RAIN != 0] = np.nan
    SSHA[SURF_TYPE != 0] = np.nan
    SSHA[FLAGS != 0] = np.nan
    SSHA[SSHA<-10] = np.nan

    print(np.nanmin(SWH),np.nanmax(SWH))
    print(np.nanmin(SSHA),np.nanmax(SSHA))
    
    fig1 = plt.figure(figsize=(21, 20), dpi=300)
    plt.rc('font', size=fsz)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    gs = gridspec.GridSpec(2, 1, height_ratios=[200,5])
    m = plt.subplot(gs[0,0], projection=ccrs.PlateCarree())

    #norm
    SWH = (SWH-SWH_min)/(SWH_max-SWH_min)
    SSHA = (SSHA-SSHA_min)/(SSHA_max-SSHA_min)

    print(np.nanmin(SWH),np.nanmax(SWH))
    print(np.nanmin(SSHA),np.nanmax(SSHA))
    
    if raster_type == 'OLCI_L2':
        plot1 = m.pcolormesh(raster_lon,raster_lat,raster_field, zorder=0,vmin=-2,vmax=2)
    elif raster_type == 'SLSTR_L1':
        plot1 = m.pcolormesh(raster_lon,raster_lat,raster_field, zorder=0,cmap=plt.cm.bone_r)
    else:
        plot1 = m.pcolormesh(raster_lon,raster_lat, red * np.nan, color=colorTuple ** Brightness, clip_on = True,
                         edgecolors=None, zorder=0, linewidth=0.6, transform=ccrs.PlateCarree())

    lin = 8000
    nlin = 3

    #m.plot(LON, LAT, c='r', linewidth=1.0, zorder=1)
    m.scatter(LON, LAT, s=SWH**nlin*lin, c=SSHA, zorder=10, cmap=plt.cm.YlOrRd, vmin=0, vmax=1, alpha=0.25)

    # key
    m.scatter(-71.6, 22.5, s=(1/0.9624)**nlin*lin, c=plt.cm.YlOrRd(0.6), zorder=2, edgecolor='k')
    m.scatter(-71.6, 23.5, s=(0.80/0.9624)**nlin*lin, c=plt.cm.YlOrRd(0.6), zorder=2, edgecolor='k')
    m.scatter(-71.6, 24.3, s=(0.60/0.9624)**nlin*lin, c=plt.cm.YlOrRd(0.6), zorder=2, edgecolor='k')
    m.scatter(-71.6, 24.9, s=(0.40/0.9624)**nlin*lin, c=plt.cm.YlOrRd(0.6), zorder=2, edgecolor='k')
    m.scatter(-71.6, 25.3, s=(0.20/0.9624)**nlin*lin, c=plt.cm.YlOrRd(0.6), zorder=2, edgecolor='k')

    noff=-0.1

    txt = plt.annotate(str(1.0*SWH_max/0.9624)+' m', xy=(-71.0, 22.5+noff), \
                   xycoords='data', size=fsz/1.5,\
                   color='#8B1A1A', zorder=100, annotation_clip=False)

    txt = plt.annotate(str(0.8*SWH_max/0.9624)+' m', xy=(-71.0, 23.5+noff), \
                   xycoords='data', size=fsz/1.5,\
                   color='#8B1A1A', zorder=100, annotation_clip=False)

    txt = plt.annotate(str(0.6*SWH_max/0.9624)+' m', xy=(-71.0, 24.3+noff), \
                   xycoords='data', size=fsz/1.5,\
                   color='#8B1A1A', zorder=100, annotation_clip=False)

    txt = plt.annotate(str(0.4*SWH_max/0.9624)+' m', xy=(-71.0, 24.9+noff), \
                   xycoords='data', size=fsz/1.5,\
                   color='#8B1A1A', zorder=100, annotation_clip=False)

    txt = plt.annotate(str(0.2*SWH_max/0.9624)+' m', xy=(-71.0, 25.3+noff), \
                   xycoords='data', size=fsz/1.5,\
                   color='#8B1A1A', zorder=100, annotation_clip=False)

    # scale labels

    txt = plt.annotate('SSHA [m]', xy=(-71.6, 30.5), \
                   xycoords='data', size=fsz/1.5,\
                   color='#FFB90F', zorder=100, annotation_clip=False)
    txt.set_path_effects([path_effects.Stroke(linewidth=5, foreground='black'), path_effects.Normal()])

    txt = plt.annotate('SWH [m]', xy=(-71.6, 25.7), \
                   xycoords='data', size=fsz/1.5,\
                   color='#FFB90F', zorder=100, annotation_clip=False)
    txt.set_path_effects([path_effects.Stroke(linewidth=5, foreground='black'), path_effects.Normal()])

    m.coastlines(resolution=land_resolution, color='k', linewidth=0.5, zorder=6)
    m.add_feature(land_poly, color='0.5', zorder=5)

    m.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())
    g1 = m.gridlines(draw_labels = True, zorder=20, color='0.5', linestyle='--',linewidth=0.5)

    g1.xlocator = mticker.FixedLocator([-86,-82,-78,-74,-70])
    g1.ylocator = mticker.FixedLocator([18,22,26,30,34])

    g1.xlabels_top = False
    g1.ylabels_right = False
    g1.xlabel_style = {'size': fsz, 'color': 'black'}
    g1.ylabel_style = {'size': fsz, 'color': 'black'}

    txt = plt.annotate(label, xy=(0.005, 0.005), xycoords='axes fraction', size=fsz/1.5,
                   color='#FFB90F', zorder=100, annotation_clip=False)
    txt.set_path_effects([path_effects.Stroke(linewidth=5, foreground='black'), path_effects.Normal()])

    ax = plt.gca()
    cbaxes = inset_axes(ax, width="100%", height="100%", loc='upper left',
                   bbox_to_anchor=(0.615,1-0.35,.5,.3), bbox_transform=ax.transAxes)

    cbaxes.yaxis.tick_right()

    cbaxes.tick_params(axis='y', colors='#8B1A1A')

    x = np.asarray([np.linspace(SSHA_min,SSHA_max,100),np.linspace(SSHA_min,SSHA_max,100)\
                   ,np.linspace(SSHA_min,SSHA_max,100),np.linspace(SSHA_min,SSHA_max,100)]).astype(float).T
  
    color_map = cbaxes.imshow(x)
    color_map.set_cmap(plt.cm.YlOrRd)

    color_map.axes.get_xaxis().set_ticks([])
    color_map.axes.get_xaxis().set_ticklabels([])

    # need to map 0 to 100 to SSHA_min to SSHA_max
    ticks = []
    tick_labels = []

    for ii in [-0.1,0,0.1,0.2,0.3,0.4,0.5]:
        ticks.append((ii - SSHA_min)/(SSHA_max - SSHA_min)*100)

    for item in ([cbaxes.title, cbaxes.xaxis.label, cbaxes.yaxis.label] +
             cbaxes.get_xticklabels() + cbaxes.get_yticklabels()):
        item.set_fontsize(fsz/1.5)

    color_map.axes.set_yticks(ticks)
    color_map.axes.set_yticklabels([-0.1,0,0.1,0.2,0.3,0.4,0.5])

    if raster_type == 'OLCI_L2':
        plt.savefig(fname, bbox_inches='tight')

        plt.gca().set_visible(False)

        ax = plt.subplot(gs[1,0])
        cbar = plt.colorbar(plot1, cax=ax, orientation='horizontal')
        cbar.set_label('Chlorophyll [mg.m$^{-3}$]', fontsize=fsz)
        ticks = [-2,-1,0,1,2]
        labels=[10**-2,10**-1,10**0,10**1,10**2]
        cbar.set_ticks(ticks)
        cbar.ax.set_xticklabels(labels)
        plt.savefig(fname+'_colbar', bbox_inches='tight')
    else:
        plt.savefig(fname, bbox_inches='tight')
 
if __name__ == '__main__':

    SWH_file = 'standard_measurement_20190902T025631.nc'
    RASTER_dir = 'S3B_SL_1_RBT____20190902T025313_20190902T025613_20190902T035352_0179_029_232_0360_MAR_O_NR_003.SEN3/'
    label='Sentinel-3B SLSTR\nS9 Brightness Temp.\n02/09/2019 (night)\nSRAL L2 SWH Overlaid'
    run(SWH_file, RASTER_dir, 'SLSTR_L1','SWH_map1', label)
    
    SWH_file = 'standard_measurement_20190902T145920.nc'
    RASTER_dir = 'S3B_OL_2_WFR____20190902T151607_20190902T151907_20190904T003825_0179_029_239_2520_MAR_O_NT_002.SEN3/'
    label='Sentinel-3B OLCI\nCHL_NN\n02/09/2019 (day)\nSRAL L2 SWH Overlaid'
    run(SWH_file, RASTER_dir, 'OLCI_L2','SWH_map2', label)

    SWH_file = 'standard_measurement_20190903T030245.nc'
    RASTER_dir = 'S3A_SL_1_RBT____20190903T030633_20190903T030933_20190903T040924_0180_049_004_0360_MAR_O_NR_003.SEN3/'
    label='Sentinel-3A SLSTR\nS9 Brightness Temp.\n03/09/2019 (night)\nSRAL L2 SWH Overlaid'
    run(SWH_file, RASTER_dir, 'SLSTR_L1','SWH_map3', label)

    SWH_file = 'standard_measurement_20190903T152903.nc'
    RASTER_dir = 'S3A_OL_1_EFR____20190903T152927_20190903T153227_20190903T173906_0179_049_011_2520_MAR_O_NR_002.SEN3/'
    label='Sentinel-3A OLCI\nTOA (L1) RGB\n03/09/2019 (day)\nSRAL L2 SWH Overlaid'
    run(SWH_file, RASTER_dir, 'OLCI_L1','SWH_map4', label)
