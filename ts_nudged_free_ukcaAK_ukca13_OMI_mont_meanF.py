# import iris
# import matplotlib
# import sys
# import csv
# import iris.plot as iplt
# import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from iris.analysis.interpolate import extract_nearest_neighbour
# from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset, num2date
# import pandas as pd
# import iris.pandas
# import glob
import datetime
import netCDF4 as nc4
import matplotlib.cm as cm
import math
import scipy.io as sio
from scipy.interpolate import interp2d  as interp2d
# cmap = matplotlib.cm.get_cmap('brewer_RdBu_11')
cmap = cm.jet

Today_date=datetime.datetime.now().strftime("%Y%m%d")

def discrete_cmap(N, base_cmap=None):
	"""Create an N-bin discrete colormap from the specified input map"""
	# Note that if base_cmap is a string or None, you can simply do
	#    return plt.cm.get_cmap(base_cmap, N)
	# The following works for string, None, or a colormap instance:
	base = plt.cm.get_cmap(base_cmap)
	color_list = base(np.linspace(0, 1, N))
	cmap_name = base.name + str(N)
	return base.from_list(cmap_name, color_list, N)

def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# print np.nansum(np.nansum(area,axis=1),axis=0)
	return area

def spatial_figure(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True): #c_bad,c_under,c_over,c_number=20,
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
		lons and lats are 1-d array while data is 2-D array
		colorbar_min,colorbar_max specifies the minimum and maximum value you want to show with the hottest and coldest color respectively
		tb_lef and tb_bot specifies if you want to have axis labels, True fro yes
	output : a spatial map of the data
	"""
	lons[lons>180]-=360; 
	# lon_b = np.min(lons); lon_e = np.max(lons)
	# lon_b = -180; lon_e = 180
	lon_b = 65; lon_e = 100
	# lat_b = np.min(lats); lat_e = np.max(lats)	

	# lat_b = -90; lat_e = 90
	lat_b = 5; lat_e = 38
	# lon_bin = 60; lat_bin = 30
	lon_bin = 10; lat_bin = 10
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)

	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=16)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=16)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines(); map.drawcountries() #map.drawstates(); # draw border lines, here only coast lines
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	#masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(10,colormap)   #  use 20 color bins, this can be changed
	cmap.set_bad([1,1,1],alpha = 1.0);
	#cmap.set_under /cmap.set_over
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	return colormesh

def box_clip(lon_s,lon_e,lat_s,lat_e,lon,lat,mask):
	"""
	fill the region outside the box with 0
	"""
	lon = np.array(lon)
	lat = np.array(lat)
	colum_s = [index for index in range(len(lon)) if np.abs(lon-lon_s)[index] == np.min(np.abs(lon-lon_s))][0]
	colum_e = [index for index in range(len(lon)) if np.abs(lon-lon_e)[index] == np.min(np.abs(lon-lon_e))][0]
	row_s = [index for index in range(len(lat)) if np.abs(lat-lat_s)[index] == np.min(np.abs(lat-lat_s))][0]
	row_e = [index for index in range(len(lat)) if np.abs(lat-lat_e)[index] == np.min(np.abs(lat-lat_e))][0]
	if (colum_s> colum_e):
		cache = colum_e; colum_e = colum_s; colum_s = cache;
	if (row_s> row_e):
		cache = row_e; row_e = row_s; row_s = cache;
	mask[:,0:colum_s] =0; mask[:,colum_e:-1] =0
	# plt.imshow(mask,origin='lower');plt.show()
	mask[0:row_s,:] =0; mask[row_e:-1,:] =0
	# plt.imshow(mask,origin='lower');plt.show()
	return mask	

def mask_weight(region_key,lon,lat,return_option,reverse=False):
	"""
	Read in the country mask
	interpolate it to the required resolution grids with lon_interp,lat_interp 
	crop the sepecified region eithr as a box or an administrtative polygon
	input: 
		region_ky: region name, say, India
		lon and lat of your data
	output: depent on output_option
		if output_option == 'mask': output mask (1 for mask and nan for others)
		elif output_option == 'area': output area of a mask
		elif output_option == 'area_weight': output weight of area against the total area of the mask, this is useful when you do an area-weighted mean
	"""
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	lons,lats = np.meshgrid(lon,lat)
	area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)

	##OCEAN_MASKS FOR COUNTRIES
	# ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_USA_AUS_BRICS_STA_720_360.mat')  ## change this accordingly
	ocean_mask = sio.loadmat('/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/Mat_File/Euro_USA_AUS_BRICS_STA_720_360.mat')  ## change this accordingly
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	## define your regions here
	box_region_dic={'North_India':[75,80,25,30],'South_India':[75,80,15,20],'East_China':[115,130,30,40],'West_China':[85,105,30,40],'All':[0,360,-90,90],'ASIA':[65,145,5,45],'US':[240,290,30,50],'ARCTIC':[0,360,60,90],'TROPICS':[0,360,-28,28],'EUROPE':[0,40,30,70],}
	if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'GloLand'):
		mask= ocean_mask[region_key][:]
	elif  region_key in box_region_dic:
		mask= ocean_mask['All'][:]
		box = box_region_dic[region_key]
		mask = box_clip(box[0],box[1],box[2],box[3],lon_mask,lat_mask,mask)
	else:
		print "error region name"
	
	# interpolate from 360*720 to your grids
	mask[np.isnan(mask)]=0;	mask[mask>0]=1;
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); mask = f(lon, lat);
	# plt.imshow(mask,origin='lower');plt.show()
	mask[mask >= 1] = 1;mask[mask < 1] = 0;
	# weight each grid cell by its area weight against the total area
	if reverse:    ## note this Flase by default, but can be used to exclude the specified  region from a larger map
		mask=1-mask
	mask[mask==0] = np.nan
	grid_area=np.multiply(mask,area); 
	mask_weighted = np.divide(grid_area,np.nansum(np.nansum(grid_area,axis=1),axis=0))
	if return_option == 'mask': return mask
	elif return_option == 'area': return grid_area
	elif return_option == 'area_weight': return mask_weighted

def nudged_run_mon():
	NA=6.022e23   #molecules/mol
	mo3=48.0      #g(o3)/mol
	mNO2=46.0 	  #g(NO2)/mol
	mair=28.97    #g(air)/mol
	DU=2.69e16    #molecules(o3)/cm2
	# Area per grid N96
	areaf = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/Models/gridarea_N96_luke.nc'
	print areaf
	areaf1 = nc4.Dataset(areaf,mode='r')
	lat = areaf1.variables['lat'][0:144]
	lon = areaf1.variables['lon'][:]
	cube_area = areaf1.variables['area'][0:144,:] #in km2

	Alldata = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/11_merge_nudged_monthly.nc'   
	Alldata1 = nc4.Dataset(Alldata,mode='r')
	Alldata_NO2 = Alldata1.variables['UM_m01s34i996_vn1100'][:]
	Alldata_AM = Alldata1.variables['UM_m01s50i063_vn1100'][:]
	Alldata_TropMask = Alldata1.variables['UM_m01s50i062_vn1100'][:]

	print Alldata_NO2.shape, 'Alldata_NO2.shape'
	print Alldata_AM.shape,
	print Alldata_TropMask.shape,'Alldata_TropMask.shape' 

	am_NO2 = Alldata_AM*Alldata_NO2 #Air Mass * NO2
	am_NO2_tropmask =  am_NO2*Alldata_TropMask

	am_NO2_tropmask[am_NO2_tropmask<=0] = 'NaN'
	print am_NO2_tropmask.shape,'am_NO2_tropmask.shape' 

	am_NO2_tropmaskA = np.nansum(am_NO2_tropmask, axis=1)
	print am_NO2_tropmaskA.shape, 'am_NO2_tropmaskA.shape'

	# am_NO2_tropmaskA = am_NO2_tropmaskA[0,:,:]
	am_NO2_tropmaskB = am_NO2_tropmaskA/(cube_area*1.0e10) #kg(NO2)/cm2
	print am_NO2_tropmaskB.shape, 'am_NO2_tropmaskB.shape'

	am_NO2_tropmaskC = am_NO2_tropmaskB/(mNO2*1.0e-3) #moles(NO2)/cm2
	print am_NO2_tropmaskC.shape, 'am_NO2_tropmaskC.shape'

	am_NO2_tropmaskD =am_NO2_tropmaskC*NA   # molecules/cm2
	am_NO2_tropmaskE = am_NO2_tropmaskD[0:-1]

	return lat, lon, am_NO2_tropmaskE

def Free_run_mon():
	NA=6.022e23   #molecules/mol
	mo3=48.0      #g(o3)/mol
	mNO2=46.0 	  #g(NO2)/mol
	mair=28.97    #g(air)/mol
	DU=2.69e16    #molecules(o3)/cm2
	# Area per grid N96
	areaf = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/Models/gridarea_N96_luke.nc'
	print areaf
	areaf1 = nc4.Dataset(areaf,mode='r')
	lat_f = areaf1.variables['lat'][0:144]
	lon_f = areaf1.variables['lon'][:]
	cube_area = areaf1.variables['area'][0:144,:] #in km2

	freerun = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/5_merged_freerun_monthly.nc'   
	freerun1 = nc4.Dataset(freerun,mode='r')
	freerun_NO2 = freerun1.variables['UM_m01s34i996_vn1100'][:]
	freerun_AM = freerun1.variables['UM_m01s50i063_vn1100'][:]
	freerun_TropMask = freerun1.variables['UM_m01s50i062_vn1100'][:]

	print freerun_NO2.shape, 'freerun_NO2.shape'
	print freerun_AM.shape,
	print freerun_TropMask.shape,'freerun_TropMask.shape' 

	free_am_NO2 = freerun_AM*freerun_NO2 #Air Mass * NO2
	free_am_NO2_tropmask =  free_am_NO2*freerun_TropMask

	free_am_NO2_tropmask[free_am_NO2_tropmask<=0] = 'NaN'
	print free_am_NO2_tropmask.shape,'free_am_NO2_tropmask.shape' 

	free_am_NO2_tropmaskA = np.nansum(free_am_NO2_tropmask, axis=1)
	print free_am_NO2_tropmaskA.shape, 'free_am_NO2_tropmaskA.shape'

	free_am_NO2_tropmaskB = free_am_NO2_tropmaskA/(cube_area*1.0e10) #kg(NO2)/cm2
	print free_am_NO2_tropmaskB.shape, 'free_am_NO2_tropmaskB.shape'

	free_am_NO2_tropmaskC = free_am_NO2_tropmaskB/(mNO2*1.0e-3) #moles(NO2)/cm2
	print free_am_NO2_tropmaskC.shape, 'free_am_NO2_tropmaskC.shape'

	free_am_NO2_tropmaskD =free_am_NO2_tropmaskC*NA   # molecules/cm2
	free_am_NO2_tropmaskE = free_am_NO2_tropmaskD[0:-1]

	
	return lat_f, lon_f, free_am_NO2_tropmaskE

def hour1330_mon():
	hour1330_data = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/11_hour1330_month2005_2015.nc'
	hour1330 = nc4.Dataset(hour1330_data,mode='r')
	lat_h = hour1330.variables['lat'][:]
	lon_h = hour1330.variables['lon'][:]
	hour1330_no2 = hour1330.variables['ukca_no2'][:]

	print hour1330_no2.shape, 'hour1330_no2.shape', np.nanmax(hour1330_no2), np.nanmean(hour1330_no2)	
	
	return lat_h, lon_h, hour1330_no2

def ukca_AK_OMI_mon():
	def netcdf_read(year):
		file_name = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_'+str(year)+'_real_comparison.nc'
		ncf = nc4.Dataset(file_name,mode='r')
		lat_AK = ncf.variables['lat'][:]
		lon_AK = ncf.variables['lon'][:]
		NO2_OMI = ncf.variables['omi_no2'][:]
		NO2_UKCA_AK = ncf.variables['ukca_no2'][:]
		NO2_OMI[NO2_OMI <= 0.1e15] = np.nan
		NO2_UKCA_AK[NO2_UKCA_AK <= 0.1e15] = np.nan

		return lon_AK,lat_AK,NO2_OMI,NO2_UKCA_AK
	
		
	def monthly_mean(year,data):
		if year in [2008,2012,2016]:
			c_days = np.array([0,31,60,91,121,152,182,213,244,274,305,335,366])
		else:
			c_days = np.array([0,31,59,90,120,151,181,212,243,273,304,334,365])
		
		mon_mean = np.empty((12,41,50))
		mon_mean[:]=np.nan
		
		for i in range(0,12):
			mon_mean[i,:,:]=np.nanmean(data[c_days[i]:c_days[i+1],:,:],axis=0)
		return mon_mean
		
		
	start_year = 2005; 
	end_year   = 2016;
	years = range(start_year,end_year);
	mon_mean_UKCA  = np.empty((12*(end_year-start_year),41,50));
	mon_mean_UKCA[:]= np.nan
	mon_mean_OMI  = np.empty((12*(end_year-start_year),41,50));
	mon_mean_OMI[:]= np.nan
	for iyear in years:
		print iyear
		lon_AK,lat_AK,NO2_OMI,NO2_UKCA_AK = netcdf_read(iyear);
		mon_mean_OMI[(iyear-start_year)*12:(iyear-start_year+1)*12,:,:]  = monthly_mean(iyear,NO2_OMI)
		mon_mean_UKCA[(iyear-start_year)*12:(iyear-start_year+1)*12,:,:] = monthly_mean(iyear,NO2_UKCA_AK)
	print mon_mean_OMI.shape, 'mon_mean_OMI.shape'
	return lon_AK, lat_AK, mon_mean_OMI, mon_mean_UKCA

title_list  = 'OMI UKCA NO2 comparison'

lat, lon, mon_mean_nudged = nudged_run_mon()
lat_f, lon_f, mon_mean_freerun = Free_run_mon()
lat_h, lon_h, hour1330_no2 = hour1330_mon()
lon_AK, lat_AK, mon_mean_OMI, mon_mean_UKCA_ak = ukca_AK_OMI_mon()

mon_mean_OMI = mon_mean_OMI[0:-1]
mon_mean_UKCA_ak = mon_mean_UKCA_ak[0:-1]

def india_ts_mon():
	region_key='India';return_option='area_weight'
	mask_IndiaA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	mask_IndiaB=mask_weight(region_key,lon_f.copy(),lat_f.copy(),return_option,reverse=False);
	mask_IndiaC=mask_weight(region_key,lon_h.copy(),lat_h.copy(),return_option,reverse=False);
	mask_IndiaD=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);

	## weighted average of your data
	India_mean_UKCA_nudged 		= np.nansum(np.nansum(np.multiply(mask_IndiaA,mon_mean_nudged),axis=2),axis=1)
	India_mean_UKCA_freerun 	= np.nansum(np.nansum(np.multiply(mask_IndiaB,mon_mean_freerun),axis=2),axis=1)
	India_mean_UKCA_hour1330 	= np.nansum(np.nansum(np.multiply(mask_IndiaC,hour1330_no2),axis=2),axis=1)
	India_mean_UKCA_ak          = np.nansum(np.nansum(np.multiply(mask_IndiaD,mon_mean_UKCA_ak),axis=2),axis=1)
	India_mean_OMI              = np.nansum(np.nansum(np.multiply(mask_IndiaD,mon_mean_OMI),axis=2),axis=1)
	return India_mean_UKCA_nudged, India_mean_UKCA_freerun, India_mean_UKCA_hour1330, India_mean_UKCA_ak, India_mean_OMI

def china_ts_mon():
	region_key='China';return_option='area_weight'
	mask_ChinaA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	mask_ChinaB=mask_weight(region_key,lon_f.copy(),lat_f.copy(),return_option,reverse=False);
	mask_ChinaC=mask_weight(region_key,lon_h.copy(),lat_h.copy(),return_option,reverse=False);
	mask_ChinaD=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);

	## weighted average of your data
	China_mean_UKCA_nudged 		= np.nansum(np.nansum(np.multiply(mask_ChinaA,mon_mean_nudged),axis=2),axis=1)
	China_mean_UKCA_freerun 	= np.nansum(np.nansum(np.multiply(mask_ChinaB,mon_mean_freerun),axis=2),axis=1)
	China_mean_UKCA_hour1330 	= np.nansum(np.nansum(np.multiply(mask_ChinaC,hour1330_no2),axis=2),axis=1)
	China_mean_UKCA_ak          = np.nansum(np.nansum(np.multiply(mask_ChinaD,mon_mean_UKCA_ak),axis=2),axis=1)
	China_mean_OMI              = np.nansum(np.nansum(np.multiply(mask_ChinaD,mon_mean_OMI),axis=2),axis=1)
	return China_mean_UKCA_nudged, China_mean_UKCA_freerun, China_mean_UKCA_hour1330, China_mean_UKCA_ak, China_mean_OMI
	
def north_india_ts_mon():
	region_key='North_India';return_option='area_weight'
	mask_North_IndiaA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	mask_North_IndiaB=mask_weight(region_key,lon_f.copy(),lat_f.copy(),return_option,reverse=False);
	mask_North_IndiaC=mask_weight(region_key,lon_h.copy(),lat_h.copy(),return_option,reverse=False);
	mask_North_IndiaD=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);

	## weighted average of your data
	North_India_mean_UKCA_nudged 		= np.nansum(np.nansum(np.multiply(mask_North_IndiaA,mon_mean_nudged),axis=2),axis=1)
	North_India_mean_UKCA_freerun 	= np.nansum(np.nansum(np.multiply(mask_North_IndiaB,mon_mean_freerun),axis=2),axis=1)
	North_India_mean_UKCA_hour1330 	= np.nansum(np.nansum(np.multiply(mask_North_IndiaC,hour1330_no2),axis=2),axis=1)
	North_India_mean_UKCA_ak          = np.nansum(np.nansum(np.multiply(mask_North_IndiaD,mon_mean_UKCA_ak),axis=2),axis=1)
	North_India_mean_OMI              = np.nansum(np.nansum(np.multiply(mask_North_IndiaD,mon_mean_OMI),axis=2),axis=1)
	return North_India_mean_UKCA_nudged, North_India_mean_UKCA_freerun, North_India_mean_UKCA_hour1330, North_India_mean_UKCA_ak, North_India_mean_OMI

def south_india_ts_mon():
	region_key='South_India';return_option='area_weight'
	mask_South_IndiaA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	mask_South_IndiaB=mask_weight(region_key,lon_f.copy(),lat_f.copy(),return_option,reverse=False);
	mask_South_IndiaC=mask_weight(region_key,lon_h.copy(),lat_h.copy(),return_option,reverse=False);
	mask_South_IndiaD=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);

	## weighted average of your data
	South_India_mean_UKCA_nudged 		= np.nansum(np.nansum(np.multiply(mask_South_IndiaA,mon_mean_nudged),axis=2),axis=1)
	South_India_mean_UKCA_freerun 	= np.nansum(np.nansum(np.multiply(mask_South_IndiaB,mon_mean_freerun),axis=2),axis=1)
	South_India_mean_UKCA_hour1330 	= np.nansum(np.nansum(np.multiply(mask_South_IndiaC,hour1330_no2),axis=2),axis=1)
	South_India_mean_UKCA_ak          = np.nansum(np.nansum(np.multiply(mask_South_IndiaD,mon_mean_UKCA_ak),axis=2),axis=1)
	South_India_mean_OMI              = np.nansum(np.nansum(np.multiply(mask_South_IndiaD,mon_mean_OMI),axis=2),axis=1)
	return South_India_mean_UKCA_nudged, South_India_mean_UKCA_freerun, South_India_mean_UKCA_hour1330, South_India_mean_UKCA_ak, South_India_mean_OMI

def east_china_ts_mon():
	region_key='East_China';return_option='area_weight'
	mask_East_ChinaA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	mask_East_ChinaB=mask_weight(region_key,lon_f.copy(),lat_f.copy(),return_option,reverse=False);
	mask_East_ChinaC=mask_weight(region_key,lon_h.copy(),lat_h.copy(),return_option,reverse=False);
	mask_East_ChinaD=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);

	## weighted average of your data
	East_China_mean_UKCA_nudged 		= np.nansum(np.nansum(np.multiply(mask_East_ChinaA,mon_mean_nudged),axis=2),axis=1)
	East_China_mean_UKCA_freerun 	= np.nansum(np.nansum(np.multiply(mask_East_ChinaB,mon_mean_freerun),axis=2),axis=1)
	East_China_mean_UKCA_hour1330 	= np.nansum(np.nansum(np.multiply(mask_East_ChinaC,hour1330_no2),axis=2),axis=1)
	East_China_mean_UKCA_ak          = np.nansum(np.nansum(np.multiply(mask_East_ChinaD,mon_mean_UKCA_ak),axis=2),axis=1)
	East_China_mean_OMI              = np.nansum(np.nansum(np.multiply(mask_East_ChinaD,mon_mean_OMI),axis=2),axis=1)
	return East_China_mean_UKCA_nudged, East_China_mean_UKCA_freerun, East_China_mean_UKCA_hour1330, East_China_mean_UKCA_ak, East_China_mean_OMI
	
def west_china_ts_mon():
	region_key='West_China';return_option='area_weight'
	mask_West_ChinaA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	mask_West_ChinaB=mask_weight(region_key,lon_f.copy(),lat_f.copy(),return_option,reverse=False);
	mask_West_ChinaC=mask_weight(region_key,lon_h.copy(),lat_h.copy(),return_option,reverse=False);
	mask_West_ChinaD=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);

	## weighted average of your data
	West_China_mean_UKCA_nudged 		= np.nansum(np.nansum(np.multiply(mask_West_ChinaA,mon_mean_nudged),axis=2),axis=1)
	West_China_mean_UKCA_freerun 	= np.nansum(np.nansum(np.multiply(mask_West_ChinaB,mon_mean_freerun),axis=2),axis=1)
	West_China_mean_UKCA_hour1330 	= np.nansum(np.nansum(np.multiply(mask_West_ChinaC,hour1330_no2),axis=2),axis=1)
	West_China_mean_UKCA_ak          = np.nansum(np.nansum(np.multiply(mask_West_ChinaD,mon_mean_UKCA_ak),axis=2),axis=1)
	West_China_mean_OMI              = np.nansum(np.nansum(np.multiply(mask_West_ChinaD,mon_mean_OMI),axis=2),axis=1)
	return West_China_mean_UKCA_nudged, West_China_mean_UKCA_freerun, West_China_mean_UKCA_hour1330, West_China_mean_UKCA_ak, West_China_mean_OMI

India_mean_UKCA_nudged, India_mean_UKCA_freerun, India_mean_UKCA_hour1330, India_mean_UKCA_ak, India_mean_OMI = india_ts_mon()
China_mean_UKCA_nudged, China_mean_UKCA_freerun, China_mean_UKCA_hour1330, China_mean_UKCA_ak, China_mean_OMI = china_ts_mon()
North_India_mean_UKCA_nudged, North_India_mean_UKCA_freerun, North_India_mean_UKCA_hour1330, North_India_mean_UKCA_ak, North_India_mean_OMI = north_india_ts_mon()
South_India_mean_UKCA_nudged, South_India_mean_UKCA_freerun, South_India_mean_UKCA_hour1330, South_India_mean_UKCA_ak, South_India_mean_OMI = south_india_ts_mon()
East_China_mean_UKCA_nudged, East_China_mean_UKCA_freerun, East_China_mean_UKCA_hour1330, East_China_mean_UKCA_ak, East_China_mean_OMI = east_china_ts_mon()
West_China_mean_UKCA_nudged, West_China_mean_UKCA_freerun, West_China_mean_UKCA_hour1330, West_China_mean_UKCA_ak, West_China_mean_OMI = west_china_ts_mon()	
	
X1 = range(0,215)
X2 = range(84,215)

#plot
fig = plt.figure(facecolor='White',figsize=[60,35]);pad= 0.5; 
# plt.suptitle(title_list, fontsize = 25, y=0.95)

ax = fig.add_subplot(321)
ax.plot(X1,India_mean_UKCA_nudged,c='b', linewidth=2);
ax.plot(X1,India_mean_UKCA_freerun,c='r',linewidth=2);
ax.plot(X2,India_mean_UKCA_hour1330,c='c',linewidth=2);
ax.plot(X2,India_mean_UKCA_ak,c='m',linestyle='-.',linewidth=3);
ax.plot(X2,India_mean_OMI,c='g',linestyle='-.',linewidth=3);

tick_locs = [0,12,24,36,48,60,72,84,96,108,120,132,144,156,168,180,192,204,215]
tick_lbls = ['Jan-98','Jan-99', 'Jan-00', 'Jan-01','Jan-02', 'Jan-03','Jan-04', 'Jan-05','Jan-06', 'Jan-07','Jan-08', 'Jan-09','Jan-10', 'Jan-11','Jan-12','Jan-13', 'Jan-14','Jan-15','Dec-15']
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
plt.ylim(1e14, 2e16)
plt.yticks(fontsize = 25)
ax.yaxis.get_offset_text().set_size(25)

# plt.xlabel("Year 1998-2015", fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)

plt.legend(('India_UKCA_nudged' ,'India_UKCA_freerun','India_UKCA_hour1330','India_UKCA_ak', 'India_sat_OMI'),fontsize=30,loc='upper left')

ax = fig.add_subplot(322)
ax.plot(X1,China_mean_UKCA_nudged,c='b', linewidth=2);
ax.plot(X1,China_mean_UKCA_freerun,c='r',linewidth=2);
ax.plot(X2,China_mean_UKCA_hour1330,c='c',linewidth=2);
ax.plot(X2,China_mean_UKCA_ak,c='m',linestyle='-.',linewidth=3);
ax.plot(X2,China_mean_OMI,c='g',linestyle='-.',linewidth=3);
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
plt.ylim(1e14, 2e16)
plt.yticks(fontsize = 25)
ax.yaxis.get_offset_text().set_size(25)
# plt.xlabel("Year 1998-2015", fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)
plt.legend(('China_UKCA_nudged' ,'China_UKCA_freerun','China_UKCA_hour1330','China_UKCA_ak', 'China_sat_OMI'),fontsize=30,loc='upper left')

ax = fig.add_subplot(323)
ax.plot(X1,North_India_mean_UKCA_nudged,c='b', linewidth=2);
ax.plot(X1,North_India_mean_UKCA_freerun,c='r',linewidth=2);
ax.plot(X2,North_India_mean_UKCA_hour1330,c='c',linewidth=2);
ax.plot(X2,North_India_mean_UKCA_ak,c='m',linestyle='-.',linewidth=3);
ax.plot(X2,North_India_mean_OMI,c='g',linestyle='-.',linewidth=3);
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
plt.ylim(1e14, 5e16)
plt.yticks(fontsize = 25)
ax.yaxis.get_offset_text().set_size(25)
# plt.xlabel("Year 1998-2015", fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)
plt.legend(('North_India_UKCA_nudged' ,'North_India_UKCA_freerun','North_India_UKCA_hour1330','North_India_UKCA_ak', 'North_India_sat_OMI'),fontsize=30,loc='upper left')

ax = fig.add_subplot(324)
ax.plot(X1,East_China_mean_UKCA_nudged,c='b', linewidth=2);
ax.plot(X1,East_China_mean_UKCA_freerun,c='r',linewidth=2);
ax.plot(X2,East_China_mean_UKCA_hour1330,c='c',linewidth=2);
ax.plot(X2,East_China_mean_UKCA_ak,c='m',linestyle='-.',linewidth=3);
ax.plot(X2,East_China_mean_OMI,c='g',linestyle='-.',linewidth=3);
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
plt.ylim(1e14, 5e16)
plt.yticks(fontsize = 25)
ax.yaxis.get_offset_text().set_size(25)
# plt.xlabel("Year 1998-2015", fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)
plt.legend(('East_China_UKCA_nudged' ,'East_China_UKCA_freerun','East_China_UKCA_hour1330','East_China_UKCA_ak', 'East_China_sat_OMI'),fontsize=30,loc='upper left')

ax = fig.add_subplot(325)
ax.plot(X1,South_India_mean_UKCA_nudged,c='b', linewidth=2);
ax.plot(X1,South_India_mean_UKCA_freerun,c='r',linewidth=2);
ax.plot(X2,South_India_mean_UKCA_hour1330,c='c',linewidth=2);
ax.plot(X2,South_India_mean_UKCA_ak,c='m',linestyle='-.',linewidth=3);
ax.plot(X2,South_India_mean_OMI,c='g',linestyle='-.',linewidth=3);
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
plt.ylim(1e14, 6e15)
plt.yticks(fontsize = 25)
ax.yaxis.get_offset_text().set_size(25)
# plt.xlabel("Year 1998-2015", fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)
plt.legend(('South_India_UKCA_nudged' ,'South_India_UKCA_freerun','South_India_UKCA_hour1330','South_India_UKCA_ak', 'South_India_sat_OMI'),fontsize=30,loc='upper left')

ax = fig.add_subplot(326)
ax.plot(X1,West_China_mean_UKCA_nudged,c='b', linewidth=2);
ax.plot(X1,West_China_mean_UKCA_freerun,c='r',linewidth=2);
ax.plot(X2,West_China_mean_UKCA_hour1330,c='c',linewidth=2);
ax.plot(X2,West_China_mean_UKCA_ak,c='m',linestyle='-.',linewidth=3);
ax.plot(X2,West_China_mean_OMI,c='g',linestyle='-.',linewidth=3);
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
plt.ylim(1e14, 6e15)
plt.yticks(fontsize = 25)
ax.yaxis.get_offset_text().set_size(25)
# plt.xlabel("Year 1998-2015", fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)
plt.legend(('West_China_UKCA_nudged' ,'West_China_UKCA_freerun','West_China_UKCA_hour1330','West_China_UKCA_ak', 'West_China_sat_OMI'),fontsize=30,loc='upper left')

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.11, hspace=0.11);

plt.savefig('/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/ESM2019/'+Today_date+'ts_mon_1998-2015_free_nud_hour_modak_omi1.png', dpi=100)   ##########	
plt.show()