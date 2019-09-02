import netCDF4 as nc4
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import cartopy.feature as cfeature
# from mpl_toolkits.basemap import Basemap, shiftgrid
import cartopy.crs as ccrs
import matplotlib.cm as cm
cmap = cm.jet
import datetime
import math
import scipy.io as sio
from scipy.interpolate import interp2d  as interp2d

# from datetime import datetime
# import cf, cfplot as cfp

Today_date=datetime.datetime.now().strftime("%Y%m%d")

def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# print np.nansum(np.nansum(area,axis=1),axis=0)
	return area
	
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
	
# Data after applying AK  
data_AK2005  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2005_real_comparison.nc'
data_AK2006  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2006_real_comparison.nc'
data_AK2007  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2007_real_comparison.nc'
data_AK2008  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2008_real_comparison.nc'
data_AK2009  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2009_real_comparison.nc'
data_AK2010  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2010_real_comparison.nc'
data_AK2011  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2011_real_comparison.nc'
data_AK2012  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2012_real_comparison.nc'
data_AK2013  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2013_real_comparison.nc'
data_AK2014  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2014_real_comparison.nc'
data_AK2015  = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/NO2_OMI_UKCA_Comparison/ukca_omi_no2omi_ukca_2015_real_comparison.nc'
# print data_AK2006  

ncfile2005 = nc4.Dataset(data_AK2005,mode='r')
lat_AK = ncfile2005.variables['lat'][:]
lon_AK = ncfile2005.variables['lon'][:]
NO2_OMI2005 = ncfile2005.variables['omi_no2'][:]
NO2_UKCA_AK2005 = ncfile2005.variables['ukca_no2'][:]

ncfile2006 = nc4.Dataset(data_AK2006,mode='r')
NO2_OMI2006 = ncfile2006.variables['omi_no2'][:]
NO2_UKCA_AK2006 = ncfile2006.variables['ukca_no2'][:]

ncfile2007 = nc4.Dataset(data_AK2007,mode='r')
NO2_OMI2007 = ncfile2007.variables['omi_no2'][:]
NO2_UKCA_AK2007 = ncfile2007.variables['ukca_no2'][:]

ncfile2008 = nc4.Dataset(data_AK2008,mode='r')
NO2_OMI2008 = ncfile2008.variables['omi_no2'][:]
NO2_UKCA_AK2008 = ncfile2008.variables['ukca_no2'][:]

ncfile2009 = nc4.Dataset(data_AK2009,mode='r')
NO2_OMI2009 = ncfile2009.variables['omi_no2'][:]
NO2_UKCA_AK2009 = ncfile2009.variables['ukca_no2'][:]

ncfile2010 = nc4.Dataset(data_AK2010,mode='r')
NO2_OMI2010 = ncfile2010.variables['omi_no2'][:]
NO2_UKCA_AK2010 = ncfile2010.variables['ukca_no2'][:]

ncfile2011 = nc4.Dataset(data_AK2011,mode='r')
NO2_OMI2011 = ncfile2011.variables['omi_no2'][:]
NO2_UKCA_AK2011 = ncfile2011.variables['ukca_no2'][:]

ncfile2012 = nc4.Dataset(data_AK2012,mode='r')
NO2_OMI2012 = ncfile2012.variables['omi_no2'][:]
NO2_UKCA_AK2012 = ncfile2012.variables['ukca_no2'][:]

ncfile2013 = nc4.Dataset(data_AK2013,mode='r')
NO2_OMI2013 = ncfile2013.variables['omi_no2'][:]
NO2_UKCA_AK2013 = ncfile2013.variables['ukca_no2'][:]

ncfile2014 = nc4.Dataset(data_AK2014,mode='r')
NO2_OMI2014 = ncfile2014.variables['omi_no2'][:]
NO2_UKCA_AK2014 = ncfile2014.variables['ukca_no2'][:]

ncfile2015 = nc4.Dataset(data_AK2015,mode='r')
NO2_OMI2015 = ncfile2015.variables['omi_no2'][:]
NO2_UKCA_AK2015 = ncfile2015.variables['ukca_no2'][:]

NO2_OMI = np.concatenate((NO2_OMI2005, NO2_OMI2006, NO2_OMI2007, NO2_OMI2008, NO2_OMI2009, NO2_OMI2010, NO2_OMI2011, NO2_OMI2012, NO2_OMI2013, NO2_OMI2014, NO2_OMI2015), axis=0)
NO2_UKCA_AK = np.concatenate((NO2_UKCA_AK2005, NO2_UKCA_AK2006, NO2_UKCA_AK2007, NO2_UKCA_AK2008, NO2_UKCA_AK2009, NO2_UKCA_AK2010, NO2_UKCA_AK2011, NO2_UKCA_AK2012, NO2_UKCA_AK2013, NO2_UKCA_AK2014, NO2_UKCA_AK2015), axis=0)

print NO2_OMI2005.shape, 'NO2_OMI2006.shape'
print NO2_OMI.shape, 'NO2_OMI.shape'
print NO2_UKCA_AK.shape, 'NO2_UKCA_AK.shape'

NO2_OMI[NO2_OMI <= 0.1e15] = np.nan
NO2_UKCA_AK[NO2_UKCA_AK <= 0.1e15] = np.nan

# NO2_OMI[NO2_OMI > 5e16] = np.nan
# NO2_UKCA_AK[NO2_UKCA_AK > 5e16] = np.nan

print NO2_OMI.shape, np.nanmax(NO2_OMI), np.nanmin(NO2_OMI)

def india_ts():
	region_key='India';return_option='area_weight'
	mask_India=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);
	## weighted average of your data
	India_mean_sat= np.nansum(np.nansum(np.multiply(mask_India,NO2_OMI),axis=2),axis=1)
	India_mean_mod= np.nansum(np.nansum(np.multiply(mask_India,NO2_UKCA_AK),axis=2),axis=1)
	India_mean_sat[India_mean_sat <=0.5e15]= np.nan 
	India_mean_mod[India_mean_mod <=0.5e15]= np.nan 
	print India_mean_sat.shape, 'India_mean_sat.shape'
	print np.nanmax(India_mean_sat),np.nanmin(India_mean_sat)
	print np.nanmax(India_mean_mod),np.nanmin(India_mean_mod)
	return India_mean_sat, India_mean_mod
	
def china_ts():
	region_key='China';return_option='area_weight'
	mask_China=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);
	## weighted average of your data
	China_mean_sat= np.nansum(np.nansum(np.multiply(mask_China,NO2_OMI),axis=2),axis=1)
	China_mean_mod= np.nansum(np.nansum(np.multiply(mask_China,NO2_UKCA_AK),axis=2),axis=1)
	China_mean_sat[China_mean_sat <=0.5e15]= np.nan 
	China_mean_mod[China_mean_mod <=0.5e15]= np.nan 
	print China_mean_sat.shape, 'China_mean_sat.shape'
	print np.nanmax(China_mean_sat),np.nanmin(China_mean_sat)
	print np.nanmax(China_mean_mod),np.nanmin(China_mean_mod)
	return China_mean_sat, China_mean_mod

def north_india_ts():
	region_key='North_India';return_option='area_weight'
	mask_North_India=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);
	## weighted average of your data
	North_India_mean_sat= np.nansum(np.nansum(np.multiply(mask_North_India,NO2_OMI),axis=2),axis=1)
	North_India_mean_mod= np.nansum(np.nansum(np.multiply(mask_North_India,NO2_UKCA_AK),axis=2),axis=1)
	North_India_mean_sat[North_India_mean_sat <=0.5e15]= np.nan 
	North_India_mean_mod[North_India_mean_mod <=0.5e15]= np.nan 
	print North_India_mean_sat.shape, 'North_India_mean_sat.shape'
	print np.nanmax(North_India_mean_sat),np.nanmin(North_India_mean_sat)
	print np.nanmax(North_India_mean_mod),np.nanmin(North_India_mean_mod)
	return North_India_mean_sat, North_India_mean_mod
	
def south_india_ts():
	region_key='South_India';return_option='area_weight'
	mask_South_India=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);
	## weighted average of your data
	South_India_mean_sat= np.nansum(np.nansum(np.multiply(mask_South_India,NO2_OMI),axis=2),axis=1)
	South_India_mean_mod= np.nansum(np.nansum(np.multiply(mask_South_India,NO2_UKCA_AK),axis=2),axis=1)
	South_India_mean_sat[South_India_mean_sat <=0.5e15]= np.nan 
	South_India_mean_mod[South_India_mean_mod <=0.5e15]= np.nan 
	print South_India_mean_sat.shape, 'South_India_mean_sat.shape'
	print np.nanmax(South_India_mean_sat),np.nanmin(South_India_mean_sat)
	print np.nanmax(South_India_mean_mod),np.nanmin(South_India_mean_mod)
	return South_India_mean_sat, South_India_mean_mod

def east_china_ts():
	region_key='East_China';return_option='area_weight'
	mask_East_China=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);
	## weighted average of your data
	East_China_mean_sat= np.nansum(np.nansum(np.multiply(mask_East_China,NO2_OMI),axis=2),axis=1)
	East_China_mean_mod= np.nansum(np.nansum(np.multiply(mask_East_China,NO2_UKCA_AK),axis=2),axis=1)
	East_China_mean_sat[East_China_mean_sat <=0.5e15]= np.nan 
	East_China_mean_mod[East_China_mean_mod <=0.5e15]= np.nan 
	print East_China_mean_sat.shape, 'East_China_mean_sat.shape'
	print np.nanmax(East_China_mean_sat),np.nanmin(East_China_mean_sat)
	print np.nanmax(East_China_mean_mod),np.nanmin(East_China_mean_mod)
	return East_China_mean_sat, East_China_mean_mod
	
def west_china_ts():
	region_key='West_China';return_option='area_weight'
	mask_West_China=mask_weight(region_key,lon_AK.copy(),lat_AK.copy(),return_option,reverse=False);
	## weighted average of your data
	West_China_mean_sat= np.nansum(np.nansum(np.multiply(mask_West_China,NO2_OMI),axis=2),axis=1)
	West_China_mean_mod= np.nansum(np.nansum(np.multiply(mask_West_China,NO2_UKCA_AK),axis=2),axis=1)
	West_China_mean_sat[West_China_mean_sat <=0.5e15]= np.nan 
	West_China_mean_mod[West_China_mean_mod <=0.5e15]= np.nan 
	print West_China_mean_sat.shape, 'West_China_mean_sat.shape'
	print np.nanmax(West_China_mean_sat),np.nanmin(West_China_mean_sat)
	print np.nanmax(West_China_mean_mod),np.nanmin(West_China_mean_mod)
	return West_China_mean_sat, West_China_mean_mod
	
India_mean_sat, India_mean_mod = india_ts()
China_mean_sat, China_mean_mod = china_ts()
North_India_mean_sat, North_India_mean_mod = north_india_ts()
South_India_mean_sat, South_India_mean_mod = south_india_ts()
East_China_mean_sat, East_China_mean_mod = east_china_ts()
West_China_mean_sat, West_China_mean_mod = west_china_ts()

# #plot
# fig = plt.figure(facecolor='White',figsize=[41,25]);pad= 1; 
# #plt.suptitle(title_list, fontsize = 25, y=0.95)

# ax1 = fig.add_subplot(321)
# ax1.plot(India_mean_sat,c='b',linewidth=2);
# ax1.plot(India_mean_mod,c='m',linestyle='-.',linewidth=2);
# # tick_locs =np.linspace(0, 4380, num=12, endpoint=False)
# xtick_locs = [0,  366,  731, 1096, 1462, 1827, 2192, 2557, 2923, 3283, 3653, 4015]
# xtick_lbls = ['Jan-05', 'Jan-06', 'Jan-07','Jan-08', 'Jan-09','Jan-10', 'Jan-11','Jan-12','Jan-13', 'Jan-14', 'Jan-15','Dec-15']
# ytick_locs = np.linspace(1.25e15, 1.5e16, num=8)

# plt.xticks(xtick_locs, xtick_lbls, fontsize = 24)
# ax1.yaxis.get_offset_text().set_size(25)
# plt.yticks(ytick_locs, fontsize = 24)
# plt.legend(('India OMI' ,'India UKCA' ),fontsize = 24)
# # plt.xlabel("Days 2005-2015", fontsize = 24)
# plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 24)


#plot
fig = plt.figure(facecolor='White',figsize=[41,25]);pad= 1; 
#plt.suptitle(title_list, fontsize = 25, y=0.95)

ax1 = fig.add_subplot(321)
ax1.plot(India_mean_sat,c='b',linewidth=2);
ax1.plot(India_mean_mod,c='m',linestyle='-.',linewidth=2);
# tick_locs =np.linspace(0, 4380, num=12, endpoint=False)
xtick_locs = [0,  366,  731, 1096, 1462, 1827, 2192, 2557, 2923, 3283, 3653, 4015]
xtick_lbls = ['Jan-05', 'Jan-06', 'Jan-07','Jan-08', 'Jan-09','Jan-10', 'Jan-11','Jan-12','Jan-13', 'Jan-14', 'Jan-15','Dec-15']
# ytick_locs = np.linspace(1.25e15, 1.5e16, num=8)
plt.ylim(1e14, 1.6e16)
plt.xticks(xtick_locs, xtick_lbls, fontsize = 24)
ax1.yaxis.get_offset_text().set_size(25)
plt.yticks( fontsize = 24)
plt.legend(('India OMI' ,'India UKCA' ),fontsize = 24,loc='upper right')
# plt.xlabel("Days 2005-2015", fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)

ax2 = fig.add_subplot(322)
ax2.plot(China_mean_sat,c='b',linewidth=2);
ax2.plot(China_mean_mod,c='m',linestyle='-.',linewidth=2);
plt.xticks(xtick_locs, xtick_lbls, fontsize = 24)
plt.legend(('China OMI' ,'China UKCA' ),fontsize = 24,loc='upper right')
plt.ylim(1e14, 1.6e16)
ax2.yaxis.get_offset_text().set_size(25)
plt.yticks(fontsize = 24)
# plt.xlabel("Days 2005-2015", fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)



# ytick_locs = np.linspace(1.25e15, 6e16, num=8)
ax1 = fig.add_subplot(323)
ax1.plot(North_India_mean_sat,c='b',linewidth=2);
ax1.plot(North_India_mean_mod,c='m',linestyle='-.',linewidth=2);
plt.xticks(xtick_locs, xtick_lbls, fontsize = 24)
plt.legend(('North_India OMI' ,'North_India UKCA' ),fontsize = 24,loc='upper right')
plt.ylim(1e14, 6e16)
ax1.yaxis.get_offset_text().set_size(25)
plt.yticks(fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)

ax1 = fig.add_subplot(324)
ax1.plot(East_China_mean_sat,c='b',linewidth=2);
ax1.plot(East_China_mean_mod,c='m',linestyle='-.',linewidth=2);
plt.xticks(xtick_locs, xtick_lbls, fontsize = 24)
plt.legend(('East_China OMI' ,'East_China UKCA' ),fontsize = 24,loc='upper right')
plt.ylim(1e14, 6e16)
ax1.yaxis.get_offset_text().set_size(25)
plt.yticks(fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)

# ytick_locs = np.linspace(1.25e15, 6.0e15, num=8)
ax2 = fig.add_subplot(325)
ax2.plot(South_India_mean_sat,c='b',linewidth=2);
ax2.plot(South_India_mean_mod,c='m',linestyle='-.',linewidth=2);
plt.xticks(xtick_locs, xtick_lbls, fontsize = 24)
plt.legend(('South_India OMI' ,'South_India UKCA' ),fontsize = 24,loc='upper right')
plt.ylim(1e14, 6e15)
ax2.yaxis.get_offset_text().set_size(25)
plt.yticks(fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)

ax2 = fig.add_subplot(326)
ax2.plot(West_China_mean_sat,c='b',linewidth=2);
ax2.plot(West_China_mean_mod,c='m',linestyle='-.',linewidth=2);
plt.xticks(xtick_locs, xtick_lbls, fontsize = 24)
plt.legend(('West_China OMI' ,'West_China UKCA' ),fontsize = 24,loc='upper right')
plt.ylim(1e14, 6e15)
ax2.yaxis.get_offset_text().set_size(25)
plt.yticks(fontsize = 24)
plt.ylabel("NO$_2$ (molecules/cm$^2$)", fontsize = 30)




plt.savefig('/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/ESM2019/'+Today_date+'TS_UKCAak_OMI_daily.png', dpi=300)   ##########	
plt.show()


