import netCDF4 as nc4
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import cartopy.feature as cfeature
from mpl_toolkits.basemap import Basemap, shiftgrid
import cartopy.crs as ccrs
import matplotlib.cm as cm
cmap = cm.jet
import datetime
# from datetime import datetime
# import cf, cfplot as cfp

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

def spatial_figure(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True): #c_bad,c_under,c_over,c_number=20,
    """
	input : all parameters and data rel;ated to the figure you want to plot_title
		lons and lats are 1-d array while data is 2-D array
		colorbar_min,colorbar_max specifies the minimum and maximum value you want to show with the hottest and coldest color respectively
		tb_lef and tb_bot specifies if you want to have axis labels, True fro yes
	output : a spatial map of the data
	"""
    lons[lons>180]-=360; 
    #lon_b = np.min(lons); lon_e = np.max(lons)
    lon_b = -180; lon_e = 180
    #lat_b = np.min(lats); lat_e = np.max(lats)	
	
    lat_b = -90; lat_e = 90
    lon_bin = 60; lat_bin = 30
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
	

#Getting data from netcdf file
file_path = '/home/v1apande/group_datastore/Research_Work/Python_notepad/IAGOS_AKP.2000.2.nc'
nc_file = nc4.Dataset(file_path,mode='r')
print nc_file.file_format
print nc_file

print len(nc_file.variables)
#print nc_file.Conventions
for i in nc_file.variables:
    #print [i]
    print [i, nc_file.variables[i].shape]
for attr in nc_file.ncattrs():
    print attr, '=', getattr(nc_file, attr)

print nc_file.dimensions.keys()
print nc_file.dimensions['lat']

print nc_file.variables.keys()
# print nc_file.variables
# print nc_file.variables['ozone']

lats = nc_file.variables['lat'][:]
lons = nc_file.variables['lon'][:]
time = nc_file.variables['time'][:]
height =nc_file.variables['height'][:]

O3 = nc_file.variables['ozone'][:]
O3units=nc_file.variables['ozone'].units
print np.nanmax(O3), np.nanmin(O3)
print O3.shape
nc_file.close()

'''
O3_A = O3[:,1,:,:]

print np.nanmax(O3_A), np.nanmin(O3_A)

print O3_A.shape
# print O3_A


fig = plt.figure(facecolor='White',figsize=[10,6]);pad= 2; 

colormap='RdBu_r'; colorbar_min=0;colorbar_max=80 ## change this accordingly
lons = lons ## 1-D
lats= lats##1-D
data=O3_A ##2-D  You do not need to shift the data as my funciton will do that for you in line 32

ax = plt.subplot(1,1,1);
colormesh_1 = spatial_figure(ax,data,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# plt.suptitle ('UKCA Ground Level Ozone', fontsize = 20, y=1)
# plt.title (''+title_list+' 2005', fontsize = 16, y=1)
plt.show()
# O3 = O3a *0.60 *1e9
'''

# time1 = time.tolist()#convet float to list
index = [ix for ix ,val in enumerate(time) if val != np.nonzero]
# print index

# index1 = index[51:65]
for idx in range(525,555): #modify as per requirement
	# print idx
	O3_1 = O3[idx,1,:,:]
	print O3_1.shape
	print np.nanmax(O3_1),np.nanmin(O3_1)
	title_list = str(idx+1)
	print title_list

	### use the funciton

	fig = plt.figure(facecolor='White',figsize=[10,6]);pad= 2; 
	colormap='RdBu_r'; colorbar_min=0;colorbar_max=300 ## change this accordingly
	
	lons = lons ## 1-D
	lats= lats##1-D
	data=O3_1 ##2-D  You do not need to shift the data as my funciton will do that for you in line 32

	ax = plt.subplot(1,1,1);
	colormesh_1 = spatial_figure(ax,data,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
	# plt.suptitle ('UKCA Ground Level Ozone', fontsize = 20, y=1)
	# plt.title (''+title_list+' 2005', fontsize = 16, y=1)
		
## Add coloabar
	cbar_ax = fig.add_axes([0.11, 0.035, 0.80, 0.03])  # position and size of the colorbar
	char = fig.colorbar(colormesh_1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
	cbar_ax.annotate(r'Ozone (ppb)',xy=(0.5,0.9), xytext=(0, pad),
			xycoords='axes fraction', textcoords='offset points',
			ha='center', va='bottom',rotation='horizontal',fontsize=15)
	plt.show()
    #plt.savefig('/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/ncfile_test/'+Today_date+ 'UKCA_nc as_'+title_list+'.png', dpi=300)"""
