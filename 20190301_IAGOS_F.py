import netCDF4 as nc4
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from netCDF4 import MFDataset
import datetime as dt
import numpy as np
import matplotlib.cm as cm
from scipy import stats
import pandas as pd
import math
import glob
import os
import matplotlib.pyplot as plt
import time,datetime,calendar
from time import clock

def reading_filelist(file_path):
	nc_files = []
	for root, dirs, files in os.walk(file_path):
		for f in files:
			if f.endswith('.nc'):
				nc_files.append(file_path+f)
	# print nc_files
	return nc_files

def read_file(file_name):
	print file_name
	nc_fid = nc4.Dataset(file_name,mode='r')
	UTCtime =nc_fid.variables['UTC_time'][:]
	# time_ref='2000-01-01 00:00:00'
	time_ref=nc_fid.variables['UTC_time'].units
	# print time_ref
	time_ref_a= time_ref[14:]	
	time_rel2000 = UTCtime+calendar.timegm(time.strptime(time_ref_a,"%Y-%m-%d %H:%M:%S"))-calendar.timegm(time.strptime('2000-01-01 00:00:00',"%Y-%m-%d %H:%M:%S"))  ## Here I am forcing all your time to be seconds since 2000-01-01 00:00:00
	lat= nc_fid.variables['lat'][:]
	lon= nc_fid.variables['lon'][:]
	height=nc_fid.variables['baro_alt_AC'][:]
	ozone= nc_fid.variables['O3_PM'][:]
	missing_value = nc_fid.variables['O3_PM'].missing_value
	fill_value = nc_fid.variables['O3_PM'].missing_value
	ozone[ozone==missing_value]=np.nan
	ozone[ozone==fill_value]=np.nan
	nc_fid.close()
	return time_rel2000,lon,lat,height,ozone
	
def regrid_data(Year,month):
	lon_grid =np.arange(0.,361.,1.875)
	lat_grid =np.arange(-90.,91.,1.25)
	height_grid = np.arange(0,12500.1,500.)
	if month in (1, 3, 5, 7, 8, 10, 12):days=31 # Add +1 day because last date flight is flying to the next month 1st day as well or delete the last date from the working folder
	elif month ==2:days=28
	else:days=30
	time_grid = np.arange(0,days*24*60*60,60*60)+calendar.timegm(time.strptime(str(Year)+'-'+str(month)+'-01 00:00:00',"%Y-%m-%d %H:%M:%S"))-calendar.timegm(time.strptime('2000-01-01 00:00:00',"%Y-%m-%d %H:%M:%S"))    ## seconds since 2000 to be output  and saved in NetCDF with an increasement of an hour (60*60 sec), I am here only initializing a monthly array, as I got memory issues with an annual array which is too large
	# time_grid = np.arange(0,720*60*60,60*60)+calendar.timegm(time.strptime(str(Year)+'-'+str(month)+'-01 00:00:00',"%Y-%m-%d %H:%M:%S"))-calendar.timegm(time.strptime('2000-01-01 00:00:00',"%Y-%m-%d %H:%M:%S","%Y-%m-%d %H:%M:%S"))    ## seconds since 2000 to be output  and saved in NetCDF with an increasement of an hour (60*60 sec), I am here only initializing a monthly array, as I got memory issues with an annual array which is too large
	# time_grid = np.arange(0,720*60*60,60*60)+calendar.timegm(time.strptime(str(Year)+'-01-01 00:00:00',"%Y-%m-%d %H:%M:%S"))-calendar.timegm(time.strptime('2000-01-01 00:00:00',"%Y-%m-%d %H:%M:%S"))    ## seconds since 2000 to be output  and saved in NetCDF with an increasement of an hour (60*60 sec), I am here only initializing a monthly array, as I got memory issues with an annual array which is too large
	time_units ='Seconds since 2000-01-01 00:00:00'
	ozone_grid=np.empty((len(time_grid),len(height_grid),len(lat_grid),len(lon_grid)));ozone_grid[:]=np.nan
	MeaTag=np.empty((len(time_grid),len(height_grid),len(lat_grid),len(lon_grid)));MeaTag[:]=0  ## Tag to count the no of measurements
	ConfiTag=np.empty((len(time_grid),len(height_grid),len(lat_grid),len(lon_grid)));ConfiTag[:]=0  ## Tag indicating if the average is confident 1-Y, 0-N
	def produce_file_name(Year,month):
		if month<=9:
			month='0'+str(month);
		else:
			month =str(month);
		file_path =  r'/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/observations/IAGOS/'+str(Year)+'/'+str(month)+'/'
		# print file_path
		return file_path
	file_path=produce_file_name(Year,month)
	# print 'file_path', file_path	
	files = reading_filelist(file_path)
	#files=nc_files ##  Please change this into a list containing all files of the specific Year, I believe you should know how to do this
	for ifile in files:
		time_rel2000,lon,lat,height,ozone = read_file(ifile);
		second_offset = calendar.timegm(time.strptime('2000-01-01 00:00:00',"%Y-%m-%d %H:%M:%S"))-calendar.timegm(time.strptime(str(Year)+'-'+str(month)+'-01 00:00:00',"%Y-%m-%d %H:%M:%S"))
		for item in range(len(time_rel2000)):
			time_index = int(np.ceil((time_rel2000[item]+second_offset)/3600)); #print time_index
			height_index = int(np.ceil(height[item]/500)); #print height_index
			lon_index = int(np.ceil(lon[item]/1.875));#print lon_index
			lat_index = int(np.ceil((lat[item]+90)/1.25));#print lat_index
			ExistingMeasueNo = MeaTag[time_index,height_index,lat_index,lon_index]	
			if ExistingMeasueNo == 0:
				# print '1d',ozone[item]
				ozone_grid[time_index,height_index,lat_index,lon_index]=ozone[item]
				# print  '4d',ozone_grid[time_index,height_index,lon_index,lat_index]
				MeaTag[time_index,height_index,lat_index,lon_index]=MeaTag[time_index,height_index,lon_index,lat_index]+1	
			else:
				ReAvg = (ExistingMeasueNo*ozone_grid[time_index,height_index,lat_index,lon_index]+ozone[item])/(ExistingMeasueNo+1)	
				ozone_grid[time_index,height_index,lat_index,lon_index]=ReAvg
				MeaTag[time_index,height_index,lat_index,lon_index]=MeaTag[time_index,height_index,lat_index,lon_index]+1
			#print ozone_grid[time_index,height_index,lon_index,lat_index]
	ConfiTag[MeaTag>=1] = 1; ## Change this to 5 for your case, I am using 1 as there are only two fake flights in the demo. However, I am not sure if judging if the estimate is rliable based on the number of measurements being avergaed is sensible, I will talk about this with you. 
	Tag = ConfiTag+0.1*MeaTag
	#ozone_grid[ozone_grid==0]=np.nan
	return time_grid,time_units,height_grid,lat_grid,lon_grid,ozone_grid,Tag

def netcdf4_write(file_name,time_series,lat,lon,height,ozone):
	
	"""
	Here is how you would normally create and store data in a netCDF file:
    1) Open/create a netCDF dataset.
    2) Define the dimensions of the data.
    3) Construct netCDF variables using the defined dimensions.
    4) Pass data into the netCDF variables.
    5) Add attributes to the variables and dataset (optional but recommended).
    6) Close the netCDF dataset.
	"""
	# 1) Open/create a netCDF dataset.
	f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write	
	"""
	The above line creates a netCDF file called "file_name" in the filepath folder.
	f is a netCDF Dataset object that provides methods for storing data to the file. 
	f also doubles as the root group. A netCDF group is basically a directory or 
	folder within the netCDF dataset. This allows you to organize data as you would 
	in a unix file system. Let's create a group for the heck of it
	"""
	
	# 2) Define the dimensions of the data.
	"""
	netCDF defines the sizes of all variables in terms of dimensions, so before any 
	variables can be created the dimen sions they use must be created first. A special 
	case, not often used in practice, is that of a scalar variable, which has no dimensions.
	A dimension is created using the Dataset.createDimension method of a Dataset. A Python 
	string is used to set the name of the dimension, and an integer value is used to set 
	the size. To create an unlimited dimension (a dimension that can be appended to), t
	he size value is set to None or 0. In this example, the time dimension is unlimited. 
	In netCDF4 files you can have more than one unlimited dimension, in netCDF 3 files 
	there may be only one, and it must be the first (leftmost) dimension of the variable
	"""
	f.createDimension('time', len(time_series))
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	f.createDimension('height', len(height))
	
	#3) Construct netCDF variables using the defined dimensions.
	
	times = f.createVariable('time',np.float64, ('time'))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	heights = f.createVariable('height',np.float32, ('height'))
	ozones = f.createVariable('ozone',np.float32,('time','height','lat','lon'))
	'''
	rx5days = f.createVariable('rx5day',np.float32,('time','lat','lon'))
	sdiis = f.createVariable('sdii',np.float32,('time','lat','lon'))
	r10s = f.createVariable('r10',np.float32,('time','lat','lon'))
	r20s = f.createVariable('r20',np.float32,('time','lat','lon'))
	rnms = f.createVariable('rnm',np.float32,('time','lat','lon'))
	cdds = f.createVariable('cdd',np.float32,('time','lat','lon'))
	cwds = f.createVariable('cwd',np.float32,('time','lat','lon'))
	r95ps = f.createVariable('r95p',np.float32,('time','lat','lon'))
	r99ps = f.createVariable('r99p',np.float32,('time','lat','lon'))
	prcptots = f.createVariable('precptot',np.float32,('time','lat','lon'))
	total_precips = f.createVariable('total_precip',np.float32,('time','lat','lon'))
	mean_preps = f.createVariable('mean_precip',np.float32,('time','lat','lon'))
	std_preps = f.createVariable('std_precip',np.float32,('time','lat','lon'))
	MPIs = f.createVariable('MPI',np.float32,('time','lat','lon'))'''
	
	#4) Passing data into variables
	times[:] = time_series
	latitudes[:] = lat
	longitudes[:] = lon
	heights[:] = height
	ozones[:] = ozone
	'''
	rx5days[:] = rx5day
	sdiis[:] = sdii
	r10s[:] = r10
	r20s[:] = r20
	rnms[:] = rnm
	cdds[:] = cdd
	cwds[:] = cwd
	r95ps[:] = r95p
	r99ps[:] =  r99p
	prcptots[:] =  precptot
	total_precips[:] = total_precip 
	mean_preps[:] = mean_prep
	std_preps[:] = std_prep
	MPIs[:] = MPI'''
	
	# 5) Add attributes to the variables and dataset (optional but recommended).
	"""
	There are two types of attributes in a netCDF file, global and variable. Global attributes provide information
	about the entire dataset, as a whole. Variable attributes provide information about one of the variables. Global
	attributes are set by assigning values to Dataset instance variables. Variable attributes are set by assigning
	values to Variable instances variables. Attributes can be strings, numbers or sequences. 
	"""
	
	# Global Attributes
	"""
	title       : Whats in the file
	institution : Where it was produced
	source      : How it was produced e.g. model version, instrument type
	history     : Audit trail of processing operations
	references  : Pointers to publications or web documentation
	comment     : Miscellaneous
	"""
	f.description = 'IAGOS Dataset'
	# f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alok Pandey at the University of Edinburgh'
	
	# Variable Attributes
	"""
	units               : mandatory for all variables containing data other than dimensionless numbers. 
	standard_name       : identifies the quantity. Units must be consistent with standard name and any 
						  statistical processing.Standard name does not include coordinate or processing information.
	long_name           : not standardised.
	ancillary variables : a pointer to variables providing metadata about the individual data values e.g. standard error 
						  or data quality information. Numeric data variables may have 
						  ##############################################################################
						  # FillValue, missing_value, valid_max, valid_min, valid_range. missing_value #
						  ##############################################################################
	"""
	times.long_name = 'Year'
	times.units = '1'
	
	latitudes.long_name = 'latitude'
	latitudes.units = 'degree_north'
	longitudes.long_name = 'longitude'
	longitudes.units = 'degree_east'
	heights.units = 'meters'
	
	ozones.standard_name = 'Ozone concentration'
	ozones.long_name = 'Ozone hourly concentration'
	ozones.units = 'ppb'
	
	'''
	rx5days.standard_name = 'maximum 5-day precip amount'
	rx5days.long_name = 'maximum 5-day precipitation amount'
	rx5days.units = 'mm'
	
	sdiis.standard_name = 'simple daily intensity index'
	sdiis.long_name = 'annual precip amount devided by the number of wet days (> 1mm)'
	sdiis.units = 'm/s/day'
	
	r10s.standard_name = 'number of heavy precip days'
	r10s.long_name = 'nnual count of days woth precipitation >= 10mm'
	r10s.units = 'day'
	
	r20s.standard_name = 'number of very heavy precip days'
	r20s.long_name = 'nnual count of days woth precipitation >= 20mm'
	r20s.units = 'day'	

	rnms.standard_name = 'number of days above nn mm here defines as 30 mm'
	rnms.long_name = 'annual count of days woth precipitation >= 30mm'
	rnms.units = 'day'	
	
	cdds.standard_name = 'consective dry days'
	cdds.long_name = 'maximum consecive days with precipitation<1mm'
	cdds.units = 'day'	
	
	cwds.standard_name = 'consective wet days'
	cwds.long_name = 'maximum consecive days with precipitation>=1mm'
	cwds.units = 'day'	
	
	r95ps.standard_name = 'very wet days'
	r95ps.long_name = 'annual total precipitation >95th percetile'
	r95ps.units = 'day'	
	
	r99ps.standard_name = 'extremely wet days'
	r99ps.long_name = 'annual total precipitation >99th percetile'
	r99ps.units = 'day'	

	prcptots.standard_name = 'annual total wet day prrecip'
	prcptots.long_name = 'annual total precipitation in wet days (rr>1mm)'
	prcptots.units = 'mm'		
	
	total_precips.standard_name = 'annual total prrecip'
	total_precips.long_name = 'annual total precipitation'
	total_precips.units = 'mm
	
	MPIs.standard_name = 'Monsoon Precipitation Index'
	MPIs.Long_name = 'Ratio of local summer(MJJAS)-winter(NDJFM) to the annual total for the North Hemisphere'''
	
	f.close()
	print 'Done!!!'
	return 1
	
def main(Years, imonth):
	for Year in Years:
		for imonth in range(1,13):
			time_series,time_units,height,lat,lon,ozone,Tag = regrid_data(Year,imonth)
			print np.nanmax(ozone),np.nanmin(ozone)
			# return time,time_units,height,lon,lat,ozone,Tag
			file_name='IAGOS_AKP.'+str(Year)+'.'+str(imonth)+'.nc'
			netcdf4_write(file_name,time_series,lat,lon,height,ozone)
		### Now save these variables in NetCDFs.....
		### ....
	return 1
	
Years = np.arange(2000,2001);
imonth =np.arange(1,13);
main(Years, imonth)