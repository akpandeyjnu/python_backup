import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import numpy as np
import matplotlib
from matplotlib.colors import LogNorm
# from iris.analysis.interpolate import extract_nearest_neighbour
import matplotlib.dates as mdates
import iris.analysis.cartography
import pandas
import matplotlib.colors as colors
import iris.coords
import pylab as pl
from datetime import timedelta, datetime
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
import iris.analysis.cartography as iac
cmap = matplotlib.cm.get_cmap('brewer_RdBu_11')

# Iris uses a strange coordinate system that make it difficult to extract regions, so create a function which will chnage the cooridnates of your date. Run this function on your imported data.
#---------------------------
def roll_cube(cube):
    """Convert lomgitudinal coordinate system from 0-360 to -180-180 - this will be needed for extraction of regions covering the origin"""

    lon = cube.coord('longitude')

    cube.data = np.roll(cube.data, len(lon.points) // 2)
    lon.points -= 180.
    if lon.bounds is not None:
        lon.bounds -= 180.
    return cube
#----------------------------
#titles_list = []


# Import the data, run the 'roll_cube' function, extract the first model level
filename = '/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/Model_Evaluation/ppfile/test/xmlxpa.pm*.pp'
stash_constraint_O3 = iris.AttributeConstraint(STASH='m01s34i001')    

cube = iris.load(filename,stash_constraint_O3)[0]
print cube.shape


rolled_cube = roll_cube(cube)
print rolled_cube

"""
#for psedo level1---------------------------
surface_cube1 = rolled_cube
print surface_cube1

surface_cube= surface_cube1*0.60* 1e9 


surface_cube.coord('latitude').guess_bounds()
surface_cube.coord('longitude').guess_bounds()



# Now extract the region accoridning to cooridnates you desire - here I'm extracting IGP W
lat_extract_W=surface_cube.extract(iris.Constraint(latitude=lambda lat: 27.5 <= lat < 32))
lon_extract_W=lat_extract_W.extract(iris.Constraint(longitude=lambda lon: 71.5 <= lon < 78))
# Now the area has been extracted, you need to average over the region. however, the gridd cells are NOT the same size, so you need an area weigted average.
aw_array_W = iris.analysis.cartography.area_weights(lon_extract_W)
averaged_area_W = lon_extract_W.collapsed(['latitude', 'longitude'],iris.analysis.MEAN,weights=aw_array_W)
#IGP C
lat_extract_C=surface_cube.extract(iris.Constraint(latitude=lambda lat: 24.5 <= lat < 29))
lon_extract_C=lat_extract_C.extract(iris.Constraint(longitude=lambda lon: 78 <= lon < 84.5))
# Now the area has been extracted, you need to average over the region. however, the gridd cells are NOT the same size, so you need an area weigted average.
aw_array_C = iris.analysis.cartography.area_weights(lon_extract_C)
averaged_area_C = lon_extract_C.collapsed(['latitude', 'longitude'],iris.analysis.MEAN,weights=aw_array_C)
#IGP E
lat_extract_E=surface_cube.extract(iris.Constraint(latitude=lambda lat: 22 <= lat < 26.5))
lon_extract_E=lat_extract_E.extract(iris.Constraint(longitude=lambda lon: 84.5<= lon < 91))
# Now the area has been extracted, you need to average over the region. however, the gridd cells are NOT the same size, so you need an area weigted average.
aw_array_E = iris.analysis.cartography.area_weights(lon_extract_E)
averaged_area_E = lon_extract_E.collapsed(['latitude', 'longitude'],iris.analysis.MEAN,weights=aw_array_E)


A = averaged_area_W.data
print A
len(A)
A.shape

A1 = A[:,:85]
A1.shape

AA =A1.transpose()
print AA 

B = averaged_area_C.data
print B

B1 = B[:,:85]
B1.shape

BB =B1.transpose()
print BB 

C = averaged_area_E.data
print C

C1 = C[:,:85]
C1.shape

CC =C1.transpose()
print CC 

'''#export as csv not working export individualy if needed
W=AA,BB,CC
print W

header = 'Lahore, Amritsar, Chandigah, NewDelhi, Agra, Pantnagar, Kanpur, Lucknow, Allahabad, Gandhi_College_Ballia, Varanashi, Patna, Durgapur, Kolkata,  Dhaka_University'
np.savetxt("/exports/csce/datastore/geos/users/s1583507/Python/CSV_Files_F/20160604_vertical_timeseries_O3.csv", W,  delimiter=",", fmt='%s', header=header, footer='IGPW, IGPC, IGPE, IGP')
np.savetxt("/exports/csce/datastore/geos/users/s1583507/Python/CSV_Files_F/20160604_vertical_timeseries_O3_1.csv",np.transpose([AA,BB,CC]), delimiter=",", fmt='%s', header=header)
'''

# Now plot the data
Y = [20, 53, 100, 160, 233, 320, 420, 533, 660, 800, 
     953, 1120, 1300, 1493, 1700, 1920, 2153, 2400, 2660, 2933,
     3220, 3520, 3833, 4160, 4500, 4853, 5220, 5600, 5993, 6400,
     6820, 7253, 7700, 8160, 8633, 9120, 9620, 10133, 10660, 11200,
     11754, 12321, 12901, 13495, 14102, 14724, 15359, 16009, 16673, 17352,
     18046, 18757, 19484, 20229, 20993, 21777, 22582, 23412, 24268, 25153,
     26071, 27024, 28018, 29058, 30150, 31301, 32518, 33811, 35190,	36666,
     38254	,39968, 41825, 43844, 46046, 48456, 51099, 54006, 57210, 60747, 64657, 68986, 73782, 79100,	85000]
    
print Y

X= np.linspace(1,24,24)
print X

levels = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000]

plt.figure(1)
plt.figure(figsize=(10,15))
plt.suptitle ('Monthly Ozone variation in IGP', fontsize = 20, y=0.95)

ax = plt.subplot(311)
#plt.contourf(AA)
plt.contourf(AA, levels, vmax=11000, vmin=0)

plt.title('O3 IGP W',fontsize = 18, y=1.01)
plt.ylabel('Height (meters)', fontsize = 16, y =0.5)
#plt.xlabel('Months')
cb1 = plt.colorbar(aspect = 25, ticks=levels, orientation='vertical',fraction=0.1,pad=0.01)
#cb1 = plt.colorbar(aspect = 25,  orientation='vertical',fraction=0.1,pad=0.01)
cb1.set_label('O3 (ppbv)', fontsize = 16)

tick_locs = [0,5,11,17,23]
tick_lbls = ['Jan-00','Jun-00', 'Dec-00','Jun-01', 'Dec-01']
plt.xticks(tick_locs, tick_lbls, fontsize = 12)

tick_locs1 = [1,10,20,30,40,50,60,70,80,84]
tick_lbls1 = ['20','800', '2933','6400', '11200','17352','25153','36666','85000']
plt.yticks(tick_locs1, tick_lbls1, fontsize = 12, y=1.01)

#-------------------------------------------------------------------------------------------------------
    
bx = plt.subplot(312)
plt.contourf(BB, levels, vmax=11000, vmin=0)

plt.title('O3 IGP C',fontsize = 18, y=1.01)
plt.ylabel('Height (meters)',fontsize = 16, y=0.5)
#plt.xlabel('Months',fontsize = 11, y=0.02)
#cb1 = plt.colorbar(aspect = 25, orientation='vertical',fraction=0.1,pad=0.01)
cb1 = plt.colorbar(aspect = 25, ticks=levels, orientation='vertical',fraction=0.1,pad=0.01)
cb1.set_label('O3 (ppbv)', fontsize = 16)

tick_locs = [0,5,11,17,23]
tick_lbls = ['Jan-00','Jun-00', 'Dec-00','Jun-01', 'Dec-01']
plt.xticks(tick_locs, tick_lbls, fontsize = 12)

tick_locs1 = [1,10,20,30,40,50,60,70,80,84]
tick_lbls1 = ['20','800', '2933','6400', '11200','17352','25153','36666','85000']
plt.yticks(tick_locs1, tick_lbls1, fontsize = 12, y=1.01)


cx = plt.subplot(313)
plt.contourf(CC, levels, vmax=11000, vmin=0)

plt.title('O3 IGP E',fontsize = 18, y=1.01)
plt.ylabel('Height (meters)',fontsize = 16, y=0.5)
plt.xlabel('Months', fontsize = 12, y=0.01)
#cb1 = plt.colorbar(aspect = 25, orientation='vertical',fraction=0.1,pad=0.01)
cb1 = plt.colorbar(aspect = 25, ticks=levels, orientation='vertical',fraction=0.1,pad=0.01)
cb1.set_label('O3 (ppbv)', fontsize = 16)

tick_locs = [0,5,11,17,23]
tick_lbls = ['Jan-00','Jun-00', 'Dec-00','Jun-01', 'Dec-01']
plt.xticks(tick_locs, tick_lbls, fontsize = 12)

tick_locs1 = [1,10,20,30,40,50,60,70,80,84]
tick_lbls1 = ['20','800', '2933','6400', '11200','17352','25153','36666','85000']
plt.yticks(tick_locs1, tick_lbls1, fontsize = 12, y=1.01)

plt.show()

# plt.savefig('/exports/csce/datastore/geos/users/s1583507/Python/Images_F/20160614_vertical_timeseries_O3.png', dpi=300)"""

