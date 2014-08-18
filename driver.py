from SSE import *
import datetime
import cPickle as pickle
import sys
sys.path.append("/home/ice/nchaney/UTILS/PYTHON")
import IO
import gdal

#Script to drive the comparison of methods to merge stations and gridded data
station_coords = '/home/raid20/nchaney/Data/Met_Data/Mesonet_Test/Mesonet_Stations.txt'
station_data_file = '/home/ice/nchaney/PROJECTS/SM_validation/SENSITIVITY/SM_LR_2001-2009_1hr.pck'#'/home/ice/nchaney/Dropbox/33883094.csv'
file_in = 'original.tif'
file_out = 'corrected.tif'
file_station = 'station_data.txt'
minlat = 31.45
minlon = -83.78
nlat = 40
nlon = 35
res = 2*0.00416666666666570
maxlat = minlat + res*(nlat-1) - res/2
maxlon = minlon + res*(nlon-1) - res/2
undef = -9.99e+08
dims = {}
dims['minlat'] = minlat #-89.8750
dims['minlon'] = minlon #0.1250
dims['nlat'] = nlat #720
dims['nlon'] = nlon #1440
dims['res'] = res
dims['maxlat'] = dims['minlat'] + dims['res']*(dims['nlat']-1)
dims['maxlon'] = dims['minlon'] + dims['res']*(dims['nlon']-1)
dims['undef'] = -9.99e+08
time = datetime.datetime(2004,1,5,18)
#ftime = datetime.datetime(2004,1,31,23)
dt = datetime.timedelta(hours=1)

#Read in all the station data
station_data = pickle.load(open(station_data_file))

#Read in the metadata
metadata = retrieve_metadata(file_in)

#Read in the gridded data
original = np.flipud(read_raster(file_in))

#Extract station data
time_station = time - datetime.timedelta(hours = 5)
var = 'apcpsfc'
data = {'id':[],var:[]}
st_coords = {'id':[],'lon':[],'lat':[]}
for station in station_data:
 if station == 'ID': continue
 if station_data[station]['LON'] >= maxlon or station_data[station]['LON'] <= minlon: continue
 if station_data[station]['LAT'] >= maxlat or station_data[station]['LAT'] <= minlat: continue
 idx = list(station_data[station]['DATE']).index(time_station)
 if time_station != datetime.datetime(time.year,1,1,0):
  value = station_data[station]['PREC'][idx] - station_data[station]['PREC'][idx-1]
 else: 
  value = station_data[station]['PREC'][idx]
 if np.isnan(value) == 1: continue
 data['id'].append(station)
 data[var].append(value)
 st_coords['id'].append(station)
 st_coords['lat'].append(station_data[station]['LAT'])
 st_coords['lon'].append(station_data[station]['LON'])
data['id'] = np.array(data['id'])
data[var] = np.array(data[var])
st_coords['id'] = np.array(st_coords['id'])
st_coords['lat'] = np.array(st_coords['lat'])
st_coords['lon'] = np.array(st_coords['lon'])
print st_coords
print data

#Write to file
fp = open('station_data.txt','wb')
fp.write('id lat lon data\n')
for id in st_coords['id']:
 mask = st_coords['id'] == id
 lat = st_coords['lat'][mask]
 lon = st_coords['lon'][mask]
 mask = data['id'] == id
 fp.write('%d %.5f %.5f %.5f\n' % (id,lat,lon,25.4*data['apcpsfc'][mask]))
fp.close()

#Read in the station data
station_data = np.loadtxt(file_station,skiprows=1)
st_coords = {'id':station_data[:,0],'lat':station_data[:,1],'lon':station_data[:,2]}
data = {'id':station_data[:,0],'data':station_data[:,3]}

#Extract all the info for the stations
data_measurements = Extract_Station_Data(st_coords,data,original,minlat,minlon,res,undef,'data')

#Compute the sample semivariogram
Compute_Sample_Semivariogram(data_measurements)

#Fit a model semivariogram
r("vnew1<-vnew")
vg_fit = r("v.fit <- fit.variogram(vnew1,model=vgm(model='Sph',range=1000.0),fit.method=1)")

#Assimilate the measurements
corrected = Assimilate_Measurements(data_measurements,original,vg_fit,nlat,nlon,undef,res,minlat,minlon,original)

#Plot the comparison
Plot_Comparison(minlat,minlon,maxlat,maxlon,nlat,nlon,res,original,data_measurements,corrected,undef)

#Write data
write_raster(np.flipud(corrected[0]),metadata,file_out)
 
