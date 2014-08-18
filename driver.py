from SSE import *

#Change these files to your data. The formats need to match.
file_in = 'example/original.tif'
file_out = 'example/corrected.tif'
file_station = 'example/station_data.txt'
undef = -9.99e+08

#Read in the metadata
metadata = retrieve_metadata(file_in)
res = metadata['resx']
minlat = metadata['miny'] + res/2
minlon = metadata['minx'] + res/2
maxlat = metadata['miny'] - res/2
maxlon = metadata['maxy'] - res/2
nlat = metadata['ny']
nlon = metadata['nx']

#Read in the gridded data
original = np.flipud(read_raster(file_in))

#Read in the station data
station_data = np.loadtxt(file_station,skiprows=1)
st_coords = {'id':station_data[:,0],'lat':station_data[:,1],'lon':station_data[:,2]}
data = {'id':station_data[:,0],'data':station_data[:,3]}

#Extract all the info for the stations
data_measurements = Extract_Station_Data(st_coords,data,original,minlat,minlon,res,undef,'data')

#Compute the sample semivariogram
Compute_Sample_Semivariogram(data_measurements)

#Fit a model semivariogram
vg_fit = r("v.fit <- fit.variogram(vnew,model=vgm(model='Sph',range=1000.0),fit.method=1)")

#Assimilate the measurements
corrected = Assimilate_Measurements(data_measurements,original,vg_fit,nlat,nlon,undef,res,minlat,minlon,original)

#Plot the comparison
Plot_Comparison(minlat,minlon,maxlat,maxlon,nlat,nlon,res,original,data_measurements,corrected,undef)

#Write data
write_raster(np.flipud(corrected[0]),metadata,file_out)
 
