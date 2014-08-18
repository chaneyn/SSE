from SSE import *
import datetime
import cPickle as pickle
import sys
sys.path.append("/home/ice/nchaney/UTILS/PYTHON")
import IO

#Script to drive the comparison of methods to merge stations and gridded data
stageiv_ctl = '/home/raid19/nchaney/StageIV_NLDAS/stageiv.ctl'
station_coords = '/home/raid20/nchaney/Data/Met_Data/Mesonet_Test/Mesonet_Stations.txt'
station_data_file = '/home/ice/nchaney/PROJECTS/SM_validation/SENSITIVITY/SM_LR_2001-2009_1hr.pck'#'/home/ice/nchaney/Dropbox/33883094.csv'
file_out = 'LR_2004-2009.nc'
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
itime = datetime.datetime(2004,1,5,18)
ftime = datetime.datetime(2004,1,5,18)
#ftime = datetime.datetime(2004,1,31,23)
dt = datetime.timedelta(hours=1)

#Define the variable to extract
var = 'apcpsfc'

#Create the output variable
data_output = Create_Data_Output_Variable()

#Initialize the output dataset
vars = [var,]
vars_info = [var,]
nt = ((ftime - itime).days + 1)*24
tstep = 'hours'
fp_out = IO.Create_NETCDF_File(dims,file_out,vars,vars_info,itime,tstep,nt)

#Open gridded file
ga.open(stageiv_ctl)
ga("set lat %f %f" % (minlat,maxlat))
ga("set lon %f %f" % (minlon,maxlon))

#Read in all the station data
station_data = pickle.load(open(station_data_file))

time = itime
t = 0
while time <= ftime:
 
 print time

 #Read in gridded data
 ga("set time %s" % time.strftime('%HZ%d%b%Y'))
 ga("data=re(%s,%d,linear,%f,%f,%d,linear,%f,%f,bl)" % (var,nlon,minlon,res,nlat,minlat,res))
 ga("set geotiff original.tif")
 ga("set gxout geotiff")
 ga("d data")
 exit()
 bdata = np.ma.getdata(ga.exp("data"))
 bdata[bdata<0] = 0.0

 #Extract station data
 time_station = time - datetime.timedelta(hours = 5)
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

 #If the mean of both is zero then no need to merge anything...
 if np.mean(bdata) == 0 and np.nanmean(data[var]) == 0: 
  #Add the data to the file
  fp_out.variables[var][t] = bdata
  #Update the time step
  t = t + 1
  time = time + dt
  continue

 #Extract all the info for the stations
 data_measurements = Extract_Station_Data(st_coords,data,bdata,minlat,minlon,res,undef,var)

 #Compute the sample semivariogram
 Compute_Sample_Semivariogram(data_measurements)

 #Fit a model semivariogram
 r("vnew1<-vnew")
 vg_fit = r("v.fit <- fit.variogram(vnew1,model=vgm(model='Sph',range=1000.0),fit.method=1)")

 #Assimilate the measurements
 data_corrected = Assimilate_Measurements(data_measurements,bdata,vg_fit,nlat,nlon,undef,res,minlat,minlon,bdata)

 #Add the data to the file
 output = data_corrected[0]
 output[output<0] = 0.0
 fp_out.variables[var][t] = data_corrected[0]
 
 #Plot the comparison
 Plot_Comparison(minlat,minlon,maxlat,maxlon,nlat,nlon,res,bdata,data_measurements,data_corrected,undef)
 
 #Update the time step
 time = time + dt
 t = t + 1

#Close the output netcdf file
fp_out.close()
