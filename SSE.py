import grads
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap,pyproj
import math
from scipy import interpolate
import numpy as np
from rpy2.robjects import r
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
import numpy.linalg as la
import datetime
import time
import pickle
import random
import gdal
import osgeo
from osgeo import osr

#Load R libraries
r('library("gstat")')
r('library("sp")')
#Start GrADS
ga = grads.GrADS(Bin='grads',Window=False,Echo=False)

def Compute_nRMSD(data_corrected,data_baseline,undef):
 
 #Compute the normalized root mean squared deviation
 idx = np.array(data_corrected > undef)
 n = idx.size
 RMSD = (np.sum((data_corrected[idx] - data_baseline[idx])**2)/n)**0.5
 if np.max(data_baseline[idx]) - np.min(data_baseline[idx]) != 0:
  nRMSD = 100*RMSD/(np.max(data_baseline[idx]) - np.min(data_baseline[idx]))
 else:
  nRMSD = 0.0

 return nRMSD

def Fit_Semivariogram_Climatology(vg_clim,data_output):

 #Fit a model semivariogram to each month
 for i in xrange(0,12):
  vg = r.assign("vnew",vg_clim[i+1]['svg'])
  vg_fit = r("v.fit <- fit.variogram(vnew,model=vgm(model='Exp',range=1000.0),fit.method=1)")
  h = np.linspace(0,np.max(np.array(vg[1])),100)
  sv = calculate_semivariogram(vg_fit,h)
  data_output["vg_clim"][i+1] = {}
  data_output["vg_clim"][i+1]["svg"] = vg_clim[i+1]['svg']
  data_output["vg_clim"][i+1]["vgfit"] = vg_fit
  plt.subplot(3,4,i+1)
  plt.plot(np.array(vg[1]),np.array(vg[2]),'bo',alpha=0.1)
  plt.plot(h,sv,'r',lw=2)
 plt.show()
 return data_output

def Plot_Comparison(minlat,minlon,maxlat,maxlon,nlat,nlon,res,bdata2,data_measurements,data_corrected,undef):
 #Plot the comparison
 llcrnrlat=minlat-res/2
 urcrnrlat=minlat+res*(nlat-1)+res/2
 llcrnrlon=minlon-res/2
 urcrnrlon=minlon+res*(nlon-1)+res/2
 fig = plt.figure(figsize=(15,4))
 ax = fig.add_subplot(132)
 vmin = np.min(bdata2)
 vmax = np.max(bdata2)
 print vmin,vmax
 ax.m = Basemap(projection='cyl',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='c')
 x,y = ax.m(data_measurements["lon"],data_measurements["lat"])
 h = ax.m.scatter(x,y,c=data_measurements["data"],s=25,vmin=vmin,vmax=vmax)
 plt.colorbar()
 data_corrected = np.array(data_corrected)
 idx = np.where(data_corrected < undef)
 data_corrected[idx] = np.nan
 ax = fig.add_subplot(131)
 ax.m = Basemap(projection='cyl',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='c')
 lons,lats = ax.m.makegrid(nlon,nlat)
 x,y = ax.m(lons,lats)
 #cs = ax.m.contourf(x,y,bdata2,vmin=vmin,vmax=vmax)
 cs = ax.m.imshow(bdata2,vmin=vmin,vmax=vmax,interpolation='nearest')
 plt.colorbar()
 plt.subplot(133)
 ax.m = Basemap(projection='cyl',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='c')
 lons,lats = ax.m.makegrid(nlon,nlat)
 x,y = ax.m(lons,lats)
 cs = ax.m.imshow(data_corrected[0,:,:],vmin=vmin,vmax=vmax,interpolation='nearest')
 #plt.imshow(data_corrected[0,:,:],vmin=vmin,vmax=vmax)
 plt.colorbar()
 plt.show()
 return

def Add_Semivariogram_Info(vg_clim,rtime):
 #Add semivariogram information to the climatology
 if vg_clim[rtime.month]['svg'] != 0:
  vg = vg_clim[rtime.month]['svg']
  vg_clim[rtime.month]['count'] = vg_clim[rtime.month]['count'] + 1.0
  r.assign("vold",vg)
  r("vnew<-rbind(vold,vnew)")
  #w2 = 1/vg_clim[rtime.month]['count']
  #w1 = 1-w2
  #r('vnew["gamma"] = ' + str(w1) + '*vold["gamma"] + ' + str(w2) + '*vnew["gamma"]')
 #r("vold<-vnew")
 vg = r("vnew")
 vg_clim[rtime.month]['svg'] = vg
 return vg_clim

def Compute_Sample_Semivariogram(data_measurements):
 
 #Compute the sample variogram
 r.assign('x',data_measurements["x"])
 r.assign('y',data_measurements["y"])
 r.assign('departure',data_measurements["departure"])
 r("d = data.frame(x=x,y=y,departures=departure)")
 r("coordinates(d)<-~x+y")
 r('g<-gstat(id="departure",formula=departure~1,data=d)')
 #boundaries = np.linspace(0.0,750.0,20)
 #r.assign('boundaries',boundaries)
 r("vnew<-variogram(g,cutoff=700.0)")#,boundaries=boundaries)")#,cutoff=1000.0,width=750.0/20.0)")
 return

def Extract_Station_Data(st_coords,data,bdata,minlat,minlon,res,undef,var):
 #Extract info for the stations

 data_measurements = {}
 data_measurements["data"] = []
 data_measurements["lon"] = []
 data_measurements["lat"] = []
 for istation in xrange(0,len(data["id"])):
  idx = np.where(st_coords["id"] == data["id"][istation])
  if len(idx[0]) > 0:# and data[var][istation] >= undef:
   data_measurements["data"].append(data[var][istation])
   data_measurements["lon"].append(st_coords["lon"][idx[0]])
   data_measurements["lat"].append(st_coords["lat"][idx[0]])
 #Convert data to the units in NLDAS2
 data_measurements["data"] = np.array(data_measurements["data"])
 idx_valid = np.where(data_measurements["data"] > undef)
 #data_measurements["data"] = Convert_Units(data_measurements["data"],var)
 #data_measurements["data"] = 5.0/9.0*(np.array(data_measurements["data"]) - 32.0)
 data_measurements["lon"] = np.array(data_measurements["lon"])
 data_measurements["lat"] = np.array(data_measurements["lat"])
 #Extract only the valid cells
 data_measurements["data"] = data_measurements["data"][idx_valid]
 data_measurements["lon"] = data_measurements["lon"][idx_valid]
 data_measurements["lat"] = data_measurements["lat"][idx_valid]

 #Departures
 #1. Find the collocated cells and the corresponding departure
 idx_ccells = []
 nstations = len(data_measurements["lon"])
 data_measurements["ilat"] = np.zeros((nstations,))
 data_measurements["ilon"] = np.zeros((nstations,))
 data_measurements["bdata"] = np.zeros((nstations,))
 data_measurements["departure"] = np.zeros((nstations,))
 data_measurements["x"] = np.zeros((nstations,))
 data_measurements["y"] = np.zeros((nstations,))

 p = pyproj.Proj(proj='utm',zone=14,ellps='WGS84')
 for istation in xrange(0,nstations):
  ilat = int(np.round((data_measurements["lat"][istation] - minlat)/res))
  ilon = int(np.round((data_measurements["lon"][istation] - minlon)/res))
  data_measurements["ilat"][istation] = ilat
  data_measurements["ilon"][istation] = ilon
  data_measurements["bdata"][istation] = bdata[ilat,ilon]
  data_measurements["departure"][istation] = data_measurements["data"][istation] - data_measurements["bdata"][istation]
  x,y = p(data_measurements["lon"][istation],data_measurements["lat"][istation])
  data_measurements["x"][istation] = x/1000.0
  data_measurements["y"][istation] = y/1000.0

 return data_measurements

def Assimilate_Measurements_IDW(data_measurements,grid_background,nlat,nlon,undef,res,minlat,minlon):

 data_background = {}
 data_background["data"] = []
 data_background["lon"] = []
 data_background["lat"] = []
 data_background["ilon"] = []
 data_background["ilat"] = []
 data_background["x"] = []
 data_background["y"] = []
 data_background["departure"] = []

 #Convert the backgr2und data into a 1-d array and memorize the position on the grid
 p = pyproj.Proj(proj='utm',zone=14,ellps='WGS84')
 for i in xrange(0,nlat):
  for j in xrange(0,nlon):
   if grid_background[i,j] > undef:
    lat = minlat + i*res
    lon = minlon + j*res
    data_background["data"].append(grid_background[i,j])
    data_background["lon"].append(lon)
    data_background["lat"].append(lat)
    data_background["ilon"].append(j)
    data_background["ilat"].append(i)
    x,y = p(lon,lat)
    data_background["x"].append(x/1000.0)
    data_background["y"].append(y/1000.0)
 #Project coordinates
 data_background["lon"] = np.array(data_background["lon"])
 data_background["lat"] = np.array(data_background["lat"])
 data_background["x"] = np.array(data_background["x"])
 data_background["y"] = np.array(data_background["y"])

 #Interpolate the departures using IDW
 nbackground = len(data_background["ilat"])
 nstations = len(data_measurements["ilat"])
 ncells = nbackground + nstations
 x = data_background["x"]
 y = data_background["y"]
 p = 2.0
 xa = np.zeros((nbackground))
 d = np.zeros((nbackground))
 for i in xrange(0,nbackground):
   h = ((data_measurements["x"][:]-x[i])**2 + (data_measurements["y"][:]-y[i])**2)**0.5
   w = h**(-p)/np.sum(h**(-p))
   d[i] = np.sum(w*data_measurements["departure"])
   xa[i] = data_background["data"][i] + d[i]

 #Place the updated data into the matrix for the region
 data_corrected = np.copy(grid_background)
 dep = np.copy(grid_background)
 for i in xrange(0,nbackground):
  ilat = data_background["ilat"][i]
  ilon = data_background["ilon"][i]
  data_corrected[ilat,ilon] = xa[i]
  dep[ilat,ilon] = d[i] 

 return [data_corrected,dep]

def Assimilate_Measurements(data_measurements,grid_background,vg_fit,nlat,nlon,undef,res,minlat,minlon,bdata2):
 
 data_background = {}
 data_background["data"] = []
 data_background["lon"] = []
 data_background["lat"] = []
 data_background["ilon"] = []
 data_background["ilat"] = []
 data_background["x"] = []
 data_background["y"] = []

 #Convert the background data into a 1-d array and memorize the position on the grid
 import time
 tic = time.clock()
 p = pyproj.Proj(proj='utm',zone=14,ellps='WGS84')
 for i in xrange(0,nlat):
  for j in xrange(0,nlon):
   if grid_background[i,j] > undef:
    lat = minlat + i*res
    lon = minlon + j*res
    data_background["data"].append(grid_background[i,j])
    data_background["lon"].append(lon)
    data_background["lat"].append(lat)
    data_background["ilon"].append(j)
    data_background["ilat"].append(i)
    x,y = p(lon,lat)
    data_background["x"].append(x/1000.0)
    data_background["y"].append(y/1000.0)
 data_background["lon"] = np.array(data_background["lon"])
 data_background["lat"] = np.array(data_background["lat"])
 data_background["x"] = np.array(data_background["x"])
 data_background["y"] = np.array(data_background["y"])

 #2.1 Create the distane matrix
 nbackground = len(data_background["ilat"])
 nstations = len(data_measurements["ilat"])
 ncells = nbackground + nstations
 distances = np.zeros((ncells,ncells))
 x = np.concatenate((data_measurements["x"][:],data_background["x"]),axis=0)
 y = np.concatenate((data_measurements["y"][:],data_background["y"]),axis=0)
 for i in xrange(0,len(x)):
   distances[i,:] = ((x-x[i])**2 + (y-y[i])**2)**0.5

 #2.2 Create the covariance matrix
 C = calculate_covariance(vg_fit,distances)

 #2.3 Construct the measurement matrix
 H = np.zeros((nstations,ncells))
 for i in xrange(0,nstations):
  H[i,i] = 1.0

 #2.3 Calculate the Kalman Gain
 CHt = np.dot(C,H.T)
 HCHt = np.dot(np.dot(H,C),H.T)
 #U,s,Vh = la.svd(HCHt,full_matrices=False)
 #singular = s < 10**-1
 #invS = 1/s
 #invS[singular] = 0.0
 #iHCHt = np.dot(np.dot(U,np.diag(invS)),Vh)
 iHCHt = la.pinv(HCHt)
 K = np.dot(CHt,iHCHt)

 #2.4 update the background
 xb = np.concatenate((data_measurements["bdata"],data_background["data"]),axis=0)
 d = data_measurements["departure"] #data_measurements["data"] - np.dot(h,xb)
 d = np.dot(K,d)
 xa = xb + d

 #2.5 Place the updated data into the matrix for the region
 data_corrected = np.copy(grid_background)  
 dep = np.copy(grid_background)
 for i in xrange(nstations,ncells):
  ilat = data_background["ilat"][i-nstations]
  ilon = data_background["ilon"][i-nstations]
  data_corrected[ilat,ilon] = xa[i]
  dep[ilat,ilon] = d[i]

 return [data_corrected,dep]

def Assimilate_Measurements_OK(data_measurements,grid_background,vg_fit,nlat,nlon,undef,res,minlat,minlon,bdata2):

 data_background = {}
 data_background["data"] = []
 data_background["lon"] = []
 data_background["lat"] = []
 data_background["ilon"] = []
 data_background["ilat"] = []
 data_background["x"] = []
 data_background["y"] = []

 #Convert the background data into a 1-d array and memorize the position on the grid
 p = pyproj.Proj(proj='utm',zone=14,ellps='WGS84')
 for i in xrange(0,nlat):
  for j in xrange(0,nlon):
   if grid_background[i,j] > undef:
    lat = minlat + i*res
    lon = minlon + j*res
    data_background["data"].append(grid_background[i,j])
    data_background["lon"].append(lon)
    data_background["lat"].append(lat)
    data_background["ilon"].append(j)
    data_background["ilat"].append(i)
    x,y = p(lon,lat)
    data_background["x"].append(x/1000.0)
    data_background["y"].append(y/1000.0)
 data_background["lon"] = np.array(data_background["lon"])
 data_background["lat"] = np.array(data_background["lat"])
 data_background["x"] = np.array(data_background["x"])
 data_background["y"] = np.array(data_background["y"])

 #2.1 Create the distane matrix
 nbackground = len(data_background["ilat"])
 nstations = len(data_measurements["ilat"])
 ncells = nbackground + nstations
 distances = np.zeros((ncells-1,ncells-1))
 distancesiN = np.zeros((ncells-1,ncells-1))
 distancesjN = np.zeros((ncells-1,ncells-1))
 x = np.concatenate((data_measurements["x"][:-1],data_background["x"]),axis=0)
 y = np.concatenate((data_measurements["y"][:-1],data_background["y"]),axis=0)
 for i in xrange(0,len(x)):
   distances[i,:] = ((x-x[i])**2 + (y-y[i])**2)**0.5
   distancesiN[i,:] = ((x-data_measurements["x"][-1])**2 + (y-data_measurements["y"][-1])**2)**0.5
   distancesjN[:,i] = ((x-data_measurements["x"][-1])**2 + (y-data_measurements["y"][-1])**2)**0.5

 #2.2 Create the generalized covariance matrix
 C = calculate_generalized_covariance(vg_fit,distances,distancesiN,distancesjN)

 #2.3 Construct the measurement matrix
 H = np.zeros((nstations-1,ncells-1))
 for i in xrange(0,nstations-1):
  H[i,i] = 1.0

 #2.3 Calculate the Kalman Gain
 CHt = np.dot(C,H.T)
 HCHt = np.dot(np.dot(H,C),H.T)
 #U,s,Vh = la.svd(HCHt,full_matrices=False)
 #singular = s < 10**-1
 #invS = 1/s
 #invS[singular] = 0.0
 #iHCHt = np.dot(np.dot(U,np.diag(invS)),Vh)
 iHCHt = la.pinv(HCHt)
 K = np.dot(CHt,iHCHt)

 #2.5 Place the updated data into the matrix for the region
 data_corrected = np.copy(grid_background)
 dep = np.copy(grid_background)
 for i in xrange(0,nbackground): 
  ilat = data_background["ilat"][i]
  ilon = data_background["ilon"][i]
  weights = K[i+nstations-1,:]
  dep[ilat,ilon] = np.dot(weights,data_measurements["departure"][:-1]) + (1-np.sum(weights))*data_measurements["departure"][-1]
  data_corrected[ilat,ilon] = data_background["data"][i] + dep[ilat,ilon]

 return [data_corrected,dep]

def calculate_generalized_covariance(vg_fit,h,hiN,hjN):
 #Create the generalzied covariance matrix Cij = gamma(i,n)+gamma(j,n)-gamma(i,j)
 p1 = calculate_semivariogram(vg_fit,hiN)#calculate_semivariogram(a,c0,hiN)
 p2 = calculate_semivariogram(vg_fit,hjN)#calculate_semivariogram(a,c0,hjN)
 p3 = calculate_semivariogram(vg_fit,h)#calculate_semivariogram(a,c0,h)
 C = p1 + p2 - p3
 return C
 
def calculate_covariance(vg_fit,h):
  #C = np.zeros((h.shape[0],h.shape[1]))
  #C = np.exp(-h**2/a**2)
  #idx = np.where(h < a)
  #C[idx] = 1.0 - (3.0/2.0)*h[idx]/a + (1.0/2.0)*(h[idx]/a)**3
  #idx = np.where(h >= a)
  #C[idx] = 0.0
  C = vg_fit[1][0] - calculate_semivariogram(vg_fit,h)
  return C

#def calculate_semivariogram(a,c0,h):
def calculate_semivariogram(vg_fit,h):
  #dist_vector = np.reshape(h,h.shape[0]*h.shape[1])
  #r.assign("dist_vector",dist_vector)
  #r.assign("v.fit",vg_fit)
  #vg_est = r("variogramLine(v.fit,dist_vector=dist_vector)")
  #V = np.reshape(np.array(vg_est[1]),(h.shape[0],h.shape[1]))
  #Exponential
  V = vg_fit[1][0]*(1 - np.exp(-h/vg_fit[2][0]))
  #Spherical
  #V = np.zeros((h.shape[0],))
  #V = np.copy(h)
  #idx = np.where(h <= a)
  #V[idx] = c0*((3.0/2.0)*(h[idx]/a) - (1.0/2.0)*(h[idx]/a)**3)
  #idx = np.where(h > a)
  #V[idx] = c0
  #V = c0*(1 - np.exp(-h**2/a**2))
  return V

def Derive_Semivariogram_Climatology(icount,fcount,tinitial,dt,itime,dtime,station_data,data_output,st_coords,minlat,minlon,res,undef,maxlat,maxlon,var,nstations):

 new_network = {}

 #Constructe climatology array
 vg_clim = {}
 for i in xrange(0,12):
  vg_clim[i+1] = {}
  vg_clim[i+1]['svg'] = 0
  vg_clim[i+1]['count'] = 0.0

 #Derive the semivariogram climatology
 for count in xrange(icount,fcount):
  i = count - icount
  rtime = itime + count*dtime
  t1 = tinitial + count*dt
  t2 = tinitial + count*dt + 23
  print rtime.year,rtime.month,rtime.day

  #Define a new network
  new_network['sample'] = random.sample(range(len(st_coords['id'])),nstations)

  #Extract the info for these stations
  new_network['st_coords'] = st_coords[new_network['sample']]

  #Read in gridded data
  ga("data = ave(" + var + ".2,t="+str(t1)+",t="+str(t2)+")")
  bdata = ga.exp("data")

  #Extract station data
  idx = np.where((station_data["year"] == rtime.year) & (station_data["month"] == rtime.month) & (station_data["day"] == rtime.day))
  data = station_data[idx]

  #Extract all the info for the stations
  data_measurements = Extract_Station_Data(new_network['st_coords'],data,bdata,minlat,minlon,res,undef,var)

  #Compute the sample semivariogram
  Compute_Sample_Semivariogram(data_measurements)

  #Add semivariogram information to the climatology
  vg_clim = Add_Semivariogram_Info(vg_clim,rtime)
  vg = vg_clim[rtime.month]['svg']

  #Pass time step output to the outgoing data
  data_output["vg"].append(vg)
  data_output["time"].append(rtime)

 #Fit a model semivariogram to each month
 Fit_Semivariogram_Climatology(vg_clim,data_output)

 return data_output

def Create_Data_Output_Variable():
 data_output = {}
 data_output["corrected"] = []
 data_output["vgfit"] = []
 data_output["vg"] = []
 data_output["time"] = []
 data_output["vg_clim"] = {}
 data_output["original"] = []
 return data_output

def Convert_Units(data,var):
 #Convert units of the data
 if var == 'tair':
  data = 5.0/9.0*(np.array(data) - 32.0) #Average C
 if var == 'rh': 
  data = np.array(data) #Average %
 if var == 'wind':
  data = np.array(data)*1609.34/3600.0 #Average m/s
 if var == 'prec':
  data = np.array(data)*25.4/24.0 #Avearege m/s
 if var == 'pres':
  data = np.array(data)*3386.0/100 #Average hPa
 if var == 'apcpsfc':
  data = np.array(data)*25.4 #Avearege mm/hr
 return data

def Convert_NLDAS(data,var):
 #Convert units of the data
 if var == 'prec':
  data = np.array(data)*3600.0*1000.0 #Avearege m/s
 return data

def Covariance_Matrix(data_measurements,grid_background,vg_fit,nlat,nlon,undef,res,minlat,minlon,bdata2):

 data_background = {}
 data_background["data"] = []
 data_background["lon"] = []
 data_background["lat"] = []
 data_background["ilon"] = []
 data_background["ilat"] = []
 data_background["x"] = []
 data_background["y"] = []

 #Convert the background data into a 1-d array and memorize the position on the grid
 p = pyproj.Proj(proj='utm',zone=14,ellps='WGS84')
 for i in xrange(0,nlat):
  for j in xrange(0,nlon):
   if grid_background[i,j] > undef:
    lat = minlat + i*res
    lon = minlon + j*res
    data_background["data"].append(grid_background[i,j])
    data_background["lon"].append(lon)
    data_background["lat"].append(lat)
    data_background["ilon"].append(j)
    data_background["ilat"].append(i)
    x,y = p(lon,lat)
    data_background["x"].append(x/1000.0)
    data_background["y"].append(y/1000.0)
 data_background["lon"] = np.array(data_background["lon"])
 data_background["lat"] = np.array(data_background["lat"])
 data_background["x"] = np.array(data_background["x"])
 data_background["y"] = np.array(data_background["y"])

 #2.1 Create the distane matrix
 nbackground = len(data_background["ilat"])
 nstations = len(data_measurements["ilat"])
 ncells = nbackground
 distances = np.zeros((ncells,ncells))
 x = data_background["x"]
 y = data_background["y"]
 for i in xrange(0,len(x)):
   distances[i,:] = ((x-x[i])**2 + (y-y[i])**2)**0.5

 #2.2 Create the covariance matrix
 C = calculate_covariance(vg_fit,distances)

 return C

def retrieve_metadata(raster):

 metadata = {}
 #Extract coordinates and projection
 ds = gdal.Open(raster)
 gt = ds.GetGeoTransform()
 cols = ds.RasterXSize
 rows = ds.RasterYSize
 srs = osgeo.osr.SpatialReference()
 srs.ImportFromWkt(ds.GetProjection())
 metadata['proj4'] = srs.ExportToProj4()
 metadata['minx'] = gt[0]
 metadata['miny'] = gt[3]+rows*gt[5]
 metadata['maxx'] = gt[0]+cols*gt[1]
 metadata['maxy'] = gt[3]
 metadata['resx'] = gt[1]
 metadata['resy'] = gt[5]
 metadata['gt'] = gt
 metadata['nx'] = cols
 metadata['ny'] = rows
 metadata['projection'] = ds.GetProjection()

 return metadata

def read_raster(file):

 #Read in the raster
 dataset = gdal.Open(file)

 #Get dimensons
 nx = dataset.RasterXSize
 ny = dataset.RasterYSize

 #Retrieve band
 band = dataset.GetRasterBand(1)

 #Convert to numpy array
 data = band.ReadAsArray(0,0,nx,ny).astype(np.float32)

 return data

def write_raster(data,metadata,file):

 #Create file
 driver = gdal.GetDriverByName('GTiff')
 ds = driver.Create(file,metadata['nx'],metadata['ny'],1,gdal.GDT_Float32)

 #Set geo information
 ds.SetGeoTransform(metadata['gt'])
 ds.SetProjection(metadata['projection'])
 outband = ds.GetRasterBand(1)
 outband.WriteArray(data,0,0)
 ds = None

 return 
