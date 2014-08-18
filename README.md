This software implements the state-space estimation procedure described in Chaney et al.,2014

You will need to install a number of python libraries before it will work. Below are the steps to follow to install the libraries. All the following steps are done at the command line (linux).

0. You should be logged in as root. Make sure that your python version is 2.7. If you don't then install it. You can find the software package here "https://www.python.org/download/releases/2.7/"

1. Install pip (http://pip.readthedocs.org/en/latest/installing.html). First download the get-pip.py file then run: "python2.7 get-pip.py"

2. Install numpy "pip install numpy"

3. Install pyproj "pip install pyproj"

4. Install rpy2 "pip instal rpy2" (Make sure that R 3.0.2 or higher is installed)

5. Install gdal. This might take a bit. Download this package ftp://ftp.remotesensing.org/gdal/1.11.0/gdal-1.11.0.tar.gz.
   Untar, configure, make, and make install. This will most likely be the most challenging part. 

6. Install python gdal library "pip install gdal"

Once you have installed all the libraries, you can run python2.7 driver.py. This will run the example and creates the corrected.tif file which you can read into other GIS packages.

NOTE: If when running the program you get:

ImportError: /usr/local/lib/python2.7/site-packages/osgeo/_gdal.so: undefined symbol: GDALRasterBandGetVirtualMem

Then just type this:

export LD_PRELOAD=/usr/local/lib/libgdal.so.1

