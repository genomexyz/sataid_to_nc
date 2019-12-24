#!/usr/bin/python3

from datetime import datetime,timedelta
from struct import unpack
from numpy import array,asarray,arange,flipud,dtype
from netCDF4 import Dataset,date2num
from sys import argv
from os import popen


def readZ(pathZ):
	fi = open(pathZ,'rb')
	recl = unpack('I',fi.read(4))
	chan = unpack('c'*8,fi.read(8))	 #Channel Name (8)
	sate = unpack('c'*8,fi.read(8))	 #Satellite Name (8)
	skip = unpack('I'*1,fi.read(4*1))
	ftim = unpack('I'*8,fi.read(4*8))   #Start Time Filming (8)
	etim = unpack('I'*8,fi.read(4*8))   #End Time Filming (8)
	calb = unpack('I'*1,fi.read(4*1))   #Coordinate Convertion (1)
	fint = unpack('I'*2,fi.read(4*2))   #Resolution before Convertion (2)
	eres = unpack('f'*2,fi.read(4*2))   #Interval Lon Lat (2)
	eint = unpack('I'*2,fi.read(4*2))   #Grid Lon Lat (2)
	nrec = unpack('I'*2,fi.read(4*2))   #Number Records and Byte Length (2)
	cord = unpack('f'*8,fi.read(4*8))   #lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4 (8)
	ncal = unpack('I'*3,fi.read(4*3))   #Number, start, end of Calibration (3)
	skip = unpack('c'*24,fi.read(1*24))
	asat = unpack('f'*6,fi.read(4*6))   #Roll, Pitch, Yaw, Lat, Lon, Altitude of Satellite (6)
	skip = unpack('c'*32,fi.read(1*32))
	vers = unpack('c'*4,fi.read(1*4))   #Sataid version
	recl = unpack('I'*1,fi.read(4*1))
	nbyt = unpack('I'*1,fi.read(4*1))
	cal = array(unpack('f'*(nbyt[0]//4-2),fi.read(4*(nbyt[0]//4-2))))   #Calibration Table
	nbyt = unpack('I'*1,fi.read(4*1))
	data = []
	if nrec[1]==2:
		#HimawariCloud 2 Byte, HimawariWIS 2 Byte
		for i in range(eint[1]):
			nbyt = unpack('I'*1,fi.read(4*1))
			line = unpack('H'*(eint[0]),fi.read(eint[0]*2))
			data.append(line[0:eint[0]])
			unpack('B'*(nbyt[0]-eint[0]*2-8),fi.read(nbyt[0]-eint[0]*2-8))
			nbyt = unpack('I'*1,fi.read(4*1))
	elif nrec[1]==1:
		#HimawariCloud 1 Byte
		for i in range(eint[1]):
			nbyt = unpack('I'*1,fi.read(4*1))
			line = unpack('B'*((nbyt[0]-8)),fi.read(((nbyt[0]-8))))
			data.append(line[0:eint[0]])
			nbyt = unpack('I'*1,fi.read(4*1))
	cal = asarray(cal)
	fi.close()
	return sate,chan,ftim,etim,eres,eint,cord,asat,cal,data

def readZWIS(pathZ):
	popen('expand '+pathZ+' '+pathZ+'.dat')
	fi = open(pathZ+'.dat','rb')
	recl = unpack('I',fi.read(4))	   
	chan = unpack('c'*8,fi.read(8))	 #Channel Name (8)
	sate = unpack('c'*8,fi.read(8))	 #Satellite Name (8)
	skip = unpack('I'*1,fi.read(4*1))
	ftim = unpack('I'*8,fi.read(4*8))   #Start Time Filming (8)
	etim = unpack('I'*8,fi.read(4*8))   #End Time Filming (8)
	calb = unpack('I'*1,fi.read(4*1))   #Coordinate Convertion (1)
	fint = unpack('I'*2,fi.read(4*2))   #Resolution before Convertion (2)
	eres = unpack('f'*2,fi.read(4*2))   #Interval Lon Lat (2)
	eint = unpack('I'*2,fi.read(4*2))   #Grid Lon Lat (2)
	nrec = unpack('I'*2,fi.read(4*2))   #Number Records and Byte Length (2)
	cord = unpack('f'*8,fi.read(4*8))   #lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4 (8)
	ncal = unpack('I'*3,fi.read(4*3))   #Number, start, end of Calibration (3)
	skip = unpack('c'*24,fi.read(1*24)) 
	asat = unpack('f'*6,fi.read(4*6))   #Roll, Pitch, Yaw, Lat, Lon, Altitude of Satellite (6)
	skip = unpack('c'*32,fi.read(1*32)) 
	vers = unpack('c'*4,fi.read(1*4))   #Sataid version
	recl = unpack('I'*1,fi.read(4*1))   
	nbyt = unpack('I'*1,fi.read(4*1))
	cal = array(unpack('f'*1022,fi.read(4*1022)))   #Calibration Table
	nbyt = unpack('I'*1,fi.read(4*1))
	data = []
	for i in range(eint[1]):
		nbyt = unpack('I'*1,fi.read(4*1))
		data.append(unpack('H'*eint[0],fi.read(2*eint[0])))
		unpack('H'*((int(nbyt[0])-int(eint[0])*2-8)/2),fi.read(int(nbyt[0])-int(eint[0])*2-8))
		nbyt = unpack('I'*1,fi.read(4*1))
	data = array(data)
	fi.close()
	popen('del '+pathZ+'.dat')
	return sate,chan,ftim,etim,eres,eint,cord,asat,cal,data

def calibrating(lut,data):
	data = [lut[i-1] for i in data]
	return data

def cropData(data,lats,lons,ul,dl,ll,rl):
	res = lats[1]-lats[0]
	uld = int((ul-lats[0])/res)
	lld = int((ll-lons[0])/res)
	dld = int((dl-lats[0])/res)
	rld = int((rl-lons[0])/res)
	data = data[dld+1:uld+2,lld:rld+1]
	return data

def convert(fsataid,fnetcdf):
	try:
		sate,chan,ftim,etim,eres,eint,cord,asat,cal,data = readZ(fsataid)
	except:
		sate,chan,ftim,etim,eres,eint,cord,asat,cal,data = readZWIS(fsataid)
	data = flipud(asarray(data))
	ul = cord[4]
	dl = cord[0]
	ll = cord[1]
	rl = cord[3]
	res = eres[0]
	lats_data = arange(cord[4],cord[0]-res,res)
	lons_data = arange(cord[1],cord[3]-res,res)
	data = data[0:len(lats_data),0:len(lons_data)]
	datatbb = calibrating(cal,data)
	ftime = datetime(ftim[0]*100+ftim[1],ftim[2],ftim[3],ftim[4],ftim[5])
	time = datetime(ftim[0]*100+ftim[1],ftim[2],ftim[3],ftim[4],int(ftim[5]/10*10))+timedelta(minutes=10)
	file = Dataset(fnetcdf,'w')
	file.title = ''.join('%c' %c.decode() for c in chan[0:2])+' Indonesia '+time.strftime('%Y%m%d-%H%M UTC')
	file.institution = 'Badan Meteorologi Klimatologi dan Geofisika'
	file.createDimension('calibration',len(cal))
	file.createDimension('latitude',len(lats_data))
	file.createDimension('longitude',len(lons_data))
	file.createDimension('time',1)
	calc = file.createVariable('calibration',dtype('f4').char,('calibration',))
	lats = file.createVariable('latitude',dtype('f4').char,('latitude',))
	lons = file.createVariable('longitude',dtype('f4').char,('longitude',))
	times = file.createVariable('time',dtype('f8').char,('time',))
	lats.units = 'degrees_north'
	lons.units = 'degrees_east'
	times.units = 'minutes since 0001-01-01 00:01:00.0'
	times.calendar = 'gregorian'
	calc[:] = cal
	times[:] = date2num(ftime,units=times.units,calendar=times.calendar)
	lons[:] = lons_data
	lats[:] = lats_data
	data_out = file.createVariable(chan[0]+chan[1],dtype('f4').char,('time','latitude','longitude'),zlib=True,least_significant_digit=3)
	data_out[:] = datatbb
	file.close()
	return fnetcdf

def main(fsataid):	
	fnetcdf = fsataid+'.nc'
	convert(fsataid,fnetcdf)

if len(argv)==2:
	fname = argv[1]
	try:
		main(fname)
	except:
		print('Not working : wrong format file')
else:
	print('Not working')
	print('Usage : sataid2nc.exe inputfile')
	print('created by Andersen 2016')
	print('Revised (to make it work in python3) by Muhammad Ryan 2019')
