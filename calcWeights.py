#!/usr/bin/python

######################################
# Loading packages 		     #
######################################




from pupynere import netcdf_file
import numpy as np
import pdb
import matplotlib.dates
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import math as m
#from scipy.interpolate import RegularGridInterpolator
from scipy import interpolate

import os, sys

sys.path.append('/Library/Python/2.7/site-packages/pupynere-1.0.15')
sys.path.append('/Library/Python/2.7/site-packages/PyNGL-1.4.0.macos-10.6-x86_64-py271-numpy160/lib/python2.7/site-packages/PyNGL')
os.environ["PYNGL_NCARG"] = "/Library/Python/2.7/site-packages/PyNGL-1.4.0.macos-10.6-x86_64-py271-numpy160/lib/python2.7/site-packages/PyNGL/ncarg"


import Ngl
from eofs.standard import Eof
from pylab import demean
import math
import time

def calcLon(x,y):
    '''Converts from x,y cartesian to longitude.
    x- x coordinate of cartesian system
    y- y coordinate of cartesian system
    
    returns: longitude
    '''
    rad = m.atan2(y,x)
    lon = rad * (180/m.pi)
    return lon

def calcLat(z):
    '''Converts from z cartesian to latitude.
    z- z coordinate of cartesian system
    
    returns: latitude
    '''
    rad = m.asin(z/6371000)
    lat = rad * (180/m.pi)
    return lat

    
def calcDist(gTot, n1, n2):
    '''Finds great circle distance in nmi
    gTot: graph used to find shortest path
    n1: start node of edge
    n2: end node of edge
    
    returns: great circle distance between n1 and n2
    '''
    r = m.pi/180
    delSigma = m.acos(m.sin(gTot.node[n1]['lat']*r)*m.sin(gTot.node[n2]['lat']*r) + m.cos(gTot.node[n1]['lat']*r)*m.cos(gTot.node[n2]['lat']*r)*m.cos((gTot.node[n2]['lon']*r)-(gTot.node[n1]['lon']*r)))
    dist = (6371000 * delSigma)/1850
    return dist    

def calcDistPreproc(lat1, lon1, lat2, lon2):
    '''Finds great circle distance in nmi before node properties are added to graph
    lat1: latitude of start node in edge
    lon1: longitude of start node in edge
    lat2: latitude of target node in edge
    lon2: longitude of target node in edge
    
    returns: great circle distance between a start and end node
    '''    
    r = m.pi/180
    delSigma = m.acos(m.sin(lat1*r)*m.sin(lat2*r) + m.cos(lat1*r)*m.cos(lat2*r)*m.cos((lon2*r)-(lon1*r)))
    dist = (6371000 * delSigma)/1850
    return dist 




def calcEdgeBearing(gTot, n1, n2): 
    '''relevant url for commented-out COG formula: http://www.movable-type.co.uk/scripts/latlong.html'''
    '''Calculates edge bearing angle between two nodes n1 and n2
    gTot: graph used to find shortest path
    n1: start node of edge
    n2: end node of edge
    
    returns: psi, the edge bearing angle
    '''
    
    deg2rad= m.pi/180    
    lon1= deg2rad * gTot.node[n1]['lon']     
    lon2= deg2rad * gTot.node[n2]['lon']
    lat1= deg2rad * gTot.node[n1]['lat']     
    lat2= deg2rad * gTot.node[n2]['lat']    
 
    psi = calcCOG1(lat1, lon1, lat2, lon2)
 
    #psi= psi / deg2rad
    return psi


def calcCOG1(lat1, lon1, lat2, lon2):
    '''Calculates COG between two nodes
    lat1: latitude of start node in edge
    lon1: longitude of start node in edge
    lat2: latitude of target node in edge
    lon2: longitude of target node in edge
    
    returns: psi, the edge bearing angle
    '''
       
    epsilon = m.atan2( lat2-lat1,   lon2-lon1 ) 
    psi = m.pi/2 - epsilon 
    if psi<0: 
        psi = psi + 2* m.pi       
       
    return psi


def calcFlow(psi, uFlow, vFlow):
    '''Calculates parallel and perpendicular components of current
    psi: COG angle
    uFlow: u component of current flow
    vFlow: v component of current flow
    
    returns:
    wParallel: parallel component of w vector (current)
    wPerp: perpendicular component of w vector (current)
    '''
    sP= m.sin(psi)
    cP= m.cos(psi)
    wParallel = uFlow * sP + vFlow * cP
    wPerp     = vFlow * sP - uFlow * cP
    return wParallel, wPerp

def calcSOG(gTot, n1, n2, vesselSpeedMax, fuDict, fvDict):
    '''Calculates the speed over ground SOG between two nodes
    gTot: graph used to find shortest path
    n1: start node ID
    n2: end node ID
    vesselSpeedMax: user specified maximum vessel speed
    fuDict: dictionary, key = node, value = u component of current interpolated at the node
    fvDict: dictionary, key = node, value = u component of current interpolated at the node
    
    returns: SOG
    '''
    psi = calcEdgeBearing(gTot, n1, n2)
    uFlow1 = fuDict[n1]
    vFlow1 = fvDict[n1]
    uFlow2 = fuDict[n2]
    vFlow2 = fvDict[n2]
    wParallel1, wPerp1 = calcFlow(psi, uFlow1, vFlow1)
    wParallel2, wPerp2 = calcFlow(psi, uFlow2, vFlow2)
    wParallel = (wParallel1 + wParallel2)/2 *1.94384 #kts
    wPerp = (wParallel1 + wParallel2)/2 *1.94384 #kts
    w = np.sqrt(wPerp**2 + wParallel**2) #kts
    delta = m.asin(wPerp/vesselSpeedMax)
    heading = psi + delta
    sog = wParallel + m.sqrt(vesselSpeedMax*vesselSpeedMax - wPerp*wPerp)
    return sog

def calcEdgeWeight(gTot, n1, n2, vesselSpeedMax, fuDict, fvDict):
    '''Calculates edge weight based on current magnitude data
    gTot: graph used to find shortest path
    n1: start node ID
    n2: end node ID
    vesselSpeedMax: user specified maximum vessel speed
    fuDict: dictionary, key = node, value = u component of current interpolated at the node
    fvDict: dictionary, key = node, value = u component of current interpolated at the node
    
    returns: edge weight based on current magnitudes
    '''
    dist = calcDist(gTot, n1, n2)
    
    '''code for time dependence. feed calcSOG the fuDict/fvDict for one individual timestep, get the weight value, and add to list where
    list index corresponds to a timestep. uncomment extra code to enable this function.'''
    dt_list  = []
    nTsteps = 1 #change
    for tstep in range(nTsteps): #depending on how you found the number of timesteps, this may need to be range(nTsteps - 1) or something
        sog = calcSOG(gTot, n1, n2, vesselSpeedMax, fuDict, fvDict)
        dt = dist/sog
        #dt_list.append(dt)
        
    #return dt_list
    return dt




def getFlowArrays(nameDict):
    '''Retrieves vectors from current data
    nameDict: csv file dictionary of user specifications/flags
    
    returns: list of information relevant for interpolating current field at each mesh node. 
    flowArrays[0] = list of current vector longitude coordinates
    flowArrays[1] = list of current vector latitude coordinates
    flowArrays[2] = list of current vector u component magnitudes
    flowArrays[3] = list of current vector v component magnitudes
    flowArrays[4] = list of vector magnitudes (module of u and v components) at each coordinate
    '''
    fU=netcdf_file(nameDict['current_db'][0],"r")
    
    iLon=fU.variables["lon"][:,:] # longitude
    iLat=fU.variables["lat"][:,:] # latitude
    
    '''in the future, write so that it looks at the database to find number of timesteps. I also don't know how to 
    access one timestep in the whole database, so I don't know how the for loop will work. The important thing is
    that once a single timestep of data is extracted, add it to a list where the index of the list corresponds to a timestep.'''
    nTsteps = 1 
    
    #iU_tot = []
    #iV_tot = []
    #iModule_tot = []
    
    for tstep in range(nTsteps): #depending on how you found the number of timesteps, this may need to be range(nTsteps - 1) or something
        
    
        iU=fU.variables["sozocrtx"][:,:,:] # u
        iV=fU.variables["somecrty"][:,:,:] # v 
        
        iU_ = np.copy(iU.squeeze())
        iV_ = np.copy(iV.squeeze())  
        
        iU_[np.where(iU_>100)]=np.nan
        iV_[np.where(iV_>100)]=np.nan    
        fU.close()
        
        iModule=np.sqrt(iU_**2 + iV_**2)*1.94384 # use built in from 2D vector magnitude; put conversion m/s->ft into a constant
        
        #iU_tot.append(iU_)
        #iV_tot.append(iV_)
        #iModule_tot.append(iModule)
    
    
    flowArrays = []
    flowArrays.append(iLon)
    flowArrays.append(iLat)
    flowArrays.append(iU_)
    flowArrays.append(iV_)    
    flowArrays.append(iModule)
    
    #flowArrays = []
    #flowArrays.append(iLon)
    #flowArrays.append(iLat)
    #flowArrays.append(iU_tot)
    #flowArrays.append(iV_tot)    
    #flowArrays.append(iModule_tot)    
 
    
    return flowArrays

def zipFunction(names, function):
    '''Associates a function value with the node that has the value.
    names: list of node IDs in graph
    function: interpolation function
    
    returns: dictionary. key=node ID, value= value of function evaluated at that node
    '''
    fDict = {}
    i = 0
    for name in names:
        fDict[name] = function[i]
        i = i + 1
    return fDict

def reshapeMags(mag,dim):
    '''Reshapes the ndarray of current magnitudes into correct format
    mag: nd array of vector magnitudes
    dim: number of data points (rows*columns of lat/lon arrays)
    
    returns: magnitude array with shape of dim rows, 1 column
    '''
    k = 0
    mArray = np.ndarray(shape=(dim,))
    for i in range(len(mag)):
        for j in range(len(mag[i])):
            mArray[k] = mag[i][j]
            k = k + 1
    return mArray

def createLatLonArray(lon,lat):
    '''Creates an array of lats and lons.
    lon: array of longitudes of shape: # rows = total # points, 1 column
    lat: array of latitudes of shape: # rows = total # points, 1 column
    
    returns: array with 2 columns (lat, lon). Number of rows = rows*columns of lat/lon arrays.
    '''
    dim = len(lon)
    llArray = np.ndarray(shape=(dim,2))
    for i in range(len(lon)):
        llArray[i][0] = lat[i]
        llArray[i][1] = lon[i]
    return llArray 

def createCurArray(iLon,iLat,dim): 
    '''Creates coordinate array from iLon, iLat arrays.
    iLon: longitude array taken directly from current .nc file. shape = array of arrays
    iLat: latitude array taken directly from current .nc file. shape = array of arrays
    dim: dimensions of iLon/iLat arrays. equal to # rows * # columns
    
    returns: array. number of rows = dim, two columns (lat, lon)
    '''
    k = 0
    cArray = np.ndarray(shape=(dim,2))
    for i in range(len(iLon)):
        for j in range(len(iLon[0])):
            cArray[k][0] = iLat[i][j]
            cArray[k][1] = iLon[i][j]
            k = k + 1
    #print cArray[:10]
    return cArray



def interpField(nProps, flowArrays, files): #not complete, development in progress
    '''Creates interpolation functions based on current flow arrays and creates dictionaries that associate an interpolated current magnitude value at each mesh node.
    nProps: node property list. contains list of mesh node longitudes, mesh node latitudes, and mesh node IDs
    flowArrays: list of lists.
    flowArrays[0] = list of current vector longitude coordinates
    flowArrays[1] = list of current vector latitude coordinates
    flowArrays[2] = list of current vector u component magnitudes
    flowArrays[3] = list of current vector v component magnitudes
    flowArrays[4] = list of vector magnitudes (module of u and v components) at each coordinate
    files: list of file input names and directories for data retreival. 
    
    returns:
    fuDict: dictionary of current u components interpolated at each mesh node. key=mesh node ID, value= interpolated u component magnitude
    fvDict: dictionary of current v components interpolated at each mesh node. key=mesh node ID, value= interpolated v component magnitude
    '''
    iLon = flowArrays[0]
    iLat = flowArrays[1]
    
    '''In the future, edit for loop (i.e. uncomment commented code and change any necessary syntax) so that this method will return a list
    of fuDicts/fvDicts for each timestep'''
    nTsteps = 1
    #numTimesteps = len(iU_tot)
    
    all_fuDicts = []
    all_fvDicts = []
    
    for tstep in range(nTsteps): #depending on how you found the number of timesteps, this may need to be range(nTsteps - 1) or something
    
        iU_ = flowArrays[2]#[tstep]
        iV_ = flowArrays[3]#[tstep]
        
        nLons = nProps[0]
        nLats = nProps[1]
        nNames = nProps[2]
    
        dim = ((iLon.shape[0]) * (iLon.shape[1]))
        
        llArray = createLatLonArray(nLons,nLats)
        cArray = createCurArray(iLon,iLat,dim)
        uArray = reshapeMags(iU_,dim)
        vArray = reshapeMags(iV_,dim)
    
        #ufunct = RegularGridInterpolator(cArray, uArray, llArray, method='linear') #still developing
        #vfunct = RegularGridInterpolator(cArray, vArray, llArray, method='linear') #still developing
        
        uFunct = interpolate.griddata(cArray, uArray, llArray, method='linear')
        vFunct = interpolate.griddata(cArray, vArray, llArray, method='linear')    
        
        #vizInterp(nLons, nLats, ufunct, vfunct, files) #still developing
        
        fuDict = zipFunction(nNames, uFunct)
        fvDict = zipFunction(nNames, vFunct)
        
        #all_fuDicts.append(fuDict)
        #all_fvDicts.append(fvDict)
        
    #return all_fuDicts, all_fvDicts
    return fuDict, fvDict


def vizInterp(nLons, nLats, ufunct, vfunct, files):
    '''Creates a visualization of the interpolation function
    nLons: list of mesh node logitudes
    nLats: list of mesh node latitudes
    uFunct: interpolation function of interpolated current u component magnitudes at each mesh node
    vFunct: interpolation function of interpolated current v component magnitudes at each mesh node
    files: list of file input names and directories for data retreival. 
    '''
    
    iModule=np.sqrt(ufunct**2 + vfunct**2)*1.94384
    
    m = Basemap(llcrnrlon=-100.,llcrnrlat=0.,urcrnrlon=30.,urcrnrlat=70, projection='merc', resolution ='l')
    x, y = m(np.array(nLons),np.array(nLats))
    
    print x[:10]
    print y[:10]
    print iModule.squeeze()[:10]
    #print len(iModule.squeeze())
    #print len(x)
    sys.exit()
    m.drawcoastlines()
    #m.fillcontinents(color='coral',lake_color='white')
    m.drawparallels(np.arange(0.,81.,20.),labels=[True,False,False,False])
    m.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,True])
    #m.drawmapboundary(fill_color='aqua')    
    m.pcolormesh(x, y, iModule.squeeze(), cmap='jet', vmin=0, vmax=2)
    #m.quiver(x[::n, ::n], y[::n, ::n], ufunct[::n, ::n], vfunct[::n, ::n], scale=10, pivot='mid')
    m.colorbar(ticks=np.arange(0,2.1,.2)) 
    
    plt.savefig(files[8], bbox_inches='tight')
    plt.show()    