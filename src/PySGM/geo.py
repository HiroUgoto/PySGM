import numpy as np
import pyproj
import sys

def distance(lon0,lat0,lon1,lat1):
    wgs84 = pyproj.Geod(ellps='WGS84')
    azimuth,azimuth_inv,dist = wgs84.inv(lon0,lat0,lon1,lat1)
    return dist*0.001

def distance_azimuth(lon0,lat0,lon1,lat1):
    wgs84 = pyproj.Geod(ellps='WGS84')
    azimuth,azimuth_inv,dist = wgs84.inv(lon0,lat0,lon1,lat1)
    return dist*0.001, azimuth

def lon_lat(lon0,lat0,azimuth,dist):
    dist = dist*1000.0
    wgs84 = pyproj.Geod(ellps='WGS84')
    lon,lat,azimuth_inv = wgs84.fwd(lon0,lat0,azimuth,dist)
    return lon,lat

def set_site(lon,lat):
    s = site(lon,lat)
    return s

######################################
##      Site and Fault class        ##
######################################
class site:
    def __init__(self,lon,lat,depth=0.0):
        self.lon = lon
        self.lat = lat
        self.depth = depth

    def distance(self,s):
        dist = distance(self.lon,self.lat,s.lon,s.lat)
        dist = np.sqrt(dist**2+(self.depth-s.depth)**2)
        return dist

    def distance_azimuth(self,s):
        dist,azimuth = distance_azimuth(self.lon,self.lat,s.lon,s.lat)
        dist = np.sqrt(dist**2+(self.depth-s.depth)**2)
        return dist,azimuth

    def to_xy(self,s):
        dist,azimuth = distance_azimuth(self.lon,self.lat,s.lon,s.lat)
        x = dist*np.cos(azimuth/180*np.pi)
        y = dist*np.sin(azimuth/180*np.pi)
        return x,y

    def lon_lat(self,azimuth,dist,depth=0.0):
        lon,lat = lon_lat(self.lon,self.lat,azimuth,dist)
        s = site(lon,lat,self.depth+depth)
        return s

class fault:
    def __init__(self,s0,angle,dimension,origin=(0.0,0.0)):
        self.s0 = s0
        self.length,self.width = dimension  # (length,width)[km]
        self.strike,self.dip = angle    # (strike,dip)
        self.origin = origin            # (origin_length,origin_width)
        self.subfault = []

    def set_subfault(self,Nl=10,Nw=10):
        dl,dw = 1.0/Nl,1.0/Nw
        l = np.linspace(-self.origin[0]+dl/2,1.0-self.origin[0]-dl/2,Nl)
        w = np.linspace(-self.origin[1]+dw/2,1.0-self.origin[1]-dw/2,Nw)
        grid_length,grid_width = np.meshgrid(l,w)
        self.dl,self.dw = dl,dw

        sd = np.sin(np.deg2rad(self.dip))
        cd = np.cos(np.deg2rad(self.dip))
        ss = np.sin(np.deg2rad(self.strike))
        cs = np.cos(np.deg2rad(self.strike))

        depth = self.width*grid_width*sd
        length_s = self.length*grid_length
        length_d = self.width*grid_width*cd
        dist = np.sqrt(length_s**2+length_d**2)
        rot = np.rad2deg(np.arctan2(length_d,length_s))
        azimuth = self.strike + rot

        self.subfault = []
        for a,d,dep in zip(azimuth.flatten(),dist.flatten(),depth.flatten()):
            self.subfault += [self.s0.lon_lat(a,d,dep)]

    def set_faultframe(self):
        l = np.array([-self.origin[0],1.0-self.origin[0],1.0-self.origin[0],-self.origin[0]])
        w = np.array([-self.origin[1],-self.origin[1],1.0-self.origin[1],1.0-self.origin[1]])

        sd = np.sin(np.deg2rad(self.dip))
        cd = np.cos(np.deg2rad(self.dip))
        ss = np.sin(np.deg2rad(self.strike))
        cs = np.cos(np.deg2rad(self.strike))

        depth = self.width*w*sd
        length_s = self.length*l
        length_d = self.width*w*cd
        dist = np.sqrt(length_s**2+length_d**2)
        rot = np.rad2deg(np.arctan2(length_d,length_s))
        azimuth = self.strike + rot

        self.faultframe = []
        for a,d,dep in zip(azimuth.flatten(),dist.flatten(),depth.flatten()):
            self.faultframe += [self.s0.lon_lat(a,d,dep)]

    def distance(self,s):
        if not self.subfault:
            self.set_subfault()

        dist = self.s0.distance(s)
        for sub in self.subfault:
            dist_sub = sub.distance(s)
            dist = min(dist,dist_sub)

        return dist
