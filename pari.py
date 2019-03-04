import numpy as np
from . import vector
import sys

#/ Parse function set /#
def parse(file_name,code,record_time,lat,lon):

    p = pari()
    p.parse(file_name)
    p.set_header(code,record_time,lat,lon)
    v = p.to_vectors()

    return v


#######################################
##          PARI       class         ##
#######################################
class pari:

    def to_vectors(self):
        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v

    def parse(self,file_name):
        try:
            with open(file_name,'r') as file:
                lines = list(file)

            header = ""
            ew,ns,ud = np.loadtxt(lines[1:], delimiter=",", usecols=(0,1,2), unpack=True)

            self.ew = np.array(ew)
            self.ns = np.array(ns)
            self.ud = np.array(ud)

        except IOError as e:
            print("File IO Error: ",e.strerror)


    def set_header(self,code,record_time,lat,lon):

        ntim = len(self.ew)

        self.header = {'code': code,
                    'record_time': record_time,
                    'lat': lat,'lon': lon,
                    'ntim':ntim}

        self.tim = np.linspace(0.0,0.01*ntim,ntim,endpoint=False)
