import numpy as np
from . import vector
import sys

#/ Parse function set /#
def parse(file_name):

    g = gns()
    g.parse(file_name)
    v = g.to_vectors()

    return v

#######################################
##       GNS (GeoNet) class          ##
#######################################
class gns:

    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):
        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v

    ### Parse JMA data ###
    def parse(self,file_name):
        try:
            file = open(file_name)
            try:
                datalines = file.readlines()
                self.read_gns_data(datalines,file_name)
            finally:
                file.close()
        except IOError as e:
            print("File IO Error: ",e.strerror)

    # --------------------------------------------------------------------------- #
    #   Parse related methods
    # --------------------------------------------------------------------------- #
    def read_gns_data(self,datalines,file_name):
#        code = file_name.strip().split("_")[-2]

        nline = 0
        code,record_time,lat,lon,dt,rot1,ndat,ndat_acc = self.read_header(datalines,nline)
        ew_body,nline = self.read_body(datalines,nline,ndat,ndat_acc)

        _,record_time,lat,lon,dt,rot2,ndat,ndat_acc = self.read_header(datalines,nline)
        ns_body,nline = self.read_body(datalines,nline,ndat,ndat_acc)

        _,record_time,lat,lon,dt,rot3,ndat,ndat_acc = self.read_header(datalines,nline)
        ud_body,nline = self.read_body(datalines,nline,ndat,ndat_acc)

        ntim = ndat_acc
        self.ew = np.array(ew_body)*0.1
        self.ns = np.array(ns_body)*0.1
        self.ud = np.array(ud_body)*0.1

        self.header = {'code': code,
                       'record_time': record_time,
                       'lat': lat,'lon': lon,
                       'ntim':ntim}

        self.tim = np.linspace(0.0,dt*ntim,ntim,endpoint=False)

    def read_header(self,datalines,nline):
        code = datalines[nline+1].strip().split()[1]

        list0 = datalines[nline+16].strip().split()
        list1 = datalines[nline+17].strip().split()
        list2 = datalines[nline+18].strip().split()
        list3 = datalines[nline+19].strip().split()

        record_time = '{0:04}/{1:02}/{2:02} {3:02}:{4:02}'.\
            format(int(list0[8]),int(list0[9]),int(list1[8]),int(list1[9]),int(list3[8]))

        lat = '-{:7.4f}'.format(int(list2[0])+(int(list2[1])+int(list2[2])/60)/60)
        lon = '{:8.4f}'.format(int(list2[3])+(int(list2[4])+int(list2[5])/60)/60)

        rot = int(list2[7])
        ndat = int(list3[0])
        ndat_acc = int(list3[3])

        listf = datalines[nline+22].strip().split()
        dt = float(listf[5])

        return code,record_time,lat,lon,dt,rot,ndat,ndat_acc

    def read_body(self,datalines,nline,ndat,ndat_acc):

        w = []
        for line in datalines[nline+26:]:
            lists = [float(line[i:i+8]) for i in range(0,len(line)-1,8)]
            w += lists
            if len(w) >= ndat_acc:
                break

        return w, nline+27+int((ndat_acc-1)/10)
