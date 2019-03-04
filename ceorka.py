import numpy as np
from . import vector
import sys

#/ Parse function set /#
def parse(file_name):

    c = ceorka()
    c.parse(file_name)
    v = c.to_vectors()

    return v

#######################################
##       CEORKA       class          ##
#######################################
class ceorka:

    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):

        vtmp = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        vtmp.trend_removal()
        v = vtmp.rotation(-self.rot)

        return v

    ### Parse CEORKA data ###
    def parse(self,file_name):
        try:
            ceo_file = open(file_name)
            try:
                datalines = ceo_file.readlines()
                self.read_ceorka_data(datalines,file_name)

            finally:
                ceo_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)


    def read_ceorka_data(self,datalines,file_name):

        code = file_name[-12:-9]
        record_time = datalines[0][25:44]

        index_lat = datalines[8].find('SITE LOCATION=')
        lat = float(datalines[8][index_lat+15:index_lat+23])
        lon = float(datalines[8][index_lat+25:index_lat+33])

        ew_body,ns_body,ud_body,self.rot = self.read_body(datalines[16:])

        self.ew = np.array(ew_body)
        self.ns = np.array(ns_body)
        self.ud = np.array(ud_body)

        ntim = len(ew_body)

        self.header = {'code': code,
                       'record_time': record_time,
                       'lat': lat,'lon': lon,
                       'ntim':ntim}

        self.tim = np.linspace(0.0,0.01*ntim,ntim,endpoint=False)

    def read_body(self,datalines):

        lists = datalines[0].strip().split()
        ndat = int(lists[0])
        nline_ns = int(ndat / 10) + 1
        sf_ns = float(datalines[0][40:48])
        rot = int(lists[6][1:4])

        ns = []
        for line in datalines[1:nline_ns]:
            s = line.rstrip()
            lists = [s[i: i+8] for i in range(0, len(line)-1, 8)]
            for w in lists:
                ns.append(float(w) * sf_ns)


        lists = datalines[nline_ns].strip().split()
        ndat = int(lists[0])
        nline_ew = int(ndat / 10) + 1
        sf_ew = float(datalines[nline_ns][40:48])

        ew = []
        for line in datalines[nline_ns+1:nline_ns+nline_ew]:
            s = line.rstrip()
            lists = [s[i: i+8] for i in range(0, len(line)-1, 8)]
            for w in lists:
                ew.append(float(w) * sf_ew)

        lists = datalines[nline_ns+nline_ew].strip().split()
        ndat = int(lists[0])
        nline = int(ndat / 10)
        sf_ud = float(datalines[nline_ns+nline_ew][40:48])

        ud = []
        for line in datalines[nline_ns+nline_ew+1:]:
            s = line.rstrip()
            lists = [s[i: i+8] for i in range(0, len(line)-1, 8)]
            for w in lists:
                ud.append(float(w) * sf_ud)

        return ew, ns, ud, rot
