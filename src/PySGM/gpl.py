# -- coding: utf-8 --
import struct
import subprocess
import numpy as np
from . import vector
import math
import sys

#/ Parse function set /#
def parse(file_name,gain):

    g = gpl()
    g.parse(file_name,gain)
    v = g.to_vectors()

    return v

#######################################
##          GPL       class          ##
#######################################
class gpl:

    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):

        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v

    ### Parse WIN data ###
    def parse(self,file_name,gain):

        try:
            file = open(file_name,'rb')

            try:
                datalines = file.read()
                self.read_data(datalines,gain)

            finally:
                file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)


    # --------------------------------------------------------------------------- #
    #   Parse related methods
    # --------------------------------------------------------------------------- #
    def read_data(self,datalines,gain):

        self.read_header(datalines)
        ew_body,ns_body,ud_body = self.read_body(datalines,gain)

        lat = 0.0
        lon = 0.0
        code = ""
        record_time = ""

        self.ew = np.array(ew_body)
        self.ns = np.array(ns_body)
        self.ud = np.array(ud_body)

        ntim = len(self.ew)

        self.header = {'code': code,
                       'record_time': record_time,
                       'lat': lat,'lon': lon,
                       'ntim':ntim}

        self.tim = np.linspace(0.0,self.dt*ntim,ntim,endpoint=False)

    def read_header(self,datalines):

        lists = datalines[0:40].decode().splitlines()
        self.header_byte = int(lists[0][-5:-1]) + 12

        lists = datalines[0:self.header_byte].decode().splitlines()

        interval = int(lists[3].split(',')[-1])
        self.dt = 0.001*interval

        volt_range = lists[4].split(',')[1].split()[0:3]
        self.ud_range = int(volt_range[0])
        self.ns_range = int(volt_range[1])
        self.ew_range = int(volt_range[2])

        self.ndat = int(int(lists[18].split(',')[-1])/3)
        self.sf = float(lists[20].split(',')[-1])


    def read_body(self,datalines,gain):

        ew_digit = []
        ns_digit = []
        ud_digit = []

        scale_factor_ew = gain * self.ew_range / 2147483647 / self.sf
        scale_factor_ns = gain * self.ns_range / 2147483647 / self.sf
        scale_factor_ud = gain * self.ud_range / 2147483647 / self.sf

        offset = self.header_byte
        for index in range(0,self.ndat):
            ud_tmp,ns_tmp,ew_tmp = struct.unpack_from('>3i',datalines,offset=offset)

            ew_digit += [ew_tmp]
            ns_digit += [ns_tmp]
            ud_digit += [ud_tmp]

            offset += 12

        ew = [float(w)*scale_factor_ew for w in ew_digit]
        ns = [float(w)*scale_factor_ns for w in ns_digit]
        ud = [float(w)*scale_factor_ud for w in ud_digit]

        return ew, ns, ud
