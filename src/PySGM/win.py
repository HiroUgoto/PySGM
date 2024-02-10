# -- coding: utf-8 --
import struct
import numpy as np
from . import vector
import math
import sys

#/ Parse function set /#
def parse(file_name,gain,bit):

#    gain = 0.2 # [G]
#    bit = 24
    w = win()
    w.parse(file_name,gain,bit)
    v = w.to_vectors()

    return v



#######################################
##          WIN       class          ##
#######################################
class win:

    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):

        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v


    ### Parse WIN data ###
    def parse(self,file_name,gain,bit):

        try:
            win_file = open(file_name,'rb')

            try:
                datalines = win_file.read()
                self.read_win_data(datalines,gain,bit)

            finally:
                win_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)


    # --------------------------------------------------------------------------- #
    #   Parse related methods
    # --------------------------------------------------------------------------- #
    def read_win_data(self,datalines,gain,bit):

        code = ""
        data_time = self.read_data_time(datalines)
        record_time = data_time[0:2]+"/"+data_time[2:4]+"/"+data_time[4:6]+ \
                    " "+data_time[6:8]+":"+data_time[8:10]+":"+data_time[10:12]

        lat = 0.0
        lon = 0.0

        ew_body,ns_body,ud_body = self.read_body(datalines,gain,bit)

        self.ew = np.array(ew_body)
        self.ns = np.array(ns_body)
        self.ud = np.array(ud_body)

        ntim = len(ew_body)

        self.header = {'code': code,
                       'record_time': record_time,
                       'lat': lat,'lon': lon,
                       'ntim':ntim}

        self.tim = np.linspace(0.0,0.01*ntim,ntim,endpoint=False)


    def read_data_time(self,datalines):

        items = struct.unpack_from('4s6s',datalines,offset=0)

        nsize = int(items[0].hex(),16)
        data_time = items[1].hex()

        return data_time

    def read_body(self,datalines,gain,bit):

        ew_digit = []
        ns_digit = []
        ud_digit = []

        scale_factor = gain*980 / float(math.pow(2.0,bit-1))

        offset = 0
        for ibrock in range(0,60):

            ew_tmp,ns_tmp,ud_tmp,offset,sec = self.hex_to_data(datalines,offset)

            ew_digit += ew_tmp
            ns_digit += ns_tmp
            ud_digit += ud_tmp

            if not self.data_decode_check(datalines,offset):
                break

        ew = [float(w)*scale_factor for w in ew_digit]
        ns = [float(w)*scale_factor for w in ns_digit]
        ud = [float(w)*scale_factor for w in ud_digit]

        return ew, ns, ud


    def hex_to_data(self,line,offset):

        items = struct.unpack_from('4s6s',line,offset=offset)

        nsize = int(items[0].hex(),16)
        data_time = int(items[1].hex())
        sec = int(items[1][-1:].hex())

        data_dict = {}

        i0 = 10
        i1, data1 = self.data_decode(line,offset+i0)
        data_dict.update(data1)

        i0 += i1
        i2, data2 = self.data_decode(line,offset+i0)
        data_dict.update(data2)

        i0 += i2
        i3, data3 = self.data_decode(line,offset+i0)
        data_dict.update(data3)

        ns_data = data_dict['0']
        ew_data = data_dict['1']
        ud_data = data_dict['2']

        offset += i0 + i3

        return ew_data, ns_data, ud_data, offset, sec

    def data_decode(self,data,i):

        def s16(val,db=4):
            if db == 1:
                return -(val & 0x8) | (val & 0x7)
            elif db == 2:
                return -(val & 0x80) | (val & 0x7f)
            elif db == 4:
                return -(val & 0x8000) | (val & 0x7fff)
            elif db == 6:
                return -(val & 0x800000) | (val & 0x7fffff)
            elif db == 8:
                return -(val & 0x80000000) | (val & 0x7fffffff)


        items = struct.unpack_from('4s',data,offset=i)
        items_data = items[0].hex()

        nc = items_data[0:4]
        nb = int(items_data[4:5],16)
        nsample = int(items_data[5:8],16)

        if nb == 0:
            db = 1    # 4bit => 1txt
            dsample = int(nsample/2)
        elif nb == 1:
            db = 2    # 8bit => 2txt
            dsample = nsample-1
        elif nb == 2:
            db = 4    # 16bit => 4txt
            dsample = 2*(nsample-1)
        elif nb == 3:
            db = 6    # 24bit => 6txt
            dsample = 3*(nsample-1)
        elif nb == 4:
            db = 8    # 32bit => 8txt
            dsample = 4*(nsample-1)

        items = struct.unpack_from('4s',data,offset=i+4)
        data_list = [s16(int(items[0].hex(),16),db=8)]

        i0 = 0
        items = struct.unpack_from(str(dsample)+'s',data,offset=i+8)
        ihex = items[0].hex()

        for i in range(0,nsample-1):
            data_list.append(data_list[-1] + s16(int(ihex[i0:i0+db],16),db=db))
            i0 += db

        return dsample+8, {str(nc)[-1]:data_list}


    def data_decode_check(self,data,i):

        try:
            item = struct.unpack_from('4s',data,offset=i)
            return True
        except:
            return False
