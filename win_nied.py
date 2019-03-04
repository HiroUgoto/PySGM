# -- config: utf-8 --
import struct
import numpy as np
from wave import vector
from wave import win
import math
import sys

#/ Parse function set /#
def win_nied_parse(file_name):

    wn = win_nied()
    wn.parse(file_name)
    v = wn.to_vectors()

    return v


###############################################
##          Kyoshin WIN       class          ##
###############################################
class win_nied:

    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):

        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v

    ### Parse WIN data ###
    def parse(self,file_name):

        try:
            win_file = open(file_name,'rb')

            try:
                datalines = win_file.read()
                self.read_win_data(datalines,file_name)

            finally:
                win_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)

    # --------------------------------------------------------------------------- #
    #   Parse related methods
    # --------------------------------------------------------------------------- #
    def read_win_data(self,datalines,file_name):

        id, lat, lon, data_time, scale_ew, scale_ns, scale_ud, ew_code, ns_code, ud_code, offset = self.read_info(datalines)
        code = '{0:03d}'.format(id)
        record_time = data_time[0:4]+"/"+data_time[4:6]+"/"+data_time[6:8]+ \
                    " "+data_time[8:10]+":"+data_time[10:12]+":"+data_time[12:14]

        ew_body,ns_body,ud_body = self.read_body(datalines,ew_code,ns_code,ud_code,offset)

        self.ew = np.array(ew_body)*scale_ew
        self.ns = np.array(ns_body)*scale_ns
        self.ud = np.array(ud_body)*scale_ud

        ntim = len(ew_body)
        self.header = {'code': code, 'record_time': record_time,
                    'lat': lat,'lon': lon, 'ntim':ntim}

        self.tim = np.linspace(0.0,0.01*ntim,ntim,endpoint=False)

    def read_info(self,datalines):
        items = struct.unpack_from('2s4s',datalines,offset=10)

        code = int(items[0].hex(),16)
        nsize = int(items[1].hex(),16)

        items = struct.unpack_from('4s4s',datalines,offset=20)
        clat = items[0].hex()
        clon = items[1].hex()
        clat_index = clat.find("e")
        clon_index = clon.find("e")
        if clat_index == -1:
            clat_index = 8
        if clon_index == -1:
            clon_index = 8
        lat = float(clat[0:3]) + float("0."+ clat[3:clat_index])
        lon = float(clon[0:3]) + float("0."+ clon[3:clon_index])

        items = struct.unpack_from('8s',datalines,offset=44)
        data_time = items[0].hex()

        items = struct.unpack_from('2s',datalines,offset=74)
        ns_code = str(items[0].hex())[-1]
        items = struct.unpack_from('2s1s',datalines,offset=76)
        su = int(items[0].hex(),16)
        gain = int(items[1].hex(),16)
        items = struct.unpack_from('4s',datalines,offset=80)
        sd = int(items[0].hex(),16)
        scale_ns = float(su)/float(sd)/float(gain)

        items = struct.unpack_from('2s',datalines,offset=94)
        ew_code = str(items[0].hex())[-1]
        items = struct.unpack_from('2s1s',datalines,offset=96)
        su = int(items[0].hex(),16)
        gain = int(items[1].hex(),16)
        items = struct.unpack_from('4s',datalines,offset=100)
        sd = int(items[0].hex(),16)
        scale_ew = float(su)/float(sd)/float(gain)

        items = struct.unpack_from('2s',datalines,offset=114)
        ud_code = str(items[0].hex())[-1]
        items = struct.unpack_from('2s1s',datalines,offset=116)
        su = int(items[0].hex(),16)
        gain = int(items[1].hex(),16)
        items = struct.unpack_from('4s',datalines,offset=120)
        sd = int(items[0].hex(),16)
        scale_ud = float(su)/float(sd)/float(gain)

        offset = 16 + nsize

        return code, lat, lon, data_time, scale_ew, scale_ns, scale_ud, ew_code, ns_code, ud_code, offset

    def read_body(self,datalines,ew_code,ns_code,ud_code,offset):

        ew_digit = []
        ns_digit = []
        ud_digit = []

        for ibrock in range(0,120):

            ew_tmp,ns_tmp,ud_tmp,offset = self.hex_to_data(datalines,ew_code,ns_code,ud_code,offset)

            ew_digit += ew_tmp
            ns_digit += ns_tmp
            ud_digit += ud_tmp

            if not self.data_decode_check(datalines,offset):
                break

        ew = [float(w) for w in ew_digit]
        ns = [float(w) for w in ns_digit]
        ud = [float(w) for w in ud_digit]

        return ew, ns, ud

    def hex_to_data(self,line,ew_code,ns_code,ud_code,offset):

        items = struct.unpack_from('8s4s4s',line,offset=offset)
        data_time = items[0].hex()
        nsize = int(items[2].hex(),16)

        data_dict = {}

        w = win.win()

        i0 = 18
        i1, data1 = w.data_decode(line,offset+i0)
        data_dict.update(data1)

        i0 += i1 + 2
        i2, data2 = w.data_decode(line,offset+i0)
        data_dict.update(data2)

        i0 += i2 + 2
        i3, data3 = w.data_decode(line,offset+i0)
        data_dict.update(data3)

        ew_data = data_dict[ew_code]
        ns_data = data_dict[ns_code]
        ud_data = data_dict[ud_code]

        offset += i0 + i3

        return ew_data, ns_data, ud_data, offset


    def data_decode_check(self,data,i):

        try:
            item = struct.unpack_from('8s',data,offset=i)
            return True
        except:
            return False
