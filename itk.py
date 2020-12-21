import struct
import datetime
import numpy as np
from . import vector


#/ Parse function sets /#
def parse(input_dir,file_basename,start_time,end_time, \
                    fmt='%Y%m%d%H%M', sc=2.5*980/(2.0**23) ):

    i = itk()
    i.parse(input_dir,file_basename,start_time,end_time,fmt,sc)
    v = i.to_vectors()

    return v

def parse_SD(input_dir,start_time,end_time,fmt='%Y%m%d%H%M',sc=2.5*980/(2.0**23)):

    i = itk()
    i.parse_SD(input_dir,start_time,end_time,fmt,sc)
    v = i.to_vectors()

    return v

#######################################
##          ITK       class          ##
#######################################
class itk:

    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):

        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v


    ### Parse ITK data ###
    def parse(self,input_dir,file_basename,start_time,end_time,fmt,sc):
        st = datetime.datetime.strptime(start_time,fmt)
        et = datetime.datetime.strptime(end_time,fmt)
        m = (et-st).seconds // 60 + 1

        datalines_list = []
        for im in range(0,m):
            ct = st + datetime.timedelta(minutes=im)
            current_time = ct.strftime('%Y%m%d%H%M')
            file_name = input_dir + current_time + "_" + file_basename + ".raw"

            try:
                itk_file = open(file_name,'rb')
                try:
                    datalines = itk_file.read()
                    datalines_list += [datalines]
                finally:
                    itk_file.close()

            except IOError as e:
                print("File IO Error: ",e.strerror)

        self.read_itk_data(st,datalines_list,sc)


    def parse_SD(self,input_dir,start_time,end_time,fmt,sc):
        st = datetime.datetime.strptime(start_time,fmt)
        et = datetime.datetime.strptime(end_time,fmt)
        m = (et-st).seconds // 60 + 1

        datalines_list = []
        for im in range(0,m):
            ct = st + datetime.timedelta(minutes=im)
            current_time = ct.strftime('%M')
            file_name = input_dir + current_time + ".RAW"

            try:
                itk_file = open(file_name,'rb')
                try:
                    datalines = itk_file.read()
                    datalines_list += [datalines]
                finally:
                    itk_file.close()

            except IOError as e:
                print("File IO Error: ",e.strerror)

        self.read_itk_data(st,datalines_list,sc)

    # --------------------------------------------------------------------------- #
    #   Parse related methods
    # --------------------------------------------------------------------------- #
    def read_itk_data(self,st,datalines_list,sc):

        ew_list = []
        ns_list = []
        ud_list = []

        for dl in datalines_list:
            ew_t,ns_t,ud_t = self.read_body(dl)
            ew_list += ew_t
            ns_list += ns_t
            ud_list += ud_t

        self.ew = np.array(ew_list) * sc
        self.ns = np.array(ns_list) * sc
        self.ud = np.array(ud_list) * sc

        ntim = len(self.ew)
        self.tim = np.linspace(0.0,0.01*ntim,ntim,endpoint=False)

        code = struct.unpack_from('6s',datalines_list[0])[0].hex()
        record_time = st.strftime('%Y/%m/%d %H:%M')
        lat = ""
        lon = ""

        self.header = {'code': code,
                       'record_time': record_time,
                       'lat': lat,'lon': lon,
                       'ntim':ntim}



    def read_body(self,datalines):

        offset_ew = 20
        offset_ns = offset_ew + 6000*4 + 4
        offset_ud = offset_ns + 6000*4 + 4

        items_ew = struct.unpack_from('6000i',datalines,offset=offset_ew)
        items_ns = struct.unpack_from('6000i',datalines,offset=offset_ns)
        items_ud = struct.unpack_from('6000i',datalines,offset=offset_ud)

        return items_ew, items_ns, items_ud
