import numpy as np
import datetime
from . import vector
import sys

#/ Parse function set /#
def parse(file_name):

    j = jma()
    j.parse(file_name)
    v = j.to_vectors()

    return v

def csv_parse(file_name):

    j = jma()
    j.csv_parse(file_name)
    v = j.to_vectors()

    return v


#######################################
##          JMA       class          ##
#######################################
class jma:

    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):

        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v


    ### Parse JMA data ###
    def parse(self,file_name):

        try:
            jma_file = open(file_name)

            try:
                datalines = jma_file.readlines()
                self.read_jma_data(datalines,file_name)

            finally:
                jma_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)


    def csv_parse(self,file_name):

        try:
            jma_file = open(file_name, encoding='shift_jis')

            try:
                datalines = jma_file.readlines()
                self.read_jma_csv_data(datalines,file_name)

            finally:
                jma_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)


    # --------------------------------------------------------------------------- #
    #   Parse related methods
    # --------------------------------------------------------------------------- #
    def read_jma_data(self,datalines,file_name):

        code = file_name[-3:]
        data_time = self.read_data_time(datalines)
        record_time = data_time[0:2]+"/"+data_time[2:4]+"/"+data_time[4:6]+ \
            " "+data_time[6:8]+":"+data_time[8:10]+":"+data_time[10:12]

        lat = 0.0
        lon = 0.0

        ew_body,ns_body,ud_body = self.read_body(datalines)

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

        for i in range(0,len(datalines)):
            if datalines[i].startswith("02/"):
                line = datalines[i+1]
                data_time = line[8:20]

                return data_time

    def read_body(self,datalines):

        ew_digit = []
        ns_digit = []
        ud_digit = []

        scale_factor = 0.3906e-3

        ibrock = 0
        for index in range(2,62):
            cindex = "%02d"%index + "/61"

            for i in range(ibrock,len(datalines)):
                if datalines[i].startswith(cindex):
                    ibrock = i+1
                    break

            ew_tmp,ns_tmp,ud_tmp = self.hex_to_data(datalines[ibrock])

            ew_digit += ew_tmp
            ns_digit += ns_tmp
            ud_digit += ud_tmp

        ew = [float(w)*scale_factor for w in ew_digit]
        ns = [float(w)*scale_factor for w in ns_digit]
        ud = [float(w)*scale_factor for w in ud_digit]

        return ew, ns, ud


    def hex_to_data(self,line):

        nsize = int(line[0:8],16)
        data_time = line[8:20]

        data_dict = {}

        i0 = 20
        i1, data1 = self.data_decode(line[i0:])
        data_dict.update(data1)

        i0 += i1
        i2, data2 = self.data_decode(line[i0:])
        data_dict.update(data2)

        i0 += i2
        i3, data3 = self.data_decode(line[i0:])
        data_dict.update(data3)

        ns_data = data_dict['0']
        ew_data = data_dict['1']
        ud_data = data_dict['2']

        return ew_data, ns_data, ud_data


    def data_decode(self,data):

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


        nc = int(data[0:4],16)
        nb = int(data[4:5],16)

        if nb == 0:
            db = 1    # 4bit => 1txt
        elif nb == 1:
            db = 2    # 8bit => 2txt
        elif nb == 2:
            db = 4    # 16bit => 4txt
        elif nb == 3:
            db = 6    # 24bit => 6txt
        elif nb == 4:
            db = 8    # 32bit => 8txt


        data_list = [s16(int(data[8:16],16),db=8)]
        i0 = 16
        for i in range(1,100):
            data_list.append(data_list[-1] + s16(int(data[i0:i0+db],16),db=db))
            i0 += db

        if nb == 0:
            i0 += db

        return i0, {str(nc):data_list}

#------------------------------------------------------------#
    def read_jma_csv_data(self,datalines,file_name):

        code = datalines[0][11:-1].replace(' ','')
#        record_time = file_name[-22:-8]
        record_time = datalines[5][14:].strip()
        parse_time = datetime.datetime.strptime(record_time,"%Y %m %d %H %M %S")
        record_time = parse_time.strftime('%Y/%m/%d %H:%M:%S')

        slat = datalines[1][-8:]
        slon = datalines[2][-9:]

        lat = float(slat)
        lon = float(slon)

        ew_body,ns_body,ud_body = self.read_csv_body(datalines)

        self.ew = np.array(ew_body)
        self.ns = np.array(ns_body)
        self.ud = np.array(ud_body)

        ntim = len(ew_body)

        self.header = {'code': code,
                       'record_time': record_time,
                       'lat': lat,'lon': lon,
                       'ntim':ntim}

        self.tim = np.linspace(0.0,0.01*ntim,ntim,endpoint=False)


    def read_csv_body(self,datalines):

        ew = []
        ns = []
        ud = []

        for line in datalines[7:]:
            lists = line.strip().split(',')
#            lists = line.strip().split()

            ns.append(float(lists[0]))
            ew.append(float(lists[1]))
            ud.append(float(lists[2]))

        return ew,ns,ud
