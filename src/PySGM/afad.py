import numpy as np
import datetime
from . import vector

#/ Parse function set /#
def parse(file_name_E,file_name_N,file_name_U):
    header_num=64
    
    af = afad()
    af.parse(file_name_E,file_name_N,file_name_U,header_num)
    v = af.to_vectors()

    return v

#######################################
##          AFAD      class          ##
#######################################
class afad:
    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):
        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v

    # --------------------------------------------------------------------------- #
    #   Parse methods
    # --------------------------------------------------------------------------- #
    def parse(self,file_name_E,file_name_N,file_name_U,header_num=64):
        try:
            ew_file = open(file_name_E)
            ns_file = open(file_name_N)
            ud_file = open(file_name_U)
            
            try:    
                ew_datalines = ew_file.readlines()
                ns_datalines = ns_file.readlines()
                ud_datalines = ud_file.readlines()
            
                self.parse_data(ew_datalines,ns_datalines,ud_datalines,header_num)

            finally:
                ew_file.close()
                ns_file.close()
                ud_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)


    def parse_data(self,ew_datalines,ns_datalines,ud_datalines,header_num=64):
        ew_header = self.parse_header(ew_datalines,header_num)
        ns_header = self.parse_header(ns_datalines,header_num)
        ud_header = self.parse_header(ud_datalines,header_num)

        self.ntim = min([ew_header['ntim'],ns_header['ntim'],ud_header['ntim']])
        self.header = ew_header.copy()

        self.ew,self.ns,self.ud = [],[],[]
        for i in range(self.roll,self.ntim):
            self.ew.append(float(ew_datalines[header_num+i].strip()))
            self.ns.append(float(ns_datalines[header_num+i].strip()))
            self.ud.append(float(ud_datalines[header_num+i].strip()))

        self.ntim = len(self.ew)
        self.header['ntim'] = self.ntim
        self.tim = np.linspace(0.0,self.ntim*self.dt,self.ntim,endpoint=False)


    def parse_header(self,datalines,header_num=64):
        header_keys = ['STATION_CODE:','STATION_LATITUDE_DEGREE:','STATION_LONGITUDE_DEGREE:','NDATA:','SAMPLING_INTERVAL_S:']
        header_keys_long = ['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS:',]

        header_dict = {}
        for i in range(header_num):
            items = datalines[i].strip().split()
            if items[0] in header_keys:
                header_dict[items[0]] = items[1]
            elif items[0] in header_keys_long:
                if len(items) <= 2:
                    header_dict[items[0]] = items[1]
                    try:
                        record_time = datetime.datetime.strptime(header_dict['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS:'],"%Y%m%d_%H%M%S")
                    except:
                        record_time = datetime.datetime.strptime(header_dict['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS:'],"%d%m%Y_%H%M%S")
                    triger_msec = 0.0
                else:
                    header_dict[items[0]] = items[1] + " " + items[2]
                    record_time = datetime.datetime.strptime(header_dict['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS:'][:19],"%Y/%m/%d %H:%M:%S")
                    record_time += datetime.timedelta(seconds=1)
                    triger_msec = float(header_dict['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS:'][19:])

        self.dt = float(header_dict['SAMPLING_INTERVAL_S:'])
        self.roll = int((1.0-triger_msec)/self.dt)          # number of trim data 

        header = {'code':header_dict['STATION_CODE:'],'record_time':record_time.strftime('%Y/%m/%d %H:%M:%S'),
                'lat':float(header_dict['STATION_LATITUDE_DEGREE:']),
                'lon':float(header_dict['STATION_LONGITUDE_DEGREE:']),
                'ntim':int(header_dict['NDATA:'])}

        return header
