# -- coding: utf-8 --
import struct
import numpy as np
from . import vector
import datetime
import sys

#/ Parse function sets /#
def parse(file_basename_x,file_basename_y,file_basename_z):
    sac_data = sac()
    sac_data.parse(file_basename_x,file_basename_y,file_basename_z)
    v = sac_data.to_vectors()
    return v


#######################################
##          SAC       class          ##
#######################################
class sac:

    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):

        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v

    ### Parse SAC data ###
    def parse(self,file_basename_x,file_basename_y,file_basename_z):

        try:
            x_file = open(file_basename_x,'rb')
            y_file = open(file_basename_y,'rb')
            z_file = open(file_basename_z,'rb')

            try:
                x_datalines = x_file.read()
                y_datalines = y_file.read()
                z_datalines = z_file.read()

                self.read_sac_data(x_datalines,y_datalines,z_datalines)

            finally:
                x_file.close()
                y_file.close()
                z_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)


    # --------------------------------------------------------------------------- #
    #   Parse related methods
    # --------------------------------------------------------------------------- #
    def read_sac_data(self,x_datalines,y_datalines,z_datalines):

        x_header,x_body = self.read_header_body(x_datalines)
        y_header,y_body = self.read_header_body(y_datalines)
        z_header,z_body = self.read_header_body(z_datalines)

        x,y,z,header = self.adjust_3comp(x_header,y_header,z_header,x_body,y_body,z_body)

        self.ns =  np.array(x) * 1.e-7
        self.ew =  np.array(y) * 1.e-7
        self.ud =  np.array(z) * 1.e-7

        ntim = header['ntim']

        dt = header['DELTA']
        duration = ntim*dt

        self.header = {'code':header['KSTNM'],
                       'record_time':header['record_time'],
                       'lat':header['STLA'],'lon':header['STLO'],
                       'ntim':ntim}

        self.tim = np.linspace(0.0,duration,ntim,endpoint=False)


    def read_header_body(self,datalines):

        header = self.read_header(datalines)
        body = self.read_body(datalines,header)
        return header, body


    def read_header(self,datalines):

        keys_float = ['DELTA','DEPMIN','DEPMAX','SCALE','ODELTA','B','E','O','A','FINTERNAL1',
                      'T0','T1','T2','T3','T4','T5','T6','T7','T8','T9',
                      'F','RESP0','RESP1','RESP2','RESP3','RESP4','RESP5','RESP6','RESP7','RESP8',
                      'RESP9','STLA','STLO','STEL','STDP','EVLA','EVLO','EVEL','EVDP','MAG',
                      'USER0','USER1','USER2','USER3','USER4','USER5','USER6','USER7','USER8','USER9',
                      'DIST','AZ','BAZ','GCARC','FINTERNAL2','FINTERNAL3','DEPMIN','CMPAZ','CMPINC','XMINIMUM',
                      'XMAXIMUM','YMINIMUM','YMAXIMUM','FUNUSED10','FUNUSED11',
                      'FUNUSED12','FUNUSED13','FUNUSED14','FUNUSED15','FUNUSED16']
        keys_int = ['NZYEAR','NZJDAY','NZHOUR','NZMIN','NZSEC','NZMSEC','NVHDR','NORID','NEVID','NPTS',
                    'IINTERNAL','NWFID','NXSIZE','NYSIZE','UNUSED','IFTYPE','IDEP','IZTYPE','IUNUSEDx','IINST',
                    'ISTREG','IEVREG','IEVTYP','IQUAL','ISYNTH','IMAGTYP','IMAGSRC','IUNUSED0','IUNUSED1','IUNUSED2',
                    'IUNUSED3','IUNUSED4','IUNUSED5','IUNUSED6','IUNUSED7']
        keys_bool = ['LEVEN','LPSPOL','LOVROK','LCALDA','LUNUSED']
        keys_string = ['KSTNM','KEVNM','KHOLE','KO','KA',
                       'KT0','KT1','KT2','KT3','KT4','KT5','KT6','KT7','KT8','KT9',
                       'KF','KUSER0','KUSER1','KUSER2','KCMPNM','KNETWK','KDATRD','KINST']

        items_float = struct.unpack_from('70f',datalines,offset=0)
        items_int = struct.unpack_from('35i',datalines,offset=280)

        items_logical = struct.unpack_from('5i',datalines,offset=420)
        items_bool = [bool(l) for l in items_logical]

        items_byte_string = struct.unpack_from('8s16s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s8s',
                                               datalines,offset=440)
        items_string = [s.decode('utf-8') for s in items_byte_string]


        header = dict(zip(keys_float,items_float))
        header.update(dict(zip(keys_int,items_int)))
        header.update(dict(zip(keys_bool,items_bool)))
        header.update(dict(zip(keys_string,items_string)))

        record_time = str(header['NZYEAR'])+"/"+str(header['NZJDAY'])+" "+str(header['NZHOUR'])+":"+str(header['NZMIN']) +":"+str(header['NZSEC'])
        header.update({'record_time':record_time})

        return header


    def read_body(self,datalines,header):

        ntim = header['NPTS']
        unpack_format = str(ntim) + 'f'

        wave = struct.unpack_from(unpack_format,datalines,offset=632)

        return wave

    def adjust_3comp(self,x_header,y_header,z_header,x_body,y_body,z_body):

        x_time = self.set_time(x_header)
        y_time = self.set_time(y_header)
        z_time = self.set_time(z_header)
        max_time = max(x_time,y_time,z_time)

        x_time_delta = max_time - x_time
        y_time_delta = max_time - y_time
        z_time_delta = max_time - z_time

        x_roll = round(x_time_delta.total_seconds()/x_header['DELTA'])
        y_roll = round(y_time_delta.total_seconds()/y_header['DELTA'])
        z_roll = round(z_time_delta.total_seconds()/z_header['DELTA'])

        nx = len(x_body[x_roll:])
        ny = len(y_body[y_roll:])
        nz = len(z_body[z_roll:])
        ntim = min(nx,ny,nz)

        x = x_body[x_roll:x_roll+ntim]
        y = y_body[y_roll:y_roll+ntim]
        z = z_body[z_roll:z_roll+ntim]

        header = x_header
        header['ntim'] = ntim
        header['record_time']= max_time.strftime('%Y/%m/%d %H:%M:%S.%f')

        return x,y,z,header

    def set_time(self,header):
        time_string = str(header['NZYEAR'])+"/"+str(header['NZJDAY'])+" "+str(header['NZHOUR'])+":"+str(header['NZMIN']) +":"+str(header['NZSEC']) +","+str(header['NZMSEC'])
        time = datetime.datetime.strptime(time_string, "%Y/%j %H:%M:%S,%f")
        return time
