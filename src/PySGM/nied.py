import numpy as np
import datetime
from . import vector


#/ Parse function sets /#
def knet_parse(file_basename):

    n = nied()
    n.knet_parse(file_basename)
    v = n.to_vectors()

    return v

def kik_surface_parse(file_basename):

    n = nied()
    n.kik_surface_parse(file_basename)
    v = n.to_vectors()

    return v

def kik_borehole_parse(file_basename):

    n = nied()
    n.kik_borehole_parse(file_basename)
    v = n.to_vectors()

    return v


def parse(file_basename, ext):

    if ext[-1] == '2':
        v = kik_surface_parse(file_basename)
    elif ext[-1] == '1':
        v = kik_borehole_parse(file_basename)
    else:
        v = knet_parse(file_basename)

    return v

#######################################
##          NIED      class          ##
#######################################
class nied:

    # --------------------------------------------------------------------------- #
    #   Convert to vectors class
    # --------------------------------------------------------------------------- #
    def to_vectors(self):

        v = vector.vectors(self.header,self.tim,self.ew,self.ns,self.ud)
        return v


    # --------------------------------------------------------------------------- #
    #   Parse related methods
    # --------------------------------------------------------------------------- #
    ### Parse K-NET data ###
    def knet_parse(self,file_basename):

        try:
            ew_file = open(file_basename + ".EW")
            ns_file = open(file_basename + ".NS")
            ud_file = open(file_basename + ".UD")

            try:
                ew_datalines = ew_file.readlines()
                ns_datalines = ns_file.readlines()
                ud_datalines = ud_file.readlines()

#                print(file_basename)
                nied.read_knet_data(self,ew_datalines,ns_datalines,ud_datalines)

            finally:
                ew_file.close()
                ns_file.close()
                ud_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)

    ### Parse KiK-net surface data ###
    def kik_surface_parse(self,file_basename):

        try:
            ew_file = open(file_basename + ".EW2")
            ns_file = open(file_basename + ".NS2")
            ud_file = open(file_basename + ".UD2")

            try:
                ew_datalines = ew_file.readlines()
                ns_datalines = ns_file.readlines()
                ud_datalines = ud_file.readlines()

                nied.read_knet_data(self,ew_datalines,ns_datalines,ud_datalines)

            finally:
                ew_file.close()
                ns_file.close()
                ud_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)

    ### Parse KiK-net surface data ###
    def kik_borehole_parse(self,file_basename):

        try:
            ew_file = open(file_basename + ".EW1")
            ns_file = open(file_basename + ".NS1")
            ud_file = open(file_basename + ".UD1")

            try:
                ew_datalines = ew_file.readlines()
                ns_datalines = ns_file.readlines()
                ud_datalines = ud_file.readlines()

                nied.read_knet_data(self,ew_datalines,ns_datalines,ud_datalines)

            finally:
                ew_file.close()
                ns_file.close()
                ud_file.close()

        except IOError as e:
            print("File IO Error: ",e.strerror)


    def read_knet_data(self,ew_datalines,ns_datalines,ud_datalines):

        ew_header,ew_body = self.read_header_body(ew_datalines)
        ns_header,ns_body = self.read_header_body(ns_datalines)
        ud_header,ud_body = self.read_header_body(ud_datalines)

        self.ew = np.array(ew_body)
        self.ns = np.array(ns_body)
        self.ud = np.array(ud_body)

        ntim = len(ew_body)

        trigger_time = datetime.datetime.strptime(ew_header['record_time'],"%Y/%m/%d %H:%M:%S")
        record_time = trigger_time - datetime.timedelta(seconds=15)

        self.header = {'code':ew_header['code'],
                       'record_time':record_time.strftime('%Y/%m/%d %H:%M:%S'),
                       'lat':ew_header['lat'],'lon':ew_header['lon'],
                       'ntim':ntim}

        self.tim = np.linspace(0.0,ew_header['duration'],ntim,endpoint=False)

    def read_header_body(self,datalines):

        header = self.read_header(datalines)
        body = self.read_body(datalines,header)
        return header, body


    def read_header(self,datalines):

        header_keys = ['event_time','event_lat','event_lon','event_depth','event_mg',
                       'code','lat','lon','height','record_time','sampling','duration',
                   'azimuth','scale_factor','PGA','correction_time']

        items = []
        for line in datalines[0:len(header_keys)]:
            item = line[18:-1]
            items.append(item)

        header = dict(zip(header_keys,items))

        header['lat'] = float(header['lat'])
        header['lon'] = float(header['lon'])
        header['height'] = float(header['height'])

        header['duration'] = float(header['duration'])
        header['PGA'] = float(header['PGA'])

        su,sd = header['scale_factor'].split("(gal)/")
        header['scale_factor'] = float(su)/float(sd)

        sp = header['sampling'].split("Hz")[0]
        header['sampling'] = float(sp)

        return header


    def read_body(self,datalines,header):

        waves = []
        for line in datalines[17:]:
            wave = line[:-1].split()
            for w in wave:
                waves.append(float(w)*header['scale_factor'])
        return waves
    # ---------------------------------------------------------------------------#
