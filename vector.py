import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import datetime
import math
import copy
import sys

from . import response
from . import spectrum
from . import jsi
from . import realtime_jsi

#/ Parse function sets /#
def parse(input_file,noheader=False):
    v = vectors.input(input_file,noheader)
    return v

#######################################
##          Vector    class          ##
#######################################
class vector:

    def __init__(self,header,tim,wave):
        self.header = header
        self.tim = tim
        self.wave = wave
        self.dt = tim[1] - tim[0]

    def copy(self):
        v2 = copy.deepcopy(self)
        return v2

    #----------------------------------------------#
    #  Wave Analysis
    #----------------------------------------------#
    def integration(wave,dt,low,high):
        w = np.fft.fft(wave)
        freq = np.fft.fftfreq(len(wave),d=dt)
        df = freq[1] - freq[0]

        nt = 10
        low0  = max(0,low - nt*df)
        high0 = high + nt*df

        w2 = np.ones_like(w)
        for (i,f) in enumerate(freq):
            coef = (0.0+1.0j)*2.0*math.pi*f

            if abs(f) < low0:
                w2[i] = 0.0 + 0.0j
            elif abs(f) < low:
                w2[i] = w[i] * (abs(f)-low0)/(low-low0) / coef
            elif abs(f) <= high:
                w2[i] = w[i] / coef
            elif abs(f) <= high0:
                w2[i] = w[i] * (abs(f)-high0)/(high-high0) / coef
            else:
                w2[i] = 0.0 + 0.0j

        wave_int = np.real(np.fft.ifft(w2))
        return wave_int

    def differentiation(wave,dt,low,high):
        w = np.fft.fft(wave)
        freq = np.fft.fftfreq(len(wave),d=dt)
        df = freq[1] - freq[0]

        nt = 10
        low0  = max(0,low - nt*df)
        high0 = high + nt*df

        w2 = np.ones_like(w)
        for (i,f) in enumerate(freq):
            coef = (0.0+1.0j)*2.0*math.pi*f

            if abs(f) < low0:
                w2[i] = 0.0 + 0.0j
            elif abs(f) < low:
                w2[i] = w[i] * (abs(f)-low0)/(low-low0) * coef
            elif abs(f) <= high:
                w2[i] = w[i] * coef
            elif abs(f) <= high0:
                w2[i] = w[i] * (abs(f)-high0)/(high-high0) * coef
            else:
                w2[i] = 0.0 + 0.0j

        wave_int = np.real(np.fft.ifft(w2))
        return wave_int

    def bandpass(wave,dt,low,high):
        w = np.fft.fft(wave)
        freq = np.fft.fftfreq(len(wave),d=dt)
        df = freq[1] - freq[0]

        nt = 10
        low0  = max(0,low - nt*df)
        high0 = high + nt*df

        w2 = np.ones_like(w)
        for (i,f) in enumerate(freq):
            if abs(f) < low0:
                w2[i] = 0.0 + 0.0j
            elif abs(f) < low:
                w2[i] = w[i] * (abs(f)-low0)/(low-low0)
            elif abs(f) <= high:
                w2[i] = w[i]
            elif abs(f) <= high0:
                w2[i] = w[i] * (abs(f)-high0)/(high-high0)
            else:
                w2[i] = 0.0 + 0.0j

        wave_int = np.real(np.fft.ifft(w2))
        return wave_int

    def bandpass_butter(wave,dt,low,high,order=5):
        fs = 1.0/dt
        nyq = 0.5*fs
        low_cut = low/nyq
        high_cut = high/nyq
        b,a = signal.butter(order,[low_cut,high_cut],btype='band')
        wave2 = signal.lfilter(b,a,wave)
        return wave2

    def envelope_vector(wave):
        wave_ana = signal.hilbert(wave)
        wave_env = np.abs(wave_ana)
        return wave_env

    def envelope(self):
        v2 = copy.deepcopy(self)
        wave_ana = signal.hilbert(self.wave)
        v2.wave = np.abs(wave_ana)
        return v2

    def power_spectrum(wave,fs):
        freq, P = signal.periodogram(wave,fs)
        return freq, P

    def fourier_spectrum(wave,fs):
        def set_parameters(wave,dt):
            n = len(wave) / 2
            ntw = 2
            while ntw < n:
                ntw = ntw*2
            freq = abs(np.fft.fftfreq(ntw,d=dt))
            return freq, ntw

        dt = 1.0/float(fs)
        freq, ntw = set_parameters(wave,dt)
        w = np.fft.fft(wave[0:ntw])

        return freq, w

    def peak_ground(wave,print_result=True):
        PG = math.sqrt(np.max(wave**2))

        if print_result:
            print("PGX: ",PG)

        return PG

    def normalize(wave,scale=1.0):
        PG = math.sqrt(np.max(wave**2))
        return wave/PG*scale

    def roll(self,shift):
        nshift = int(shift/self.dt)
        v2 = copy.deepcopy(self)
        v2.wave = np.roll(self.wave,nshift)
        return v2

    def amplify(self,amp):
        v2 = copy.deepcopy(self)
        v2.wave = self.wave*amp
        return v2

    def sum(self,v1):
        v2 = copy.deepcopy(self)
        v2.wave += v1.wave
        return v2

    def extract(self,start,end):
        ns = int(start/self.dt)
        ne = int(end/self.dt)
        ntim = ne-ns
        v2 = copy.deepcopy(self)
        v2.wave = self.wave[ns:ne]
        v2.tim = self.tim[ns:ne]
        v2.header['ntim'] = ntim
        return v2

    #----------------------------------------------#
    #  Residual functions
    #----------------------------------------------#
    def nonlinear_residual(self,v):
        v0 = np.sum(self.wave**2)
        v1 = np.sum(v.wave**2)
        v01 = np.sum((self.wave-v.wave)**2)
        return v01*v01/v0/v1

    def linear_residual(self,v):
        v0 = np.sum(self.wave**2)
        v01 = np.sum((self.wave-v.wave)**2)
        return v01/v0

    #----------------------------------------------#
    #  Plot functions
    #----------------------------------------------#
    def plot(self):

        start = self.tim[0]
        end = self.tim[-1]

        ns = int(start/self.dt)
        ne = int(end/self.dt)

        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax1 = fig.add_subplot(111)

        if 'code' in self.header and 'record_time' in self.header:
            ax1.set_title(self.header['code']+" "+self.header['record_time'])

        ax1.plot(self.tim,self.wave,color='k',lw=1)

        amax = np.amax(np.abs(self.wave[ns:ne]))

        ax1.set_xlim(start,end)
        ax1.set_ylim(-1.2*amax,1.2*amax)
        ax1.set_xlabel("time [s]")

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
        plt.show()

    def plot_with(self,v):
        start = self.tim[0]
        end = min(self.tim[-1],v.tim[-1])

        ns = 0
        ne = int((end-start)/self.dt)

        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        if 'code' in self.header and 'record_time' in self.header:
            ax1.set_title(self.header['code']+" "+self.header['record_time'])

        ax1.plot(self.tim,self.wave,color='k',lw=1)
        ax2.plot(v.tim,v.wave,color='r',lw=1)

        amax = np.amax(np.abs(self.wave[ns:ne]))

        ax1.set_xlim(start,end)
        ax2.set_xlim(start,end)

        ax1.set_ylim(-1.2*amax,1.2*amax)
        ax2.set_ylim(-1.2*amax,1.2*amax)

        ax2.set_xlabel("time [s]")

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
        plt.show()

#######################################
##          Vectors   class          ##
#######################################
class vectors(vector):

    def __init__(self,header,tim,ew,ns,ud):
        self.header = header
        self.tim = tim
        self.ew = ew
        self.ns = ns
        self.ud = ud
        self.dt = tim[1] - tim[0]

    def copy(self):
        v2 = copy.deepcopy(self)
        return v2

    #----------------------------------------------#
    #  Construction functions
    #----------------------------------------------#
    def append(self,v):
        v2 = self
        v2.header = self.header
        v2.ew = np.hstack((self.ew,v.ew))
        v2.ns = np.hstack((self.ns,v.ns))
        v2.ud = np.hstack((self.ud,v.ud))

        dt = self.tim[1] - self.tim[0]
        shift = np.full_like(v.tim,self.tim[-1]+dt)
        tim = v.tim + shift

        v2.tim = np.hstack((self.tim,tim))
        v2.header['ntim'] = len(v2.tim)

        return v2

    def to_vector(self,comp):
        if comp == "ew":
            v = vector(self.header,self.tim,self.ew)
        elif comp == "ns":
            v = vector(self.header,self.tim,self.ns)
        elif comp == "ud":
            v = vector(self.header,self.tim,self.ud)
        return v

    #----------------------------------------------#
    #  Correction functions
    #----------------------------------------------#
    def trend_removal(self):
        self.ew = self.ew - np.average(self.ew)
        self.ns = self.ns - np.average(self.ns)
        self.ud = self.ud - np.average(self.ud)

    def rotation(self,rot,deg='deg'):

        if deg == 'deg':
            rad = rot/180.0 * math.pi
        else:
            rad = rot

        p = self.ew * math.sin(rad) + self.ns * math.cos(rad)
        n = self.ew * math.cos(rad) - self.ns * math.sin(rad)

        v2 = copy.deepcopy(self)
        v2.ew = n
        v2.ns = p

        return v2

    def zero_padding(self,ntim):
        self.trend_removal()
        ntim0 = len(self.tim)

        if ntim0 < ntim:
            zero = np.zeros(ntim-ntim0,dtype=float)

            v2 = self
            v2.header = self.header
            v2.ew = np.hstack((self.ew,zero))
            v2.ns = np.hstack((self.ns,zero))
            v2.ud = np.hstack((self.ud,zero))

            dt = self.tim[1] - self.tim[0]
            v2.tim = np.linspace(0,ntim*dt,ntim,endpoint=False)
            v2.header['ntim'] = len(v2.tim)

        else:
            v2 = copy.deepcopy(self)

        return v2

    def time_shift(self,shift):
        nshift = int(shift/self.dt)
        v2 = copy.deepcopy(self)
        v2.ew = np.roll(self.ew,-nshift)
        v2.ns = np.roll(self.ns,-nshift)
        v2.ud = np.roll(self.ud,-nshift)
        return v2

    def normalize(self,scale=1.0):
        v2 = copy.deepcopy(self)
        v2.ew = vector.normalize(self.ew,scale)
        v2.ns = vector.normalize(self.ns,scale)
        v2.ud = vector.normalize(self.ud,scale)
        return v2

    def level_shift(self,level):
        v2 = copy.deepcopy(self)
        v2.ew = self.ew + level
        v2.ns = self.ns + level
        v2.ud = self.ud + level
        return v2

    def trim(self,ntim):
        v2 = copy.deepcopy(self)
        v2.ew = self.ew[0:ntim]
        v2.ns = self.ns[0:ntim]
        v2.ud = self.ud[0:ntim]
        v2.tim = self.tim[0:ntim]
        v2.header['ntim'] = ntim

        return v2

    def adjust(self,v):
        ntim = v.header['ntim']

        if self.header['ntim'] > ntim:
            v2 = self.trim(ntim)
        else:
            v2 = self.zero_padding(ntim)

        return v2

    def every(self,n):
        v2 = copy.deepcopy(self)
        v2.ew = self.ew[0::n]
        v2.ns = self.ns[0::n]
        v2.ud = self.ud[0::n]
        v2.tim = self.tim[0::n]
        v2.dt = v2.tim[1] - v2.tim[0]

        return v2

    def resampling(self,fs):
        f = int(1.0/self.dt)
        if f > fs:
            n = int(f/fs)
            v2 = self.every(n)
        else:
            v2 = self

        return v2

    def amplify(self,amp):
        v2 = copy.deepcopy(self)
        v2.ew = self.ew*amp
        v2.ns = self.ns*amp
        v2.ud = self.ud*amp
        return v2

    def extract(self,start,end=60,to_end=True):
        ns = int(start/self.dt)
        if to_end:
            ne = self.header['ntim']
        else:
            ne = int(end/self.dt)
        ntim = ne-ns

        v2 = copy.deepcopy(self)
        v2.ew = self.ew[ns:ne]
        v2.ns = self.ns[ns:ne]
        v2.ud = self.ud[ns:ne]
        v2.tim = self.tim[ns:ne]
        v2.header['ntim'] = ntim
        return v2

    #----------------------------------------------#
    #  Wave Analysis
    #----------------------------------------------#
    def integration(self,low=0.2,high=50):
        v2 = copy.deepcopy(self)
        v2.ew = vector.integration(self.ew,self.dt,low,high)
        v2.ns = vector.integration(self.ns,self.dt,low,high)
        v2.ud = vector.integration(self.ud,self.dt,low,high)

        return v2

    def differentiation(self,low=0.05,high=50):
        v2 = copy.deepcopy(self)
        v2.ew = vector.differentiation(self.ew,self.dt,low,high)
        v2.ns = vector.differentiation(self.ns,self.dt,low,high)
        v2.ud = vector.differentiation(self.ud,self.dt,low,high)

        return v2

    def bandpass(self,low=0.05,high=10):
        v2 = copy.deepcopy(self)
        v2.ew = vector.bandpass(self.ew,self.dt,low,high)
        v2.ns = vector.bandpass(self.ns,self.dt,low,high)
        v2.ud = vector.bandpass(self.ud,self.dt,low,high)

        return v2

    def bandpass_butter(self,low=0.05,high=10,order=5):
        v2 = copy.deepcopy(self)
        v2.ew = vector.bandpass_butter(self.ew,self.dt,low,high,order)
        v2.ns = vector.bandpass_butter(self.ns,self.dt,low,high,order)
        v2.ud = vector.bandpass_butter(self.ud,self.dt,low,high,order)
        return v2

    def envelope(self):
        v2 = copy.deepcopy(self)
        v2.ew = vector.envelope_vector(self.ew)
        v2.ns = vector.envelope_vector(self.ns)
        v2.ud = vector.envelope_vector(self.ud)
        return v2

    def power_spectrum(self):
        fs = 1.0/self.dt
        freq, Pew = vector.power_spectrum(self.ew,fs)
        freq, Pns = vector.power_spectrum(self.ns,fs)
        freq, Pud = vector.power_spectrum(self.ud,fs)

        self.ps = spectrum.spectrum(self.header,freq,Pew,Pns,Pud)
        return self.ps

    def fourier_spectrum(self):
        fs = 1.0/self.dt
        freq, Few = vector.fourier_spectrum(self.ew,fs)
        freq, Fns = vector.fourier_spectrum(self.ns,fs)
        freq, Fud = vector.fourier_spectrum(self.ud,fs)

        self.fs = spectrum.spectrum(self.header,freq,Few,Fns,Fud)
        return self.fs

    def fourier_amp_spectrum(self):
        fs = 1.0/self.dt
        freq, Few = vector.fourier_spectrum(self.ew,fs)
        freq, Fns = vector.fourier_spectrum(self.ns,fs)
        freq, Fud = vector.fourier_spectrum(self.ud,fs)

        self.fs = spectrum.spectrum(self.header,freq,np.abs(Few),np.abs(Fns),np.abs(Fud))
        return self.fs

    def spectrum_ratio(self,ps0):
        fs = 1.0/self.dt
        freq, Pew = vector.power_spectrum(self.ew,fs)
        freq, Pns = vector.power_spectrum(self.ns,fs)
        freq, Pud = vector.power_spectrum(self.ud,fs)

        nf = min(len(freq),len(ps0.freq))

        sr_ew = np.sqrt(Pew[0:nf] / ps0.ew[0:nf])
        sr_ns = np.sqrt(Pns[0:nf] / ps0.ns[0:nf])
        sr_ud = np.sqrt(Pud[0:nf] / ps0.ud[0:nf])

        sr = spectrum.spectrum(self.header,freq,sr_ew,sr_ns,sr_ud)
        return sr

    def response_spectrum(self):
        self.period = np.logspace(-1,1,100)

        self.Sa_ew, self.Sv_ew, self.Sd_ew = response.response_spectrum(self.ew,self.period,self.dt)
        self.Sa_ns, self.Sv_ns, self.Sd_ns = response.response_spectrum(self.ns,self.period,self.dt)
        self.Sa_ud, self.Sv_ud, self.Sd_ud = response.response_spectrum(self.ud,self.period,self.dt)

    def response_spectrum_max(self):
        self.period = np.logspace(-1,1,100)

        self.Sa, self.Sv, self.Sd, self.Sa_rot, self.Sv_rot, self.Sd_rot = response.response_spectrum_max(self.ew,self.ns,self.period,self.dt)

    def response_spectrum_peak_period(self):
        self.period = np.logspace(-1,1,100)
        pSv,peak_pSv,peak_period = response.response_spectrum_pseudo(self.ew,self.ns,self.period,self.dt)

        return peak_pSv, peak_period

    def calc_SI(self):
        SI = response.calc_SI(self.ew,self.ns,self.dt)
        return SI

    def peak_ground_3d(self,print_result=True):
        PG = math.sqrt(np.max(self.ew**2 + self.ns**2 + self.ud**2))

        if print_result:
            print("PGX (3d): ",PG)

        return PG

    def peak_ground_2d(self,print_result=True):
        PG = math.sqrt(np.max(self.ew**2 + self.ns**2))

        if print_result:
            print("PGX (2d): ",PG)

        return PG

    def peak_ground_all(self,print_result=True):
        PG_ew = math.sqrt(np.max(self.ew**2))
        PG_ns = math.sqrt(np.max(self.ns**2))
        PG_ud = math.sqrt(np.max(self.ud**2))

        if print_result:
            print("PGX (ew): ",PG_ew)
            print("PGX (ns): ",PG_ns)
            print("PGX (ud): ",PG_ud)

        return PG_ew, PG_ns, PG_ud

    def jma_seismic_intensity(self,print_result=True):
        si = jsi.jsi(self.ew,self.ns,self.ud,self.dt)

        if print_result:
            print("JMA Seismic Intensity: ", si)

        return si

    def realtime_seismic_intensity(self):
        rsi = realtime_jsi.realtime_jsi(self.ew,self.ns,self.ud,self.dt)
        return rsi

    def peak_period_Sv(self):
        period = np.logspace(-1,1,100)
        rot_list = np.linspace(0,180,num=8)
        Sv_max = 0.0

        for rot in rot_list:
            w = self.rotation(rot)
            Sa, Sv, Sd = response.response_spectrum(w.ns,period,self.dt)

            if Sv_max < np.max(Sv):
                Sv_max = np.max(Sv)
                peak_period = period[np.argmax(Sv)]

        return peak_period, Sv_max

    def hv_spectrum(self,nt=4096,start=60,ncut=10,window=0.2):

        fsample = 1.0/self.dt
        nseg = int(((len(self.tim) - 30) - start)/(nt/2))

        index = start
        for i in range(0,nseg):
            ew = self.ew[index:index+nt]
            ns = self.ns[index:index+nt]
            ud = self.ud[index:index+nt]
            rms = math.sqrt(np.var(ew) + np.var(ns) + np.var(ud))

            try:
                rms_list += [rms]
                index_list += [index]
            except:
                rms_list = [rms]
                index_list = [index]

            index += int(nt/2)

        r = np.array(rms_list)
        sorted_index = np.argsort(r)

        for i in sorted_index[0:10]:
            index = index_list[i]

            ew = self.ew[index:index+nt]
            ns = self.ns[index:index+nt]
            ud = self.ud[index:index+nt]

            ew -= np.average(ew)
            ns -= np.average(ns)
            ud -= np.average(ud)

            freq, Few = vector.fourier_spectrum(ew,fsample)
            freq, Fns = vector.fourier_spectrum(ns,fsample)
            freq, Fud = vector.fourier_spectrum(ud,fsample)

            fs = spectrum.spectrum(self.header,freq,Few,Fns,Fud)
            fs_smooth = fs.smoothing(window)

            try:
                fs_list.append(fs_smooth)
            except:
                fs_list = spectrum.spectrum_list(fs_smooth)

        fs_list.average(window)

        hv_list = []
        for s in fs_list.sr:
            hv_tmp = np.sqrt(np.abs(s.ew)**2 + np.abs(s.ns)**2) / np.abs(s.ud)
            hv_list += [hv_tmp]

        hv = np.sqrt(fs_list.ave.ew**2 + fs_list.ave.ns**2) / fs_list.ave.ud

        return fs_list.freq, hv, hv_list

    #----------------------------------------------#
    #  I/O functions
    #----------------------------------------------#
    def input(input_file,noheader=False):
        with open(input_file,'r') as file:
            lines = list(file)

        if noheader:
            header = ""
            tim,ew,ns,ud = np.loadtxt(lines[0:], usecols=(0,1,2,3), unpack=True)
        else:
            h = lines[0].strip('#').rstrip()
            header = eval(h)
            tim,ew,ns,ud = np.loadtxt(lines[1:], usecols=(0,1,2,3), unpack=True)

        v = vectors(header,tim,ew,ns,ud)
        return v


    def output(self,output_file,fmt="%15.7f"):
        header = str(self.header)
        output = np.c_[self.tim,self.ew,self.ns,self.ud]
        np.savetxt(output_file,output,fmt=fmt,header=header,comments="#")


    def output_rs(self,output_file,fmt="%15.7f"):
        header = str(self.header)

        output = np.c_[self.period,self.Sa_ew,self.Sa_ns,self.Sa_ud]
        np.savetxt(output_file+".sa",output,fmt=fmt,header=header,comments="#")

        output = np.c_[self.period,self.Sv_ew,self.Sv_ns,self.Sv_ud]
        np.savetxt(output_file+".sv",output,fmt=fmt,header=header,comments="#")

        output = np.c_[self.period,self.Sd_ew,self.Sd_ns,self.Sd_ud]
        np.savetxt(output_file+".sd",output,fmt=fmt,header=header,comments="#")

    def output_rs_max(self,output_file,fmt="%15.7f"):
        header = str(self.header)

        output = np.c_[self.period,self.Sa,self.Sa_rot]
        np.savetxt(output_file+".samax",output,fmt=fmt,header=header,comments="#")

        output = np.c_[self.period,self.Sv,self.Sv_rot]
        np.savetxt(output_file+".svmax",output,fmt=fmt,header=header,comments="#")

        output = np.c_[self.period,self.Sd,self.Sd_rot]
        np.savetxt(output_file+".sdmax",output,fmt=fmt,header=header,comments="#")


    def output_ps(self,output_file,fmt="%15.7f"):
        self.ps.output(output_file,fmt=fmt)

    #----------------------------------------------#
    #  Plot functions
    #----------------------------------------------#
    def plot_all(self,start=0,end=60,to_end=False):

        if to_end:
            start = self.tim[0]
            end = self.tim[-1]
            amax = np.amax(np.r_[np.abs(self.ew),np.abs(self.ns),np.abs(self.ud)])
        else:
            ns = int(start/self.dt)
            ne = int(end/self.dt)
            amax = np.amax(np.r_[np.abs(self.ew[ns:ne]),np.abs(self.ns[ns:ne]),np.abs(self.ud[ns:ne])])

        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        if 'code' in self.header and 'record_time' in self.header:
            ax1.set_title(self.header['code']+" "+self.header['record_time'])


        ax1.plot(self.tim,self.ew,color='k',lw=1)
        ax2.plot(self.tim,self.ns,color='k',lw=1)
        ax3.plot(self.tim,self.ud,color='k',lw=1)



        ax1.set_xlim(start,end)
        ax2.set_xlim(start,end)
        ax3.set_xlim(start,end)

        ax1.set_ylim(-1.2*amax,1.2*amax)
        ax2.set_ylim(-1.2*amax,1.2*amax)
        ax3.set_ylim(-1.2*amax,1.2*amax)

        ax1.set_xticklabels([])
        ax2.set_xticklabels([])

        ax3.set_xlabel("time [s]")

        ax1.text(0.02,0.85,"EW",transform=ax1.transAxes)
        ax2.text(0.02,0.85,"NS",transform=ax2.transAxes)
        ax3.text(0.02,0.85,"UD",transform=ax3.transAxes)

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
        plt.show()


    def plot_rs_all(self):
        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        plt.subplot(331)
        plt.title("Sa")
        plt.plot(self.period,self.Sa_ew,color='maroon',lw=1.5)
        plt.ylabel("EW")

        plt.subplot(332)
        plt.title("Sv")
        plt.plot(self.period,self.Sv_ew,color='darkorange',lw=1.5)

        plt.subplot(333)
        plt.title("Sd")
        plt.plot(self.period,self.Sd_ew,color='darkgreen',lw=1.5)


        plt.subplot(334)
        plt.plot(self.period,self.Sa_ns,color='maroon',lw=1.5)
        plt.ylabel("NS")

        plt.subplot(335)
        plt.plot(self.period,self.Sv_ns,color='darkorange',lw=1.5)

        plt.subplot(336)
        plt.plot(self.period,self.Sd_ns,color='darkgreen',lw=1.5)


        plt.subplot(337)
        plt.plot(self.period,self.Sa_ud,color='maroon',lw=1.5)
        plt.ylabel("UD")
        plt.xlabel("period [s]")

        plt.subplot(338)
        plt.plot(self.period,self.Sv_ud,color='darkorange',lw=1.5)
        plt.xlabel("period [s]")

        plt.subplot(339)
        plt.plot(self.period,self.Sd_ud,color='darkgreen',lw=1.5)
        plt.xlabel("period [s]")

        axs = plt.gcf().get_axes()
        for ax in axs:
            plt.axes(ax)
            plt.xscale("log")
            plt.yscale("log")
            plt.xlim(0.1,10)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        plt.show()

    def plot_rs_max(self):
        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        plt.subplot(231)
        plt.xscale("log")
        plt.yscale("log")
        plt.title("Sa")
        plt.plot(self.period,self.Sa,color='maroon',lw=1.5)
        plt.ylabel("Sa")

        plt.subplot(232)
        plt.xscale("log")
        plt.yscale("log")
        plt.title("Sv")
        plt.plot(self.period,self.Sv,color='darkorange',lw=1.5)
        plt.ylabel("Sv")

        plt.subplot(233)
        plt.xscale("log")
        plt.yscale("log")
        plt.title("Sd")
        plt.plot(self.period,self.Sd,color='darkgreen',lw=1.5)
        plt.ylabel("Sd")


        plt.subplot(234)
        plt.xscale("log")
        plt.ylim(-90,90)
        plt.plot(self.period,self.Sa_rot,color='maroon',lw=1.5)
        plt.ylabel("azimth")
        plt.xlabel("period [s]")

        plt.subplot(235)
        plt.xscale("log")
        plt.ylim(-90,90)
        plt.plot(self.period,self.Sv_rot,color='darkorange',lw=1.5)
        plt.xlabel("period [s]")

        plt.subplot(236)
        plt.xscale("log")
        plt.ylim(-90,90)
        plt.plot(self.period,self.Sd_rot,color='darkgreen',lw=1.5)
        plt.xlabel("period [s]")

        plt.show()

    def plot_ps_all(self):
        self.ps.plot_all()
