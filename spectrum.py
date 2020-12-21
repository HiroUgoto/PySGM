import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import math
import copy
import sys

#######################################
##        Spectrum    class          ##
#######################################
class spectrum:

    def __init__(self,header,freq,ew,ns,ud):
        self.header = header
        self.freq = freq
        self.ew = ew
        self.ns = ns
        self.ud = ud
        self.df = freq[1] - freq[0]

    #----------------------------------------------#
    #  Spectrum Analysis
    #----------------------------------------------#
    def smoothing(self,window,abs=True):
        s = copy.deepcopy(self)
        nw = int(window/s.df)

        if nw > 0:
            if nw%2 == 0:
                nw += 1

            nw2 = int((nw-1)/2)
            w = signal.parzen(nw)

            a = np.r_[self.ew[nw2:0:-1],self.ew,self.ew[0],self.ew[-1:-nw2:-1]]

            if abs:
                s.ew = np.convolve(w/w.sum(),np.abs(a),mode='valid')
            else:
                s.ew = np.convolve(w/w.sum(),a,mode='valid')

            a = np.r_[self.ns[nw2:0:-1],self.ns,self.ns[0],self.ns[-1:-nw2:-1]]

            if abs:
                s.ns = np.convolve(w/w.sum(),np.abs(a),mode='valid')
            else:
                s.ns = np.convolve(w/w.sum(),a,mode='valid')

            a = np.r_[self.ud[nw2:0:-1],self.ud,self.ud[0],self.ud[-1:-nw2:-1]]

            if abs:
                s.ud = np.convolve(w/w.sum(),np.abs(a),mode='valid')
            else:
                s.ud = np.convolve(w/w.sum(),a,mode='valid')

        return s

    def spectrum_average(s_list,window):
        n = len(s_list)

        if n > 0:
            s = copy.deepcopy(s_list[0])

            s_ew = np.zeros_like(s_list[0].ew)
            s_ns = np.zeros_like(s_list[0].ns)
            s_ud = np.zeros_like(s_list[0].ud)

            for s_tmp in s_list:
                s.ew += s_tmp.ew
                s.ns += s_tmp.ns
                s.ud += s_tmp.ud

            s.ew /= float(n)
            s.ns /= float(n)
            s.ud /= float(n)

            s2 = s.smoothing(window)

            return s2

    def spectrum_ratio(self,s0):
        sr = copy.deepcopy(self)

        sr.ew = self.ew/s0.ew
        sr.ns = self.ns/s0.ns
        sr.ud = self.ud/s0.ud

        return sr

    #----------------------------------------------#
    #  I/O functions
    #----------------------------------------------#
    def output(self,output_file,fmt="%15.7f"):
        header = str(self.header)
        output = np.c_[self.freq,self.ew,self.ns,self.ud]
        np.savetxt(output_file,output,fmt=fmt,header=header,comments="#")


    #----------------------------------------------#
    #  Plot functions
    #----------------------------------------------#
    def plot_fs_all(self):

        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        if 'code' in self.header and 'record_time' in self.header:
            ax1.set_title(self.header['code']+" "+self.header['record_time'])

        ax1.plot(self.freq,10*np.log10(np.abs(self.ew)),color='k',lw=1)
        ax2.plot(self.freq,10*np.log10(np.abs(self.ns)),color='k',lw=1)
        ax3.plot(self.freq,10*np.log10(np.abs(self.ud)),color='k',lw=1)

        amax = np.amax(np.r_[np.abs(self.ew),np.abs(self.ns),np.abs(self.ud)])

        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax3.set_xscale("log")

        ax1.set_xlim(0.05,20)
        ax2.set_xlim(0.05,20)
        ax3.set_xlim(0.05,20)

#        ax1.set_ylim(-50,0)
#        ax2.set_ylim(-50,0)
#        ax3.set_ylim(-50,0)

        ax3.set_xlabel("frequency [Hz]")

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
        plt.show()


    def plot_ps_all(self):

        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        ax1.set_title(self.header['code']+" "+self.header['record_time'])

        ax1.plot(self.freq,10*np.log10(self.ew),color='k',lw=1)
        ax2.plot(self.freq,10*np.log10(self.ns),color='k',lw=1)
        ax3.plot(self.freq,10*np.log10(self.ud),color='k',lw=1)

        amax = np.amax(np.r_[np.abs(self.ew),np.abs(self.ns),np.abs(self.ud)])

        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax3.set_xscale("log")

        ax1.set_xlim(0.05,50)
        ax2.set_xlim(0.05,50)
        ax3.set_xlim(0.05,50)

#        ax1.set_ylim(-150,0)
#        ax2.set_ylim(-150,0)
#        ax3.set_ylim(-150,0)

#        ax1.set_ylim(-100,50)
#        ax2.set_ylim(-100,50)
#        ax3.set_ylim(-100,50)

        ax3.set_xlabel("frequency [Hz]")

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
        plt.show()


    def plot_sr_all(self):

        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        ax1.set_title(self.header['code']+" "+self.header['record_time'])

        ax1.plot(self.freq,np.abs(self.ew),color='k',lw=1)
        ax2.plot(self.freq,np.abs(self.ns),color='k',lw=1)
        ax3.plot(self.freq,np.abs(self.ud),color='k',lw=1)

        amax = np.amax(np.r_[np.abs(self.ew),np.abs(self.ns),np.abs(self.ud)])

        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax3.set_xscale("log")

        ax1.set_xlim(0.05,20)
        ax2.set_xlim(0.05,20)
        ax3.set_xlim(0.05,20)

        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax3.set_yscale("log")

        ax1.set_ylim(0.05,20)
        ax2.set_ylim(0.05,20)
        ax3.set_ylim(0.05,20)

        ax3.set_xlabel("frequency [Hz]")

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
        plt.show()


#######################################
##        Spectrum list class        ##
#######################################
class spectrum_list:

    def __init__(self,sr):
        self.sr = [sr]
        self.header = sr.header
        self.freq = sr.freq

    def append(self,sr):
        self.sr += [sr]

    def average(self,window):
        n = len(self.sr)

        if n > 0:
            s = copy.deepcopy(self.sr[0])

            s.ew = np.zeros_like(self.sr[0].ew)
            s.ns = np.zeros_like(self.sr[0].ns)
            s.ud = np.zeros_like(self.sr[0].ud)

            for s_tmp in self.sr:
                s.ew += np.abs(s_tmp.ew)
                s.ns += np.abs(s_tmp.ns)
                s.ud += np.abs(s_tmp.ud)

            s.ew /= float(n)
            s.ns /= float(n)
            s.ud /= float(n)

            self.ave = s.smoothing(window)

    def statistics_horizontal(self,nw=12):
        n = len(self.sr)
        stat_freq_bounds = np.logspace(-1,1.2,nw)
        stat_freq_width = np.diff(stat_freq_bounds)
        self.stat_freq = (stat_freq_bounds[0:-2] + stat_freq_bounds[1:-1])/2

        freq_window = []
        for i,fc in enumerate(self.stat_freq):
            freq_window += [np.where((self.freq > stat_freq_bounds[i]) & (self.freq < stat_freq_bounds[i+1]))]

        samples = [[] for i in range(len(self.stat_freq))]
        for s in self.sr:
            for i,fc in enumerate(self.stat_freq):
                try:
                    samples[i] = np.append(samples[i],[s.ew[freq_window[i]],s.ns[freq_window[i]]])
                except:
                    samples[i] = np.append(s.ew[freq_window[i]],s.ns[freq_window[i]])

        self.mean = []
        self.mean_nstd,self.mean_pstd = [],[]
        for sample_data in samples:
            m = np.mean(np.log10(sample_data))
            std = np.std(np.log10(sample_data),ddof=1)

            self.mean += [10**m]
            self.mean_pstd += [10**(m+std)]
            self.mean_nstd += [10**(m-std)]


    def output_average(self,output_file,fmt="%15.7f"):
        header = str(self.header)
        output = np.c_[self.freq,self.ave.ew,self.ave.ns,self.ave.ud]
        np.savetxt(output_file,output,fmt=fmt,header=header,comments="#")


    def output_statistics(self,output_file,fmt="%15.7f"):
        header = str(self.header)
        output = np.c_[self.stat_freq,self.mean,self.mean_nstd,self.mean_pstd]
        np.savetxt(output_file,output,fmt=fmt,header=header,comments="#")


    def plot_all(self):

        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        ax1.set_title(self.header['code'])

        for sr in self.sr:
            ax1.plot(self.freq,np.abs(sr.ew),color='k',lw=1)
            ax2.plot(self.freq,np.abs(sr.ns),color='k',lw=1)
            ax3.plot(self.freq,np.abs(sr.ud),color='k',lw=1)

        ax1.plot(self.freq,np.abs(self.ave.ew),color='r',lw=2)
        ax2.plot(self.freq,np.abs(self.ave.ns),color='r',lw=2)
        ax3.plot(self.freq,np.abs(self.ave.ud),color='r',lw=2)

        amax = np.amax(np.r_[np.abs(self.ave.ew),np.abs(self.ave.ns),np.abs(self.ave.ud)])

        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax3.set_xscale("log")

        ax1.set_xlim(0.1,10)
        ax2.set_xlim(0.1,10)
        ax3.set_xlim(0.1,10)

        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax3.set_yscale("log")

        # ax1.set_ylim(0.01,100)
        # ax2.set_ylim(0.01,100)
        # ax3.set_ylim(0.01,100)

        ax3.set_xlabel("frequency [Hz]")

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
        plt.show()
