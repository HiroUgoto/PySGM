import numpy as np
from scipy.signal import lfilter
import matplotlib.pyplot as plt


def realtime_jsi(ew,ns,ud,dt):
    r = rsi()
    r.calc_realtime_jsi(ew,ns,ud,dt)
    return r

# ============================== #
class rsi():
    def calc_realtime_jsi(self,ew,ns,ud,dt):
        self.dt = dt
        self.ntim = len(ew)
        self.tim = np.linspace(0,self.ntim*dt,self.ntim)

        fs = 1.0/dt

        ew_filtered = self.jsi_filter(ew,dt)
        ns_filtered = self.jsi_filter(ns,dt)
        ud_filtered = self.jsi_filter(ud,dt)

        wave = np.sqrt(ew_filtered**2 + ns_filtered**2 + ud_filtered**2)
        pm = self.pseudo_max(wave)
        self.I = 2*np.log10(pm) + 0.94

    def pseudo_max(self,x):
        pm_list = np.zeros(self.ntim)
        pm_stock = np.zeros(30)
        for i in range(0,self.ntim):
            pm_stock_31 = sorted(np.append(pm_stock,x[i]),reverse=True)
            pm_stock = pm_stock_31[0:-1]
            pm_list[i] = max(pm_stock[-1],x[0])

        return pm_list

    def jsi_filter(self,x,dt):
        y1 = self.iir_1st( x,0.0,1.0,0.45,dt)
        y2 = self.iir_1st(y1,1.0,2.0,7.00,dt)
        y3 = self.iir_1st(y2,4.0,8.0,7.00,dt)
        y4 = self.iir_1st(y3,0.25,0.5,7.00,dt)
        y5 = self.iir_2nd(y4,0.9,11.0,dt)
        y6 = 1.409*y5
        return y6

    def iir_1st(self,x,a,b,f,dt):
        w = 2*np.pi*f
        a0 = w+2*b/dt
        a1 = w-2*b/dt
        b0 = w*a+2/dt
        b1 = w*a-2/dt

        y = lfilter([b0,b1],[a0,a1],x)
        return y

    def iir_2nd(self,x,h,f,dt):
        w = 2*np.pi*f
        a0 = 12/dt**2 +12*h*w/dt + w**2
        a1 = 10*w**2 - 24/dt**2
        a2 = 12/dt**2 -12*h*w/dt + w**2
        b0 = w**2
        b1 = 10*w**2
        b2 = w**2

        y = lfilter([b0,b1,b2],[a0,a1,a2],x)
        return y

    def plot_all(self,start=0,end=60,to_end=False):
        if to_end:
            start = self.tim[0]
            end = self.tim[-1]
            amax = np.amax(self.I)
        else:
            ns = int(start/self.dt)
            ne = int(end/self.dt)
            amax = np.amax(self.I[ns:ne])

        fig = plt.figure()

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax1 = fig.add_subplot(111)
        ax1.plot(self.tim,self.I,color='k',lw=1)
        ax1.set_xlim(start,end)
        ax1.set_ylim(0,1.2*amax)
        ax1.set_xlabel("time [s]")
        ax1.text(0.02,0.85,"Ir",transform=ax1.transAxes)
        plt.show()
