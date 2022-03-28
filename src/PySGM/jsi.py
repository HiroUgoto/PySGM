import numpy as np
import math
from scipy.fftpack import fft, ifft, fftfreq

#######################################
##  JMA Seismic intensity            ##
#######################################
def jsi(ew,ns,ud,dt):

    ntim = len(ew)
    fs = 1.0/dt

    ew_filtered = jsi_filter(ew,dt)
    ns_filtered = jsi_filter(ns,dt)
    ud_filtered = jsi_filter(ud,dt)

    wave = np.abs(np.power(ew_filtered,2) + np.power(ns_filtered,2) + np.power(ud_filtered,2))
    wave = sorted(wave, reverse=True)

    I = math.log10(wave[int(0.3*fs)-1]) + 0.94

    return I


def jsi_filter(wave,dt):

    ntim = len(wave)

    freq = np.abs(fftfreq(ntim,d=dt))
    wave_f = fft(wave)

    w1 = np.sqrt(1.0 / freq[1:])
    w1 = np.concatenate(([0],w1))

    x = freq / 10.0
    w2 = 1.0 / np.sqrt(1.0 + 0.694*np.power(x,2) + 0.241*np.power(x,4) + 0.0557*np.power(x,6) \
                    + 0.009664*np.power(x,8) + 0.00134*np.power(x,10) + 0.000155*np.power(x,12))

    w3 = np.sqrt(1.0 - np.exp(-np.power(freq/0.5,3)))

    wave_f = w1*w2*w3 * wave_f
    wave_filtered = ifft(wave_f)

    return wave_filtered
