import numpy as np
import copy
from . import vector


def sin_wave(dt,nt,fp,t0=0.0):
    tim = np.linspace(0,dt*nt,num=nt,endpoint=False)
    sig = np.sin(2.0*np.pi*fp*(tim-t0))
    h0 = np.heaviside(tim-t0,1)
    h1 = np.heaviside(tim-(t0+1.0/fp),1)
    wave = sig*(h0-h1)

    header = {'code':'sin_wave',        'record_time':'NA','lat':0.0,'lon':0.0,'ntim':nt}

    v = vector.vector(header,tim,wave)

    return v

def ricker_wave(dt,nt,fp,t0=0.0):
    tim = np.linspace(0,dt*nt,num=nt,endpoint=False)
    tp = 1.0/fp+t0
    wave = -(1.0-2.0*(np.pi*fp)**2*(tim-tp)**2)*np.exp(-(tim-tp)**2*(np.pi*fp)**2)

    header = {'code':'sin_wave',        'record_time':'NA','lat':0.0,'lon':0.0,'ntim':nt}

    v = vector.vector(header,tim,wave)

    return v

def zero_wave(dt,nt):
    tim = np.linspace(0,dt*nt,num=nt,endpoint=False)
    wave = np.zeros_like(tim)
    header = {'code':'',        'record_time':'NA','lat':0.0,'lon':0.0,'ntim':nt}
    v = vector.vector(header,tim,wave)
    return v

def zeros_like(vector):
    v = zero_wave(vector.dt,vector.header["ntim"])
    return v
