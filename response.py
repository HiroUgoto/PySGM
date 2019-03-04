# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
import math

try:
    from numba import jit
except ImportError:
    def jit(*args, **_kwargs):
        if len(args) > 0 and hasattr(args[0], "__call__"):
            return args[0]
        else:
            def _(func):
                return func
            return _

#######################################
##          Resp    Method           ##
#######################################
def response_spectrum(wave,period,dt):

    sa_list = []
    sv_list = []
    sd_list = []
    for p in period:
        sa,sv,sd = response_1dof(wave,p,dt)
        sa_list += [sa]
        sv_list += [sv]
        sd_list += [sd]

    return np.array(sa_list), np.array(sv_list), np.array(sd_list)

@jit
def response_spectrum_pseudo(ew,ns,period,dt):
    ntim = min(len(ew),30000)
    psv_list = []

    for p in period:
        sa = 0.0
        sa_rot = 0

        acc_ew, vel_ew, disp_ew = response_1dof_full(ew[0:ntim],p,dt)
        acc_ns, vel_ns, disp_ns = response_1dof_full(ns[0:ntim],p,dt)

        acc_index = np.argmax(acc_ew**2 + acc_ns**2)
        sa = math.sqrt(acc_ew[acc_index]**2 + acc_ns[acc_index]**2)
        psv = sa*p/(2.0*np.pi)

        psv_list += [psv]

    peak_index = np.argmax(np.array(psv_list))
    peak_psv = psv_list[peak_index]
    peak_period = period[peak_index]

    return np.array(psv_list),peak_psv,peak_period

@jit
def response_spectrum_max(ew,ns,period,dt):

    ntim = min(len(ew),30000)

    sa_list = []
    sv_list = []
    sd_list = []

    sa_rot_list = []
    sv_rot_list = []
    sd_rot_list = []

    for p in period:
        sa = 0.0
        sv = 0.0
        sd = 0.0
        sa_rot = 0
        sv_rot = 0
        sd_rot = 0

        acc_ew, vel_ew, disp_ew = response_1dof_full(ew[0:ntim],p,dt)
        acc_ns, vel_ns, disp_ns = response_1dof_full(ns[0:ntim],p,dt)

        acc_index = np.argmax(acc_ew**2 + acc_ns**2)
        sa = math.sqrt(acc_ew[acc_index]**2 + acc_ns[acc_index]**2)
        sa_rot = math.atan(acc_ew[acc_index]/acc_ns[acc_index])/math.pi*180

        vel_index = np.argmax(vel_ew**2 + vel_ns**2)
        sv = math.sqrt(vel_ew[vel_index]**2 + vel_ns[vel_index]**2)
        sv_rot = math.atan(vel_ew[vel_index]/vel_ns[vel_index])/math.pi*180

        disp_index = np.argmax(disp_ew**2 + disp_ns**2)
        sd = math.sqrt(disp_ew[disp_index]**2 + disp_ns[disp_index]**2)
        sd_rot = math.atan(disp_ew[disp_index]/disp_ns[disp_index])/math.pi*180

        sa_list += [sa]
        sv_list += [sv]
        sd_list += [sd]

        sa_rot_list += [sa_rot]
        sv_rot_list += [sv_rot]
        sd_rot_list += [sd_rot]

    return np.array(sa_list), np.array(sv_list), np.array(sd_list), \
        np.array(sa_rot_list), np.array(sv_rot_list), np.array(sd_rot_list)


@jit
def calc_SI(ew,ns,dt):
    ntim = min(len(ew),30000)
    period,dp = np.linspace(0.1,2.5,25,retstep=True)
    rot = np.linspace(0,np.pi,8,endpoint=False)
    SI_list = np.zeros_like(rot)

    for p in period:
        acc_ew, vel_ew, disp_ew = response_1dof_full(ew[0:ntim],p,dt,h=0.2)
        acc_ns, vel_ns, disp_ns = response_1dof_full(ns[0:ntim],p,dt,h=0.2)

        if p == period[0]:
            d = 0.5*dp
        elif p == period[-1]:
            d = 0.5*dp
        else:
            d = dp

        for index,r in enumerate(rot):
            sv = np.max(np.abs(vel_ns*np.cos(r) + vel_ew*np.sin(r)))
            SI_list[index] += sv*d

    SI = np.max(SI_list)/2.4
    return SI


@jit
def response_1dof(wave,period,dt,h=0.05):
    ntim = len(wave)
    beta = 1.0 / 6.0

    omega = 2.0*math.pi / period

    acc = 0.0
    vel = 0.0
    disp = 0.0

    sa = 0.0
    sv = 0.0
    sd = 0.0

    dt2 = dt*dt
    o2 = omega*omega
    oh = omega*h
    ohdt = oh*dt
    odt2 = o2*dt2

    for it in range(1,ntim):
        acc1 = -(wave[it] + 2.0*oh*(vel+0.5*acc*dt) \
                        + o2*(disp+vel*dt+(0.5-beta)*acc*dt2)) \
                        / (1.0 + ohdt + beta*odt2)
        vel1 = vel + 0.5*(acc+acc1)*dt
        disp1 = disp + vel*dt + (0.5-beta)*acc*dt2 + beta*acc1*dt2

        acc_abs = acc1 + wave[it]

        acc = acc1
        vel = vel1
        disp = disp1

        sa = max(sa,abs(acc_abs))
        sv = max(sv,abs(vel))
        sd = max(sd,abs(disp))

    return sa, sv, sd

@jit
def response_1dof_full(wave,period,dt,h=0.05):
    ntim = len(wave)
    beta = 1.0 / 6.0

    omega = 2.0*math.pi / period

    acc = 0.0
    vel = 0.0
    disp = 0.0

    sa = 0.0
    sv = 0.0
    sd = 0.0

    dt2 = dt*dt
    o2 = omega*omega
    oh = omega*h
    ohdt = oh*dt
    odt2 = o2*dt2

    acc_list = [acc]
    vel_list = [vel]
    disp_list = [disp]

    for it in range(1,ntim):
        acc1 = -(wave[it] + 2.0*oh*(vel+0.5*acc*dt) \
                        + o2*(disp+vel*dt+(0.5-beta)*acc*dt2)) \
                        / (1.0 + ohdt + beta*odt2)
        vel1 = vel + 0.5*(acc+acc1)*dt
        disp1 = disp + vel*dt + (0.5-beta)*acc*dt2 + beta*acc1*dt2

        acc_abs = acc1 + wave[it]

        acc = acc1
        vel = vel1
        disp = disp1

        acc_list += [acc_abs]
        vel_list += [vel1]
        disp_list += [disp1]

    return np.array(acc_list), np.array(vel_list), np.array(disp_list)
