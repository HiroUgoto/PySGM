import numpy as np
import matplotlib.pyplot as plt
import math

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

def response_spectrum_FD(wave,period,dt):
    h = 0.05
    img = 0.0 + 1.0j

    ntim = len(wave)
    sa_list = []
    sv_list = []
    sd_list = []

    w = np.fft.fft(wave)/(ntim//2)
    freq = np.fft.fftfreq(ntim,d=dt)
    omega = 2.0*np.pi*freq[1:ntim//2]

    acc_FD = np.zeros_like(w)
    vel_FD = np.zeros_like(w)
    disp_FD = np.zeros_like(w)

    for p in period:
        omega0 = 2*np.pi/p

        resp = omega**2 / (-omega**2 + 2*img*h*omega0*omega + omega0**2)
        coef = img*omega

        acc_FD[1:ntim//2] = w[1:ntim//2] * resp
        vel_FD[1:ntim//2] = acc_FD[1:ntim//2] / coef
        disp_FD[1:ntim//2] = vel_FD[1:ntim//2] / coef

        acc = np.real(np.fft.ifft(acc_FD))*ntim + wave
        vel = np.real(np.fft.ifft(vel_FD))*ntim
        disp = np.real(np.fft.ifft(disp_FD))*ntim

        sa = np.max(np.abs(acc))
        sv = np.max(np.abs(vel))
        sd = np.max(np.abs(disp))

        sa_list += [sa]
        sv_list += [sv]
        sd_list += [sd]

    return np.array(sa_list), np.array(sv_list), np.array(sd_list)

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

def response_spectrum_pseudo_FD(ew,ns,period,dt):
    h = 0.05
    img = 0.0 + 1.0j

    ntim = min(len(ew),30000)
    psv_list = []

    ew_FD = np.fft.fft(ew)/(ntim//2)
    ns_FD = np.fft.fft(ns)/(ntim//2)
    freq = np.fft.fftfreq(ntim,d=dt)
    omega = 2.0*np.pi*freq[1:ntim//2]

    acc_ew_FD = np.zeros_like(ew_FD)
    acc_ns_FD = np.zeros_like(ns_FD)

    for p in period:
        omega0 = 2*np.pi/p

        resp = omega**2 / (-omega**2 + 2*img*h*omega0*omega + omega0**2)
        acc_ew_FD[1:ntim//2] = ew_FD[1:ntim//2] * resp
        acc_ns_FD[1:ntim//2] = ns_FD[1:ntim//2] * resp

        acc_ew = np.real(np.fft.ifft(acc_ew_FD))*ntim + ew
        acc_ns = np.real(np.fft.ifft(acc_ns_FD))*ntim + ns

        acc_index = np.argmax(acc_ew**2 + acc_ns**2)
        sa = math.sqrt(acc_ew[acc_index]**2 + acc_ns[acc_index]**2)
        psv = sa*p/(2.0*np.pi)

        psv_list += [psv]

    peak_index = np.argmax(np.array(psv_list))
    peak_psv = psv_list[peak_index]
    peak_period = period[peak_index]

    return np.array(psv_list),peak_psv,peak_period

def response_spectrum_peak_period_sv(ew,ns,period,dt):
    ntim = min(len(ew),30000)
    sv_list = []

    for p in period:
        sv = 0.0

        acc_ew, vel_ew, disp_ew = response_1dof_full(ew[0:ntim],p,dt)
        acc_ns, vel_ns, disp_ns = response_1dof_full(ns[0:ntim],p,dt)

        vel_index = np.argmax(vel_ew**2 + vel_ns**2)
        sv = math.sqrt(vel_ew[vel_index]**2 + vel_ns[vel_index]**2)

        sv_list += [sv]

    peak_index = np.argmax(np.array(sv_list))
    peak = sv_list[peak_index]
    peak_period = period[peak_index]

    return peak,peak_period

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

def calc_SI_FD(ew,ns,dt):
    h = 0.20
    img = 0.0 + 1.0j

    ntim = min(len(ew),30000)
    period,dp = np.linspace(0.1,2.5,25,retstep=True)
    rot = np.linspace(0,np.pi,8,endpoint=False)
    SI_list = np.zeros_like(rot)

    ew_FD = np.fft.fft(ew)/(ntim//2)
    ns_FD = np.fft.fft(ns)/(ntim//2)
    freq = np.fft.fftfreq(ntim,d=dt)
    omega = 2.0*np.pi*freq[1:ntim//2]

    vel_ew_FD = np.zeros_like(ew_FD)
    vel_ns_FD = np.zeros_like(ns_FD)

    for p in period:
        omega0 = 2*np.pi/p

        resp = -img*omega / (-omega**2 + 2*img*h*omega0*omega + omega0**2)
        vel_ew_FD[1:ntim//2] = ew_FD[1:ntim//2] * resp
        vel_ns_FD[1:ntim//2] = ns_FD[1:ntim//2] * resp

        vel_ew = np.real(np.fft.ifft(vel_ew_FD))*ntim
        vel_ns = np.real(np.fft.ifft(vel_ns_FD))*ntim

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

def response_1dof_full_FD(wave,period,dt,h=0.05):
    n = len(wave)
    w = np.fft.fft(wave)/(n/2)
    freq = np.fft.fftfreq(len(wave),d=dt)
    omega0 = 2*np.pi/period

    acc_FD = np.zeros_like(w)
    vel_FD = np.zeros_like(w)
    disp_FD = np.zeros_like(w)
    for (i,f) in enumerate(freq):
        if f > 0.0:
            omega = 2.0*np.pi*f
            img = 0.0 + 1.0j
            resp = omega**2 / (-omega**2 + 2*img*h*omega0*omega + omega0**2)
            coef = img*omega

            acc_FD[i] = w[i] * resp
            vel_FD[i] = acc_FD[i] / coef
            disp_FD[i] = vel_FD[i] / coef

    acc = np.real(np.fft.ifft(acc_FD))*n
    vel = np.real(np.fft.ifft(vel_FD))*n
    disp = np.real(np.fft.ifft(disp_FD))*n

    return acc+wave, vel, disp
