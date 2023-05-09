import numpy as np
import math

def demand_curve(wave,period,ductility,dt):
    kh_min = 20.0
    kh_max = 0.01

    demand_kh = []
    for p in period:
        duc_min = elasto_plastic_response(wave,p,kh_min,dt)
        duc_max = elasto_plastic_response(wave,p,kh_max,dt)

        dc_kh = dsc_binary_search(wave,p,ductility,kh_min,duc_min,kh_max,duc_max,dt)
        demand_kh += [dc_kh]
        print("   ",p,dc_kh)

    return demand_kh


def dsc_binary_search(wave,period,ductility,kh_min,duc_min,kh_max,duc_max,dt):
    eps = 1.e-4
    kh_mid = (kh_min + kh_max) * 0.5
    duc_mid = elasto_plastic_response(wave,period,kh_mid,dt)
    
    if abs(ductility - duc_mid) > eps:
        if (abs(duc_max - duc_mid) > eps) and (abs(duc_mid - duc_min) > eps):
            if duc_mid < ductility:
                # print("   * kh:",kh_mid,"(",duc_mid,"-",ductility,"-",duc_max,")")
                kh_mid = dsc_binary_search(wave,period,ductility,kh_mid,duc_mid,kh_max,duc_max,dt)
            else:
                # print("   * kh:",kh_mid,"(",duc_min,"-",ductility,"-",duc_mid,")")
                kh_mid = dsc_binary_search(wave,period,ductility,kh_min,duc_min,kh_mid,duc_mid,dt)

    return kh_mid

# --- Wrapper --- #
def elasto_plastic_response(wave,period,kh,dt):
    # return elasto_plastic_response_bilinear(wave,period,kh,dt)
    return elasto_plastic_response_clough(wave,period,kh,dt)

# --- Elasto-plastic response of 1DOF system --- #
def elasto_plastic_response_bilinear(wave,period,kh,dt):
    ntim = len(wave)

    # data refinement
    nref = 5
    ddt = dt / nref 

    wave_norm = wave / (nref-1)

    # Set parameters
    h = 0.10
    kr = 0.0

    omega = 2.0*math.pi / period
    k1 = np.power(omega,2)
    k2 = k1*kr
    g = 980.0
    uy = g*kh/k1

    c0 = -np.power(ddt,2) / (1. + ddt*omega*h)
    c1 = 2. / (1. + ddt*omega*h)
    c2 = - (1. - ddt*omega*h) / (1. + ddt*omega*h)

    disp = 0.0
    disp1 = 0.0
    disp2 = 0.0
    q = 0.0

    umax =  uy
    umin = -uy
    qp =  (k1-k2)*uy
    qm = -(k1-k2)*uy
    
    disp_max = 0.0
    for it in range(1,ntim-1):
        for j in range(nref):
            fine_wave = wave_norm[it]*(nref-1-j) + wave_norm[it+1]*j

            disp = c0*(fine_wave + q) + c1*disp1 + c2*disp2

            du = disp - disp1
            q_pre = q + k1*du
            if du > 0:
                q_sur = qp + k2*disp
                if q_pre > q_sur:
                    q = q_sur
                else:
                    q = q_pre
                    
            else:
                q_sur = qm + k2*disp
                if q_pre < q_sur:
                    q = q_sur
                else:
                    q = q_pre

            disp2 = disp1
            disp1 = disp

            disp_max = max(disp_max,abs(disp))
            # print(disp,q)

    ductility = disp_max / uy
    return ductility


def elasto_plastic_response_clough(wave,period,kh,dt):
    ntim = len(wave)

    # data refinement
    nref = 5
    ddt = dt / nref 

    wave_norm = wave / (nref-1)

    # Set parameters
    h = 0.10
    kr = 0.0

    omega = 2.0*math.pi / period
    k1 = np.power(omega,2)
    k2 = k1*kr
    g = 980.0
    uy = g*kh/k1

    c0 = -np.power(ddt,2) / (1. + ddt*omega*h)
    c1 = 2. / (1. + ddt*omega*h)
    c2 = - (1. - ddt*omega*h) / (1. + ddt*omega*h)

    disp = 0.0
    disp1 = 0.0
    disp2 = 0.0
    vel1 = 0.0
    acc1 = 0.0
    q = 0.0

    qp =  (k1-k2)*uy
    qm = -(k1-k2)*uy

    umax =  uy
    umin = -uy
    qmax =  g*kh
    qmin = -g*kh
    u0 = 0.0

    disp_max = 0.0
    for it in range(1,ntim-1):
        for j in range(nref):
            fine_wave = wave_norm[it]*(nref-1-j) + wave_norm[it+1]*j

            disp = c0*(fine_wave + q) + c1*disp1 + c2*disp2

            du = disp - disp1
            q_pre = q + k1*du
            if du > 0:
                if q_pre <= 0.:
                    q = q_pre
                    u0 = disp1 - q/k1
                elif disp < umax:
                    kc = qmax / (umax - u0)
                    q = kc*(disp - u0)
                else:
                    q = qp + k2*disp
                    umax = disp
                    qmax = q
            else:
                if q_pre >= 0.:
                    q = q_pre
                    u0 = disp1 - q/k1
                elif disp > umin:
                    kc = qmin / (umin - u0)
                    q = kc*(disp - u0)
                else:
                    q = qm + k2*disp
                    umin = disp
                    qmin = q

            disp2 = disp1
            disp1 = disp

            disp_max = max(disp_max,abs(disp))
            # print(disp,q)
    
    ductility = disp_max / uy
    return ductility
