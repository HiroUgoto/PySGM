import numpy as np

def Si_Midorikawa_1999(R,Mw,D,type,output="PGV"):
    # R: fault distance [km]
    # Mw: Moment magnitude
    # D: hypocenter depth [km]
    # type: EQ type; crust 0, inter-plate 1, intra-plate 2

    if output == "PGA":
        c = 0.0055 * 10**(0.50*Mw)
        if type == 0:
            d = 0.0
        elif type == 1:
            d = 0.01
        else:
            d = 0.22
        b = 0.50*Mw + 0.0043*D + d +0.61

        A = 10**(b - np.log10(R+c) - 0.003*R)
        std = 10**0.27
        A_sp = A*std
        A_sn = A/std

    elif output == "PGV":
        c = 0.0028 * 10**(0.50*Mw)
        if type == 0:
            d = 0.0
        elif type == 1:
            d = -0.02
        else:
            d = 0.12
        b = 0.58*Mw + 0.0038*D + d -1.29
        A = 10**(b - np.log10(R+c) - 0.002*R)
        std = 10**0.23
        A_sp = A*std
        A_sn = A/std

    else:
        print("Input Error in Si_Midorikawa_1999")

    return A,A_sn,A_sp

def Morikawa_Fujiwara_2013(R,Mw,type,Vs30=350,D1400=250,output="PGV"):
    # R: fault distance [km]
    # Mw: Moment magnitude
    # Vs30: Averaged S-wave velocity to 30m [m/s]
    # D1400: Top depth of Vs=1400m/s [m]
    # type: EQ type; crust 0, inter-plate 1, intra-plate 2

    if output == "PGA":
        Mw0 = min(Mw,8.2)
        Mw1 = 16.0
        a = -0.0321
        d = 0.011641
        e = 0.5

        if type == 0:
            b = -0.005315
            c = 7.0830
        elif type == 1:
            b = -0.005042
            c = 7.1181
        else:
            b = -0.005605
            c = 7.5035

        pd = 0.0663
        Dmin = 100.00
        Gd = 10**(pd*np.log10(max(Dmin,D1400)/250))

        ps = -0.3709
        Vsmax = 1950.0
        Gs = 10**(ps*np.log10(min(Vsmax,Vs30)/350))

        A = 10**(a*(Mw0 - Mw1)**2 +b*R + c - np.log10(R+d*10**(e*Mw0)))
        std = 10**0.3761

        A = A*Gd*Gs
        A_sp = A*std
        A_sn = A/std

    elif output == "PGV":
        Mw0 = min(Mw,8.2)
        Mw1 = 16.0
        a = -0.0325
        d = 0.002266
        e = 0.5
        if type == 0:
            b = -0.002654
            c = 5.6952
        elif type == 1:
            b = -0.002408
            c = 5.6026
        else:
            b = -0.003451
            c = 6.0030

        pd = 0.2317
        Dmin = 60.00
        Gd = 10**(pd*np.log10(max(Dmin,D1400)/250))

        ps = -0.5546
        Vsmax = 1100.0
        Gs = 10**(ps*np.log10(min(Vsmax,Vs30)/350))

        A = 10**(a*(Mw0 - Mw1)**2 +b*R + c - np.log10(R+d*10**(e*Mw0)))
        std = 10**0.3399

        A = A*Gd*Gs
        A_sp = A*std
        A_sn = A/std

    else:
        print("Input Error in Morikawa_Fujiwara_2013")

    return A,A_sn,A_sp
