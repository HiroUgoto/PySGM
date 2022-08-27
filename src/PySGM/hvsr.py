import numpy as np
import scipy.signal
import matplotlib.pyplot as plt

# from . import vector

#-----------------------------------------------------------------#
def detrend(data):
    v = scipy.signal.detrend(data)
    return v

def cos_taper(data):
    taper_data = []
    for line in data:
        w = scipy.signal.tukey(len(line),0.05)
        taper_data += [w*line]
    return np.array(taper_data)

def rms(data):
    r = np.mean(data**2,axis=-1)
    return r

def smoothing(fourier_data,nw):
    if nw%2 == 0:
        nwt = nw + 1
    else:
        nwt = nw

    nw2 = int((nwt-1)/2)
    w = scipy.signal.parzen(nwt)

    if fourier_data.ndim == 1:
        line = fourier_data
        a = np.r_[line[nw2:0:-1],line,line[0],line[-1:-nw2:-1]]
        smooth_data = np.convolve(w/w.sum(),a,mode='valid')
    else:
        smooth_data = []
        for line in fourier_data:
            a = np.r_[line[nw2:0:-1],line,line[0],line[-1:-nw2:-1]]
            smooth_data += [np.convolve(w/w.sum(),a,mode='valid')]

    return np.array(smooth_data)

#-----------------------------------------------------------------#
def segment_selection(vector,tseg=40.96,mseg=10,print_flag=True,plot_flag=False):

    vector_bp = vector.bandpass(0.05,50)

    ud = vector_bp.ud
    ns = vector_bp.ns
    ew = vector_bp.ew

    fs = 1.0/vector.dt
    nt = int(tseg*fs)
    nseg = len(ud)//nt
    ns = 3

    v_list = []
    rms_stack = np.zeros(2*nseg)

    for d in [ud,ns,ew]:
        v_raw = np.resize(d,[nseg,nt]).copy()
        v_raw2 = np.resize(np.roll(d,int(nt/2)),[nseg,nt]).copy()
        v = detrend(np.vstack((v_raw,v_raw2)))
        v_list += [cos_taper(v)]
        rms_stack += rms(v)

    nseg = nseg*2
    total_rms = rms_stack.mean()/ns
    rms_index = np.argsort(rms_stack)

    total_rms = rms_stack.mean()/ns
    selected_rms = rms_stack[rms_index[0:mseg]].mean()/ns

    segment_data = []
    for v in v_list:
        sv = v[rms_index[0:mseg],:].copy()
        segment_data += [sv]

    if print_flag:
        print("------------------------------------")
        print("segment selection")
        print("+ selected segment :",mseg,"/",nseg)
        print("+ rms ratio        :",selected_rms,"/",total_rms)
        print("------------------------------------")

    if plot_flag:
        fig = plt.figure()
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        ax1.plot(segment_data[0][0],color='k',lw=1)
        ax2.plot(segment_data[0][1],color='k',lw=1)
        ax3.plot(segment_data[0][2],color='k',lw=1)

        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        ax3.set_xlabel("time samples")

        ax1.text(0.02,0.85,"UD",transform=ax1.transAxes)
        ax2.text(0.02,0.85,"H1",transform=ax2.transAxes)
        ax3.text(0.02,0.85,"H2",transform=ax3.transAxes)

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
        plt.show()

    return segment_data

#-----------------------------------------------------------------#
def hv_spectra(vector,tseg=40.96,mseg=10,band_width=0.2,plot_flag=True):
    segment_data = segment_selection(vector,tseg,mseg,print_flag=False)

    fs = 1.0/vector.dt
    ns = len(segment_data)

    nt = len(segment_data[0][0,:])
    nyq = int(nt/2)
    freq = np.fft.fftfreq(nt,d=1/fs)

    df = freq[1]-freq[0]
    nw = int(band_width/df)

    X = []
    for v in segment_data:
        F = np.fft.fft(v)
        X += [smoothing(F,nw)]

    UD = np.mean(np.abs(np.conj(X[0])*X[0]),axis=0)
    H1 = np.mean(np.abs(np.conj(X[1])*X[1]),axis=0)
    H2 = np.mean(np.abs(np.conj(X[2])*X[2]),axis=0)

    hv = smoothing(np.sqrt(H1+H2)/np.sqrt(UD),nw)

    if plot_flag:
        ud_seg = np.abs(np.conj(X[0])*X[0])
        h1_seg = np.abs(np.conj(X[1])*X[1])
        h2_seg = np.abs(np.conj(X[2])*X[2])
        hv_seg = []
        for i in range(mseg):
            hv_seg += [smoothing(np.sqrt(h1_seg[i,:]+h2_seg[i,:])/np.sqrt(ud_seg[i,:]),nw)]

        plt.figure()
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("frequency (Hz)")
        plt.ylabel("H/V spectrum")

        for i in range(mseg):
            plt.plot(freq[0:nyq],hv_seg[i][0:nyq],lw=1,color="gray")

        plt.plot(freq[0:nyq],hv[0:nyq],color="red")
        plt.grid()
        plt.show()


    return freq[0:nyq],hv[0:nyq]
