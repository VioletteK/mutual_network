"""
Copyright (c) 2016, Domenico GUARINO
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GUARINO BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import os
import glob
import gc
import json
import pickle
import sys
import resource
import collections
import random
import warnings
import itertools
from functools import cmp_to_key
from itertools import zip_longest # analysis
import numpy as np
from scipy.signal import savgol_filter

################################
import matplotlib
matplotlib.use('Agg') # to be used when DISPLAY is undefined
################################
# matplotlib.rc('image', cmap='viridis')
# matplotlib.rc('image', cmap='Reds')
matplotlib.rc('image', cmap='summer')
################################

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.collections import LineCollection
import matplotlib.colors as mpcolors
import matplotlib.cm as mpcm

from neo.core import AnalogSignal # analysis
import quantities as pq

# SciPy related
import scipy.signal as signal
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist, jaccard, hamming
import scipy.stats as stats
from scipy.stats import ttest_ind
from scipy.stats import shapiro
from scipy.stats import hypergeom
from scipy.stats import entropy
from scipy import sparse
from scipy.sparse.linalg import spsolve
import scipy
import scipy.linalg
import scipy.signal as signal
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.optimize import curve_fit


from pyNN.utility.plotting import Figure, Panel # analysis
from sklearn.metrics import mutual_info_score


# ------------------------------------------------------------------------------

class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))

# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------


def analyse(params, folder, addon='', removeDataFile=False):
    print("\nAnalysing data...")

    # populations key-recorders match
    populations = {}
    for popKey,popVal in params['Populations'].items():
        # if popKey != 'ext': # if is not the drive, to be dropped
        #     if popKey in params['Recorders']:
        #         populations[popKey] = list(params['Recorders'][popKey].keys())
        # we do the analysis on what we recorded
        if popKey in params['Recorders']:
            populations[popKey] = list(params['Recorders'][popKey].keys())

    ###################################
    # iteration over populations and selective plotting based on params and available recorders
    for key,rec in populations.items():
        print("\n\nAnalysis for:",key)

        # assuming 2D structure to compute the edge N
        n = 0
        if isinstance(params['Populations'][key]['n'], dict):
            n = int(params['Populations'][params['Populations'][key]['n']['ref']]['n'] * params['Populations'][key]['n']['ratio'])
        else:
            n = int(params['Populations'][key]['n'])
        edge = int(np.sqrt(n))

        # trials
        for trial_id,trial in enumerate(params['trials']):
            print("\n"+trial['name'])

            for itrial in range(trial['count']):
                print("trial #",itrial)
                timeslice_start = params['run_time'] * trial_id + params['Analysis']['transient'] # to avoid initial transient
                timeslice_end   = params['run_time'] * (trial_id+itrial+1)
                print("trial-based slicing window (ms):", timeslice_start, timeslice_end)

                # get data
                print("from file:",key+addon+'_'+trial['name']+str(itrial))
                neo = pickle.load( open(folder+'/'+key+addon+'_'+trial['name']+str(itrial)+'.pkl', "rb") )
                data = neo.segments[0]

                # getting and slicing data
                # continuous variables (Vm, Gsyn, W) are sliced just according to the dt
                if 'w' in rec:
                    w = data.filter(name = 'w')[0]#[timeslice_start:timeslice_end]
                if 'v' in rec:
                    vm = data.filter(name = 'v')[0]#[timeslice_start:timeslice_end]
                if 'gsyn_exc' in rec:
                    gexc = data.filter(name = 'gsyn_exc')[0]#[timeslice_start:timeslice_end]
                if 'gsyn_inh' in rec:
                    ginh = data.filter(name = 'gsyn_inh')[0]#[timeslice_start:timeslice_end]
                # discrete variables (spiketrains) are sliced according to their time signature
                if 'spikes' in rec:
                    # spiketrains = data.spiketrains[ (data.spiketrains[:]>=timeslice_start) & (data.spiketrains[:]<=timeslice_end) ]
                    spiketrains = []
                    for spiketrain in data.spiketrains:
                        spiketrains.append(spiketrain[ (spiketrain>=timeslice_start) & (spiketrain<=timeslice_end) ])

                # Get cell indexes and ids
                cell_coords = []
                cell_indexes = []
                cell_ids = []
                with open(folder+'/'+key+addon+'_'+trial['name']+str(itrial)+'_positions.txt', 'r') as posfile:
                    print("... getting cell indexes and ids")
                    lines = posfile.readlines()
                    posfile.close()
                    for line in lines:
                        cell_coords.append( [int(float(i)) for i in line.split(' ')[:4]] ) # not including 4
                    # print(cell_coords) # id, idx, x, y
                    #[ ..., [13287, 4064, 63, 32], [13288, 4065, 63, 33], ... ]
                    cell_coords = np.array(cell_coords)
                    cell_indexes = (cell_coords[:,1]).tolist()
                    cell_ids = cell_coords[:,0]

                panels = []
                # return list for param search
                # scores = []
                # scores.append(0)   # 0. Spike count
                # scores.append(0.0) # 1. Inter-Spike Interval
                # scores.append(0.0) # 2. Coefficient of Variation
                # scores.append("")  # 3. additional label for adapting ISI
                # scores.append(0.0) # 4. Firing rate
                # scores.append(0.0) # 5. mean power in the Delta range (0.1-4 Hz)
                # scores.append(0.0) # 6. mean power in the Spindles range (7-14 Hz)
                # scores.append(0.0) # 7. Cross-Correlation


                ###################################
                if 'Vm' in params['Analysis'] and params['Analysis']['Vm'] and key in params['Analysis']['Vm']:
                    print('Vm')
                    dt = params['dt']
                    run_time = params['run_time']
                    size = np.sqrt(params['Populations']['py']['n'])
                    injection_start,injection_end = params['Injections']['py']['start']
                    interval = 50



                    x = params['Recorders']['py']['v']['x']
                    y = params['Recorders']['py']['v']['y']
                    window = params['Recorders']['py']['v']['size']+1
                    X, Y = np.meshgrid([i for i in range(window)],[i for i in range(window)])
                    list_coord = [(x+i)*size+(y+j) for j in range(window) for i in range(window)]
                    #Here is a list of the ancient coord at the index of their new
                    #to obtain the new coord of a neurone : list_coord.index(coord_neurone)

                    injection_points=params['Injections']['py']['cellidx']
                    min_point, max_point = min(injection_points)//size,max(injection_points)//size

                    #the 'center' of the annulus taken as the central point of all the cells
                    ref_neurone = [np.floor(np.mean([min_point,max_point])) for i in range(2)]

                    #Gaussian Blur
                    Kernel1 = 1/16*np.array([
                    [1,2,1],
                    [2,4,2],
                    [1,2,1]
                    ])

                    def accentuation1(x):
                        ''' this function amplifies the high values and decreases the low ones '''
                        return 0.2*(np.tanh(10*x-2)+1)


                    def accentuation2(x):
                        return 0.2*(np.tanh((10*x-2)**3)+1)


                    Recorded_cell = np.zeros((window, window))
                    #The MI matrix

                    Vm_t = np.zeros((window,window))
                    #The Vm matrix


                    max_time = run_time-injection_start-interval

                    Time_delay = np.arange(-100,800,10)
                    #the windows where we calculate the MI have to intersect
                    V=len(vm)
                    vm=vm.T

                    indice = list_coord.index(ref_neurone[0]*size+ref_neurone[1])
                    vm_base = vm[indice][int(injection_start/dt):int((injection_start+interval)/dt)]
                    #list of the vm values of the central neuron beetween the beginning of the injection and beginnin+interval

                    c_X,xedges = np.histogram(vm_base,500,range=(-90.,-40.))
                    #every interval (ms) we calculate the MI
                    for time_delay in Time_delay:
                        for i in range(window**2):
                            vm_neurone = vm[i][int((injection_start+time_delay)/dt):int((injection_start+time_delay+interval)/dt)]
                            c_Y,xedges = np.histogram(vm_neurone,500,range=(-90.,-40.))#,density=True)
                            Recorded_cell[i%window][i//window]= mutual_info_score(c_X,c_Y)
                            Vm_t[i%window][i//window] = np.mean(vm_neurone)


                        Recorded_cell[indice%window][indice//window]=0.4
                        #Set up this value to see the central neuron

                        Recorded_cell1=signal.convolve2d(accentuation1(Recorded_cell),Kernel1, mode='same')
                        #This is the colormap of accentuated and then blurred MI to obtain a good contouring

                        fig=plt.figure()
                        plt.imshow(Recorded_cell, cmap = 'inferno',interpolation='none')
                        plt.clim([0,0.4])
                        plt.contour(Recorded_cell1, 5)
                        plt.colorbar()
                        tmin=float(i*dt)
                        tmax=float((i+interval)*dt)
                        plt.title('MI '+str(injection_start+time_delay)+','+str(injection_start+time_delay+interval))
                        fig.savefig(folder+'/Tau='+str(params['Populations']['py']['cellparams']['tau_w'])+'Convolution'+str(time_delay)+'.png')
                        plt.close()
                        fig.clf()


















                # for systems with low memory :)
                if removeDataFile:
                    os.remove(folder+'/'+key+addon+'_'+trial['name']+str(itrial)+'.pkl')

                #print("scores",key,":",scores)

    #return scores # to fix: is returning only the last scores!




###############################################
# ADDITIONAL FUNCTIONS

# objective function to be fit using point data from spiketrains (5th order polynomial)
def objective(x, a, b, c, d, e, f):
	return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f


from itertools import islice
def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


# Finds baseline in the firing rate
# Params:
#   l for smoothness (λ)
#   p for asymmetry
# Both have to be tuned to the data at hand.
# We found that generally is a good choice (for a signal with positive peaks):
#   10^2 ≤ l ≤ 10^9
#   0.001 ≤ p ≤ 0.1
# but exceptions may occur.
# In any case one should vary l on a grid that is approximately linear for log l
def baseline(y, l, p, niter=10):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = l * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w) # Do not create a new matrix, just update diagonal values
        Z = W + D
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z


def select_spikelist( spiketrains, edge=None, limits=None ):
    new_spiketrains = []
    for i,st in enumerate(spiketrains):

        # reject cells outside limits
        if limits:
            # use st entries as indexes to put ones
            x = int(i % edge)
            y = int(i / edge)
            # print(i, x,y)
            # if (x>10 and x<50) and (y>10 and y<54):
            if (x>limits[0][0] and x<limits[0][1]) and (y>limits[1][0] and y<limits[1][1]):
                # print('taken')
                new_spiketrains.append( st )
        else:
            new_spiketrains.append( st )

    return new_spiketrains


from scipy.signal import hilbert, butter, filtfilt
"""
from https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
Regrading low and high band pass limits:
For digital filters, Wn are in the same units as fs.
By default, fs is 2 half-cycles/sample, so these are normalized from 0 to 1, where 1 is the Nyquist frequency.
(Wn is thus in half-cycles / sample.)
"""
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs # fs is 10 samples per ms (dt=0.1) therefore fs = 10 kHz; its Nyquist freq = fs/2
    low = lowcut / nyq # 5000 / 5000 = 1 (lowest period)
    high = highcut / nyq # 200 / 5000 = 0.04 (highest period)
    b, a = butter(order, [low, high], btype='band')
    return b, a
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y
# convolution window
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def phase_coherence(s1, s2, folder, addon, lowcut=5000., highcut=200., fs=10000, order=5):
    """
    from:
    http://jinhyuncheong.com/jekyll/update/2017/12/10/Timeseries_synchrony_tutorial_and_simulations.html
    https://towardsdatascience.com/four-ways-to-quantify-synchrony-between-time-series-data-b99136c4a9c9

    Phase coherence is here the instantaneous phase synchrony measuring the phase similarities between two signals at each timepoint.
    The phase is to the angle of the signal when it is mapped between -pi to pi degrees.
    When two signals line up in phase their angular difference becomes zero.
    The angles are calculated through the hilbert transform of the signal.
    Phase coherence can be quantified by subtracting the angular difference from 1.
    """
    s1 = butter_bandpass_filter(s1,lowcut=lowcut,highcut=highcut,fs=fs,order=order)
    s2 = butter_bandpass_filter(s2,lowcut=lowcut,highcut=highcut,fs=fs,order=order)
    # s1 = butter_bandpass_filter(s1+1.,lowcut=lowcut,highcut=highcut,fs=fs,order=order)
    # s2 = butter_bandpass_filter(s2+1.,lowcut=lowcut,highcut=highcut,fs=fs,order=order)
    # s1 = butter_bandpass_filter(s1[100:]/np.max(s1[100:]),lowcut=lowcut,highcut=highcut,fs=fs,order=order)
    # s2 = butter_bandpass_filter(s2[100:]/np.max(s2[100:]),lowcut=lowcut,highcut=highcut,fs=fs,order=order)

    angle1 = np.angle(hilbert(s1),deg=False)
    angle2 = np.angle(hilbert(s2),deg=False)
    phase_coherence = 1-np.sin(np.abs(angle1-angle2)/2)
    N = len(angle1)

    # Plot results
    fig,ax = plt.subplots(3,1,figsize=(14,7),sharex=True)
    ax[0].plt(s1,color='r',label='fr1')
    ax[0].plt(s2,color='b',label='fr2')
    ax[0].legend(bbox_to_anchor=(0., 1.02, 1., .102),ncol=2)
    ax[0].set(xlim=[0,N], title='Band-passed firing rate')
    ax[1].plt(angle1,color='r')
    ax[1].plt(angle2,color='b')
    ax[1].set(ylabel='Angle',title='Angle at each Timepoint',xlim=[0,N])
    ax[2].plt(phase_coherence)
    ax[2].set(ylim=[0,1.1],xlim=[0,N],title='Instantaneous Phase Coherence',xlabel='Time',ylabel='Coherence')
    plt.tight_layout()
    fig.savefig(folder+'/PhaseCoherence_'+addon+'.svg', transparent=True)
    plt.close()
    fig.clear()


def LFP(data):
    v = data.filter(name="v")[0]
    g = data.filter(name="gsyn_exc")[0]
    # g = data.filter(name="gsyn_inh")[0]
    # We produce the current for each cell for this time interval, with the Ohm law:
    # I = g(V-E), where E is the equilibrium for exc, which usually is 0.0 (we can change it)
    # (and we also have to consider inhibitory condictances)
    iex = g*(v) #AMPA
    # the LFP is the result of cells' currents, with the Coulomb law:
    avg_i_by_t = np.sum(i,axis=1)/i.shape[0] # no distance involved for the moment
    sigma = 0.1 # [0.1, 0.01] # Dobiszewski_et_al2012.pdf
    lfp = (1/(4*np.pi*sigma)) * avg_i_by_t
    return lfp



def adaptation_index(data):
    # from NaudMarcilleClopathGerstner2008
    k = 2
    st = data.spiketrains
    if st == []:
        return None
    # ISI
    isi = np.diff(st)
    running_sum = 0
    for i,interval in enumerate(isi):
        if i < k:
            continue
        print(i, interval)
        running_sum = running_sum + ( (isi[i]-isi[i-1]) / (isi[i]+isi[i-1]) )
    return running_sum / len(isi)-k-1



def CC( duration, spiketrains, bin_size=10, auto=False ):
    """
    Binned-time Cross-correlation
    """
    print("CC")
    if spiketrains == [] :
        return NaN
    # create bin edges based on number of times and bin size
    # binning absolute time, and counting the number of spike times in each bin
    bin_edges = np.arange( 0, duration, bin_size )
    # print("bin_edges", bin_edges.shape, bin_edges)
    counts = []
    for spike_times in spiketrains:
        counts.append( np.histogram( spike_times, bin_edges )[0] )# spike count of bin_size bins
    counts = np.array(counts)
    CC = np.corrcoef(counts)
    # print(CC.shape, CC)
    return np.nanmean(CC)


def aCC( duration, spiketrains, bin_size=10, auto=False ):
    """
    Binned-time Cross-correlation
    """
    print("use CC")
    # if spiketrains == [] :
    #     return NaN
    # # create bin edges based on number of times and bin size
    # # binning absolute time, and counting the number of spike times in each bin
    # bin_edges = np.arange( 0, duration, bin_size )
    # for i, spike_times_i in enumerate(spiketrains):
    #     # print(i)
    #     for j, spike_times_j in enumerate(spiketrains):
    #         if auto and i!=j:
    #             continue
    #         itrain = np.histogram( spike_times_i, bin_edges )[0] # spike count of bin_size bins
    #         jtrain = np.histogram( spike_times_j, bin_edges )[0]

    #         if auto:
    #             CC.append( np.correlate(itrain, jtrain, mode="full") )
    #         else:
    #             CC.append( np.corrcoef(itrain, jtrain)[0,1] )
    # if auto:
    #     # return np.sum(CC, axis=0)
    #     return CC
    # else:
    #     return np.nanmean(CC)



def firingrate( start, end, spiketrains, bin_size=10 ):
    """
    Population rate
    as in https://neuronaldynamics.epfl.ch/online/Ch7.S2.html
    """
    if spiketrains == [] :
        return NaN
    # create bin edges based on start and end of slices and bin size
    bin_edges = np.arange( start, end, bin_size )
    print("coucou")
    # print("bin_edges",bin_edges.shape)
    # binning total time, and counting the number of spike times in each bin
    hist = np.zeros( bin_edges.shape[0]-1 )
    for spike_times in spiketrains:
        hist = hist + np.histogram( spike_times, bin_edges )[0]
    return ((hist / len(spiketrains) ) / bin_size ) * 1000 # average over population; result in ms *1000 to have it in sp/s



def isi( spiketrains, bins ):
    """
    Mean Inter-Spike Intervals for all spiketrains
    """
    isih = np.zeros(bins)
    for st in spiketrains:
        # print("st diff (int)", np.diff( st.magnitude ).astype(int) )
        isih = isih + np.histogram( np.diff( st.magnitude ).astype(int), len(isih) )[0]
    return isih



def cv( spiketrains, bins ):
    """
    Coefficient of variation
    """
    ii = isi(spiketrains, bins)
    return np.std(ii) / np.mean(ii)



# # ----------------------------------
# # -----      Phase Space      ------
# # ----------------------------------
def _v(x,y,p,I):
    gL = p['cm'] / p['tau_m'] # tau_m = cm / gL
    return ( -gL*(x-p['v_rest']) + gL*p['delta_T']*np.exp((x-p['v_thresh'])/p['delta_T']) -y +I ) / p['cm'] # Brette Gerstner 2005

def v_nullcline(x,p,I):
    gL = p['cm'] / p['tau_m']
    return -gL*(x-p['v_rest']) + gL*p['delta_T']*np.exp((x-p['v_thresh'])/p['delta_T']) + I # Touboul Brette 2008

def _w(x,y,p):
    return ( p['a']*(x-p['v_rest']) ) / p['tau_w'] # Brette Gerstner 2005

def w_nullcline(x,p,I):
    return p['a']*(x-p['v_rest']) # Touboul Brette 2008

def nullcline(f, params, I, limits, steps):
    fn = []
    c = np.linspace( limits[0], limits[1], steps ) # linearly spaced numbers
    for i in c:
        fn.append( f(i,params,I) )
    return c, fn

def Euler( f1, f2, params, Input, iv, dt, time):
    x = np.zeros(time)
    y = np.zeros(time)
    # initial values:
    x[0] = iv[0]
    y[0] = iv[1]
    I = 0 # init
    # compute and fill lists
    i=1
    while i < time:
        # integrating
        x[i] = x[i-1] + ( f1(x[i-1],y[i-1],params,I) )*dt
        y[i] = y[i-1] + ( f2(x[i-1],y[i-1],params) )*dt
        # discontinuity
        if x[i] >= params['v_spike']:
            x[i-1] = params['v_spike']
            x[i] = params['v_spike']
            y[i] = y[i] + params['b']
        i = i+1 # iterate
    return x, y

# # ----------------------------------
# # ----------------------------------

def H(params, limits, steps):
    p = params['Populations']['cell']['cellparams']
    # gL = p['cm'] / p['tau_m']
    cm = p['cm'] * pq.nF
    tau_m = p['tau_m']*pq.mS
    gL = cm / tau_m
    gc = gL/cm
    delta_T = p['delta_T']*pq.mV
    v_thresh = p['v_thresh']*pq.mV
    tau_w = p['tau_w']*pq.mS

    fn = []
    Vm = np.linspace( limits[0], limits[1], steps ) # linearly spaced numbers
    for v in Vm:
        v = v*pq.mV
        h = gc * (np.exp((v-v_thresh)/delta_T) -1.) - (1/tau_w)
        fn.append( h )
    return Vm, fn
