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
import matplotlib.pyplot as plot
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


def shan_entropy(c):
    #c_normalized = c
    c_normalized = c / float(np.sum(c))
    c_normalized = c_normalized[np.nonzero(c_normalized)]
    H = -sum(c_normalized* np.log2(c_normalized))
    return H


# ------------------------------------------------------------------------------


def analyse(params, folder, addon='', removeDataFile=False):
    print("\nAnalysing data...")
    print("----------------coucou-------------")
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
                scores = []
                scores.append(0)   # 0. Spike count
                scores.append(0.0) # 1. Inter-Spike Interval
                scores.append(0.0) # 2. Coefficient of Variation
                scores.append("")  # 3. additional label for adapting ISI
                scores.append(0.0) # 4. Firing rate
                scores.append(0.0) # 5. mean power in the Delta range (0.1-4 Hz)
                scores.append(0.0) # 6. mean power in the Spindles range (7-14 Hz)
                scores.append(0.0) # 7. Cross-Correlation


                ###################################

                if 'Vm' in params['Analysis'] and params['Analysis']['Vm'] and key in params['Analysis']['Vm']:
                    print('Vm')
                    vm=vm.T
                    print(vm.shape)
                    for i in range(len(vm)):#[3601]:#""[480,472]:#
                        fig = plot.figure()
                        #plot.plot(vm[i],color="red",linewidth=1)
                        plot.plot(vm[i],linewidth=1)
                        plot.ylim([-100,0.0]) # GENERIC
                        plot.ylabel('Membrane Potential (mV)')
                        plot.xlabel('Time (dt='+str(params['dt'])+' ms)')
                        #plot.xlim([7000,17000])#ZZOOOM
                        print(i)
                        fig.savefig(folder+'/vm_cell_'+key+str(i)+'__0.4.png', transparent=True)
                        # elif compteurbis==2:
                        #     print('coucou2')
                        #     fig.savefig(folder+'/vm_cell_'+key+str(i)+'__0.25''.png', transparent=True)
                        # else :
                        #     print('coucou3')
                        #     fig.savefig(folder+'/vm_cell_'+key+str(i)+'__0.4''.png', transparent=True)
                        plot.close()
                        fig.clf()

                    # Histo=[0 for i in range(30)]

                    ###### Tracer le Vm pour chaque cellule ########


                        ###### Tracer l'histogramme qui correspond à la moyenne des vm sur le temps ########

                        # h,xedges=np.histogram(vm[i],bins=30,range=(-70.,-40.))
                        # for yh in range(len(Histo)):
                        #     Histo[yh]+=h[yh]
                        # X,Y=np.meshgrid(xedges,yedges)
                        # plot.pcolormesh(X,Y,h)



                    ########## Commande pour tracer le vm #####################
                                # X=vm[60,15000:45000]
                                # Y=vm[48,15000:45000]
                                # bins=200
                                # interval=10
                                # MI=[]
                                #
                                # for dt in range(0,30000,interval):
                                #     Xdt=X[dt:dt+interval]
                                #     Ydt=Y[dt:dt+interval]
                                #     c_XY,xedges,yedges = np.histogram2d(Xdt,Ydt,bins,range=[[-90.,-40.],[-90.,-40.]])#,normed=True)
                                #     c_X,xedges = np.histogram(Xdt,bins,range=(-90.,-40.))#,density=True)#[8:]
                                #     c_Y,xedges = np.histogram(Ydt,bins,range=(-90.,-40.))#,density=True)
                                #     #print(c_XY,c_X,c_Y)
                                #     #print(xedges,yedges)
                                #
                                #
                                #     H_X = shan_entropy(c_X)
                                #     H_Y = shan_entropy(c_Y)
                                #     H_XY = shan_entropy(c_XY)
                                #     #print(H_X,H_Y,H_XY)
                                #     MIdt = H_X + H_Y - H_XY
                                #
                                #
                                #     MI.append(MIdt)
                                #
                                #
                                #
                                # MI_filtered=savgol_filter(MI, 51, 3)
                                # fig=plot.figure()
                                # plot.plot(MI)
                                # plot.plot(MI_filtered)
                                # plot.tight_layout()
                                # fig.savefig(folder+'/MI_test_fig_withoutdensity_2'+addon+'.png', transparent=True,dpi=300)
                                # plot.close()
                                # fig.clf()
                                #





                    # for yh in range(len(Histo)):
                    #     Histo[yh]=Histo[yh]/(len(vm)*60000)
                    #
                    # fig=plot.figure()
                    # plot.plot(Histo[8:])
                    # plot.xticks(np.arange(len(xedges[8:])),xedges[8:],rotation=90)
                    # plot.tight_layout()
                    # fig.savefig(folder+'/vm_cell_hist.png', transparent=True)
                    # plot.close()
                    # fig.clf()
                    # 0/0


                    # Mutual Info
                    # MIneighbors = np.mutual(vm[60], vm[61])
                    # MIforein = np.mutual(vm[60], vm[101])
                    # fig = plot.figure()
                    # plot.plot(MIneighbors,linewidth=1)
                    # plot.plot(MIforein,linewidth=1)
                    # plot.ylim([-100,0.0]) # GENERIC
                    # plot.ylabel('Mut')
                    # plot.xlabel('Time (dt='+str(params['dt'])+' ms)')
                    # fig.savefig(folder+'/MI_'+key+str(i)+'.png', transparent=True)
                    # plot.close()
                    # fig.clf()

                    # plot.plot(vm,linewidth=2)
                    # plot.ylim([-100,0.0]) # GENERIC
                    #################
                    # plot.xlim([9900,12500]) # E/IPSP single pulse (0.1 Hz)
                    #################
                    # plot.xlim([10000,60000]) # E/IPSP single pulse (0.1 Hz)
                    # plot.ylim([-66,-54]) # Control EPSP on RE
                    #################
                    # plot.ylim([-75.5,-71.5]) # Control EPSP on RS
                    # plot.ylim([-79.,-74.5]) # ACh EPSP on RS
                    # plot.ylim([-79.,-74.5]) # Control IPSP on RS
                    # plot.ylim([-79.,-74.5]) # ACh IPSP on RS
                    #################
                    # plot.ylim([-64.5,-60]) # Control EPSP on FS
                    # plot.ylim([-51.5,-47.]) # ACh EPSP on FS
                    # plot.ylim([-67.5,-63.5]) # Control IPSP on FS
                    # plot.ylim([-54.5,-50.5]) # ACh IPSP on FS
                    #################
                    # all Vms
                    # fig.savefig(folder+'/vm_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    # plot.ylabel('Membrane Potential (mV)')
                    # plot.xlabel('Time (dt='+str(params['dt'])+' ms)')
                    # plot.close()
                    # fig.clf()
                    # #################
                    # # Average Vms
                    # fig = plot.figure()
                    # plot.plot(np.mean(vm,axis=1),linewidth=2)
                    # plot.ylim([-100,0.0]) # GENERIC
                    # fig.savefig(folder+'/avg_vm_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    # plot.close()
                    # fig.clf()
                    # #################
                    # # Vm histogram
                    # fig = plot.figure()
                    # ylabel = key
                    # n,bins,patches = plot.hist(np.mean(vm,1), bins=50, normed=True) # 50*dt = 5ms bin
                    # fig.savefig(folder+'/Vm_histogram_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    # plot.close()
                    # fig.clear()



                ###################################
                if 'Rasterplot' in params['Analysis'] and key in params['Analysis']['Rasterplot']:
                    print('Rasterplot')

                    # cell selection
                    if key in params['Analysis']['Rasterplot']:
                        if params['Analysis']['Rasterplot'][key]['limits'] == 'all':
                            spikelist = select_spikelist( spiketrains=spiketrains, edge=edge )
                        else:
                            spikelist = select_spikelist( spiketrains=spiketrains, edge=edge, limits=params['Analysis']['Rasterplot'][key]['limits'] )

                    local_addon = addon
                    if params['Analysis']['Rasterplot']['interval']:
                        local_addon = local_addon +'_zoom'

                    fig = plot.figure()
                    for row,st in enumerate(spikelist):
                        train = st
                        if params['Analysis']['Rasterplot']['interval']:
                            train = train[ train>params['Analysis']['Rasterplot']['interval'][0] ]
                            train = train[ train<params['Analysis']['Rasterplot']['interval'][1] ]
                        plot.scatter( train, [row]*len(train), marker='o', edgecolors='none', s=0.2, c=params['Analysis']['Rasterplot'][key]['color'] )
                    #fig.savefig(folder+'/spikes_'+key+local_addon+'_'+trial['name']+str(itrial)+params['Analysis']['Rasterplot']['type'], transparent=True, dpi=params['Analysis']['Rasterplot']['dpi'])
                    plot.close()
                    fig.clf()



                ###################################
                if 'ISI' in params['Analysis'] and key in params['Analysis']['ISI'] and 'spikes' in rec:
                    print('ISI')
                    # Spike Count
                    if hasattr(spiketrains[0], "__len__"):
                        scores[0] = len(spiketrains[0])

                    # cell selection
                    spikelist = select_spikelist( spiketrains=spiketrains, edge=edge, limits=params['Analysis']['ISI'][key]['limits'] )

                    # ISI
                    print("time:", params['run_time'], "bins:", params['run_time']/params['Analysis']['ISI'][key]['bin'] )
                    isitot = isi(spikelist, int(params['run_time']/params['Analysis']['ISI'][key]['bin']) )
                    # print("ISI", isitot, isitot.shape)
                    if isinstance(isitot, (np.ndarray)):
                        if len(isitot)>1:
                            scores[1] = 0.0 #
                            scores[2] = 0.0 #

                            # ISI histogram
                            fig = plot.figure()
                            plot.semilogy(range(len(isitot)), isitot)
                            # plot.plot(range(len(isitot)), isitot)
                            plot.title("mean:"+str(scores[1])+" CV:"+str(scores[2]))
                            plot.ylabel('count')
                            plot.xlabel('ISI (bin=50ms)')
                            fig.savefig(folder+'/ISI_histogram_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                            plot.close()
                            fig.clear()

                            if 'ISI#' in params['Analysis'] and params['Analysis']['ISI#']:
                                # if strictly increasing, then spiking is adapting
                                # but check that there are no spikes beyond stimulation
                                # if spiketrains[0][-1] < params['Injections']['cell']['start'][-1] and all(x<y for x, y in zip(isitot, isitot[1:])):
                                if spiketrains[0][-1] < params['run_time']:
                                    if all(x<y for x, y in zip(isitot, isitot[1:])):
                                        scores[3] = 'adapting'
                                        # ISIs plotted against spike interval position
                                        fig = plot.figure()
                                        plot.plot(isitot,linewidth=2)
                                        plot.title("CV:"+str(scores[2])+" "+str(addon))
                                        # plot.xlim([0,10])
                                        plot.ylim([2,12.])
                                        plot.xlabel('Spike Interval #')
                                        plot.ylabel('ISI (ms)')
                                        fig.savefig(folder+'/ISI_interval_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                                        plot.close()
                                        fig.clf()



                ###################################
                # if 'FiringRate' in params['Analysis'] and key in params['Analysis']['FiringRate'] and 'spikes' in rec:
                #     print('FiringRate')
                #
                #     # spikelist = select_spikelist( spiketrains=spiketrains, edge=edge, limits=params['Analysis']['FiringRate'][key]['limits'] )
                #
                #     # for row,st in enumerate(spikelist):
                #     #     train = st
                #     #     if params['Analysis']['FiringRate']['interval']:
                #     #         train = train[ train>params['Analysis']['FiringRate']['interval'][0] ]
                #     #         train = train[ train<params['Analysis']['FiringRate']['interval'][1] ]
                #
                #     # firing rate
                #     #fr = firingrate(timeslice_start, timeslice_end, spiketrains, bin_size=10) # ms
                #     bin_size=10
                #     bin_edges = np.arange(timeslice_start, timeslice_end, bin_size )
                #     liste_histo=[]
                #     for sidx,spike_times in enumerate(spiketrains):
                #
                #         """
                #         Population rate
                #         as in https://neuronaldynamics.epfl.ch/online/Ch7.S2.html
                #         """
                #         # create bin edges based on start and end of slices and bin size
                #
                #         # print("bin_edges",bin_edges.shape)
                #         # binning total time, and counting the number of spike times in each bin
                #         #hist = np.zeros( bin_edges.shape[0]-1 )
                #
                #         hist = np.histogram( spike_times, bin_edges )[0]
                #         liste_histo.append(hist)
                #             #return ((hist / len(spiketrains) ) / bin_size ) * 1000 # average over population; result in ms *1000 to have it in sp/s
                #         #fr = firingrate(params, spikelist, bin_size=10) # ms
                #         fr=hist
                #         scores[4] = fr.mean()
                #         cvtot = cv(spiketrains, int(params['run_time']/params['Analysis']['FiringRate']['bin']) )
                #         fig = plot.figure()
                #         plot.plot(fr,linewidth=0.5)
                #         #plot.title("mean rate: %.2f$\pm$%.2f sp/s " )
                #         plot.ylim(params['Analysis']['FiringRate'][key]['firing'])
                #         fig.savefig(folder+'/firingrate_'+key+str(sidx)+'.png', transparent=True)
                #         # plot.xlim([3000.,4000.])
                #         # fig.savefig(folder+'/zoom_firingrate_'+key+addon+'.svg')
                #         plot.close()
                #         fig.clf()
                #     # liste_histo=np.array(liste_histo)
                #     # fig = plot.figure()
                #     # plot.pcolormesh(liste_histo)
                #     # plot.colorbar()
                #     # fig.savefig(folder+'/liste_histo'+key+addon+'.png', transparent=True)
                #     # # fig.savefig(folder+'/Matrix_'+key+addon+'.svg', transparent=True)
                #     plot.close()
                #     fig.clear()
                #     fig.clf()

                    # ###################
                    # ## spectrum
                    # fig = plot.figure()
                    # Fs = 1 / params['dt']  # sampling frequency
                    # sr = Fs**2 # sample rate
                    # # Compute the power spectrum 'classical way', with 2sec temporal window and 1sec overlap
                    # freq, P = signal.welch(fr, sr, window='hamming')
                    # # plot different spectrum types:
                    # sp = plot.semilogx(freq, P, color = 'r')
                    # delta = sp[0].get_ydata()[1:12] # 0.1-4 Hz interval power values
                    # spindle = sp[0].get_ydata()[18:35] # 7-14 Hz interval power values
                    # scores[5] = delta.mean()
                    # scores[6] = spindle.mean()
                    # plot.xlabel('Frequency (Hz)')
                    # plot.ylabel('Power spectrum (µV**2)')
                    # fig.savefig(folder+'/FR_Spectrum_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    # plot.close()
                    # fig.clear()



                ###################################
                if 'CrossCorrelation' in params['Analysis'] and params['Analysis']['CrossCorrelation'] and key in params['Analysis']['CrossCorrelation'] and 'spikes' in rec:
                    print('CrossCorrelation')
                    # cross-correlation
                    scores[7] = aCC(params, spiketrains, bin_size=params['Analysis']['CrossCorrelation'][key]['bin_size'])
                    print(scores[7])
                    print("CC:", scores[7])



                ###################################
                if 'LFP' in params['Analysis'] and params['Analysis']['LFP'] and 'v' in rec and 'gsyn_exc' in rec:
                    print('LFP')
                    lfp = LFP(data)
                    fig = plot.figure()
                    plot.plot(lfp)
                    fig.savefig(folder+'/LFP_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    plot.close()
                    fig.clear()
                    # # spectrum
                    # fig = plot.figure()
                    # Fs = 1 / params['dt']  # sampling frequency
                    # plot.title("Spectrum")
                    # # plot.magnitude_spectrum(lfp, Fs=Fs, scale='dB', color='red')
                    # plot.psd(lfp, Fs=Fs)
                    # fig.savefig(folder+'/Spectrum_'+key+addon+'_'+str(trial_name)+str(trial)+'.png')
                    # fig.clear()



                ###################################
                if 'ConductanceBalance' in params['Analysis'] and params['Analysis']['ConductanceBalance'] and 'gsyn_exc' in rec and 'gsyn_inh' in rec:
                    if not key in params['Analysis']['ConductanceBalance']:
                        continue
                    if trial['name'] not in params['Analysis']['ConductanceBalance'][key]['trials']:
                        continue
                    print('Conductance Balance')

                    # Average conductances
                    avg_gexc = np.mean(gexc, axis=1)
                    avg_ginh = np.mean(ginh, axis=1)
                    # conductances
                    fig = plot.figure()
                    plot.plot(avg_gexc,linewidth=2, color='red')
                    plot.plot(avg_ginh,linewidth=2, color='blue')
                    plot.xlabel('Time (s)')
                    plot.ylabel('Conductance (µS)')
                    fig.savefig(folder+'/avg_conductance_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    plot.close()
                    fig.clf()
                    # Conductance balance — a measure of contrast between excitation and inhibition
                    avg_gbalance = avg_gexc / (avg_gexc+avg_ginh)
                    fig = plot.figure()
                    plot.plot(avg_gbalance,linewidth=2, color='black')
                    plot.title("Mean conductance ratio: %.2f" % (avg_gbalance.nanmean()) )
                    plot.xlabel('Time (s)')
                    plot.ylabel('Conductance ratio')
                    fig.savefig(folder+'/avg_conductanceratio_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    plot.close()
                    fig.clf()





                # for systems with low memory :)
                if removeDataFile:
                    os.remove(folder+'/'+key+addon+'_'+trial['name']+str(itrial)+'.pkl')

                print("scores",key,":",scores)

    return scores # to fix: is returning only the last scores!




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
    fig,ax = plot.subplots(3,1,figsize=(14,7),sharex=True)
    ax[0].plot(s1,color='r',label='fr1')
    ax[0].plot(s2,color='b',label='fr2')
    ax[0].legend(bbox_to_anchor=(0., 1.02, 1., .102),ncol=2)
    ax[0].set(xlim=[0,N], title='Band-passed firing rate')
    ax[1].plot(angle1,color='r')
    ax[1].plot(angle2,color='b')
    ax[1].set(ylabel='Angle',title='Angle at each Timepoint',xlim=[0,N])
    ax[2].plot(phase_coherence)
    ax[2].set(ylim=[0,1.1],xlim=[0,N],title='Instantaneous Phase Coherence',xlabel='Time',ylabel='Coherence')
    plot.tight_layout()
    fig.savefig(folder+'/PhaseCoherence_'+addon+'.svg', transparent=True)
    plot.close()
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
