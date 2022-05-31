"""
Copyright (c) 2016-2021, Domenico GUARINO
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

################################
import matplotlib
matplotlib.use('Agg') # to be used when DISPLAY is undefined
################################
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
from scipy import sparse
from scipy.sparse.linalg import spsolve
import scipy
import scipy.linalg
import scipy.signal as signal
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.optimize import curve_fit



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
    # transient
    if not 'transient' in params['Analysis']:
       params['Analysis']['transient'] = 1000 # ms, INIT in case of older param files


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
                if 'scores' in params['Analysis'] and key in params['Analysis']['scores']:
                    scores.append(0)   # 0. Spike count
                    scores.append(0.0) # 1. Inter-Spike Interval
                    scores.append(0.0) # 2. Coefficient of Variation
                    scores.append("")  # 3. additional label for adapting ISI
                    scores.append(0.0) # 4. Firing rate
                    scores.append(0.0) # 5. mean power in the Delta range (0.1-4 Hz)
                    scores.append(0.0) # 6. mean power in the Spindles range (7-14 Hz)
                    scores.append(0.0) # 7. Cross-Correlation
                    scores.append(0.0) # 8. Conductance Balance


                ###################################
                if 'Vm' in params['Analysis'] and params['Analysis']['Vm'] and key in params['Analysis']['Vm']:
                    print('Vm')
                    # print(vm)
                    fig = plt.figure()
                    plt.plot(vm,linewidth=2)
                    plt.ylim([-100,0.0]) # GENERIC
                    #################
                    # plt.xlim([9900,12500]) # E/IPSP single pulse (0.1 Hz)
                    #################
                    # plt.xlim([10000,60000]) # E/IPSP single pulse (0.1 Hz)
                    # plt.ylim([-66,-54]) # Control EPSP on RE
                    #################
                    # plt.ylim([-75.5,-71.5]) # Control EPSP on RS
                    # plt.ylim([-79.,-74.5]) # ACh EPSP on RS
                    # plt.ylim([-79.,-74.5]) # Control IPSP on RS
                    # plt.ylim([-79.,-74.5]) # ACh IPSP on RS
                    #################
                    # plt.ylim([-64.5,-60]) # Control EPSP on FS
                    # plt.ylim([-51.5,-47.]) # ACh EPSP on FS
                    # plt.ylim([-67.5,-63.5]) # Control IPSP on FS
                    # plt.ylim([-54.5,-50.5]) # ACh IPSP on FS
                    #################
                    # all Vms
                    fig.savefig(folder+'/vm_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    plt.ylabel('Membrane Potential (mV)')
                    plt.xlabel('Time (dt='+str(params['dt'])+' ms)')
                    plt.close()
                    fig.clf()
                    #################
                    # Average Vms
                    fig = plt.figure()
                    plt.plot(np.mean(vm,axis=1),linewidth=2)
                    plt.ylim([-100,0.0]) # GENERIC
                    fig.savefig(folder+'/avg_vm_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    plt.close()
                    fig.clf()
                    #################
                    # Vm histogram
                    fig = plt.figure()
                    ylabel = key
                    n,bins,patches = plt.hist(np.mean(vm,1), bins=50, normed=True) # 50*dt = 5ms bin
                    fig.savefig(folder+'/Vm_histogram_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    plt.close()
                    fig.clear()



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

                    fig = plt.figure()
                    for row,st in enumerate(spikelist):
                        train = st
                        if params['Analysis']['Rasterplot']['interval']:
                            train = train[ train>params['Analysis']['Rasterplot']['interval'][0] ]
                            train = train[ train<params['Analysis']['Rasterplot']['interval'][1] ]
                        plt.scatter( train, [row]*len(train), marker='o', edgecolors='none', s=0.02, c=params['Analysis']['Rasterplot'][key]['color'] )
                    fig.savefig(folder+'/spikes_'+key+local_addon+'_'+trial['name']+str(itrial)+params['Analysis']['Rasterplot']['type'], transparent=True, dpi=params['Analysis']['Rasterplot']['dpi'])
                    plt.close()
                    fig.clf()


                ###################################
                if 'FiringRate' in params['Analysis'] and key in params['Analysis']['FiringRate'] and 'spikes' in rec:
                    print('FiringRate')
                    fr = firingrate(timeslice_start, timeslice_end, spiketrains, bin_size=params['Analysis']['FiringRate']['bin']) # ms
                    fr_mean = fr.mean()
                    fr_cc = CC(timeslice_end-timeslice_start, spiketrains, bin_size=params['Analysis']['FiringRate']['bin'])
                    if 'scores' in params['Analysis'] and key in params['Analysis']['scores']:
                        scores[4] = fr_mean
                        scores[7] = fr_cc
                    cvtot = cv(spiketrains, int(params['run_time']/params['Analysis']['FiringRate']['bin']) )
                    fig = plt.figure()
                    plt.plot(fr,linewidth=0.5)
                    plt.title("mean rate: %.2f$\pm$%.2f sp/s - CV: %.2f - CC: %.3f" % (fr_mean, fr.std(), cvtot, fr_cc) )
                    plt.ylim(params['Analysis']['FiringRate'][key]['firing'])
                    fig.savefig(folder+'/firingrate_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    # plt.xlim([3000.,4000.])
                    # fig.savefig(folder+'/zoom_firingrate_'+key+addon+'.svg')
                    plt.close()
                    fig.clf()
                    ###################
                    ## spectrum
                    fig = plt.figure()
                    Fs = 1 / params['dt']  # sampling frequency
                    sr = Fs**2 # sample rate
                    # Compute the power spectrum 'classical way', with 2sec temporal window and 1sec overlap
                    freq, P = signal.welch(fr, sr, window='hamming')
                    # plot different spectrum types:
                    sp = plt.semilogx(freq, P, color = 'r')
                    delta = sp[0].get_ydata()[1:12] # 0.1-4 Hz interval power values
                    spindle = sp[0].get_ydata()[18:35] # 7-14 Hz interval power values
                    if 'scores' in params['Analysis'] and key in params['Analysis']['scores']:
                        scores[5] = delta.mean()
                        scores[6] = spindle.mean()
                    plt.xlabel('Frequency (Hz)')
                    plt.ylabel('Power spectrum (µV**2)')
                    fig.savefig(folder+'/FR_Spectrum_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    plt.close()
                    fig.clear()
                    ###################
                    if 'stats' in params['Analysis']['FiringRate'] and 'stats'==True:
                        print("spikes stats with bin =", params['Analysis']['FiringRate']['bin'])
                        # as for dF/F
                        sp_bins = int((timeslice_end-timeslice_start)/params['Analysis']['FiringRate']['bin']) # in ms
                        spikes_isih = []
                        spikes_rate = []
                        # for each cell
                        for st in spiketrains:
                            # per cell firingrate
                            spikes_isih.extend( np.histogram( np.diff( st.magnitude ).astype(int), sp_bins )[0] ) # ms, grouped by bin size
                            spikes_rate.extend( np.histogram( st, sp_bins )[0] )
                        # plot ISI
                        spikes_isihs = np.histogram(spikes_isih,4000)[0]
                        for dfri,dffrtis in enumerate(spikes_isihs): # correct for 0 hist
                            if dfri>0:
                                if dffrtis<1:
                                    spikes_isihs[dfri] = spikes_isihs[dfri-1]
                        x_isi = np.arange(len(spikes_isihs))
                        fig = plt.figure()
                        plt.plot( x_isi, spikes_isihs )
                        plt.ylabel('occurrences')
                        plt.yscale('log')
                        plt.xscale('log')
                        plt.xlim([1,4000])
                        plt.xlabel('spikes ISI')
                        fig.savefig(folder+'/spikes_isi_'+key+addon+'_'+trial['name']+str(itrial)+'.png', transparent=True)
                        fig.savefig(folder+'/spikes_isi_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                        plt.close()
                        fig.clear()
                        fig.clf()
                        # plot Rates
                        spikes_rates = np.histogram(spikes_rate,1000)[0]
                        for dfri,dffrtis in enumerate(spikes_rates): # correct for 0 hist
                            if dfri>0:
                                if dffrtis<1:
                                    spikes_rates[dfri] = spikes_rates[dfri-1]
                        x_rate = np.arange(len(spikes_rates))
                        fig = plt.figure()
                        plt.plot(x_rate, spikes_rates)
                        plt.ylabel('occurrences')
                        plt.yscale('log')
                        plt.xscale('log')
                        plt.xlim([1,50])
                        plt.xlabel('spikes/s')
                        fig.savefig(folder+'/spikes_rate_'+key+addon+'_'+trial['name']+str(itrial)+'.png', transparent=True)
                        fig.savefig(folder+'/spikes_rate_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                        plt.close()
                        fig.clear()
                        fig.clf()


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
                    fig = plt.figure()
                    plt.plot(avg_gexc,linewidth=2, color='red')
                    plt.plot(avg_ginh,linewidth=2, color='blue')
                    plt.xlabel('Time (s)')
                    plt.ylabel('Conductance (µS)')
                    fig.savefig(folder+'/avg_conductance_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    plt.close()
                    fig.clf()
                    # Conductance balance — a measure of contrast between excitation and inhibition
                    avg_gbalance = avg_gexc / (avg_gexc+avg_ginh)
                    if 'scores' in params['Analysis'] and key in params['Analysis']['scores']:
                        scores[8] = avg_gbalance.nanmean()
                    fig = plt.figure()
                    plt.plot(avg_gbalance,linewidth=2, color='black')
                    plt.title("Mean conductance ratio: %.2f" % (avg_gbalance.nanmean()) )
                    plt.xlabel('Time (s)')
                    plt.ylabel('Conductance ratio')
                    fig.savefig(folder+'/avg_conductanceratio_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    plt.close()
                    fig.clf()


                # for systems with low memory :)
                if removeDataFile:
                    os.remove(folder+'/'+key+addon+'_'+trial['name']+str(itrial)+'.pkl')

                print("scores",key,":",scores)
                # end of analysis level
            # end of trial level
        # end of population level
    return scores




###############################################
# ADDITIONAL FUNCTIONS

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


def firingrate( start, end, spiketrains, bin_size=10 ):
    """
    Population rate
    as in https://neuronaldynamics.epfl.ch/online/Ch7.S2.html
    """
    if spiketrains == [] :
        return NaN
    # create bin edges based on start and end of slices and bin size
    bin_edges = np.arange( start, end, bin_size )
    # print("bin_edges",bin_edges.shape)
    # binning total time, and counting the number of spike times in each bin
    hist = np.zeros( bin_edges.shape[0]-1 )
    for spike_times in spiketrains:
        hist += np.histogram( spike_times, bin_edges )[0]
    return ((hist / len(spiketrains) ) / bin_size ) * 1000 # average over population; result in ms *1000 to have it in sp/s


def isi( spiketrains, bins ):
    """
    Mean Inter-Spike Intervals for all spiketrains
    """
    isih = np.zeros(bins)
    for st in spiketrains:
        # print("st diff (int)", np.diff( st.magnitude ).astype(int) )
        isih += np.histogram( np.diff( st.magnitude ).astype(int), bins )[0]
    return isih


def cv( spiketrains, bins ):
    """
    Coefficient of variation
    """
    ii = isi(spiketrains, bins)
    return np.std(ii) / np.mean(ii)
