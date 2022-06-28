"""
Copyright (c) 2016-2022, Domenico GUARINO
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Universite Paris Saclay nor the
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
import scipy
import scipy.linalg



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
                    if 'subsampling' in params['Analysis'] and params['Analysis']['subsampling']:
                        choice_indices = np.random.choice(len(data.spiketrains), params['Analysis']['subsampling'], replace=False)
                        for spiketrain in [data.spiketrains[i] for i in choice_indices]:
                            spiketrains.append(spiketrain[ (spiketrain>=timeslice_start) & (spiketrain<=timeslice_end) ])
                    else:
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
                # if 'scores' in params['Analysis'] and key in params['Analysis']['scores']:
                #     scores.append(0)   # 0. Spike count
                #     scores.append(0.0) # 1. Inter-Spike Interval
                #     scores.append(0.0) # 2. Coefficient of Variation
                #     scores.append("")  # 3. additional label for adapting ISI
                #     scores.append(0.0) # 4. Firing rate
                #     scores.append(0.0) # 5. mean power in the Delta range (0.1-4 Hz)
                #     scores.append(0.0) # 6. mean power in the Spindles range (7-14 Hz)
                #     scores.append(0.0) # 7. Cross-Correlation
                #     scores.append(0.0) # 8. Conductance Balance


                ###################################
                if 'Vm' in params['Analysis'] and params['Analysis']['Vm'] and key in params['Analysis']['Vm']:
                    print('Vm')
                    ##### COLORMAP
                    size = params['Recorders']['py']['v']['size']+1
                    dt = params['dt']
                    run_time = params['run_time']
                    interval = 10

                    for i in range(10000,int(run_time//dt),interval):
                    #for i in range(int(injection_start//dt),int((injection_end+100)//dt),interval):
                        #a colormap each interval*dt from the injection's beginning
                        Vm_i=np.zeros((size,size)) #our future colormap
                        for j in range(size**2):
                            Vm_i[j%size][j//size] = np.mean([vm[k][j] for k in range(i,min(i+interval,len(vm)))])
                            # Vm_i[j%size][j//size] = vm[i][j]
                        fig = plt.figure()
                        plt.imshow(Vm_i, cmap=matplotlib.cm.get_cmap('RdBu_r'),interpolation='none',vmin=-80,vmax=-50)
                        plt.colorbar()
                        tmin=float(i*dt)
                        tmax=float((i+interval)*dt)
                        plt.title('window ['+str(round(tmin,3))+':'+str(round(tmax,3))+'] ms')
                        fig.savefig(folder+'/tau='+str(params['Populations']['py']['cellparams']['tau_w'])+'Colormap'+str(i)+'.png', transparent=True)
                        plt.close()
                        fig.clf()
                    ##### Vm plot
                    # vm = vm.transpose()
                    # for iv,v in enumerate(vm):
                    # ###### tracer tous les vms en svg
                    #     fig = plt.figure()
                    #     plt.plot(v,linewidth=2)
                    #     plt.plot([i for i in range(len(v))],[-50 for i in range(len(v))],'.')
                    #     #plt.ylim([-100,0.0]) # GENERIC
                    #     #################
                    #     # plt.xlim([9900,12500]) # E/IPSP single pulse (0.1 Hz)
                    #     #################
                    #     # plt.xlim([10000,60000]) # E/IPSP single pulse (0.1 Hz)
                    #     # plt.ylim([-66,-54]) # Control EPSP on RE
                    #     #################
                    #     # plt.ylim([-75.5,-71.5]) # Control EPSP on RS
                    #     # plt.ylim([-79.,-74.5]) # ACh EPSP on RS
                    #     # plt.ylim([-79.,-74.5]) # Control IPSP on RS
                    #     # plt.ylim([-79.,-74.5]) # ACh IPSP on RS
                    #     #################
                    #     # plt.ylim([-64.5,-60]) # Control EPSP on FS
                    #     # plt.ylim([-51.5,-47.]) # ACh EPSP on FS
                    #     # plt.ylim([-67.5,-63.5]) # Control IPSP on FS
                    #     # plt.ylim([-54.5,-50.5]) # ACh IPSP on FS
                    #     #################
                    #     # each Vms
                    #     fig.savefig(folder+'/vm_'+'neurone\n'+str(iv)+'.svg', transparent=True)
                    #     plt.ylabel('Membrane Potential (mV)')
                    #     plt.xlabel('Time (dt='+str(params['dt'])+' ms)')
                    #     plt.ylim([-80,-40])
                    #     plt.close()
                    #     fig.clf()
                    #################
                    # # Average Vms
                    # fig = plt.figure()
                    # plt.plot(np.mean(vm,axis=1),linewidth=2)
                    # plt.ylim([-100,0.0]) # GENERIC
                    # fig.savefig(folder+'/avg_vm_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    # plt.close()
                    # fig.clf()
                    # #################
                    # # Vm histogram
                    # fig = plt.figure()
                    # ylabel = key
                    # n,bins,patches = plt.hist(np.mean(vm,1), bins=50, normed=True) # 50*dt = 5ms bin
                    # fig.savefig(folder+'/Vm_histogram_'+key+addon+'_'+trial['name']+str(itrial)+'.svg', transparent=True)
                    # plt.close()
                    # fig.clear()








                # for systems with low memory :)
                if removeDataFile:
                    os.remove(folder+'/'+key+addon+'_'+trial['name']+str(itrial)+'.pkl')

    #             print("scores",key,":",scores)
    #             # end of analysis level
    #         # end of trial level
    #     # end of population level
    # return scores
