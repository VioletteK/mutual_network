import os
import numpy as np
import scipy
import scipy.signal as signal

################################
import matplotlib
matplotlib.use('Agg') # to be used when DISPLAY is undefined
################################
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.signal import savgol_filter
from sklearn.metrics import mutual_info_score

RdBu = cm.get_cmap('RdBu_r', 256)


# # folder = '_Nostim_ trials with only spontaneous activity, second set of 10 trials/'
# folder = '_Nostim_ trials with only spontaneous activity, first 10 trials/'
# #for i in range(11,21):
# for i in range(1,11):
#     file = 'Nostim_'+str(i)
#     L=[[[0 for k in range(100)] for l in range(100)]for i in range(511)]
#     with open(folder+file, 'r') as fr:
#         lines = fr.readlines()
#         for t in range(511):
#             for i in range(100):
#                 for j in range(100):
#                     L[t][-i+99][j]=float(lines[t*10000+i*100+j])
#folder = 'Mouses_3-5-6/20160912/'
#folder = 'Mouses_3-5-6/20160914/'
folder = 'Mouses_3-5-6/20160916/'
for i in range(1,21):
    file = 'C2_'+str(i)+'.txt'
    L=[[[0 for k in range(100)] for l in range(100)]for i in range(511)]
    with open(folder+file, 'r') as fr:
        lines = fr.readlines()
        for t in range(511):
            for i in range(100):
                    L[t][i]=[float(lines[t*100+j].split('\t')[i]) for j in range(100)]
    dt = 1
    run_time = 511
    size = 100


    injection_start,injection_end = 100,511
    x = 0
    y = 0
    window = 100
    list_coord = [(x+i)*size+(y+j) for j in range(window) for i in range(window)]
    injection_points=[5050]
    min_point, max_point = 50,50


    ref_neurone = [50,50]
    looked_neurone = [40,60]

    Time_delay = np.arange(-50,100,1)
    interval = 10

    header='/home/margauxvrech/mutual_network/DATAVSD/'
    newheader=folder+file +'data'
    if not os.path.exists(header+newheader):
        os.makedirs(header+newheader)
    print('\nAnalysing data of '+file)
    print('\nData will be saved in '+newheader)


    V = len(L)

    vm_base=[L[t][ref_neurone[0]][ref_neurone[1]] for t in range(injection_start,min(injection_start+interval,injection_end))]
    c_X,xedges = np.histogram(vm_base,500,range=(-0.1,0.02))
    MI_delay=[]
    for time_delay in Time_delay :
        vm_neurone = [L[t][looked_neurone[0]][looked_neurone[1]] for t in range(min(injection_start+time_delay,510),min(injection_start+time_delay+interval,510))]
        c_Y,xedges = np.histogram(vm_neurone,500,range=(-0.1,0.02))
        MI_delay.append(mutual_info_score(c_X,c_Y))

    fig=plt.figure()
    fig.add_subplot(2,1,1)
    plt.plot(Time_delay,MI_delay)
    plt.title('Mutual information des neurones '+str(looked_neurone)+' et '+str(ref_neurone))
    plt.xlabel('Time delay')
    plt.ylabel("MI")
    plt.legend()

    fig.add_subplot(2,1,2)
    v = [L[injection_start+t][looked_neurone[0]][looked_neurone[1]] for t in Time_delay]


    plt.plot(Time_delay,v,linewidth=2, label = 'looked_'+str(looked_neurone))
    v =np.array([L[injection_start+t][ref_neurone[0]][ref_neurone[1]] for t in Time_delay])
    plt.plot(Time_delay,v+0.01,linewidth=2,color = 'r', label= 'ref_'+str(ref_neurone))

    # plt.ylim([-90,0])
    plt.ylabel("Vm")
    plt.legend()

    fig.savefig(header+newheader+'/MI_'+str(looked_neurone)+str(ref_neurone)+'.png')
    #fig.savefig(folder+'/Mutual Information avec '+str(number_of_annulus)+' anneaux' +'.png')
    plt.close()
    fig.clf()
