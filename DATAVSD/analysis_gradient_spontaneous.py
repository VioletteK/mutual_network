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


# folder = '_Nostim_ trials with only spontaneous activity, second set of 10 trials/'
folder = '_Nostim_ trials with only spontaneous activity, first 10 trials/'
# for i in range(11,21):
# for i in range(1,11):
for i in [2]:
    file = 'Nostim_'+str(i)
    L=[[[0 for k in range(100)] for l in range(100)]for i in range(511)]
    with open(folder+file, 'r') as fr:
        lines = fr.readlines()
        for t in range(511):
            for i in range(100):
                for j in range(100):
                    L[t][-i+99][j]=float(lines[t*10000+i*100+j])


    dt = 1
    run_time = 511
    size = 100
    injection_start,injection_end = 366,511
    interval = 20
    x = 0
    y = 0
    window = 100
    list_coord = [(x+i)*size+(y+j) for j in range(window) for i in range(window)]
    injection_points=[8555]
    min_point, max_point = 85,55
    ref_neurone = [85,55]

    # RdBu = cm.get_cmap('RdBu', 256)
    header='/home/margauxvrech/mutual_network/DATAVSD/'
    newheader=folder+file +'data'
    if not os.path.exists(header+newheader):
        os.makedirs(header+newheader)
    print('\nAnalysing data of '+file)
    print('\nData will be saved in '+newheader)
    Time_delay = np.arange(-20,30,1)
    V=len(L)
    Recorded_cell = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
    Recorded_cell_brut = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
    Norm_vect = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]


    Kernel1 = 1/16*np.array([
    [1,2,1],
    [2,4,2],
    [1,2,1]
    ])
    Kernel2=1/256*np.array([
    [1,4,6,4,1],
    [4,16,24,16,4],
    [6,24,36,24,6],
    [4,16,24,16,4],
    [1,4,6,4,1]
    ])
    X, Y = np.meshgrid([i for i in range(window)],[i for i in range(window)])

    def accentuation1(x):
        return 0.008*(np.tanh(300*x-2)+1)

    def accentuation2(x):
        return 0.008*(np.tanh((300*x-2))+1)
    # filtered_L=[]
    # for t in range(511):
    #     filtered_L.append(signal.convolve2d(accentuation1(np.array(L[t])),Kernel1, mode='same'))
    # vm=filtered_L

    # vm_base=[vm[t][ref_neurone[0]][ref_neurone[1]] for t in range(injection_start,injection_start+interval)]
    vm_base_brut=[L[t][ref_neurone[0]][ref_neurone[0]] for t in range(injection_start,injection_start+interval)]
    # c_X,xedges = np.histogram(vm_base,200,range=(-0.1,0.02))
    c_X_brut,xedges = np.histogram(vm_base_brut,200,range=(-0.1,0.02))


    VM_moyen=[np.mean([L[injection_start+t][i] for i in range(100)]) for t in Time_delay]
    fig=plt.figure()
    plt.plot(Time_delay,VM_moyen)
    filename= '/VM_moyen'+'.png'
    fig.savefig(header+newheader+filename, transparent=True)
    plt.close()
    fig.clf()

    print('Plotting MI ...')
    for time_delay in Time_delay:
        for i in range(window):
            for j in range(window):
                # vm_neurone = [vm[t][i][j] for t in range(injection_start+time_delay,injection_start+interval+time_delay)]
                vm_neurone_brut = [L[t][i][j] for t in range(injection_start+time_delay,injection_start+interval+time_delay)]
                # c_Y,xedges = np.histogram(vm_neurone,200,range=(-0.1,0.02))#,density=True)
                c_Y_brut,xedges = np.histogram(vm_neurone_brut,200,range=(-0.1,0.02))
                # Recorded_cell[injection_start+time_delay][i][j]= mutual_info_score(c_X,c_Y)
                Recorded_cell_brut[injection_start+time_delay][i][j]= mutual_info_score(c_X_brut,c_Y_brut)


        fig=plt.figure()
        # fig.add_subplot(2,1,1)


        plt.imshow(Recorded_cell_brut[injection_start+time_delay], cmap = 'inferno',interpolation='none')
        plt.clim([0,0.3])
        plt.colorbar()
        plt.title('MI brut '+str(injection_start+time_delay))

        # fig.add_subplot(2,1,2)
        # plt.imshow(Recorded_cell[injection_start+time_delay], cmap = 'inferno',interpolation='none')
        # plt.clim([0,0.028])
        # plt.colorbar()
        # plt.title('MI filtré '+str(injection_start+time_delay))
        [U,V]=np.gradient(Recorded_cell_brut[injection_start+time_delay])
        # plt.quiver(X,Y,-V,U,color='white')

        plt.show(block=True)
        filename= '/Gradient_'+str(injection_start+time_delay)+'.png'
        fig.savefig(header+newheader+filename, transparent=True)

        plt.close()
        fig.clf()

        for i in range(window):
            for j in range(window):

                Norm_vect[injection_start+time_delay][i][j]=np.linalg.norm([U[i][j],V[i][j]])
        fig=plt.figure()
        plt.imshow(Norm_vect[injection_start+time_delay], interpolation = 'none',vmin=0,vmax= 0.1,cmap= 'YlGnBu')
        plt.colorbar()
        filename= '/Norm_vect_'+str(injection_start+time_delay)+'.png'
        fig.savefig(header+newheader+filename)
        plt.close()
        fig.clf()

    fig = plt.figure()
    plt.imshow([[np.mean([Norm_vect[injection_start+time_delay][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)], interpolation = 'none',vmin = 0,vmax=0.05, cmap = 'viridis')
    plt.colorbar()
    filename= '/Moyenne_generale_'+str(injection_start+time_delay)+'.png'
    fig.savefig(header+newheader+filename)
    plt.close()
    fig.clf()


    def weighted_mean(values,weights):
        if len(values) != len(weights):
            return 'Error, must be same length'
        else :
            return np.sum([values[i]*weights[i] for i in range(len(values))])/np.sum(weights)


    fig = plt.figure()
    plt.imshow([[weighted_mean([Norm_vect[injection_start+time_delay][i][j] for time_delay in Time_delay], np.linspace(0,10,len(Time_delay))) for j in range(window)] for i in range(window)], interpolation = 'none',vmin = 0,vmax=0.05, cmap = 'viridis')
    plt.colorbar()
    filename= '/Weighted_mean_'+str(injection_start+time_delay)+'.png'
    fig.savefig(header+newheader+filename)
    plt.close()
    fig.clf()
