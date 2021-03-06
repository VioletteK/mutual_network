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
    #where the activity starts
    interval = 20
    x = 0
    y = 0
    window = 100
    list_coord = [(x+i)*size+(y+j) for j in range(window) for i in range(window)]
    #list of the cooridinate at their index in the vm list
    ref_neurone = [85,55]


    header='/home/margauxvrech/mutual_network/DATAVSD/'
    newheader=folder+file +'data'
    if not os.path.exists(header+newheader):
        os.makedirs(header+newheader)
    print('\nAnalysing data of '+file)
    print('\nData will be saved in '+newheader)


    Time_delay = np.arange(-20,35,1)
    V=len(L)

    Recorded_cell_brut = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
    #MI
    Norm_vect = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
    #Norm of the gradient vector



    X, Y = np.meshgrid([i for i in range(window)],[i for i in range(window)])

    def accentuation1(x):
        return 0.008*(np.tanh(300*x-2)+1)




    vm_base_brut=[L[t][ref_neurone[0]][ref_neurone[0]] for t in range(injection_start,injection_start+interval)]
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

                vm_neurone_brut = [L[t][i][j] for t in range(injection_start+time_delay,injection_start+interval+time_delay)]

                c_Y_brut,xedges = np.histogram(vm_neurone_brut,200,range=(-0.1,0.02))

                Recorded_cell_brut[injection_start+time_delay][i][j]= mutual_info_score(c_X_brut,c_Y_brut)

        #Plotting the MI
        # fig=plt.figure()
        # plt.imshow(Recorded_cell_brut[injection_start+time_delay], cmap = 'inferno',interpolation='none')
        # plt.clim([0,0.3])
        # plt.colorbar()
        # plt.title('MI brut '+str(injection_start+time_delay))


        [U,V]=np.gradient(Recorded_cell_brut[injection_start+time_delay])
        # filename= '/Gradient_'+str(injection_start+time_delay)+'.png'
        # fig.savefig(header+newheader+filename, transparent=True)
        #
        # plt.close()
        # fig.clf()

        for i in range(window):
            for j in range(window):

                Norm_vect[injection_start+time_delay][i][j]=np.linalg.norm([U[i][j],V[i][j]])

        #Plotting the norm of the gradient vectors
        # fig=plt.figure()
        # plt.imshow(Norm_vect[injection_start+time_delay], interpolation = 'none',vmin=0,vmax= 0.1,cmap= 'YlGnBu')
        # plt.colorbar()
        # filename= '/Norm_vect_'+str(injection_start+time_delay)+'.png'
        # fig.savefig(header+newheader+filename)
        # plt.close()
        # fig.clf()

    #Mean norm
    fig = plt.figure()
    plt.imshow([[np.mean([Norm_vect[injection_start+time_delay][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)], interpolation = 'none', vmax = 0.07, cmap = 'viridis')
    plt.colorbar()
    filename= '/Mean_vector_norm_MI'+'.png'
    fig.savefig(header+newheader+filename)
    plt.close()
    fig.clf()



    #Mean MI
    fig = plt.figure()
    plt.imshow([[np.mean([Recorded_cell_brut[injection_start+time_delay][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)], interpolation = 'none', cmap = 'viridis')
    plt.colorbar()
    filename= '/Mean_MI'+'.png'
    fig.savefig(header+newheader+filename)
    plt.close()
    fig.clf()
