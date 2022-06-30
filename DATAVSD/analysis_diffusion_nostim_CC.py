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
#for i in range(11,21):
for i in range(1,11):
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

    x = 0
    y = 0
    window = 100
    list_coord = [(x+i)*size+(y+j) for j in range(window) for i in range(window)]
    injection_points=[5050]
    min_point, max_point = 50,50


    ref_neurone = [85,55]


    Time_delay = np.arange(-20,18,1)
    number_of_annulus = 20
    interval = 20
    r = 20
    header='/home/margauxvrech/mutual_network/DATAVSD/'
    newheader=folder+file +'data'
    if not os.path.exists(header+newheader):
        os.makedirs(header+newheader)
    print('\nAnalysing data of '+file)
    print('\nData will be saved in '+newheader)

    V=len(L)
    Recorded_cell_brut = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]




    X, Y = np.meshgrid([i for i in range(window)],[i for i in range(window)])

    vm_base_brut=[L[t][ref_neurone[0]][ref_neurone[1]] for t in range(min(injection_start,injection_end-1),min(injection_start+interval,injection_end))]

    Time_maximum1 = []
    list_rho1 = []


    Annulus = [[] for i in range(number_of_annulus)]



    #this is the ray of the disk around our central cell
    print('... the width of the annulus is '+ str(r/number_of_annulus))
    for a in range(number_of_annulus):
        for i in range(100):
            for j in range(100):
                d = np.linalg.norm([i-ref_neurone[0],j-ref_neurone[1]])
                for a in range(number_of_annulus):
                    if a* r/number_of_annulus<d<= (a+1)*r/number_of_annulus:
                        Annulus[a].append([i,j])
                        break

    Annulus = list(filter(([]).__ne__, Annulus))



    #########################################################
    #ANNULUS ANALYSIS
    ######################################################





    print('Calculating CC...')




    fig=plt.figure()
    for A in Annulus :
        MI_annulus = []
        for time_delay in Time_delay :
            MI_delay = []
            for neuron in A :
                vm_neurone_brut = [L[t][neuron[0]][neuron[1]] for t in range(min(injection_start+time_delay,injection_end-1),min(injection_start+interval+time_delay,injection_end))]
                MI_delay.append(float(signal.correlate(vm_base_brut,vm_neurone_brut, mode = 'valid')))
            MI_annulus.append(np.mean(MI_delay, dtype=np.float64))



        if Annulus.index(A)>1:
            maxi1 = max(MI_annulus)
            Time_maximum1.append(Time_delay[list(MI_annulus).index(maxi1)])
            list_rho1.append(Annulus.index(A))
            plt.plot(Time_delay,MI_annulus, color = cm.hsv(Annulus.index(A)/len(Annulus)),marker='x')#,linestyle='-.')
            # plt.plot(Time_delay,MI_annulus_filtered,label='Anneau '+str(Annulus.index(A)), color = cm.hsv(Annulus.index(A)/len(Annulus)))#,linestyle='-.')
    plt.title('CC')
    plt.xlabel('Time delay')
    plt.ylabel("CC")
    plt.legend(frameon=False)
    fig.savefig(header+newheader+'/CC_'+str(len(Annulus))+'_annulus'+str(ref_neurone)+'.png')
    plt.close()
    fig.clf()
    print('\nPlotting the Annulus')
    #plot the Annulus
    Annulus_plot = np.zeros((window, window))
    for A in Annulus:
        for neuron in A :
            Annulus_plot[neuron[0]][neuron[1]]=Annulus.index(A)+2

    fig1 = plt.figure()
    plt.imshow(Annulus_plot, cmap = 'inferno',interpolation = 'none')
    plt.colorbar()
    plt.title('Annulus')
    fig1.savefig(header+newheader+'/Annulus_CC'+'.png')
    plt.close()
    fig1.clf()

    fig2 = plt.figure()

    plt.plot(list_rho1,Time_maximum1,'x')
    lr = scipy.stats.linregress(list_rho1,Time_maximum1)
    y=[lr[0]*i+lr[1] for i in list_rho1]
    plt.plot(list_rho1,y,c='r', label="pente="+str(lr[0])+", R2="+str(lr[2]**2))
    plt.show()
    plt.legend()
    plt.title('Maximum')

    fig2.savefig(header+newheader+'/Raw_Maximum_CC'+'.png')
    plt.close()
    fig2.clf()
