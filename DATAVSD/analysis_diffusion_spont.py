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
# for i in range(1,11):
# for i in [1]:
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

    # injection_start,injection_end = 110,511
    injection_start,injection_end = 366,511
    # injection_start,injection_end = 393,511
    # injection_start,injection_end = 340,511
    # injection_start,injection_end = 168,511
    # injection_start,injection_end = 57,511
    # injection_start,injection_end =
    # injection_start,injection_end = 233,511
    # injection_start,injection_end = 400,511
    # injection_start,injection_end = 364,511
    x = 0
    y = 0
    window = 100
    list_coord = [(x+i)*size+(y+j) for j in range(window) for i in range(window)]
    injection_points=[5050]
    min_point, max_point = 50,50

    # ref_neurone = [80,50]
    ref_neurone = [85,55]
    # ref_neurone = [70,50]
    # ref_neurone = [30,5]
    # ref_neurone = [0,99]
    # ref_neurone = [12,62]
    # ref_neurone = [,]
    # ref_neurone = [99,82]
    # ref_neurone = [75,99]
    # ref_neurone = [52,99]

    Time_delay = np.arange(-20,30,1)
    number_of_annulus = 10
    interval = 20
    r = 30
    # RdBu = cm.get_cmap('RdBu', 256)
    header='/home/margauxvrech/mutual_network/DATAVSD/'
    newheader=folder+file +'data'
    if not os.path.exists(header+newheader):
        os.makedirs(header+newheader)
    print('\nAnalysing data of '+file)
    print('\nData will be saved in '+newheader)

    V=len(L)
    # Recorded_cell = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
    Recorded_cell_brut = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]



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
        return 0.014*(np.tanh(300*x-2)+1)

    def accentuation2(x):
        return 0.014*(np.tanh((300*x-2))+1)
    # filtered_L=[]
    # for t in range(511):
    #     filtered_L.append(signal.convolve2d(accentuation2(np.array(L[t])),Kernel1, mode='same'))
    # vm=filtered_L

    # vm_base=[vm[t][50][50] for t in range(injection_start,injection_start+interval)]
    vm_base_brut=[L[t][ref_neurone[0]][ref_neurone[1]] for t in range(min(injection_start,injection_end-1),min(injection_start+interval,injection_end))]
    # c_X,xedges = np.histogram(vm_base,200,range=(-0.1,0.02))
    c_X_brut,xedges = np.histogram(vm_base_brut,500,range=(-0.1,0.02))

    Time_maximum = []
    list_rho = []
    Time_maximum1 = []
    list_rho1 = []
    # Time_maximum2 = []
    # list_rho2 = []

    Annulus = [[] for i in range(number_of_annulus)]



    #this is the ray of the disk around our central cell
    print('... the width of the annulus is '+ str(r/number_of_annulus))
    #maybe no need to go to the border of the image
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





    print('Calculating MI...')




    fig=plt.figure()
    for A in Annulus :
        MI_annulus = []
        for time_delay in Time_delay :
            MI_delay = []
            for neuron in A :
                vm_neurone_brut = [L[t][neuron[0]][neuron[1]] for t in range(min(injection_start+time_delay,injection_end-1),min(injection_start+interval+time_delay,injection_end))]
                c_Y_brut,xedges = np.histogram(vm_neurone_brut,500,range=(-0.1,0.02))
                MI_delay.append(mutual_info_score(c_X_brut,c_Y_brut))
            MI_annulus.append(np.mean(MI_delay, dtype=np.float64))
        #

        MI_annulus_filtered=savgol_filter(MI_annulus, 11, 3)
        # MI_annulus_filtered=MI_annulus
        if Annulus.index(A)>1:
            maxi = max(MI_annulus_filtered)
            Time_maximum.append(Time_delay[list(MI_annulus_filtered).index(maxi)])
            list_rho.append(Annulus.index(A)* r/number_of_annulus)
            maxi1 = max(MI_annulus)
            Time_maximum1.append(Time_delay[list(MI_annulus).index(maxi1)])
            list_rho1.append(Annulus.index(A))
            plt.scatter(Time_delay,MI_annulus, color = cm.hsv(Annulus.index(A)/len(Annulus)),marker='x')#,linestyle='-.')
            plt.plot(Time_delay,MI_annulus_filtered,label='Anneau '+str(Annulus.index(A)), color = cm.hsv(Annulus.index(A)/len(Annulus)))#,linestyle='-.')
    plt.title('Mutual information')
    plt.xlabel('Time delay')
    plt.ylabel("MI")
    plt.legend(frameon=False)
    fig.savefig(header+newheader+'/MI_'+str(len(Annulus))+'_anneaux'+str(ref_neurone)+'.png')
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
    fig1.savefig(header+newheader+'/Annulus'+'.png')
    plt.close()
    fig1.clf()

    fig2 = plt.figure()

    plt.plot(list_rho1,Time_maximum1,'x')
    lr = scipy.stats.linregress(list_rho1,Time_maximum1)
    y=[lr[0]*i+lr[1] for i in list_rho1]
    plt.plot(list_rho1,y,c='r', label="pente="+str(lr[0])+", R2="+str(lr[2]**2))
    plt.show()
    plt.legend()
    plt.title('Maximum brut')


    fig2.savefig(header+newheader+'/Maximum'+'.png')
    plt.close()
    fig2.clf()

    fig2 = plt.figure()

    plt.plot(list_rho,Time_maximum,'x')
    lr = scipy.stats.linregress(list_rho,Time_maximum)
    y=[lr[0]*i+lr[1] for i in list_rho1]
    plt.plot(list_rho,y,c='r', label="pente="+str(lr[0])+", R2="+str(lr[2]**2))
    plt.show()
    plt.legend()
    plt.title('Maximum savgol')


    fig2.savefig(header+newheader+'/Maximum_savgol'+'.png')
    plt.close()
    fig2.clf()

    carac_time = np.mean([Time_maximum1[i+1]-Time_maximum1[i] for i in range(len(Time_maximum1)-1)])
    print("Le temps caracteristique est "+ str(carac_time)+ 'et la vitesse est donc :'+str(round(lr[0]**2*carac_time/100,3)))
