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
number_of_annulus = 21

# # folder = 'Single C2 whisker evoked responses, second set of 10 trials/'
# folder = 'Single_C2_whisker_evoked_responses_first_10_trials/'
# # for i in [11,12,13,14,15,16,17,18,19,20]:
# for i in range(1,11):
#     file = 'C2_'+str(i)
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
    interval = 5
    x = 0
    y = 0
    window = 100
    list_coord = [(x+i)*size+(y+j) for j in range(window) for i in range(window)]
    injection_points=[5050]
    min_point, max_point = 50,50
    ref_neurone = [50,50]


    header='/home/margauxvrech/mutual_network/DATAVSD/'
    newheader=folder+file +'data'
    if not os.path.exists(header+newheader):
        os.makedirs(header+newheader)
    print('\nAnalysing data of '+file)
    print('\nData will be saved in '+newheader)
    Time_delay = np.arange(-5,20,1)
    V=len(L)
    Recorded_cell_brut = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
    Recorded_cell_correlate = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
    VM = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]



    X, Y = np.meshgrid([i for i in range(window)],[i for i in range(window)])

    def accentuation1(x):
        return 0.014*(np.tanh(300*x-2)+1)


    vm_base_brut=[L[t][50][50] for t in range(injection_start,injection_start+interval)]
    c_X_brut,xedges = np.histogram(vm_base_brut,200,range=(-0.1,0.02))


    print('Plotting VM and MI...')
    for time_delay in Time_delay:
        for i in range(window):
            for j in range(window):

                vm_neurone_brut = [L[t][i][j] for t in range(injection_start+time_delay,injection_start+interval+time_delay)]

                c_Y_brut,xedges = np.histogram(vm_neurone_brut,200,range=(-0.1,0.02))
                Recorded_cell_brut[injection_start+time_delay][i][j]= mutual_info_score(c_X_brut,c_Y_brut)
                Recorded_cell_correlate[injection_start+time_delay][i][j]= float(signal.correlate(vm_base_brut,vm_neurone_brut,mode = 'valid'))
                VM[injection_start+time_delay][i][j] = np.mean([L[t][i][j] for t in range(injection_start+time_delay,injection_start+interval+time_delay)])

        fig=plt.figure()
        fig.add_subplot(1,3,1)
        plt.title('VM '+str(injection_start+time_delay))
        plt.imshow(VM[injection_start+time_delay], cmap = RdBu,interpolation='none', vmin =0 , vmax =0.015 )
        plt.colorbar()
        fig.add_subplot(1,3,2)
        plt.imshow(Recorded_cell_correlate[injection_start+time_delay],cmap = 'magma',interpolation='none')
        plt.clim([0,0.0001])
        plt.colorbar()
        plt.title('CC')

        fig.add_subplot(1,3,3)
        plt.imshow(Recorded_cell_brut[injection_start+time_delay], cmap = 'inferno',interpolation='none')
        plt.clim([0,0.028])
        plt.colorbar()

        plt.title('MI '+str(injection_start+time_delay))

        filename= '/VM_MI_CC_'+str(injection_start+time_delay)+'.png'
        fig.savefig(header+newheader+filename, transparent=True)
        plt.show(block=True)
        plt.close()
        fig.clf()

    #########################################################
    #ANNULUS ANALYSIS
    ######################################################

    # print('\nPlotting the Annulus')
    #
    #
    #
    # number_of_annulus = 21
    # interval = 30

    # Time_maximum1 = []
    # list_rho1 = []
    # # Time_maximum2 = []
    # # list_rho2 = []
    #
    # Annulus = [[] for i in range(number_of_annulus)]
    #
    # r = np.linalg.norm([50,50])
    #
    # #this is the ray of the disk around our central cell
    # print('... the width of the annulus is '+ str(r/number_of_annulus))
    #
    # for a in range(number_of_annulus):
    #     for i in range(100):
    #         for j in range(100):
    #             d = np.linalg.norm([i-50,j-50])
    #             for a in range(number_of_annulus):
    #                 if a* r/number_of_annulus<d<= (a+1)*r/number_of_annulus:
    #                     Annulus[a].append([i,j])
    #                     break
    #
    # Annulus = list(filter(([]).__ne__, Annulus))
    #
    #
    #
    # fig=plt.figure()
    # for A in Annulus :
    #     MI_annulus = []
    #     for time_delay in Time_delay :
    #         MI_delay = []
    #         for neuron in A :
    #             MI_delay.append(Recorded_cell_brut[injection_start+time_delay][neuron[0]][neuron[1]])
    #         MI_annulus.append(np.mean(MI_delay, dtype=np.float64))

    #
    #
    #     maxi1 = max(MI_annulus)
    #     Time_maximum1.append(Time_delay[list(MI_annulus).index(maxi1)])
    #     list_rho1.append(Annulus.index(A))
    #     plt.scatter(Time_delay,MI_annulus, color = cm.hsv(Annulus.index(A)/len(Annulus)),marker='x')#,linestyle='-.')
    #     plt.plot(Time_delay,MI_annulus_filtered,label='Anneau '+str(Annulus.index(A)), color = cm.hsv(Annulus.index(A)/len(Annulus)))#,linestyle='-.')
    # plt.title('Mutual information')
    # plt.xlabel('Time delay')
    # plt.ylabel("MI")
    # plt.legend(frameon=False)
    # fig.savefig(header+newheader+'/MI_'+str(len(Annulus))+' anneaux'+'.png')
    # plt.close()
    # fig.clf()
    # #plot the Annulus
    # Annulus_plot = np.zeros((window, window))
    # for A in Annulus:
    #     for neuron in A :
    #         Annulus_plot[neuron[0]][neuron[1]]=Annulus.index(A)+2
    #
    # fig1 = plt.figure()
    # plt.imshow(Annulus_plot, cmap = 'inferno',interpolation = 'none')
    # plt.colorbar()
    # plt.title('Annulus')
    # fig1.savefig(header+newheader+'/Annulus'+'.png')
    # plt.close()
    # fig1.clf()
    #
    # fig2 = plt.figure()
    #
    # plt.plot(list_rho1,Time_maximum1,'x')
    # lr = scipy.stats.linregress(list_rho1,Time_maximum1)
    # y=[lr[0]*i+lr[1] for i in list_rho1]
    # plt.plot(list_rho1,y,c='r', label="pente="+str(lr[0])+", R2="+str(lr[2]**2))
    # plt.show()
    # plt.legend()
    # plt.title('Maximum brut')
    #
    #
    # fig2.savefig(header+newheader+'/Maximum'+'.png')
    # plt.close()
    # fig2.clf()
