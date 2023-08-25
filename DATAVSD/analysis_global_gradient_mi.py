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
Grad_list = []
MI_list = []
#MWill contain each tab we get for each file
index = 0 # at which mouse and record are we looking at
for folder in  ['Mouses_3-5-6/20160912/','Mouses_3-5-6/20160914/','Mouses_3-5-6/20160916/']:
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
        newheader='/Analyse_globale'
        if not os.path.exists(header+newheader):
            os.makedirs(header+newheader)
        print('\nAnalysing data of '+file)
        print('\nData will be saved in '+newheader)
        Time_delay = np.arange(-5,25,1)
        V=len(L)
        Recorded_cell_brut = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
        Norm_vect = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
        vm_base_brut=[L[t][50][50] for t in range(injection_start,injection_start+interval)]
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

            #
            # fig=plt.figure()
            #
            #
            # plt.imshow(Recorded_cell_brut[injection_start+time_delay], cmap = 'inferno',interpolation='none')
            # plt.clim([0,0.028])
            # plt.colorbar()
            # plt.title('MI brut '+str(injection_start+time_delay))
            #
            [U,V]=np.gradient(Recorded_cell_brut[injection_start+time_delay])
            # # plt.quiver(X,Y,-V,U,color='white')
            #
            # plt.show(block=True)
            # filename= '/Gradient_'+str(injection_start+time_delay)+'.png'
            # fig.savefig(header+newheader+filename, transparent=True)
            #
            # plt.close()
            # fig.clf()

            for i in range(window):
                for j in range(window):

                    Norm_vect[int((injection_start+time_delay)/dt)][i][j]=np.linalg.norm([U[i][j],V[i][j]])
            # fig=plt.figure()
            # plt.imshow(Norm_vect[int((injection_start+time_delay)/dt)], interpolation = 'none',vmin=0,vmax= 0.01,cmap= 'YlGnBu')
            # plt.colorbar()
            # filename= '/Norm_vect_'+str(injection_start+time_delay)+'.png'
            # fig.savefig(header+newheader+filename)
            # plt.close()
            # fig.clf()

            #Plotting mean norm of each gradient vector

        Grad_list.append([[np.mean([Norm_vect[int((injection_start+time_delay)/dt)][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])
        MI_list.append([[np.mean([Recorded_cell_brut[int((injection_start+time_delay)/dt)][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])

        fig = plt.figure()
        plt.imshow([[MI_list[index][i][j] for j in range(window)] for i in range(window)])
        plt.colorbar()
        filename= '/Mean_gradient_MI_'+str(index)+'.png'
        fig.savefig(header+newheader+filename)
        plt.close()
        fig.clf()

        index+= 1


N = len(Grad_list)
fig=plt.figure()
plt.imshow([[np.mean([Grad_list[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
plt.clim([0,0.018])
plt.colorbar()
filename= '/Global_mean_gradient.png'
fig.savefig(header+newheader+filename)
plt.close()
fig.clf()

fig=plt.figure()
plt.imshow([[np.mean([MI_list[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
plt.colorbar()
filename= '/Global_MI_gradient.png'
fig.savefig(header+newheader+filename)
plt.close()
fig.clf()
