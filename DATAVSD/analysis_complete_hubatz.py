import os
import numpy as np
import scipy
import scipy.signal as signal
import cv2 as cv

################################
import matplotlib
matplotlib.use('Agg') # to be used when DISPLAY is undefined
import matplotlib.pyplot as plt
from matplotlib import cm


RdBu = cm.get_cmap('RdBu_r', 256)

################################

from scipy.signal import savgol_filter
from sklearn.metrics import mutual_info_score

# #  To create an ods file containing all the mean VM
# import jpype
# import asposecells
# jpype.startJVM()
# from asposecells.api import Workbook, FileFormatType

Mean_VM = True
MI = False
Annulus = False
Gradient = False
Ellipse = True
S2_PP = True
treshold = [[-1 for _ in range(5)] for _ in range(5)] #Treshold to filter the data
# Recommended : ???

Time_delay = np.arange(7,25,1)
run_time = 50 # nombre d'images
size = 100 # taille de l'image
injection_start,injection_end = 0,50 #début de la stimulation
interval = 1 # intervalle sur lequel on calcule la MI


VM_list = [] #VM
VM_Whisker_all = [[[],[],[],[],[]] for _ in range(5)]
VM_Whisker_tbt = [[[[np.zeros((100,100)) for _ in Time_delay] for _ in range(5)] for _ in range(5)] for _ in range(4)]

if MI :
    MI_list = [] #MI
if Gradient :
    VM_list_smooth = []
    if MI :
        MI_list_smooth = []

if Annulus :
    ref_neurone = [50,50] #neurone de référence
    number_of_annulus = 21



x = 0
y = 0
window = 100
list_coord = [(x+i)*size+(y+j) for j in range(window) for i in range(window)]



if Gradient :
    Grad_list_VM = [] #Gradient de la "VM"
    Grad_list_VM_smooth = []
    if MI :
        Grad_list_MI = [] #Gradient de la MI
        Grad_list_MI_smooth = []


def eq_S1S2(x):
    return 0.35*x+20
def eq_S1PP(x):
    return -0.23*x+83


if Ellipse :
    rise = 1.020
    ray = 10
    area_limit = 1000
    if S2_PP :
        area_PP_min = 200
        area_PP_max = 1000000000
        area_S2_min = 200
        area_S2_max = 1000000000
    color_ellipse = (100,100,100)
    Ellipses_mouse_VM = []
    Mean_vm_whisker = [[[[] for i in range(5)] for j in range(5)] for mouse in range(4)]
    Mask = np.zeros((100,100))
    MaskS2 = np.zeros((100,100))
    MaskPP = np.zeros((100,100))
    Dist_foyer = [[[[] for _ in Time_delay] for _ in range(5)] for _ in range(5)]
    List_thet_tbt = [[[[] for _ in Time_delay] for i in range(5)] for j in range(5)]
    List_thet_tbt_tri = [[[[] for _ in Time_delay] for i in range(5)] for j in range(5)]
    mean_dist_foyer = 19
    # Recommended : ???
    mean_dist_foyer_tbt = [[[3 for _ in Time_delay] for _ in range(5)] for _ in range(5)]
    # Recommended : ???
    for x in range(100):
        for y in range(100):
            if eq_S1PP(x) < y:
                MaskPP[y,x] = 255
            elif eq_S1PP(x) > y > eq_S1S2(x) :
                Mask[y,x] = 255
            else :
                MaskS2[y,x] = 255
    fig_whisker_foyer, axs_w_f = plt.subplots(5, 5)
    fig=plt.figure()
    plt.imshow(Mask+0.5*MaskS2+0.2*MaskPP)
    filename= '/Mask'+'.png'
    fig.savefig('/home/margauxvrech/mutual_network/DATAVSD/Analyse_complete_Hubatz_1/'+filename, transparent=True)
    plt.close()
    fig.clf()

    if MI :
         Ellipses_mouse_MI = []
    if Gradient :
        Ellipses_mouse_grad_VM = []
        if MI :
            Ellipses_mouse_grad_MI = []


header='/home/margauxvrech/mutual_network/DATAVSD/'
data_header='/Analyse_complete_Hubatz_1/Data'
if not os.path.exists(header+data_header):
    os.makedirs(header+data_header)
print('\nNumeric data will be saved in '+data_header)



data = open(header+data_header +"/complete_analysis.txt", 'w')
data.write('Complete analysis\n Analysis options :\n')
data.write('Mean_VM = ' + str(Mean_VM) + '\n' +
'MI = ' + str(MI) + '\n' +
'Annulus = ' + str(Annulus) + '\n' +
'Gradient = ' + str(Gradient) + '\n' +
'Ellipse = ' + str(Ellipse) + '\n' +
'S2_PP = ' + str(S2_PP) + '\n' +
'treshold = ' + str(treshold) + '\n')
if Annulus :
    data.write('Reference neuron : ' + str(ref_neurone) + '\n' +
    'Number of annulus = ' + str(number_of_annulus) + '\n')
if Ellipse :
    data.write('Number of rays : ' + str(ray) + '\n')
data.close()
# if Mean_VM :
#     VMbook = Workbook(FileFormatType.ODS)
#     VMsheet = VMbook.getWorksheets().get(0)

Mouse_index = 0
nb_previous_records = 0
nb_of_records = 0

fig_thet, axs_thet = plt.subplots(5,5)
fig_thet_tri, axs_thet_tri = plt.subplots(5,5)


#
for folder in [ 'Hubatz_data/E1','Hubatz_data/E2','Hubatz_data/E3','Hubatz_data/E4']:
    Mouse = Mouse_index+1
    nb_previous_records = nb_of_records+nb_previous_records
    nb_of_records = 0
    VM_Whisker = [[[],[],[],[],[]] for i in range(5)]
    if Ellipse :
        Ellipses_mean_VM = []
        if MI :
            Ellipses_mean_MI = []
        if Gradient :
            Ellipses_grad_VM = []
            if MI :
                Ellipses_grad_MI = []
    # if Mean_VM :
    #     column = str(chr(str(65+Mouse_index)))
    #     VMsheet.getCells().get(column+"1").putValue("Mouse "+ str(Mouse))
    for row in range(1,6):
        for line in range(1,6):
            mouse_header = '/Analyse_complete_Hubatz_1/Whisker_'+str(row)+'_'+str(line)+'/Mouse_'+str(Mouse)
            if row == 1 and line == 5 :
                break
            for exp in range(1,11):
                if folder[-1]=='2' and exp == 10:
                    break
                file = '/Stim_{}_{}_{}'.format(row,line,exp)
                L=[[[0 for k in range(100)] for l in range(100)] for i in range(50)]
                with open(folder+file+'.txt', 'r') as fr:
                    lines = fr.readlines()
                    for t in range(50):
                        for i in range(100):
                            L[t][i]=[float(lines[t*100+j].split('\t')[i]) for j in range(100)]

                ######################################################
                # Création d'un dossier associé à chaque souris
                ######################################################

                newheader='/Analyse_complete_Hubatz_1/Whisker_'+str(row)+'_'+str(line)+'/Mouse_'+str(Mouse)+file
                if not os.path.exists(header+newheader):
                    os.makedirs(header+newheader)
                print('\nAnalysing data of '+file)
                print('\nData will be saved in '+newheader)

                if not os.path.exists(header+data_header+'/Whisker_'+str(row)+'_'+str(line)+'/Mouse_'+str(Mouse)+file):
                    os.makedirs(header+data_header+'/Whisker_'+str(row)+'_'+str(line)+'/Mouse_'+str(Mouse)+file)

                min_min_vm, min_min_mi = 10000,10000
                max_max_vm, max_max_mi = -10000,-10000

                V=len(L)
                VM = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                VM_smooth = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                if MI :
                    Recorded_cell_brut = [[[0 for k in range(100)] for l in range(100)] for t in range(511)] #Raw MI
                    Recorded_cell_smooth = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                if Ellipse :
                    Ellipses_VM = []
                    Ellipses_VM_S2 = []
                    Ellipses_VM_PP = []
                    Nb_elps = 0
                    Nb_elps_PP = 0
                    Nb_elps_S2 = 0
                    Dist = [[] for k in range(ray)]
                    Dist_S2 = [[] for k in range(ray)]
                    Dist_PP = [[] for k in range(ray)]
                    VM_ellipse = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                    VM_ellipse_S1 = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                    VM_ellipse_S2 = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                    VM_ellipse_PP = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                    if MI :
                        Ellipses_MI = []
                        Ellipses_MI_PP = []
                        Ellipses_MI_S2 = []
                        MI_ellipse = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]


                if Gradient :
                    Norm_vect_VM = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                    Norm_vect_VM_smooth= [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                    if MI :
                        Norm_vect_MI = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
                        Norm_vect_MI_smooth = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]



                X, Y = np.meshgrid([i for i in range(window)],[i for i in range(window)])

                if MI :
                    vm_base_brut=[L[t][50][50] for t in range(injection_start,injection_start+interval)]
                    c_X_brut,xedges = np.histogram(vm_base_brut,200,range=(-0.1,0.02))

                ### Affichage du VM moyen en fonction du temps

                VM_moyen=[np.mean([L[injection_start+t][i] for i in range(100)]) for t in Time_delay]

                if Mean_VM :
                    fig=plt.figure()
                    plt.plot(Time_delay,VM_moyen)
                    filename= '/VM_moyen'+'.png'
                    fig.savefig(header+newheader+filename, transparent=True)
                    plt.close()
                    fig.clf()
                    max_mean_vm = max(VM_moyen)
                    mean_vm = np.mean(VM_moyen)
                    print("The max value of the mean VM is "+str(max_mean_vm))
                    # VMsheet.getCells().get(column+str(1+row+line+exp)).putValue(max_mean_vm)


                ##########
                #Création des anneaux
                ##########

                if Annulus :
                    Annulus = [[] for i in range(number_of_annulus)]
                    r = np.linalg.norm(ref_neurone)

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

                    #Plotting them for visualisation

                    Annulus_plot = np.zeros((window, window))
                    for A in Annulus:
                        for neuron in A :
                            Annulus_plot[neuron[0]][neuron[1]]=Annulus.index(A)+2

                    fig1 = plt.figure()
                    plt.imshow(Annulus_plot, cmap = 'inferno',interpolation = 'none')
                    plt.colorbar()
                    plt.title('Annulus')
                    fig1.savefig(header+newheader+'/Annulus_MI'+'.png')
                    plt.close()
                    fig1.clf()


                if max_mean_vm>=treshold[row-1][line-1]:
                    print("We analyse this data\n")
                    nb_of_records+=1


                    max_vm,min_vm = -1000,1000
                    max_mi,min_mi = -1000,1000

                    if Annulus :
                        vm_base_brut=[L[t][ref_neurone[0]][ref_neurone[1]] for t in range(injection_start,injection_start+interval)]

                    for count,time_delay in enumerate(Time_delay):
                        for i in range(window):
                            for j in range(window):
                                vm_neurone_brut = [L[t][i][j] for t in range(injection_start+time_delay,injection_start+interval+time_delay)]
                                VM[injection_start+time_delay][i][j] = np.mean([L[t][i][j] for t in range(injection_start+time_delay,injection_start+interval+time_delay)])

                                if VM[injection_start+time_delay][i][j] < min_vm :
                                    min_vm = VM[injection_start+time_delay][i][j]
                                if VM[injection_start+time_delay][i][j] > max_vm :
                                    max_vm = VM[injection_start+time_delay][i][j]
                                if MI :
                                    c_Y_brut,xedges = np.histogram(vm_neurone_brut,200,range=(-0.1,0.02))
                                    Recorded_cell_brut[injection_start+time_delay][i][j]= mutual_info_score(c_X_brut,c_Y_brut)

                                    if Recorded_cell_brut[injection_start+time_delay][i][j] < min_mi :
                                        min_mi = Recorded_cell_brut[injection_start+time_delay][i][j]
                                    if Recorded_cell_brut[injection_start+time_delay][i][j] > max_mi :
                                        max_mi = Recorded_cell_brut[injection_start+time_delay][i][j]
                        if Mouse_index == 1 :
                            coef = 9
                        else :
                            coef = 10
                        VM_Whisker_tbt[Mouse_index][row-1][line-1][count] += np.array(VM[injection_start+time_delay])/coef
                    if min_vm < min_min_vm :
                        min_min_vm = min_vm
                    if max_vm > max_max_vm:
                        max_max_vm = max_vm
                    if MI :
                        if min_mi < min_min_mi :
                            min_min_mi = min_mi
                        if max_mi > max_max_mi:
                            max_max_mi = max_mi

                    if MI :
                        mean_mi = np.mean([Recorded_cell_brut[injection_start+time_delay][i][j] for i in range(window) for j in range(window) for time_delay in Time_delay])
                    if Annulus :
                        print('Calculating the propagation of VM & MI...')
                        Time_maximum = []
                        list_rho = []
                        if MI :
                            Time_maximum1 = []
                            list_rho1 = []


                        fig_a, axs_1 = plt.subplots(1, 1)
                        fig_b, axs_2 = plt.subplots(1, 1)

                        for A in Annulus :
                            MI_annulus = []
                            VM_annulus = []
                            for time_delay in Time_delay :
                                MI_delay = []
                                VM_delay = []
                                for neuron in A :
                                    VM_delay.append(VM[injection_start+time_delay][neuron[0]][neuron[1]])
                                    if MI :
                                        MI_delay.append(Recorded_cell_brut[injection_start+time_delay][neuron[0]][neuron[1]])
                                VM_annulus.append(np.mean(VM_delay, dtype=np.float64))
                                if MI :
                                    MI_annulus.append(np.mean(MI_delay, dtype=np.float64))

                            maxi = max(VM_annulus)
                            Time_maximum.append(Time_delay[list(VM_annulus).index(maxi)])
                            list_rho.append(Annulus.index(A))
                            axs_1.plot(Time_delay,VM_annulus, color = cm.hsv(Annulus.index(A)/len(Annulus)),marker='x')
                            # MI
                            if MI :
                                maxi1 = max(MI_annulus)
                                Time_maximum1.append(Time_delay[list(MI_annulus).index(maxi1)])
                                list_rho1.append(Annulus.index(A))
                                axs_2.plot(Time_delay,MI_annulus, color = cm.hsv(Annulus.index(A)/len(Annulus)),marker='x')


                        fig_a.suptitle('VM')
                        fig_a.savefig(header+newheader+'/VM_'+str(len(Annulus))+'_annulus'+'.png')
                        plt.close()
                        fig_a.clf()
                        if MI :
                            fig_b.suptitle('MI')
                            fig_b.savefig(header+newheader+'/MI_'+str(len(Annulus))+'_annulus'+'.png')
                            plt.close()
                            fig_b.clf()
                        #plot the Annulus

                        #OLS
                        fig2 = plt.figure()

                        plt.plot(list_rho,Time_maximum,'x')
                        lr = scipy.stats.linregress(list_rho,Time_maximum)
                        y=[lr[0]*i+lr[1] for i in list_rho]
                        plt.plot(list_rho,y,c='r', label="pente="+str(lr[0])+", R2="+str(lr[2]**2))
                        plt.show()
                        plt.legend()
                        plt.title('Raw Maximum VM')


                        fig2.savefig(header+newheader+'/Raw_maximum_VM'+'.png')
                        plt.close()
                        fig2.clf()


                        #plot the Annulus

                        #OLS
                        if MI :
                            fig2 = plt.figure()

                            plt.plot(list_rho1,Time_maximum1,'x')
                            lr = scipy.stats.linregress(list_rho1,Time_maximum1)
                            y=[lr[0]*i+lr[1] for i in list_rho1]
                            plt.plot(list_rho1,y,c='r', label="pente="+str(lr[0])+", R2="+str(lr[2]**2))
                            plt.show()
                            plt.legend()
                            plt.title('Raw Maximum MI')


                            fig2.savefig(header+newheader+'/Raw_maximum_MI'+'.png')
                            plt.close()
                            fig2.clf()


                    print('Plotting VM and MI...')
                    first_ellipse = False
                    first_ellipse_PP = False
                    first_ellipse_S2 = False
                    S1_stop = False
                    S2_stop = False
                    PP_stop = False
                    abscisse_thet = [0 for _ in Time_delay]


                    for index,time_delay in enumerate(Time_delay) :
                    ###########
                    # AFFICHAGE VM ET MI
                    ###########
                        if MI :
                            fig, axs = plt.subplots(1, 2)

                            axs[0].imshow(VM[injection_start+time_delay], cmap = RdBu,interpolation='none', vmin =min_vm , vmax =max_vm)
                            axs[0].set_title('VM '+str(injection_start+time_delay))
                            axs[1].imshow(Recorded_cell_brut[injection_start+time_delay], cmap = 'inferno',interpolation='none', vmin = min_mi, vmax = max_mi)
                            axs[1].set_title('MI '+str(injection_start+time_delay))
                            filename= '/VM_MI_'+str(injection_start+time_delay)+'.png'
                        else :
                            fig, axs = plt.subplots(1, 1)
                            axs.imshow(VM[injection_start+time_delay], cmap = RdBu,interpolation='none', vmin =min_vm , vmax =max_vm)
                            axs.set_title('VM '+str(injection_start+time_delay))

                            filename= '/VM_'+str(injection_start+time_delay)+'.png'
                        fig.savefig(header+newheader+filename, transparent=True)
                        plt.show(block=True)
                        plt.close()
                        fig.clf()
                        if Gradient :
                            VM_smooth[injection_start+time_delay] = [[min_vm if VM[injection_start+time_delay][i][j]<mean_vm else max_vm for j in range(100)] for i in range(100)]
                            if MI :
                                Recorded_cell_smooth[injection_start+time_delay] = [[min_mi if Recorded_cell_brut[injection_start+time_delay][i][j]<mean_mi else max_mi for j in range(100)] for i in range(100)]

                        ###########
                        #CALCUL GRADIENT VM ET MI
                        ###########

                        if Gradient :
                            [U,V]=np.gradient(VM[injection_start+time_delay])
                            [U2,V2]=np.gradient(VM_smooth[injection_start+time_delay])
                            if MI :
                                [U1,V1]=np.gradient(Recorded_cell_brut[injection_start+time_delay])
                                [U3,V3]=np.gradient(Recorded_cell_smooth[injection_start+time_delay])

                            for i in range(window):
                                for j in range(window):
                                    Norm_vect_VM[int((injection_start+time_delay))][i][j]=np.linalg.norm([U1[i][j],V1[i][j]])
                                    Norm_vect_VM_smooth[int((injection_start+time_delay))][i][j]=np.linalg.norm([U2[i][j],V2[i][j]])
                                    if MI :
                                        Norm_vect_MI[int((injection_start+time_delay))][i][j]=np.linalg.norm([U1[i][j],V1[i][j]])
                                        Norm_vect_MI_smooth[int((injection_start+time_delay))][i][j]=np.linalg.norm([U3[i][j],V3[i][j]])
                        if Ellipse :

                            Dist_foyer[row-1][line-1][index].append(0)

                            VM_ellipse[injection_start+time_delay] = np.array( (VM[injection_start+time_delay] - min_vm ) * 255/(max_vm-min_vm) , dtype = np.uint8)
                            blur_vm = cv.GaussianBlur(VM_ellipse[injection_start+time_delay],(9,9),0)
                            ret1,VM_ellipse[injection_start+time_delay] = cv.threshold(blur_vm,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                            ret1,VM_ellipse[injection_start+time_delay] = cv.threshold(blur_vm,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)

                            np.savetxt(header+newheader+'/VM_'+str(time_delay)+'.txt',VM_ellipse[injection_start+time_delay])
                            if MI :
                                MI_ellipse[injection_start+time_delay] = np.array( (Recorded_cell_brut[injection_start+time_delay] - min_MI ) * 255/(max_MI-min_MI) , dtype = np.uint8)
                                blur_mi = cv.GaussianBlur(MI_ellipse[injection_start+time_delay],(9,9),0)
                                #ret2,MI_smooth_treshold = cv.threshold(blur_mi,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                                ret2,MI_ellipse[injection_start+time_delay] = cv.threshold(blur_mi,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                                ret2,MI_ellipse = cv.threshold(blur_mi,int(np.floor(ret2*rise)),255,cv.THRESH_BINARY)
                                np.savetxt(header+newheader+'/MI'+str(time_delay)+'.txt',MI_ellipse[injection_start+time_delay])
                            for i in range(100):
                                for j in range(100):
                                    VM_ellipse_S1[injection_start+time_delay][i][j] = min(VM_ellipse[injection_start+time_delay][i][j],Mask[i][j])
                                    VM_ellipse_S2[injection_start+time_delay][i][j] = min(VM_ellipse[injection_start+time_delay][i][j],MaskS2[i][j])
                                    VM_ellipse_PP[injection_start+time_delay][i][j] = min(VM_ellipse[injection_start+time_delay][i][j],MaskPP[i][j])

                            contours1,hierarchy1 = cv.findContours(np.array(VM_ellipse_S1[injection_start+time_delay]).astype(np.uint8), 1, 2)
                            if S2_PP :
                                contoursS2,hierarchyS2 = cv.findContours(np.array(VM_ellipse_S2[injection_start+time_delay]).astype(np.uint8), 1, 2)
                                contoursPP,hierarchyPP = cv.findContours(np.array(VM_ellipse_PP[injection_start+time_delay]).astype(np.uint8), 1, 2)
                            if MI :
                                contours2,hierarchy2 = cv.findContours(MI_ellipse[injection_start+time_delay], 1, 2)

                            if len(contours1) != 0:
                                cont = max(contours1,key = len)
                                if len(cont) >= 5:
                                    elps = cv.fitEllipse(cont)
                                    area = np.pi*elps[1][0]*elps[1][1]/4
                                    if area < area_limit and not S1_stop :
                                        if not first_ellipse :
                                            first_ellipse =True
                                            centre_S1 = (int(round(elps[0][0])),int(round(elps[0][1])))
                                            theta_S1 = elps[2]
                                        Nb_elps+=1
                                        abscisse_thet[index]+=1
                                        Ellipses_VM.append(elps)
                                        Ellipses_mouse_VM.append(elps)
                                        if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                            cv.ellipse(VM_ellipse[injection_start+time_delay], elps, color_ellipse)
                                            center = (int(round(elps[0][0])),int(round(elps[0][1])))
                                            axes = (int(round(elps[1][0]/2)),int(round(elps[1][1]/2)))
                                            theta = elps[2]
                                            cv.ellipse(VM_ellipse[injection_start+time_delay],center,axes,theta ,0. ,3.,(170, 170, 170), 1)
                                            ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                                            bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                                            ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                                            by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                                            for count,angle in enumerate([k*2*np.pi/ray for k in range(ray)]):
                                                point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))))
                                                # Distance par rapport au centre de chaque ellipse
                                                # cv.line(VM_ellipse[injection_start+time_delay][:][20:80],(int(round(elps[0][0])),int(round(elps[0][1]))),point, color_ray)
                                                # Dist[count].append(np.linalg.norm([int(round(elps[0][0])),int(round(elps[0][1]))]-np.array(point)))
                                                cv.line(VM_ellipse[injection_start+time_delay],centre_S1,point, (255-255/ray*count,255-255/ray*count,255-255/ray*count))
                                                Dist[count].append(np.linalg.norm(np.array(centre_S1)-np.array(point)))
                                    else :
                                        Ellipses_VM.append(0)
                                        if Nb_elps != 0 :
                                            S1_stop = True
                                else:
                                    Ellipses_VM.append(0)
                            else :
                                Ellipses_VM.append(0)
                            if S2_PP :
                                if len(contoursS2) != 0:
                                    cont = max(contoursS2,key = len)
                                    if len(cont) >= 5:
                                        elps = cv.fitEllipse(cont)
                                        area_S2 = np.pi*elps[1][0]*elps[1][1]/4
                                        if area_S2_min < area_S2 < area_S2_max and not S2_stop:
                                            if not first_ellipse_S2 :
                                                first_ellipse_S2 = True
                                                centre_S2 = (int(round(elps[0][0])),int(round(elps[0][1])))
                                                theta_S2 = elps[2]*2*np.pi/360
                                            Nb_elps_S2+=1
                                            Ellipses_VM_S2.append(elps)
                                            if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                                cv.ellipse(VM_ellipse[injection_start+time_delay], elps, color_ellipse)
                                                center = (int(round(elps[0][0])),int(round(elps[0][1])))
                                                axes = (int(round(elps[1][0]/2)),int(round(elps[1][1]/2)))
                                                theta = elps[2]
                                                cv.ellipse(VM_ellipse[injection_start+time_delay],center,axes,theta ,0. ,3.,(170, 170, 170), 1)
                                                ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                                                bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                                                ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                                                by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                                                for count,angle in enumerate([k*2*np.pi/ray for k in range(ray)]):
                                                    point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_S2)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_S2)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_S2)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_S2)*2*np.pi/360))))
                                                    cv.line(VM_ellipse[injection_start+time_delay],centre_S2,point, (255-255/ray*count,255-255/ray*count,255-255/ray*count))
                                                    Dist_S2[count].append(np.linalg.norm(np.array(centre_S2)-np.array(point)))
                                        else :
                                            if Nb_elps_S2 != 0 :
                                                S2_stop = True
                                if len(contoursPP) != 0:
                                    cont = max(contoursPP,key = len)
                                    if len(cont) >= 5:
                                        elps = cv.fitEllipse(cont)
                                        area_PP = np.pi*elps[1][0]*elps[1][1]/4
                                        if area_PP_min < area_PP < area_PP_max and not PP_stop:
                                            if not first_ellipse_PP :
                                                first_ellipse_PP = True
                                                centre_PP = (int(round(elps[0][0])),int(round(elps[0][1])))
                                                theta_PP = elps[2]*2*np.pi/360
                                            Nb_elps_PP+=1
                                            Ellipses_VM_PP.append(elps)
                                            if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                                cv.ellipse(VM_ellipse[injection_start+time_delay], elps, color_ellipse)
                                                center = (int(round(elps[0][0])),int(round(elps[0][1])))
                                                axes = (int(round(elps[1][0]/2)),int(round(elps[1][1]/2)))
                                                theta = elps[2]
                                                cv.ellipse(VM_ellipse[injection_start+time_delay],center,axes,theta ,0. ,3.,(170, 170, 170), 1)
                                                ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                                                bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                                                ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                                                by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                                                for count,angle in enumerate([k*360/ray for k in range(ray)]):
                                                    point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_PP)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_PP)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_PP)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_PP)*2*np.pi/360))))
                                                    cv.line(VM_ellipse[injection_start+time_delay],centre_PP,point, (255-255/ray*count,255-255/ray*count,255-255/ray*count))
                                                    Dist_PP[count].append(np.linalg.norm(np.array(centre_PP)-np.array(point)))
                                            else :
                                                if Nb_elps_PP != 0 :
                                                    PP_stop = True
                            if MI :
                                print("LA VERSION MI N EST PAS CODEE ")
                                if len(contours2) != 0:
                                    cont = max(contours2,key = len)
                                    if len(cont) >= 5:
                                        elps = cv.fitEllipse(cont)
                                        Ellipses_MI.append(elps)
                                        if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                            cv.ellipse(MI_ellipse[injection_start+time_delay], elps, color_ellipse, 2 )
                                fig, axs = plt.subplots(1, 2)

                                axs[0].imshow(VM_ellipse[injection_start+time_delay], cmap = RdBu,interpolation='none')
                                axs[0].set_title('VM treshold')
                                axs[1].imshow(MI_ellipse[injection_start+time_delay], cmap = 'inferno',interpolation='none')
                                axs[1].set_title('MI treshold')


                                filename= '/VM_MI_ellipse_'+str(injection_start+time_delay)+'.png'
                            else :
                                fig, axs = plt.subplots(1, 1)

                                axs.imshow(VM_ellipse[injection_start+time_delay], cmap = RdBu,interpolation='none')
                                axs.set_title('VM tresh Aire centrale')
                                filename= '/VM_ellipse_'+str(injection_start+time_delay)+'.png'
                            fig.savefig(header+newheader+filename, transparent=True)
                            plt.show(block=True)
                            plt.close()
                            fig.clf()

                    for index_1,time_delay in enumerate(Time_delay) :
                        if abscisse_thet[index_1]==1 :
                            ind = 1
                            while Ellipses_VM[index_1-ind] == 0 :
                                ind -= 1
                            Dist_foyer[row-1][line-1][index_1][-1] += np.linalg.norm(np.array([Ellipses_VM[index_1][0][0],Ellipses_VM[index_1][0][1]])-np.array([Ellipses_VM[index_1-ind][0][0],Ellipses_VM[index_1-ind][0][1]]))

                    fig,axs = plt.subplots(1, 1)
                    fig_foyer,axs_foyer = plt.subplots(1, 1)
                    for index_1,time_delay in enumerate(Time_delay) :
                        if abscisse_thet[index_1]==1 :
                            axs.scatter(time_delay,90-Ellipses_VM[index_1][2])
                            axs_foyer.scatter(Ellipses_VM[index_1][0][0],Ellipses_VM[index_1][0][1], color = cm.hsv(index_1/len(Time_delay)))
                            axs_w_f[row-1][line-1].scatter(Ellipses_VM[index_1][0][0],Ellipses_VM[index_1][0][1], color = cm.hsv(index_1/len(Time_delay)))
                            List_thet_tbt[row-1][line-1][index_1].append(90-Ellipses_VM[index_1][2])
                            if Dist_foyer[row-1][line-1][index_1][-1] < mean_dist_foyer_tbt[row-1][line-1][index_1] :
                                axs_thet[row-1][line-1].scatter(time_delay,90-Ellipses_VM[index_1][2], marker = 'x', color = 'b')
                                axs_thet_tri[row-1][line-1].scatter(time_delay,90-Ellipses_VM[index_1][2], marker = 'x', color = 'b')
                                List_thet_tbt_tri[row-1][line-1][index_1].append(90-Ellipses_VM[index_1][2])
                            else :
                                axs_thet[row-1][line-1].scatter(time_delay,90-Ellipses_VM[index_1][2], marker = 'x', color = 'g')

                    axs.set_ylim(-90,90)
                    axs.set_title('Evolution theta')
                    axs_w_f[row-1][line-1].set_title(str(row) + '_' + str(line))
                    axs_w_f[row-1][line-1].set_xlim(0,100)
                    axs_w_f[row-1][line-1].set_ylim(0,100)
                    axs_foyer.set_xlim(0,100)
                    axs_foyer.set_ylim(0,100)
                    axs_foyer.set_title('Evolution position des foyers: distance totale = ' + str(np.sum(Dist_foyer[row-1][line-1][-1])))
                    filename= '/Evol_theta.png'
                    filename_foyer = '/Evol_foyer.png'
                    fig.savefig(header+newheader+filename)
                    fig_foyer.savefig(header+newheader+filename_foyer)
                    plt.close()
                    fig.clf()
                    fig_foyer.clf()



                    VM_list.append([[np.mean([VM[int(injection_start+time_delay)][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])
                    VM_Whisker[row-1][line-1].append([[np.mean([VM[int(injection_start+time_delay)][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])
                    if MI :
                        MI_list.append([[np.mean([Recorded_cell_brut[int((injection_start+time_delay))][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])

                    vm = np.array([[VM_list[-1][i][j] for j in range(window)] for i in range(window)])
                    fig = plt.figure()
                    plt.imshow(vm, vmin =0 , vmax =0.009)
                    plt.colorbar()
                    filename= '/Mean_VM.png'
                    fig.savefig(header+newheader+filename)
                    plt.close()
                    fig.clf()
                    np.savetxt(header+data_header+'/Whisker_'+str(row)+'_'+str(line)+'/Mouse_'+str(Mouse)+file+'/Mean_VM.txt',vm)
                    if MI :
                        mi = np.array([[MI_list[-1][i][j] for j in range(window)] for i in range(window)])
                        fig = plt.figure()
                        plt.imshow(mi, vmin =0 , vmax =0.08 )
                        plt.colorbar()
                        filename= '/Mean_MI.png'
                        fig.savefig(header+newheader+filename)
                        plt.close()
                        fig.clf()
                        np.savetxt(header+data_header+'/Whisker_'+str(row)+'_'+str(line)+'/Mouse_'+str(Mouse)+file+'/Mean_MI.txt',mi)
                    any_ellipse = False
                    if Ellipse :
                        for elps in Ellipses_VM :
                            if elps != 0 :
                                cv.ellipse(vm, elps, color_ellipse, 1 )
                                any_ellipse = True
                        if any_ellipse:
                            elps  = 0
                            i = 1
                            while elps == 0 :
                                elps = Ellipses_VM[-i]
                                i+=1
                            ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                            bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                            ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                            by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                            for count,angle in enumerate([k*2*np.pi/ray for k in range(ray)]):
                                # point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))))
                                # cv.line(vm[:][20:80],(int(round(elps[0][0])),int(round(elps[0][1]))),point, color_ray)
                                point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))))
                                cv.line(vm,centre_S1,point, (0.009/ray*count,0.009/ray*count,0.009/ray*count,1))

                        fig = plt.figure()
                        plt.imshow(vm, vmin =0 , vmax =0.009)
                        plt.colorbar()
                        filename= '/Mean_VM_ellipses.png'
                        fig.savefig(header+newheader+filename)
                        plt.close()
                        fig.clf()
                        if MI :
                            for elps in Ellipses_MI :
                                cv.ellipse(mi[:][27:72], elps, color_ellipse, 1 )
                            fig = plt.figure()
                            plt.imshow(mi, vmin =0 , vmax =0.008)
                            plt.colorbar()
                            filename= '/Mean_MI_ellipses.png'
                            fig.savefig(header+newheader+filename)
                            plt.close()
                            fig.clf()

                        VM_list_tresh = np.array( (VM_list[-1] - min_min_vm ) * 255/(max_max_vm-min_min_vm) , dtype = np.uint8)
                        blur_vm = cv.GaussianBlur(VM_list_tresh,(9,9),0)
                        ret1,VM_list_treshold= cv.threshold(blur_vm,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                        ret1,VM_list_treshold= cv.threshold(blur_vm,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)

                        fig = plt.figure()
                        for k in range(ray):
                            plt.plot(Time_delay[:Nb_elps], Dist[k], label = 'Angle = ' + str(k*360/ray), color = cm.hsv(k/ray))
                        plt.legend(loc = 'upper left')
                        filename= '/Distance_VM_ellipse_S1_rays_'+str(ray)+'.png'
                        fig.savefig(header+newheader+filename)
                        plt.close()
                        fig.clf()
                        if S2_PP :
                            fig = plt.figure()
                            for k in range(ray):
                                plt.plot(Time_delay[:Nb_elps_S2], Dist_S2[k], label = 'Angle = ' + str(k*360/ray), color = cm.hsv(k/ray))
                            plt.legend(loc = 'upper left')
                            filename= '/Distance_VM_ellipse_S2_rays_'+str(ray)+'.png'
                            fig.savefig(header+newheader+filename)
                            plt.close()
                            fig.clf()
                            fig = plt.figure()
                            for k in range(ray):
                                plt.plot(Time_delay[:Nb_elps_PP], Dist_PP[k], label = 'Angle = ' + str(k*360/ray), color = cm.hsv(k/ray))
                            plt.legend(loc = 'upper left')
                            filename= '/Distance_VM_ellipse_PP_rays_'+str(ray)+'.png'
                            fig.savefig(header+newheader+filename)
                            plt.close()
                            fig.clf()
                        if MI :
                            MI_list_tresh = np.array( (MI_list[-1] - min_min_mi ) * 255/(max_max_mi-min_min_mi) , dtype = np.uint8)
                            blur_mi = cv.GaussianBlur(MI_list_tresh,(9,9),0)
                            ret1,MI_list_treshold = cv.threshold(blur_mi,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                            ret1,MI_list_treshold = cv.threshold(blur_mi,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)


                        contours1,hierarchy1 = cv.findContours(VM_list_treshold, 1, 2)
                        if MI :
                            contours2,hierarchy2 = cv.findContours(MI_list_treshold, 1, 2)
                        if len(contours1) != 0:
                            cont = max(contours1,key = len)
                            if len(cont) >= 5:
                                elps = cv.fitEllipse(cont)
                                Ellipses_mean_VM.append(elps)
                                if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                    cv.ellipse(VM_list_treshold, elps, color_ellipse, 2 )
                        if MI :
                            if len(contours2) != 0:
                                cont = max(contours2,key = len)
                                if len(cont) >= 5:
                                    elps = cv.fitEllipse(cont)
                                    Ellipses_mean_MI.append(elps)
                                    if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                        cv.ellipse(MI_list_treshold, elps, color_ellipse, 2 )
                        if MI :
                            fig, axs = plt.subplots(1, 2)
                            axs[0].imshow(VM_list_treshold,interpolation='none', vmin =0 , vmax =255)
                            axs[0].set_title('Mean VM treshold')
                            axs[1].imshow(MI_list_treshold, cmap = 'inferno',interpolation='none', vmin = 0, vmax = 255)
                            axs[1].set_title('Mean MI treshold')

                            filename= '/Mean VM_MI treshold.png'
                        else :
                            fig, axs = plt.subplots(1, 1)
                            axs.imshow(VM_list_treshold,interpolation='none', vmin =0 , vmax =255)
                            axs.set_title('Mean VM treshold')
                            filename= '/Mean VM treshold.png'
                        fig.savefig(header+newheader+filename, transparent=True)
                        plt.show(block=True)
                        plt.close()
                        fig.clf()


                    if Gradient :
                        Grad_list_VM.append([[np.mean([Norm_vect_VM[int(injection_start+time_delay)][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])
                        if MI :
                            Grad_list_MI.append([[np.mean([Norm_vect_MI[int(injection_start+time_delay)][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])

                        VM_list_smooth.append([[np.mean([VM_smooth[int((injection_start+time_delay))][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])
                        Grad_list_VM_smooth.append([[np.mean([Norm_vect_VM_smooth[int((injection_start+time_delay))][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])
                        if MI :
                            MI_list_smooth.append([[np.mean([Recorded_cell_smooth[int((injection_start+time_delay))][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])
                            Grad_list_MI_smooth.append([[np.mean([Norm_vect_MI_smooth[int((injection_start+time_delay))][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])


                        #Gradient VM & MI
                        grad_vm = np.array([[Grad_list_VM[-1][i][j] for j in range(window)] for i in range(window)])
                        fig = plt.figure()
                        plt.imshow(grad_vm)#, vmin =0 , vmax =0.003 )
                        plt.colorbar()
                        filename= '/Mean_gradient_Norm_VM.png'
                        fig.savefig(header+newheader+filename)
                        plt.close()
                        np.savetxt(header+data_header+'/Whisker_'+str(row)+'_'+str(line)+'/Mouse_'+str(Mouse)+file+'/Mean_gradient_Norm_VM.txt',grad_vm)
                        if MI :
                            grad_mi = np.array([[Grad_list_MI[-1][i][j] for j in range(window)] for i in range(window)])
                            fig = plt.figure()
                            plt.imshow(grad_mi)#, vmin =0 , vmax =0.03 )
                            plt.colorbar()
                            filename= '/Mean_gradient_Norm_MI.png'
                            fig.savefig(header+newheader+filename)
                            plt.close()
                            fig.clf()
                            np.savetxt(header+data_header+'/Whisker_'+str(row)+'_'+str(line)+'/Mouse_'+str(Mouse)+file+'/Mean_gradient_Norm_MI.txt',grad_mi)

                        #Gradient VM & MI smooth
                        grad_vm_smooth = np.array([[Grad_list_VM_smooth[-1][i][j] for j in range(window)] for i in range(window)])
                        fig = plt.figure()
                        plt.imshow(grad_vm_smooth)#, vmin =0 , vmax =0.003 )
                        plt.colorbar()
                        filename= '/Mean_gradient_Norm_VM_smooth.png'
                        fig.savefig(header+newheader+filename)
                        plt.close()
                        fig.clf()
                        np.savetxt(header+newheader+'/Mean_gradient_Norm_VM_smooth.txt',grad_vm_smooth)
                        if MI :
                            grad_mi_smooth = np.array([[Grad_list_MI_smooth[-1][i][j] for j in range(window)] for i in range(window)])
                            fig = plt.figure()
                            plt.imshow(grad_mi_smooth)#, vmin =0 , vmax =0.03 )
                            plt.colorbar()
                            filename= '/Mean_gradient_Norm_MI_smooth.png'
                            fig.savefig(header+newheader+filename)
                            plt.close()
                            fig.clf()
                            np.savetxt(header+newheader+'/Mean_gradient_Norm_MI_smooth.txt',grad_mi_smooth)

                        if Ellipse :
                            max_grad_vm = max([grad_vm_smooth[i][j] for i in range(window) for j in range(window)])
                            min_grad_vm = min([grad_vm_smooth[i][j] for i in range(window) for j in range(window)])
                            grad_vm_smooth = np.array((grad_vm_smooth - min_grad_vm ) * 255/(max_grad_vm-min_grad_vm) , dtype = np.uint8)
                            blur_grad_vm = cv.GaussianBlur(grad_vm_smooth,(9,9),0)
                            ret1,Grad_vm_treshold= cv.threshold(blur_grad_vm,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                            ret1,Grad_vm_treshold= cv.threshold(blur_grad_vm,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)
                            if MI :
                                max_grad_mi = max([grad_mi_smooth[i][j] for i in range(window) for j in range(window)])
                                min_grad_mi = min([grad_mi_smooth[i][j] for i in range(window) for j in range(window)])
                                grad_mi_smooth = np.array((grad_mi_smooth - min_grad_mi ) * 255/(max_grad_mi-min_grad_mi) , dtype = np.uint8)
                                blur_grad_mi = cv.GaussianBlur(grad_mi_smooth,(9,9),0)
                                ret1,Grad_mi_treshold= cv.threshold(blur_grad_mi,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                                ret1,Grad_mi_treshold= cv.threshold(blur_grad_mi,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)


                            contours1,hierarchy1 = cv.findContours(Grad_vm_treshold, 1, 2)
                            if MI :
                                contours2,hierarchy2 = cv.findContours(Grad_mi_treshold, 1, 2)
                            if len(contours1) != 0:
                                cont = max(contours1,key = len)
                                if len(cont) >= 5:
                                    elps = cv.fitEllipse(cont)
                                    Ellipses_grad_VM.append(elps)
                                    if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                        cv.ellipse(Grad_vm_treshold, elps, color_ellipse, 2 )
                            if MI :
                                if len(contours2) != 0:
                                    cont = max(contours2,key = len)
                                    if len(cont) >= 5:
                                        elps = cv.fitEllipse(cont)
                                        Ellipses_grad_MI.append(elps)
                                        if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                            cv.ellipse(Grad_mi_treshold, elps, color_ellipse, 2 )
                            if MI :
                                fig, axs = plt.subplots(1, 2)
                                axs[0].imshow(Grad_vm_treshold,interpolation='none', vmin =0 , vmax =255)
                                axs[0].set_title('Mean Gradient Norm VM treshold')
                                axs[1].imshow(Grad_mi_treshold, cmap = 'inferno',interpolation='none', vmin = 0, vmax = 255)
                                axs[1].set_title('Mean Gradient Norm MI treshold')

                                filename= '/Mean Gradient Norm VM_MI treshold.png'
                            else :
                                fig, axs = plt.subplots(1, 1)
                                axs.imshow(Grad_vm_treshold,interpolation='none', vmin =0 , vmax =255)
                                axs.set_title('Mean Gradient Norm VM treshold')
                                filename= '/Mean Gradient Norm VM treshold.png'

                            fig.savefig(header+newheader+filename, transparent=True)
                            plt.show(block=True)
                            plt.close()
                            fig.clf()
    # Ici on est a la fin dun dossier

    for row in range(5):
        for line in range(5):
            if row == 0 and line == 4 :
                break
            first_ellipse = False
            first_ellipse_PP = False
            first_ellipse_S2 = False
            S1_stop = False
            S2_stop = False
            PP_stop = False
            abscisse_thet = [0 for _ in Time_delay]
            Ellipses_VM_tbt = []
            Nb_elps = 0
            Nb_elps_PP = 0
            Nb_elps_S2 = 0
            Dist = [[] for k in range(ray)]
            Dist_S2 = [[] for k in range(ray)]
            Dist_PP = [[] for k in range(ray)]
            mean_vm_tbt_S1 = np.zeros((100,100))
            mean_vm_tbt_S2 = np.zeros((100,100))
            mean_vm_tbt_PP = np.zeros((100,100))

            for count, time_delay in enumerate(Time_delay) :
                mean_vm_tbt = np.array( (VM_Whisker_tbt[Mouse_index][row][line][count] - min_min_vm ) * 255/(max_max_vm-min_min_vm) , dtype = np.uint8)
                blur_vm = cv.GaussianBlur(mean_vm_tbt,(9,9),0)
                ret1,mean_vm_tbt = cv.threshold(blur_vm,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                ret1,mean_vm_tbt = cv.threshold(blur_vm,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)
                for i in range(100):
                    for j in range(100):
                        mean_vm_tbt_S1[i][j] = min(mean_vm_tbt[i][j],Mask[i][j])
                        mean_vm_tbt_S2[i][j] = min(mean_vm_tbt[i][j],MaskS2[i][j])
                        mean_vm_tbt_PP[i][j] = min(mean_vm_tbt[i][j],MaskPP[i][j])

                contours1,hierarchy1 = cv.findContours(np.array(mean_vm_tbt_S1).astype(np.uint8), 1, 2)
                if S2_PP :
                    contoursS2,hierarchyS2 = cv.findContours(np.array(mean_vm_tbt_S2).astype(np.uint8), 1, 2)
                    contoursPP,hierarchyPP = cv.findContours(np.array(mean_vm_tbt_PP).astype(np.uint8), 1, 2)


                if len(contours1) != 0:
                    cont = max(contours1,key = len)
                    if len(cont) >= 5:
                        elps = cv.fitEllipse(cont)
                        area = np.pi*elps[1][0]*elps[1][1]/4
                        if area < area_limit and not S1_stop :
                            if not first_ellipse :
                                first_ellipse =True
                                centre_S1 = (int(round(elps[0][0])),int(round(elps[0][1])))
                                theta_S1 = elps[2]
                            Nb_elps+=1
                            abscisse_thet[count]+=1
                            Ellipses_VM_tbt.append(elps)
                            if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                cv.ellipse(mean_vm_tbt, elps, color_ellipse)
                                center = (int(round(elps[0][0])),int(round(elps[0][1])))
                                axes = (int(round(elps[1][0]/2)),int(round(elps[1][1]/2)))
                                theta = elps[2]
                                cv.ellipse(VM_ellipse[injection_start+time_delay],center,axes,theta ,0. ,3.,(170, 170, 170), 1)
                                ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                                bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                                ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                                by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                                for count,angle in enumerate([k*2*np.pi/ray for k in range(ray)]):
                                    point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))))
                                    # Distance par rapport au centre de chaque ellipse
                                    # cv.line(VM_ellipse[injection_start+time_delay][:][20:80],(int(round(elps[0][0])),int(round(elps[0][1]))),point, color_ray)
                                    # Dist[count].append(np.linalg.norm([int(round(elps[0][0])),int(round(elps[0][1]))]-np.array(point)))
                                    cv.line(mean_vm_tbt,centre_S1,point, (255-255/ray*count,255-255/ray*count,255-255/ray*count))
                                    Dist[count].append(np.linalg.norm(np.array(centre_S1)-np.array(point)))
                        else :
                            Ellipses_VM_tbt.append(0)
                            if Nb_elps != 0 :
                                S1_stop = True
                    else:
                        Ellipses_VM_tbt.append(0)
                else :
                    Ellipses_VM_tbt.append(0)

                fig, axs = plt.subplots(1, 1)
                axs.imshow(mean_vm_tbt, cmap = RdBu,interpolation='none')
                axs.set_title('Mean VM tresh Aire centrale')
                head ='/Analyse_complete_Hubatz_1/Whisker_'+str(row+1)+'_'+str(line+1)+'/Mouse_'+str(Mouse)
                filename= '/Mean_VM_ellipse_'+str(injection_start+time_delay)+'_Mouse_'+str(Mouse)+'.png'
                fig.savefig(header+head+filename, transparent=True)
                plt.show(block=True)
                plt.close()
                fig.clf()

            fig = plt.figure()
            for index_1,time_delay in enumerate(Time_delay) :
                if abscisse_thet[index_1]==1 :
                    plt.scatter(time_delay,90-Ellipses_VM_tbt[index_1][2])
            plt.ylim(-90,90)
            plt.title('Evolution theta moyen selon le temps souris '+str(Mouse))
            filename= '/Evol_theta_mean_VM.png'
            fig.savefig(header+head+filename)
            plt.close()
            fig.clf()

            fig = plt.figure()
            for k in range(ray):
                plt.plot(Time_delay[:Nb_elps], Dist[k], label = 'Angle = ' + str(k*360/ray), color = cm.hsv(k/ray))
            plt.legend(loc = 'upper left')
            filename= '/Distance_mean_VM_ellipse_S1_rays_'+str(ray)+'.png'
            fig.savefig(header+head+filename)
            plt.close()
            fig.clf()

            Mean_vm_whisker[Mouse_index][row][line] = np.array([[np.mean([VM_Whisker[row][line][k][i][j] for k in range(len(VM_Whisker[row][line]))]) for j in range(window)] for i in range(window)])

    Mean_vm = np.array([[np.mean([VM_list[k][i][j] for k in range(nb_previous_records,nb_previous_records+nb_of_records)]) for j in range(window)] for i in range(window)])
    for row in range(5):
        for line in range(5):
            if row == 0 and line == 4 :
                break
            fig=plt.figure()
            plt.imshow(Mean_vm_whisker[Mouse_index][row][line])
            plt.colorbar()
            head ='/Analyse_complete_Hubatz_1/Whisker_'+str(row+1)+'_'+str(line+1)
            filename= '/Whisker_'+str(row+1)+'_'+str(line+1)+'_mouse_'+str(Mouse_index+1)+'_mean_VM_treshold='+str(treshold)+'.png'
            fig.savefig(header+head+filename)
            plt.close()
            fig.clf()



    if Gradient :
        grad_vm_smooth = np.array([[np.mean([Grad_list_VM_smooth[k][i][j] for k in range(nb_previous_records,nb_previous_records+nb_of_records)]) for j in range(window)] for i in range(window)])
        if MI :
            grad_mi_smooth = np.array([[np.mean([Grad_list_MI_smooth[k][i][j] for k in range(nb_previous_records,nb_previous_records+nb_of_records)]) for j in range(window)] for i in range(window)])





    if Ellipse :
        # Minor/Major axis VM & MI
        fig = plt.figure()
        plt.scatter([Ellipses_mean_VM[i][1][0] for i in range(len(Ellipses_mean_VM))],[Ellipses_mean_VM[i][1][1] for i in range(len(Ellipses_mean_VM))], c = [i for i in range(len(Ellipses_mean_VM))])
        plt.xlabel('Demi-grand axe')
        plt.ylabel('Demi-petit axe')
        plt.title('Repartition VM souris '+str(Mouse_index+2))
        filename= '/Petit-Grand Axe VM souris '+str(Mouse_index+2)+'.png'
        fig.savefig(header+filename)
        plt.close()
        fig.clf()
        if MI:
            fig = plt.figure()
            plt.scatter([Ellipses_mean_MI[i][1][0] for i in range(len(Ellipses_mean_MI))],[Ellipses_mean_MI[i][1][1] for i in range(len(Ellipses_mean_MI))], c = [i for i in range(len(Ellipses_mean_MI))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Demi-petit axe')
            plt.title('Repartition MI souris '+str(Mouse_index+2))
            filename= '/Petit-Grand Axe MI souris '+str(Mouse_index+2)+'.png'
            fig.savefig(header+filename)
            plt.close()
            fig.clf()

        # Major axis THETA VM & MI
        fig = plt.figure()
        plt.scatter([Ellipses_mean_VM[i][1][0] for i in range(len(Ellipses_mean_VM))],[Ellipses_mean_VM[i][-1] for i in range(len(Ellipses_mean_VM))], c = [i for i in range(len(Ellipses_mean_VM))])
        plt.xlabel('Demi-grand axe')
        plt.ylabel('Theta (degr)')
        plt.title('Repartition VM souris '+str(Mouse_index+2))
        filename= '/GrandAxe-Theta_VM souris'+str(Mouse_index+2)+'.png'
        fig.savefig(header+filename)
        plt.close()
        fig.clf()
        if MI:
            fig = plt.figure()
            plt.scatter([Ellipses_mean_MI[i][1][0] for i in range(len(Ellipses_mean_MI))],[Ellipses_mean_MI[i][-1] for i in range(len(Ellipses_mean_MI))], c = [i for i in range(len(Ellipses_mean_MI))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Theta (degr)')
            plt.title('Repartition MI souris '+str(Mouse_index+2))
            filename= '/GrandAxe-Theta_MI souris'+str(Mouse_index+2)+'.png'
            fig.savefig(header+filename)
            plt.close()

        if Gradient :
            # Minor/Major axis GRAD SMOOTh VM & MI
            fig = plt.figure()
            plt.scatter([Ellipses_grad_VM[i][1][0] for i in range(len(Ellipses_grad_VM))],[Ellipses_grad_VM[i][1][1] for i in range(len(Ellipses_grad_VM))], c = [i for i in range(len(Ellipses_grad_VM))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Demi-petit axe')
            plt.title('Repartition VM souris '+str(Mouse_index+2))
            filename= '/Minor-Major Axis Gradient VM Smooth souris'+str(Mouse_index+2)+'.png'
            fig.savefig(header+filename)
            plt.close()
            fig.clf()
            if MI :
                fig = plt.figure()
                plt.scatter([Ellipses_grad_MI[i][1][0] for i in range(len(Ellipses_grad_MI))],[Ellipses_grad_MI[i][1][1] for i in range(len(Ellipses_grad_MI))], c = [i for i in range(len(Ellipses_grad_MI))])
                plt.xlabel('Demi-grand axe')
                plt.ylabel('Demi-petit axe')
                plt.title('Repartition MI souris '+str(Mouse_index+2))
                filename= '/Minor-Major Axis Gradient MI Smooth souris'+str(Mouse_index+2)+'.png'
                fig.savefig(header+filename)
                plt.close()
                fig.clf()

            # Major axis THETA GRAD SMOOTH VM & MI
            fig = plt.figure()
            plt.scatter([Ellipses_grad_VM[i][1][0] for i in range(len(Ellipses_grad_VM))],[Ellipses_grad_VM[i][-1] for i in range(len(Ellipses_grad_VM))], c = [i for i in range(len(Ellipses_grad_VM))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Theta (degr)')
            plt.title('Repartition VM souris '+str(Mouse_index+2))
            filename= '/MajorAxis-Theta_Grad_VM souris'+str(Mouse_index+2)+'.png'
            fig.savefig(header+filename)
            plt.close()
            fig.clf()
            if MI :
                fig = plt.figure()
                plt.scatter([Ellipses_grad_MI[i][1][0] for i in range(len(Ellipses_grad_MI))],[Ellipses_grad_MI[i][-1] for i in range(len(Ellipses_grad_MI))], c = [i for i in range(len(Ellipses_grad_MI))])
                plt.xlabel('Demi-grand axe')
                plt.ylabel('Theta (degr)')
                plt.title('Repartition MI souris '+str(Mouse_index+2))
                filename= '/MajorAxis-Theta_Grad_MI souris'+str(Mouse_index+2)+'.png'
                fig.savefig(header+filename)
                plt.close()
    Mouse_index+=1

if Ellipse :

    mean_dist_foyer = [[[] for _ in range(5)] for _ in range(5)]
    Mean_thet_tbt = [[[] for _ in range(5)] for _ in range(5)]
    Var_thet_tbt = [[[] for _ in range(5)] for _ in range(5)]
    Mean_thet_tbt_tri = [[[] for _ in range(5)] for _ in range(5)]
    Var_thet_tbt_tri = [[[] for _ in range(5)] for _ in range(5)]
    for row in range(5):
        for line in range(5):
            for i in range(len(Time_delay)):
                mean_dist_foyer[row][line].append(np.mean(Dist_foyer[row][line][i]))
                Mean_thet_tbt[row][line].append(np.mean(List_thet_tbt[row][line][i]))
                Var_thet_tbt[row][line].append(np.var(List_thet_tbt[row][line][i]))
                Mean_thet_tbt_tri[row][line].append(np.mean(List_thet_tbt_tri[row][line][i]))
                Var_thet_tbt_tri[row][line].append(np.var(List_thet_tbt_tri[row][line][i]))
            axs_thet[row][line].set_ylim(-90,90)
            axs_thet[row][line].plot(Time_delay,Mean_thet_tbt[row][line], color = 'r', marker = 'o')
            axs_thet_tri[row][line].set_ylim(-90,90)
            axs_thet_tri[row][line].plot(Time_delay,Mean_thet_tbt_tri[row][line], color = 'r', marker = 'o')
    for ax in axs_thet.flat:
        ax.label_outer()
    for ax in axs_thet_tri.flat:
        ax.label_outer()
    print('The mean distance made by the focus is '+str(mean_dist_foyer))
    fig_thet.suptitle('Evolution globale des thetas')
    filename= '/Global_evol_theta.png'
    fig_thet.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
    fig_thet.clf()
    fig_thet_tri.suptitle('Evolution globale des thetas filtre selon la distance parcourue')
    filename= '/Global_evol_theta_dist_treshold.png'
    fig_thet_tri.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
    fig_thet_tri.clf()

    data = open(header+data_header +"/complete_analysis.txt", 'w')
    data.write('\n Mean Dist   = ' + str(mean_dist_foyer)+ '\n')
    data.write('\n Variance of the thetas :\n')
    data.write('Var_all = ' + str(Var_thet_tbt) + '\n' +
    'Var_tri = ' + str(Var_thet_tbt_tri))
    data.close()


    for compt, ax in enumerate(axs_w_f.flat):
        ax.set(xlabel=compt, ylabel=compt)
    fig_whisker_foyer.suptitle('Evolution des foyers')
    filename_foyer = '/Evol_foyer_all_whiskers.png'
    fig_whisker_foyer.savefig(header+'/Analyse_complete_Hubatz_1'+filename_foyer)
    plt.close()
    fig_whisker_foyer.clf()


N = len(VM_list)

Mean_vm = np.array([[np.mean([VM_list[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
fig=plt.figure()
plt.imshow(Mean_vm)
plt.colorbar()
filename= '/Global_mean_VM_treshold='+str(treshold)+'.png'
fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
plt.close()
fig.clf()
np.savetxt(header+data_header+'/Global_mean_VM_treshold='+str(treshold)+'.txt',Mean_vm)
if MI :
    Mean_mi = np.array([[np.mean([MI_list[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
    fig=plt.figure()
    plt.imshow(Mean_mi)
    plt.colorbar()
    filename= '/Global_mean_MI_treshold='+str(treshold)+'.png'
    fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
    plt.close()
    fig.clf()
    np.savetxt(header+data_header+'/Global_mean_MI_treshold='+str(treshold)+'.txt',Mean_mi)





if Gradient :
    #Gradient VM
    Mean_grad_vm = np.array([[np.mean([Grad_list_VM[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
    fig=plt.figure()
    plt.imshow(Mean_grad_vm)
    plt.clim([0,0.004])
    plt.colorbar()
    filename= '/Global_mean_gradient_norm_VM_='+str(treshold)+'.png'
    fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
    plt.close()
    fig.clf()
    np.savetxt(header+data_header+'/Global_mean_gradient_norm_VM_='+str(treshold)+'.txt',Mean_grad_vm)
    if MI :
        #Gradient MI
        Mean_grad_mi = np.array([[np.mean([Grad_list_MI[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
        fig=plt.figure()
        plt.imshow(Mean_grad_mi)
        plt.clim([0.014,0.05])
        plt.colorbar()
        filename= '/Global_mean_gradient_norm_MI_='+str(treshold)+'.png'
        fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
        plt.close()
        fig.clf()
        np.savetxt(header+data_header+'/Global_mean_gradient_norm_MI_='+str(treshold)+'.txt',Mean_grad_mi)

        grad_mi_smooth = np.array([[np.mean([Grad_list_MI_smooth[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
        fig=plt.figure()
        plt.imshow(grad_mi_smooth)
        #plt.clim([0.014,0.05])
        plt.colorbar()
        filename= '/Global_mean_gradient_norm_MI_smooth.png'
        fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
        plt.close()
        fig.clf()
        np.savetxt(header+'/Analyse_complete_Hubatz_1'+'/Global_mean_gradient_norm_MI_smooth.txt',grad_mi_smooth)


    grad_vm_smooth = np.array([[np.mean([Grad_list_VM_smooth[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
    fig=plt.figure()
    plt.imshow(grad_vm_smooth)
    #plt.clim([0,0.004])
    plt.colorbar()
    filename= '/Global_mean_gradient_norm_VM_smooth.png'
    fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
    plt.close()
    fig.clf()
    np.savetxt(header+'/Analyse_complete_Hubatz_1'+'/Global_mean_gradient_VM_smooth.txt',grad_vm_smooth)


if Ellipse :
    vm = Mean_vm
    minimum_vm = min([vm[i][j] for i in range(window) for j in range(window)])
    maximum_vm = max([vm[i][j] for i in range(window) for j in range(window)])
    if MI :
        mi = Mean_mi
        minimum_mi = min([mi[i][j] for i in range(window) for j in range(window)])
        maximum_mi = max([mi[i][j] for i in range(window) for j in range(window)])
    VM = np.array( (vm - minimum_vm ) * 255/(maximum_vm-minimum_vm) , dtype = np.uint8)
    blur_vm = cv.GaussianBlur(VM,(9,9),0)
    ret1,VM = cv.threshold(blur_vm,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
    ret1,VM = cv.threshold(blur_vm,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)
    if MI :
        MI_map = np.array( (mi - minimum_mi ) * 255/(maximum_mi-minimum_mi) , dtype = np.uint8)
        blur_mi = cv.GaussianBlur(MI_map,(9,9),0)
        ret1,MI_map = cv.threshold(blur_mi,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
        ret1,MI_map = cv.threshold(blur_mi,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)



    contours1,hierarchy1 = cv.findContours(VM, 1, 2)
    if MI :
        contours2,hierarchy2 = cv.findContours(MI_map, 1, 2)

    if len(contours1) != 0:
        cont = max(contours1,key = len)
        if len(cont) >= 5:
            elps = cv.fitEllipse(cont)
            if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                cv.ellipse(VM, elps, color_ellipse, 2 )
    if MI :
        if len(contours2) != 0:
            cont = max(contours2,key = len)
            if len(cont) >= 5:
                elps = cv.fitEllipse(cont)
                if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                    cv.ellipse(MI_map, elps, color_ellipse, 2 )
    if MI :
        fig, axs = plt.subplots(1, 2)

        axs[0].imshow(VM,interpolation='none', vmin =0 , vmax =255)
        axs[0].set_title('Global Mean VM treshold')
        axs[1].imshow(MI_list_treshold, cmap = 'inferno',interpolation='none', vmin = 0, vmax = 255)
        axs[1].set_title('Global Mean MI treshold')

        filename= '/Global Mean VM_MI treshold.png'
    else :
        fig, axs = plt.subplots(1, 1)
        axs.imshow(VM,interpolation='none', vmin =0 , vmax =255)
        axs.set_title('Global Mean VM treshold')
        filename= '/Global Mean VM treshold.png'
    fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename, transparent=True)
    plt.show(block=True)
    plt.close()
    fig.clf()
    # Minor/Major axis VM & MI
    fig = plt.figure()
    plt.scatter([Ellipses_mouse_VM[i][1][0] for i in range(len(Ellipses_mouse_VM))],[Ellipses_mouse_VM[i][1][1] for i in range(len(Ellipses_mouse_VM))], c = [i for i in range(len(Ellipses_mouse_VM))])
    plt.xlabel('Demi-grand axe')
    plt.ylabel('Demi-petit axe')
    plt.title('Repartition VM toutes souris')
    filename= '/Petit-Grand Axe VM.png'
    fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
    plt.close()
    fig.clf()
    if MI :
        fig = plt.figure()
        plt.scatter([Ellipses_mouse_MI[i][1][0] for i in range(len(Ellipses_mouse_MI))],[Ellipses_mouse_MI[i][1][1] for i in range(len(Ellipses_mouse_MI))], c = [i for i in range(len(Ellipses_mouse_MI))])
        plt.xlabel('Demi-grand axe')
        plt.ylabel('Demi-petit axe')
        plt.title('Repartition MI toutes souris ')
        filename= '/Petit-Grand Axe MI.png'
        fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
        plt.close()
        fig.clf()

    # Major axis THETA VM & MI
    fig = plt.figure()
    plt.scatter([Ellipses_mouse_VM[i][1][0] for i in range(len(Ellipses_mouse_VM))],[Ellipses_mouse_VM[i][2] for i in range(len(Ellipses_mouse_VM))], c = [i for i in range(len(Ellipses_mouse_VM))])
    plt.xlabel('Demi-grand axe')
    plt.ylabel('Theta (degr)')
    plt.title('Repartition VM souris '+str(Mouse_index))
    filename= '/GrandAxe-Theta_VM.png'
    fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
    plt.close()
    fig.clf()
    if MI :
        fig = plt.figure()
        plt.scatter([Ellipses_mouse_MI[i][1][0] for i in range(len(Ellipses_mouse_MI))],[Ellipses_mouse_MI[i][-1] for i in range(len(Ellipses_mouse_MI))], c = [i for i in range(len(Ellipses_mouse_MI))])
        plt.xlabel('Demi-grand axe')
        plt.ylabel('Theta (degr)')
        plt.title('Repartition MI toutes souris')
        filename= '/GrandAxe-Theta_MI.png'
        fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
        plt.close()


    if Gradient :
        minimum_vm = min([grad_vm_smooth[i][j] for i in range(window) for j in range(window)])
        maximum_vm = max([grad_vm_smooth[i][j] for i in range(window) for j in range(window)])
        if MI :
            minimum_mi = min([grad_mi_smooth[i][j] for i in range(window) for j in range(window)])
            maximum_mi = max([grad_mi_smooth[i][j] for i in range(window) for j in range(window)])

        VM = np.array( (grad_vm_smooth - minimum_vm ) * 255/(maximum_vm-minimum_vm) , dtype = np.uint8)
        blur_vm = cv.GaussianBlur(VM,(9,9),0)
        ret1,VM = cv.threshold(blur_vm,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
        ret1,VM = cv.threshold(blur_vm,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)
        if MI :
            MI_map = np.array( (grad_mi_smooth - minimum_mi ) * 255/(maximum_mi-minimum_mi) , dtype = np.uint8)
            blur_mi = cv.GaussianBlur(MI_map,(9,9),0)
            ret1,MI_map = cv.threshold(blur_mi,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
            ret1,MI_map = cv.threshold(blur_mi,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)


        contours1,hierarchy1 = cv.findContours(VM, 1, 2)
        if MI :
            contours2,hierarchy2 = cv.findContours(MI_map, 1, 2)
        if len(contours1) != 0:
            cont = max(contours1,key = len)
            if len(cont) >= 5:
                elps = cv.fitEllipse(cont)
                if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                    cv.ellipse(VM, elps, color_ellipse, 2 )
        if MI :
            if len(contours2) != 0:
                cont = max(contours2,key = len)
                if len(cont) >= 5:
                    elps = cv.fitEllipse(cont)
                    if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                        cv.ellipse(MI_map, elps, color_ellipse, 2 )
        if MI :
            fig, axs = plt.subplots(1, 2)

            axs[0].imshow(VM,interpolation='none', vmin =0 , vmax =255)
            axs[0].set_title('Global Mean Gradient Norm VM treshold')
            axs[1].imshow(MI_list_treshold, cmap = 'inferno',interpolation='none', vmin = 0, vmax = 255)
            axs[1].set_title('Global Mean  Gradient Norm MI treshold')
            filename= '/Global Mean Gradient Norm VM_MI treshold.png'
        else :
            fig, axs = plt.subplots(1, 1)
            axs.imshow(VM,interpolation='none', vmin =0 , vmax =255)
            axs.set_title('Global Mean Gradient Norm VM treshold')
            filename= '/Global Mean Gradient Norm VM treshold.png'
        fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename, transparent=True)
        plt.show(block=True)
        plt.close()
        fig.clf()

                # Minor/Major axis GRAD SMOOTh VM & MI
        fig = plt.figure()
        plt.scatter([Ellipses_mouse_grad_VM[i][1][0] for i in range(len(Ellipses_mouse_grad_VM))],[Ellipses_mouse_grad_VM[i][1][1] for i in range(len(Ellipses_mouse_grad_VM))], c = [i for i in range(len(Ellipses_mouse_grad_MI))])
        plt.xlabel('Demi-grand axe')
        plt.ylabel('Demi-petit axe')
        plt.title('Repartition VM toutes souris ')
        filename= '/Minor-Major Axis Gradient VM Smooth.png'
        fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
        plt.close()
        fig.clf()
        if MI :
            fig = plt.figure()
            plt.scatter([Ellipses_mouse_grad_MI[i][1][0] for i in range(len(Ellipses_mouse_grad_MI))],[Ellipses_mouse_grad_MI[i][1][1] for i in range(len(Ellipses_mouse_grad_MI))], c = [i for i in range(len(Ellipses_mouse_grad_MI))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Demi-petit axe')
            plt.title('Repartition MI totues souris ')
            filename= '/Minor-Major Axis Gradient MI Smooth.png'
            fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
            plt.close()
            fig.clf()

        # Major axis THETA GRAD SMOOTH VM & MI
        fig = plt.figure()
        plt.scatter([Ellipses_mouse_grad_VM[i][1][0] for i in range(len(Ellipses_mouse_grad_VM))],[Ellipses_mouse_grad_VM[i][-1] for i in range(len(Ellipses_mouse_grad_VM))], c = [i for i in range(len(Ellipses_mouse_grad_VM))])
        plt.xlabel('Demi-grand axe')
        plt.ylabel('Theta (degr)')
        plt.title('Repartition VM toutes souris ')
        filename= '/MajorAxis-Theta_Grad_VM.png'
        fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
        plt.close()
        fig.clf()
        if MI :
            fig = plt.figure()
            plt.scatter([Ellipses_mouse_grad_MI[i][1][0] for i in range(len(Ellipses_mouse_grad_MI))],[Ellipses_mouse_grad_MI[i][-1] for i in range(len(Ellipses_mouse_grad_MI))], c = [i for i in range(len(Ellipses_mouse_grad_MI))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Theta (degr)')
            plt.title('Repartition MI toutes souris ')
            filename= '/MajorAxis-Theta_Grad_MI.png'
            fig.savefig(header+'/Analyse_complete_Hubatz_1'+filename)
            plt.close()
for row in range(5):
    for line in range(5):
        if row == 0 and line == 4 :
            break
        fig=plt.figure()
        plt.imshow([[np.mean([Mean_vm_whisker[Mouse_index][row][line][i][j] for Mouse_index in range(4)]) for j in range(window)] for i in range(window)])
        plt.colorbar()
        head ='/Analyse_complete_Hubatz_1/'
        filename= '/Whisker_'+str(row+1)+'_'+str(line+1)+'_all mice_-mean_VM_treshold='+str(treshold)+'.png'
        fig.savefig(header+head+filename)
        plt.close()
        fig.clf()
