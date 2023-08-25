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


Mean_VM = True
MI = False
Annulus = False
Gradient = False
Ellipse = True
S2 = False
PP = False
treshold = 0.0015 #Treshold to filter the data
# Recommended : 0.0015

Time_delay = np.arange(7,16,1)
run_time = 511 # nombre d'images
size = 100 # taille de l'image
injection_start,injection_end = 100,511 #début de la stimulation
interval = 1 # intervalle sur lequel on calcule la MI


for rise in [1.01,1.015,1.02,1.3,1.35,1.4]:
        
    header='/home/margauxvrech/mutual_network/DATAVSD/'
    data_header='/Analyse_complete_2/Data'
    if not os.path.exists(header+data_header):
        os.makedirs(header+data_header)
    print('\nNumeric data will be saved in '+data_header)

    VM_list = [] #VM
    VM_each_mice = [[[[0 for k in range(100)] for l in range(100)] for t in range(511)] for _ in range(4)]
    if MI :
        MI_list = [] #MI


    if Annulus :
        ref_neurone = [50,50] #neurone de référence
        number_of_annulus = 21


    x = 0
    y = 0
    window = 100
    list_coord = [(x+i)*size+(y+j) for j in range(window) for i in range(window)]

    VM_list_smooth = []
    if MI :
        MI_list_smooth = []

    if Gradient :
        Grad_list_VM = [] #Gradient de la "VM"
        Grad_list_VM_smooth = []
        if MI :
            Grad_list_MI = [] #Gradient de la MI
            Grad_list_MI_smooth = []



    def eq_S1S2(x):
        return 20
    def eq_S1PP(x):
        return 0.35*x+60
    def eq_S1PP_bis(x):
        return -0.66*x+126

    if Ellipse :
        rise = 1.025
        ray = 10
        area_limit = 6000
        area_PP_min = 400
        area_PP_max = 1000000000
        area_S2_min = 200
        area_S2_max = 1000000000
        color_ellipse = (100,100,100)
        Ellipses_mouse_VM = []
        Mask = np.zeros((100,100))
        MaskS2 = np.zeros((100,100))
        MaskPP = np.zeros((100,100))
        List_thet_tbt = [[] for _ in Time_delay]
        List_thet_tbt_tri = [[] for _ in Time_delay]
        Dist_foyer = [[] for _ in Time_delay]
        mean_dist_foyer_tbt= [12, 1, 1, 1, 1, 2,2, 4, 8]
        # Recommended : [12, 1, 1, 1, 1, 2,2, 4, 8] Must be same lenght as Time_delay

        for x in range(100):
            for y in range(100):
                if eq_S1PP(x) < y or eq_S1PP_bis(x) < y:
                    MaskPP[y,x] = 255
                elif eq_S1PP(x) > y > eq_S1S2(x) and y < eq_S1PP_bis(x):
                    Mask[y,x] = 255
                else :
                    MaskS2[y,x] = 255

        fig=plt.figure()
        plt.imshow(Mask+0.5*MaskS2+0.2*MaskPP)
        filename= '/Mask'+'.png'
        fig.savefig('/home/margauxvrech/mutual_network/DATAVSD/Analyse_complete_2/'+filename, transparent=True)
        plt.close()
        fig.clf()

        if MI :
             Ellipses_mouse_MI = []
        if Gradient :
            Ellipses_mouse_grad_VM = []
            if MI :
                Ellipses_mouse_grad_MI = []







    data = open(header+data_header +"/complete_analysis.txt", 'w')
    data.write('Complete analysis\n Analysis options :\n')
    data.write('Mean_VM = ' + str(Mean_VM) + '\n' +
    'MI = ' + str(MI) + '\n' +
    'Annulus = ' + str(Annulus) + '\n' +
    'Gradient = ' + str(Gradient) + '\n' +
    'Ellipse = ' + str(Ellipse) + '\n' +
    'treshold = ' + str(treshold) + '\n' +
    'S2 = ' + str(S2) + '\n' +
    'PP = ' + str(PP) + '\n')
    if Annulus :
        data.write('Reference neuron : ' + str(ref_neurone) + '\n' +
        'Number of annulus = ' + str(number_of_annulus) + '\n')
    if Ellipse :
        data.write('Number of rays : ' + str(ray) + '\n')
    data.close()

    Mouse_index = 0
    mouse_list = [3,4,4,5,6]

    nb_previous_records = 0
    nb_of_records = 0
    fig_thet, axs_thet = plt.subplots(1, 1)
    fig_thet_tri, axs_thet_tri = plt.subplots(1, 1)


    #
    for folder in ['Mouses_3-5-6/20160912','Single_C2_whisker_evoked_responses_first_10_trials','Single C2 whisker evoked responses, second set of 10 trials','Mouses_3-5-6/20160914/','Mouses_3-5-6/20160916/']:
        Mouse = mouse_list[Mouse_index]
        nb_previous_records = nb_of_records+nb_previous_records
        nb_of_records = 0
        mouse_header = '/Analyse_complete_2/Mouse_'+str(Mouse)
        if Ellipse :


            Ellipses_mean_VM = []
            if MI :
                Ellipses_mean_MI = []
            if Gradient :
                Ellipses_grad_VM = []
                if MI :
                    Ellipses_grad_MI = []

        for a in range(1,21):
            if folder[0]=='M':
                file = '/C2_'+str(a)+'.txt'
                L=[[[0 for k in range(100)] for l in range(100)]for i in range(511)]
                with open(folder+file, 'r') as fr:
                    lines = fr.readlines()
                    for t in range(511):
                        for i in range(100):
                            L[t][i]=[float(lines[t*100+j].split('\t')[i]) for j in range(100)]
                file = '/C2_'+str(a)
            elif folder[36]!='s':
                if a>=11:
                    break
                else:
                    file = '/C2_'+str(a)
                    L=[[[0 for k in range(100)] for l in range(100)]for i in range(511)]
                    with open(folder+file, 'r') as fr:
                        lines = fr.readlines()
                        for t in range(511):
                            for i in range(100):
                                for j in range(100):
                                    L[t][-i+99][j]=float(lines[t*10000+i*100+j])
            else :
                if a>=11:
                    break
                else:
                    file = '/C2_'+str(10+a)
                    L=[[[0 for k in range(100)] for l in range(100)]for i in range(511)]
                    with open(folder+file, 'r') as fr:
                        lines = fr.readlines()
                        for t in range(511):
                            for i in range(100):
                                for j in range(100):
                                    L[t][-i+99][j]=float(lines[t*10000+i*100+j])


            ######################################################
            # Création d'un dossier associé à chaque souris
            ######################################################

            newheader='/Analyse_complete_2/Mouse_'+str(Mouse)+file
            if not os.path.exists(header+newheader):
                os.makedirs(header+newheader)
            print('\nAnalysing data of '+file)
            print('\nData will be saved in '+newheader)

            if not os.path.exists(header+data_header+'/Mouse_'+str(Mouse)+file):
                os.makedirs(header+data_header+'/Mouse_'+str(Mouse)+file)

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


            if max_mean_vm>=treshold:
                print("We analyse this data\n")
                nb_of_records+=1


                max_vm,min_vm = -1000,1000
                max_mi,min_mi = -1000,1000

                if Annulus :
                    vm_base_brut=[L[t][ref_neurone[0]][ref_neurone[1]] for t in range(injection_start,injection_start+interval)]

                for time_delay in Time_delay:
                    for i in range(window):
                        for j in range(window):

                            vm_neurone_brut = [L[t][i][j] for t in range(injection_start+time_delay,injection_start+interval+time_delay)]
                            a = np.mean([L[t][i][j] for t in range(injection_start+time_delay,injection_start+interval+time_delay)])
                            VM[injection_start+time_delay][i][j] = a
                            VM_each_mice[Mouse-3][injection_start+time_delay][i][j] += a/20

                            if a < min_vm :
                                min_vm = a
                            if a > max_vm :
                                max_vm = a
                            if MI :
                                c_Y_brut,xedges = np.histogram(vm_neurone_brut,200,range=(-0.1,0.02))
                                Recorded_cell_brut[injection_start+time_delay][i][j]= mutual_info_score(c_X_brut,c_Y_brut)

                                if Recorded_cell_brut[injection_start+time_delay][i][j] < min_mi :
                                    min_mi = Recorded_cell_brut[injection_start+time_delay][i][j]
                                if Recorded_cell_brut[injection_start+time_delay][i][j] > max_mi :
                                    max_mi = Recorded_cell_brut[injection_start+time_delay][i][j]

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
                thet = 0
                first_thet = True
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
                        Dist_foyer[index].append(0)
                        VM_ellipse[injection_start+time_delay] = np.array( (VM[injection_start+time_delay] - min_vm ) * 255/(max_vm-min_vm) , dtype = np.uint8)
                        blur_vm = cv.GaussianBlur(VM_ellipse[injection_start+time_delay],(9,9),0)
                        ret1,VM_ellipse[injection_start+time_delay] = cv.threshold(blur_vm,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                        ret1,VM_ellipse[injection_start+time_delay] = cv.threshold(blur_vm,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)

                        np.savetxt(header+newheader+'/VM'+str(time_delay)+'.txt',VM_ellipse[injection_start+time_delay])
                        if MI :
                            MI_ellipse[injection_start+time_delay] = np.array( (Recorded_cell_brut[injection_start+time_delay] - min_MI ) * 255/(max_MI-min_MI) , dtype = np.uint8)
                            blur_mi = cv.GaussianBlur(MI_ellipse[injection_start+time_delay],(9,9),0)
                            #ret2,MI_smooth_treshold = cv.threshold(blur_mi,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                            ret2,MI_ellipse[injection_start+time_delay] = cv.threshold(blur_mi,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
                            ret2,MI_ellipse = cv.threshold(blur_mi,int(np.floor(ret2*1.005)),255,cv.THRESH_BINARY)
                            np.savetxt(header+newheader+'/MI'+str(time_delay)+'.txt',MI_ellipse[injection_start+time_delay])


                        for i in range(100):
                            for j in range(100):
                                VM_ellipse_S1[injection_start+time_delay][i][j] = min(VM_ellipse[injection_start+time_delay][i][j],Mask[i][j])
                                VM_ellipse_S2[injection_start+time_delay][i][j] = min(VM_ellipse[injection_start+time_delay][i][j],MaskS2[i][j])
                                VM_ellipse_PP[injection_start+time_delay][i][j] = min(VM_ellipse[injection_start+time_delay][i][j],MaskPP[i][j])


                        contours1,hierarchy1 = cv.findContours(np.array(VM_ellipse_S1[injection_start+time_delay]).astype(np.uint8), 1, 2)
                        if S2 :
                            contoursS2,hierarchyS2 = cv.findContours(np.array(VM_ellipse_S2[injection_start+time_delay]).astype(np.uint8), 1, 2)
                        if PP :
                            contoursPP,hierarchyPP = cv.findContours(np.array(VM_ellipse_PP[injection_start+time_delay]).astype(np.uint8), 1, 2)
                        if MI :
                            contours2,hierarchy2 = cv.findContours(MI_ellipse[injection_start+time_delay][:][20:70], 1, 2)
                        if len(contours1) != 0:
                            cont = max(contours1,key = len)
                            if len(cont) >= 5:
                                elps = cv.fitEllipse(cont)
                                area = np.pi*elps[1][0]*elps[1][1]/4
                                ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                                bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                                ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                                by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                                out = 0
                                for t in range(360):
                                    ord = elps[0][1]+ay*np.cos(t*2*np.pi/360)+by * np.sin(t*2*np.pi/360)
                                    if ord > 70 or ord < 20 :
                                        out+=1
                                if out <200  and 70>elps[0][1]>20 : #  and (np.abs(elps[2]-thet) < 30 or np.abs(elps[2]-thet-360) < 30 or first_thet )
                                    first_thet = False
                                    thet = elps[2]
                                    if area < area_limit and not S1_stop:
                                        if not first_ellipse :
                                            first_ellipse =True
                                            centre_S1 = (int(round(elps[0][0])),int(round(elps[0][1])))
                                            theta_S1 = elps[2]
                                        Nb_elps+=1
                                        Ellipses_VM.append(elps)
                                        abscisse_thet[index]+=1
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
                                                # cv.line(VM_ellipse[injection_start+time_delay][:][20:70],(int(round(elps[0][0])),int(round(elps[0][1]))),point, color_ray)
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
                        else :
                            Ellipses_VM.append(0)
                        if S2 :
                            if len(contoursS2) != 0:
                                cont = max(contoursS2,key = len)
                                if len(cont) >= 5:
                                    elps = cv.fitEllipse(cont)
                                    area_S2 = np.pi*elps[1][0]*elps[1][1]/4
                                    ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                                    bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                                    ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                                    by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                                    out = 0
                                    for t in range(360):
                                        ord = elps[0][1]+ay*np.cos(t*2*np.pi/360)+by * np.sin(t*2*np.pi/360)
                                        if ord > 20 :
                                            out+=1
                                    if  out <= 100 and -10 < elps[0][1]< 20:
                                        if area_S2_min < area_S2 < area_S2_max and not S2_stop :
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

                                                for count,angle in enumerate([k*2*np.pi/ray for k in range(ray)]):
                                                    point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_S2)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_S2)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_S2)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_S2)*2*np.pi/360))))
                                                    cv.line(VM_ellipse[injection_start+time_delay],centre_S2,point, (255-255/ray*count,255-255/ray*count,255-255/ray*count))
                                                    Dist_S2[count].append(np.linalg.norm(np.array(centre_S2)-np.array(point)))
                                        else :
                                            if Nb_elps_S2 != 0 :
                                                S2_stop = True
                        if PP :
                            if len(contoursPP) != 0:
                                cont = max(contoursPP,key = len)
                                if len(cont) >= 5:
                                    elps = cv.fitEllipse(cont)
                                    area_PP = np.pi*elps[1][0]*elps[1][1]/4
                                    ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                                    bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                                    ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                                    by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                                    out = 0
                                    for t in range(360):
                                        ord = elps[0][1]+ay*np.cos(t*2*np.pi/360)+by * np.sin(t*2*np.pi/360)
                                        if ord < 70 :
                                            out+=1
                                    if out < 250 and 100 > elps[0][1]>60:
                                        if area_PP_min < area_PP < area_PP_max and not PP_stop :
                                            if not first_ellipse_PP :
                                                first_ellipse_PP = True
                                                centre_PP = (int(round(elps[0][0])),int(round(elps[0][1])))
                                                theta_PP = elps[2]*2*np.pi/360
                                            Nb_elps_PP+=1
                                            Ellipses_VM_PP.append(elps)
                                            if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                                cv.ellipse(VM_ellipse[injection_start+time_delay], elps, color_ellipse)
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

                                        for angle in range(360):
                                            cv.ellipse(MI_ellipse[injection_start+time_delay][:][20:70], elps,angle,angle+1, color_ellipse, 2 )
                            fig, axs = plt.subplots(1, 2)

                            axs[0].imshow(VM_ellipse[injection_start+time_delay], cmap = RdBu,interpolation='none')
                            axs[0].set_title('VM tresh Aire centrale = '+str(round(area)) + ' PP ' + str(round(area_PP))+ ' S2 ' + str(round(area_S2)))
                            axs[1].imshow(MI_ellipse[injection_start+time_delay], cmap = 'inferno',interpolation='none')
                            axs[1].set_title('MI treshold')


                            filename= '/VM_MI_ellipse_'+str(injection_start+time_delay)+'.png'
                        else :
                            fig, axs = plt.subplots(1, 1)

                            axs.imshow(VM_ellipse[injection_start+time_delay], cmap = RdBu,interpolation='none')

                            axs.set_title('VM tresh')
                            filename= '/VM_ellipse_'+str(injection_start+time_delay)+'.png'
                        fig.savefig(header+newheader+filename, transparent=True)
                        plt.show(block=True)
                        plt.close()
                        fig.clf()

                for index_1,time_delay in enumerate(Time_delay) :
                    if abscisse_thet[index_1]==1 :
                        ind = 1
                        while Ellipses_VM[index_1-ind] == 0 :
                            ind += 1
                        Dist_foyer[index_1][-1] += np.linalg.norm(np.array([Ellipses_VM[index_1][0][0],Ellipses_VM[index_1][0][1]])-np.array([Ellipses_VM[index_1-ind][0][0],Ellipses_VM[index_1-ind][0][1]]))

                fig,axs = plt.subplots(1, 1)
                fig_foyer,axs_foyer = plt.subplots(1, 1)
                for index_1,time_delay in enumerate(Time_delay) :
                    if abscisse_thet[index_1]==1 :
                        List_thet_tbt[index_1].append(90-Ellipses_VM[index_1][2])
                        axs.scatter(time_delay,90-Ellipses_VM[index_1][2])
                        axs_foyer.scatter(Ellipses_VM[index_1][0][0],Ellipses_VM[index_1][0][1], color = cm.hsv(index_1/len(Time_delay)))
                        if Dist_foyer[index_1][-1] < mean_dist_foyer_tbt[index_1] :
                            axs_thet.scatter(time_delay,90-Ellipses_VM[index_1][2], marker = 'x', color = 'b')
                            axs_thet_tri.scatter(time_delay,90-Ellipses_VM[index_1][2], marker = 'x', color = 'b')
                            List_thet_tbt_tri[index_1].append(90-Ellipses_VM[index_1][2])
                        else :
                            axs_thet.scatter(time_delay,90-Ellipses_VM[index_1][2], marker = 'x', color = 'g')

                axs.set_ylim(-90,90)
                axs.set_title('Evolution theta')
                axs_foyer.set_xlim(0,100)
                axs_foyer.set_ylim(0,100)
                axs_foyer.set_title('Evolution position des foyers: distance totale = ' + str(np.sum(Dist_foyer[-1])))
                filename= '/Evol_theta.png'
                filename_foyer = '/Evol_foyer.png'
                fig.savefig(header+newheader+filename)
                fig_foyer.savefig(header+newheader+filename_foyer)
                plt.close()
                fig.clf()
                fig_foyer.clf()


                VM_list.append([[np.mean([VM[int(injection_start+time_delay)][i][j] for time_delay in Time_delay]) for j in range(window)] for i in range(window)])
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
                np.savetxt(header+data_header+'/Mouse_'+str(Mouse)+file+'/Mean_VM.txt',vm)
                if MI :
                    mi = np.array([[MI_list[-1][i][j] for j in range(window)] for i in range(window)])
                    fig = plt.figure()
                    plt.imshow(mi, vmin =0 , vmax =0.008 )
                    plt.colorbar()
                    filename= '/Mean_MI.png'
                    fig.savefig(header+newheader+filename)
                    plt.close()
                    fig.clf()
                    np.savetxt(header+data_header+'/Mouse_'+str(Mouse)+file+'/Mean_MI.txt',mi)

                any_ellipse = False
                if Ellipse :
                    for elps in Ellipses_VM :
                        if elps != 0 :
                            cv.ellipse(vm, elps, color_ellipse, 1 )
                            any_ellipse = True
                    if any_ellipse :
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
                            # point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_S2)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_S2)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_S2)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_S2)*2*np.pi/360))))
                            # cv.line(vm[:][20:70],(int(round(elps[0][0])),int(round(elps[0][1]))),point, color_ray)
                            point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_S1)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_S1)*2*np.pi/360))))
                            cv.line(vm,centre_S1,point, (0.009/ray*(ray-count),0.009/ray*(ray-count),0.009/ray*(ray-count),1))

                    fig = plt.figure()
                    plt.imshow(vm, vmin =0 , vmax =0.009)
                    plt.colorbar()
                    filename= '/Mean_VM_ellipses.png'
                    fig.savefig(header+newheader+filename)
                    plt.close()
                    fig.clf()
                    if MI :
                        for elps in Ellipses_MI :
                            cv.ellipse(mi[:][20:70], elps, color_ellipse, 1 )
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
                    if S2 :
                        fig = plt.figure()
                        for k in range(ray):
                            plt.plot(Time_delay[:Nb_elps_S2], Dist_S2[k], label = 'Angle = ' + str(k*360/ray), color = cm.hsv(k/ray))
                        plt.legend(loc = 'upper left')
                        filename= '/Distance_VM_ellipse_S2_rays_'+str(ray)+'.png'
                        fig.savefig(header+newheader+filename)
                        plt.close()
                        fig.clf()
                    if PP :
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
                    np.savetxt(header+data_header+'/Mouse_'+str(Mouse)+file+'/Mean_gradient_Norm_VM.txt',grad_vm)
                    if MI :
                        grad_mi = np.array([[Grad_list_MI[-1][i][j] for j in range(window)] for i in range(window)])
                        fig = plt.figure()
                        plt.imshow(grad_mi)#, vmin =0 , vmax =0.03 )
                        plt.colorbar()
                        filename= '/Mean_gradient_Norm_MI.png'
                        fig.savefig(header+newheader+filename)
                        plt.close()
                        fig.clf()
                        np.savetxt(header+data_header+'/Mouse_'+str(Mouse)+file+'/Mean_gradient_Norm_MI.txt',grad_mi)

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
        Mean_vm = np.array([[np.mean([VM_list[k][i][j] for k in range(nb_previous_records,nb_previous_records+nb_of_records)]) for j in range(window)] for i in range(window)])
        fig=plt.figure()
        plt.imshow(Mean_vm)
        plt.colorbar()
        filename= '/Mouse_mean_VM_treshold='+str(treshold)+'.png'
        fig.savefig(header+mouse_header+filename)
        plt.close()
        fig.clf()
        np.savetxt(header+mouse_header+'/Mouse_mean_VM_treshold='+str(treshold)+'.txt',Mean_vm)
        if MI :
            Mean_mi = np.array([[np.mean([MI_list[k][i][j] for k in range(nb_previous_records,nb_previous_records+nb_of_records)]) for j in range(window)] for i in range(window)])
            fig=plt.figure()
            plt.imshow(Mean_mi)
            plt.colorbar()
            filename= '/Mouse_mean_MI_treshold='+str(treshold)+'.png'
            fig.savefig(header+mouse_header+filename)
            plt.close()
            fig.clf()
            np.savetxt(header+mouse_header+'/Mouse_mean_MI_treshold='+str(treshold)+'.txt',Mean_mi)

        if Gradient :
            if MI :
                grad_mi_smooth = np.array([[np.mean([Grad_list_MI_smooth[k][i][j] for k in range(nb_previous_records,nb_previous_records+nb_of_records)]) for j in range(window)] for i in range(window)])
                fig=plt.figure()
                plt.imshow(grad_mi_smooth)
                #plt.clim([0.014,0.05])
                plt.colorbar()
                filename= '/Mouse_mean_gradient_norm_MI_smooth.png'
                fig.savefig(header+mouse_header+filename)
                plt.close()
                fig.clf()
                np.savetxt(header+mouse_header+'/Mouse_mean_gradient_norm_MI_smooth.txt',grad_mi_smooth)


            grad_vm_smooth = np.array([[np.mean([Grad_list_VM_smooth[k][i][j] for k in range(nb_previous_records,nb_previous_records+nb_of_records)]) for j in range(window)] for i in range(window)])
            fig=plt.figure()
            plt.imshow(grad_vm_smooth)
            #plt.clim([0,0.004])
            plt.colorbar()
            filename= '/Mouse_mean_gradient_norm_VM_smooth.png'
            fig.savefig(header+mouse_header+filename)
            plt.close()
            fig.clf()
            np.savetxt(header+mouse_header+'/Mouse_mean_gradient_VM_smooth.txt',grad_vm_smooth)


        if Ellipse :
            vm = Mean_vm
            if MI :
                mi = Mean_mi
                minimum_mi = min([mi[i][j] for i in range(window) for j in range(window)])
                maximum_mi = max([mi[i][j] for i in range(window) for j in range(window)])
            minimum_vm = min([vm[i][j] for i in range(window) for j in range(window)])
            maximum_vm = max([vm[i][j] for i in range(window) for j in range(window)])
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
                    Ellipses_mouse_VM.append(elps)
                    if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                        cv.ellipse(VM, elps, color_ellipse, 2 )
            if MI :
                if len(contours2) != 0:
                    cont = max(contours2,key = len)
                    if len(cont) >= 5:
                        elps = cv.fitEllipse(cont)
                        Ellipses_mouse_MI.append(elps)
                        if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                            cv.ellipse(MI_map, elps, color_ellipse, 2 )
            if MI :
                fig, axs = plt.subplots(1, 2)

                axs[0].imshow(VM,interpolation='none', vmin =0 , vmax =255)
                axs[0].set_title('Mouse Mean VM treshold')
                axs[1].imshow(MI_list_treshold, cmap = 'inferno',interpolation='none', vmin = 0, vmax = 255)
                axs[1].set_title('Mouse Mean MI treshold')

                filename= '/Mouse Mean VM_MI treshold.png'
            else :
                fig, axs = plt.subplots(1, 1)
                axs.imshow(VM,interpolation='none', vmin =0 , vmax =255)
                axs.set_title('Mouse Mean VM treshold')

                filename= '/Mouse Mean VM treshold.png'
            fig.savefig(header+mouse_header+filename, transparent=True)
            plt.show(block=True)
            plt.close()
            fig.clf()

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
                        Ellipses_mouse_grad_VM.append(elps)
                        if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                            cv.ellipse(VM, elps, color_ellipse, 2 )
                if MI :
                    if len(contours2) != 0:
                        cont = max(contours2,key = len)
                        if len(cont) >= 5:
                            elps = cv.fitEllipse(cont)
                            Ellipses_mouse_grad_MI.append(elps)
                            if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                cv.ellipse(MI_map, elps, color_ellipse, 2 )
                if MI :
                    fig, axs = plt.subplots(1, 2)

                    axs[0].imshow(VM,interpolation='none', vmin =0 , vmax =255)
                    axs[0].set_title('Mouse Mean Gradient Norm VM treshold')
                    axs[1].imshow(MI_list_treshold, cmap = 'inferno',interpolation='none', vmin = 0, vmax = 255)
                    axs[1].set_title('Mouse Mean  Gradient Norm MI treshold')

                    filename= '/Mouse Mean Gradient Norm VM_MI treshold.png'
                else :
                    fig, axs = plt.subplots(1, 1)

                    axs.imshow(VM,interpolation='none', vmin =0 , vmax =255)
                    axs.set_title('Mouse Mean Gradient Norm VM treshold')
                    filename= '/Mouse Mean Gradient Norm VM treshold.png'
                fig.savefig(header+mouse_header+filename, transparent=True)
                plt.show(block=True)
                plt.close()
                fig.clf()




            # Minor/Major axis VM & MI
            fig = plt.figure()
            plt.scatter([Ellipses_mean_VM[i][1][0] for i in range(len(Ellipses_mean_VM))],[Ellipses_mean_VM[i][1][1] for i in range(len(Ellipses_mean_VM))], c = [i for i in range(len(Ellipses_mean_VM))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Demi-petit axe')
            plt.title('Repartition VM souris '+str(Mouse_index+2))
            filename= '/Petit-Grand Axe VM.png'
            fig.savefig(header+mouse_header+filename)
            plt.close()
            fig.clf()
            if MI:
                fig = plt.figure()
                plt.scatter([Ellipses_mean_MI[i][1][0] for i in range(len(Ellipses_mean_MI))],[Ellipses_mean_MI[i][1][1] for i in range(len(Ellipses_mean_MI))], c = [i for i in range(len(Ellipses_mean_MI))])
                plt.xlabel('Demi-grand axe')
                plt.ylabel('Demi-petit axe')
                plt.title('Repartition MI souris '+str(Mouse_index+2))
                filename= '/Petit-Grand Axe MI.png'
                fig.savefig(header+mouse_header+filename)
                plt.close()
                fig.clf()

            # Major axis THETA VM & MI
            fig = plt.figure()
            plt.scatter([Ellipses_mean_VM[i][1][0] for i in range(len(Ellipses_mean_VM))],[Ellipses_mean_VM[i][-1] for i in range(len(Ellipses_mean_VM))], c = [i for i in range(len(Ellipses_mean_VM))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Theta (degr)')
            plt.title('Repartition VM souris '+str(Mouse_index+2))
            filename= '/GrandAxe-Theta_VM.png'
            fig.savefig(header+mouse_header+filename)
            plt.close()
            fig.clf()
            if MI:
                fig = plt.figure()
                plt.scatter([Ellipses_mean_MI[i][1][0] for i in range(len(Ellipses_mean_MI))],[Ellipses_mean_MI[i][-1] for i in range(len(Ellipses_mean_MI))], c = [i for i in range(len(Ellipses_mean_MI))])
                plt.xlabel('Demi-grand axe')
                plt.ylabel('Theta (degr)')
                plt.title('Repartition MI souris '+str(Mouse_index+2))
                filename= '/GrandAxe-Theta_MI.png'
                fig.savefig(header+mouse_header+filename)
                plt.close()

            if Gradient :
                # Minor/Major axis GRAD SMOOTh VM & MI
                fig = plt.figure()
                plt.scatter([Ellipses_grad_VM[i][1][0] for i in range(len(Ellipses_grad_VM))],[Ellipses_grad_VM[i][1][1] for i in range(len(Ellipses_grad_VM))], c = [i for i in range(len(Ellipses_grad_VM))])
                plt.xlabel('Demi-grand axe')
                plt.ylabel('Demi-petit axe')
                plt.title('Repartition VM souris '+str(Mouse_index+2))
                filename= '/Minor-Major Axis Gradient VM Smooth.png'
                fig.savefig(header+mouse_header+filename)
                plt.close()
                fig.clf()
                if MI :
                    fig = plt.figure()
                    plt.scatter([Ellipses_grad_MI[i][1][0] for i in range(len(Ellipses_grad_MI))],[Ellipses_grad_MI[i][1][1] for i in range(len(Ellipses_grad_MI))], c = [i for i in range(len(Ellipses_grad_MI))])
                    plt.xlabel('Demi-grand axe')
                    plt.ylabel('Demi-petit axe')
                    plt.title('Repartition MI souris '+str(Mouse_index+2))
                    filename= '/Minor-Major Axis Gradient MI Smooth.png'
                    fig.savefig(header+mouse_header+filename)
                    plt.close()
                    fig.clf()

                # Major axis THETA GRAD SMOOTH VM & MI
                fig = plt.figure()
                plt.scatter([Ellipses_grad_VM[i][1][0] for i in range(len(Ellipses_grad_VM))],[Ellipses_grad_VM[i][-1] for i in range(len(Ellipses_grad_VM))], c = [i for i in range(len(Ellipses_grad_VM))])
                plt.xlabel('Demi-grand axe')
                plt.ylabel('Theta (degr)')
                plt.title('Repartition VM souris '+str(Mouse_index+2))
                filename= '/MajorAxis-Theta_Grad_VM.png'
                fig.savefig(header+mouse_header+filename)
                plt.close()
                fig.clf()
                if MI :
                    fig = plt.figure()
                    plt.scatter([Ellipses_grad_MI[i][1][0] for i in range(len(Ellipses_grad_MI))],[Ellipses_grad_MI[i][-1] for i in range(len(Ellipses_grad_MI))], c = [i for i in range(len(Ellipses_grad_MI))])
                    plt.xlabel('Demi-grand axe')
                    plt.ylabel('Theta (degr)')
                    plt.title('Repartition MI souris '+str(Mouse_index+2))
                    filename= '/MajorAxis-Theta_Grad_MI.png'
                    fig.savefig(header+mouse_header+filename)
                    plt.close()




        Mouse_index+=1

    if Ellipse :
        mean_dist_foyer = []
        Mean_thet_tbt = []
        Var_thet_tbt = []
        Mean_thet_tbt_tri = []
        Var_thet_tbt_tri = []
        for i in range(len(Time_delay)):
            mean_dist_foyer.append(np.mean(Dist_foyer[i]))
            Mean_thet_tbt.append(np.mean(List_thet_tbt[i]))
            Var_thet_tbt.append(np.var(List_thet_tbt[i]))
            Mean_thet_tbt_tri.append(np.mean(List_thet_tbt_tri[i]))
            Var_thet_tbt_tri.append(np.var(List_thet_tbt_tri[i]))
        print('The mean distance made by the focus is '+str(mean_dist_foyer))
        axs_thet.set_ylim(-90,90)
        axs_thet.plot(Time_delay,Mean_thet_tbt, color = 'r', marker = 'o')
        fig_thet.suptitle('Evolution globale des thetas')
        filename= '/Global_evol_theta'+str(rise)+'.png'
        fig_thet.savefig(header+'/Analyse_complete_2'+filename)
        fig_thet.clf()
        axs_thet_tri.set_ylim(-90,90)
        axs_thet_tri.plot(Time_delay,Mean_thet_tbt_tri, color = 'r', marker = 'o')
        fig_thet_tri.suptitle('Evolution globale des thetas filtre selon la distance parcourue')
        filename= '/Global_evol_theta_dist_treshold'+str(rise)+'.png'
        fig_thet_tri.savefig(header+'/Analyse_complete_2'+filename)
        fig_thet_tri.clf()

        print('Variance of each records : ' + str(Var_thet_tbt))
        print('Variance of records where the focus did not change much : ' + str(Var_thet_tbt_tri))

        data = open(header+data_header +"/complete_analysis.txt", 'a')
        data.write('\n Mean Dist   = ' + str(mean_dist_foyer)+ '\n')
        data.write('\n Variance of the thetas :\n')
        data.write('Var_all = ' + str(Var_thet_tbt) + '\n' +
        'Var_tri = ' + str(Var_thet_tbt_tri))
        data.close()


    for mouse in range(4):
        newheader = '/Analyse_complete_2/Mouse_'+str(mouse+3)
        Ellipses_VM = []
        abscisse_thet = [0 for _ in Time_delay]
        Ellipses_VM_S2 = []
        Ellipses_VM_PP = []
        Nb_elps = 0
        Nb_elps_PP = 0
        Nb_elps_S2 = 0
        Dist = [[] for k in range(ray)]
        Dist_S2 = [[] for k in range(ray)]
        Dist_PP = [[] for k in range(ray)]
        first_ellipse = False
        first_ellipse_PP = False
        first_ellipse_S2 = False
        S1_stop = False
        S2_stop = False
        PP_stop = False
        thet = 0
        first_thet = True
        VM_ellipse = [[[0 for k in range(100)] for l in range(100)] for t in range(511)]
        min_vm = min([VM_each_mice[mouse][i][j][k] for i in range(511) for j in range(100) for k in range(100)])
        max_vm = max([VM_each_mice[mouse][i][j][k] for i in range(511) for j in range(100) for k in range(100)])
        for index, time_delay in enumerate(Time_delay) :
            fig, axs = plt.subplots(1, 1)
            axs.imshow(VM_each_mice[mouse][injection_start+time_delay], cmap = RdBu,interpolation='none', vmin =min_vm , vmax =max_vm)
            axs.set_title('Mean VM '+str(injection_start+time_delay))
            filename= '/Mouse_'+str(mouse+3)+'_mean_VM_'+str(injection_start+time_delay)+'.png'
            fig.savefig(header+newheader+filename, transparent=True)
            plt.show(block=True)
            plt.close()
            fig.clf()

            VM_ellipse[injection_start+time_delay] = np.array( (np.array(VM_each_mice[mouse][injection_start+time_delay]) - min_vm ) * 255/(max_vm-min_vm) , dtype = np.uint8)
            blur_vm = cv.GaussianBlur(VM_ellipse[injection_start+time_delay],(9,9),0)
            ret1,VM_ellipse[injection_start+time_delay] = cv.threshold(blur_vm,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
            ret1,VM_ellipse[injection_start+time_delay] = cv.threshold(blur_vm,int(np.floor(ret1*rise)),255,cv.THRESH_BINARY)

            for i in range(100):
                for j in range(100):
                    VM_ellipse_S1[injection_start+time_delay][i][j] = min(VM_ellipse[injection_start+time_delay][i][j],Mask[i][j])
                    VM_ellipse_S2[injection_start+time_delay][i][j] = min(VM_ellipse[injection_start+time_delay][i][j],MaskS2[i][j])
                    VM_ellipse_PP[injection_start+time_delay][i][j] = min(VM_ellipse[injection_start+time_delay][i][j],MaskPP[i][j])


            contours1,hierarchy1 = cv.findContours(np.array(VM_ellipse_S1[injection_start+time_delay]).astype(np.uint8), 1, 2)
            if S2 :
                contoursS2,hierarchyS2 = cv.findContours(np.array(VM_ellipse_S2[injection_start+time_delay]).astype(np.uint8), 1, 2)
            if PP :
                contoursPP,hierarchyPP = cv.findContours(np.array(VM_ellipse_PP[injection_start+time_delay]).astype(np.uint8), 1, 2)

            if len(contours1) != 0:
                cont = max(contours1,key = len)
                if len(cont) >= 5:
                    elps = cv.fitEllipse(cont)
                    area = np.pi*elps[1][0]*elps[1][1]/4
                    ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                    bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                    ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                    by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                    out = 0
                    for t in range(360):
                        ord = elps[0][1]+ay*np.cos(t*2*np.pi/360)+by * np.sin(t*2*np.pi/360)
                        if ord > 70 or ord < 20 :
                            out+=1
                    if out <200  and 20<elps[0][1]<70 : #  and (np.abs(elps[2]-thet) < 30 or first_thet )
                        first_thet = False
                        thet = elps[2]
                        if area < area_limit and not S1_stop:
                            if not first_ellipse :
                                first_ellipse =True
                                centre_S1 = (int(round(elps[0][0])),int(round(elps[0][1])))
                                theta_S1 = elps[2]
                            Nb_elps+=1
                            abscisse_thet[index]+= 1
                            Ellipses_VM.append(elps)
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
                                    # cv.line(VM_ellipse[injection_start+time_delay][:][20:70],(int(round(elps[0][0])),int(round(elps[0][1]))),point, color_ray)
                                    # Dist[count].append(np.linalg.norm([int(round(elps[0][0])),int(round(elps[0][1]))]-np.array(point)))
                                    cv.line(VM_ellipse[injection_start+time_delay],centre_S1,point, (255-255/ray*count,255-255/ray*count,255-255/ray*count))
                                    Dist[count].append(np.linalg.norm(np.array(centre_S1)-np.array(point)))
                        else :
                            Ellipses_VM.append(0)
                            if Nb_elps != 0 :
                                S1_stop = True
                    else :
                        Ellipses_VM.append(0)
                else :
                    Ellipses_VM.append(0)
            else :
                Ellipses_VM.append(0)
            if S2 :
                if len(contoursS2) != 0:
                    cont = max(contoursS2,key = len)
                    if len(cont) >= 5:
                        elps = cv.fitEllipse(cont)
                        area_S2 = np.pi*elps[1][0]*elps[1][1]/4
                        ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                        bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                        ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                        by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                        out = 0
                        for t in range(360):
                            ord = elps[0][1]+ay*np.cos(t*2*np.pi/360)+by * np.sin(t*2*np.pi/360)
                            if ord > 20 :
                                out+=1
                        if  out <= 100 :
                            if area_S2_min < area_S2 < area_S2_max and not S2_stop :
                                if not first_ellipse_S2 :
                                    first_ellipse_S2 = True
                                    centre_S2 = (int(round(elps[0][0])),int(round(elps[0][1])))
                                    theta_S2 = elps[2]*2*np.pi/360
                                Nb_elps_S2+=1
                                Ellipses_VM_S2.append(elps)
                                if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                    cv.ellipse(VM_ellipse[injection_start+time_delay], elps, color_ellipse)


                                    for count,angle in enumerate([k*2*np.pi/ray for k in range(ray)]):
                                        point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_S2)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_S2)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_S2)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_S2)*2*np.pi/360))))
                                        cv.line(VM_ellipse[injection_start+time_delay],centre_S2,point, (255-255/ray*count,255-255/ray*count,255-255/ray*count))
                                        Dist_S2[count].append(np.linalg.norm(np.array(centre_S2)-np.array(point)))
                            else :
                                if Nb_elps_S2 != 0 :
                                    S2_stop = True
            if PP :
                if len(contoursPP) != 0:
                    cont = max(contoursPP,key = len)
                    if len(cont) >= 5:
                        elps = cv.fitEllipse(cont)
                        area_PP = np.pi*elps[1][0]*elps[1][1]/4
                        ax = elps[1][0]/2*np.cos(elps[2]*2*np.pi/360)
                        bx = -elps[1][1]/2*np.sin(elps[2]*2*np.pi/360)
                        ay = elps[1][0]/2*np.sin(elps[2]*2*np.pi/360)
                        by = elps[1][1]/2*np.cos(elps[2]*2*np.pi/360)
                        out = 0
                        for t in range(360):
                            ord = elps[0][1]+ay*np.cos(t*2*np.pi/360)+by * np.sin(t*2*np.pi/360)
                            if ord < 70 :
                                out+=1
                        if out < 250 :
                            if area_PP_min < area_PP < area_PP_max and not PP_stop :
                                if not first_ellipse_PP :
                                    first_ellipse_PP = True
                                    centre_PP = (int(round(elps[0][0])),int(round(elps[0][1])))
                                    theta_PP = elps[2]*2*np.pi/360
                                Nb_elps_PP+=1
                                Ellipses_VM_PP.append(elps)
                                if not any([np.isnan(elps[k][l]) for k in range(2) for l in range(2)]) or np.isnan(elps[2]) :
                                    cv.ellipse(VM_ellipse[injection_start+time_delay], elps, color_ellipse)
                                    for count,angle in enumerate([k*360/ray for k in range(ray)]):
                                        point = (int(round(elps[0][0]+ax*np.cos(angle-(elps[2]-theta_PP)*2*np.pi/360)+bx*np.sin(angle-(elps[2]-theta_PP)*2*np.pi/360))),int(round(elps[0][1]+ay*np.cos(angle-(elps[2]-theta_PP)*2*np.pi/360)+by*np.sin(angle-(elps[2]-theta_PP)*2*np.pi/360))))
                                        cv.line(VM_ellipse[injection_start+time_delay],centre_PP,point, (255-255/ray*count,255-255/ray*count,255-255/ray*count))
                                        Dist_PP[count].append(np.linalg.norm(np.array(centre_PP)-np.array(point)))
                                else :
                                    if Nb_elps_PP != 0 :
                                        PP_stop = True

                fig, axs = plt.subplots(1, 1)

                axs.imshow(VM_ellipse[injection_start+time_delay], cmap = RdBu,interpolation='none')
                if not np.isnan(area):
                    if PP and S2:
                        if not np.isnan(area_PP):
                            if not np.isnan(area_S2):
                                axs.set_title('VM tresh Aire centrale = '+str(round(area)) + ' PP ' + str(round(area_PP))+ ' S2 ' + str(round(area_S2)))
                            else :
                                axs.set_title('VM tresh Aire centrale = '+str(round(area)) + ' PP ' + str(round(area_PP))+ ' S2 NaN')
                        else :
                            if not np.isnan(area_S2):
                                axs.set_title('VM tresh Aire centrale = '+str(round(area)) + ' PP NaN S2 ' + str(round(area_S2)))
                            else :
                                axs.set_title('VM tresh Aire centrale = '+str(round(area)) + ' PP NaN S2 NaN')
                    else :
                        axs.set_title('VM tresh Aire centrale = '+str(round(area)))
                else :
                    if PP and S2 :
                        if not np.isnan(area_PP):
                            if not np.isnan(area_S2):
                                axs.set_title('VM tresh Aire centrale = NaN PP ' + str(round(area_PP))+ ' S2 ' + str(round(area_S2)))
                            else :
                                axs.set_title('VM tresh Aire centrale = NaN PP ' + str(round(area_PP))+ ' S2 NaN')
                        else :
                            if not np.isnan(area_S2):
                                axs.set_title('VM tresh Aire centrale = NaN PP NaN S2 ' + str(round(area_S2)))
                            else :
                                axs.set_title('VM tresh Aire centrale = NaN PP NaN S2 NaN')
                    else :
                        axs.set_title('VM tresh Aire centrale = NaN')
                filename= '/Mouse_'+str(mouse+3)+'_mean_VM_ellipse_'+str(injection_start+time_delay)+'.png'
                fig.savefig(header+newheader+filename, transparent=True)
                plt.show(block=True)
                plt.close()
                fig.clf()
        fig = plt.figure()
        for k in range(ray):
            plt.plot(Time_delay[:Nb_elps], Dist[k], label = 'Angle = ' + str(k*360/ray), color = cm.hsv(k/ray))
        plt.legend(loc = 'upper left')
        filename= '/Distance_mean_VM_ellipse_S1_rays_'+str(ray)+'.png'
        fig.savefig(header+newheader+filename)
        plt.close()
        fig.clf()
        if S2 :
            fig = plt.figure()
            for k in range(ray):
                plt.plot(Time_delay[:Nb_elps_S2], Dist_S2[k], label = 'Angle = ' + str(k*360/ray), color = cm.hsv(k/ray))
            plt.legend(loc = 'upper left')
            filename= '/Distance_mean_VM_ellipse_S2_rays_'+str(ray)+'.png'
            fig.savefig(header+newheader+filename)
            plt.close()
            fig.clf()
        if PP :
            fig = plt.figure()
            for k in range(ray):
                plt.plot(Time_delay[:Nb_elps_PP], Dist_PP[k], label = 'Angle = ' + str(k*360/ray), color = cm.hsv(k/ray))
            plt.legend(loc = 'upper left')
            filename= '/Distance_mean_VM_ellipse_PP_rays_'+str(ray)+'.png'
            fig.savefig(header+newheader+filename)
            plt.close()
            fig.clf()


        fig = plt.figure()
        for index,time_delay in enumerate(Time_delay) :
            if abscisse_thet[index]==1 :
                plt.scatter(time_delay,90-Ellipses_VM[index][2])
        plt.ylim(-90,90)
        plt.title('Evolution theta selon le temps souris '+str(mouse+3))
        filename= '/Evol_theta.png'
        fig.savefig(header+newheader+filename)
        plt.close()
        fig.clf()

        fig = plt.figure()
        for index,time_delay in enumerate(Time_delay) :
            if abscisse_thet[index]==1 :
                plt.scatter(Ellipses_VM[index][0][0],Ellipses_VM[index][0][1], color = cm.hsv(index/len(Time_delay)))# , c = abscisse_thet[index])
        plt.xlim(0,100)
        plt.ylim(0,100)
        plt.title('Foyer des ellipses')
        filename= '/Foyer_ellipses.png'
        fig.savefig(header+newheader+filename)
        plt.close()
        fig.clf()

    N = len(VM_list)

    Mean_vm = np.array([[np.mean([VM_list[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
    fig=plt.figure()
    plt.imshow(Mean_vm)
    plt.colorbar()
    filename= '/Global_mean_VM_treshold='+str(treshold)+'.png'
    fig.savefig(header+'/Analyse_complete_2'+filename)
    plt.close()
    fig.clf()
    np.savetxt(header+data_header+'/Global_mean_VM_treshold='+str(treshold)+'.txt',Mean_vm)
    if MI :
        Mean_mi = np.array([[np.mean([MI_list[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
        fig=plt.figure()
        plt.imshow(Mean_mi)
        plt.colorbar()
        filename= '/Global_mean_MI_treshold='+str(treshold)+'.png'
        fig.savefig(header+'/Analyse_complete_2'+filename)
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
        fig.savefig(header+'/Analyse_complete_2'+filename)
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
            fig.savefig(header+'/Analyse_complete_2'+filename)
            plt.close()
            fig.clf()
            np.savetxt(header+data_header+'/Global_mean_gradient_norm_MI_='+str(treshold)+'.txt',Mean_grad_mi)

            grad_mi_smooth = np.array([[np.mean([Grad_list_MI_smooth[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
            fig=plt.figure()
            plt.imshow(grad_mi_smooth)
            #plt.clim([0.014,0.05])
            plt.colorbar()
            filename= '/Global_mean_gradient_norm_MI_smooth.png'
            fig.savefig(header+'/Analyse_complete_2'+filename)
            plt.close()
            fig.clf()
            np.savetxt(header+'/Analyse_complete_2'+'/Global_mean_gradient_norm_MI_smooth.txt',grad_mi_smooth)


        grad_vm_smooth = np.array([[np.mean([Grad_list_VM_smooth[k][i][j] for k in range(N)]) for j in range(window)] for i in range(window)])
        fig=plt.figure()
        plt.imshow(grad_vm_smooth)
        #plt.clim([0,0.004])
        plt.colorbar()
        filename= '/Global_mean_gradient_norm_VM_smooth.png'
        fig.savefig(header+'/Analyse_complete_2'+filename)
        plt.close()
        fig.clf()
        np.savetxt(header+'/Analyse_complete_2'+'/Global_mean_gradient_VM_smooth.txt',grad_vm_smooth)


    if Ellipse :
        vm = Mean_vm
        if MI :
            mi = Mean_mi
            minimum_mi = min([mi[i][j] for i in range(window) for j in range(window)])
            minimum_vm = min([vm[i][j] for i in range(window) for j in range(window)])
            maximum_mi = max([mi[i][j] for i in range(window) for j in range(window)])
            maximum_vm = max([vm[i][j] for i in range(window) for j in range(window)])
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
        fig.savefig(header+'/Analyse_complete_2'+filename, transparent=True)
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
        fig.savefig(header+'/Analyse_complete_2'+filename)
        plt.close()
        fig.clf()
        if MI :
            fig = plt.figure()
            plt.scatter([Ellipses_mouse_MI[i][1][0] for i in range(len(Ellipses_mouse_MI))],[Ellipses_mouse_MI[i][1][1] for i in range(len(Ellipses_mouse_MI))], c = [i for i in range(len(Ellipses_mouse_MI))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Demi-petit axe')
            plt.title('Repartition MI toutes souris ')
            filename= '/Petit-Grand Axe MI.png'
            fig.savefig(header+'/Analyse_complete_2'+filename)
            plt.close()
            fig.clf()

        # Major axis THETA VM & MI
        fig = plt.figure()
        plt.scatter([Ellipses_mouse_VM[i][1][0] for i in range(len(Ellipses_mouse_VM))],[Ellipses_mouse_VM[i][-1] for i in range(len(Ellipses_mouse_VM))], c = [i for i in range(len(Ellipses_mouse_VM))])
        plt.xlabel('Demi-grand axe')
        plt.ylabel('Theta (degr)')
        plt.title('Repartition VM souris '+str(Mouse_index))
        filename= '/GrandAxe-Theta_VM.png'
        fig.savefig(header+'/Analyse_complete_2'+filename)
        plt.close()
        fig.clf()
        if MI :
            fig = plt.figure()
            plt.scatter([Ellipses_mouse_MI[i][1][0] for i in range(len(Ellipses_mouse_MI))],[Ellipses_mouse_MI[i][-1] for i in range(len(Ellipses_mouse_MI))], c = [i for i in range(len(Ellipses_mouse_MI))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Theta (degr)')
            plt.title('Repartition MI toutes souris')
            filename= '/GrandAxe-Theta_MI.png'
            fig.savefig(header+'/Analyse_complete_2'+filename)
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
                ret1,MI_map = cv.threshold(blur_mi,ret1+10,255,cv.THRESH_BINARY)


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
            fig.savefig(header+'/Analyse_complete_2'+filename, transparent=True)
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
            fig.savefig(header+'/Analyse_complete_2'+filename)
            plt.close()
            fig.clf()
            if MI :
                fig = plt.figure()
                plt.scatter([Ellipses_mouse_grad_MI[i][1][0] for i in range(len(Ellipses_mouse_grad_MI))],[Ellipses_mouse_grad_MI[i][1][1] for i in range(len(Ellipses_mouse_grad_MI))], c = [i for i in range(len(Ellipses_mouse_grad_MI))])
                plt.xlabel('Demi-grand axe')
                plt.ylabel('Demi-petit axe')
                plt.title('Repartition MI totues souris ')
                filename= '/Minor-Major Axis Gradient MI Smooth.png'
                fig.savefig(header+'/Analyse_complete_2'+filename)
                plt.close()
                fig.clf()

            # Major axis THETA GRAD SMOOTH VM & MI
            fig = plt.figure()
            plt.scatter([Ellipses_mouse_grad_VM[i][1][0] for i in range(len(Ellipses_mouse_grad_VM))],[Ellipses_mouse_grad_VM[i][-1] for i in range(len(Ellipses_mouse_grad_VM))], c = [i for i in range(len(Ellipses_mouse_grad_VM))])
            plt.xlabel('Demi-grand axe')
            plt.ylabel('Theta (degr)')
            plt.title('Repartition VM toutes souris ')
            filename= '/MajorAxis-Theta_Grad_VM.png'
            fig.savefig(header+'/Analyse_complete_2'+filename)
            plt.close()
            fig.clf()
            if MI :
                fig = plt.figure()
                plt.scatter([Ellipses_mouse_grad_MI[i][1][0] for i in range(len(Ellipses_mouse_grad_MI))],[Ellipses_mouse_grad_MI[i][-1] for i in range(len(Ellipses_mouse_grad_MI))], c = [i for i in range(len(Ellipses_mouse_grad_MI))])
                plt.xlabel('Demi-grand axe')
                plt.ylabel('Theta (degr)')
                plt.title('Repartition MI toutes souris ')
                filename= '/MajorAxis-Theta_Grad_MI.png'
                fig.savefig(header+'/Analyse_complete_2'+filename)
                plt.close()
