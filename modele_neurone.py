import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
folder = 'modele_tau'
##########################
#Integrate and Fire model#
##########################

# time step for integration
DT=0.01
DTlist=np.linspace(0,0.10,10)
# number of iterations
Max_It=1000

#Parameters
Cm=1.
I=1.
v=0
#initial condition
V=0.
#Threshold for spike detection
Threshold= 0.8
#Create a list to save the value of the variable at each time step
Vm=[[] for i in range(10)]

# iterate
for j in range(10):
    DT = DTlist[j]
    for i in range(Max_It):

        # calculate the variation of V at each time step
        dVdt = (-V+I)/Cm
        # calculate the new value of V at each time step
        V = V+dVdt*DT
        if V>Threshold:
            V=0
        # Save the value of V into a list
        Vm[j].append(V)

    # create the time vector
    TimeVec=[i*DT for i in range(Max_It)]
    # Plot the Vm over time
    plt.plot(TimeVec, Vm[j])
    plt.show()


# time step for integration
DT=0.1
# number of iterations
Max_It=10000

#Parameters
Cm=0.15
I=4.
gL=0.01
El=-65.
deltaT=2.
Vt=-50.
a=1.
b=0.02
tau = 500
Treshold=-50
refractorytime=5

Refrac_It = refractorytime/DT

#initial condition
V=-65.
w= 0
#Create a list to save the value of the variable at each time step
Vm=[]
W=[]
def f(V):
    return -gL*(V-El)+gL*deltaT*np.exp((V-Vt)/deltaT)

# iterate

Spike = Refrac_It

for tau in np.arange(400,600,25):
    Vm = []
    for i in range(Max_It):
        if Spike>=Refrac_It:
            if 1000<i<1200:
                I = 1.
            else:
                I =0.
            # calculate the variation of V at each time step
            dVdt = (f(V)+I-w)/Cm
            dwdt = (a*(V-El)-w)/tau
            # calculate the new value of V at each time step
            V = V+dVdt*DT
            w = w+dwdt*DT
            if V>Treshold:
                V=-65.
                w+=b
                Spike=1
            # Save the value of V into a list
        else :
            Spike+=1
            V=-65
        Vm.append(V)
        W.append(w)
    # create the time vector
    TimeVec=[i*DT for i in range(Max_It)]

    # Plot the Vm over time
    fig =plt.figure()
    plt.plot(TimeVec, Vm)
    plt.title('Modele pour tau = '+str(tau))
    fig.savefig(folder+'/Influence_tau ='+str(tau)+'.png')
    plt.close()
    fig.clf()
