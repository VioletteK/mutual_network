import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
Cm = 0.15
gl = 0.01
El = -65
deltat = 2
Vt = -50
a = 4
b = 0.02
I = 0.5

t = 0
T = 500
y = np.array([-65.0,0])
h = 0.1

folder = 'modele_tau'

def pointMilieuMultiD(fun, t0, T, y0, h):

    """ fun est le second membre de l'equation y' = fun(t,y)
        [t0,t0+T] est l'intervalle de temps sur lequel on resout
        y0 est un numpy array : c'est la condition initiale y(t0)=y0
        h est le pas de temps"""
    """euler renvoie
        * tps un numpy array contenant la discretisation de l'intervalle de temps
        * sol un numpy ndarray de taille (tps.size, y0.size) """

    tps = np.arange(t0, t0+T+h, h)
    d = y0.size #
    N = tps.size
    sol = np.zeros((N, d)) #
    sol[0]=y0
    for i in range(N-1):
        p1 = fun(tps[i],sol[i])
        sol[i+1] = sol[i] + h*fun(tps[i]+0.5*h, sol[i]+0.5*h*p1)

    return [tps, sol]


for tau in np.arange(400,575,25):
    def Adex(t,y):
        V = y[0]
        w = y[1]
        if V >= -50 :
            return np.array([-65,w+b])
        else :
            return np.array([
            (-gl/Cm)*((V-El)-deltat*np.exp((V-Vt)/deltat))-w/Cm+I,
            1/tau*(a*(V-El)-w)
            ])
    X,Y = pointMilieuMultiD(Adex, t, T, y, h )
    fig =plt.figure()
    plt.plot(X,Y[:,1])
    plt.title('Modele pour tau = '+str(tau))
    fig.savefig(folder+'/Influence_tau ='+str(tau)+'.png')
    plt.close()
    fig.clf()
