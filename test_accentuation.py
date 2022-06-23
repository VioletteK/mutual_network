import numpy as np
import matplotlib
matplotlib.use('Agg') # to be used when DISPLAY is undefined

import matplotlib.pyplot as plt

def f(x):
    return 0.008*(np.tanh(300*x-2)+1)

fig = plt.figure()
X = np.linspace(0,0.016,100)
plt.plot(X,f(X))
plt.show()
fig.savefig('Salut')
