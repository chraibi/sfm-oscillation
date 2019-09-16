"""
plot trajectories
input *txt-file containing the trajectories
"""

import matplotlib.pyplot as plt
from sys import *
from numpy import *
ms = 10
fs = 20
lw = 2
file1 = argv[1]
print('load file ',file1)
A1 = loadtxt(file1)#, delimiter=",")
print('finished loading')
ids = unique(A1[:,0]).astype(int)
frames = A1[:,1]
print("ids=", ids)
Length = 200
figname = file1.split(".txt")[0]+".pdf"
for i in ids[::]:
    #print("i=",i)
    p = A1[ A1[:,0] == i ]
    x1 = p[:,1] #time
    y1 = fmod(p[:,2], Length) #x
    abs_d_data = abs(diff(y1))
    abs_mean = abs_d_data.mean()
    abs_std = abs_d_data.std()
    if abs_std <=0.5*abs_mean:
        T = []
        plt.plot(x1, y1,'-r',lw=lw)
    else:
        T = nonzero(abs_d_data > abs_mean + 3*abs_std)[0] 
        start = 0
        for t in T:
            print("start=",start, "t=",t)
            #plt.plot(y1[start:t], x1[start:t],'-k',ms=lw, lw=lw)
            plt.plot(x1[start:t],y1[start:t],'-k',ms=lw, lw=lw)
            start = t+1
        plt.plot(x1[start:],y1[start:],'-k',ms=lw, lw=lw)
        #plt.plot(y1[start:], x1[start:],'-k',ms=lw, lw=lw)
plt.ylabel(r'$x_n\; \rm{[m]}$', size=20)
plt.xlabel(r'$t\; \rm{[s]}$', size=20)

plt.savefig(figname)
print('fgname: ', figname)
plt.show()
