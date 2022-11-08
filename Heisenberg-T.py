# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

siteNum = 3*2

paraNum = 8

field = ['E','Esq','M','Msq','Px','Pxsq','Py','Pysq']

title = ["E per site ~ T", "C per site ~ T", 
         "M per site ~ T", r"$\chi_m$ per site ~ T",
         r"$P_x$ per site ~ T",r"$\chi_{e,x}$ per site ~ T", 
         r"$P_y$ per site ~ T",r"$\chi_{e,y}$ per site ~ T"]

superTitle = r"Heisenberg, $s_\theta=\pi/4 \quad s_\phi=\pi/5 \quad \tau_\theta = \pi/6 \quad \tau_\phi = \pi/7$"

# theoretical value
# theory = pd.read_csv("theoretical.csv",header=None)
# theory = theory.dropna(axis=1, how='all')
# theorySize = theory.shape
# theory = theory.to_numpy()
# theory /= siteNum

temperature = pd.read_csv("temperature.csv",header=None)
temperature = temperature.dropna(axis=1, how='all')
temperature = temperature.to_numpy()[0]

size = pd.read_csv("resultE.csv",header=None).dropna(axis=1, how='all').shape

result = np.zeros((paraNum,size[0],size[1]),dtype='float')

for i in range(paraNum):
    result[i] = pd.read_csv("result"+field[i]+".csv",header=None).dropna(axis=1, how='all').to_numpy()

for i in range(size[0]):
    for j in range(size[1]):
        result[1][i][j] = (result[1][i][j] - result[0][i][j] ** 2) / temperature[j] ** 2

for para in range(3,8,2):
    for i in range(size[0]):
        for j in range(size[1]):
            result[para][i][j] = (result[para][i][j] - result[para-1][i][j] ** 2) / temperature[j]

resultMean = np.zeros((paraNum,size[1]),dtype='float')
resultDev = np.zeros((paraNum,size[1]),dtype='float')

for para in range(paraNum):
    resultMean[para] = np.mean(result[para]/siteNum,axis=0)
    resultDev[para] = np.sqrt(np.var(result[para]/siteNum,axis=0)/(size[0]-1))

baseline = np.zeros(size[1],dtype='float')

fig = plt.figure(figsize=(16,16),dpi=150)

for i in range(4):
    for j in range(2):
        plotAxes1 = fig.add_subplot(4,2,i*2+j+1)
        plotAxes1.errorbar(temperature, resultMean[i*2+j], yerr=resultDev[i*2+j], marker='o',
                ms=2, mew=4, capsize=2, capthick=1)
        # plotAxes1.plot(temperature,theory[i*2+j],linewidth=3, alpha=0.2)
        # plotAxes2 = plotAxes1.twinx()
        # plotAxes2.errorbar(temperature, resultMean[i*2+j]-theory[i*2+j], yerr=resultDev[i*2+j],
        #                    marker='o',ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
        # simuDevFromTheo = np.abs(resultMean[i*2+j]-theory[i*2+j]).mean() * 2
        # plotAxes2.set_ylim(-simuDevFromTheo,simuDevFromTheo)
        plotAxes1.plot(temperature, baseline, linestyle='dashed')
        plotAxes1.set_title(title[i*2+j])

plt.suptitle(superTitle,fontsize = 'xx-large',weight = 'extra bold', y=0.92)

plt.savefig("Heisenberg-T_2D.png")