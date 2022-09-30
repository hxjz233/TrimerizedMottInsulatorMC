# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

temperature = pd.read_csv("temperature.csv",header=None)
temperature = temperature.dropna(axis=1, how='all')
temperature = temperature.to_numpy()[0]

resultE = pd.read_csv("resultE.csv",header=None)
resultE = resultE.dropna(axis=1, how='all')
size = resultE.shape
resultE = resultE.to_numpy()

resultEsq = pd.read_csv("resultEsq.csv",header=None)
resultEsq = resultEsq.dropna(axis=1, how='all')
resultEsq = resultEsq.to_numpy()

resultM = pd.read_csv("resultM.csv",header=None)
resultM = resultM.dropna(axis=1, how='all')
resultM = resultM.to_numpy()

resultMsq = pd.read_csv("resultMsq.csv",header=None)
resultMsq = resultMsq.dropna(axis=1, how='all')
resultMsq = resultMsq.to_numpy()

specificHeat = np.zeros((size[0],size[1]),dtype='float')
for i in range(size[0]):
    for j in range(size[1]):
        specificHeat[i][j] = (resultEsq[i][j] - resultE[i][j] ** 2) / temperature[j] ** 2
        
susceptibility = np.zeros((size[0],size[1]),dtype='float')
for i in range(size[0]):
    for j in range(size[1]):
        susceptibility[i][j] = (resultMsq[i][j] - resultM[i][j] ** 2) / temperature[j]

resultC = np.mean(specificHeat,axis=0) / 64
deviationC = np.std(specificHeat,axis=0) / 64
resultChi = np.mean(susceptibility,axis=0)
deviationChi = np.std(susceptibility,axis=0)
# print(resultC)
# print(deviationC)
# print(resultChi)
# print(deviationChi)
# %% 
plt.figure(figsize=(8,8),dpi=200)
plt.subplot(2,1,1)
plt.errorbar(temperature, resultC, yerr=deviationC, marker='s',
         ms=2, mew=4, capsize=2, capthick=1)
plt.title("C per site ~ T")

plt.subplot(2,1,2)
plt.errorbar(temperature, resultChi, yerr=deviationChi, marker='s',
             ms=2, mew=4, capsize=2, capthick=1)
plt.title(r"$\chi$ per site ~ T")
plt.savefig("isingErrorbar.png")
