# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

siteNum = 3*3

# theoretical value
theory = pd.read_csv("theoretical.csv",header=None)
theory = theory.dropna(axis=1, how='all')
theoryU = theory.to_numpy()[0] / siteNum
theoryC = theory.to_numpy()[1] / siteNum
theoryM = theory.to_numpy()[2] / siteNum
theoryChi = theory.to_numpy()[3] / siteNum

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

resultU = np.mean(resultE/siteNum,axis=0)
deviationU = np.sqrt(np.var(resultE/siteNum,axis=0)/(size[0]-1))
resultC = np.mean(specificHeat/siteNum,axis=0)
deviationC = np.sqrt(np.var(specificHeat/siteNum,axis=0)/(size[0]-1))
resultaveM = np.mean(resultM/siteNum,axis=0)
deviationM = np.sqrt(np.var(resultM/siteNum,axis=0)/(size[0]-1))
resultChi = np.mean(susceptibility/siteNum,axis=0)
deviationChi = np.sqrt(np.var(susceptibility/siteNum,axis=0)/(size[0]-1))

# %% Draw Physical Quantities
baseline = np.zeros(size[1],dtype='float')

fig = plt.figure(figsize=(16,8),dpi=150)
ax1 = fig.add_subplot(221)
ax1.errorbar(temperature, resultU, yerr=deviationU, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
ax1.plot(temperature,theoryU,linewidth=3, alpha=0.2)
ax2 = ax1.twinx()
ax2.errorbar(temperature, resultU-theoryU, yerr=deviationU, marker='o',
          ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
ax2.plot(temperature, baseline, linestyle='dashed')
ax1.set_title("E per site ~ T")


specificHeatax1 = fig.add_subplot(222)
specificHeatax1.errorbar(temperature, resultC, yerr=deviationC, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
specificHeatax1.plot(temperature,theoryC,linewidth=3, alpha=0.2)
specificHeatax2 = specificHeatax1.twinx()
specificHeatax2.errorbar(temperature, resultC-theoryC, yerr=deviationC, marker='o',
          ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
specificHeatax2.set_ylim(-0.01,0.01)
specificHeatax2.plot(temperature, baseline, linestyle='dashed')
specificHeatax1.set_title("C per site ~ T")

magax1 = fig.add_subplot(223)
magax1.errorbar(temperature, resultaveM, yerr=deviationM, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
magax1.plot(temperature,theoryM,linewidth=3, alpha=0.2)
magax2 = magax1.twinx()
magax2.errorbar(temperature, resultaveM-theoryM, yerr=deviationM, marker='o',
          ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
magax2.set_ylim(-0.02,0.02)
magax2.plot(temperature, baseline, linestyle='dashed')
magax1.set_title("M per site ~ T")

susceptibilityax1 = fig.add_subplot(224)
susceptibilityax1.errorbar(temperature, resultChi, yerr=deviationChi, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
susceptibilityax1.plot(temperature,theoryChi,linewidth=3, alpha=0.2)
susceptibilityax2 = susceptibilityax1.twinx()
susceptibilityax2.errorbar(temperature, resultChi-theoryChi, yerr=deviationChi, marker='o',
          ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
susceptibilityax2.set_ylim(-0.02,0.02)
susceptibilityax2.plot(temperature, baseline, linestyle='dashed')
susceptibilityax1.set_title(r"$\chi$ per site ~ T")

plt.savefig("IsingVerify.png")