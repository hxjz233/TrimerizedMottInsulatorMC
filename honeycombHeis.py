# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# theoretical value for 2*1 case
# theoryU = [-1.95, -1.9, -1.85, -1.80009, -1.75067, -1.70255, -1.65662, -1.61357, \
# -1.57377, -1.53731, -1.50412, -1.47399, -1.44666, -1.42186, -1.39934, \
# -1.37885, -1.36016, -1.34308, -1.32742, -1.31304, -1.29978, -1.28755, \
# -1.27622, -1.26571, -1.25594, -1.24683, -1.23833, -1.23037, -1.22291, \
# -1.2159, -1.20931, -1.2031, -1.19724, -1.1917, -1.18645, -1.18148, \
# -1.17677, -1.17228, -1.16802, -1.16395, -1.16008, -1.15638, -1.15285, \
# -1.14947, -1.14623, -1.14313, -1.14016, -1.13731, -1.13457]
# theoryC = [0.5, 0.5, 0.499856, 0.49773, 0.489258, 0.471647, 0.44579, 0.414629, \
# 0.381241, 0.347956, 0.316243, 0.286876, 0.260168, 0.236144, 0.214674, \
# 0.195555, 0.178552, 0.16343, 0.14997, 0.137969, 0.127248, 0.117649, \
# 0.109033, 0.10128, 0.0942856, 0.0879603, 0.0822253, 0.0770127, \
# 0.0722637, 0.0679269, 0.0639576, 0.0603169, 0.0569704, 0.0538882, \
# 0.0510438, 0.0484141, 0.0459784, 0.0437185, 0.0416183, 0.0396632, \
# 0.0378405, 0.0361387, 0.0345475, 0.0330577, 0.0316609, 0.0303497, \
# 0.0291172, 0.0279575, 0.026865]
theoryU = np.zeros(49)
theoryC = np.zeros(49)
    
siteNum = 1

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

# resultM = pd.read_csv("resultM.csv",header=None)
# resultM = resultM.dropna(axis=1, how='all')
# resultM = resultM.to_numpy()

# resultMsq = pd.read_csv("resultMsq.csv",header=None)
# resultMsq = resultMsq.dropna(axis=1, how='all')
# resultMsq = resultMsq.to_numpy()


specificHeat = np.zeros((size[0],size[1]),dtype='float')
for i in range(size[0]):
    for j in range(size[1]):
        specificHeat[i][j] = (resultEsq[i][j] - resultE[i][j] ** 2) / temperature[j] ** 2
        
# susceptibility = np.zeros((size[0],size[1]),dtype='float')
# for i in range(size[0]):
#     for j in range(size[1]):
#         susceptibility[i][j] = (resultMsq[i][j] - resultM[i][j] ** 2) / temperature[j]

resultU = np.mean(resultE/siteNum,axis=0)
deviationU = np.sqrt(np.var(resultE/siteNum,axis=0)/(size[0]-1))
resultC = np.mean(specificHeat/siteNum,axis=0)
deviationC = np.sqrt(np.var(specificHeat/siteNum,axis=0)/(size[0]-1))
# resultaveM = np.mean(resultM/siteNum,axis=0)
# deviationM = np.sqrt(np.var(resultM/siteNum,axis=0)/(size[0]-1))
# resultChi = np.mean(susceptibility/siteNum,axis=0)
# deviationChi = np.sqrt(np.var(susceptibility/siteNum,axis=0)/(size[0]-1))

# %% Draw Physical Quantities
baseline = np.zeros(size[1],dtype='float')

fig = plt.figure(figsize=(8,8),dpi=150)
ax1 = fig.add_subplot(211)
ax1.errorbar(temperature, resultU, yerr=deviationU, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
ax1.plot(temperature,theoryU,linewidth=3, alpha=0.2)
ax2 = ax1.twinx()
ax2.errorbar(temperature, resultU-theoryU, yerr=deviationU, marker='o',
          ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
ax2.plot(temperature, baseline, linestyle='dashed')
ax1.set_title("E per site ~ T")


specificHeatax1 = fig.add_subplot(212)
specificHeatax1.errorbar(temperature, resultC, yerr=deviationC, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
specificHeatax1.plot(temperature,theoryC,linewidth=3, alpha=0.2)
specificHeatax2 = specificHeatax1.twinx()
specificHeatax2.errorbar(temperature, resultC-theoryC, yerr=deviationC, marker='o',
          ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
specificHeatax2.set_ylim(-0.01,0.01)
specificHeatax2.plot(temperature, baseline, linestyle='dashed')
specificHeatax1.set_title("C per site ~ T")

# magax1 = fig.add_subplot(223)
# magax1.errorbar(temperature, resultaveM, yerr=deviationM, marker='o',
#         ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
# magax1.plot(temperature,theoryM,linewidth=3, alpha=0.2)
# magax2 = magax1.twinx()
# magax2.errorbar(temperature, resultaveM-theoryM, yerr=deviationM, marker='o',
#           ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
# magax2.plot(temperature, baseline, linestyle='dashed')
# magax1.set_title("M per site ~ T")

# susceptibilityax1 = fig.add_subplot(224)
# susceptibilityax1.errorbar(temperature, resultChi, yerr=deviationChi, marker='o',
#         ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
# susceptibilityax1.plot(temperature,theoryChi,linewidth=3, alpha=0.2)
# susceptibilityax2 = susceptibilityax1.twinx()
# susceptibilityax2.errorbar(temperature, resultChi-theoryChi, yerr=deviationChi, marker='o',
#           ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
# susceptibilityax2.plot(temperature, baseline, linestyle='dashed')
# susceptibilityax1.set_title(r"$\chi$ per site ~ T")

plt.savefig("Heisenberg2by1")