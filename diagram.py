# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# theoretical value for 2*2 case
theoryU = [-2.9, -2.80009, -2.70255, -2.61357, -2.53731, -2.47399, -2.42186, \
-2.37885, -2.34308, -2.31304, -2.28755, -2.26571, -2.24683, -2.23037, \
-2.2159, -2.2031, -2.1917, -2.18148, -2.17228, -2.16395, -2.15638, \
-2.14947, -2.14313, -2.13731, -2.13193, -2.12696, -2.12234, -2.11805, \
-2.11404, -2.1103, -2.10679, -2.10349, -2.1004, -2.09748, -2.09472, \
-2.09212, -2.08965, -2.08732, -2.0851, -2.08299, -2.08098, -2.07907, \
-2.07724, -2.0755, -2.07383, -2.07224, -2.07071, -2.06924, -2.06784]
theoryC = [0.999999, 0.99546, 0.943294, 0.829258, 0.695913, 0.573753, 0.472287, \
0.39111, 0.326861, 0.275938, 0.235297, 0.202559, 0.175921, 0.154025, \
0.135854, 0.120634, 0.107776, 0.0968282, 0.087437, 0.0793264, \
0.0722774, 0.0661153, 0.0606993, 0.055915, 0.0516691, 0.0478844, \
0.0444971, 0.041454, 0.0387103, 0.0362283, 0.033976, 0.031926, \
0.0300551, 0.028343, 0.0267723, 0.0253281, 0.023997, 0.0227678, \
0.0216302, 0.0205755, 0.0195958, 0.0186841, 0.0178344, 0.0170412, \
0.0162996, 0.0156052, 0.0149541, 0.0143429, 0.0137682]
theoryM = [0.9, 0.800091, 0.702549, 0.613567, 0.537315, 0.473987, 0.421864, \
0.378851, 0.343078, 0.313035, 0.287548, 0.265713, 0.246834, 0.230371, \
0.215905, 0.203102, 0.191699, 0.181483, 0.172282, 0.163953, 0.156381, \
0.149468, 0.143133, 0.137308, 0.131932, 0.126958, 0.122342, 0.118047, \
0.114042, 0.110297, 0.106788, 0.103495, 0.100397, 0.0974784, \
0.0947238, 0.0921198, 0.0896544, 0.087317, 0.0850978, 0.0829882, \
0.0809802, 0.0790667, 0.0772413, 0.075498, 0.0738314, 0.0722365, \
0.0707089, 0.0692443, 0.0678391]
theoryChi = [0.0999999, 0.199092, 0.282988, 0.331703, 0.347956, 0.344252, \
0.330601, 0.312888, 0.294175, 0.275938, 0.258827, 0.243071, 0.228697, \
0.215636, 0.203781, 0.193014, 0.18322, 0.174291, 0.16613, 0.158653, \
0.151783, 0.145454, 0.139608, 0.134196, 0.129173, 0.124499, 0.120142, \
0.116071, 0.11226, 0.108685, 0.105326, 0.102163, 0.0991817, 0.096366, \
0.0937031, 0.091181, 0.0887891, 0.0865176, 0.0843578, 0.0823019, \
0.0803426, 0.0784733, 0.076688, 0.0749813, 0.0733482, 0.071784, \
0.0702845, 0.0688458, 0.0674643]
    
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

fig = plt.figure(figsize=(8,16),dpi=150)
ax1 = fig.add_subplot(411)
ax1.errorbar(temperature, resultU, yerr=deviationU, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
ax1.plot(temperature,theoryU,linewidth=3, alpha=0.2)
ax2 = ax1.twinx()
ax2.errorbar(temperature, resultU-theoryU, yerr=deviationU, marker='o',
         ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
ax2.plot(temperature, baseline, linestyle='dashed')
ax1.set_title("E per site ~ T")


specificHeatax1 = fig.add_subplot(412)
specificHeatax1.errorbar(temperature, resultC, yerr=deviationC, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
specificHeatax1.plot(temperature,theoryC,linewidth=3, alpha=0.2)
specificHeatax2 = specificHeatax1.twinx()
specificHeatax2.errorbar(temperature, resultC-theoryC, yerr=deviationC, marker='o',
          ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
specificHeatax2.plot(temperature, baseline, linestyle='dashed')
specificHeatax1.set_title("C per site ~ T")

magax1 = fig.add_subplot(413)
magax1.errorbar(temperature, resultaveM, yerr=deviationM, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
magax1.plot(temperature,theoryM,linewidth=3, alpha=0.2)
magax2 = magax1.twinx()
magax2.errorbar(temperature, resultaveM-theoryM, yerr=deviationM, marker='o',
          ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
magax2.plot(temperature, baseline, linestyle='dashed')
magax1.set_title("M per site ~ T")

susceptibilityax1 = fig.add_subplot(414)
susceptibilityax1.errorbar(temperature, resultChi, yerr=deviationChi, marker='o',
        ms=2, mew=4, capsize=2, capthick=1, alpha=0.2)
susceptibilityax1.plot(temperature,theoryChi,linewidth=3, alpha=0.2)
susceptibilityax2 = susceptibilityax1.twinx()
susceptibilityax2.errorbar(temperature, resultChi-theoryChi, yerr=deviationChi, marker='o',
          ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
susceptibilityax2.plot(temperature, baseline, linestyle='dashed')
susceptibilityax1.set_title(r"$\chi$ per site ~ T")

plt.savefig("Heisenberg1by1ExternalField")

# %% thermalization

valid = pd.read_csv("validateEquilibrium.csv",header=None)
valid = valid.dropna(axis=1, how='all')
sizeValid = valid.shape
valid = valid.to_numpy()
meanforValid = np.mean(valid/siteNum,axis=0)
devforValid = np.sqrt(np.var(valid/siteNum,axis=0)/(size[0]-1))

avgE = []

for i in range(0,50):
    avgE.append(np.mean(meanforValid[0:i+50]))

for i in range(50,sizeValid[1]-50):
    avgE.append(np.mean(meanforValid[i-50:i+50]))

# valid = np.arange(sizeValid[1])
# plt.errorbar(np.arange(100), meanforValid[0:100], yerr=devforValid[0:100], marker='s',
#         ms=2, mew=4, capsize=2, capthick=1)
# plt.title("E ~ MC Step")
# plt.savefig("HeisenbergThermoE")
plt.plot(np.arange(9000),avgE[0:9000])
plt.title("E' per site ~ MC Step")
plt.savefig("HeisenbergThermoEprime")