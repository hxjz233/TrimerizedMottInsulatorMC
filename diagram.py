# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# theoretical value for 2*2 case
theoryC = [0.0000432135, 0.000431885, 0.00213136, 0.00680649, 0.0163197,
0.0320823, 0.0546476, 0.083632, 0.117872, 0.15569, 0.195179,
0.234457, 0.271847, 0.305999, 0.335946, 0.361096, 0.381203, 0.396304,
0.406649, 0.412638, 0.414759, 0.413539, 0.40951, 0.403179, 0.395012,
0.385426, 0.374787, 0.363406, 0.351544, 0.339417, 0.3272, 0.315033,
0.303024, 0.291259, 0.2798, 0.268693, 0.257967, 0.247644, 0.237733,
0.228239, 0.219159, 0.210486, 0.202212, 0.194325, 0.18681]
theoryChi = [8., 6.66661, 5.71397, 4.99887, 4.44138, 3.9933, 3.62379, 3.31228, \
3.04462, 2.81091, 2.60405, 2.41894, 2.25186, 2.10006, 1.96144, \
1.83442, 1.71775, 1.61039, 1.51149, 1.42031, 1.3362, 1.25858, \
1.18692, 1.12072, 1.05955, 1.00299, 0.950652, 0.90219, 0.85728, \
0.815624, 0.77695, 0.741008, 0.707571, 0.676429, 0.647394, 0.620291, \
0.594963, 0.571266, 0.54907, 0.528256, 0.508714, 0.490347, 0.473064, \
0.456783, 0.44143]
    
siteNum = 4

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

resultE = np.mean(resultE/siteNum,axis=0)
deviationE = np.sqrt(np.var(resultE/siteNum,axis=0)/(size[0]-1))

resultC = np.mean(specificHeat/siteNum,axis=0)
deviationC = np.sqrt(np.var(specificHeat/siteNum,axis=0)/(size[0]-1))
# resultChi = np.mean(susceptibility/siteNum,axis=0)
# deviationChi = np.sqrt(np.var(susceptibility/siteNum,axis=0)/(size[0]-1))

# print(resultC)
# print(deviationC)
# print(resultChi)
# print(deviationChi)
# # %% 
# baseline = np.zeros(size[1],dtype='float')

# plt.figure(figsize=(16,8),dpi=200)
# plt.subplot(2,2,1)
# plt.errorbar(temperature, resultC, yerr=deviationC, marker='s',
#          ms=2, mew=4, capsize=2, capthick=1)
# plt.plot(temperature,theoryC,linewidth=3)
# plt.title("C per site ~ T")
plt.errorbar(temperature, resultE, yerr=deviationE, marker='s',
        ms=2, mew=4, capsize=2, capthick=1)
#plt.plot(temperature,theoryC,linewidth=3)
plt.title("E per site ~ T")
plt.savefig("Heisenberg2by2E-T")

# plt.subplot(2,2,2)
# plt.errorbar(temperature, resultChi, yerr=deviationChi, marker='s',
#               ms=2, mew=4, capsize=2, capthick=1)
# plt.plot(temperature,theoryChi,linewidth=3)
# plt.title(r"$\chi$ per site ~ T")

# plt.subplot(2,2,3)
# plt.errorbar(temperature, resultC-theoryC, yerr=deviationC, marker='s',
#          ms=2, mew=4, capsize=2, capthick=1)
# plt.plot(temperature,baseline)
# plt.title("C per site, MC-ED ~ T")

# plt.subplot(2,2,4)
# plt.errorbar(temperature, resultChi-theoryChi, yerr=deviationChi, marker='s',
#          ms=2, mew=4, capsize=2, capthick=1)
# plt.plot(temperature,baseline)
# plt.ylim(-0.01,0.01)
# plt.title(r"$\chi$ per site, MC-ED ~ T")

# plt.savefig("isingErrorbarN=2.png")
# %% 
baseline = np.zeros(size[1],dtype='float')

# plt.figure(figsize=(8,8),dpi=200)
# plt.subplot(2,1,1)
plt.errorbar(temperature, resultC, yerr=deviationC, marker='s',
          ms=2, mew=4, capsize=2, capthick=1)
# plt.plot(temperature,theoryC,linewidth=3)
plt.title("C per site ~ T")
plt.savefig("Heisenberg2by2C-T")

# plt.subplot(2,1,2)
# plt.errorbar(temperature, resultC-theoryC, yerr=deviationC, marker='s',
#          ms=2, mew=4, capsize=2, capthick=1)
# plt.plot(temperature,baseline)
# plt.title("C per site, MC-ED ~ T")
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
plt.scatter(np.arange(2000),avgE[0:2000])
plt.title("E' per site ~ MC Step")
plt.savefig("HeisenbergThermoEprime")