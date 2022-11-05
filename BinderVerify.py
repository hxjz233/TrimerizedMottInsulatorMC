import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from matplotlib import cm

BinderNum = 3
simuCount = 10
temperatureCount = 49
qCount = 3 * 2

superTitle = r"Deduced Ising, $s_\theta=\pi/4 \quad s_\phi=\pi/5 \quad \tau_\theta = \pi/6 \quad \tau_\phi = \pi/7$, "+"3*2 triangular lattice"

temperature = pd.read_csv("temperature.csv",header=None)
temperature = temperature.dropna(axis=1, how='all')
temperature = temperature.to_numpy()[0]

simufilename='./resultBinder.csv'
simudata = np.zeros((BinderNum,simuCount,temperatureCount,qCount),dtype='float')
predata = []
with open(simufilename) as csvfile:
    csv_reader = csv.reader(csvfile)
    for row in csv_reader:
        if row == []:
            continue
        else:
            predata.append(row)
rowNum = 0
for binderIndex in range(BinderNum):
    for simuIndex in range(simuCount):
        for tempIndex in range(temperatureCount):
            for qIndex in range(qCount):
                if predata[rowNum][qIndex] != '':
                    simudata[binderIndex][simuIndex][tempIndex][qIndex] = predata[rowNum][qIndex]
            rowNum += 1

theofilename='./theoreticalBinder.csv'
theodata = np.zeros((BinderNum,temperatureCount,qCount),dtype='float')
predata = []
with open(theofilename) as csvfile:
    csv_reader = csv.reader(csvfile)
    for row in csv_reader:
        if row == []:
            continue
        else:
            predata.append(row)
rowNum = 0
for binderIndex in range(BinderNum):
    for tempIndex in range(temperatureCount):
        for qIndex in range(qCount):
            if predata[rowNum][qIndex] != '':
                theodata[binderIndex][tempIndex][qIndex] = predata[rowNum][qIndex]
        rowNum += 1

# %%
relaMean = np.zeros((BinderNum,temperatureCount,qCount),dtype='float')
relaDev = np.zeros((BinderNum,temperatureCount,qCount),dtype='float')
for binderIndex in range(BinderNum):
        for tempIndex in range(temperatureCount):
            for qIndex in range(qCount):
                relaMean[binderIndex,tempIndex,qIndex] = np.mean(simudata[binderIndex,:,tempIndex,qIndex])
                relaDev[binderIndex,tempIndex,qIndex] = np.sqrt(np.var(simudata[binderIndex,:,tempIndex,qIndex])/(simuCount-1))
                
relaMean -= theodata
relaMean /= theodata
relaDev /= theodata

# # %%
# X = np.arange(qCount)
# Y = temperature
# X,Y = np.meshgrid(X,Y)
# base = np.zeros((temperatureCount,qCount),dtype='float')
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(8,8),dpi=150)

# surfMean = ax.scatter(X, Y, relaMean[0], color='m',
#                        linewidth=0, antialiased=False, alpha=0.25)
# surfDevUp = ax.plot_surface(X, Y, relaMean[0]+relaDev[0], color='r',
#                        linewidth=0, antialiased=False, alpha=0.25)
# surfDevDown = ax.plot_surface(X, Y, relaMean[0]-relaDev[0], color='b',
#                        linewidth=0, antialiased=False, alpha=0.25)
# surfBase = ax.plot_surface(X, Y, base,
#                        linewidth=0, antialiased=False, alpha=0.25)
# ax.set(zlim=(-0.02,0.02))
# ax.view_init(15,30,45)
# plt.show()
# %%
X = np.arange(qCount)
base = np.zeros((qCount),dtype='float')
fig = plt.figure(figsize=(18,12),dpi=150)

xlabel = [r"$T=0.1$",r"$T=1.3$",r"$T=2.5$",r"$T=3.7$"]
ylabel = [r"$U_M^z(q)$",r"$U_P^x(q)$",r"$U_P^y(q)$"]
qflattened = [r"q=(0,0)",r"q=(0,1)",r"q=(1,0)",r"q=(1,1)",r"q=(2,0)",r"q=(2,1)"]

for i in range(3):
    for j in range(4):
        plotAxes = fig.add_subplot(3,4,i*4+j+1)
        plotAxes.errorbar(X,relaMean[i,j*12],yerr=relaDev[i,j*12],
                          marker='o',ms=1, mew=4, capsize=4, capthick=1, linestyle='none')
        plotAxes.set_xticks(X,qflattened,rotation=30)
        # plotAxes.set_yticklabels(rotation=60)
        if i == 2:
            plotAxes.set_xlabel(xlabel[j])
        if j == 0:
            plotAxes.set_ylabel(ylabel[i])
        plotAxes.plot(X, base, linestyle='dashed')
        
plt.suptitle(superTitle,fontsize = 'xx-large',weight = 'extra bold', y=0.92)
plt.savefig("VerifyBinder.png")