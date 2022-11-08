# -*- coding: utf-8 -*-
import numpy as np
import re
import csv

field = ['latSizeXMax','latSizeYMax','latType','interactionJprime',
         'interactionJprimeprime','externalZMageticField',
         'interactionType','temperature']
size = np.logspace(2,4,base=2,num=3,endpoint=True)
latType = ['Triangle']
interactionJprime = [-1.0]
interactionJprimeprime = [0]
magfield = np.linspace(0,1,21,endpoint=True)
interactionType = ['Heisenberg']
temperature = np.linspace(0.1,5,50,endpoint=True)

with open('./input.csv', 'w', newline='') as f:
    write = csv.writer(f)
    write.writerow(field)
    for para1 in size:
        for para2 in latType:
            for para3 in interactionJprime:
                for para4 in interactionJprimeprime:
                    for para5 in magfield:
                        for para6 in interactionType:
                            for para7 in temperature:
                                write.writerow([para1,para1,para2,para3,para4,para5,para6,para7])