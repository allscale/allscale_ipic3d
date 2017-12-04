#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py

from pylab import genfromtxt;  

mat = genfromtxt("../results/density.tiny.1204.txt");

plt.plot(mat[:,0], mat[:,1], 'ro')
#plt.hist(mat[:,3], 20, facecolor='green')
plt.title("Density")
plt.xlabel("Value")
plt.ylabel("Frequency")

plt.show()

#fig = plt.gcf()

#plot_url = py.plot_mpl(fig, filename='mpl-basic-histogram')
