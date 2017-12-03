#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py

from pylab import genfromtxt;  

mat = genfromtxt("../results/particles.test.1201.txt");

#plt.plot(mat[:,3], mat[:,4], 'ro')
plt.hist(mat[:,3], 20, facecolor='green')
plt.title("Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")

plt.show()

#fig = plt.gcf()

#plot_url = py.plot_mpl(fig, filename='mpl-basic-histogram')
