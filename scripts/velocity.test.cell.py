#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py
# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

from pylab import genfromtxt;  

mat = genfromtxt("../results/particles.test.cell.000.txt");
#pyplot.plot(mat0[:,0], mat0[:,1], label = "data0");

#gaussian_numbers = np.random.randn(1000)
plt.hist(mat[:,5], 20, facecolor='green')
plt.title("Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")

plt.show()

#fig = plt.gcf()

#plot_url = py.plot_mpl(fig, filename='mpl-basic-histogram')
