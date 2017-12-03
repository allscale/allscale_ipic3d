#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

import plotly.plotly as py
# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

from pylab import genfromtxt;  

mat = genfromtxt("../results/particles.test.1202.txt");

plt.plot(mat[:,0], mat[:,1], 'ro')
plt.xlabel("Value")
plt.ylabel("Frequency")

plt.show()

#fig = plt.gcf()

#plot_url = py.plot_mpl(fig, filename='mpl-basic-histogram')
