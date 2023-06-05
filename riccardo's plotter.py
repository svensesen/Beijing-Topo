import pandas as pd
import numpy as np
import os 
import matplotlib.plot as plt

def plotter(data:dict, marker = 'o-'):
  for key in data.keys():
    for point in data[key]:
      plt.scatter(point[0], point[1], marker = marker)

def data_plot(folder:None):

  
  if folder:
    path = os.listdir(folder)
  else:
    path = os.listdir()
  pt_dict = {}
  for fileplt in path:
    if fileplt.endswith('.plt'):
      data_i = np.genfromtxt(path + '/' + fileplt, delimiter=',', skip_header=6)
      pt_dict[fileplt] = [list(data[:3]) for data in data_i]

  plotter(pt_dict)

