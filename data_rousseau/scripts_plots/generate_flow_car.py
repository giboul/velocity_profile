#%%
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import tools as tl
from matplotlib.ticker import FuncFormatter
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
from matplotlib import rc


data=tl.loadData()

cases=['A1','A2','A3','A4','B1','B2','B3','B4','B5']

for c in cases:
    print(data[c]['flow_characteristics']['uet']*data[c]['flow_characteristics']['d_p']/(data['physico-chemical_properties']['isoMix']['visc']*1e-6))

cases=['L12','L14','L15']


for c in cases:
    d=data[c]['modellings']['LiAndSawamoto']['dictModel']['d']

    print(data[c]['uetV']*d/(1e-6))







# data=tl.convtolist(data)
# f = open('data/data_Rousseau_JFM3.json', "w")
# json.dump(data, f, indent=4)
# f.close()

# f = open('data/data_Rousseau_JFM3.json')
# data = json.load(f)
# f.close()
# data = tl.convNp(data)




#%%