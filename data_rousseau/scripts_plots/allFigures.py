

#%%
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))


from figure4 import figure4
from figure5 import figure5
from figure6 import figure6
from figure7 import figure7
from figure8 import figure8
from figure9 import figure9
from figure10 import figure10
from figure11 import figure11
from figures12_13_14_simu import figures12_13_14
import matplotlib.pyplot as plt
# plt.ion()
plt.ioff() #set interactive off to save the figure on your computer 
from matplotlib import rc
rc('text', usetex=True) #set to True if you have latex availaible on your computer
# import matplotlib
# matplotlib.use('Qt5Agg')
## 
frmt='eps'
#%% plot paper figures
figure4(format=frmt)
figure5(format=frmt)
figure6(format=frmt)
figure7(format=frmt)
figure8(format=frmt)
figure9(format=frmt)
figure10(format=frmt)
figure11(format=frmt)
figures12_13_14('A2',format=frmt)
figures12_13_14('B4',format=frmt)
figures12_13_14('L12',format=frmt)

#%%
# plot all simulation
# cases=['A1','A2','A3','A4','B1','B2','B3','B4','B5','L11','L12','L13','L14','L15']
# for c in cases:
#         figures12_13_14(selected_Case=c,format=frmt)



 
#%%
# from figures12_13_14_simu_contrib import figures_contrib
# cases = ['A1','A2','A3','A4','B1','B2','B3','B4','B5','L11','L12','L13','L14','L15']
# for c in cases:
#     for k in ['LiAndSawamoto', 'vD_LiAndSawamoto', 'dispersion_Li_vD']:
#         figures_contrib(selected_Case=c, closures=k)

# %%
