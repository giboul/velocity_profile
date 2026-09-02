# -*- coding: utf-8 -*-
"""
Created on 4 Feb  2021

@author: Gauthier Rousseau

"""
# %%
import sys
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import custom_cmap
import tools as tl
from matplotlib.ticker import FuncFormatter
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
from matplotlib import rc

rc('text', usetex=True)
# matplotlib.use('Qt5Agg')

# make sure that you  type the correct path

data=tl.loadData()
major_formatter = FuncFormatter(tl.my_formatter)

# to display the possible keys type :
# data['A1'].keys()
# # to call the porosity profile for instance type
# data['A1']['profiles']['Porosity']
#%%
# %matplotlib qt5
plt.close('all')

def figure5(format='pdf'):
    cmap1=custom_cmap.make_cmap_customized(Palette='mountain')
    cmap2=custom_cmap.make_cmap_customized(Palette='green')
    fig, Axs=plt.subplots(1,2,figsize=(4.88, 2.5))

    for ax in Axs: 
        ax.xaxis.set_major_formatter(major_formatter)
        ax.yaxis.set_major_formatter(major_formatter)
        Pos=ax.get_position().bounds
        ax.set_position([Pos[0]-0.04,Pos[1]+0.05,Pos[2]*1.05,Pos[3]*1.05])
        ax.set_ylim([-1.7, 2.7])
        ax.set_yticks([-1, 0, 1, 2])
    Axs[1].set_yticklabels([]) 
    Pos=Axs[1].get_position().bounds
    Axs[1].set_position([Pos[0]+0.04,Pos[1],Pos[2],Pos[3]]) 
    Axs[0].set_ylabel(tl.rsd)



    Exp=['L11','L12','L13','L14','L15']

    for i in [1, 3, 4]:
        ca=Exp[i]
        dp=0.025
        uet=data[ca]['uetV']
        indup=0
        Axs[0].plot(data[ca]['profiles']['ux'][indup:],data[ca]['profiles']['z_80'][indup:]/dp,'-.',c=cmap2(0.1*i+0.2),label=Exp[i],linewidth=1)
        Axs[1].plot(data[ca]['profiles']['ux'][indup:]/uet,data[ca]['profiles']['z_80'][indup:]/dp,'-.',c=cmap2(0.1*i+0.2),label=Exp[i],linewidth=1)
        

    Exp=['A1','A2','A3','A4','B1','B2','B3','B4','B5']


    for i in [0, 3, 5]:
        ca=Exp[i]
        dp=data[ca]['flow_characteristics']['d_p']
        uet=data[ca]['flow_characteristics']['uet']
        indup=int(data[ca]['flow_characteristics']['free_surface_index'])+2
        Axs[0].plot(data[ca]['profiles']['ux'][indup:],data[ca]['profiles']['z_80'][indup:]/dp,'-',c=cmap1(0.1*i+0.2),label=Exp[i],linewidth=1)
        Axs[1].plot(data[ca]['profiles']['ux'][indup:]/uet,data[ca]['profiles']['z_80'][indup:]/dp,'-',c=cmap1(0.1*i+0.2),label=Exp[i],linewidth=1)


    Axs[0].set_xlabel(r'$U_x$ in [$\mathrm{m~s}^{-1}$]')
    Axs[1].set_xlabel(r'$U_x/u_*$')
    Axs[1].legend()
    # fig.text(0.04, 0.9, '(a)')
    # fig.text(0.5,0.9,'(b)')

    fig.text(0.04,0.95,'(a)',fontsize=9)
    fig.text(0.5,0.95,'(b)',fontsize=9)

    # plt.show()
    fig.savefig('figure05.'+format,dpi=1200)
    print('Figure 5 plotted')
if __name__ == "__main__":
    print("File one executed when ran directly")
    figure5()
else:
   print("File one executed when imported")

#%%