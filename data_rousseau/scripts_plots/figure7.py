"""
Created on 4 Feb  2021

@author: Gauthier Rousseau

"""
# %%
import sys
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import json
import custom_cmap
import tools as tl
from matplotlib.ticker import FuncFormatter
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
from matplotlib import rc

rc('text', usetex=True)
# matplotlib.use('Qt4Agg')


# make sure that you  type the correct path
data=tl.loadData()
# to display the possible keys type :
# data['A1'].keys()
# # to call the porosity profile for instance type
# data['A1']['profiles']['Porosity']
plt.close('all')
major_formatter = FuncFormatter(tl.my_formatter)

def flip(items, ncol):
    import itertools
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])


# %matplotlib qt5
#%%

def figure7(format='pdf'):

    cmap=custom_cmap.make_cmap_customized(Palette='mountain')
    cmap2=custom_cmap.make_cmap_customized(Palette='green')
    MT=['o','D','^','s','*']

    Tt=[]
    i=1
    ProfT=[]
    fig,[ax1,ax2]=plt.subplots(1,2,figsize=(4.88, 3.))
    Axs=[ax1,ax2]
    for ax in Axs: 
        ax.xaxis.set_major_formatter(major_formatter)
        ax.yaxis.set_major_formatter(major_formatter)
        Pos=ax.get_position().bounds
        ax.set_position([Pos[0]-0.04,Pos[1]+0.05,Pos[2]*1.05,Pos[3]*1.05])
        ax.set_ylim([-0.75,1])

    Xlims=[0.005,20]
    Pos=Axs[1].get_position().bounds
    Axs[1].set_position([Pos[0]+0.04,Pos[1],Pos[2],Pos[3]]) 

    Exps=[['A1','A2','A3','A4'],['B1','B2','B3','B4','B5']]
    for Exp,ax in zip(Exps,Axs):
        i=0
        for ca in Exp:  
            fvecZ=data[ca]['profiles']['z_80']
            Vel=data[ca]['profiles']['ux']
            uetp=data[ca]['flow_characteristics']['uetp']
            dp=data[ca]['flow_characteristics']['d_p']
            PR=data[ca]['profiles']['porosity']
            taut=data[ca]['profiles']['tautxz'] 
            indup=int(data[ca]['flow_characteristics']['free_surface_index'])
            indRC=int(data[ca]['flow_characteristics']['roughness_crest_index'])

            ax.semilogx(Vel[indup:]/uetp,fvecZ[indup:]/dp,color=cmap((i+1)/len(Exp)/1.5),linewidth=1.,label=Exp[i])
            i+=1

        ax.set_xlim(Xlims)
        ax.plot(Xlims,np.array([fvecZ[indRC],fvecZ[indRC]])/dp,linestyle='--',linewidth=1,color='k') 
        ax.plot(Xlims,np.array([0,0])/dp,linewidth=1,color='k')   
        ax.set_ylim([-1.5,1.8])

        ax.set_xlabel(r'$U_x /  u_{p}$')
        ax.legend(fontsize=8,loc=2,ncol=1)
        ax.yaxis.set_label_coords(-0.1,0.5)

    ax1.set_ylabel(tl.rsd)
    Axs[1].set_yticklabels([]) 
    dxz=30
    Axs[0].text(Xlims[1]+dxz, 0.28, '$z_{rc,A}$', fontsize=11, horizontalalignment='center',verticalalignment='center')
    Axs[1].text(Xlims[1]+dxz, 0.14, '$z_{rc,B}$', fontsize=11, horizontalalignment='center',verticalalignment='center')

    ax1.plot([0.6,0.6],[-0.7,-1.1],'-',color='k')
    ax1.text(1,-0.95,'$U_{x,SSL}$')

    ax2.plot([1.8,1.8],[-0.6,-1.],'-',color='k')
    ax2.text(3,-0.85,'$U_{x,SSL}$')
    # fig.text(0.01,-1.3,'(a)')
    # fig.text(0.01,-1.3,'(b)')

    fig.text(0.04,0.95,'(a)',fontsize=9)
    fig.text(0.5,0.95,'(b)',fontsize=9)

    fig.savefig('figure07.'+format,dpi=1200)
    print('Figure 7 plotted')
if __name__ == "__main__":
    print("File one executed when ran directly")
    figure7()
else:
   print("File one executed when imported")

