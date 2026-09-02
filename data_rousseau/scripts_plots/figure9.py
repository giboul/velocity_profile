# -*- coding: utf-8 -*-
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

#%% figure 9

def figure9(format='pdf'):
    cmap=custom_cmap.make_cmap_customized(Palette='mountain')
    cmap2=custom_cmap.make_cmap_customized(Palette='green')

    Tt=[]
    i=1
    ProfT=[]
    fig,[ax1,ax2]=plt.subplots(1,2,figsize=(4.88, 2.5))
    Axs=[ax1,ax2]
    Ylims=[-0.9,0.35]
    for ax in Axs: 
        ax.xaxis.set_major_formatter(major_formatter)
        ax.yaxis.set_major_formatter(major_formatter)
        Pos=ax.get_position().bounds
        ax.set_position([Pos[0]-0.04,Pos[1]+0.1,Pos[2]*1.05,Pos[3]*0.9])
        ax.set_ylim(Ylims)
    Axs[1].set_yticklabels([]) 
    Pos=Axs[1].get_position().bounds
    Axs[1].set_position([Pos[0]+0.04,Pos[1],Pos[2],Pos[3]]) 
    PhyCh=data['physico-chemical_properties']

    Exps=[['A1','A2','A3','A4'],['B1','B2','B3','B4']]
    for Exp,ax in zip(Exps,Axs):
        i=0
        
        for ca in Exp:   
            fvecZ=data[ca]['profiles']['z_80']
            Vel=data[ca]['profiles']['ux']
            Por=data[ca]['profiles']['porosity']
            taud=data[ca]['profiles']['tausxz']
            indup=int(data[ca]['flow_characteristics']['free_surface_index'])
            dz=fvecZ[0]-fvecZ[1]    
            dudz=np.abs(np.diff(Vel))/dz
            zlm=fvecZ[1:]+dz/2
            dp=data[ca]['flow_characteristics']['d_p']
            H=data[ca]['flow_characteristics']['H']
            sinbeta=np.sqrt((1-Por)/(1-np.min(Por)))
            cosbeta = np.sqrt(1 - (1 - Por) / (1 - np.min(Por)))
            lambdap=0.015
            dv=-np.diff(Vel)/np.absolute(fvecZ[0]-fvecZ[1])
            dvM=np.mean(np.array([dv[1:], dv[:-1]]), axis=0)
            dvM2 = np.zeros_like(Vel)
            dvM2[1:-1] = dvM
            TheoDisp = -Vel * PhyCh['isoMix']['rho'] * 1000 * sinbeta ** (2) * lambdap * dp * dvM2
            # ax.plot(lm[indup::3]/lpores[indup::3],zlm[indup::3]/dp,MT[i],color=cmap((i+1)/len(Exp)/1.5),label=Exp[i],ms=3)
            ax.plot(-taud[indup:],fvecZ[indup:]/dp,color=cmap((i+1)/len(Exp)/1.5),linewidth=1,label=Exp[i])
            ax.plot(-TheoDisp[indup:], fvecZ[indup:] / dp, label='(2.13)', linestyle='dotted', color=cmap((i + 1) / len(Exp) / 1.5), linewidth=1)
        
            i+=1
        Xlims= [-0.01,0.5]   
        ax.set_xlim(Xlims)
        indRC=int(data[ca]['flow_characteristics']['roughness_crest_index'])
        ax.plot(Xlims,np.array([fvecZ[indRC],fvecZ[indRC]])/dp,linestyle='--',color='k') 
        ax.plot(Xlims,np.array([0,0])/H,color='k') 
        ax.plot(Xlims,np.array([-0.5*dp,-0.5*dp])/dp,linestyle='dotted',color='k') 




        uet=data[ca]['flow_characteristics']['uet']

        #Axs[1].set_xlim([-0.,1.5])
        # ax.plot(zlm*0.4/0.001,zlm*1000,'--',color=cmap2(0.3),label=r'$0.4\, z$')
        ax.set_xlabel(r'$\tau_d ~[\mathrm{N~m^{-2}}]$')



        #plt.plot(PR,fvecZ/dic
        # tExp['Dbeads'],color=cmap((i+1)/len(Cpl)/1.5),linestyle='dotted')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(flip(handles, 2), flip(labels, 2),fontsize=8,loc=4,ncol=2)
        ax.set_xlim(Xlims)
        indRC=int(data[ca]['flow_characteristics']['roughness_crest_index'])
        ax.plot(Xlims,np.array([fvecZ[indRC],fvecZ[indRC]])/dp,linestyle='--',linewidth=1,color='k') 
    Axs[0].yaxis.set_label_coords(-0.12,0.5)
    Axs[0].set_ylabel(tl.rsd)
    dxz=0.06
    Axs[0].text(Xlims[1]+dxz, 0.28, '$z_{rc,A}$', fontsize=11, horizontalalignment='center',verticalalignment='center')
    Axs[1].text(Xlims[1]+dxz, 0.14, '$z_{rc,B}$', fontsize=11, horizontalalignment='center',verticalalignment='center')

    xlamb=0.
    Axs[0].plot([xlamb,xlamb],[-5,5],'-',c='k')
    Axs[1].plot([xlamb,xlamb],[-5,5],'-',c='k')

    fig.text(0.04,0.94,'(a)',fontsize=9)
    fig.text(0.5,0.94,'(b)',fontsize=9)
    fig.savefig('figure09.'+format,dpi=1200)
    print('Figure 9 plotted')
if __name__ == "__main__":
    print("File one executed when ran directly")
    figure9()
else:
   print("File one executed when imported")


#%%



# %%
