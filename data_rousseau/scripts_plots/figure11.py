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

# %matplotlib qt5

#%% figure 11

def figure11(format='pdf'):

    cmap=custom_cmap.make_cmap_customized(Palette='mountain')
    cmap2=custom_cmap.make_cmap_customized(Palette='green')
    Exp=['A1','A2','A3','A4']
    MT=['o','D','^','s','*']

    Tt=[]
    i=1
    ProfT=[]
    fig,[ax1,ax2]=plt.subplots(1,2,figsize=(4.88, 2.6))
    Axs=[ax1,ax2]
    Ylims=[-0.5,0.35]
    Xlims = [0.,4*(0.02)**0.5]
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

    Exps=[['A1','A2','A3','A4'],['B1','B2','B3','B4','B5']]
    for Exp,ax in zip(Exps,Axs):
        i=0
        
        for ca in Exp:   
            fvecZ=data[ca]['profiles']['z_80']
            Vel=data[ca]['profiles']['ux']
            PR=data[ca]['profiles']['porosity']
            taut=data[ca]['profiles']['tautxz']
            indup=int(data[ca]['flow_characteristics']['free_surface_index'])
            dz=fvecZ[0]-fvecZ[1]    
            dudz=np.abs(np.diff(Vel))/dz
            zlm=fvecZ[1:]+dz/2
            porMexp=np.mean(np.array([PR[1:],PR[:-1]]),axis=0)
            dp=data[ca]['flow_characteristics']['d_p']
            lpores = dp * ((1 - porMexp) /(1-0.38))**0.5
            tautexpM=np.mean([taut[1:],taut[:-1]],axis=0) 
            lm=np.sqrt((-tautexpM/(PhyCh['isoMix']['rho']*1000)))/(dudz)
            H=data[ca]['flow_characteristics']['H']
            ax.plot(lm[indup::3]/lpores[indup::3],zlm[indup::3]/dp,MT[i],color=cmap((i+1)/len(Exp)/1.5),label=Exp[i],ms=3)
            i+=1
        ax.set_xlim([-0.,0.4])
        ax.set_xlim([-0.,0.4])


        ax.plot([-0.,3],np.array([0,0]),linewidth=1,color='k') 

        uet=data[ca]['flow_characteristics']['uet']

        #Axs[1].set_xlim([-0.,1.5])
        # ax.plot(zlm*0.4/0.001,zlm*1000,'--',color=cmap2(0.3),label=r'$0.4\, z$')
        ax.set_xlabel(r'$\ell_t/(d_p \sqrt{\frac{1-\epsilon}{1-\epsilon_b}})$')

        ax.yaxis.set_label_coords(-0.1,0.5)

        #plt.plot(PR,fvecZ/dic
        # tExp['Dbeads'],color=cmap((i+1)/len(Cpl)/1.5),linestyle='dotted')
        ax.legend(fontsize=8,loc=4)
        ax.set_xlim(Xlims)
        indRC=int(data[ca]['flow_characteristics']['roughness_crest_index'])
        ax.plot(Xlims,np.array([fvecZ[indRC],fvecZ[indRC]])/dp,linestyle='--',linewidth=1,color='k') 
    Axs[0].set_ylabel(tl.rsd)

    Axs[0].text(Xlims[1]+0.06, 0.28, '$z_{rc,A}$', fontsize=11, horizontalalignment='center',verticalalignment='center')
    Axs[1].text(Xlims[1]+0.06, 0.14, '$z_{rc,B}$', fontsize=11, horizontalalignment='center',verticalalignment='center')

    xlamb=0.1
    Axs[0].plot([xlamb,xlamb],[-5,5],'-',c='k')
    Axs[1].plot([xlamb,xlamb],[-5,5],'-',c='k')


    Axs[0].text(xlamb, -0.57, r'$\lambda^{\prime}=0.1$', fontsize=8, horizontalalignment='center',verticalalignment='center')
    Axs[1].text(xlamb, -0.57, r'$\lambda^{\prime}=0.1$', fontsize=8, horizontalalignment='center',verticalalignment='center')

    fig.text(0.04,0.94,'(a)',fontsize=9)
    fig.text(0.5,0.94,'(b)',fontsize=9)
    fig.savefig('figure11.'+format,dpi=1200)
    print('Figure 11 plotted')
if __name__ == "__main__":
    print("File one executed when ran directly")
    figure11()
else:
   print("File one executed when imported")


#%%



# %%
