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

# %% Load data

# make sure that you  type the correct path
f = open('../data_Rousseau_JFM.json')
data = json.load(f)
cases = ['A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3', 'B4', 'B5']
for ca in cases:
    for k in data[ca]['profiles'].keys():
        data[ca]['profiles'][k] = np.array(data[ca]['profiles'][k])
f.close()
# to display the possible keys type :
# data['A1'].keys()
# # to call the porosity profile for instance type
# data['A1']['profiles']['Porosity']
plt.close('all')
major_formatter = FuncFormatter(tl.my_formatter)

# %matplotlib qt5

#%% figure 10


def figure10(format='pdf'):

    cmap=custom_cmap.make_cmap_customized(Palette='mountain')
    cmap2=custom_cmap.make_cmap_customized(Palette='green')
    MT=['o','D','^','s','*']

    Tt=[]
    i=1
    ProfT=[]
    fig,[ax1,ax2]=plt.subplots(1,2,figsize=(4.88, 4))
    Axs=[ax1,ax2]
    for ax in Axs: 
        ax.xaxis.set_major_formatter(major_formatter)
        ax.yaxis.set_major_formatter(major_formatter)
        Pos=ax.get_position().bounds
        ax.set_position([Pos[0]-0.04,Pos[1],Pos[2]*1.05,Pos[3]*1.05])
        ax.set_ylim([-0.75,1.1])
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
            tautexpM=np.mean([taut[1:],taut[:-1]],axis=0) 
            lm=np.sqrt((-tautexpM/(PhyCh['isoMix']['rho']*1000)))/(dudz)
            H = data[ca]['flow_characteristics']['H']
            ax.plot(lm[indup::3]/H,zlm[indup::3]/H,MT[i],color=cmap((i+1)/len(Exp)/1.5),label=Exp[i],ms=3)
            Xlims= [-0.01,0.5]
            # indRC=int(data[ca]['flow_characteristics']['roughness_crest_index'])

            # ax.plot(Xlims,np.array([fvecZ[indRC],fvecZ[indRC]])/H,linestyle='--',linewidth=1,color='k') 
            i+=1
        ax.set_xlim([-0.,0.4])
        ax.set_xlim([-0.,0.4])


        ax.plot([-0.,3],np.array([0,0]),linewidth=1,color='k') 
        ksi=np.arange(-0.,1,0.001)
        Cole=0.1
        vD=1
        uet=data[ca]['flow_characteristics']['uet']
        vD=1-np.exp(-uet*ksi*H/(PhyCh['isoMix']['visc']*10**-6)/26)
        ColeCurve0=0.41*(1-ksi)**0.5*(1/ksi+np.pi*0*np.sin(np.pi*ksi))**(-1)
        ColeCurve1=0.41*(1-ksi)**0.5*(1/ksi+np.pi*Cole*np.sin(np.pi*ksi))**(-1)
        ColeCurve2=0.41*(1-ksi)**0.5*(1/ksi+np.pi*Cole*np.sin(np.pi*ksi))**(-1)*vD
        ax.plot(ColeCurve0,ksi,'--',linewidth=1,color='k',label='Coles - $\Pi$=0') 
        ax.plot(ColeCurve1,ksi,'-.',linewidth=1,color='k',label='Coles - $\Pi$=0.1')  
        ax.plot(ColeCurve2,ksi,linewidth=1,color='k',label=r'(4.7)') 
        indA=300
        ax.annotate(r'$\Pi = 0.1 $ and $ \Gamma$', xy=(ColeCurve2[indA], ksi[indA]), xytext=(ColeCurve1[indA]+0.1, ksi[indA]+0.1),
                    arrowprops=dict(color='k',arrowstyle="->", connectionstyle="arc3"))
        indA=500
        ax.annotate(r'$\Pi = 0.1 $', xy=(ColeCurve1[indA], ksi[indA]), xytext=(ColeCurve1[indA]+0.14, ksi[indA]+0.1),
                    arrowprops=dict(color='k',arrowstyle="->", connectionstyle="arc3"))

        indA=800
        ax.annotate(r'$\Pi = 0$', xy=(ColeCurve0[indA], ksi[indA]), xytext=(ColeCurve1[indA]+0.1, ksi[indA]+0.1),
                    arrowprops=dict(color='k',arrowstyle="->", connectionstyle="arc3"))
        #Axs[1].set_xlim([-0.,1.5])
        ax.plot(zlm*0.4/0.001,zlm*1000,'--',color=cmap2(0.3),label=r'$0.4\, z$')
        ax.set_xlabel(r'$\ell_t/h_f$')

        ax.yaxis.set_label_coords(-0.1,0.5)

        #plt.plot(PR,fvecZ/dic
        # tExp['Dbeads'],color=cmap((i+1)/len(Cpl)/1.5),linestyle='dotted')
 
        ax.legend(fontsize=8,loc=4)

        
    fig.text(0.04,0.94,'(a)',fontsize=9)
    fig.text(0.5,0.94,'(b)',fontsize=9)

    Axs[0].set_ylabel(tl.rsh)
    Axs[0].yaxis.set_label_coords(-0.13,0.4)
    fig.savefig('figure10.'+format,dpi=1200)
    print('Figure 10 plotted')
if __name__ == "__main__":
    print("File one executed when ran directly")
    figure10()
else:
   print("File one executed when imported")

#%%



# %%
