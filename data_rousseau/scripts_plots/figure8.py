
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

rc('text.latex', preamble=r'\usepackage{amsmath}')
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

def figure8(format='pdf'):
        
    cmap=custom_cmap.make_cmap_customized(Palette='mountain')
    cmap2=custom_cmap.make_cmap_customized(Palette='green')

    fig,[[ax1,ax2,ax3],[ax4,ax5,ax6]]=plt.subplots(2,3,figsize=(4.88, 6.))
    Axs=[ax1,ax2,ax3,ax4,ax5,ax6]
    for ax in Axs: 
        ax.xaxis.set_major_formatter(major_formatter)
        ax.yaxis.set_major_formatter(major_formatter)
        Pos=ax.get_position().bounds
        ax.set_position([Pos[0]+0.0,Pos[1]+0.05,0.97*Pos[2],1.05*Pos[3]])
        ax.set_ylim([-1.1,1.8])
    for ax in Axs[3:6]: 
        Pos=ax.get_position().bounds
        ax.set_position([Pos[0],Pos[1]-0.06,Pos[2],Pos[3]])
        
        ax.set_ylim([-1.1,0.55])
    Axs[1].set_yticklabels([])    

    Axs[2].set_yticklabels([])  
    Axs[4].set_yticklabels([])    
    Axs[5].set_yticklabels([])    

    PhyCh=data['physico-chemical_properties']

    Exp=['A1','A4','B1','B4']
    i=0
    for ca,kca in zip(Exp,range(len(Exp))):   
    #for i in [4,5,6,7]:  
    #for i in [0,1,2,3]:

        fvecZ=data[ca]['profiles']['z_80']
        Vel=data[ca]['profiles']['ux']
        uet=data[ca]['flow_characteristics']['uet']
        uetp=data[ca]['flow_characteristics']['uetp']
        dp=data[ca]['flow_characteristics']['d_p']
        H=data[ca]['flow_characteristics']['H']
        PR=data[ca]['profiles']['porosity']
        slope=data[ca]['flow_characteristics']['i']
        taut=data[ca]['profiles']['tautxz'] 
        indup=int(data[ca]['flow_characteristics']['free_surface_index'])
        indRC=int(data[ca]['flow_characteristics']['roughness_crest_index'])

    
        Axs[0].plot(data[ca]['profiles']['utx'][indup:]/(uet),fvecZ[indup:]/dp,color=cmap((i+1)/len(Exp)/1.5),linewidth=1.)
        Axs[0].plot(2.3*np.exp(-fvecZ[indup:]/H),fvecZ[indup:]/dp,'--',color=cmap((i+1)/len(Exp)/1.5),linewidth=1.,label=Exp[i]+r' - Eq. 5.2')    
        
        Axs[1].plot(data[ca]['profiles']['utz'][indup:] / (uet), fvecZ[indup:] / dp, color=cmap((i + 1) / len(Exp) / 1.5), linewidth=1.)
        Axs[1].plot(1.27*np.exp(-fvecZ[indup:]/H),fvecZ[indup:]/dp,'--',color=cmap((i+1)/len(Exp)/1.5),linewidth=1.,label=Exp[i]+r' - Eq. 5.3')    
    
        Axs[2].plot(data[ca]['profiles']['tautxz'][indup:]/PR[indup:]/(PhyCh['isoMix']['rho']*1000*uet**2),
        
        
        fvecZ[indup:]/dp,color=cmap((i+1)/len(Exp)/1.5),linewidth=1.)
        intgravity=np.zeros_like(Vel)
        dz=fvecZ[0]-fvecZ[1] 
        
        for j in range(indup,len(PR)):
            intgravity[j]=np.sum(9.81*slope*PR[indup:j])*dz  
        ax3.plot(-intgravity[indup:]*PhyCh['isoMix']['rho']*1000/(PhyCh['isoMix']['rho']*1000*uet**2),fvecZ[indup:]/dp,linestyle='dotted',color=cmap((i+1)/len(Exp)/1.5),label=Exp[i]+r' - $G$',linewidth=1.)

        Axs[3].plot(data[ca]['profiles']['udx'][indup:]/(uet),fvecZ[indup:]/dp,color=cmap((i+1)/len(Exp)/1.5),linewidth=1.)
        Axs[4].plot(data[ca]['profiles']['udz'][indup:]/(uet),fvecZ[indup:]/dp,color=cmap((i+1)/len(Exp)/1.5),linewidth=1.)
        Axs[5].plot(data[ca]['profiles']['tausxz'][indup:]/PR[indup:]/(PhyCh['isoMix']['rho']*1000*uet**2),fvecZ[indup:]/dp,color=cmap((i+1)/len(Exp)/1.5),label=Exp[i],linewidth=1.)
        

    #    Axs[5].plot(Prof[1][:,-1][indup:]/(PhyCh['isoMix']['rho']*1000*uet**2),fvecZ[indup:]/dp,color=cmap((i+1)/len(Exp)/1.5),label=Exp[i],linewidth=1.5)
    #     pormin=0.4
    #     sinbeta=((1-PR)/(1-pormin))**0.5
    #     cosbeta=(1-(1-PR)/(1-pormin))**0.5
    #     lambdap=0.5
    #     alpha=1/(lambdap*cosbeta+1-lambdap)
    #     tauttheo=PhyCh['isoMix']['rho']*1000*PR*sinbeta*alpha*lambdap*(1-alpha*cosbeta)*Vel**2
    # #    plt.plot(tauttheo,fvecZ/dp,'--',color=cmap((i+1)/len(Exp)/1.5))
    #     tauttheo2=PhyCh['isoMix']['rho']*1000*0.5*lambdap*(1-lambdap)*sinbeta**3*Vel**2

        i += 1
        if ca==Exp[0]:
            zrcA=fvecZ[indRC]/dp
        if ca==Exp[3]:
            zrcB=fvecZ[indRC]/dp
                    

    Axs[0].set_xlim([0.,3])
    for ix in [0,1,2,3,4,5]:
        Axs[ix].plot([-3, 3],[zrcA,zrcA], linestyle='--', linewidth=1, color='k', alpha=0.5)
        Axs[ix].plot([-3, 3], [zrcB,zrcB], linestyle='--', linewidth=1, color='k',alpha=0.5)
        Axs[ix].plot([-3,3],np.array([-0.7*dp,-0.7*dp])/dp,linestyle='dotted',linewidth=2,color='k') 
        Axs[ix].plot([-3,3],np.array([0,0])/dp,linewidth=1,color='k')   

    ExpV=['L12', 'L15']
    rhof =1.77e3
    i=0
    for ca in ExpV:
        rhof =1.77e3

        taudv=data[ca]['profiles']['tausxz']
        Zv=data[ca]['profiles']['z_80']
        uetV=data[ca]['uetV']
        Axs[0].plot(data[ca]['profiles']['utx'][:]/uetV,Zv/ 25e-3,'-.',color=cmap2(i*0.2+0.3),label=ca,linewidth=1.2,alpha=0.6) 
        Axs[1].plot(data[ca]['profiles']['utz'][:]/uetV,Zv/ 25e-3,'-.',color=cmap2(i*0.2+0.3),label=ca,linewidth=1.2,alpha=0.6) 
        Axs[2].plot(data[ca]['profiles']['tautxz'][:]/(rhof*uetV**2),Zv/ 25e-3,'-.',color=cmap2(i*0.2+0.3),label=ca,linewidth=1.2,alpha=0.6) 
        Axs[3].plot(data[ca]['profiles']['udx'][:]/uetV,Zv/ 25e-3,'-.',color=cmap2(i*0.2+0.3),label=ca,linewidth=1.2,alpha=0.6)  
        Axs[4].plot(data[ca]['profiles']['udz'][:]/uetV,Zv/ 25e-3,'-.',color=cmap2(i*0.2+0.3),label=ca,linewidth=1.2,alpha=0.6) 
        Axs[5].plot(data[ca]['profiles']['tausxz'][:]/(rhof*uetV**2),Zv/ 25e-3,'-.',color=cmap2(i*0.2+0.3),label=ca,linewidth=1.2,alpha=0.6) 
        i += 1

    Axs[1].set_xlim([-0.,1.5])

    xlim_max_2=0.1
    Axs[2].set_xlim([-1.,xlim_max_2])


    Axs[3].set_xlim([-0.,1.5])

    Axs[4].set_xlim([-0.,0.6])
    xlim_max_5=0.08
    Axs[5].set_xlim([-0.3,xlim_max_5])



    Axs[0].set_xlabel(r"$\sigma_x /  u_{*}$")
    Axs[0].set_ylabel(tl.rsd)

    Axs[3].set_ylabel(tl.rsd)
    Axs[1].set_xlabel(r"$\sigma_z /  u_*$")
    Axs[2].set_xlabel(r"$\langle \overline{u^\prime_x u^\prime_z} \rangle/  u_{*}^2$  ")
    Axs[3].set_xlabel(r"$\sqrt{\langle  \tilde{u}^2_x  \rangle}/  u_*$")
    Axs[4].set_xlabel(r"$\sqrt{\langle   \tilde{u}^2_z   \rangle}/  u_*$")
    Axs[5].set_xlabel(r"$\langle \tilde{u}_x \tilde{u}_z \rangle/  u_*^2$")
    Axs[0].yaxis.set_label_coords(-0.2, 0.5)
    Axs[3].yaxis.set_label_coords(-0.2, 0.5)
    #plt.plot(PR,fvecZ/dp,color=cmap((i+1)/len(Exp)/1.5),linestyle='dotted")
    # Axs[2].legend(fontsize=8,loc=3)
    Axs[5].legend(fontsize = 8, loc = 3)
    # Axs[0].legend(fontsize=8,loc=4)
    # Axs[1].legend(fontsize=8,loc=4)
    Dx, Dy = -0.07, 1.04



    Axs[0].text(Dx, Dy, '(a)', fontsize=9, horizontalalignment='center',verticalalignment='center', transform=Axs[0].transAxes)
    Axs[1].text(Dx, Dy, '(b)', fontsize=9, horizontalalignment='center',verticalalignment='center', transform=Axs[1].transAxes)
    Axs[2].text(Dx, Dy,'(c)', fontsize=9, horizontalalignment='center',verticalalignment='center', transform=Axs[2].transAxes)
    Axs[3].text(Dx, Dy, '(d)', fontsize=9, horizontalalignment='center',verticalalignment='center', transform=Axs[3].transAxes)
    Axs[4].text(Dx, Dy, '(e)', fontsize=9, horizontalalignment='center',verticalalignment='center', transform=Axs[4].transAxes)
    Axs[5].text(Dx, Dy,'(f)', fontsize=9, horizontalalignment='center',verticalalignment='center', transform=Axs[5].transAxes)





    Axs[2].text(xlim_max_2 + 0.16, zrcA, r'$z_{rc_A}$', fontsize=10, horizontalalignment='center', verticalalignment='center')
    Axs[2].text(xlim_max_2+0.16, zrcB-0.01, r'$z_{rc_B}$', fontsize=10, horizontalalignment='center',verticalalignment='center')
    Axs[5].text(xlim_max_5+0.07, zrcA, r'$z_{rc_A}$', fontsize=10, horizontalalignment='center', verticalalignment='center')
    Axs[5].text(xlim_max_5+0.07, zrcB, r'$z_{rc_B}$', fontsize=10, horizontalalignment='center',verticalalignment='center')
    Axs[2].text(xlim_max_2 + 0.16, -0.7, r'$z_{t}$', fontsize=10, horizontalalignment='center',verticalalignment='center')
    Axs[5].text(xlim_max_5 + 0.07,-0.7, '$z_{t}$', fontsize=10, horizontalalignment='center',verticalalignment='center')

    fig.savefig('figure08.'+format,dpi=1200)
    print('Figure 8 plotted')
if __name__ == "__main__":
    print("File one executed when ran directly")
    figure8()
else:
   print("File one executed when imported")


    # %%
