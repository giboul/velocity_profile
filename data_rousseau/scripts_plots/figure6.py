
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

def figure6(format='pdf'):

    slopes=np.array([0.5,
    0.5,
    0.5,
    0.5,
    1,
    1,
    2,
    2,
    4,
    4,
    8,
    8])/100
    por1=0.38
    por2=0.33




    #% Second for PhD
    subVel1=np.array([0.0016632619247300405,
    0.0034833368152520743,
    0.004163071563776872,
    0.004299888242200031,
    0.0015324195272295373,
    0.004080125947974336,
    0.004658626694567437,
    0.007157885529367379,
    0.009898724991223169,
    0.00901799790586682,
    0.011635966707738843,
    0.014515363631529396,
    0.011277696871258926,
    0.0013800074685492683,
    0.0030054662554341526,
    0.00735352612734195,
    0.01324132360723663])

    slopes=np.array([0.005,
    0.005,
    0.005,
    0.01,
    0.01,
    0.01,
    0.01,
    0.02,
    0.02,
    0.02,
    0.04,
    0.04,
    0.04,
    0.005,
    0.01,
    0.02,
    0.04])

    ind=np.arange(len(slopes))

    np.std([0.0016632619247300405,
    0.0034833368152520743,
    0.004163071563776872,
    0.0013800074685492683])


    fig,[ax1,ax2]=plt.subplots(1,2,figsize=(4.88, 3.))
    Axs=[ax1,ax2]
    for ax in Axs: 
        ax.xaxis.set_major_formatter(major_formatter)
        ax.yaxis.set_major_formatter(major_formatter)
        Pos=ax.get_position().bounds
        ax.set_position([Pos[0]-0.02,Pos[1]-0.03,Pos[2]*1.05,Pos[3]*1.05])

    Pos=Axs[1].get_position().bounds
    Axs[1].set_position([Pos[0]+0.04,Pos[1],Pos[2],Pos[3]]) 



    cmap=custom_cmap.make_cmap_customized(Palette='mountain')
    cmap2=custom_cmap.make_cmap_customized(Palette='green')
    MT=['o','D','^','s','*']

    uvec=np.arange(0,0.04,0.0005)

    # ax3 = Axs[0].twiny()

    Axs[0].errorbar(slopes[ind],subVel1[ind], yerr=0.002*np.ones(len(slopes)),linewidth=0,marker='o', mfc='k',
            mec='k', ms=3, mew=1,elinewidth=0.8,ecolor=cmap2(0.7),capsize=4,label='Measured $U_{x,SSL}$')




    Axs[0].set_ylabel(r'$U_{x,SSL} ~[\mathrm{m/s}]$')
    Axs[0].set_xlabel(r'$i$ [-]')

    nu=3e-6/0.95
    por1=0.39
    por2=0.35
    dp1=0.008
    dp2=0.015
    g=9.81
    Ae=180




    Axs[0].plot(Ae*(1-por1)**2*nu/(por1**2*dp1**2*g)*uvec+1.75*(1-por1)/(por1*dp1*g)*uvec**2, uvec,'-.',color=cmap2(0.8),label=r'Ergun equation')

    Axs[0].plot(Ae*(1-por1)**2*nu/(por1**2*dp1**2*g)*uvec,uvec,'--',color=cmap2(0.3),label='Kozeny-Carman',lw=2)




    Axs[0].legend(fontsize=8,framealpha=0.4,loc=2)
    # ax3.set_ylim([0,0.018*0.008/3e-6])
    # ax3.grid()
    # ax3.set_xlabel(r'$Re_p$')
    Axs[0].set_ylim([0,0.022])
    Axs[0].set_xlim([0,0.05])

    Axs[0].xaxis.set_major_formatter(major_formatter)
    Axs[0].yaxis.set_major_formatter(major_formatter)

    subVel1=np.array([0.004578756167502515,
    0.006010113536591663,
    0.004992477092281067,
    0.004384371434580384,
    0.006515899101549208,
    0.012933242069416083,
    0.021647428157208232,
    0.037608783879505336,
    0.037929150123693786,
    0.007415051539814834])

    slopes=np.array([0.005, 0.005, 0.005, 0.005, 0.005, 0.02, 0.04, 0.08, 0.08, 0.01])

    ind=np.arange(len(slopes))




    uvec=np.arange(0,0.06,0.0005)

    # ax3 = Axs[1].twiny()



    Axs[1].errorbar(slopes[ind],subVel1[ind], yerr=0.002*np.ones(len(slopes)),linewidth=0,marker='o', mfc='k',
            mec='k', ms=3, mew=1,elinewidth=0.8,ecolor=cmap2(0.7),capsize=3,label='Measured $U_{x,SSL}$')



    Axs[1].set_ylabel(r'$U_{x,SSL}~ [\mathrm{m/s}]$')
    Axs[1].set_xlabel(r'$i$ [-]')

    nu=3e-6/0.95
    por1=0.35
    por2=0.35
    dp1=0.013
    dp2=0.015
    g=9.81
    Ae=180




    dp1=0.015
    dp2=0.013



    dp1=0.014
    # dp2=0.012

    Axs[1].plot(Ae*(1-por1)**2*nu/(por1**2*dp1**2*g)*uvec+1.75*(1-por1)/(por1*dp1*g)*uvec**2,uvec,'-.',color=cmap2(0.8),label=r'Ergun equation')

    Axs[1].plot(Ae*(1-por1)**2*nu/(por1**2*dp1**2*g)*uvec,uvec,'--',color=cmap2(0.3),label='Kozeny-Carman ',lw=2)



    Axs[1].legend(fontsize=8,framealpha=0.4,loc=2)
    # ax3.set_xlim([0,0.04*0.014/3e-6])
    # ax3.grid()
    # ax3.set_xlabel(r'$Re_p$')
    Axs[1].set_ylim([0,0.06])
    Axs[1].set_xlim([0,0.1])

    Axs[1].xaxis.set_major_formatter(major_formatter)
    Axs[1].yaxis.set_major_formatter(major_formatter)
    #AX[1].plot(subVel1[ind],slopes[ind],'x',ms=8,label='7-9 mm',color='k')
    # fig.set_size_inches(6.5, 3.5)
    fig.text(0.03,0.92,'(a)',fontsize=9)
    fig.text(0.51,0.92,'(b)',fontsize=9)

    for a in Axs:
        Pos=a.get_position().bounds
        a.set_position([Pos[0], Pos[1] + 0.1, Pos[2], Pos[3] - 0.12])
    Pos=a.get_position().bounds
    a.set_position([Pos[0]+0.02, Pos[1], Pos[2] , Pos[3]])
    fig.savefig('figure06.'+format,dpi=1200)
    print('Figure 6 plotted')

if __name__ == "__main__":
    print("File one executed when ran directly")
    figure6()
else:
   print("File one executed when imported")

# %%
