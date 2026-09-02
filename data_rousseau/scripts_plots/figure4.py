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

# %% Figure 2
def figure4(format='pdf'):
    # plotted runs
    cases = [r'A2', r'A3']
    cmap1 = custom_cmap.make_cmap_customized(Palette='green')  # optional
    cmap2 = custom_cmap.make_cmap_customized(Palette='green')
    colors = [cmap1(0.4), cmap1(0.8)]
    ic = 0
    lws = [1, 1]
    fig, ax1, ax2, ax3, ax5, ax6 = tl.setPLotProfilFig2(
        namefig='profil'+str(0), major_formatter=major_formatter)

    PhyCh = data['physico-chemical_properties']
    for ca in cases:
        cl = colors[ic]
        data_case = data[ca]
        Vel = data_case['profiles']['ux']
        Z = data_case['profiles']['z_80']
        d = data_case['flow_characteristics']['d_p']
        slope = data_case['flow_characteristics']['i']
        indup = int(data_case['flow_characteristics']['free_surface_index'])+4
        indRC = int(data_case['flow_characteristics']['roughness_crest_index'])
        ax1.plot(Vel[indup:], Z[indup:]/d, '-', color=cl,
                label=ca, ms=5, linewidth=lws[ic])
        ax5.semilogx(Vel[indup:], Z[indup:]/d, '-', color=cl,
                    label=ca, ms=5, linewidth=lws[ic])
    #    ax5.plot(Vel,Z/d,'-',color=cl,label=Exps[i],ms=5,linewidth=lws[i])
        tautexp = data_case['profiles']['tautxz']
        taudexp = data_case['profiles']['tausxz']
        PR = data_case['profiles']['porosity']
        dz = Z[0]-Z[1]
        dudz = np.abs(np.diff(Vel))/dz
        porMexp = np.mean(np.array([PR[1:], PR[:-1]]), axis=0)
        zlm = Z[1:]+dz/2
        if ic == 1:
            ax2.plot(PR, Z/d, '-', color=cl, label=ca, ms=5, linewidth=lws[ic])
            ax6.plot(-tautexp[indup:], Z[indup:]/d, '--', color=cl,
                    label=ca+r' - $\tau_{t,exp}$', ms=5, linewidth=lws[ic])
            ax6.plot(-taudexp[indup:], Z[indup:] / d, '-', color=cl,
                    label=ca + r' - $\tau_{d,exp}$', ms=5, linewidth=lws[ic])
            ax6.plot(PhyCh['isoMix']['visc']*1e-3*dudz[indup:]*porMexp[indup:], zlm[indup:] /
                    d, '-.', color=cl, label=ca+r' - $\tau_{v,exp}$', ms=5, linewidth=lws[ic])
        if ic == 0:
            ax3.plot(-tautexp[indup:], Z[indup:]/d, '--', color=cl,
                    label=ca+r' - $\tau_{t,exp}$', ms=5, linewidth=lws[ic])
            ax3.plot(-taudexp[indup:], Z[indup:] / d, '-', color=cl,
                    label=ca + r' - $\tau_{d,exp}$', ms=5, linewidth=lws[ic])
            ax3.plot(PhyCh['isoMix']['visc']*1e-3*dudz[indup:]*porMexp[indup:], zlm[indup:] /
                    d, '-.', color=cl, label=ca+r' - $\tau_{v,exp}$', ms=5, linewidth=lws[ic])

        intgravity = np.zeros_like(Vel)
        for j in range(indup, len(PR)):
            intgravity[j] = np.sum(9.81*slope*PR[indup:j])*dz
        tautexpM = np.mean(np.array([tautexp[1:], tautexp[:-1]]), axis=0)
        ax1.errorbar(Vel[indup], Z[indup]/d,  xerr=0.05, color=(188/255, 44/255, 44/255), linewidth=0, marker='+', mfc=cmap2(0.4),
                    mec=(188 / 255, 44 / 255, 44 / 255), ms=5, mew=1, elinewidth=1, ecolor=(188 / 255, 44 / 255, 44 / 255), capsize=4)
        if ic == 0:
            ax3.plot(PhyCh['isoMix']['visc']*1e-3*dudz[indup-1:]*porMexp[indup-1:]-tautexp[indup:]-taudexp[indup:],
                    Z[indup:]/d, '-', color=(188/255, 44/255, 44/255), label=ca+r' - $\tau_{Tot}$', ms=5, linewidth=1.5)
        if ic == 1:
            ax6.plot(PhyCh['isoMix']['visc']*1e-3*dudz[indup-1:]*porMexp[indup-1:]-tautexp[indup:]-taudexp[indup:],
                    Z[indup:]/d, '-', color=(188/255, 44/255, 44/255), label=ca+r' - $\tau_{Tot}$', ms=5, linewidth=1.5)
        # k+=1
        if ic == 0:
            ax3.plot(intgravity[indup:] * PhyCh['isoMix']['rho'] * 1000, Z[indup:] / d,
                    linestyle='dotted', color=cl, label=ca + r' - $G$', linewidth=lws[ic])
        if ic == 1:
            ax6.plot(intgravity[indup:] * PhyCh['isoMix']['rho'] * 1000, Z[indup:] / d,
                    linestyle='dotted', color=cl, label=ca + r' - $G$', linewidth=lws[ic])

        ic += 1

    ax3.legend( loc=4, ncol=1)
    ax6.legend( loc=4, ncol=1)

    ax1.legend( loc=2)
    xM = 0.48
    ax1.set_xlim([-0.015, xM])
    ax1.plot([-0.015, xM], np.array([Z[indRC], Z[indRC]]) /
            d, linestyle='--', linewidth=0.6, color='k')
    ax1.plot([-0.015, xM], np.array([0, 0])/d, linewidth=0.4, color='k')

    ax1.plot([-0.015, 0.5], np.array([-0.7, -0.7]),
            linestyle='--', linewidth=0.6, color='k')

    # ax2.plot([-0.1, 1.1], np.array([Z[indRC], Z[indRC]]) /
    #         d, linestyle='--', linewidth=0.4, color='k')
        
    # ax2.plot([-0.1, 1.1], np.array([-0.7, -0.7]), linestyle='--', linewidth=0.4, color='k')

    xMb = 1.6
    for ax in [ax3, ax2 , ax6]:
        ax.set_xlim([-0.1, xMb])
        ax.plot([-0.1, xMb], np.array([Z[indRC], Z[indRC]]) /
                d, linestyle='--', linewidth=0.6, color='k')
        ax.plot([-0.1, xMb], np.array([-0.7, -0.7]),
                linestyle='--', linewidth=0.6, color='k')
        ax.plot([-0.1, xMb], np.array([0, 0])/d, linewidth=0.6, color='k')
        ax.plot([0, 0], np.array([-10, 10]) / d, linewidth=0.6, color='k')
    ax2.set_xlim([0., 1.1])


    ax6.text(xMb + 0.01, Z[indRC] / d, r'$z_{rc}$')
    ax6.text(xMb + 0.01, -0.7, r'$z_{t}$')
    ax5.yaxis.set_label_coords(-0.2, 0.4)
    fig.set_size_inches([4.88, 4.])
    fig.savefig('figure04.'+format,dpi=1200)
    print('Figure 4 plotted')
if __name__ == "__main__":
    print("File one executed when ran directly")
    figure4()
else:
   print("File one executed when imported")


# %%
