
# -*- coding: utf-8 -*-
"""
Created on 4 Feb  2021

@author: Gauthier Rousseau

"""
# %%
import sys
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import custom_cmap as cm
import tools as tl
from matplotlib.ticker import FuncFormatter
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
from matplotlib import rc

rc('text', usetex=True)
# matplotlib.use('Qt4Agg')

# make sure that you  type the correct path

data=tl.loadData(jsonF='../data_Rousseau_JFM.json')
# data=tl.loadData()
data=tl.convNp(data)
major_formatter = FuncFormatter(tl.my_formatter)

# to display the possible keys type :
# data['A1'].keys()
# # to call the porosity profile for instance type
# data['A1']['profiles']['Porosity']
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
#%%
# %matplotlib qt5
plt.close('all')
def figures_contrib(selected_Case='B5', closures='dispersion_Li_vD', format='pdf'):
    # cases='rousseau'
    sC=selected_Case

    if sC == 'L11' or sC == 'L12' or sC == 'L13' or sC == 'L14' or sC == 'L15':
        cases='voermans'
    else:
        cases='rousseau'

    print('Cases : ' + cases)
    print('Plot case : ' + sC)

    rc('text', usetex=True)
    z80 = data[sC]['profiles']['z_80']
    # for k in data.keys():
    #     print(data[k].keys())
    i_bottom =0
   
    # i_mixinglength = 4
    # i_trubstress = 3
    # i_dispstress = 2
    i_forces = 1
    i_velocity = 0

    dictModel = tl.dictModelGen('BAE', data)

    z80 = data[sC]['profiles']['z_80']
    if cases=='rousseau':
        Por_raw = data[sC]['profiles']['porosity']
        flowC = data[sC]['flow_characteristics']
        d = flowC['d_p']
        indup = int(flowC['free_surface_index'])
        indRC = int(flowC['roughness_crest_index'])
        Por = tl.porGen(z80, Por_raw, indRC, sC, d)
        zsurf = z80[indup]
        uetp = data[sC]['flow_characteristics']['uetp']
        rhof = data['physico-chemical_properties']['isoMix']['rho']*1e3
        nu = 3.0137130801687765e-06


    elif cases=='voermans':
        Por = data[sC]['profiles']['porosity']
        indup=0
        zsurf = z80[indup]
        d = 0.025
        i=3.08e-05
        uetp =(9.81*d*i)**0.5
        rhof =1.77e3 
        nu = 1.35e-6
        data[sC]['flow_characteristics']={}
        data[sC]['flow_characteristics']['i'] = i
        data[sC]['flow_characteristics']['uetp']=uetp
        data[sC]['flow_characteristics']['d_p'] = d
        data[sC]['flow_characteristics']['rhof'] = rhof
        flowC = data[sC]['flow_characteristics']
        uetdepth = (9.81*zsurf*d*flowC['i'])**0.5


    # % plots

    from matplotlib.ticker import FuncFormatter
    def my_formatter(x, pos):
        """Format 1 as 1, 0 as 0, and all values whose absolute values is between
        0 and 1 without the leading "0." (e.g., 0.7 is formatted as .7 and -0.4 is
        formatted as -.4)."""
        val_str = '{:g}'.format(x)
        if np.abs(x) > 0 and np.abs(x) < 1:
            return val_str.replace("0", "", 1)
        else:
            return val_str
    major_formatter = FuncFormatter(my_formatter)
    def setPLotProfilComp(namefig='profil'):
        major_formatter = FuncFormatter(my_formatter)
        exty=0.19
        extx=0.8
        my, sy = 0.1, 0.08
        mx = 0.12
        Nplots = 2
        exty=0.18*4/Nplots
        axiswindows = []
        AX = []
        mxlet=mx/4
        fig = plt.figure(namefig, figsize=(4.88, 8*Nplots/4))
        letters = ['(a)', '(b)', '(c)', '(d)', '(e)']
        for k in range(Nplots):
            aw=[mx, my + (k ) * exty + (k ) * sy, extx, exty]
            axiswindows.append(aw)
            a=plt.Axes(fig, aw)
            AX.append(a)
            fig.add_axes(a)
            a.xaxis.set_major_formatter(major_formatter)
            a.yaxis.set_major_formatter(major_formatter)
            fig.text(mxlet, my+(k+1)*exty+(k)*sy,   letters[Nplots-k-1])
            # Pos=ax.get_position().bounds
        # ax.set_position([Pos[0]-0.04,Pos[1],Pos[2]*1.05,Pos[3]*1.05])
        # for a in AX[:3]:
        #     a.set_xticklabels([])  
        # AX[i_trubstress].set_ylabel(r'$-\langle u^\prime_x u^\prime_z \rangle/  u_p^2$  ')
        # AX[i_mixinglength].set_ylabel(r'$\ell_t /d_p$')
        AX[i_forces].set_ylabel(r'$f_{contrib} /(\epsilon \rho g i)$')
        # AX[i_dispstress].set_ylabel(r'$-\langle \tilde{u}_x \tilde{u}_z \rangle/  u_p^2$')
        AX[i_velocity].set_ylabel(r'$ U_x / u_p$')
        AX[i_bottom].set_xlabel(tl.rsd)

        
        
        # fig.text(mxlet, my+2*exty+1*sy, '(c)')
        # fig.text(mxlet, my+1*exty, '(d)')
                    
        return fig, AX



    fig, AX = setPLotProfilComp(sC)

    porC=0.8

    indzinf=0
    indzsup = 2

    for k in data[sC]['modellings'].keys():
        z = data[sC]['modellings'][k]['z']

    for a in AX:
        a.set_xlim([-0.8,z80[indup]/d+0.2])
    cm1 = cm.make_cmap_customized(Palette='mountain')
    cm2 = cm.make_cmap_customized(Palette='green')

    ### plot expe

    scatter =AX[i_velocity].scatter(z80[indup::2]/d,data[sC]['profiles']['ux'][indup::2]/uetp,s=5,color=cm2(0.7))
    # AX[1].plot(z80[indup:]/d,-data[sC]['profiles']['tautxz'][indup:],'-.',linewidth=1,c=cm2(0.7))
    # l0=AX[i_trubstress].scatter(z80[indup::2]/d,-data[sC]['profiles']['tautxz'][indup::2]/(rhof*uetp**2),s=5,color=cm2(0.7))

    # AX[i_dispstress].plot(z80[indup:] / d, -data[sC]['profiles']['tausxz'][indup:], '-.', linewidth=1, c=cm2(0.7))
    # scatter =AX[i_dispstress].scatter(z80[indup::2] / d, -data[sC]['profiles']['tausxz'][indup::2]/(rhof*uetp**2),s=5,color=cm2(0.7))

    dz=z80[0]-z80[1]
    Vel = data[sC]['profiles']['ux']
        
    dudz=np.abs(np.diff(Vel))/dz

    porMexp=np.mean(np.array([Por[1:],Por[:-1]]),axis=0)
    #    ViscTerm=dictModel['nu']*porMexp*dudz 
    intgravity=np.zeros_like(porMexp)
    zlm=z80[1:]+dz/2
    for j in range(indup,len(porMexp)):
        intgravity[j]=np.sum(9.81*flowC['i']*porMexp[indup:j])*dz
    tautexp=data[sC]['profiles']['tautxz']
    tautexpM=np.mean(np.array([tautexp[1:],tautexp[:-1]]),axis=0)
    lm=np.sqrt((-tautexpM/dictModel['rhof'])/porMexp)/(dudz)
    # AX[2].plot(zlm[indup::2] / d, lm[indup::2], '+', color=cm2(0.7), label=r'$l_{m,exp}$', ms=5)
    # scatter = AX[i_mixinglength].scatter(zlm[indup::2] / d, lm[indup::2]/ d, s=5, color=cm2(0.7))


    #### plot Model ###
    h1, l1 = [scatter], ["Exp. - "+sC]
    ic = 0
    lw = [.8, 0.9, 0.7]
    ls = ['--', 'dotted', '-']

    for k in [ closures]:
        z = data[sC]['modellings'][k]['z']
        lz=len(z)
        Indz=range(indzinf+1,lz-2-indzsup)
        por = data[sC]['modellings'][k]['por']
        dz=np.diff(z)[0]
        hsup = z[np.where(por > porC)[0][0]]
        Zr = (z[Indz] - hsup) / d
        v=data[sC]['modellings'][k]['v']
        dv=np.diff(v)/dz
        dvM=np.mean(np.array([dv[1:],dv[:-1]]),axis=0)
        AX[i_velocity].plot(Zr, v[Indz]/uetp, c=cm1((ic+1) / 3),lw=lw[ic],ls=ls[ic])
        T=data[sC]['modellings'][k]['T']
        taut = rhof * T[1:-1] * (dvM ** 2)
        dtaut=np.diff(taut)/dz
        dtaut=(dtaut[:-1]+dtaut[1:])/2
        dtaut=smooth(dtaut,6)
        Afi=data[sC]['modellings'][k]['Afi']
        AfiM=np.mean(np.array([Afi[1:],Afi[:-1]]),axis=0) 
        taufi = rhof * Afi[1:-1] * v[1:-1] * dvM
        # dtaufi=np.diff(taufi)/dz
        dAfi=np.diff(Afi)/dz
        vM=np.mean(np.array([v[1:],v[:-1]]),axis=0)
        elem1=(dAfi * vM + dv * AfiM)*dv
        dtaufi=rhof*(elem1[:-1]+Afi[1:-1]*v[1:-1]*np.diff(dv)/dz)
        dtaufi=smooth(dtaufi,6)
        # dtaufi=rhof*(Afi[1:-1]*v[1:-1]*np.diff(dv)/dz)
        # dtaufi=rhof*(elem1[:-1])
        # dtaufi=(dtaufi[:-1]+dtaufi[1:])/2
        Indz2=range(indzinf,lz-3-indzsup)
        # AX[i_trubstress].plot(Zr, taut[Indz2]/(rhof*uetp**2), c=cm1((ic+1) / 3),lw=lw[ic],ls=ls[ic])
        Indz3=range(indzinf,lz-5-indzsup)
        ZrM=np.mean(np.array([Zr[1:],Zr[:-1]]),axis=0)
        G=data[sC]['modellings'][k]['G']*dictModel['rhof']
        [line2D]=AX[i_forces].plot(Zr[1:-1],(dtaut[Indz3])/G[Indz3],color=cm1((ic+1) / 3),label=r'$\frac{\mathrm{d} \tau_t}{\mathrm{d} z} $',lw=lw[ic],ls=ls[ic])
        [line2D]=AX[i_forces].plot(Zr[1:-1],(dtaufi[Indz3])/G[Indz3],color=cm1((ic+1) / 3),label=r'$\frac{\mathrm{d} \tau_d}{\mathrm{d} z} $',lw=lw[ic],ls='-')
        tau_v=nu*rhof*np.diff(por*v)/dz
        tau_v=(tau_v[:-1]+tau_v[1:])/2
        dtau_v=np.diff(tau_v)/dz

        dtau_vM=(dtau_v[:-1]+dtau_v[1:])/2
        dtau_vM=smooth(dtau_vM,6)
        [line2D]=AX[i_forces].plot(Zr[1:-1],(dtau_vM[Indz3])/G[Indz3],color=cm1((ic+1) / 3),lw=lw[ic],ls=':',label=r'$\frac{\mathrm{d} \tau_v}{\mathrm{d} z} $',)
        dpor=np.diff(por)/dz
        dporM=np.mean(np.array([dpor[1:],dpor[:-1]]),axis=0)    
        SecondBrickCorr=nu*rhof*dvM*dporM
        [line2D]=AX[i_forces].plot(Zr,SecondBrickCorr[Indz2]/G[Indz2],color='b',label=r'$-\varrho \nu \frac{\mathrm{d} U_x}{\mathrm{d} z} \frac{\mathrm{d} \epsilon}{\mathrm{d} z} $',lw=lw[ic],ls=':')
        lM=data[sC]['modellings'][k]['lM']

        fp=data[sC]['modellings'][k]['fp']
        fv=data[sC]['modellings'][k]['fv']
        f=rhof*(-(dictModel['AEr']*nu*(1-por[1:-1])**2/(por[1:-1]*d**2)*v[1:-1]+dictModel['BEr']*(1-por[1:-1])/(d)*v[1:-1]**2)+nu*np.diff(dpor)/dz*v[1:-1])
        [line2D]=AX[i_forces].plot(Zr,(fp[Indz]+fv[Indz])/G[Indz],color='k',label=r'$ \epsilon \bar{f}$',lw=lw[ic],ls='dotted')
        # [line2D]=AX[i_forces].plot(Zr,f[Indz]/G[Indz],color='k',label=r'$ f$',lw=lw[ic],ls='dotted')

        #plot the sum
        # sumF=fp[2:-2]+fv[2:-2]+dtau_vM+dtaut+SecondBrickCorr[1:-1]
        sumF=fp[2:-2]+fv[2:-2]+dtaut+dtaufi[1:-1]+dtau_vM-SecondBrickCorr[1:-1]
    
        # [line2D]=AX[i_forces].plot(Zr[1:-1],(sumF[Indz3])/G[Indz3],color='k',label=r'$\Sigma$',lw=1,ls='-')

        # [line2D]=AX[i_mixinglength].plot(Zr,lM[Indz]/d,color=cm1((ic+1) / 3),label=r'$ l_{m,mod}$',lw=lw[ic],ls=ls[ic])

        # [line2D]=AX[i_dispstress].plot(Zr, taufi[Indz2]/(rhof*uetp**2), c=cm1((ic + 1) / 3), lw=lw[ic], ls=ls[ic])  
        # h, l = [line2D],[leg[ic]]
        # h1, l1 = h1+h, l1+l  
        ic += 1

    # AX[0].plot(z80[indup::2] / d, data[sC]['profiles']['ux'][indup::2], '-.', linewidth=1, c=cm2(0.7))

    # l1b=AX[1].plot(z/d,z*0.4/d,'-',label=r'$0.4\, z$',lw=0.6,ls='-',c='k')

        
        fig.legend(h1 , l1, fontsize=7,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=4,bbox_to_anchor=(0.1, 0.95, 0.4, 0.2))
        # AX[i_forces].set_ylim([0, 5])
        # AX[i_mixinglength].set_ylim([0, 0.004/d])
        # if cases=='voermans':
            # AX[i_mixinglength].set_ylim([0, 0.01/d])
        for a in AX:
            a.yaxis.set_label_coords(-0.07, 0.5)

        # plt.pause(0.1)
        AX[i_forces].legend()

        # fig.savefig('test2.pdf')
        AX[i_forces].plot([-10,10],[0,0],c=(0.1,0.1,0.1),lw=0.4,zorder=-1)

        leg = [r'L \& S', r'vD + L \& S', r' Disp. + vD + L \& S']
        if closures=='LiAndSawamoto': 
            fig.suptitle(leg[0])
        elif closures=='vD_LiAndSawamoto': 
            fig.suptitle(leg[1])
        elif closures=='dispersion_Li_vD':
            fig.suptitle(leg[2])
        fig.savefig(sC+'_'+closures+'_compContrib.'+format,dpi=1200)
if __name__ == "__main__":
    print("File one executed when ran directly")
    figures_contrib('B5')
    plt.show()
else:
   print("File one executed when imported")

# %%
