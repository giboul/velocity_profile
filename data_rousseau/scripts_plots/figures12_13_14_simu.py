
# -*- coding: utf-8 -*-
"""
Created on 4 Feb  2021

@author: Gauthier Rousseau

"""
# %%
# %matplotlib qt5
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

data=tl.loadData()
data=tl.convNp(data)
major_formatter = FuncFormatter(tl.my_formatter)

# to display the possible keys type :
# data['A1'].keys()
# # to call the porosity profile for instance type
# data['A1']['profiles']['Porosity']


plt.close('all')


def figures12_13_14(selected_Case='A1',format='pdf'):
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
   
    i_mixinglength = 5
    i_trubstress = 4
    i_dispstress = 3
    i_forces = 2
    i_velocity = 1
    i_comp = 0

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
        exty=0.16
        extx=0.8
        my, sy = 0.06, 0.025
        mx = 0.12
        Nplots = 6
        exty=0.18*4/Nplots
        axiswindows = []
        AX = []
        mxlet=mx/4
        fig = plt.figure(namefig, figsize=(4.88, 4.5*Nplots/4))
        letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
        for k in range(Nplots):
            aw=[mx, my + (k ) * exty + (k ) * sy, extx, exty]
            axiswindows.append(aw)
            a=plt.Axes(fig, aw)
            AX.append(a)
            fig.add_axes(a)
            a.xaxis.set_major_formatter(major_formatter)
            a.yaxis.set_major_formatter(major_formatter)
            fig.text(mxlet, my+(k+1)*exty+(k)*sy+0.035,   letters[Nplots-k-1])
            # 

        for a in AX[1:Nplots]:
            a.set_xticklabels([]) 
            Pos=a.get_position().bounds 
            a.set_position([Pos[0],Pos[1]+0.035,Pos[2],Pos[3]])
        Pos=AX[0].get_position().bounds 
        AX[0].set_position([Pos[0],Pos[1],Pos[2],Pos[3]+0.035])
        AX[i_trubstress].set_ylabel(r'$-\langle \overline{u^\prime_x u^\prime_z} \rangle/  u_p^2$  ')
        AX[i_mixinglength].set_ylabel(r'$\ell_t /d_p$')
        AX[i_forces].set_ylabel(r'$\frac{f}{\epsilon \rho g i}$')
        AX[i_dispstress].set_ylabel(r'$-\langle \tilde{u}_x \tilde{u}_z \rangle/  u_p^2$')
        AX[i_velocity].set_ylabel(r'$ U_x / u_p$')
        AX[i_comp].set_ylabel(r'$\frac{f_{contrib}}{\epsilon \rho g i}$')
        AX[i_bottom].set_xlabel(tl.rsd)

        
        
        # fig.text(mxlet, my+2*exty+1*sy, '(c)')
        # fig.text(mxlet, my+1*exty, '(d)')
                    
        return fig, AX



    fig, AX = setPLotProfilComp(sC)

    porC=0.8

    indzinf=0
    indzsup = 0
    mod='modellings'
    for k in data[sC][mod].keys():
        z = data[sC][mod][k]['z']

    for a in AX:
        a.set_xlim([-0.8,z80[indup]/d+0.2])
    cm1 = cm.make_cmap_customized(Palette='mountain')
    cm2 = cm.make_cmap_customized(Palette='green')

    ### plot expe
    por_exp=data[sC]['profiles']['porosity']
    AX[i_velocity].scatter(z80[indup::2]/d,data[sC]['profiles']['ux'][indup::2]/uetp,s=5,color=cm2(0.7))
    # AX[1].plot(z80[indup:]/d,-data[sC]['profiles']['tautxz'][indup:],'-.',linewidth=1,c=cm2(0.7))
    l0=AX[i_trubstress].scatter(z80[indup::2]/d,-data[sC]['profiles']['tautxz'][indup::2]/(por_exp[indup::2]*rhof*uetp**2),s=5,color=cm2(0.7))

    # AX[i_dispstress].plot(z80[indup:] / d, -data[sC]['profiles']['tausxz'][indup:], '-.', linewidth=1, c=cm2(0.7))
    if sC=='L12':
        scatter =AX[i_dispstress].scatter(z80[indup::2] / d, -data[sC]['profiles']['tausxz'][indup::2]/(rhof*uetp**2),s=5,color=cm2(0.7))
    else:
        scatter =AX[i_dispstress].scatter(z80[indup::2] / d, -data[sC]['profiles']['tausxz'][indup::2]/(por_exp[indup::2]*rhof*uetp**2),s=5,color=cm2(0.7))

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
    scatter = AX[i_mixinglength].scatter(zlm[indup::2] / d, lm[indup::2]/ d, s=5, color=cm2(0.7))


    #### plot Model ###
    h1, l1 = [scatter], ["Exp. - "+sC]
    ic = 0
    lw = [.8, 0.9, 0.7]
    ls = ['--', 'dotted', '-']
    leg = [r'L \& S', r'damp. + L \& S', r' disp. + damp. + L \& S']
    mod='modellings'
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    for k in ['LiAndSawamoto', 'vD_LiAndSawamoto','dispersion_Li_vD']:
        z = data[sC][mod][k]['z']
        lz=len(z)
        Indz=range(indzinf+1,lz-2-indzsup)
        por = data[sC][mod][k]['por']
        porM=por[1:-1]
        dz=np.diff(z)[0]
        hsup = z[np.where(por > porC)[0][0]]
        Zr = (z[Indz] - hsup) / d
        v=data[sC][mod][k]['v']
        dv=np.diff(v)/dz
        dvM=np.mean(np.array([dv[1:],dv[:-1]]),axis=0)
        porM=np.mean(np.array([porM[1:],porM[:-1]]),axis=0)
        AX[i_velocity].plot(Zr, v[Indz]/uetp, c=cm1((ic+1) / 3),lw=lw[ic],ls=ls[ic])
        T=data[sC][mod][k]['T']
        taut =  T[1:-1] * (dvM ** 2)
        Indz2=range(indzinf,lz-3-indzsup)
        AX[i_trubstress].plot(Zr, taut[Indz2]/(porM*uetp**2), c=cm1((ic+1) / 3),lw=lw[ic],ls=ls[ic])
        lM=data[sC][mod][k]['lM']
        G=data[sC][mod][k]['G']
        fp=data[sC][mod][k]['fp']
        fv=data[sC][mod][k]['fv']
        [line2D]=AX[i_forces].plot(Zr,(fp[Indz]+fv[Indz])/rhof/G[Indz],color=cm1((ic+1) / 3),label=r'$ $',lw=lw[ic],ls=ls[ic])
        [line2D]=AX[i_mixinglength].plot(Zr,lM[Indz]/d,color=cm1((ic+1) / 3),label=r'$ l_{m,mod}$',lw=lw[ic],ls=ls[ic])
        Afi=data[sC][mod][k]['Afi']
        taufi =  Afi[1:-1] * v[1:-1] * dvM

        [line2D]=AX[i_dispstress].plot(Zr, taufi[Indz2]/(porM*uetp**2), c=cm1((ic + 1) / 3), lw=lw[ic], ls=ls[ic])  
        h, l = [line2D],[leg[ic]]

        h1, l1 = h1+h, l1+l  
        if k=='dispersion_Li_vD':
        # plot forces contrib comparison
            dtaut=np.diff(taut)/dz
            dtaut=(dtaut[:-1]+dtaut[1:])/2
            dtaut=smooth(dtaut,10)

            vM=np.mean(np.array([v[1:],v[:-1]]),axis=0)
            AfiM=np.mean(np.array([Afi[1:],Afi[:-1]]),axis=0) 
            # taufi =  Afi[1:-1] * v[1:-1] * dvM
            dAfi=np.diff(Afi)/dz
            elem1=(dAfi * vM + dv * AfiM)*dv
            dtaufi=(elem1[:-1]+Afi[1:-1]*v[1:-1]*np.diff(dv)/dz)
            Indz3=range(indzinf,lz-6-indzsup)
            [line2D]=AX[i_comp].plot(Zr[1:-2],(dtaut[Indz3])/G[Indz3],color=cm1((ic+1) / 3),label=r'$\frac{\mathrm{d} \tau_t}{\mathrm{d} z} $',lw=1,ls=':')
            [line2D]=AX[i_comp].plot(Zr[1:-2],(dtaufi[Indz3])/G[Indz3],color=cm1((ic+1) / 3),label=r'$\frac{\mathrm{d} \tau_d}{\mathrm{d} z} $',lw=lw[ic],ls='--')

            tau_v=nu*np.diff(por*v)/dz
            tau_v=(tau_v[:-1]+tau_v[1:])/2
            dtau_v=np.diff(tau_v)/dz
            dtau_vM=(dtau_v[:-1]+dtau_v[1:])/2
            dtau_vM=smooth(dtau_vM,6)
            [line2D]=AX[i_comp].plot(Zr[1:-2],(dtau_vM[Indz3])/G[Indz3],color=cm1((ic+1) / 3),lw=lw[ic],ls=':',label=r'$\frac{\mathrm{d} \tau_v}{\mathrm{d} z} $',)
        
            dpor=np.diff(por)/dz
            dporM=np.mean(np.array([dpor[1:],dpor[:-1]]),axis=0)    
            SecondBrickCorr=nu*dvM*dporM

            [line2D]=AX[i_comp].plot(Zr,SecondBrickCorr[Indz2]/G[Indz2],color='b',label=r'$-\varrho \nu \frac{\mathrm{d} U_x}{\mathrm{d} z} \frac{\mathrm{d} \epsilon}{\mathrm{d} z} $',lw=lw[ic],ls=':')

            [line2D]=AX[i_comp].plot(Zr,(fp[Indz]+fv[Indz])/rhof/G[Indz],color=cm1((ic+1) / 3),label=r'$\epsilon \bar{f}$',lw=lw[ic],ls=ls[ic])

            sumF=(fp[2:-2]+fv[2:-2])/rhof+dtaut+dtaufi[1:-1]+dtau_vM-SecondBrickCorr[1:-1]
            # [line2D]=AX[i_comp].plot(Zr[1:-2],sumF[Indz3]/G[Indz3],color='g',label=r'$ f$',lw=lw[ic],ls=ls[ic])
            AX[i_comp].legend(fontsize=7)

            # fig.savefig('test2.pdf')
            AX[i_comp].plot([-10,10],[0,0],c=(0.1,0.1,0.1),lw=0.4,zorder=-1)

            # fig.text(0.5,0.205,r' Disp. + vD + L \& S - force contributions',fontsize=8,ha='center')
        ic += 1

    # AX[0].plot(z80[indup::2] / d, data[sC]['profiles']['ux'][indup::2], '-.', linewidth=1, c=cm2(0.7))

    # l1b=AX[1].plot(z/d,z*0.4/d,'-',label=r'$0.4\, z$',lw=0.6,ls='-',c='k')

        
    fig.legend(h1 , l1, fontsize=6,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=4,bbox_to_anchor=(0.1, 0.95, 0.4, 0.2))
    # AX[i_forces].set_ylim([0, 5])
    AX[i_mixinglength].set_ylim([0, 0.004/d])
    if cases=='voermans':
        AX[i_mixinglength].set_ylim([0, 0.01/d])
    for a in AX:
        a.yaxis.set_label_coords(-0.07, 0.5)

    # plt.pause(0.1)


    # fig.savefig('test2.pdf')

    fig.savefig(sC+'.'+format,dpi=1200)


if __name__ == "__main__":
    print("File one executed when ran directly")
    figures12_13_14('L12')
    plt.show()
else:
   print("File one executed when imported")

# %%
