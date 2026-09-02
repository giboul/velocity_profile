# -*- coding: utf-8 -*-
"""
Created on 4 Feb  2021

@author: Gauthier Rousseau

"""
import matplotlib
import os


import matplotlib.pyplot as plt
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import json


plt.style.use('seaborn-paper')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rc('font', size=9)
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.minor.width'] = 0.5
plt.rcParams['ytick.minor.size'] = 4.
plt.rcParams['ytick.minor.width'] = 0.5
plt.rcParams['xtick.major.size'] = 4.
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.size'] = 4.
plt.rcParams['ytick.major.width'] = 0.5

plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams["legend.fancybox"] =False
plt.rcParams['legend.labelspacing']= 0.5
plt.rcParams['legend.fontsize'] =  'small'
plt.rcParams['legend.borderpad'] = 0.3
plt.rcParams["legend.markerscale"] = 0.9
plt.rcParams["legend.framealpha"] = 1

rsd=r'$z^{\prime}/d_p$'
rsh=r'$z^{\prime}/h_f$'
# plt.rc('grid', linestyle="-", color='grey',linewidth=0.3)
# plt.rcParams['axes.grid'] = True
# plt.rc('text', usetex=True) may be activated if you prefer a Latex rendering font

def my_formatter(x, pos):
    """Format 1 as 1, 0 as 0, and all values whose absolute values is between
    0 and 1 without the leading "0." (e.g., 0.7 is formatted as .7 and -0.4 is
    formatted as -.4)."""
    val_str = '{:g}'.format(x)
    if np.abs(x) > 0 and np.abs(x) < 1:
        return val_str.replace("0", "", 1)
    else:
        return val_str


def loadData(jsonF='../data_Rousseau_JFM.json'):

    f = open(jsonF)
    data = json.load(f)
    cases = ['A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3', 'B4', 'B5',
    'L11','L12','L13','L14','L15']
    for ca in cases:
        for k in data[ca]['profiles'].keys():
            data[ca]['profiles'][k] = np.array(data[ca]['profiles'][k])
    f.close()
    return data 

def setPLotProfilFig2(namefig='profil', major_formatter=None):
    exty = 0.8
    extx = 0.7
    Ltot = 14
    mx = 1 / Ltot
    sx= 0.5 /Ltot
    shifty = 0.14
    xPor=1./Ltot
    axiswindow1 = [mx, shifty, 3.5/Ltot, exty]
    axiswindow2 = [mx+3.5/Ltot+sx, shifty, xPor, exty]
    axiswindow3 = [mx + 3.5 / Ltot + 2 * sx + xPor, shifty, 3.25/Ltot, exty]
    axiswindow4 = [mx+6.75/Ltot+3*sx+xPor, shifty, 3.25/Ltot, exty]

    # axiswindow3=[7./Ltot, shifty, 6./Ltot, exty]
#    axiswindow4=[11/Ltot, 0.1, 2.5/Ltot, exty]
    fig = plt.figure(namefig, figsize=(5.5, 3.5))
    ypos = 0.96
    fig.text(mx-sx/2, ypos, r'(a)')
    fig.text(mx+3.5/Ltot+sx-sx/2, ypos, r'(b)')
    fig.text(mx + 3.5 / Ltot + 2 * sx + xPor-sx/2, ypos, r'(c)')
    fig.text(mx+6.75/Ltot+3*sx+xPor-sx/2, ypos, r'(d)')

    ax1 = plt.Axes(fig, axiswindow1)
    ax2 = plt.Axes(fig, axiswindow2)
    ax3 = plt.Axes(fig, axiswindow3)
    ax4 = plt.Axes(fig, axiswindow4)
    fig.add_axes(ax1)
    fig.add_axes(ax2)
    fig.add_axes(ax3)
    fig.add_axes(ax4)

    ax2.set_yticklabels([])
    ax3.set_yticklabels([])
    ax1.set_xlabel(r' $ U_x ~[\mathrm{m~s^{-1}}]$')
    ax1.set_ylabel(rsd)
    ax2.set_xlabel(r'$\epsilon(z)$')
    ax3.set_xlabel(r'$\tau ~{[\mathrm{Pa}]}$')
    ax1.xaxis.set_major_formatter(major_formatter)
    ax1.yaxis.set_major_formatter(major_formatter)
    ax2.xaxis.set_major_formatter(major_formatter)
    ax3.xaxis.set_major_formatter(major_formatter)

    ax5 = plt.Axes(fig, [0.2, 0.2, 0.2, 0.2])
    fig.add_axes(ax5)
    ax5.set_position([2.08/Ltot, shifty+0.12, 2.24/Ltot, 0.25])

    ax4.set_xlabel(r'$\tau~{[\mathrm{Pa}]}$')
    ax4.set_yticklabels([])
    ax4.set_xlim([-0.05, 0.4])
    ax5.set_ylim(-1.5, 1.8)
    ax3.set_ylim(-1.9, 1.9)
    ax4.set_ylim(-1.9, 1.9)
    ax2.set_ylim(-1.9, 1.9)
    ax1.set_ylim(-1.9, 1.9)
    ax5.set_xticks([0.001, 0.01, 0.1])
    ax2.set_xticks([0, 0.5, 1])
    ax3.set_xticks([0, 1])
    ax5.set_ylabel(rsd)
    ax5.yaxis.set_label_coords(-0.17, 0.5)
    ax5.xaxis.set_label_coords(0.5, -0.2)
    ax1.yaxis.set_label_coords(-0.15, 0.5)
    ax5.set_xlabel(r'$ U_x~[\mathrm{m \cdotp s^{-1}}]$')


    ax2.set_xlim([0, 1.2])

    return fig, ax1, ax2, ax3, ax5, ax4


def dictModelGen(fluid, data=None):

    dictModel = {}
    dictModel['g'] = 9.81
    dictModel['mu'] = 0.67
    dictModel['kappa'] = 0.41
    dictModel['type'] = 'Ergun-vanDriest'
    T = 273.15+20
    dictModel['T'] = T
    dictModel['Pressure'] = 10**5
    dictModel['rhos'] = 2.43e3

    if fluid == 'BAE':
        # Attention : for T=293.15K=20C
        dictModel['xBA'] = 0.6
        dictModel['rhof'] = data['physico-chemical_properties']['isoMix']['rho']*1e3
        dictModel['dyn_visc'] = data['physico-chemical_properties']['isoMix']['visc']*1e-3
        dictModel['nu'] = dictModel['dyn_visc']/dictModel['rhof']
        dictModel['fluid'] = fluid
    if fluid == 'NaI':
        # Attention : for T=293.15K=20C
        dictModel['rhof']=1.77e3  
        dictModel['nu'] = 1.35e-6
        dictModel['dyn_visc']=dictModel['nu']*dictModel['rhof']
        dictModel['fluid'] = fluid
    # 

    dictModel['AEr'] = 180
    dictModel['BEr'] = 1.75

    return dictModel


def porGen(z80, Por_raw, indRC, sC, d):
    z100p = (z80[indRC]+z80[indRC-1])/2
    z100 = z80 - z100p
    if sC == 'A1':
        indysub = np.where(z100 < -d*0.8)
        avP = np.mean(Por_raw[indysub[0][:-4]])
        indT = np.where(Por_raw < avP)[0][0]
        Por = np.copy(Por_raw)
        Por[indT:] = avP
    elif sC == 'A2' or sC == 'A3' or sC == 'A4':
        indysub = np.where(z100 < -d*0.5)
        avP = np.mean(Por_raw[indysub[0][:-4]])
        indT = np.where(Por_raw < avP)[0][0]
        Por = np.copy(Por_raw)
        Por[indT:] = avP
    elif sC == 'B1' or sC == 'B2' or sC == 'B3' or sC == 'B4' or sC == 'B5':
        indysub = np.where(z80 < -d*0.4)
        avP = np.mean(Por_raw[indysub[0][:-20]])
        indT = np.where(Por_raw < avP)[0][0]
        Por = np.copy(Por_raw)
        Por[indT:] = avP
    return Por


def setPLotProfilComp(namefig='profil'):
    exty=0.8
    extx=0.7
    Ltot=14
    axiswindow1=[1./Ltot, 0.1, 3/Ltot, exty]
    axiswindow2=[4.5/Ltot, 0.1,2/Ltot, exty]
    axiswindow3=[7/Ltot, 0.1, 3/Ltot, exty]
    axiswindow4=[10.5/Ltot, 0.1, 3/Ltot, exty]
    fig = plt.figure(namefig,figsize=(5.5,3.5)) 
    ax1 = plt.Axes(fig, axiswindow1)
    ax2 = plt.Axes(fig, axiswindow2)
    ax3 = plt.Axes(fig, axiswindow3)
    ax4 = plt.Axes(fig, axiswindow4)
    fig.add_axes(ax1)
    fig.add_axes(ax2)
    fig.add_axes(ax3)
    fig.add_axes(ax4)
    ax2.set_yticklabels([])  
    ax3.set_yticklabels([])
    ax4.set_yticklabels([])
    ax1.set_xlabel(r'$\small  U_x ~[\mathrm{m \cdotp s^{-1}}]$',fontsize=10)
    ax1.set_ylabel(r'$z^{\prime}_{\epsilon=0.8}/d_p$',fontsize=10)
    ax2.set_xlabel(r'$\epsilon(z) $ ~[-]',fontsize=10)
    ax3.set_xlabel(r' $ \small \tau ~\small{[\mathrm{N~m^{-2}}]}$',fontsize=10)
    ax4.set_xlabel(r'$l_m \small{[\mathrm{m}]}$',fontsize=10)
    ax2.set_xlim([0.,1.])
    ax1.xaxis.set_major_formatter(major_formatter)
    ax1.yaxis.set_major_formatter(major_formatter)
    ax2.xaxis.set_major_formatter(major_formatter)
    ax3.xaxis.set_major_formatter(major_formatter)
    ax4.xaxis.set_major_formatter(major_formatter)
    ax5=ax4.twiny()
    ax5.xaxis.set_major_formatter(major_formatter)
    plt.show()
    
    return fig, ax1, ax2, ax3, ax4, ax5

def convtolist(data):
    cases = ['A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3', 'B4', 'B5','L11','L12','L13','L14','L15']
    for ca in cases:
        if ca in data.keys():
            for k in data[ca]['profiles'].keys():
                data[ca]['profiles'][k] = data[ca]['profiles'][k].tolist()
            if 'modellings' in data[ca].keys():
                keys = ['z', 'v', 'por', 'Afi', 'T', 'G', 'lM', 'fp', 'fv']
                dic = data[ca]['modellings']
                for kmod in dic.keys():
                    for k in keys:
                        if k in dic[kmod].keys():
                            dic[kmod][k] = dic[kmod][k].tolist()
                        data[ca]['modellings'][kmod]=dic[kmod]

    return data


def convNp(data):
    cases = ['A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3', 'B4', 'B5','L11', 'L12', 'L13', 'L14', 'L15']
    for ca in cases:
        if ca in data.keys():
            for k in data[ca]['profiles'].keys():
                data[ca]['profiles'][k] = np.array(data[ca]['profiles'][k])
            if 'modellings' in data[ca].keys():
                keys = ['z', 'v', 'por', 'Afi', 'T', 'G', 'lM', 'fp', 'fv']
                dic = data[ca]['modellings']
                for kmod in dic.keys():
                    for k in keys:
                        if k in dic[kmod].keys():
                            dic[kmod][k] = np.array(dic[kmod][k])
                        data[ca]['modellings'][kmod]=dic[kmod]

                                               

    return data