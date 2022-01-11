# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:04:35 2021

@author: CU
"""
import matplotlib.pyplot as plt
import numpy as np
from pandas import read_excel
from scipy.integrate import simpson

def eps2_viewer(bondingtype, txtfile, buildplot = True, latexfont=False, textsize = 20):    
        font = {'size'   : textsize}
        
        plt.rcParams.update({
        "text.usetex": latexfont,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        })
        
        data = np.loadtxt('./dataORDF/{0}/{1}'.format(bondingtype,txtfile), comments = '----')
        solid = txtfile[:-4]
        
        
        E = data[:,0]
        eps2_1 = data[:,1]
        eps2_2 = data[:,2]
        eps2_3 = data[:,3]
        
        
        eps2_winnderindex = np.argmax([max(eps2_1), max(eps2_2), max(eps2_3)])
        eps2_ofinterest = data[:,1+eps2_winnderindex]
        
        
        #Hier wird Metadata retrieved aus eps2maxDS.xlsx
    
        excel = 'appendix_metadata_DV.xlsx'
        df = read_excel(excel)
        #ET = df.loc[df['solid'] == '{}'.format(solid), 'ET'].iloc[0]
        #ES = df.loc[df['solid'] == '{}'.format(solid), 'ES'].iloc[0]
        bond = df.loc[df['solid_DF'] == '{}'.format(solid), 'bonding'].iloc[0]
        solid_alt = df.loc[df['solid_DF'] == '{}'.format(solid), 'solid'].iloc[0]
        structure = df.loc[df['solid_DF'] == '{}'.format(solid), 'structure'].iloc[0]
        structure_output = df.loc[df['solid_DF'] == '{}'.format(solid), 'structure_DV'].iloc[0]
        
        
        
        A = simpson(eps2_ofinterest, E)
        E_bar = simpson(eps2_ofinterest*E, E)/A
        eps2_bar = simpson(eps2_ofinterest*eps2_ofinterest, E)/(2*A)
        
        A_DELDIS = simpson(eps2_ofinterest[:1000], E[:1000])
        E_bar_DELDIS = simpson((eps2_ofinterest*E)[:1000], E[:1000])/A_DELDIS
        eps2_bar_DELDIS = simpson((eps2_ofinterest*eps2_ofinterest)[:1000], E[:1000])/(2*A_DELDIS)
        
        eps2max = max(eps2_ofinterest)
        Emax = E[np.argmax(eps2_ofinterest)]
        eps2_ofinterest_NORMALIZED = eps2_ofinterest/eps2max
        fwhm_indexlist = [i for i,element in enumerate(eps2_ofinterest_NORMALIZED) if element > eps2_bar/eps2max]
        FWHM = E[fwhm_indexlist[-1]] - E[fwhm_indexlist[0]]
        
        if structure == 'Fm3m':
            structurelabel = r'Fm$\bar{3}$m'
        elif structure == 'Pnma':
            structurelabel = 'Pnma'
        elif structure == 'R3m':
            structurelabel = 'R3m'
        elif structure == 'R-3m':
            structurelabel = r'R$\bar{3}$m'
        elif structure == 'P4nmm':
            structurelabel = 'P4/nmm'
        elif structure == 'P4mmm':
            structurelabel = 'P4/mmm'
        elif structure == 'C2m':
            structurelabel = 'C2/m'
        elif structure == 'Pm3m':
            structurelabel = r'Pm$\bar{3}$m'
        elif structure == 'Fd3m':
            structurelabel = r'Fd$\bar{3}$m'
        elif structure == 'F43m':
            structurelabel = r'F$\bar{4}3$m'
        elif structure == 'P63mc':
            structurelabel = 'P6$_3$mc'
        elif structure == 'Pmc21':
            structurelabel = 'Pmc2$_1$'
        else:
            structurelabel = ''
        
        boxdata = (r'$\varepsilon_{2}^\mathrm{max}$ = %s ' % round(eps2max,3),
                   r'$E_\mathrm{max}$ = %seV' % round(Emax,3),
                   '',
                   r'$\varepsilon_{2}^\mathrm{cen}$ = %s' % round(eps2_bar,3),
                   r'$E_\mathrm{cen}$ = %seV' % round(E_bar,3),
                   '',
                   r'$\Delta E$ = %seV' % round(FWHM,3))
        textstr = ' \n '.join(boxdata)
        
        fig = plt.figure(figsize = (12.8,7.2))
        if buildplot == True:
            if bond == 1:
                color = 'black'
                name = 'Ionic'
            elif bond == 0:
                color = 'red'
                name = 'Covalent'
            elif bond == 3:
                color = 'green'
                name = 'Metavalent'
                
        
            ax = fig.add_subplot(111)
            ax.plot(E, eps2_ofinterest, color = color)
            ax.vlines(Emax, ymin = 0, ymax = max(eps2_ofinterest), color = '#00549f')
            ax.hlines(eps2_bar, xmin = E[fwhm_indexlist[0]], xmax = E[fwhm_indexlist[-1]], color = '#00549f')
            ax.scatter(E_bar_DELDIS, eps2_bar_DELDIS, color = 'orange')
            ax.scatter(E_bar, eps2_bar, color = '#00549f')
            #ax.annotate(r'($E_\mathrm{cen}, \varepsilon_{2}^\mathrm{cen}$)', (E_bar, eps2_bar), fontsize = 'xx-large')
            
            ax.set_xlim(0,35)
            ax.set_ylim(0)
            ax.set_title(solid_alt+' ('+structurelabel+') '+' - Imaginary part of dielectric function', **font)
            ax.set_xlabel('Energy (eV)', **font)
            ax.set_ylabel(r'$\varepsilon_2$', **font)
            #ax.legend() # <---handles = legend_elements (stand vorher mal drin)
            ax.tick_params(axis='both', labelsize=16)
            
            if latexfont == True:
                x_axtext = 0.825
                y_axtext = 0.666
            else:
                x_axtext = 0.8
                y_axtext = 0.6
            
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            ax.text(x_axtext, y_axtext, s=textstr, transform = ax.transAxes, fontsize = 'xx-large', bbox=props)
            
        #plt.xticks(**font)
        plt.tight_layout()
        plt.show(block=False)
        
        EpsilonMAP = eps2max/(Emax+eps2max)
        
        print('____________________________')
        print(solid_alt+' '+structure_output, '\n')
        
        print('eps2_max', eps2max)
        print('E_max: ', Emax, 'eV', '\n')
        
        print('eps2_cen ', eps2_bar)
        print('E_cen ', E_bar, 'eV', '\n')
        
        print('DeltaE: ', FWHM, 'eV')
        print('____________________________')
        
        return eps2max, Emax, FWHM, EpsilonMAP, eps2_ofinterest, fig


from tkinter import *
import os
df = read_excel('appendix_metadata_DV.xlsx')
y = []
root = Tk()
root.title("DieViz")
try:
    root.iconbitmap('iconDieViz.ico')
except:
    None



def ionic():
    path = './dataORDF/ionic'
    solids = os.listdir(path)
    for i in y:
        i.destroy()
    l=1
    n=0
    for i in solids:
        solid_alt = df.loc[df['solid_DF'] == '{}'.format(i[:-4]), 'solid_DV'].iloc[0]
        structure = df.loc[df['solid_DF'] == '{}'.format(i[:-4]), 'structure_DV'].iloc[0]
        
        solidBut = Button(root, text = solid_alt+'\n'+structure, fg = 'white', bg = 'black', command =lambda i=i: eps2_viewer(bondingtype='ionic', txtfile=i))
        solidBut.grid(row=l, column=n)
        n+=1
        if n % 10 == 0:
            l+=1
            n=0
            
        y.append(solidBut)



def cova():
    path = './dataORDF/covalent'
    solids = os.listdir(path)
    for i in y:
        i.destroy()
    l=1
    n=0
    for i in solids:
        solid_alt = df.loc[df['solid_DF'] == '{}'.format(i[:-4]), 'solid_DV'].iloc[0]
        structure = df.loc[df['solid_DF'] == '{}'.format(i[:-4]), 'structure_DV'].iloc[0]
        
        solidBut = Button(root, text = solid_alt+'\n'+structure, fg = 'white', bg = 'red', command =lambda i=i: eps2_viewer(bondingtype='covalent', txtfile=i))
        solidBut.grid(row=l, column=n)
        n+=1
        if n % 10 == 0:
            l+=1
            n=0
        
        y.append(solidBut)




def mv():
    path = './dataORDF/metavalent'
    solids = os.listdir(path)
    for i in y:
        i.destroy()
    l=1
    n=0
    for i in solids:
        solid_alt = df.loc[df['solid_DF'] == '{}'.format(i[:-4]), 'solid_DV'].iloc[0]
        structure = df.loc[df['solid_DF'] == '{}'.format(i[:-4]), 'structure_DV'].iloc[0]
        
        solidBut = Button(root, text = solid_alt+'\n'+structure, fg = 'white', bg = 'green', command =lambda i=i: eps2_viewer(bondingtype='metavalent', txtfile=i))
        solidBut.grid(row=l, column=n)
        n+=1
        if n % 10 == 0:
            l+=1
            n=0
        
        y.append(solidBut)





ionButton = Button(root, text='Ionic', fg = 'white', bg = 'black', command = ionic)
covButton = Button(root, text='Covalent', fg = 'white', bg = 'red', command = cova)
mvButton = Button(root, text='Metavalent', fg = 'white', bg = 'green', command = mv)



ionButton.grid(row=0, column=0, padx=10)
covButton.grid(row=0, column=1, padx=10)
mvButton.grid(row=0, column=2, padx=10)

root.mainloop()