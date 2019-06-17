import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import astropy 
from astropy import cosmology
from astropy.cosmology import WMAP9 as cosmo
import math as mt
import os, sys
import matplotlib.lines as mlines

#import pyfits 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import astropy 
from astropy import cosmology 
import math as mt
from scipy.integrate import quad
from scipy.stats import chi2_contingency
from pylab import *
from scipy.optimize import curve_fit
import scipy as sp
import scipy.special
import scipy.stats as stats
import seaborn as sns 
import pandas as pd
from astropy.modeling.models import Sersic1D
from numpy import * 
import bces.bces
import nmmn.stats



def uvj(cf_Ha,gf_Ha,crf_VJ,crf_UV,crf_m,cz):
    plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    f = open('/Users/jennifercooper/Projects/thesis/gprior/Ha_catalogs/all_core_sn.txt', 'r')
    lines = f.readlines()
    f.close()
    
    #create arrays 
    cf_Ha    = [] 
    crf_VJ   = []
    crf_UV   = []
    crf_m    = []
    cz       = []
  
    
    #pull array column 
    for line in lines: 
        a = line.split()
        cf_Ha.append(float(a[4]))
        crf_VJ.append(float(a[7]))
        crf_UV.append(float(a[6]))
        crf_m.append(float(a[8]))
        cz.append(float(a[3]))
        
    cf_Ha = np.array(cf_Ha)
    cf_Ha_log = np.log10(cf_Ha)
    crf_VJ = np.array(crf_VJ)
    crf_UV = np.array(crf_UV)
    crf_m  = np.array(crf_m)
    cz     = np.array(cz)
    cD_l_mpc = np.array(cosmo.luminosity_distance(cz))
    cD_l_cm = cD_l_mpc*3.08567758128E+24
    clogLHa = np.log10(4*np.pi*cf_Ha*cD_l_cm**2)
    cSFR = clogLHa - 41.27
    #print(np.mean(cSFR))
    #print(np.std(cSFR))


    f = open('/Users/jennifercooper/Projects/thesis/gprior/Ha_catalogs/all_group_sn.txt', 'r')
    lines = f.readlines()
    f.close()
    
    #create arrays 
    gf_Ha    = [] 
    grf_VJ   = []
    grf_UV   = []
    grf_m    = []
    gz       = []
  
    
    #pull array column 
    for line in lines: 
        a = line.split()
        gf_Ha.append(float(a[4]))
        grf_VJ.append(float(a[7]))
        grf_UV.append(float(a[6]))
        grf_m.append(float(a[8]))
        gz.append(float(a[3]))
        
    gf_Ha = np.array(gf_Ha)
    gf_Ha_log = np.log10(gf_Ha)
    grf_VJ = np.array(grf_VJ)
    grf_UV = np.array(grf_UV)
    grf_m  = np.array(grf_m)
    gz     = np.array(gz)
    gD_l_mpc = np.array(cosmo.luminosity_distance(gz))
    gD_l_cm = gD_l_mpc*3.08567758128E+24
    glogLHa = np.log10(4*np.pi*gf_Ha*gD_l_cm**2)
    gSFR = glogLHa - 41.27
    #print(np.mean(gSFR))
    #print(np.std(gSFR))

    f = open('/Users/jennifercooper/Projects/thesis/gprior/Ha_catalogs/all_field_sn.txt', 'r')
    lines = f.readlines()
    f.close()
    
    #create arrays 
    f_Ha    = [] 
    rf_VJ   = []
    rf_UV   = []
    rf_m    = []
    z       = []
  
    
    #pull array column 
    for line in lines:
        a = line.split()
        f_Ha.append(float(a[4]))
        rf_VJ.append(float(a[7]))
        rf_UV.append(float(a[6]))
        rf_m.append(float(a[8]))
        z.append(float(a[3]))
        
    f_Ha = np.array(f_Ha)
    f_Ha_log = np.log10(f_Ha)
    rf_VJ = np.array(rf_VJ)
    rf_UV = np.array(rf_UV)
    rf_m  = np.array(rf_m)
    z     = np.array(z)
    D_l_mpc = np.array(cosmo.luminosity_distance(z))
    D_l_cm = D_l_mpc*3.08567758128E+24
    logLHa = np.log10(4*np.pi*f_Ha*D_l_cm**2)
    SFR = logLHa - 41.27
    #print(np.mean(SFR))
    #print(np.std(SFR))
    

 
    
    x1 = [-5,0.85]
    x2 = [0.85,1.6]
    x3 = [1.6,1.6]
    y1 = [1.3,1.3]
    y2 = [1.3,2]
    y3 = [2,5]

    f = open('/Users/jennifercooper/Projects/thesis/WFI_catalogs/1301.txt', 'r')
    lines = f.readlines()[1:]
    f.close()
    
    #create arrays 
    LDP_UV    = [] 
    LDP_VJ   = []
    
    
    #pull array column 
    for line in lines: 
        a = line.split()
        LDP_UV.append(float(a[5]))
        LDP_VJ.append(float(a[6]))
        
        
    LDP_UV = np.array(LDP_UV)
    LDP_VJ = np.array(LDP_VJ)

    
    plt.subplot(3,2,1)
    #plt.hist2d(LDP_VJ, LDP_UV, bins=(250, 150), cmap=plt.cm.Greys)
    plt.scatter(LDP_VJ, LDP_UV, color = 'grey', alpha = 0.3)
    plt.scatter(crf_VJ,crf_UV,c=cSFR,s=100)
    cbar = plt.colorbar()
    #cbar.set_clim(-0.5,1.2)
    cbar.set_label("log(SFR)")
    #plt.colorbar()

    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster core UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,2.5)
    
    plt.subplot(3,2,2)
    #plt.hist2d(LDP_VJ, LDP_UV, bins=(150, 150), cmap=plt.cm.Greys)
    plt.scatter(LDP_VJ, LDP_UV, color = 'grey', alpha = 0.3)
    plt.scatter(grf_VJ,grf_UV,c=gSFR,s=100)
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    #cbar.set_clim(-0.5,1.2)
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster group UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,2.5)

    plt.subplot(3,2,3)
    #plt.hist2d(LDP_VJ, LDP_UV, bins=(150, 150), cmap=plt.cm.Greys)
    plt.scatter(LDP_VJ, LDP_UV, color = 'grey', alpha = 0.3)
    plt.scatter(rf_VJ,rf_UV,c=SFR,s=100)
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    #cbar.set_clim(-0.5,1.2)
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Field UVJ')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,2.5)
   
    
    plt.subplot(3,2,4)
    plt.scatter(crf_m,cSFR,color='blue',alpha=0.7,marker = 'o',s=100)
    plt.scatter(grf_m,gSFR,color='blue',alpha=0.5,marker = 'x',s=100)
    plt.scatter(rf_m,SFR,c='red',alpha=0.5,s=100) 
    plt.xlim(8,11)
    #plt.ylim(-1,3)
    plt.xlabel('logM')
    plt.ylabel('logSFR')
    #plt.title('All pointings matched with LPD Halpha>0')
    red_patch = mpatches.Patch(color='red', label='Field', alpha = 0.5)
    blue_patch = mpatches.Patch(color='blue', label='Core', alpha = 0.5)
    blue_star = mlines.Line2D([], [], color='blue', marker='x', linestyle='None',
                          markersize=10, label='Groups')
    plt.legend(handles=[red_patch,blue_patch,blue_star],loc=1)
    
    
    #zcl1       = cz[np.where((cz>0.4628)&(cz<0.5028))]
    
   # zf1       = z[np.where(z<0.4628)]   
   # zf2       = z[np.where((z>0.5028))]

   # x1 = zc
   # x2 = zf1
   # x3 = zf2
   # 
    
   # kwargs1 = dict(histtype='stepfilled', alpha=0.3, normed=False, bins=2)
   # kwargs2 = dict(histtype='stepfilled', alpha=0.3, normed=False, bins=5)
   ## kwargs3 = dict(histtype='stepfilled', alpha=0.3, normed=False, bins=5)



    
   # plt.hist(x1,color='blue', **kwargs1)
   # plt.hist(x2,color='red', **kwargs2)
   # plt.hist(x3,color='red', **kwargs3);
   # plt.xlabel('z')


    plt.show()

def boot_SFR_M(crf_m,cz):
    plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    f = open('/Users/jennifercooper/Projects/thesis/gprior/Ha_catalogs/all_core.txt', 'r')
    lines = f.readlines()
    f.close()
    
    #create arrays 
    cf_Ha    = [] 
    crf_m    = []
    cz       = []
  
    
    #pull array column 
    for line in lines: 
        a = line.split()
        cf_Ha.append(float(a[4]))
        crf_m.append(float(a[8]))
        cz.append(float(a[3]))
        
    cf_Ha = np.array(cf_Ha)
    cf_Ha_log = np.log10(cf_Ha)
    crf_m  = np.array(crf_m)
    cz     = np.array(cz)
    cD_l_mpc = np.array(cosmo.luminosity_distance(cz))
    cD_l_cm = cD_l_mpc*3.08567758128E+24
    clogLHa = np.log10(4*np.pi*cf_Ha*cD_l_cm**2)
    cSFR = clogLHa - 41.27


    f = open('/Users/jennifercooper/Projects/thesis/gprior/Ha_catalogs/all_group_sn.txt', 'r')
    lines = f.readlines()
    f.close()
    
    #create arrays 
    gf_Ha    = [] 
    grf_m    = []
    gz       = []
  
    
    #pull array column 
    for line in lines: 
        a = line.split()
        gf_Ha.append(float(a[4]))
        grf_m.append(float(a[8]))
        gz.append(float(a[3]))
        
    gf_Ha = np.array(gf_Ha)
    gf_Ha_log = np.log10(gf_Ha)
    grf_m  = np.array(grf_m)
    gz     = np.array(gz)
    gD_l_mpc = np.array(cosmo.luminosity_distance(gz))
    gD_l_cm = gD_l_mpc*3.08567758128E+24
    glogLHa = np.log10(4*np.pi*gf_Ha*gD_l_cm**2)
    gSFR = glogLHa - 41.27


    f = open('/Users/jennifercooper/Projects/thesis/gprior/Ha_catalogs/all_field_sn.txt', 'r')
    lines = f.readlines()
    f.close()
    
    #create arrays 
    f_Ha    = [] 
    rf_m    = []
    z       = []
  
    
    #pull array column 
    for line in lines: 
        a = line.split()
        f_Ha.append(float(a[4]))
        rf_m.append(float(a[8]))
        z.append(float(a[3]))
        
    f_Ha = np.array(f_Ha)
    f_Ha_log = np.log10(f_Ha)
    rf_m  = np.array(rf_m)
    z     = np.array(z)
    D_l_mpc = np.array(cosmo.luminosity_distance(z))
    D_l_cm = D_l_mpc*3.08567758128E+24
    logLHa = np.log10(4*np.pi*f_Ha*D_l_cm**2)
    SFR = logLHa - 41.27

    x = 0
    y = 0
    yer = 0
    lcb1 = 0
    lcb2 = 0
    lcb3 = 0
    ucb1 = 0
    ucb2 = 0
    ucb3 = 0
    xcb1 = 0
    import warnings
    warnings.filterwarnings("ignore")

    ax=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    ax = sns.set(style="white", color_codes=True)
    ax = subplot(2,2,1)
    ax = xlim(8,11)
    ax = ylim(-1,2)
    x, y = pd.Series(crf_m, name="Mass Core"), pd.Series(cSFR, name="SFR Core")
    sort = np.argsort(x)
    x = x[sort]
    x = np.array(x)
    y = y[sort] 
    y = np.array(y)
    yer = zeros(len(x))
    xer=zeros(len(x))
    cov=zeros(len(x))   # no correlation between error measurements
    i = 0
    nboot=10000   # number of bootstrapping trials
    def func(x): return x[1]*x[0]+x[2]
    a,b,aerr,berr,covab=bces.bces.bcesp(x,xer,y,yer,cov,nboot)
    ybces=a[3]*x+b[3]  # the integer corresponds to the desired BCES method for plotting (3-ort, 0-y|x, 1-x|y, *don't use bissector*)
    # array with best-fit parameters
    fitm=np.array([ a[i],b[i] ])
    # covariance matrix of parameter uncertainties
    covm=np.array([ (aerr[i]**2,covab[i]), (covab[i],berr[i]**2) ])
    # Gets lower and upper bounds on the confidence band
    print(lcb1,lcb2)
    lcb1,ucb1,x=nmmn.stats.confband(x, y, a[i], b[i], 0.68, x)
    lcb2,ucb2,x2=nmmn.stats.confband(x, y, a[i], b[i], 0.95, x)
    lcb3,ucb3,x3=nmmn.stats.confband(x, y, a[i], b[i], 0.997, x)
    #errorbar(x,y,yerr=yer,fmt='o')
    ax = plot(x,a[i]*x+b[i],'-k')
    ax = fill_between(x, lcb1, ucb1, alpha=0.6, facecolor='purple')
    ax = fill_between(x, lcb2, ucb2, alpha=0.3, facecolor='blue')
    ax = fill_between(x, lcb3, ucb3, alpha=0.4, facecolor='grey')
    #ax = bar(x,y,yerr=[f_dk_p[np.where(f_flag>0.9)], f_dk_n[np.where(f_flag>0.9)]], facecolor='none')
    ax = legend(loc='best')
    ax = xlabel('LogM')
    ax = ylabel('Log(SFR)')



    

