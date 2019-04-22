import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import astropy 
from astropy import cosmology
from astropy.cosmology import WMAP9 as cosmo
import math as mt
import os, sys
import matplotlib.lines as mlines




def c1301(cf_Ha,gf_Ha,crf_VJ,crf_UV,crf_m,cz):
    plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/c-1301-gprior.txt', 'r')
    lines = f.readlines()[1:]
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
        crf_VJ.append(float(a[6]))
        crf_UV.append(float(a[7]))
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


    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/g-1301-gprior.txt', 'r')
    lines = f.readlines()[1:]
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
        grf_VJ.append(float(a[6]))
        grf_UV.append(float(a[7]))
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
    
    x1 = [-5,0.85]
    x2 = [0.85,1.6]
    x3 = [1.6,1.6]
    y1 = [1.3,1.3]
    y2 = [1.3,2]
    y3 = [2,5]
    
    plt.subplot(3,2,1)
    plt.scatter(crf_VJ[np.where(np.logical_and(cz>0.4628,cz<0.5028))],crf_UV[np.where(np.logical_and(cz>0.4628,cz<0.5028))],c=cSFR[np.where(np.logical_and(cz>0.4628,cz<0.5028))])
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster core UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)
    
    plt.subplot(3,2,2)
    #plt.scatter(rf_VJ[np.where(np.logical_and(z<0.4628,z>0.5028))],rf_UV[np.where(np.logical_and(z<0.4628,z>0.5028))],c=SFR[np.where(np.logical_and(z<0.4628,z>0.5028))])
    plt.scatter(crf_VJ[np.where(cz<0.4628)],crf_UV[np.where(cz<0.4628)],c=cSFR[np.where(cz<0.4628)])
    plt.scatter(crf_VJ[np.where(cz>0.5028)],crf_UV[np.where(cz>0.5028)],c=cSFR[np.where(cz>0.5028)])
    plt.scatter(grf_VJ[np.where(gz<0.4628)],grf_UV[np.where(gz<0.4628)],c=gSFR[np.where(gz<0.4628)])
    plt.scatter(grf_VJ[np.where(gz>0.5028)],grf_UV[np.where(gz>0.5028)],c=gSFR[np.where(gz>0.5028)])
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    cbar.set_clim(-0.2,0.8)
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Field UVJ')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)

    plt.subplot(3,2,3)
    plt.scatter(grf_VJ[np.where(np.logical_and(gz>0.4628,gz<0.5028))],grf_UV[np.where(np.logical_and(gz>0.4628,gz<0.5028))],c=gSFR[np.where(np.logical_and(gz>0.4628,gz<0.5028))])
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster groups UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)
   
    
    plt.subplot(3,2,4)
    plt.scatter(crf_m[np.where((cz>0.4628)&(cz<0.5028))],cSFR[np.where((cz>0.4628)&(cz<0.5028))],color='blue',alpha=0.7,marker = 'o')
    plt.scatter(grf_m[np.where((gz>0.4628)&(gz<0.5028))],gSFR[np.where((gz>0.4628)&(gz<0.5028))],color='blue',alpha=0.5,marker = 'x')
    plt.scatter(crf_m[np.where(cz<0.4628)],cSFR[np.where(cz<0.4628)],c='red',alpha=0.5)
    plt.scatter(crf_m[np.where(cz>0.5028)],cSFR[np.where(cz>0.5028)],c='red',alpha=0.5)
    plt.scatter(grf_m[np.where(gz<0.4628)],gSFR[np.where(gz<0.4628)],c='red',alpha=0.5)
    plt.scatter(grf_m[np.where(gz>0.5028)],gSFR[np.where(gz>0.5028)],c='red',alpha=0.5) 
    plt.xlim(6.8,11)
    #plt.ylim(-1,3)
    plt.xlabel('logM')
    plt.ylabel('logSFR')
    #plt.title('All pointings matched with LPD Halpha>0')
    red_patch = mpatches.Patch(color='red', label='Field', alpha = 0.5)
    blue_patch = mpatches.Patch(color='blue', label='Core', alpha = 0.5)
    blue_star = mlines.Line2D([], [], color='blue', marker='x', linestyle='None',
                          markersize=10, label='Groups')
    plt.legend(handles=[red_patch,blue_patch,blue_star],loc=2)
    
    
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
    


def c1138(cf_Ha,gf_Ha,crf_VJ,crf_UV,crf_m,cz):
    plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/c-1138-gprior.txt', 'r')
    lines = f.readlines()[1:]
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
        crf_VJ.append(float(a[6]))
        crf_UV.append(float(a[7]))
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


    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/g-1138-gprior.txt', 'r')
    lines = f.readlines()[1:]
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
        grf_VJ.append(float(a[6]))
        grf_UV.append(float(a[7]))
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

    x1 = [-5,0.85]
    x2 = [0.85,1.6]
    x3 = [1.6,1.6]
    y1 = [1.3,1.3]
    y2 = [1.3,2]
    y3 = [2,5]

    plt.subplot(3,2,1)
    plt.scatter(crf_VJ[np.where(np.logical_and(cz<0.4996,cz>0.4596))],crf_UV[np.where(np.logical_and(cz<0.4996,cz>0.4596))],c=cSFR[np.where(np.logical_and(cz<0.4996,cz>0.4596))])
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster core UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)
    
    plt.subplot(3,2,2)
    #plt.scatter(rf_VJ[np.where(np.logical_and(z<0.4628,z>0.5028))],rf_UV[np.where(np.logical_and(z<0.4628,z>0.5028))],c=SFR[np.where(np.logical_and(z<0.4628,z>0.5028))])
    plt.scatter(crf_VJ[np.where(cz<0.4596)],crf_UV[np.where(cz<0.4628)],c=cSFR[np.where(cz<0.4596)])
    plt.scatter(crf_VJ[np.where(cz>0.4996)],crf_UV[np.where(cz>0.4996)],c=cSFR[np.where(cz>0.4996)])
    plt.scatter(grf_VJ[np.where(gz<0.4596)],grf_UV[np.where(gz<0.4628)],c=gSFR[np.where(gz<0.4596)])
    plt.scatter(grf_VJ[np.where(gz>0.4996)],grf_UV[np.where(gz>0.4996)],c=gSFR[np.where(gz>0.4996)])
    cbar = plt.colorbar()
    cbar.set_clim(0,1)
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Field UVJ')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)
    
    plt.subplot(3,2,3)
    plt.scatter(grf_VJ[np.where(np.logical_and(gz<0.4996,gz>0.4596))],grf_UV[np.where(np.logical_and(gz<0.4996,gz>0.4596))],c=gSFR[np.where(np.logical_and(gz<0.4996,gz>0.4596))])
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster groups UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)
    

    
    plt.subplot(3,2,4)
    plt.scatter(crf_m[np.where((cz>0.4596)&(cz<0.4996))],cSFR[np.where((cz>0.4596)&(cz<0.4996))],color='blue',alpha=0.7,marker = 'o')
    plt.scatter(grf_m[np.where((gz>0.4596)&(gz<0.4996))],gSFR[np.where((gz>0.4596)&(gz<0.4996))],color='blue',alpha=0.5,marker = 'x')
    plt.scatter(crf_m[np.where(cz<0.4596)],cSFR[np.where(cz<0.4596)],c='red',alpha=0.5)
    plt.scatter(crf_m[np.where(cz>0.4996)],cSFR[np.where(cz>0.4996)],c='red',alpha=0.5)
    plt.scatter(grf_m[np.where(gz<0.4596)],gSFR[np.where(gz<0.4596)],c='red',alpha=0.5)
    plt.scatter(grf_m[np.where(gz>0.4996)],gSFR[np.where(gz>0.4996)],c='red',alpha=0.5) 
    plt.xlim(7.5,11)
    #plt.ylim(-1,3)
    plt.xlabel('logM')
    plt.ylabel('logSFR')
    #plt.title('All pointings matched with LPD Halpha>0')
    red_patch = mpatches.Patch(color='red', label='Field', alpha = 0.5)
    blue_patch = mpatches.Patch(color='blue', label='Core', alpha = 0.5)
    blue_star = mlines.Line2D([], [], color='blue', marker='x', linestyle='None',
                          markersize=10, label='Groups')
    plt.legend(handles=[red_patch,blue_patch,blue_star],loc=1)

    #plt.subplot(3,2,6)
    
    #zc       = z[np.where((z>0.4596)&(z<0.4996))]
    
    #zf1       = z[np.where(z<0.4596)]   
    #zf2       = z[np.where((z>0.4996))]

    #x1 = zc
   # x2 = zf1
   # x3 = zf2
    
    
   # kwargs1 = dict(histtype='stepfilled', alpha=0.3, normed=False, bins=2)
   # kwargs2 = dict(histtype='stepfilled', alpha=0.3, normed=False, bins=5)
   # kwargs3 = dict(histtype='stepfilled', alpha=0.3, normed=False, bins=5)



   # 
  # plt.hist(x2,color='red', **kwargs2)
   # plt.hist(x3,color='red', **kwargs3);
   # plt.xlabel('z')
   # 
    plt.show()


def c1059(cf_Ha,gf_Ha,crf_VJ,crf_UV,crf_m,cz):
    plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/c-1059-gprior.txt', 'r')
    lines = f.readlines()[1:]
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
        crf_VJ.append(float(a[6]))
        crf_UV.append(float(a[7]))
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


    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/g-1059-gprior.txt', 'r')
    lines = f.readlines()[1:]
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
        grf_VJ.append(float(a[6]))
        grf_UV.append(float(a[7]))
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


    x1 = [-5,0.85]
    x2 = [0.85,1.6]
    x3 = [1.6,1.6]
    y1 = [1.3,1.3]
    y2 = [1.3,2]
    y3 = [2,5]
    
    plt.subplot(3,2,1)
    plt.scatter(crf_VJ[np.where((cz<0.4764)&(cz>0.4364))],crf_UV[np.where((cz<0.4764)&(cz>0.4364))],c=cSFR[np.where(np.logical_and(cz<0.4764,cz>0.4364))])
    cbar = plt.colorbar()
    cbar.set_clim(-1.5,1.5)
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster core UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)
    
    plt.subplot(3,2,2)
    #plt.scatter(crf_VJ[np.where((cz>0.4764)&(cz<0.4364))],crf_UV[np.where((cz>0.4764)&(cz<0.4364))],c=cSFR[np.where((cz>0.4764)&(cz<0.4364))])
    #plt.scatter(grf_VJ[np.where((gz>0.4764)&(gz<0.4364))],grf_UV[np.where((gz>0.4764)&(gz<0.4364))],c=gSFR[np.where((gz>0.4764)&(gz<0.4364))])
    plt.scatter(crf_VJ[np.where(cz<0.4364)],crf_UV[np.where(cz<0.4364)],c=cSFR[np.where(cz<0.4364)])
    plt.scatter(crf_VJ[np.where(cz>0.4764)],crf_UV[np.where(cz>0.4764)],c=cSFR[np.where(cz>0.4764)])
    plt.scatter(grf_VJ[np.where(gz<0.4364)],grf_UV[np.where(gz<0.4364)],c=gSFR[np.where(gz<0.4364)])
    plt.scatter(grf_VJ[np.where(gz>0.4764)],grf_UV[np.where(gz>0.4764)],c=gSFR[np.where(gz>0.4764)])
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    cbar.set_clim(-1.5,1.5)
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Field UVJ')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)

    plt.subplot(3,2,3)
    plt.scatter(grf_VJ[np.where((gz<0.4764)&(gz>0.4364))],grf_UV[np.where((gz<0.4764)&(gz>0.4364))],c=gSFR[np.where(np.logical_and(gz<0.4764,gz>0.4364))])
    cbar = plt.colorbar()
    cbar.set_clim(-1.5,1.5)
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster groups UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)

    plt.subplot(3,2,4)
    plt.scatter(crf_m[np.where((cz>0.4364)&(cz<0.4764))],cSFR[np.where((cz>0.4364)&(cz<0.4764))],color='blue',alpha=0.7,marker = 'o')
    plt.scatter(grf_m[np.where((gz>0.4364)&(gz<0.4764))],gSFR[np.where((gz>0.4364)&(gz<0.4764))],color='blue',alpha=0.5,marker = 'x')
    plt.scatter(crf_m[np.where(cz<0.4364)],cSFR[np.where(cz<0.4364)],c='red',alpha=0.5)
    plt.scatter(crf_m[np.where(cz>0.4764)],cSFR[np.where(cz>0.4764)],c='red',alpha=0.5)
    plt.scatter(grf_m[np.where(gz<0.4364)],gSFR[np.where(gz<0.4364)],c='red',alpha=0.5)
    plt.scatter(grf_m[np.where(gz>0.4764)],gSFR[np.where(gz>0.4764)],c='red',alpha=0.5)  
    plt.xlim(7.5,11.5)
    #plt.ylim(-1,3)
    plt.xlabel('logM')
    plt.ylabel('logSFR')
    #plt.title('All pointings matched with LPD Halpha>0')
    red_patch = mpatches.Patch(color='red', label='Field', alpha = 0.5)
    blue_patch = mpatches.Patch(color='blue', label='Core', alpha = 0.5)
    blue_star = mlines.Line2D([], [], color='blue', marker='x', linestyle='None',
                          markersize=10, label='Groups')
    plt.legend(handles=[red_patch,blue_patch,blue_star],loc=2)
    plt.show()


def c1227(cf_Ha,gf_Ha,crf_VJ,crf_UV,crf_m,cz):
    plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/c-1227-gprior.txt', 'r')
    lines = f.readlines()[1:]
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
        crf_VJ.append(float(a[6]))
        crf_UV.append(float(a[7]))
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


    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/g-1227-gprior.txt', 'r')
    lines = f.readlines()[1:]
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
        grf_VJ.append(float(a[6]))
        grf_UV.append(float(a[7]))
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

    x1 = [-5,0.85]
    x2 = [0.85,1.6]
    x3 = [1.6,1.6]
    y1 = [1.3,1.3]
    y2 = [1.3,2]
    y3 = [2,5]

    plt.subplot(3,2,1)
    plt.scatter(crf_VJ[np.where(np.logical_and(cz<0.6557,cz>0.6157))],crf_UV[np.where(np.logical_and(cz<0.6557,cz>0.6157))],c=cSFR[np.where(np.logical_and(cz<0.6557,cz>0.6157))])
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster core UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)
    
    plt.subplot(3,2,2)
    #plt.scatter(rf_VJ[np.where(np.logical_and(z<0.4628,z>0.5028))],rf_UV[np.where(np.logical_and(z<0.4628,z>0.5028))],c=SFR[np.where(np.logical_and(z<0.4628,z>0.5028))])
    plt.scatter(crf_VJ[np.where(cz<0.6157)],crf_UV[np.where(cz<0.6157)],c=cSFR[np.where(cz<0.6157)])
    plt.scatter(crf_VJ[np.where(cz>0.6557)],crf_UV[np.where(cz>0.6557)],c=cSFR[np.where(cz>0.6557)])
    plt.scatter(grf_VJ[np.where(gz<0.6157)],grf_UV[np.where(gz<0.6157)],c=gSFR[np.where(gz<0.6157)])
    plt.scatter(grf_VJ[np.where(gz>0.6557)],grf_UV[np.where(gz>0.6557)],c=gSFR[np.where(gz>0.6557)])
    cbar = plt.colorbar()
    cbar.set_clim(0,0.5)
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Field UVJ')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)

    plt.subplot(3,2,3)
    plt.scatter(grf_VJ[np.where(np.logical_and(gz<0.6557,gz>0.6157))],grf_UV[np.where(np.logical_and(gz<0.6557,gz>0.6157))],c=gSFR[np.where(np.logical_and(gz<0.6557,gz>0.6157))])
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster groups UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)

    plt.subplot(3,2,4)
    plt.scatter(crf_m[np.where((cz>0.6157)&(cz<0.6557))],cSFR[np.where((cz>0.6157)&(cz<0.6557))],color='blue',alpha=0.7,marker = 'o')
    plt.scatter(grf_m[np.where((gz>0.6157)&(gz<0.6557))],gSFR[np.where((gz>0.6157)&(gz<0.6557))],color='blue',alpha=0.5,marker = 'x')
    plt.scatter(crf_m[np.where(cz<0.6157)],cSFR[np.where(cz<0.6157)],c='red',alpha=0.5)
    plt.scatter(crf_m[np.where(cz>0.6557)],cSFR[np.where(cz>0.6557)],c='red',alpha=0.5)
    plt.scatter(grf_m[np.where(gz<0.6157)],gSFR[np.where(gz<0.6157)],c='red',alpha=0.5)
    plt.scatter(grf_m[np.where(gz>0.6557)],gSFR[np.where(gz>0.6557)],c='red',alpha=0.5)   
    plt.xlim(7.5,11)
    #plt.ylim(-1,3)
    plt.xlabel('logM')
    plt.ylabel('logSFR')
    #plt.title('All pointings matched with LPD Halpha>0')
    red_patch = mpatches.Patch(color='red', label='Field', alpha = 0.5)
    blue_patch = mpatches.Patch(color='blue', label='Core', alpha = 0.5)
    blue_star = mlines.Line2D([], [], color='blue', marker='x', linestyle='None',
                          markersize=10, label='Groups')
    plt.legend(handles=[red_patch,blue_patch,blue_star],loc=2)

    plt.show()
    
