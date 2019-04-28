import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import astropy 
from astropy import cosmology
from astropy.cosmology import WMAP9 as cosmo
import math as mt
import os, sys
import matplotlib.lines as mlines




def uvj(cf_Ha,gf_Ha,crf_VJ,crf_UV,crf_m,cz):
    plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/all_core.txt', 'r')
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
    print(np.mean(cSFR))
    print(np.std(cSFR))


    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/all_group.txt', 'r')
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
    print(np.mean(gSFR))
    print(np.std(gSFR))

    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/all_field.txt', 'r')
    lines = f.readlines()[1:]
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
        rf_VJ.append(float(a[6]))
        rf_UV.append(float(a[7]))
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
    print(np.mean(SFR))
    print(np.std(SFR))
    
    x1 = [-5,0.85]
    x2 = [0.85,1.6]
    x3 = [1.6,1.6]
    y1 = [1.3,1.3]
    y2 = [1.3,2]
    y3 = [2,5]

    f = open('/data2/jrcooper/notebooks/reduction/EDisCS/Ha_catalogs/gaussian_prior/1301_LDP.txt', 'r')
    lines = f.readlines()[1:]
    f.close()
    
    #create arrays 
    LDP_UV    = [] 
    LDP_VJ   = []
    
    
    #pull array column 
    for line in lines: 
        a = line.split()
        LDP_UV.append(float(a[3]))
        LDP_VJ.append(float(a[2]))
        
        
    LDP_UV = np.array(LDP_UV)
    LDP_VJ = np.array(LDP_VJ)

    
    plt.subplot(3,2,1)
    plt.hist2d(LDP_VJ-0.155, LDP_UV-0.306, bins=(150, 150), cmap=plt.cm.Greys)
    plt.scatter(crf_VJ,crf_UV,c=cSFR,s=100)
    cbar = plt.colorbar()
    cbar.set_clim(-0.5,1.2)
    cbar.set_label("log(SFR)")
    #plt.colorbar()

    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster core UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)
    
    plt.subplot(3,2,2)
    plt.hist2d(LDP_VJ-0.155, LDP_UV-0.306, bins=(150, 150), cmap=plt.cm.Greys)
    plt.scatter(grf_VJ,grf_UV,c=gSFR,s=100)
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    cbar.set_clim(-0.5,1.2)
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Cluster group UVJ z+/- 0.02')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)

    plt.subplot(3,2,3)
    plt.hist2d(LDP_VJ-0.155, LDP_UV-0.306, bins=(150, 150), cmap=plt.cm.Greys)
    plt.scatter(rf_VJ,rf_UV,c=SFR,s=100)
    cbar = plt.colorbar()
    cbar.set_label("log(SFR)")
    cbar.set_clim(-0.5,1.2)
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    plt.title('Field UVJ')
    plt.plot(x1,y1,color='black')
    plt.plot(x2,y2,color='black')
    plt.plot(x3,y3,color='black')
    plt.xlim(-0.5,2)
    plt.ylim(0,3)
   
    
    plt.subplot(3,2,4)
    plt.scatter(crf_m,cSFR,color='blue',alpha=0.7,marker = 'o',s=100)
    plt.scatter(grf_m,gSFR,color='blue',alpha=0.5,marker = 'x',s=100)
    plt.scatter(rf_m,SFR,c='red',alpha=0.5,s=100) 
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
    


