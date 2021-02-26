# -*- coding: utf-8 -*-
"""
Created on Mon Nov 07 14:14:30 2016

@author: b_seit01
"""
print '########################################'
print 'start'
print '########################################'

import numpy as np
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
import matplotlib.cm as cm

 ###########################################################################################################
 ###########################################################################################################
 ###########################################################################################################
name1 = 'output_sigma_plot';
F = np.genfromtxt(name1,comments = '#', dtype = np.float64);
print 'read file', name1
print ''
x = F[:,0];
y = F[:,1];
sigma = F[:,2];
FTLEmaxi = np.amax(sigma);
FTLEmini = np.amin(sigma);

name2 = 'output_lcs';
G = np.genfromtxt(name2,comments = '#', dtype = np.float64);
print 'read file', name2
print ''
Rx = G[:,0];
Ry = G[:,1];

name3 = 'output_initcond_plot';
H = np.genfromtxt(name3,comments = '#', dtype = np.float64);
print 'read file', name3
print ''
lambda1 = H[:,2];
lambda1_maxi = np.amax(lambda1);
lambda1_mini = np.amin(lambda1);

print 'FTLE_max =', FTLEmaxi
print 'FTLE_min =', FTLEmini
print ''
 ###########################################################################################################
 ###########################################################################################################
 ###########################################################################################################



 ###########################################################################################################
 ###########################################################################################################
 ###########################################################################################################

NX = 2j*1000;
NY = 2j*200;

print 'interpolation'
print 'NX =', NX
print 'NY =', NY
print '' 

xmax = max(x)
xmin = min(x)
ymax = max(y)
ymin = min(y)

extent = (xmin,xmax,xmin,ymax);

xs,ys = np.mgrid[0:xmax:NX, 0:ymax:NY];
resampled = ml.griddata(x,y, sigma, xs,ys, interp = 'linear');
resampled_lambda1 = ml.griddata(x,y, lambda1, xs,ys, interp = 'linear');
 ###########################################################################################################
 ###########################################################################################################
 ###########################################################################################################



 ###########################################################################################################
 ###########################################################################################################
 ###########################################################################################################
print 'plot FTLE'
print''
# and now... plot!
plt.figure(1);

#plt.subplot(211);
plt.imshow(resampled.T, aspect= 'equal', origin='lower',vmin=FTLEmini,vmax=FTLEmaxi, cmap=cm.rainbow, extent=extent)
plt.plot(Rx,Ry,'ro',ms=2)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.title(r'$\sigma(x,y)$',fontsize=40);
plt.colorbar(orientation = 'horizontal');
plt.ylabel('y',rotation=0,fontsize=24);
plt.xlabel('x',fontsize=24);

#plt.subplot(212);
#plt.imshow(resampled_lambda1.T, aspect= 'equal', origin='lower',vmin=lambda1_mini,vmax=lambda1_maxi, cmap=cm.bone, extent=extent);
#plt.title('condition (A) and (B)')
#plt.colorbar(orientation = 'horizontal')
#plt.ylabel('z',fontsize=20)
#plt.xlabel('x',fontsize=20)

plt.savefig('lcsftle_gyre_stat.png',orientation = 'portrait',dpi = 150);
plt.show();
print 'figure saved'
 ###########################################################################################################
 ###########################################################################################################
 ###########################################################################################################



print '########################################'
print 'end'
print '########################################'

