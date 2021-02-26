# -*- coding: utf-8 -*-
"""
reated on Thu Dec 15 14:29:30 2016

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
name1 = 'output_sigma';
F = np.genfromtxt(name1,comments = '#', dtype = np.float64);
print 'read file', name1
print ''
sigma = F
FTLEmaxi = np.amax(sigma);
FTLEmini = np.amin(sigma);


name2 = 'output_lcs';
G = np.genfromtxt(name2,comments = '#', dtype = np.float64);
print 'read file', name2
print ''
Rx = G[:,0];
Ry = G[:,1];

name3 = 'output_initcond';
H = np.genfromtxt(name3,comments = '#', dtype = np.float64);
print 'read file', name3
print ''
lambda1 = H
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

xmax = 2.*np.pi
xmin = 0.
ymax = 2.*np.pi
ymin = 0.

extent = (xmin,xmax,xmin,ymax);

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
plt.imshow(sigma.T, aspect= 'equal', origin='lower',vmin=0,vmax=500000, cmap=cm.rainbow, extent=extent);
plt.plot(Rx,Ry,'ro',ms=2)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.title('FTLE-Field');
plt.colorbar(orientation = 'horizontal');
plt.ylabel('z',fontsize=20);
plt.xlabel('x',fontsize=20);

#plt.subplot(212);
#plt.imshow(lambda1.T, aspect= 'equal', origin='lower',vmin=lambda1_mini,vmax=lambda1_maxi, cmap=cm.bone, extent=extent);
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
