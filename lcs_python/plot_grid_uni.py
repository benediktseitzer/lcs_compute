# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14:06:30 2016

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
sigma = F;
sigma_maxi = np.amax(sigma);
sigma_mini = np.amin(sigma);

name2 = 'output_alpha';
A = np.genfromtxt(name2,comments = '#', dtype = np.float64);
print 'read file', name2
print ''
alpha = A;
alpha_maxi = np.amax(alpha);
alpha_mini = np.amin(alpha);


name3 = 'output_initcond';
L1 = np.genfromtxt(name3,comments = '#', dtype = np.float64);
print 'read file', name3
print ''
lambda1 = L1;
lambda1_maxi = np.amax(lambda1)
lambda1_mini = np.amin(lambda1)

name4 = 'output_lambda2';
L2 = np.genfromtxt(name4,comments = '#', dtype = np.float64);
print 'read file', name4
print ''
lambda2 = L2;
lambda2_maxi = np.amax(lambda2)
lambda2_mini = np.amin(lambda2)


print 'FTLE_max =', sigma_maxi
print 'FTLE_min =', sigma_mini
print 'Alphamaxi =', alpha_maxi
print 'Alphamini =', alpha_mini
print 'L1maxi =', lambda1_maxi
print 'L1mini =', lambda1_mini
print 'L2maxi =', lambda2_maxi
print 'L2mini =', lambda2_mini

###########################################################################################################
###########################################################################################################
###########################################################################################################




###########################################################################################################
###########################################################################################################
###########################################################################################################

xmax = 2.*np.pi
xmin = 0.*np.pi
ymax = 2.*np.pi
ymin = 0.

extent = (xmin,xmax,ymin,ymax);

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
######## SIGMA ########
plt.subplot(221);
plt.imshow(sigma.T, aspect= 'equal', origin='lower',vmin=sigma_mini,vmax=sigma_maxi, cmap=cm.rainbow, extent=extent);
plt.title('FTLE-Field');
plt.colorbar(orientation = 'horizontal');
plt.ylabel('z',fontsize=20);
plt.xlabel('x',fontsize=20);
######## ALPHA ########
plt.subplot(222);
plt.imshow(alpha.T, aspect= 'equal', origin='lower',vmin=alpha_mini,vmax=alpha_maxi, cmap=cm.hot, extent=extent);
plt.title(r'$\alpha(x,y)$')
plt.colorbar(orientation = 'horizontal')
plt.ylabel('z',fontsize=20);
plt.xlabel('x',fontsize=20);
######## conditionAB ########
plt.subplot(223)
plt.imshow(lambda1.T, aspect= 'equal', origin='lower',vmin=lambda1_mini,vmax=lambda1_maxi, cmap=cm.rainbow, extent=extent);
plt.title('condition (A) and (B)')
plt.colorbar(orientation = 'horizontal')
plt.ylabel('z',fontsize=20)
plt.xlabel('x',fontsize=20)
######## lambda2 #######
plt.subplot(224)
plt.imshow(lambda2.T, aspect= 'equal', origin='lower',vmin=lambda2_mini,vmax=lambda2_maxi, cmap=cm.bone, extent=extent);
plt.title(r'$\lambda_2$')
plt.colorbar(orientation = 'horizontal')
plt.ylabel('z',fontsize=20)
plt.xlabel('x',fontsize=20)
######## show&save ########
plt.savefig('outputparas_gyre_stat.png',orientation = 'portrait',dpi = 150);
plt.show();
print 'figure saved'
###########################################################################################################
###########################################################################################################
###########################################################################################################




print '########################################'
print 'end'
print '########################################'
