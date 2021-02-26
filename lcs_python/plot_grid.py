# -*- coding: utf-8 -*-
"""
Created on Mon May 23 11:45:30 2016

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
sigma = F[:,2];
ix = F[:,0];
iy = F[:,1];
FTLEmaxi = np.amax(sigma);
FTLEmini = np.amin(sigma);

name2 = 'output_alpha_plot';
A = np.genfromtxt(name2,comments = '#', dtype = np.float64);
print 'read file', name2
print ''
alpha = A[:,2];
Alphamaxi = np.amax(alpha);
Alphamini = np.amin(alpha);


name3 = 'output_initcond_plot';
L1 = np.genfromtxt(name3,comments = '#', dtype = np.float64);
print 'read file', name3
print ''
lambda1 = L1[:,2];
L1maxi = np.amax(lambda1)
L1mini = np.amin(lambda1)

name4 = 'output_lambda2_plot';
L2 = np.genfromtxt(name4,comments = '#', dtype = np.float64);
print 'read file', name4
print ''
lambda2 = L2[:,2];
L2maxi = np.amax(lambda2)
L2mini = np.amin(lambda2)


print 'FTLE_max =', FTLEmaxi
print 'FTLE_min =', FTLEmini
print 'Alphamaxi =', Alphamaxi
print 'Alphamini =', Alphamini
print 'L1maxi =', L1maxi
print 'L1mini =', L1mini
print 'L2maxi =', L2maxi
print 'L2mini =', L2mini

x = F[:,0];
y = F[:,1];
###########################################################################################################
###########################################################################################################
###########################################################################################################




###########################################################################################################
###########################################################################################################
###########################################################################################################
NX = 1300j;
NY = 500j;

print 'interpolation'
print 'NX =', NX
print 'NY =', NY
print '' 

xmax = max(ix)
xmin = min(ix)
ymax = max(iy)
ymin = min(iy)

extent = (xmin,xmax,ymin,ymax);

xs,ys = np.mgrid[0:xmax:NX, 0:ymax:NY];
resampled_sigma = ml.griddata(x,y, sigma, xs,ys, interp = 'linear');
resampled_alpha = ml.griddata(x,y, alpha, xs,ys, interp = 'linear');
resampled_lambda1 = ml.griddata(x,y, lambda1, xs,ys, interp = 'linear');
resampled_lambda2 = ml.griddata(x,y, lambda2, xs,ys, interp = 'linear');

resampled_sigma_maxi = np.amax(resampled_sigma)
resampled_alpha_maxi = np.amax(resampled_alpha)
resampled_lambda1_maxi = np.amax(resampled_lambda1)
resampled_lambda2_maxi = np.amax(resampled_lambda2)
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
plt.imshow(resampled_sigma.T, aspect= 'equal', origin='lower',vmin=0,vmax=resampled_sigma_maxi, cmap=cm.rainbow, extent=extent);
plt.title('FTLE-Field');
plt.colorbar(orientation = 'horizontal');
plt.ylabel('y',fontsize=20);
plt.xlabel('x',fontsize=20);
######## ALPHA ########
plt.subplot(222);
plt.imshow(resampled_alpha.T, aspect= 'equal', origin='lower',vmin=0,vmax=resampled_alpha_maxi, cmap=cm.hot, extent=extent);
plt.title(r'$\alpha(x,y)$')
plt.colorbar(orientation = 'horizontal')
plt.ylabel('y',fontsize=20);
plt.xlabel('x',fontsize=20);
######## conditionAB ########
plt.subplot(223)
plt.imshow(resampled_lambda1.T, aspect= 'equal', origin='lower',vmin=0,vmax=resampled_lambda1_maxi, cmap=cm.bone, extent=extent);
plt.title('condition (A) and (B)')
plt.colorbar(orientation = 'horizontal')
plt.ylabel('y',fontsize=20)
plt.xlabel('x',fontsize=20)
######## lambda2 #######
plt.subplot(224)
plt.imshow(resampled_lambda2.T, aspect= 'equal', origin='lower',vmin=0,vmax=resampled_lambda2_maxi, cmap=cm.rainbow, extent=extent);
plt.title(r'$\lambda_2$')
plt.colorbar(orientation = 'horizontal')
plt.ylabel('y',fontsize=20)
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
