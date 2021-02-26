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
name1 = 'output_eigv1_x_plot'
A = np.genfromtxt(name1,comments = '#', dtype = np.float64);
print 'read', name1
print '' 
x = A[:,0];
y = A[:,1];
V2X = A[:len(x),2];

name2 = 'output_eigv1_y_plot' 
B = np.genfromtxt(name2,comments = '#', dtype = np.float64);
print 'read', name2 
print ''
V2Y = B[:,2];

name3 = 'output_initcond_plot'
C = np.genfromtxt(name3,comments = '#',dtype = np.float64);
print 'read', name3 
print ''
sigma = C[:,2];

name4 = 'output_lcs'
D = np.genfromtxt(name4, comments = '#',dtype = np.float64)
print 'read', name4
print ''
Rx = D[:,0];
Ry = D[:,1];

xmax = max(x)
xmin = min(x)

ymax = max(y)
ymin = min(y)
###########################################################################################################
###########################################################################################################
###########################################################################################################




###########################################################################################################
###########################################################################################################
###########################################################################################################
# ploty plot!
plt.figure(1)
#plt.quiver(x,y,V2X,V2Y,sigma,units = 'xy',scale = 90,width = 0.0001,cmap = cm.rainbow)
plt.quiver(x,y,V2X,V2Y,sigma,units = 'xy',scale = 50,width = 0.001,cmap = cm.rainbow)
#plt.quiver(x,y,V2X,V2Y,sigma,units = 'xy',scale = 2,width = 0.04,cmap = cm.rainbow)
#plt.plot(Rx,Ry,'ro')
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.title(r'$\xi_1$',fontsize=30)
plt.ylabel('y',rotation=0,fontsize=30);
plt.xlabel('x',fontsize=24);
plt.savefig('lcsxi1_gyre_stat.png',orientation = 'portrait',dpi = 150);
plt.show()
###########################################################################################################
###########################################################################################################
###########################################################################################################




print '########################################'
print 'end'
print '########################################'
