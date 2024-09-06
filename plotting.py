import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import matplotlib.lines as lne
from matplotlib.colors import LinearSegmentedColormap

import numpy as np
import sys, os, glob,csv
import random
import colorsys

DDIR = './cfgs/'

# =============== plot the network ========================================

def plotNetwork(netfile):
    
    #try:
    network = np.loadtxt(netfile)
    L = network[0][0]
    Lx=L
    xtilt=L/2
    Ly = L*np.sqrt(3)/2.
    print(L,Ly,xtilt)
    fig=plt.figure(figsize=(Lx,Ly))
    axval= fig.add_subplot(1, 1, 1)
    network = network[1:]
    
    
    for pair in network:
        xu = pair[0]%L+0.5*(pair[0]//L)
        yu = (pair[0]//L)*np.sqrt(3)/2
        xv = pair[1]%L+0.5*(pair[1]//L)
        yv = (pair[1]//L)*np.sqrt(3)/2

        cir1=ptch.Circle((xu,yu),radius=.1,ec=None,fc='black')
        axval.add_patch(cir1)
        axval.text(xu+0.1,yu+0.1,str(pair[0].astype(int)),fontsize=8)
        cir2=ptch.Circle((xv,yv),radius=.1,ec=None,fc='black')
        axval.add_patch(cir2)
        axval.text(xv+0.1,yv+0.1,str(pair[1].astype(int)),fontsize=8)
        ####################
        if xu-xv>1 and yu-yv<1:
           xu=xu-L
        elif xv-xu>1 and yv-yu<1:
           xv=xv-Lx
        
        if yu-yv>1:
           yu=yu-Ly
           xu=xu-xtilt
        elif yv-yu>1:
           yv=yv-Ly
           xv=xv-xtilt
        
        edge=lne.Line2D([xu,xv],[yu,yv],color='black')              
        axval.add_line(edge)

    axval.axis('scaled')
    #axval.axis('off')
    return fig


# =============== configuration plotting ========================================
'''
Display the bonds saved in bondfile [save_config function in plotting.cpp]. 
You can highlight a cluster with given index (in red), and the largest cluster (in black)
'''
def plotConf(bondfile, cluster_ind=0, show_largest=True):
    
    #try:
    bonds = np.loadtxt(bondfile)
    Lx = bonds[0][0]
    Ly = Lx*np.sqrt(3)/2.
    xtilt = Lx/2.
    print("Lx, Ly, xtilt",Lx,Ly,xtilt)
    bonds=bonds[1:]

    fig=plt.figure(figsize=(Lx,Ly))
    axval= fig.add_subplot(1, 1, 1)
     
    max_index=max(bonds[:,4]).astype(int)
    rnd_cols = [ [random.random() for i in range(3)] for j in range(1,max_index+1)]

    values, counts = np.unique(bonds[:,4], return_counts=True)
    largest_index = values[counts.argmax()].astype(int) #index of the largest cluster


    for bond in bonds:
        xu,yu,xv,yv=bond[0:4]
        cir1=ptch.Circle((xu,yu),radius=.1,ec=None,fc='black')
        axval.add_patch(cir1)
        cir2=ptch.Circle((xv,yv),radius=.1,ec=None,fc='black')
        axval.add_patch(cir2)
        ####################
        if xu-xv>1 and yu-yv<1:
           xu=xu-Lx
        elif xv-xu>1 and yv-yu<1:
           xv=xv-Lx
        
        if yu-yv>1:
           yu=yu-Ly
           xu=xu-xtilt
        elif yv-yu>1:
           yv=yv-Ly
           xv=xv-xtilt
        
        
        
        if cluster_ind>0 and bond[4]==cluster_ind:
           edge=lne.Line2D([xu,xv],[yu,yv],linewidth=4,color='red') 
        elif show_largest and bond[4]==largest_index:
           edge=lne.Line2D([xu,xv],[yu,yv],linewidth=4,color='black')
        else:
           edge=lne.Line2D([xu,xv],[yu,yv],color=rnd_cols[bond[4].astype(int)-1])              
        axval.add_line(edge)

    axval.axis('scaled')
    #axval.axis('off')
    return fig




###############################################

if __name__ == "__main__":
    L = input("enter linear size ")
    k = input("enter step ")
    
    fig=plotConf(DDIR+'rigid_clusters_L'+L+'_step_'+k+'.txt', cluster_ind=147)
    fig = plotNetwork(DDIR+"network_L"+L+'_step_'+k+'.txt')
    plt.show()
    
