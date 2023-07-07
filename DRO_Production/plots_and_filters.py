#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:37:42 2021

@author: eveshalom
"""
from matplotlib import pyplot as plt
import numpy as np
from scipy import ndimage as ndimage

def checkROIs(Ea,Ev):
    plt.rcParams.update({'font.size':18})
    fig,[ax1,ax2] = plt.subplots(1,2,figsize=(20,10))
    ax1.plot(Ea)
    ax2.plot(Ev)
    ax1.set(ylabel='Signal Enhancement',title='Middle Cerebral Arteries ETC')
    ax2.set(ylabel='Signal Enhancement',title='Sagittal Sinus ETC')
    return

def MedianFilter(paramMap):
    for i in range(paramMap.shape[-1]):
        paramMap[:,:,i] = ndimage.median_filter(paramMap[:,:,i], size=(3,3))
    return paramMap

def displaySlices(data3d,outdir,fname,cbarlabel):
    plt.rcParams.update({'font.size':22})
    title = fname.replace("_", " ")
    fig, axs = plt.subplots(nrows=4,ncols=4,sharey=True,sharex=True,figsize=(22,20))
    fig.subplots_adjust(wspace=0.01)
    fig.subplots_adjust(hspace=0.1)
    vmin=np.percentile(data3d,2)
    vmax=np.percentile(data3d,98)
    quad = []
    z=0
    for ax in axs.flatten():
        im=ax.imshow(data3d[:,:,z],vmin=vmin,vmax=vmax)
        quad.append(im)
        z+=1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label('{}'.format(cbarlabel),fontsize=30)
    fig.suptitle(" \n {} ".format(title),fontsize=50)
    plt.savefig('{}/{}.png'.format(outdir,fname))
    return

def animateSingle(data,fname,fdir,datashape,dt,cbarlabel):
    from matplotlib.animation import FuncAnimation
    plt.rcParams.update({'font.size':22})
    title = fname.replace("_", " ")
    fig, axs = plt.subplots(nrows=4,ncols=4,sharey=True,sharex=True,figsize=(22,20))
    fig.subplots_adjust(wspace=0.01)
    fig.subplots_adjust(hspace=0.1)
    vmin=0
    vmax=np.percentile(data,98)#224
    quad = []
    for ax in axs.flatten():
        im=ax.imshow(data[:,:,0,0],vmin=vmin,vmax=vmax)
        quad.append(im)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label('{}'.format(cbarlabel),fontsize=30)
    fig.suptitle(" \n {} at T = 0 mins".format(title),fontsize=50)
    def init():
        z=0
        for pan in quad:
            pan.set_array(data[:,:,z,0])
            z+=1
        return quad
        
    def animate(i):
        z=0
        for pan in quad:
            pan.set_array(data[:,:,z,i])
            z+=1
        fig.suptitle(" \n {} at T = {:.2f} mins".format(title,i*dt),fontsize=50)
        return quad
    
    anim = FuncAnimation(fig=fig,  func=animate,frames=datashape[-1],init_func=init, interval=500, blit=True)
    anim.save('{}/{}.gif'.format(fdir,fname),writer='ffmpeg')
    return