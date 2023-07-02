#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:30:45 2021

@author: eveshalom
"""
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.widgets import LassoSelector
from matplotlib import path

def ROILasso(MaxSlice, datashape, z):
    
    plt.ion()
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.set_title('lasso selection:')
    ax1.imshow(MaxSlice,vmin=np.percentile(MaxSlice,70),vmax=np.percentile(MaxSlice,95))
    ax1.set_aspect('equal')
    
    # Empty array to be filled with lasso selector
    array = np.zeros_like(MaxSlice)
    ax2 = fig.add_subplot(122)
    ax2.set_title('numpy array:')
    mask = ax2.imshow(array,vmax=1, interpolation='nearest')
    
    plt.show()
    # Pixel coordinates
    pix = np.arange(MaxSlice.shape[1])
    xv, yv = np.meshgrid(pix,pix)
    pix = np.vstack( (xv.flatten(), yv.flatten()) ).T
    
    def updateArray(array, indices):
        lin = np.arange(array.size)
        newArray = array.flatten()
        newArray[lin[indices]] = 1
        return newArray.reshape(array.shape)

    def onselect(verts):
        array=mask._A._data
        p = path.Path(verts)
        ind = p.contains_points(pix, radius=1)
        array = updateArray(array, ind)
        mask.set_data(array)
        fig.canvas.draw_idle()
    
    lineprops = {'color': 'red', 'linewidth': 4, 'alpha': 0.8}
    lasso = LassoSelector(ax1, onselect,lineprops, button=1)

    return lasso,mask

def maskto4d(mask2d,datashape,z):
    mask2d=mask2d._A._data.copy()
    roivox = np.sum(mask2d)
    timemask=np.tile(mask2d[:,:,np.newaxis],(1,1,datashape[-1]))
    mask4d = np.zeros(datashape)
    mask4d[:,:,z,:]=timemask
    return mask4d, roivox


def ROISelector(Max,z,datashape,use_save,encdir):
    
    if use_save is True:
        tv = input("Which save roi would you like? Type aif or sag or tum \n")
        tv = str(tv)
        mask = np.load('{}/np_param_arrays/{}mask.npy'.format(encdir,tv))
        roivox = np.load('{}/np_param_arrays/{}maskvox.npy'.format(encdir,tv))
        return mask, roivox
    else:
        lasso_obj, maskarray = ROILasso(Max[:,:,z],datashape,z)
        mask,roivox = maskto4d(maskarray, datashape, z)
        return mask, roivox, lasso_obj

def ICfromROI(E,mask,roivox,numaxis):
    
    Eroi = (np.sum(mask*E,axis=tuple(range(0, numaxis))))/roivox # calculates average roi signal enhancement

    return Eroi