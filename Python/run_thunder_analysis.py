# -*- coding: utf-8 -*-
"""
Run thunder analysis like pca, ica and get data 
"""
#Import python libraries
import numpy as np
import scipy as sp

#Import thunder libraries
from thunder.factorization import PCA
from thunder.viz import Colorize
from thunder.utils import subset, pack

#Run PCA and pack scores into images
def run_pca(data, k):   
    #PCA with k components
    pca = PCA(k=k).fit(data)    
    #Convert scores to images
    imgs = pack(pca.scores)
        
    return pca, imgs

#Make maps of the pca scores as polar colors    
def make_pca_maps(pca, imgs, img_size_x, img_size_y, scale_num, nsamples, thresh_num, color_maps):
            
    #Get polorized maps of first two components
            
    maps = Colorize(color_maps, scale=scale_num).image(imgs)
    #Get number of planes based on map dimesnions
    if len(maps.shape)==3:
        num_planes = 1
    else:
        num_planes = np.size(maps,2)
    num_time = np.size(pca.comps.T,0)
        
    pts = subset(pca.scores, nsamples=nsamples, thresh=thresh_num)
    clrs = Colorize(color_maps, scale=scale_num).points(pts)
    
    #Get distribution and reconstruction using first two components

    recon = map(lambda x: (x[0] * pca.comps[0, :] + x[1] * pca.comps[1, :]).tolist(), pts)
    tt = range(0, np.shape(pca.comps[0,:])[0])
    
    #Get specific color matches across animals and get mean and standard deviation
    array = [map(int,single_dim) for single_dim in clrs] #Convert the colors to RGB integers
    new_array = [tuple(row) for row in array] 
    unique_clrs = list(set(new_array))    #Get unique combination of colors
    unique_clrs.remove((0,0,0))
    matches = [np.where((np.array(array) == match).all(axis=1)) for match in unique_clrs] #Match the colors with the original rows
    
    #From maps get number of pixel matches with color for each plane
    array_maps = maps.astype(np.int64)
    matched_pixels = np.zeros((np.size(unique_clrs,0),num_planes))
    if len(maps.shape) == 3:
        array_maps_plane = np.reshape(array_maps, np.size(array_maps,0)*np.size(array_maps,1),3)
    else:     
        for ii in xrange(0,num_planes):
            array_maps_plane = np.reshape(array_maps[:,:,ii,:], (np.size(array_maps,0)*np.size(array_maps,1),3))
            matched_pixels[:,ii] = [np.size(np.where((np.array(array_maps_plane) == match).all(axis=1))) for match in unique_clrs]
             
    
    #Find stats based on the color
    matched_signals = [ structtype() for i in range(np.size(matches,0)*num_planes)]
    
    mean_signal = np.zeros((np.size(matches,0), num_planes, num_time))
    sem_signal = np.zeros((np.size(matches,0), num_planes, num_time))
    for ii in xrange(0,np.size(matches,0)):
            temp_ele = np.array(matches[ii])
            matched_signals[ii].clr_grped_signal = [np.array(recon[ele]) for ele in temp_ele[0,:]] #Get signals from the reconstruction that match the colors                     
            mean_signal[ii,:] = np.mean(matched_signals[ii].clr_grped_signal,axis=0) 
            sem_signal[ii,:] = sp.stats.sem(matched_signals[ii].clr_grped_signal,axis=0)  
            
        
    return maps, pts, clrs, recon, tt, unique_clrs, matched_pixels, matched_signals, mean_signal, sem_signal
    
class structtype():
    pass