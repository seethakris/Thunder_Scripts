# -*- coding: utf-8 -*-
"""
Get user input and Run Thunder analysis
"""

#User input for experiment folder
Exp_Folder = '/Users/seetha/Desktop/Thunder/Data/Dorsal_Raphe_Gcamp_Tph2/Registered/'

#EOther experiment parameters
num_time = 285  #Number of time points
num_time_req  = 205 #Number of time points to include in analysis
img_size_x = 126 #X and Y resolution - if there are images that dont have this resolution, 
img_size_y = 260 #they will be resized
f_f_flag = 0 #0-raw data, 1-delta f/f
filename_suffix = 'Allfiles' #Suffix to save different files with
pca_components = 2 #Number of pca components to detect from files
num_pca_colors = 300 #Number of colors on the pca maps
thresh_pca = 0.0001 #Threshold above which to plot the pca components
num_samples = 10000 #number of random samples to select to do reconstruction
color_map = 'polar' #Colormap for plotting principle components

#Import python libraries
import numpy as np
import os
filesep = os.path.sep #File seperator by operating system
import time
import scipy 

#Import custom files
from create_data_for_thunder import create_textfile #To create text files
from run_thunder_analysis import run_pca, make_pca_maps
from plot_thunder_raphe import plot_pca_figures
#Import thunder libraries
from thunder.utils import ThunderContext

### ~~~~~~~~~~~~~~~~~~~ Main Script ~~~~~~~~~~~~~~~~~~~~~ ###

############### STEP 1 ######################
#Create text file 
#Check if text file already exists, else create it
txt_file = [f for f in os.listdir(Exp_Folder) if (f.endswith('.txt') and f.find(filename_suffix+'.txt')==0)]

if len(txt_file)==0: 
    start_time = time.time() 
    print 'Saving images to text on '+Exp_Folder
    exp_filenames = create_textfile(Exp_Folder, filename_suffix, num_time, num_time_req, img_size_x, img_size_y, f_f_flag)
    print 'Saving to text file took '+ str(int(time.time()-start_time)) +' seconds'
else:
    exp_filenames = [ f for f in os.listdir(Exp_Folder) if (os.path.isfile(os.path.join(Exp_Folder,f)) and f.find('.tif')>0)]


############## STEP 2 ######################
#Start Thunder Context if it doesnt already exist

print 'Starting Thunder Now. Check console for details'
tsc = ThunderContext.start(appName="pca")
time.sleep(2)
   
############## STEP 3 ######################  
#Load data from the text file using thunder context
data = tsc.loadText(Exp_Folder+filename_suffix+'.txt').cache()

############## STEP 3 ######################  
#Run PCA      
start_time = time.time()   
print 'Running pca...on '+Exp_Folder
pca, imgs_pca = run_pca(data,pca_components)                
print 'Running PCA took '+ str(int(time.time()-start_time)) +' seconds'

#Create polar maps
start_time = time.time()   
print 'Creating polar maps...on '+ Exp_Folder
maps, pts, clrs, recon, tt, unique_clrs, matched_pixels, matched_signals, mean_signal, sem_signal = make_pca_maps(pca, imgs_pca, img_size_x,\
img_size_y, num_pca_colors, num_samples, thresh_pca, color_map)
print 'Creating polar maps took '+ str(int(time.time()-start_time)) +' seconds'

#Save pca components, scores, maps and clrs
#Create directory first
Np_matfiles_Directory = Exp_Folder+'Np_mat_files'+filesep
if not os.path.exists(Np_matfiles_Directory ):
    os.makedirs(Np_matfiles_Directory) 

#Save as matfile
scipy.io.savemat(Np_matfiles_Directory+filename_suffix+'_PCA.mat', dict(comp=pca.comps.T,\
latent=pca.latent,scores=imgs_pca, color_maps=maps,pca_scatter_pts=pts, pts_clrs=clrs, \
unique_clrs=unique_clrs, matched_pixels=matched_pixels,\
mean_signal=mean_signal, sem_signal=sem_signal))
#Save as numpy array
np.savez(Np_matfiles_Directory+filename_suffix+'_PCA.npz', pca.comps.T,pca.latent,imgs_pca,maps, pts,recon, clrs,\
unique_clrs, matched_pixels, mean_signal, sem_signal)

#Plot PCA components and maps
start_time = time.time()   
print 'Plotting PCA figures...on '+Exp_Folder
plot_pca_figures(pca, maps, pts, clrs, recon,tt, unique_clrs, matched_pixels, matched_signals,Exp_Folder, filename_suffix, exp_filenames)        
print 'Plotting PCA figures took '+ str(int(time.time()-start_time)) +' seconds'


############## STEP 4 ######################  
#Run ICA
