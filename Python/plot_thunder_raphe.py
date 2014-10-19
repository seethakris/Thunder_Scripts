# -*- coding: utf-8 -*-
"""
Plot PCA and ICA components for Raphe - Gcamp6f and tph2
"""
#Import python libraries
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns #For creating nice plots
import pandas as pd

filesep = os.path.sep

#Plot various PCA figures    
def plot_pca_figures(pca, maps, pts, clrs, recon,tt,unique_clrs, matched_pixels,matched_signals, Exp_Folder, filename_suffix, Exp_Name):  

    #Plotting as pdf
    Figure_PDFDirectory = Exp_Folder+'Figures'+filesep
    if not os.path.exists(Figure_PDFDirectory):
        os.makedirs(Figure_PDFDirectory)           
    pp = PdfPages(Figure_PDFDirectory+filename_suffix+'_PCA.pdf')

    sns.set_context("poster")  
    

############ Plot Colormaps of scores ############
    
    #If there is only one stack, else plot each stack
    if len(maps.shape)==3:
        #Plot colored maps for each stack
        with sns.axes_style("white"):
            fig2 = plt.imshow(maps[:,:,:].transpose((1,0,2)))
            plt.title(Exp_Name[0])
            fig2 = plt.gcf()
            pp.savefig(fig2)
            plt.close()

    else:
        for ii in range(0,np.size(maps,2)):
            with sns.axes_style("white"):
                fig2 = plt.imshow(maps[:,:,ii,:].transpose((1,0,2)))
                plt.title(Exp_Name[ii])
                fig2 = plt.gcf()
                pp.savefig(fig2)
                plt.close()
                

    ########### Plot components ##################
    fig2 = plt.figure()
    sns.set_context("talk", font_scale=1.25)
    with sns.axes_style("dark"):
        ax1 = plt.subplot(221)
        plt.plot(pca.comps.T);
        plt.locator_params(axis = 'y', nbins = 4)
        sns.axlabel("Time (seconds)","a.u")
        A = []
        for ii in xrange(0,np.size(pca.comps.T, 0)):
            A = np.append(A, ['comp' + str(ii+1)])
        ax1.legend(A, loc=4)
        plot_vertical_lines()
        plt.axhline(y=0, linestyle='-', color='k', linewidth=1)
        
    #Plot mean signals according to color and boxplot of number of pixels in each plane
    with sns.axes_style("dark"):
        for ii in range(0,np.size(unique_clrs,0)):       
            fig2 = plt.subplot(223)
            sns.tsplot(np.array(matched_signals[ii].clr_grped_signal), linewidth=3, ci=95, err_style="ci_band", color=unique_clrs[ii])
            plt.locator_params(axis = 'y', nbins = 4)            
            sns.axlabel("Time (seconds)","a.u")            
            plot_vertical_lines()
        plt.axhline(y=0, linestyle='-', color='k', linewidth=1)

    
    with sns.axes_style("white"):
        fig2 = plt.subplot(247)
        fig2 = sns.boxplot(np.transpose(matched_pixels[:,0:6]/220),linewidth=3, widths=.5, color=unique_clrs)
        for ii in range(0,np.size(unique_clrs,0)):
            fig2 = plt.plot(np.repeat(ii+1,np.size(matched_pixels[:,0:6],1)), np.transpose(matched_pixels[ii,0:6]/220),'s', \
            color=unique_clrs[ii], markersize=5, markeredgecolor='k', markeredgewidth=2) 
            plt.locator_params(axis = 'y', nbins = 4)
        sns.axlabel("Colors", "Number of Pixels")
        plt.title('Gcamp6f')
        sns.despine(offset=10, trim=True);    
        
        fig2 = plt.subplot(248)
        fig2 = sns.boxplot(np.transpose(matched_pixels[:,6:]/80),linewidth=3, widths=.5, color=unique_clrs)
        for ii in range(0,np.size(unique_clrs,0)):
            fig2 = plt.plot(np.repeat(ii+1,np.size(matched_pixels[:,6:],1)), np.transpose(matched_pixels[ii,6:]/220),'s', \
            color=unique_clrs[ii], markersize=5, markeredgecolor='k', markeredgewidth=2) 
            plt.locator_params(axis = 'y', nbins = 4)
        plt.title('tph2')    
        sns.axlabel("Colors", "Number of Pixels")
        sns.despine(offset=10, trim=True);   
        
            #Plot mean projection        
    with sns.axes_style("white"):  
        temp = (np.max(maps[:,:,0:6,:], axis=2))
        fig2 = plt.subplot(243)
        plt.imshow(temp.astype(np.float16))
        plt.axis('off')
        plt.title('Max Gcamp6f')
        
        temp = (np.max(maps[:,:,6:,:], axis=2))
        fig2 = plt.subplot(244)
        plt.imshow(temp.astype(np.float16))
        plt.axis('off')
        plt.title('Max tph2')
        
        plt.tight_layout()
        fig2 = plt.gcf()
        pp.savefig(fig2)
        plt.close()
    
    matched_pixels1 = np.delete(matched_pixels, 3,0)    
    A = np.zeros((3,np.size(matched_pixels1,0)*np.size(matched_pixels1,1)), dtype=np.float)
    count = 0
    for ii in range(0,np.size(matched_pixels1,0)):
        for jj in range(0,np.size(matched_pixels1,1)):
            if ii == 0 or ii == 1 or ii == 4:
                A[1,count] = 1
            else:
                A[1,count] = 0     
    
            if jj < 6:
                A[2,count] = 0
                A[0,count] = matched_pixels1[ii,jj]/220
            else:
                A[2,count] = 1
                A[0,count] = matched_pixels1[ii,jj]/80
            count = count+1 
    A = np.transpose(A)        
    B = pd.DataFrame({'Cells':A[:,0], 'response':A[:,1], 'Gcamportph2':A[:,2]})
    B["Gcamptph"] = B.Gcamportph2.map({0: "Gcamp6", 1: "tph2"})
    B["Response"] = B.response.map({0: "Excitatory On", 1: "Inhibitory On"})
    fig2 = plt.figure()
    fig2 = sns.factorplot("Response", "Cells", "Gcamptph", B)
    fig2 = plt.gcf()
    pp.savefig(fig2)
    plt.close()

    pp.close()
    
def plot_vertical_lines():
    plt.axvline(x=46, linestyle='-', color='k', linewidth=1)
    plt.axvline(x=66, linestyle='--', color='k', linewidth=1)
    plt.axvline(x=86, linestyle='-', color='k', linewidth=1)
    plt.axvline(x=106, linestyle='--', color='k', linewidth=1)
    plt.axvline(x=126, linestyle='-', color='k', linewidth=1)
    plt.axvline(x=146, linestyle='--', color='k', linewidth=1)
    plt.axvline(x=166, linestyle='-', color='k', linewidth=1)
    plt.axvline(x=186, linestyle='--', color='k', linewidth=1)