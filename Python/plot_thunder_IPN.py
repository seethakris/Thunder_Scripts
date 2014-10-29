# -*- coding: utf-8 -*-
"""
Plot PCA and ICA components for after and before ablation data of AF4 in the dHb
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
    
    dIPN = [Exp_Name.index(s) for s in Exp_Name if "dIPN" in s]
    vIPN = [Exp_Name.index(s) for s in Exp_Name if "vIPN" in s]

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
        #Plot dIPN
        fig2 = plt.figure()
        for ii in range(0,np.size(dIPN,0)):           
            with sns.axes_style("white"):
                plt.subplot(2,3,ii+1)                        
                plt.imshow(maps[:,:,dIPN[ii],:])
                plt.axis('off')
                plt.title(Exp_Name[dIPN[ii]][15:18],fontsize=10)            
#        plt.tight_layout()   
        fig2 = plt.gcf()
        pp.savefig(fig2)
        plt.close()
                
        #Plot vIPN
        fig2 = plt.figure()
        for ii in range(0,np.size(vIPN,0)):           
            with sns.axes_style("white"):
                plt.subplot(2,3,ii+1)                        
                plt.imshow(maps[:,:,vIPN[ii],:])
                plt.axis('off')
                plt.title(Exp_Name[vIPN[ii]][15:18], fontsize=10)            
#        plt.tight_layout()   
        fig2 = plt.gcf()
        pp.savefig(fig2)
        plt.close()        

    ########### Plot components ##################
    fig2 = plt.figure()
    sns.set_context("talk", font_scale=1.25)
    with sns.axes_style("dark"):
        ax1 = plt.subplot(231)
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
            fig2 = plt.subplot(234)
            sns.tsplot(np.array(matched_signals[ii].clr_grped_signal), ci=95, err_style="ci_band", color=unique_clrs[ii])
            plt.locator_params(axis = 'y', nbins = 4)            
            sns.axlabel("Time (seconds)","a.u")            
            plot_vertical_lines()
        plt.axhline(y=0, linestyle='-', color='k', linewidth=1)

    matching = [Exp_Name.index(s) for s in Exp_Name if "vIPN" in s]
    temp_matched_pixels = matched_pixels[:,matching]
    with sns.axes_style("white"):
        fig2 = plt.subplot(232)
        fig2 = sns.boxplot(np.transpose(temp_matched_pixels),linewidth=2, widths=.5, color=unique_clrs)
        for ii in range(0,np.size(unique_clrs,0)):
            fig2 = plt.plot(np.repeat(ii+1,np.size(temp_matched_pixels,1)), np.transpose(temp_matched_pixels[ii,:]),'s', \
            color=unique_clrs[ii], markersize=5, markeredgecolor='k', markeredgewidth=2) 
            plt.locator_params(axis = 'y', nbins = 4)
        plt.title('vIPN')

#            plt.ylim((0, 1000))
        sns.axlabel("Colors", "Number of Pixels")
        sns.despine(offset=10, trim=True);    
        
    matching = [Exp_Name.index(s) for s in Exp_Name if "dIPN" in s]
    temp_matched_pixels = matched_pixels[:,matching]
    with sns.axes_style("white"):
        fig2 = plt.subplot(233)
        fig2 = sns.boxplot(np.transpose(temp_matched_pixels),linewidth=2, widths=.5, color=unique_clrs)
        for ii in range(0,np.size(unique_clrs,0)):
            fig2 = plt.plot(np.repeat(ii+1,np.size(temp_matched_pixels,1)), np.transpose(temp_matched_pixels[ii,:]),'s', \
            color=unique_clrs[ii], markersize=5, markeredgecolor='k', markeredgewidth=2) 
            plt.locator_params(axis = 'y', nbins = 4)
#            plt.ylim((0, 1000))
        sns.axlabel("Colors", "Number of Pixels")
        plt.title('dIPN')
        sns.despine(offset=10, trim=True);    
        
        
    #Create an lm plot seperating Before and After number of pixels
    #Make Panda data frame    
    A = np.zeros((3,np.size(matched_pixels,0)*np.size(matched_pixels,1)), dtype=np.int)
    count = 0
    for ii in range(0,np.size(matched_pixels,0)):
        for jj in range(0,np.size(matched_pixels,1)):
            A[0,count] = matched_pixels[ii,jj]
            A[1,count] = ii
            if jj in [0,1,3,5]:
                A[2,count] = 1 #vIPN
            else:
                A[2,count] = 0 #dIPN
            count = count+1
    A = np.transpose(A)
    B = pd.DataFrame({'Pixel':A[:,0], 'response':A[:,1], 'IPN':A[:,2]})
    B["Region"] = B.IPN.map({1: "vIPN", 0: "dIPN"})
    



    #Plot mean projection        
    with sns.axes_style("white"):  
        fig2 = plt.subplot(235)
        matching = [Exp_Name.index(s) for s in Exp_Name if "vIPN" in s]
        temp = np.max(maps[:,:,matching,:], axis=2)
        plt.imshow(temp.astype(np.float16).transpose((1,0,2)))
        plt.axis('off')
        plt.title('Max vIPN')
        
        fig2 = plt.subplot(236)
        matching = [Exp_Name.index(s) for s in Exp_Name if "dIPN" in s]
        temp = np.max(maps[:,:,matching,:], axis=2)
        plt.imshow(temp.astype(np.float16).transpose((1,0,2)))
        plt.axis('off')
        plt.title('Max dIPN')
        
        plt.tight_layout()
        fig2 = plt.gcf()
        pp.savefig(fig2)
        plt.close()

    with sns.axes_style("dark"):
         fig3 = plt.figure()
         g = sns.lmplot("Pixel", "IPN", B, y_jitter=.20, hue="response", fit_reg=False, palette=unique_clrs, markers="s", scatter_kws={"s": 50})
         g.set(ylim = (-0.2, 1.2), yticks=[0, 1], yticklabels=["dIPN", "vIPN"])
         plt.axhline(y=0.5, linestyle='-', color='w', linewidth=0.5)
         fig3 = plt.gcf()
         pp.savefig(fig3)
         plt.close()

    with sns.axes_style("dark"):
         fig3 = plt.figure()
         g = sns.lmplot("Pixel", "IPN", B, y_jitter=.20, fit_reg=True, palette=unique_clrs, markers="s", scatter_kws={"s": 50})
         g.set(ylim = (-0.2, 1.2), yticks=[0, 1], yticklabels=["dIPN", "vIPN"])
         plt.axhline(y=0.5, linestyle='-', color='w', linewidth=0.5)
         fig3 = plt.gcf()
         pp.savefig(fig3)
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