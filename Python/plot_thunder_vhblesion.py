# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 11:16:25 2014

Plot PCA and ICA components for after and before ablation data of vhb lesion with raphe and dhb

@author: seetha
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
    
#Pick experiemnt names and seperate to indeces for plotting     
    Fishnum = np.unique([np.int(ii[15:18]) for ii in Exp_Name])

############ Plot Colormaps of scores ############
    
    count = 0
    for ii in range(0,np.size(Fishnum,0)):
        fig2 = plt.figure()
        Num_exp = np.size([jj for jj in Exp_Name if jj[15:18] == str(Fishnum[ii])],0)
        count2 = 5
        count1 = 1
        for jj in range(0, Num_exp):
            with sns.axes_style("white"):                                          
                if "Before" in Exp_Name[count]:
                    plt.subplot(2,4,count1)                        
                    plt.imshow(maps[:,:,count,:].transpose((1,0,2)))
                    plt.axis('off')
                    plt.title(Exp_Name[count][15:-4])
                    count = count+1
                    count1 = count1+1
                else: 
                    plt.subplot(2,4,count2)
                    plt.imshow(maps[:,:,count,:].transpose((1,0,2)))
                    plt.axis('off')
                    plt.title(Exp_Name[count][15:-4])
                    count = count+1
                    count2=count2+1

            
            plt.tight_layout()   
            
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

#Boxplot of number of responses        
    matching = [Exp_Name.index(s) for s in Exp_Name if "Before" in s]
    temp_matched_pixels = matched_pixels[:,matching]
    with sns.axes_style("white"):
        fig2 = plt.subplot(232)
        fig2 = sns.boxplot(np.transpose(temp_matched_pixels),linewidth=2, widths=.5, color=unique_clrs)
        for ii in range(0,np.size(unique_clrs,0)):
            fig2 = plt.plot(np.repeat(ii+1,np.size(temp_matched_pixels,1)), np.transpose(temp_matched_pixels[ii,:]),'s', \
            color=unique_clrs[ii], markersize=5, markeredgecolor='k', markeredgewidth=2) 
            plt.locator_params(axis = 'y', nbins = 4)
            plt.ylim((0, 3000))
        sns.axlabel("Colors", "Number of Pixels")
        sns.despine(offset=10, trim=True);    
        
    matching = [Exp_Name.index(s) for s in Exp_Name if "After" in s]
    temp_matched_pixels = matched_pixels[:,matching]
    with sns.axes_style("white"):
        fig2 = plt.subplot(233)
        fig2 = sns.boxplot(np.transpose(temp_matched_pixels),linewidth=2, widths=.5, color=unique_clrs)
        for ii in range(0,np.size(unique_clrs,0)):
            fig2 = plt.plot(np.repeat(ii+1,np.size(temp_matched_pixels,1)), np.transpose(temp_matched_pixels[ii,:]),'s', \
            color=unique_clrs[ii], markersize=5, markeredgecolor='k', markeredgewidth=2) 
            plt.locator_params(axis = 'y', nbins = 4)
            plt.ylim((0, 3000))
        sns.axlabel("Colors", "Number of Pixels")
        sns.despine(offset=10, trim=True);    
        
#Plot mean projection        
    with sns.axes_style("white"):  
#        fig2 = plt.subplot(438)
#        matching = [Exp_Name.index(s) for s in Exp_Name if "dHbBefore" in s]
#        temp = np.median(maps[:,:,matching,:], axis=2)
#        plt.imshow(temp.astype(np.float16).transpose((1,0,2)))
#        plt.axis('off')
#        plt.title('Max Before')
#        
#        fig2 = plt.subplot(4,3,11)
#        matching = [Exp_Name.index(s) for s in Exp_Name if "dHbAfter" in s]
#        temp = np.median(maps[:,:,matching,:], axis=2)
#        plt.imshow(temp.astype(np.float16).transpose((1,0,2)))
#        plt.axis('off')
#        plt.title('Max After')
#        
        fig2 = plt.subplot(439)
        matching = [Exp_Name.index(s) for s in Exp_Name if "DRBefore" in s]
        temp = np.median(maps[:,:,matching,:], axis=2)
        plt.imshow(temp.astype(np.float16).transpose((1,0,2)))
        plt.axis('off')
        plt.title('Max Before')
        
        fig2 = plt.subplot(4,3,12)
        matching = [Exp_Name.index(s) for s in Exp_Name if "DRAfter" in s]
        temp = np.median(maps[:,:,matching,:], axis=2)
        plt.imshow(temp.astype(np.float16).transpose((1,0,2)))
        plt.axis('off')
        plt.title('Max After')
        
        plt.tight_layout()
        fig2 = plt.gcf()
        pp.savefig(fig2)
        plt.close()


#Create an lm plot seperating Before and After number of pixels
#Make Panda data frame    
    A = np.zeros((4,np.size(matched_pixels,0)*np.size(matched_pixels,1)), dtype=np.int)
    matching_BA= [Exp_Name.index(s) for s in Exp_Name if "Before" in s]
    matching_dhb= [Exp_Name.index(s) for s in Exp_Name if "dHb" in s]
    matching_DR= [Exp_Name.index(s) for s in Exp_Name if "DR" in s]

    count = 0
    for ii in range(0,np.size(matched_pixels,0)):
        for jj in range(0,np.size(matched_pixels,1)):
            A[0,count] = matched_pixels[ii,jj]
            A[1,count] = ii
            if jj in matching_BA:
                A[2,count] = 1 #Before Lesion
            else:
                A[2,count] = 0 #After Lesion
                
            #Check which region
            if jj in matching_dhb:
                A[3,count] = 0
            if jj in matching_DR:
                A[3,count] = 1
                
            count = count+1
    A = np.transpose(A)
    B = pd.DataFrame({'Pixel':A[:,0], 'response':A[:,1], 'Lesion':A[:,2], 'Region':A[:,3]})
    B["BeforeAfter"] = B.Lesion.map({1: "Before", 0: "After"})
    B["Region_Map"] = B.Region.map({0: "dHb", 1:"DorsalRaphe"})
    
    with sns.axes_style("dark"):
         fig3 = plt.figure()
         g = sns.lmplot("Pixel", "Lesion", B, y_jitter=.20, hue="response",col="Region_Map", fit_reg=False, palette=unique_clrs, markers="s", scatter_kws={"s": 50})
         g.set(ylim = (-0.2, 1.2), yticks=[0, 1], yticklabels=["After", "Before"])
         plt.axhline(y=0.5, linestyle='-', color='w', linewidth=0.5)
         fig3 = plt.gcf()
         pp.savefig(fig3)
         plt.close()

    with sns.axes_style("dark"):
         fig3 = plt.figure()
         g = sns.lmplot("Pixel", "Lesion", B,col="Region_Map", y_jitter=.20, fit_reg=True, palette=unique_clrs, markers="s", scatter_kws={"s": 50})
         g.set(ylim = (-0.2, 1.2), yticks=[0, 1], yticklabels=["After", "Before"])
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