# -*- coding: utf-8 -*-
"""
Convert data to text from multi tiffs
"""

#Import relevant libraries
import numpy as np #for numerical operations on arrays
import PIL as pil # for image resizing
import os  

import matplotlib.pyplot as plt #for plotting
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

from libtiff import TIFF #for reading multiTiffs

#Custom libraries
from smooth import smooth

def create_textfile(Exp_Folder, text_filename, num_time, num_time_req, img_size_x, img_size_y, f_f_flag):
    
    filesep = os.path.sep #Get fileseperator according to operating system
    
    zz = 0 #each multitiff file is considered as one stack
    Matfile_for_Thunder = None
    
    #Find tiff files in given experiment folder
    onlyfiles = [ f for f in os.listdir(Exp_Folder) if (os.path.isfile(os.path.join(Exp_Folder,f)) and f.find('.tif')>0)]
        
    #Plotting as pdf
    Figure_PDFDirectory = Exp_Folder+'Figures'+filesep
    if not os.path.exists(Figure_PDFDirectory):
        os.makedirs(Figure_PDFDirectory) 
    pp = PdfPages(Figure_PDFDirectory+text_filename+'_PreprossData.pdf')
    
    #Loop through files found
    for lst in xrange(0,np.size(onlyfiles, axis=0)): 
        tif = TIFF.open(os.path.join(Exp_Folder,onlyfiles[lst]), mode='r') #Open multitiff             
        
        #initialize temp variables    
        zz = zz+1
        data = np.zeros((img_size_x,img_size_y,num_time), dtype=np.uint8)
        ii = 0
        count = 0
    
        #Loop through each image in multitiff. Convert to uint8 if necessary
        for image in tif.iter_images():
            #Check if image is of the xy resolution specified, else resize
            if np.size(image,1)!=img_size_y or np.size(image,0)!=img_size_x:
                if ii == 0:
                    print 'Size mismatch..Resizing '+ onlyfiles[lst]
                temp_image = pil.Image.fromarray(image)    
                #Check if image is uint8, else convert
                if image.dtype == 'uint8':
                    data[:,:,ii] = np.array(temp_image.resize((img_size_y, img_size_x), pil.Image.NEAREST))
                else:
                    data[:,:,ii] = np.uint8(np.array(temp_image.resize((img_size_y, img_size_x), pil.Image.NEAREST))/256)
                ii = ii+1
            else:
                if image.dtype == 'uint8':
                    data[:,:,ii] = image
                else:
                    data[:,:,ii] = np.uint8(image/256)
                ii = ii+1   
                    
       
        #Plot summed data over time for reviewing
        with sns.axes_style("white"):
            fig1 = plt.imshow(np.sum(data[:,:,:], axis=2), cmap='jet')
            plt.title(onlyfiles[lst])
            fig1 = plt.gcf()
            pp.savefig(fig1)
            plt.close()
            
        #Get data in thunder format [xx,yy,zz,time]
        print 'Creating array from stack for ' + onlyfiles[lst]
        temp_matfile_for_thunder = np.zeros([np.size(data, axis=0)*np.size(data, axis=1),3+num_time_req+4], dtype=np.int)
        for yy in xrange(0,np.size(data, axis=1)):
            for xx in xrange(0,np.size(data, axis=0)): 
                    temp_matfile_for_thunder[count,0] = xx+1;
                    temp_matfile_for_thunder[count,1] = yy+1;
                    temp_matfile_for_thunder[count,2] = zz;
                    # Create delta f/f values if necessary
                    if f_f_flag==0:
                        temp_matfile_for_thunder[count,3:] = smooth(data[xx,yy,0:num_time_req],5,'hanning')
                    else:
                        temp_matfile_for_thunder[count,3:] = smooth((data[xx,yy,0:num_time_req])/np.mean(data[xx,yy,5:num_time_req]),5,'hanning')
                    count = count+1 
        
        #Plot data heatmaps
        with sns.axes_style("white"):
            A = temp_matfile_for_thunder[:,3:]            
            fig2 = plt.imshow(A,aspect='auto', cmap='jet')
            plot_vertical_lines()
            plt.title(onlyfiles[lst])
            plt.colorbar()
            fig2 = plt.gcf()
            pp.savefig(fig2)
            plt.close()
            A = None
        
        #Append each tiff files data to a bigger matrix
        if Matfile_for_Thunder == None:
            Matfile_for_Thunder = temp_matfile_for_thunder
        else:
            Matfile_for_Thunder = np.append(Matfile_for_Thunder,temp_matfile_for_thunder, axis=0)
    
    pp.close()        
    print 'Saving all the data in the text file '            
    np.savetxt(Exp_Folder+text_filename+'.txt', Matfile_for_Thunder, fmt='%i')#Save as text file
    return onlyfiles
                
    
def plot_vertical_lines():
    plt.axvline(x=46, linestyle='-', color='k', linewidth=1)
    plt.axvline(x=66, linestyle='--', color='k', linewidth=1)
    plt.axvline(x=86, linestyle='-', color='k', linewidth=1)
    plt.axvline(x=106, linestyle='--', color='k', linewidth=1)
    plt.axvline(x=126, linestyle='-', color='k', linewidth=1)
    plt.axvline(x=146, linestyle='--', color='k', linewidth=1)
    plt.axvline(x=166, linestyle='-', color='k', linewidth=1)
    plt.axvline(x=186, linestyle='--', color='k', linewidth=1)
    
    
        
    
