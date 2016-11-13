# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 13:22:04 2013

@author: erika
@modifier: riccardo
"""
#executes data extraction and no sample fitting

import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys, getopt, os
import re
import ConfigParser

from IPython import embed



def peak_det_hor(spect,ROI_start,ROI_end):
    """
    Search for the number of peaks in horizontal spectrum. 
    If one or two peaks are found computes and adds as decimal, the FWHM of the peak (10pixel of FWHM = 0.1)
    Returns the number of peaks n_p, eventually with the FWHM, the peaks positions p_p, and the peaks heights p_h
    
    Parameters
    ----------
    spect : 1D or 2D array
        Array or 2D array in which each row contains the signal to be analyzed
    thr_value : double default = 50000
        threshold value, discard all the data-points below this value
    leng : int default = 9
        length of the window where search for the peak

    Example
    --------
    >>>import numpy as np
    >>> a=np.asanyarray([[100000,100000,100000,250000,500000,250000,100000,100000,100000],[100000,100000,100000,250000,500000,250000,100000,100000,100000]])
    >>> n_p,p_p,p_h=peak_det_hor(a)
    >>>n_p
    array([1.3,1.3])   
    >>>p_p
    array([4,4])
    >>>p_h
    array([500000,500000])
    """

    n_s=spect.shape[0]
    area=np.zeros(n_s)
    k=0
    for r in np.arange(n_s):
        an_v=spect[r,ROI_start:ROI_end]
	ROI_ext=ROI_end-ROI_start
	#an_back=spect[r,ROI_start-(ROI_ext+1):ROI_start-1]
	#an_back=spect[r,150:200]
	#an_back2=spect[r,ROI_end+1:ROI_end+ROI_ext+1]
	#bkg_val=(an_back.mean()+an_back2.mean())/2
	bkg_val=(spect[r,500:950].mean()) +8500 #(spect[r,500:950].mean())*1.07 -8500
	#an_back2=spect[r,100:200]	
	#area[r]=an_v.sum()/float(ROI_end-ROI_start)-(an_back1.sum()/100.+an_back1.sum()/100.)/2
    	#area[r]= an_v.sum()/float(ROI_end-ROI_start)- an_back.sum()/(499.)
	area[r]= an_v.sum()/float(ROI_end-ROI_start)- bkg_val
	# x=


    return area




def normalize(x,y):
    """
    The same of numpy.divide
    Returns x/y
    """
    return np.divide(x,y)



def scatplt(x,y,labx,laby,col,size):
    """
    Creates a figure and plots a scatterplot
    
    Parameters
    ----------
    x : 1D array of double
        X-axys values
    y : 1D array of double
        Y-axys values
    labx : string
        label of the x-axis
    laby : string
        label of the y-axis
    col : char
        color of the points
    s : int
        size of the points
    """
    fig1=plt.figure()
    plt.xlabel(labx)
    plt.ylabel(laby) 
    plt.scatter (x,y,color=col,s=size)
#    plt.show()
    

def fitter(i0,i1,deg):
    """
    Fits i0 vs i1 using a polynomium of degree deg in a runtime defined range of i1.\n
    Returns an array containing the coefficients,(first coefficient is the higest degree monomium's one)  and the bounds of i1 interval inserted by the user.\n
    See http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html for more references.
    
    Parameters
    ----------
    i0 : 1D array of double
        array
    i1 : 1D array of double
        array    
    deg : int
        degree of the polynomium

    Examples
    --------
    >>> p1,i1_min,i1_max = fitter(i0,i1,3)
    """
    scatplt(i0,i1,'i0','i1','black',2)

#    i1_min=float(raw_input('i1_min value: '))
#    i1_max=float(raw_input('i1_max value: '))
    i0_min=float(raw_input('i0_min value: '))
    i0_max=float(raw_input('i0_max value: '))

#       cutting on i1
#    i0_b=i0[np.where(i1[np.where(i1>=i1_min)]<=i1_max)]    
#    i1_b=i1[np.where(i1[np.where(i1>=i1_min)]<=i1_max)]
#    p1 =np.polyfit(i1_b,i0_b,deg) 
    i0_b=i0[np.where(i0[np.where(i0>=i0_min)]<=i0_max)]    
    i1_b=i1[np.where(i0[np.where(i0>=i0_min)]<=i0_max)]
    p1 =np.polyfit(i0_b,i1_b,deg) 

    p2=np.zeros(deg+1)
    for i in range(deg+1): #need to revert the array 'cause polyfit and polynomial.polieval has the opposite convention for coefficients
        p2[i]=p1[deg-i]
    p1=p2    
    
    return p1,i0_min,i0_max    
#    return p1,i1_min,i1_max





def my_filter_multiple(fil_fin,iom_sh_a,caen_area,ov,i0,i1,ov_fin,n_var,good): 
    """
    Returns, in multiple shot mode, the values of all imported variables that satisfies to the filters conditions, eventually appended to previously filtered data.
    
    Parameters
    ----------
    fil_fin : 1D array of 0 and 1
        0 represents values to be rejected, 1 represents values to be accepted
    iom_sh_a : 1D array of double
        array of values representing the X variable    
    caen_area : 1D array of double
        array of values representing the Y variable
    ov : 2D array of double
        array cointaining all the other imported variables
    i0 : 1D array where the filtered i0 is stored
        each time the function is called append to i0 the corresponding filtered value of iom_sh_a starting from position good
    i1 : 1D array where the filtered i1 is stored
        each time the function is called append to i1 the corresponding filtered value of caen_area starting from position good
    ov_fin : 2D array where the filtered ov is stored
        each time the function is called append to ov_fin the corresponding values (row) of ov starting from position good
    n_var : int
        is the number of extra variables imported, and thus the number of columns in ov ov_fin 2D array
    good : int
        starting point for the appending

    Examples
    --------
    >>> i0,i1,ov_fin,good=my_filter_multiple(fil_fin,iom_sh_a,caen_area,ov,i0,i1,ov_fin,n_var,good)
    """

    for index in np.arange(fil_fin.shape[0],dtype='uint16'):
        if fil_fin[index]==1 :
            i0[good]=iom_sh_a[index]
            i1[good]=caen_area[index]
            if n_var>0:
                ov_fin[good,:]=ov[index,:]
            good+=1    
   
    
    return i0,i1,ov_fin,good



def my_filter_single(fil_fin,iom_sh_a,caen_area,ov,k,i0,i1,ov_fin,n_var,good): #where???
    """
    Returns, in single shot mode, the values of all imported variables that are firts shots and satisfies to the filters conditions, eventually appended to previously filtered data.
    
    Parameters
    ----------
    fil_fin : 1D array of 0 and 1
        0 represents values to be rejected, 1 represents values to be accepted
    iom_sh_a : 1D array of double
        array of values representing the X variable    
    caen_area : 1D array of double
        array of values representing the Y variable
    ov : 2D array of double
        array cointaining all the other imported variables
    k : 1D array of double
        array containing boundary condition to be satisfied for each variable
    i0 : 1D array where the filtered i0 is stored
        each time the function is called append to i0 the corresponding filtered value of iom_sh_a starting from position good
    i1 : 1D array where the filtered i1 is stored
        each time the function is called append to i1 the corresponding filtered value of caen_area starting from position good
    ov_fin : 2D array where the filtered ov is stored
        each time the function is called append to ov_fin the corresponding values (row) of ov starting from position good
    n_var : int
        is the number of extra variables imported, and thus the number of columns in ov ov_fin 2D array
    good : int
        starting point for the appending

    Examples
    --------
    >>> i0,i1,ov_fin,good=my_filter_single(fil_fin,iom_sh_a,caen_area,ov,k,i0,i1,ov_fin,n_var,good)
    """
    for index in np.arange(fil_fin.shape[0],dtype='uint16'):
        if ((iom_sh_a[index] >= k[1]) and (iom_sh_a[index]) <= k[2]) :
            if fil_fin[index]==1:
                i0[good]=iom_sh_a[index]
                i1[good]=caen_area[index]
                if n_var>0:
                    ov_fin[good,:]=ov[index,:]
                good+=1
                break 
            else:
                break
        
    return i0,i1,ov_fin,good



def build_fil(iom_sh_a,k,ov,n_var,total_fil): 
    """
    Returns the filter matrix for the current file (fil) the total filter matrix (total_fil) for the whole examined files and the fil_fin array.\n
    fil_fin contains 1 (the corresponding row of data has been validated) and 0 (the corresponding row of data has been discarded).
    
    Parameters
    ----------
    iom_sh_a : 1D array
        array of i0 value to be filtered
    k : 1D array of double
        array containing boundary condition to be satisfied for each variable    
    ov : 2D array of double
        array cointaining all the other imported variables
    n_var : int
        is the number of extra variables imported, and thus the number of columns in ov ov_fin 2D array
    total_fil : 2D array of 0 and 1 
        contains the filter matrix representative of all the imported data 

    Examples
    --------
    >>> fil,total_fil,fil_fin = build_fil(iom_sh_a,k,ov,n_var,total_fil)
    """
    fil=np.zeros((iom_sh_a.shape[0],n_var+1))
# Fill the first column of fil with 1 if the condition is satisfied and with 0 otherwise
    fil[:,0] = np.in1d(iom_sh_a,iom_sh_a[np.where(iom_sh_a<=k[2])][np.where(iom_sh_a[np.where(iom_sh_a<=k[2])]>=k[1])])
# doing the same for the other variables 
    for i in np.arange(1,n_var+1,dtype='uint8'):
        fil[:,i] = np.in1d(ov[:,i-1],ov[:,i-1][np.where(ov[:,i-1]<=k[2*i+2])][np.where(ov[:,i-1][np.where(ov[:,i-1]<=k[2*i+2])]>=k[2*i+1])])
    total_fil=np.append(total_fil,fil,axis=0)    
    fil_fin=np.zeros(iom_sh_a.shape[0])
    for j in np.arange(iom_sh_a.shape[0],dtype='uint16'):  
        if fil[j,:].sum()==n_var+1:
            fil_fin[j]=1
    

    return fil,total_fil,fil_fin




def load_var(file_name, variables,total_i0,total_i1,total_ov=0,n_var=0):
    """
    Load the specified variables from file_name, and perform a first peak analisys if the hor_spectrum or vert_spectrum is imported.    
    The variables are appended to the corresponding input array 
    Parameters
    ----------
    file_name : string 
        name of the file to open
    variables : list of string in pair
        in each pair the first value is the full h5 address of the variable to import, the second is the bkg to subtract   
    total_i0 : 1D array of double
        array cointaining all the i0 values imported from all files
    total_i1 : 1D array of double
        array cointaining all the i1 values imported from all files
    total_ov : 2D array of double
        2D array cointaining all the other values imported from all files
    n_var : int
        is the number of extra variables imported, and thus the number of columns in ov ov_fin 2D array

    Examples
    --------
    >>> total_i0,total_i1,len_iom, iom_sh_a,caen_area,ov,total_ov = load_var(file_name, variables,total_i0,total_i1,total_ov,n_var)
    """

    input_file = h5py.File(file_name,'r')
    iom_sh_a  = np.asanyarray(input_file[variables[0]])-int(variables[1])
    caen_area = np.asanyarray(input_file[variables[2]])-int(variables[3])
    total_i0=np.append(total_i0,iom_sh_a)
    total_i1=np.append(total_i1,caen_area)
    len_iom=iom_sh_a.shape[0]
#dinamically imports all the extra variables and putting into a matrix

    if n_var>0:
        ov=np.zeros((iom_sh_a.shape[0],n_var)) #viene sovrascritto ad ogni giro mentre a total_ov viene appeso        
        for i in np.arange(n_var,dtype='uint8'):
            if variables[2*i+4]=='photon_diagnostics/Spectrometer/hor_spectrum':
                d=np.asanyarray(input_file[variables[2*i+4]])-float(variables[2*i+5])
		# change when changing spectrometer roi
                ov1=peak_det_hor(d,250,400)
                #ov1=peak_det_hor(d,380,450)
		#print ov1
            else:
                ov1=np.asanyarray(input_file[variables[2*i+4]])-float(variables[2*i+5])
		ov1=ov1.ravel()
		#embed()
		#ov1.reshape((1000))
		#print ov1.shape
            ov[:,i]=ov1
        total_ov=np.append(total_ov,ov,axis=0)
    input_file.close()  
    if n_var>0:
        return total_i0,total_i1,len_iom, iom_sh_a,caen_area,ov,total_ov
    else:
        return total_i0,total_i1,len_iom, iom_sh_a,caen_area


def search(input_path):
    """
    Search for hdf in folders and subfolders
    Create a list of all subfolders and files.
    
    Parameters
    ----------
    input_path : str
        Name of main folder

    Examples
    --------
    >>> fold_list,file_list = search('inputh_path')
    """
    file_list=[]
#    file_list=np.asarray(0)
    folder_list=np.asarray(glob.glob(input_path+'/'+'*'))
    print '\n'    
    subfolders=1    
    for i in range(folder_list.shape[0]):   
        if glob.glob(folder_list[i]+'/'+'*'):
            file_list.extend(glob.glob(folder_list[i]+'/'+'*'))
         
    if not file_list:
        subfolders=0 
        file_list=folder_list

    if subfolders:
        print 'list of sub-folders:'
        print folder_list, folder_list.shape[0], '\n'
 
    print 'list of files:'
    print file_list, len(file_list), '\n'
   
    return subfolders,file_list




def main(argv):
# read configuration parameters and sets the variables   
    length=100000
    corr_file=[]
    plot0=1
    plot1=0
    fitting=0
    multiple=0
    filtered_data=0
    plt.ion() #don't wait to close the plot to go on
    print argv[0]        
    config=ConfigParser.RawConfigParser()
    config.read('timex_conf.cfg')
    sect='Section%s' %argv[0]
    print sect
    arg1 = config.get(sect,'input')
    arg2 = config.get(sect,'output')
    variables=config.get(sect,'variables').split(',')
    n_var=len(variables)/2-2 
    k=map(float,config.get(sect,'mod').split(','))
    multiple=k[0]
    if config.get(sect,'plotting')=='1':
        filtered_data=1
    if config.has_option(sect,'fit'):
        fitting=1
        deg=config.getint(sect,'fit')
        
    
    
    #print filtered_data
#    search()
    input_path =arg1
    n=re.split("/",input_path)
    name=n[-3]
    output_path =arg2
    try:
        os.makedirs(output_path)
    except:
        os.stat(output_path)
#search for all subdir and files
    subfolders,file_list=search(input_path)

    print 'additional variables: ', n_var, '\n'
#    print 'variable:', variables[4],variables[5], '\n'
        
    i0=np.zeros(length)
    i1=np.zeros(length)
    ov_fin=np.zeros((length,n_var))
#initialing the variables to fill during import    

    num=0
    total_i0=np.zeros(0)
    total_i1=np.zeros(0)
    total_ov=np.zeros((0,n_var))
    total_fil=np.zeros((0,n_var+1),dtype='uint8')
    good=0    

    black_list=['/net/online4eis/store/eis-timex/inhouse_Jul2014/Al_FEL-self-transm/17.2nm_S-foc_01/rawdata/17.2nm_S-foc_01_acq_500_shots_64297329.h5']
    for file in range(len(file_list)):
        file_name='%s' %file_list[file]
                 
#searches for variables to import                           
        try:
            if file_name not in black_list:
                total_i0,total_i1,len_iom, iom_sh_a,caen_area,ov,total_ov = load_var(file_name, variables,total_i0,total_i1,total_ov,n_var)  
       
        except IOError, ex:
            corr_file.append(file_name)
	    print ex
            continue
        except KeyError, ex:
            corr_file.append(file_name)
	    print ex
            continue
#populating the filter matrix

#creating the filter matrix and setting the value 1 for each element to be kept ,total_fil is the matrix with the filters for all the hdf file in dataset
        fil,total_fil,fil_fin = build_fil(iom_sh_a,k,ov,n_var,total_fil)
         


#filtering the data


# for first shot mode first check if filter on i0 is passed then applyes other filters 
        if multiple==0: 
            i0,i1,ov_fin,good=my_filter_single(fil_fin,iom_sh_a,caen_area,ov,k,i0,i1,ov_fin,n_var,good)
        else:
           i0,i1,ov_fin,good=my_filter_multiple(fil_fin,iom_sh_a,caen_area,ov,i0,i1,ov_fin,n_var,good) 
 

    i0=i0[0:good]
    i1=i1[0:good]
    if n_var>0:
        ov_fin=ov_fin[0:good,:]
    number_i=good
    good=0

  
    print 'total amount of good shots:', number_i, '\n'
    print 'list of corrupted files:'
    print corr_file, '\n'
   # fitting only if required
    if fitting:
	print 'deg',deg
        p1,i0_min,i0_max = fitter(i0,i1,deg)
#        p1,i1_min,i1_max = fitter(i1,i0,deg)
        plot1=1 #enables fit plotting
    plt.figure(1)
    plt.subplot(2,1,1)
    if plot0:

        filt=plt.scatter(i0, i1/i0, color='black', s=2)        
        if filtered_data:

            raw=plt.scatter(total_i0, total_i1/total_i0, color='gray', s=2) 
            filt=plt.scatter(i0, i1/i0, color='black', s=2)                        
            plt.legend((raw,filt),('unfiltered','filtered'),scatterpoints=1, loc='upper right',ncol=1, prop={'size':10}) #, fontsize=8     
        plt.xlabel('i0')
        plt.ylabel('i1 / i0')
  

    plt.subplot(2,1,2)
    plt.scatter(i0, i1,color='black', s=2)

    if filtered_data:

        raw=plt.scatter(total_i0, total_i1, color='gray', s=2) 
        filt=plt.scatter(i0, i1, color='black', s=2)                        
        plt.legend((raw,filt),('unfiltered','filtered'),scatterpoints=1, loc='upper right',ncol=1, prop={'size':10}) #, fontsize=8   

    plt.xlabel('i0')
    plt.ylabel('i1')  

    if plot1:
        l=np.linspace(i0_min,i0_max,200)
        plt.plot( l,np.polynomial.polynomial.polyval(l,p1),color='red',linewidth=2)
#        l=np.linspace(i1_min,i1_max,200)
#        plt.plot( np.polynomial.polynomial.polyval(l,p1),l,color='yellow')


    output_txt = '%s%s_F.txt' %(output_path,name)  #address of the filtered data
    output_par ='%s%s_F_par.txt' %(output_path,name)  #address of the fit parameters 
    output_UF_txt ='%s%s_UF.txt' %(output_path,name)  #address of the unfiltered data
    output_fil ='%s%s_filters.txt' %(output_path,name)  #address of the whole filters matrix



    plt.savefig('%s%s.png' %(output_path,name))

    if fitting:
        print 'output parameter file: ', output_par
#        output_txt = output_path +name+'_F.txt'

        
        try:     

            coefficienti=np.array(p1)
            degree=np.array([deg])
            tagli=np.array([i0_min,i0_max])
#            tagli=np.array([i1_min,i1_max])
            parametri=np.concatenate((degree,coefficienti,tagli))
            np.savetxt(output_par, parametri)
        except IOError, ex:
            print 'Exception',ex
    
    try:     
        
        if n_var>0:
            output=np.column_stack((i0,i1,ov_fin))
        else:
            output=np.column_stack((i0,i1))
        np.savetxt(output_txt, output,delimiter='	')

        if n_var>0:

            output_UF=np.column_stack((total_i0,total_i1,total_ov))
        else:
            output_UF=np.column_stack((total_i0,total_i1))
        np.savetxt(output_UF_txt, output_UF,delimiter='	')
        np.savetxt(output_fil,total_fil,delimiter='  ')
    except IOError, ex:
        print 'Exception',ex 
    print 'output txt file: ', output_txt, '\n'
    raw_input('Done')



if __name__=='__main__':
    main(sys.argv[1:]) 
