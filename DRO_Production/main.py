#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:55:12 2021

@author: Eve Shalom
"""
#%% Import all relevent functions

from matplotlib import pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from ROIselection import ROISelector, ICfromROI,ROILasso,maskto4d
from TKModelFunctions import modifiedToftsMurase, modifiedToftsMurase1Vox, ForwardsModTofts
from DICOM_processing import dicomdirscan, Slice_Time_Readings, SpatiotempSorter, SignalEnhancementExtract
from DICOM_processing import produceDICOM,produce_vfa_DICOM, mkdir_p, produceDICOMktrans
from plots_and_filters import checkROIs,MedianFilter, animateSingle, displaySlices
from Conc2DROSignal import createR10_withref, calcR1_R2
from Conc2DROSignal import Conc2Sig_Linear
from Conc2DROSignal import STDmap, addnoise
import pydicom

#%% DRO directory
# IMPORTANT YOU MUST EXTRACT THE ZIP FILES IN DIRECTORY 'AIFmasksforDRO' TO RUN IN CURRENT STATE #
DROnum = 2 #Choose 1 or 2
uidprefix = '1.3.6.1.4.1.9328.50.16.'

if DROnum == 1:
    dro_IDnum = '7868737135'
    dirPath = 'original_RIDER_Patients/RIDER Neuro MRI-1023805636/09-01-1904-BRAINRESEARCH-87950/13.000000-perfusion-23726/'
    dir_prot = 'original_RIDER_Patients/RIDER Neuro MRI-1023805636'
    aif_dir = 'AIFmasksforDRO/Synthetic_P1/'
    
if DROnum == 2:
    dro_IDnum = '9215224289'
    dirPath = 'original_RIDER_Patients/RIDER Neuro MRI-3369019796/03-21-1904-BRAINRESEARCH-71264/12.000000-perfusion-17557/
    dir_prot = 'original_RIDER_Patients/RIDER Neuro MRI-3369019796'
    aif_dir = 'AIFmasksforDRO/Synthetic_P2/'

Study1UID = pydicom.uid.generate_uid(prefix=uidprefix)
Study2UID = pydicom.uid.generate_uid(prefix=uidprefix)
StudyUID = [Study1UID,Study2UID]

save_new_DRO_data = False

if save_new_DRO_data == True:
    
    my_new_label = '' #Add new label to name DRO with any new additions
    dro_dir = 'newOutputs/DRO_Patient_{}_{}'.format(dro_IDnum,my_new_label)
    array_dir = dro_dir+'/np_param_arrays/'
    dicom_Ktrans = dro_dir+'/ktransDICOM/'
    dicom_save_path = dro_dir+'/paramDICOM/'
    signal_save_path = dro_dir+'/signalDICOM/'
    
    mkdir_p(dicom_Ktrans)
    mkdir_p(dicom_save_path)
    mkdir_p(array_dir)

    mkdir_p(array_dir+'Visit1/')
    mkdir_p(array_dir+'Visit2/')
    mkdir_p(dicom_save_path+'Visit1/')
    mkdir_p(dicom_save_path+'Visit2/')
    mkdir_p(signal_save_path+'Visit1/')
    mkdir_p(signal_save_path+'Visit2/')

#%% Read in DICOM and extract enhancement

# Specify directory path where .dcm files are stored
dicomRawList,dicomSliceLoc, dicomAcqTime,Refdcm = dicomdirscan(dirPath) # reads in all dicom files in folder into a list for sorting
Slices, Acqs,SliceSet = Slice_Time_Readings(dicomSliceLoc, dicomAcqTime) # Creat slice/aquasition times list and creates dictionary of correct size
Slices = SpatiotempSorter(dicomRawList,Slices,SliceSet) # sorts dictionary of acquisitions by slice and time

# gets data shape for the 4D data
datashape = [Refdcm.Rows,Refdcm.Columns,len(SliceSet),len(Acqs)]

E, S0, S = SignalEnhancementExtract(Slices,datashape,5) # returns signal enhancement (E), baseline signal (S0), and the raw signal (S)

Max = np.max(E,axis=3) # Max signal enhancement of all times for each voxel
auc = np.sum(E,axis=3) # Area under curve (SIC) for each voxel

dt = 4.8/60 # mins from the RIDER website DCE description
t = np.linspace(0,dt*S.shape[-1],S.shape[-1]) # time series points
#%%
slices_3d =[] #create dicom template for Ktrans
for i in range(0,16):
    slices_3d += [Slices['{}'.format(i+1)][0]]#use first time point in each slice
    
#%% AIF lasso slection - aim to target Middle Cerebral Arteries
# Plotting must be set to automatic in spyder or notebook in jupyter
#If not the graph will not open interactively  if choose to draw ROIs

use_saved = 1 # 1 if using pre saved defined regions; 0 if want to draw ROI by hand
z = 5
if use_saved == 0:
	aifmask,roivoxa, aiflasso = ROISelector(Max,z, datashape,False,dro_dir)
	#%% Select sagittal sinus voxels to form Ev
	# Plotting must be set to automatic in spyder or notebook in jupyter
	# If not the graph will not open interactively  if choose to draw ROIs
	sagmask,roivoxv, saglasso = ROISelector(Max,z, datashape,False,dro_dir)

	# Select tumour region
	#tummask, roivoxt, tumlasso = ROISelector(Max,6,datashape,False,dro_dir) # Only use for DRO 2 

#%% Extract MCAs and SSs enhancement time curves
#If you wish to use saved masks input use_saved = 1
#int(input('Do you want to use saved masks? (0:no, 1:yes)\n'))

if use_saved == 1:
    aifmask = np.load('{}aifmask.npy'.format(aif_dir))
    sagmask = np.load('{}sagmask.npy'.format(aif_dir))
    roivoxa = np.load('{}aifmaskvox.npy'.format(aif_dir))
    roivoxv = np.load('{}sagmaskvox.npy'.format(aif_dir))


Ea = ICfromROI(E, aifmask, roivoxa, aifmask.ndim-1)
Ev = ICfromROI(E, sagmask, roivoxv, sagmask.ndim-1)
S0ref = ICfromROI(S0[:,:,z], sagmask[:,:,z,0], roivoxv, sagmask.ndim-2)

checkROIs(Ea,Ev) #uncomment to check ROI enhancement curves (Ea and Ev)
# Voxel wise modified Tofts
# K1, k2, Vp all 3D arrays over the whole spatial domain
# K1 = Ktrans, k2 = kep, Vp = Vp

#%% Fit modified Tofts
Hct = 0.45 # Hematocrit standard value (Sourbron & Buckley, 2013)

K1, k2, Vp = modifiedToftsMurase(Ea/(1-Hct), E, dt, datashape)

# Partial Volume correction
pvc_K1, pvc_k2, pvc_Vp = modifiedToftsMurase1Vox(Ea/(1-Hct), Ev, dt, datashape)
pvc = abs((1-Hct)/pvc_Vp) # sagittal sinus Vp should be 0.55, so calc correction factor if not

# Apply correction factor to fitted parameters
cor_K1 = K1*pvc
cor_k2 = k2*pvc
cor_Vp = Vp*pvc 
# Apply Median Filter to parameters all with footprint (3,3)

mf_K1 = MedianFilter(cor_K1)
mf_k2 = MedianFilter(cor_k2)
mf_Vp = MedianFilter(cor_Vp)

# Sanity check parameter values

# No negative rates or volumes 
# volume max of 1

fmod = 1 # forward model chooser. 1: Modified Tofts, 2: 2CXM

mf_Vp[mf_Vp<=0]=0
mf_Vp[mf_Vp>1]=1.0
mf_K1[mf_K1<=0]=0
mf_k2[mf_k2<=0]=0

# evolve forwards model
aif_cor = np.max((Ea/(1-Hct))*pvc)/6 # caluclate enhancement to concentration correction
Cp = ((Ea/(1-Hct))*pvc)/aif_cor # calculate Cp using signal enhancement aif and concentration conversion factor

 
Ctiss=ForwardsModTofts(mf_K1,mf_k2, mf_Vp, Cp, dt) # run forwards model (exponential convolution approach)

# State/find MR parameters

r1      = 3.9 #longitudinal relaxivity Gd-DTPA (Hz/mM) source: (Pintaske,2006)
r2st    = 10 # transverse relaxivity Gd-DTPA (Hz/mM) roughly estimated using (Pintaske,2006) and (Siemonsen, 2008)
fa      = np.deg2rad(float(Refdcm.FlipAngle)) #flip angle (rads)

Te      = 1e-3*float(Refdcm.EchoTime) # Echo Time (s)
Tr      = 1e-3*float(Refdcm.RepetitionTime) # Repetition Time (s)
T1b     = 1.48 # T1 for blood measured in sagittal sinus @ 1.5T (s) (Zhang et al., 2013)

R10_value   = 1.18 #precontrast T1 relaxation rate (Hz) brain avg (radiopedia)
R20st_value = 17.24 #precontrast T1 relaxation rate (Hz) brain avg using T2* from (Siemonsen, 2008)
R10         = createR10_withref(S0ref,S0,Tr,fa,T1b,datashape) # precontrast R1 map (normalised to sagittal sinus)
#
sigConv = 1 # Sigconv: 1 = linear conc2sig, 2: water excahnge conc2sig

if sigConv == 1:
    
    R1, R2st    = calcR1_R2(R10,R20st_value,r1,r2st,Ctiss) # returns R10 and R2st maps
    dro_S, M    = Conc2Sig_Linear(Tr,Te,fa,R1,R2st,S,1,0) # return dro Signal scaled to range of original data set

# Add nosie
stdS = STDmap(S, t0=5) # caluclate Standard deviation for original data
dro_Snoise = addnoise(1, stdS, dro_S,Nt=datashape[-1]) # apply rician noise at SNR of original data

##% Adjust map parameters for visit 2

trans_K1=mf_K1.copy()
trans_k2=mf_k2.copy()
trans_Vp=mf_Vp.copy()

if DROnum == 1:
    vmax_K1, vmax_k2, vmax_Vp = 1, 1, 0.2
    vmin_K1,vmin_k2, vmin_Vp= 0.2, 0.2, 0.01
    lb_K1, lb_k2, lb_Vp = 0.5, 0.48, 0.52
    ub_K1, ub_k2, ub_Vp = 1.5, 1.55, 1.4
    lim_K1, lim_k2, lim_Vp = vmax_K1+0.5,vmax_k2+0.1, vmax_Vp+0.5
    ub_lim = 1.01
    
if DROnum == 2:
    vmax_K1, vmax_k2, vmax_Vp = 1, 1, 0.2
    vmin_K1,vmin_k2, vmin_Vp= 0.2, 0.2, 0.01
    lb_K1, lb_k2, lb_Vp = 0.54, 0.52, 0.49
    ub_K1, ub_k2, ub_Vp = 1.52, 1.5, 1.43
    lim_K1, lim_k2, lim_Vp = vmax_K1+0.5,vmax_k2+0.1, vmax_Vp+0.5
    ub_lim = 1.01
    

trans_K1[trans_K1 <= vmin_K1] = trans_K1[trans_K1 <= vmin_K1]*lb_K1
trans_K1[trans_K1 >= lim_K1]=trans_K1[trans_K1 >= lim_K1]*ub_lim
trans_K1[(trans_K1 >= vmax_K1) & (trans_K1 < lim_K1)] = trans_K1[(trans_K1 >= vmax_K1) & (trans_K1 < lim_K1)]*ub_K1
trans_K1[(trans_K1 > vmin_K1) & (trans_K1 < vmax_K1)]= trans_K1[(trans_K1 > vmin_K1) & (trans_K1 < vmax_K1)]*(lb_K1 + ( ( (trans_K1[(trans_K1 > vmin_K1) & (trans_K1 < vmax_K1)] - vmin_K1) / (vmax_K1-vmin_K1) ) *(ub_K1-lb_K1) ))

trans_k2[trans_k2 <= vmin_k2] = trans_k2[trans_k2 <= vmin_k2]*lb_k2
trans_k2[trans_k2 >= lim_k2]=trans_k2[trans_k2 >= lim_k2]*ub_lim
trans_k2[(trans_k2 >= vmax_k2) & (trans_k2 < lim_k2)] = trans_k2[(trans_k2 >= vmax_k2) & (trans_k2 < lim_k2)]*ub_k2
trans_k2[(trans_k2 > vmin_k2) & (trans_k2 < vmax_k2)]= trans_k2[(trans_k2 > vmin_k2) & (trans_k2 < vmax_k2)]*(lb_k2 + ( ( (trans_k2[(trans_k2 > vmin_k2) & (trans_k2 < vmax_k2)] - vmin_k2) / (vmax_k2-vmin_k2) ) *(ub_k2-lb_k2) ))

trans_Vp[trans_Vp <= vmin_Vp] = trans_Vp[trans_Vp <= vmin_Vp]*lb_Vp
trans_Vp[(trans_Vp >= lim_Vp)] = trans_Vp[trans_Vp >= lim_Vp]*ub_lim
trans_Vp[(trans_Vp >= vmax_Vp) & (trans_Vp < lim_Vp)] = trans_Vp[(trans_Vp >= vmax_Vp) & (trans_Vp < lim_Vp)]*ub_Vp
trans_Vp[(trans_Vp > vmin_Vp) & (trans_Vp < vmax_Vp)]= trans_Vp[(trans_Vp > vmin_Vp) & (trans_Vp < vmax_Vp)]*(lb_Vp + ( ( (trans_Vp[(trans_Vp > vmin_Vp) & (trans_Vp < vmax_Vp)] - vmin_Vp) / (vmax_Vp-vmin_Vp) ) *(ub_Vp-lb_Vp) ))

trans_Vp[trans_Vp>1]=1

Ctiss_tr=ForwardsModTofts(trans_K1,trans_k2, trans_Vp, Cp, dt)

sigConv = 1 # Sigconv: 1 = linear conc2sig

if sigConv == 1:
    
    R1_tr, R2st_tr = calcR1_R2(R10,R20st_value,r1,r2st,Ctiss_tr)
    dro_S_tr, M_tr = Conc2Sig_Linear(Tr,Te,fa,R1_tr,R2st_tr,S,1,M)
    
dro_Snoise_tr = addnoise(1, stdS, dro_S_tr,Nt=datashape[-1])
#%%

if save_new_DRO_data == True:
    # Produce visit 1 DICOM
    produceDICOM('{}'.format(signal_save_path),Slices,dro_Snoise,1,dro_IDnum,StudyUID[0])

    # Produce visit 2 DICOM
    produceDICOM('{}'.format(signal_save_path),Slices,dro_Snoise_tr,2,dro_IDnum,StudyUID[1])
    
    # Save DICOM Ktrans maps (*1000 as DICOM integer) parameter maps
    produceDICOMktrans('{}'.format(dicom_Ktrans),slices_3d,mf_K1,1,dro_IDnum,StudyUID[0])
    produceDICOMktrans('{}'.format(dicom_Ktrans),slices_3d,trans_K1,2,dro_IDnum,StudyUID[1])

    #Save DICOM of R10, kep, vp, Ctiss, R1, R2st (also *1000)

    produceDICOMktrans('{}/{}'.format(dicom_save_path,'R10'),slices_3d,R10,1,dro_IDnum,StudyUID[0])

    produceDICOMktrans('{}/{}/{}'.format(dicom_save_path,'Visit1','kep'),slices_3d,mf_k2,1,dro_IDnum,StudyUID[0])
    produceDICOMktrans('{}/{}/{}'.format(dicom_save_path,'Visit1','vp'),slices_3d,mf_Vp,1,dro_IDnum,StudyUID[0])

    produceDICOMktrans('{}/{}/{}'.format(dicom_save_path,'Visit2','kep'),slices_3d,trans_k2,2,dro_IDnum,StudyUID[1])
    produceDICOMktrans('{}/{}/{}'.format(dicom_save_path,'Visit2','vp'),slices_3d,trans_Vp,2,dro_IDnum,StudyUID[1])

    produceDICOM('{}/{}/{}'.format(dicom_save_path,'Visit1','Ctiss'),Slices,Ctiss,1,dro_IDnum,StudyUID[0])
    produceDICOM('{}/{}/{}'.format(dicom_save_path,'Visit2','Ctiss'),Slices,Ctiss_tr,2,dro_IDnum,StudyUID[1])

    produceDICOM('{}/{}/{}'.format(dicom_save_path,'Visit1','R1'),Slices,R1,1,dro_IDnum,StudyUID[0])
    produceDICOM('{}/{}/{}'.format(dicom_save_path,'Visit1','R2st'),Slices,R2st,1,dro_IDnum,StudyUID[0])

    produceDICOM('{}/{}/{}'.format(dicom_save_path,'Visit2','R1'),Slices,R1_tr,2,dro_IDnum,StudyUID[1])
    produceDICOM('{}/{}/{}'.format(dicom_save_path,'Visit2','R2st'),Slices,R2st_tr,2,dro_IDnum,StudyUID[1])

    # Save np parameter maps visit 1

    np.save('{}/{}/{}'.format(array_dir,'Visit1','R10.npy'),R10,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','Ea.npy'),Ea,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','Ev.npy'),Ev,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','aifmask.npy'),aifmask,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','sagmask.npy'),sagmask,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','aifmaskvox.npy'),roivoxa,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','sagmaskvox.npy'),roivoxv,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','M.npy'),M,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','R1.npy'),R1,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','R2st.npy'),R2st,allow_pickle=False)

    np.save('{}/{}/{}'.format(array_dir,'Visit1','K1.npy'),mf_K1,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','k2.npy'),mf_k2,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','Vp.npy'),mf_Vp,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','AIF.npy'),Cp,allow_pickle=False)

    np.save('{}/{}/{}'.format(array_dir,'Visit1','stdS.npy'),stdS,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','dro_S.npy'),dro_S,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','dro_Snoise.npy'),dro_Snoise,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit1','Ctiss.npy'),Ctiss,allow_pickle=False)

    np.save('{}/{}/{}'.format(array_dir,'Visit2','M_tr.npy'),M_tr,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit2','R1_tr.npy'),R1_tr,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit2','R2st_tr.npy'),R2st_tr,allow_pickle=False)

    np.save('{}/{}/{}'.format(array_dir,'Visit2','K1.npy'),trans_K1,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit2','k2.npy'),trans_k2,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit2','Vp.npy'),trans_Vp,allow_pickle=False)

    np.save('{}/{}/{}'.format(array_dir,'Visit2','dro_S.npy'),dro_S_tr,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit2','dro_Snoise.npy'),dro_Snoise_tr,allow_pickle=False)
    np.save('{}/{}/{}'.format(array_dir,'Visit2','Ctiss.npy'),Ctiss_tr,allow_pickle=False)


    # creating variable flip angle data and saving
    import os
    exclude_directories = set(['Visit2'])
    dir_vfas = []
    sdir_vfas = []
    for dname, dirs, files in os.walk(dir_prot):  #loop recursively through directories
        # will exclude dir if in the exclude list
        dirs[:] = [d for d in dirs if d not in exclude_directories]
        for fname in dirs:
            if 'flip' in fname:
                dir_vfas += [os.path.join(dname, fname)]
                sdir_vfas += [fname]

    vfa_Te = 2.1 * 1e-3 # (s) vfa echo time
    vfa_Tr = 4.43 * 1e-3 # (s) vfa repetition time
    vfas = [5,10,15,20,25,30] # flip angles (Deg)

    for vfa in vfas:
        for dname in dir_vfas:
            if 'ax {} flip'.format(str(vfa)) in dname:
                study_dir = dname
                for file in sdir_vfas:
                    if 'ax {} flip'.format(str(vfa)) in file:
                        fname = file
                
                        vfa_S_v1, M_vfa1 = Conc2Sig_Linear(vfa_Tr,vfa_Te,np.deg2rad(vfa),R10,R20st_value,S,2,M)
                        vfa_S_v2, M_vfa2 = Conc2Sig_Linear(vfa_Tr,vfa_Te,np.deg2rad(vfa),R10,R20st_value,S,2,M_tr)
                        vfa_S_v1_noise = addnoise(1, stdS, vfa_S_v1,Nt=datashape[-1])
                        vfa_S_v2_noise = addnoise(1, stdS, vfa_S_v2,Nt=datashape[-1])
                        vfa_S_v1_noise = vfa_S_v1_noise[:,:,:,0]
                        vfa_S_v2_noise = vfa_S_v2_noise[:,:,:,0]
                        produce_vfa_DICOM(dro_dir,study_dir,fname,dro_IDnum,1,vfa_S_v1_noise, StudyUID[0])
                        produce_vfa_DICOM(dro_dir,study_dir,fname,dro_IDnum,2,vfa_S_v2_noise, StudyUID[1])



    from DICOM_processing import changeDICOM
    import os

    # Script written to replace the DCE MRI in protocol with the DRO DCE-MRI
    # dicom files
    # Also to replace the patient identifiers and unique terms in the rest of the protocol
    ##### IMPORTANT: #####
    # move all scans form visit 2 in real data set in a 'Visit2' directory

    visit = [1,2] # list of visits (if want more will have to change study date allocation in changeDICOM)

    # state intended new patient ID number and the path to the patient study (the enclosing directory)

    # state new directory path to save the DRO outputs into


    dir_T1maps = [] # empty list to store T1 mapping driectory paths
    sdir_T1maps = [] # empty list to store t1 mapping directory names

    # Directory exclusion based on:
    # https://stackoverflow.com/questions/55394131/python-list-all-files-in-subdirectories-but-exclude-some-directories/55395209
    # Poster: Shanteshwar Inde, Mar 28 2019 at 10:22

    exclude_directories = set(['Visit2']) # Excludes Visit2 directory from os.walk
    for dname, dirs, files in os.walk(dir_prot):  #loop recursively through directories
        # will exclude dir if in the exclude list
        dirs[:] = [d for d in dirs if d not in exclude_directories]
        for fname in dirs:
            # allocates directory names based on scan type
            # known str segment in each one even with Series number and UID changing
            if 'flair' in fname:
                dir_flair = os.path.join(dname, fname)
                sdir_flair = fname
            if 'gre' in fname:
                dir_flash = os.path.join(dname, fname)
                sdir_flash = fname
            if 'tensor' in fname:
                dir_tensor = os.path.join(dname, fname)
                sdir_tensor = fname

    # run change dicom for all of the additional scans
    for n in range(0, len(visit)):
        study_flairUID = pydicom.uid.generate_uid(prefix=uidprefix)
        study_flashUID = pydicom.uid.generate_uid(prefix=uidprefix)
        
        changeDICOM(dro_dir,dir_flair,sdir_flair,dro_IDnum,visit[n],study_flairUID) 
        changeDICOM(dro_dir,dir_flash,sdir_flash,dro_IDnum,visit[n],study_flashUID) 
        changeDICOM(dro_dir,dir_tensor,sdir_tensor,dro_IDnum,visit[n],StudyUID[n]) 

save_anim_and_plot = False

if save_anim_and_plot == True:
    # Save animation dro and plots both visits 1
    displaySlices(trans_K1,'{}/'.format(dro_dir),'transK1','K1')
    displaySlices(trans_k2,'{}/'.format(dro_dir),'transk2','k2')
    displaySlices(trans_Vp,'{}/'.format(dro_dir),'transVp','Vp')

    displaySlices(mf_K1,'{}/'.format(dro_dir),'DRO_K1','K1')
    displaySlices(mf_k2,'{}/'.format(dro_dir),'DRO_k2','k2')
    displaySlices(mf_Vp,'{}/'.format(dro_dir),'DRO_Vp','Vp')

    animateSingle(S,'Data_Signal',dro_dir,datashape,dt,'Signal')
    animateSingle(dro_Snoise,'DRO_Signal_Noisy',dro_dir,datashape,dt,'Signal')
    animateSingle(Ctiss,'DRO_Ctiss',dro_dir,datashape,dt,'Concentration')

    animateSingle(dro_Snoise_tr,'Scan2_DRO_Signal_Noisy',dro_dir,datashape,dt,'Signal')
    animateSingle(Ctiss_tr,'Scan2_DRO_Ctiss',dro_dir,datashape,dt,'Concentration')

#%% Visualise output signal in random coordinate points

view_random_signal_points = False

if view_random_signal_points == True:
    import random
    a = random.randint(35,225)
    b = random.randint(50,200)
    c = random.randint(0,15)
    v1_region = np.mean(np.mean(dro_Snoise[125:150,150:165,5,:],axis=1),axis=0)
    v2_region = np.mean(np.mean(dro_Snoise_tr[125:150,150:165,5,:],axis=1),axis=0)

    with plt.style.context('ggplot'):

        plt.rcParams.update({'font.size':18})
        fig, axs = plt.subplots(nrows=1,ncols=1,figsize=(10,10))
        fig.subplots_adjust(wspace=0.2)
        
        axs.plot(t,dro_Snoise[a,b,c,:]-dro_S[a,b,c,0],'b',label='Visit 1')
        axs.plot(t,dro_Snoise_tr[a,b,c,:]-dro_S_tr[a,b,c,0],'c',label='Visit 2')
        axs.set(xlim=(0,None),xlabel='Time (mins)',ylabel='Signal Intensity',title='\n In voxel (x,y,z)= =({},{},{}) \n'.format(b,a,c))
        plt.legend()

