#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 09:11:26 2021

@author: pyess
"""

import numpy as np
import os
import nibabel as nib
from sigfig import round as sfround
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
import itertools
## Syntheitic_P1 is patient ID 'RIDER Neuro MRI-7868737135'
## Synthetic_P2 is is patient ID 'RIDE Neuro MRI-9215224289'

def nocanonical_get_array_from_nifti(pathname):
    
    nifti_file = nib.load(pathname) # load nifti file
    nifti_data = nifti_file.dataobj # extract array data from nifti file
    nifti_array = np.array(nifti_data) # convert nifti data into numpy array
    
    return nifti_array

plotting = 'Off'#

list_dir = 'entryDirectories' # path to enclosing directory for entries
entry_list = ['constantKtransModel'] # e.g. ['entry1', 'entry2', 'entry3'] # directory names for entries, can be single or multiple
entry_list_rep = ['constantKtransModel'] # e.g. ['entry1', 'entry2', 'entry3'] entry_list and entry_list_rep are not required to be identical if only some entries have reporduced data directory names for entries with reproduceability, can be single or multiple

labels = ['C1_v1','C1_v2','C2_v1','C2_v2','C3_v1','C3_v2','C4_v1','C4_v2','C5_v1','C5_v2','C6_v1','C6_v2','C7_v1','C7_v2','C8_v1','C8_v2','S1_v1','S1_v2','S2_v1','S2_v2']

dro_dir = 'DROKtransNifti' # path to stored ground truth maps
mask_dir = 'Masks' # path to stored mask files

scoring_type = 'mean'# set to 'median' to evaluate scores with median tumor Ktrans values  

Ktrans_mean_all = [] # for repeatability analysis
Ktrans_mask_allvox_rep = [] #voxelwise reproducibility analysis
Ktrans_synthetic_diff_all = [] #voxelwise accuracy analysis
all_repeat_scores = []

with open('scoringOutputs/OSIPI_scores.txt', 'w') as f:
    f.write('Scores for teams listed. \n \n \n')

with open('scoringOutputs/proportional_change_Ktrans_from_DRO.txt', 'w') as f:
    f.write('Team \t dK (prop) SP1 \t dK (prop) SP2 \n')

with open('scoringOutputs/OSIPI_score_tabular.txt', 'w') as f:
    f.write('Team \t Accuracy \t Repeatability \t Reproducibility \t OSIPI score silver \t OSIPI score gold \n')

with open('scoringOutputs/TMROI_Ktrans.txt', 'w') as f:
    f.write("Team \t"+"\t".join(labels)+"\n")
    
for entries in entry_list: # cycles through each emtrance directory from list above
    
    # walk through the directory and carry out extraction if finds .nii type files
    c_fnames = [] # to loop through correct clinical patients (repeatability)
    s_fnames = [] # to loop through correct synthetic patients (accuracy)
    all_fnames = []  # to loop through correct patients (reproducibility)
    
    for dirName, subdirs,files in os.walk('{}/{}'.format(list_dir,entries)):
        for file in files:
            if ".nii" in file.lower():  # check whether the file's nifti
                if "clinical" in file.lower(): # add to c_fnames list if file name contains 'clinical' string
                    c_fnames += [file]
                if "synthetic" in file.lower(): # add to c_fnames list  if file name contains 'synthetic' string
                    s_fnames += [file]
                all_fnames += [file]
                
    c_mask_fnames = [] # to loop through correct clinical masks (repeatability)
    s_mask_fnames = [] # to loop through correct synthetic masks (accuracy)
    all_mask_fnames = [] # to loop through correct masks (reproducibility)
    Ktrans_mask_mean = [] # for table output
    Ktrans_mask_std = [] # for table output
    Ktrans_mask_vox_rep = []  # for voxelwise reproducibility
    dKtrans_prop = [] # for proportional change analysis wrt DRO
    dKtrans_prop_gt = [] # for proportional in DRO
    
    for dirName, subdirs,files in os.walk('{}'.format(mask_dir)):
        for file in files:
            if ".nii" in file.lower():  # check whether the file's DICOM
                if "clinical" in file.lower():
                    c_mask_fnames += [file]
                if "synthetic" in file.lower():
                    s_mask_fnames += [file]
                all_mask_fnames += [file]
                    
    clinical_P = np.arange(0,8,1) # creates array with same entries as clinical patients
    
    rsum = 0
    
    if entries == 'Alice':
        factor = 1000 # scaling factor from SOP
    else:
        factor = 1
    
    for i in clinical_P: # calculate and sum for each clinical patient
        v1 = nocanonical_get_array_from_nifti('{}/{}/{}'.format(list_dir,entries, c_fnames[i*2]))/factor # read in submission NIFTI visit 1
        v2 = nocanonical_get_array_from_nifti('{}/{}/{}'.format(list_dir,entries, c_fnames[(i*2)+1]))/factor # read submission in NIFTI visit 2
        
        if entries == 'entry_requires_rotation': # Apply transformations to align with mask NIFTI
            v1 = np.rot90(v1,k=-1)
            v2 = np.rot90(v2,k=-1)
        if entries == 'entry_requires_flip':
            v1 = np.flip(v1,axis=1)
            v2 = np.flip(v2,axis=1)
            
        m_1 = nocanonical_get_array_from_nifti('{}/{}'.format(mask_dir, c_mask_fnames[i*2])) # read in mask NIFTI visit 2
        m_2 = nocanonical_get_array_from_nifti('{}/{}'.format(mask_dir, c_mask_fnames[(i*2)+1])) # read in mask NIFTI visit 2
        
        masked_v1 = v1 * m_1  
        masked_v2 = v2 * m_2
        
        if scoring_type == 'mean':
            mean_v1 = masked_v1[m_1==1].mean()
            mean_v2 = masked_v2[m_2==1].mean()
        elif scoring_type == 'median':
            mean_v1 = np.median(masked_v1[m_1==1])
            mean_v2 = np.median(masked_v2[m_2==1])
        else:
            print('No/invalid scoring type selected')
        
        Ktrans_mask_std += [sfround(masked_v1[m_1==1].std(), sigfigs=3)]
        Ktrans_mask_std += [sfround(masked_v2[m_2==1].std(), sigfigs=3)]
        
        Ktrans_mask_mean += [sfround(mean_v1, sigfigs=3)] # for Ktrans table and repeatability analysis
        Ktrans_mask_mean += [sfround(mean_v2, sigfigs=3)] # for Ktrans table and repeatability analysis
        
        mean = np.mean((mean_v1,mean_v2))
        stdev = np.std((mean_v1,mean_v2))
        frac_var = stdev/mean # coefficient of variation
        rsum += frac_var
    r_score = (np.exp(  - (rsum/len(clinical_P)) )) # caluclate repeatability score
    
    synthetic_P = np.arange(0,2,1) # creates array with same entries as synthetic patients
    
    asum = 0
    Ktrans_mask_mean_gt = ['NA',]*(2*len(clinical_P))
    total_diff = []

    for i in synthetic_P: # calculate and sum for every synthetic patient
        v1 = nocanonical_get_array_from_nifti('{}/{}/{}'.format(list_dir,entries, s_fnames[2*i]))/factor
        v2 = nocanonical_get_array_from_nifti('{}/{}/{}'.format(list_dir,entries, s_fnames[(2*i)+1]))/factor
        
        gt_v1 = nocanonical_get_array_from_nifti('{}/Synthetic_P{}_ktrans_from_dicom/Visit1/Synthetic_P{}_Visit1_Ktrans_aligned_from_dicom.nii'.format(dro_dir,i+1,i+1))/1000#[:,:,:,0]
        gt_v2 = nocanonical_get_array_from_nifti('{}/Synthetic_P{}_ktrans_from_dicom/Visit2/Synthetic_P{}_Visit2_Ktrans_aligned_from_dicom.nii'.format(dro_dir,i+1,i+1))/1000#[:,:,:,0]
        gt_v1 = np.flip(gt_v1,axis=0)
        gt_v2 = np.flip(gt_v2,axis=0)
            
        m_1 = nocanonical_get_array_from_nifti('{}/{}'.format(mask_dir, s_mask_fnames[i*2])) # load mask arrays
        m_2 = nocanonical_get_array_from_nifti('{}/{}'.format(mask_dir, s_mask_fnames[(i*2)+1]))

        # apply masks
        masked_v1 = v1 * m_1  
        masked_v2 = v2 * m_2
        masked_gt_v1 = gt_v1 * m_1
        masked_gt_v2 = gt_v2 * m_2
        
        if scoring_type == 'mean':
            mean_v1 = masked_v1[m_1==1].mean()
            mean_gtv1 = masked_gt_v1[m_1==1].mean()
        elif scoring_type == 'median':
            mean_v1 = np.median(masked_v1[m_1==1])
            mean_gtv1 = np.median(masked_gt_v1[m_1==1])
        else:
            print('No/invalid scoring type selected')
        
        mean = np.mean((mean_v1,mean_gtv1))
        stdev = np.std((mean_v1,mean_gtv1))
        frac_var_v1 = stdev/mean #coefficient of variation v1
        
        if scoring_type == 'mean':
            mean_v2 = masked_v2[m_2==1].mean()
            mean_gtv2 = masked_gt_v2[m_2==1].mean()
        elif scoring_type == 'median':
            mean_v2 = np.median(masked_v2[m_2==1])
            mean_gtv2 = np.median(masked_gt_v2[m_2==1])
        else:
            print('No/invalid scoring type selected')
        
        mean = np.mean((mean_v2,mean_gtv2))
        stdev = np.std((mean_v2,mean_gtv2))
        frac_var_v2 = stdev/mean #coefficient of variation v2
        
        asum += frac_var_v1 + frac_var_v2 # variation sum for score
        
        dK = (mean_v1-mean_v2)/mean_v1
        dK_gt = (mean_gtv1-mean_gtv2)/mean_gtv1
        
        dKtrans_prop += [dK,]
        dKtrans_prop_gt += [dK_gt,]

        Ktrans_mask_mean += [sfround(mean_v1, sigfigs=3)] # mean Ktrans for table and repeatability analysis
        Ktrans_mask_mean += [sfround(mean_v2, sigfigs=3)]
        Ktrans_mask_mean_gt += [sfround(mean_gtv1, sigfigs=3)]
        Ktrans_mask_mean_gt += [sfround(mean_gtv2, sigfigs=3)]
        
        v1_diff = (masked_v1[m_1==1]-masked_gt_v1[m_1==1])
        total_diff += [v1_diff.flatten()]

        v2_diff = (masked_v2[m_2==1]-masked_gt_v2[m_2==1])
        total_diff += [v2_diff.flatten(),]
        
    all_diff = np.concatenate( total_diff, axis=0 )
    Ktrans_synthetic_diff_all += [all_diff] # For voxelwise accuracy analysis
    
    a_score = (np.exp(  - (asum/(2*len(synthetic_P)) ) ))  # calculate accuracy score    
    
    if entries in entry_list_rep:
        
        all_P = np.arange(0,10,1)
        
        Ktrans_mask_mean_rep = []
        
        repro_sum =0
        for i in all_P:
            v1 = nocanonical_get_array_from_nifti('{}/{}/{}'.format(list_dir,entries, all_fnames[2*i]))/factor
            v2 = nocanonical_get_array_from_nifti('{}/{}/{}'.format(list_dir,entries, all_fnames[(2*i)+1]))/factor
            repro_v1 = nocanonical_get_array_from_nifti('{}/{}_neutral/{}'.format(list_dir,entries, all_fnames[2*i]))/factor
            repro_v2 = nocanonical_get_array_from_nifti('{}/{}_neutral/{}'.format(list_dir,entries, all_fnames[(2*i)+1]))/factor
                    
            m_1 = nocanonical_get_array_from_nifti('{}/{}'.format(mask_dir, all_mask_fnames[i*2])) # load mask arrays
            m_2 = nocanonical_get_array_from_nifti('{}/{}'.format(mask_dir, all_mask_fnames[(i*2)+1]))
 
            masked_v1 = v1 * m_1  
            masked_v2 = v2 * m_2    
            masked_repro_v1 = repro_v1 * m_1
            masked_repro_v2 = repro_v2 * m_2
        
            if scoring_type == 'mean':
                mean_v1 = masked_v1[m_1==1].mean()
                mean_repro_v1 = masked_repro_v1[m_1==1].mean()
            elif scoring_type == 'median':
                mean_v1 = np.median(masked_v1[m_1==1])
                mean_repro_v1 = np.median(masked_repro_v1[m_1==1])
            else:
                print('No/invalid scoring type selected')
            mean = np.mean((mean_v1,mean_repro_v1))
            stdev = np.std((mean_v1,mean_repro_v1))
            frac_var_v1 = stdev/mean
            
            if scoring_type == 'mean':
                mean_v2 = masked_v2[m_2==1].mean()
                mean_repro_v2 = masked_repro_v2[m_2==1].mean()
            elif scoring_type == 'median':
                mean_v2 = np.median(masked_v2[m_2==1])
                mean_repro_v2 = np.median(masked_repro_v2[m_2==1])
            else:
                print('No/invalid scoring type selected')
                
            mean = np.mean((mean_v2,mean_repro_v2))
            stdev = np.std((mean_v2,mean_repro_v2))
            frac_var_v2 = stdev/mean
            repro_sum += frac_var_v1 + frac_var_v2
            
            repro_diff_v1 = masked_v1-masked_repro_v1 # for voxelwise reproducibility
            repro_diff_v2 = masked_v2-masked_repro_v2 # for voxelwise reproducibility
            
            
            Ktrans_mask_mean_rep += [sfround(mean_repro_v1, sigfigs=3)] # for Ktrans table
            Ktrans_mask_mean_rep += [sfround(mean_repro_v2, sigfigs=3)]# for Ktrans table
            Ktrans_mask_vox_rep += [repro_diff_v1[m_1==1].flatten()] #for voxelwise reproducibility
            Ktrans_mask_vox_rep += [repro_diff_v2[m_2==1].flatten()]#for voxelwise reproducibility

        repro_score = (np.exp( - (repro_sum/(2*len(all_P)) ) ) ) # calculate reproducibility score        

    else:
        repro_score = np.nan
        Ktrans_mask_mean_rep = [] # for Ktrans table
        Ktrans_mask_mean_rep += ['NA']*20
    
    Ktrans_mask_allvox_rep += [Ktrans_mask_vox_rep] #for voxelwise reproducibility analysis
    Ktrans_mean_all += [Ktrans_mask_mean,] #For repeatability visit analysis
    all_repeat_scores += [r_score,]
    # Write output file 
    with open('scoringOutputs/OSIPI_scores.txt', 'a') as f: # append scores into file
        f.write("Entry team: "+"{}".format(entries)+"\n")
        f.write("Reproducibility Score: "+"{:.3f}".format(repro_score)+ "\n")
        f.write("Repeatability Score: "+ "{:.3f}".format(r_score)+"\n")
        f.write("Accuracy Score: " +"{:.3f}".format(a_score) + "\n")
        f.write("OSIPI Score silver: "+"{:.1f}".format((r_score * a_score )*100) +"% " +" \n")
        f.write("OSIPI Score gold: "+"{:.1f}".format((r_score * a_score*repro_score )*100) +"% " +" \n \n \n")
    with open('scoringOutputs/OSIPI_score_tabular.txt', 'a') as f:
        f.write('{} \t {:.3f} \t {:.3f} \t {:.3f} \t {:.1f} \t {:.1f} \n'.format(entries,a_score,r_score,repro_score,(r_score * a_score )*100,(r_score * a_score *repro_score)*100))
    
    with open('scoringOutputs/TMROI_Ktrans.txt', 'a') as f:
        if entries == entry_list[0]:
            f.write("gt \t"+"\t".join([str(a) for a in Ktrans_mask_mean_gt])+"\n")
        f.write("{}_entry \t".format(entries)+"\t".join([str(a) for a in Ktrans_mask_mean])+"\n")
        f.write("{}_sd \t".format(entries)+"\t".join([str(a) for a in Ktrans_mask_std])+"\n")
        f.write("{}_repro \t".format(entries)+"\t".join([str(a) for a in Ktrans_mask_mean_rep])+"\n")
    
    with open('scoringOutputs/proportional_change_Ktrans_from_DRO.txt', 'a') as f:
        if entries == entry_list[0]:
            f.write('DRO \t {:.3f} \t {:.3f} \n'.format(dKtrans_prop_gt[0], dKtrans_prop_gt[1]))
        f.write('{} \t {:.3f} \t {:.3f} \n'.format(entries,dKtrans_prop[0], dKtrans_prop[1]))

#%% Reproducibility - all data points

with open('scoringOutputs/reproducability_statistics.txt', 'w') as f:
    f.write('Team \t repro mean diff \t repro SD diff \t repro median diff \t repro LQ diff  \t  repro UQ diff \n')


for j in range(0,len(entry_list)):
    boxplot_list = list(itertools.chain(*Ktrans_mask_allvox_rep[j]))
    with open('scoringOutputs/reproducability_statistics.txt', 'a') as f:
        f.write('{} \t {:.3e} \t {:.3e} \t {:.3e} \t {:.3e} \t {:.3e}\n'.format(entry_list[j-1],np.mean(np.asarray(boxplot_list)),np.std(np.asarray(boxplot_list)),np.median(np.asarray(boxplot_list)),np.percentile(np.asarray(boxplot_list),25),np.percentile(np.asarray(boxplot_list),75)))

if plotting == 'On':
    print('Plotting enabled...')
    #%% Accuracy diffs
    colors = sns.color_palette("tab20b", n_colors=12).as_hex()
    adj_colors = np.delete(np.asarray(colors),(0,len(entry_list)+1))
    plt.rcParams.update({'font.size':26})
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,12),dpi=100,sharey=True, sharex=True)
    
    medianprops = dict(linestyle='-', linewidth=2.5, color='black')
    meanpointprops = dict(marker='o', markeredgecolor='black',markerfacecolor='firebrick')
    bplot=ax.boxplot(Ktrans_synthetic_diff_all,zorder=2,vert = False,whis=(5,95),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True, positions=np.flip(range(1,len(entry_list)+1)))
    for patch, color in zip(bplot['boxes'], adj_colors):
                patch.set_facecolor(color)
    ax.set_xlim([-0.3,0.5])
    ax.set_xlabel('$\Delta K^{trans}$')
    ax.axvline(x=0,linestyle='--',color='k',zorder=1)
    plt.yticks(np.flip(np.arange(1, len(entry_list)+1,1)), entry_list)
    fig.savefig('scoringOutputs/accuracy_voxdiff_boxplot.eps', bbox_inches='tight')
    fig.savefig('scoringOutputs/accuracy_voxdiff_boxplot.png', bbox_inches='tight')
    plt.show()
    
    #%% Repeatability
    diff_means = []
    diff_stdevs = []
    Ktrans_diff_all = []
    RCperc = []
    for j in range(0,len(entry_list)):
        Ktrans_diff = []
        wCV = []
        for i in np.arange(0, 15, 2):
            v1 = Ktrans_mean_all[j][i]
            v2 = Ktrans_mean_all[j][i+1]
            Ktrans_diff += [(v1-v2),]
            wCV += [np.std((v1,v2))**2/np.mean((v1,v2))**2,]
        current_RCperc = 2.77*(100*(np.sqrt(np.mean(np.asarray(wCV)))))
        RCperc += [current_RCperc,]
        Ktrans_diff_all +=[Ktrans_diff,]
        Ktrans_diff = np.asarray(Ktrans_diff)
        diff_means += [np.mean(Ktrans_diff),]
        diff_stdevs += [np.std(Ktrans_diff),]
        
    if len(RCperc) >= 2:
        Pcoeff = scipy.stats.pearsonr(RCperc,all_repeat_scores)
        print('r = {}, p = {}'.format(Pcoeff[0],Pcoeff[1]))
    else:
        print('Insufficient entry teams to calculate correlation.')
    
    plt.rcParams.update({'font.size':26})
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,12),dpi=100,sharey=True, sharex=True)
    
    medianprops = dict(linestyle='-', linewidth=2.5, color='black')
    meanpointprops = dict(marker='o', markeredgecolor='black',markerfacecolor='firebrick')
    bplot=ax.boxplot(Ktrans_diff_all,vert = False,whis=(5,95),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True, positions=np.flip(range(1,len(entry_list)+1)))
    for patch, color in zip(bplot['boxes'], adj_colors):
                patch.set_facecolor(color)
    ax.set_xlim([-0.1,0.1])
    ax.axvline(x=0,linestyle='--',color='k',zorder=1)
    ax.set_xlabel('$K^{trans}$ change between clinical visits',fontsize=26)
    plt.yticks(np.flipud(np.arange(1,len(entry_list)+1,1)), entry_list)
    
    fig.savefig('scoringOutputs/clinical_patients_repeatability_boxplot.eps', bbox_inches='tight')
    fig.savefig('scoringOutputs/clinical_patients_repeatability_boxplot.png', bbox_inches='tight')
    
    plt.show()

else:
    print('Finished! To display plots run with plotting to enabled.')
    
print('View scoringOutputs for all generated ouputs')