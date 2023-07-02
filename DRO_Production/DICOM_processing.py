#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:35:06 2021

@author: eveshalom
"""

import os
import pydicom
import numpy as np
from plots_and_filters import MedianFilter

def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise
    return

def dicomdirscan(dirpath):
    dicomRawList = []  # create an empty list to hold raw dcm
    dicomSliceLoc = [] # create an empty list to hold slice locations
    dicomAcqTime = [] # create an empty list to hold acquisition times
    
    # walk through the directory and carry out extraction if finds .dcm type files
    for dirName, subdir,files in os.walk(dirpath):
        for file in files:
            if ".dcm" in file.lower():  # check whether the file's DICOM
                ds = pydicom.read_file(os.path.join(dirName,file)) #pydicom used to read file
                
                dicomSliceLoc.append(float(ds.SliceLocation)) # save all slice locations
                dicomAcqTime.append(float(ds.AcquisitionTime)) # save all acquisition times
                dicomRawList.append(ds) # lsit of all read in dicom files
    
    Refdcm = dicomRawList[0].copy() # save a reference dicom (use for data shape purposes)

    
    return dicomRawList,dicomSliceLoc, dicomAcqTime,Refdcm

def Slice_Time_Readings(dicomSliceLoc, dicomAcqTime):
    
    SliceSet = list(set(dicomSliceLoc)) # list all the distinct slice location in data set
    SliceSet = sorted(SliceSet) # Slice locations sorted in asending order
    
    Acqs = list(set(dicomAcqTime)) # list all distinct acquisition times
    Acqs = sorted(Acqs) # Acquisiton times in asending order
    
    Slices={} # Create empty dictionary
    for i in range(1,len(SliceSet)+1):
        Slices.update({'{}'.format(i):[]}) # Creat dictionary enteries for each slice e.g: Slice1 ... SliceN
    
    return Slices, Acqs,SliceSet

def SpatiotempSorter(dicomRawList,Slices,SliceSet):
    for snap in dicomRawList:
        currentSlice=float(snap.SliceLocation) # finds slice location of current dicom
        i=SliceSet.index(currentSlice) # seraches for position of current slice in sliceset
        Slices['{}'.format(i+1)].append(snap) # Adds dicom to correct slice

    for entries in Slices:
        Slices['{}'.format(entries)].sort(key=lambda x: float(x.AcquisitionNumber)) # sorts each slice into time series
    
    return Slices

def SignalEnhancementExtract(dicomdata,datashape,baselinepoints):
    
    S = np.zeros(datashape) # Aloocate 4d np array of correct size
    
    #Go through Slice dictionary sequentially to fill signal array
    z=0
    for entries in dicomdata: # Every slice is a dictionary entry
        t=0
        for snap in dicomdata['{}'.format(entries)]: # Time points are sequentially stored
            S[:,:,z,t] = snap.pixel_array # output signal pixel_array into np array
            t += 1
        z += 1
    
    #Take baseline average
    S0=np.average(S[:,:,:,0:baselinepoints],axis=3) # Take baseline signal
    E=np.zeros_like(S)
    
    #Calcualte siganl enhancement
    for i in range(0,datashape[-1]):
        E[:,:,:,i]=S[:,:,:,i]-S0
        E[:,:,:,i]=MedianFilter(E[:,:,:,i]) # Median filter size (3,3)
        
    return E, S0, S

def produceDICOMktrans(encdir,ogdicom,K,visit,PatientIDNum,StudyInstanceUID):
    newDicom = ogdicom.copy()
    PatientID = 'RIDER Neuro MRI-{}'.format(PatientIDNum)
    
    uidprefix = '1.3.6.1.4.1.9328.50.16.'
    SeriesInstanceUID = pydicom.uid.generate_uid(prefix=uidprefix)
    StorageMediaFileSetUID = pydicom.uid.generate_uid(prefix=uidprefix)
    FrameOfReferenceUID = pydicom.uid.generate_uid(prefix=uidprefix)
    SeriesID = SeriesInstanceUID[-5:]
    Seriesdir = '{}.000000-perfusion-{}'.format(ogdicom[0].SeriesNumber,SeriesID)
    
    if visit == 1:
        StudyDate = '19040321'
        DateID = StudyInstanceUID[-5:]
        Datedir = '21-03-1904-BRAINRESEARCH-{}'.format(DateID)
        
    if visit == 2:
        StudyDate = '19040323'
        DateID = StudyInstanceUID[-5:]
        Datedir = '23-03-1904-BRAINRESEARCH-{}'.format(DateID)
    
    mkdir_p('{}/{}/{}/{}'.format(encdir,PatientID,Datedir,Seriesdir))

    z=0
    for snap in newDicom:
        SOPInstanceUID = pydicom.uid.generate_uid(prefix=uidprefix)
        snap.PatientID = PatientID
        snap.SOPInstanceUID = SOPInstanceUID
        snap.StudyInstanceUID = StudyInstanceUID
        snap.SeriesInstanceUID = SeriesInstanceUID
        snap.StorageMediaFileSetUID = StorageMediaFileSetUID
        snap.FrameOfReferenceUID = FrameOfReferenceUID
        snap.StudyDate = StudyDate
        snap.ContentDate = StudyDate
        Kout = abs(K[:,:,z])*1000
        Kout = Kout.astype(np.uint16)
        snap.PixelData = Kout.tobytes()
        fname = str(0+1).zfill(2)+'-'+str(z+1).zfill(4)
        snap.save_as("{}/{}/{}/{}/{}.dcm".format(encdir,PatientID,Datedir,Seriesdir,fname), write_like_original=False)
        z += 1

    return

def produceDICOM(encdir,ogdicom,S,visit,PatientIDNum,StudyInstanceUID):
    newDicom = ogdicom.copy()
    PatientID = 'RIDER Neuro MRI-{}'.format(PatientIDNum)
    
    uidprefix = '1.3.6.1.4.1.9328.50.16.'
    SeriesInstanceUID = pydicom.uid.generate_uid(prefix=uidprefix)
    StorageMediaFileSetUID = pydicom.uid.generate_uid(prefix=uidprefix)
    FrameOfReferenceUID = pydicom.uid.generate_uid(prefix=uidprefix)
    SeriesID = SeriesInstanceUID[-5:]
    Seriesdir = '{}.000000-perfusion-{}'.format(ogdicom['1'][0].SeriesNumber,SeriesID)
    if visit == 1:
        StudyDate = '19040321'
        DateID = StudyInstanceUID[-5:]
        Datedir = '21-03-1904-BRAINRESEARCH-{}'.format(DateID)
    if visit == 2:
        StudyDate = '19040323'
        DateID = StudyInstanceUID[-5:]
        Datedir = '23-03-1904-BRAINRESEARCH-{}'.format(DateID)
    
    mkdir_p('{}/{}/{}/{}'.format(encdir,PatientID,Datedir,Seriesdir))

    z=0
    for entries in newDicom:
        t=0
        for snap in newDicom['{}'.format(entries)]:
            SOPInstanceUID = pydicom.uid.generate_uid(prefix=uidprefix)
            snap.PatientID = PatientID
            snap.SOPInstanceUID = SOPInstanceUID
            snap.StudyInstanceUID = StudyInstanceUID
            snap.SeriesInstanceUID = SeriesInstanceUID
            snap.StorageMediaFileSetUID = StorageMediaFileSetUID
            snap.FrameOfReferenceUID = FrameOfReferenceUID
            snap.StudyDate = StudyDate
            snap.ContentDate = StudyDate
            Sout = abs(S[:,:,z,t])
            Sout = Sout.astype(np.uint16)
            snap.PixelData = Sout.tobytes()
            fname = str(t+1).zfill(2)+'-'+str(z+1).zfill(4)
            snap.save_as("{}/{}/{}/{}/{}.dcm".format(encdir,PatientID,Datedir,Seriesdir,fname), write_like_original=False)
            t += 1
        z += 1
    return

def produce_vfa_DICOM(dro_dir,study_dir,fname,PatientIDNum,visit,S,StudyInstanceUID):
    import random
    from DICOM_processing import mkdir_p
    import pydicom
    
    PatientID = 'RIDER Neuro MRI-{}'.format(PatientIDNum)
    uidprefix = '1.3.6.1.4.1.9328.50.16.'
    SeriesInstanceUID = pydicom.uid.generate_uid(prefix=uidprefix)
    StorageMediaFileSetUID = pydicom.uid.generate_uid(prefix=uidprefix)
    FrameOfReferenceUID = pydicom.uid.generate_uid(prefix=uidprefix)
    SeriesID = SeriesInstanceUID[-5:]
    Seriesdir = fname.replace(fname[-5:],SeriesID)
    
    if visit == 1:
        StudyDate = '19040321'
        DateID = StudyInstanceUID[-5:]
        Datedir = '21-03-1904-BRAINRESEARCH-{}'.format(DateID)
    if visit == 2:
        StudyDate = '19040323'
        DateID = StudyInstanceUID[-5:]
        Datedir = '23-03-1904-BRAINRESEARCH-{}'.format(DateID)
    
    newpath = '{}/{}/{}/{}'.format(dro_dir,PatientID,Datedir,Seriesdir)
    mkdir_p(newpath)
    
    for dirName, subdir,files in os.walk(study_dir):
        for file in files:
            if ".dcm" in file.lower():  # check whether the file's DICOM
                ds = pydicom.read_file(os.path.join(dirName,file)) #pydicom used to read file
                SOPInstanceUID = pydicom.uid.generate_uid(prefix=uidprefix)
                ds.PatientID = PatientID
                ds.SOPInstanceUID = SOPInstanceUID
                ds.StudyInstanceUID = StudyInstanceUID
                ds.SeriesInstanceUID = SeriesInstanceUID
                ds.StorageMediaFileSetUID = StorageMediaFileSetUID
                ds.FrameOfReferenceUID = FrameOfReferenceUID
                ds.StudyDate = StudyDate
                ds.ContentDate = StudyDate
                z = int(ds.InstanceNumber)-1
                Sout = abs(S[:,:,z])
                Sout = Sout.astype(np.uint16)
                ds.PixelData = Sout.tobytes()
                ds.save_as("{}/{}".format(newpath,file), write_like_original=False)
    return

def changeDICOM(dro_dir,study_dir,fname,PatientIDNum,visit, StudyInstanceUID):
    import random
    from DICOM_processing import mkdir_p
    import pydicom
    
    PatientID = 'RIDER Neuro MRI-{}'.format(PatientIDNum)
    uidprefix = '1.3.6.1.4.1.9328.50.16.'
    SeriesInstanceUID = pydicom.uid.generate_uid(prefix=uidprefix)
    StorageMediaFileSetUID = pydicom.uid.generate_uid(prefix=uidprefix)
    FrameOfReferenceUID = pydicom.uid.generate_uid(prefix=uidprefix)
    SeriesID = SeriesInstanceUID[-5:]
    Seriesdir = fname.replace(fname[-5:],SeriesID)
    
    if visit == 1:
        StudyDate = '19040321'
        DateID = StudyInstanceUID[-5:]
        Datedir = '21-03-1904-BRAINRESEARCH-{}'.format(DateID)
    if visit == 2:
        StudyDate = '19040323'
        DateID = StudyInstanceUID[-5:]
        Datedir = '23-03-1904-BRAINRESEARCH-{}'.format(DateID)
    
    newpath = '{}/{}/{}/{}'.format(dro_dir,PatientID,Datedir,Seriesdir)
    mkdir_p(newpath)
    
    for dirName, subdir,files in os.walk(study_dir):
        for file in files:
            if ".dcm" in file.lower():  # check whether the file's DICOM
                ds = pydicom.read_file(os.path.join(dirName,file)) #pydicom used to read file
                SOPInstanceUID = pydicom.uid.generate_uid(prefix=uidprefix)
                ds.PatientID = PatientID
                ds.SOPInstanceUID = SOPInstanceUID
                ds.StudyInstanceUID = StudyInstanceUID
                ds.SeriesInstanceUID = SeriesInstanceUID
                ds.StorageMediaFileSetUID = StorageMediaFileSetUID
                ds.FrameOfReferenceUID = FrameOfReferenceUID
                ds.StudyDate = StudyDate
                ds.ContentDate = StudyDate
                ds.save_as("{}/{}".format(newpath,file), write_like_original=False)
    return