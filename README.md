# ISMRM-OSIPI DCE-MRI Challenge Repository

This repository accompanies the ISMRM-OSIPI DCE-MRI Challenge and provides supplementary data from the manuscript. OSIPI TF6.2 designed this OSIPI-DCE challenge for evaluating DCE-MRI methods that estimate the volume transfer constant (Ktrans) in brain tumors.

A set of parameter maps and scripts for analyzing the challenge data are provided withi this repository. For access to the challenge guidelines and challenge DICOM data files see the OSF repository (https://osf.io/u7a6f).

## Challenge Overview

 Analysis methods were submitted by competing teams which were evaluated and compared with each other for quantification of Ktrans from DCE-MRI in terms of accuracy, repeatability, and reproducibility. The accuracy of Ktrans estimations by the participating teams were scored based on synthetic data, i.e., 2 digital reference objects (DROs) each including 2 visits, specifically designed for this challenge; repeatability was rated based on test-retest DCE-MRI scans of 8 patients from RIDER Neuro MRI Database; reproducibility was assessed based on an independent analysis by a neutral team.

The submissions were then ranked according to a global score reflecting that an ideal method should be accurate AND repeatable AND reproducible.

## Repository Overview
This data and code collection includes:
- The Digital Reference Object parameter maps used the simulate the synthetic patient data
    * NIfTI files of the Ktrans maps: DROKtransNifti
    * Compressed numpy array files: pythonArraysDRO
- Challenge DICOM data:
    * Synthetic patient data provided for the challenge: SyntheticDicom
- Scoring script used for the challenge:
    * challengeScoring.py: Just scores caluclated with plotting options.
- Nifti mask files to extract tumor regions-of-interest defined for the challenge: Masks

## Parameter Maps

All the parameter maps used in this repository can be found in the syntehtic data part of the repository [here](SyntheticData). The are serveral data types available: [NIfII](SyntheticData/NIfTI) or [numpy array files](SyntheticData/pythonArraysDRO). These maps include:

-**Ktrans**
- **Plasma Volume Fraction (vp)**
- **Leakage Rate (kep)**
- **Arterial Input Function (AIF)**
-**Precontrast Longitudinal Relaxation Rate (R10)**
- **Precontrast Transverse Relaxation Rate (R*20)**

## Mask NIfTI Files

The mask NIfTI files, which define the regions of interest for the analysis, can be found in the following file: `[path/to/mask/files]`. These mask files should be used to restrict the analysis to specific regions of interest within the DCE-MRI data.

## Scoring Script

A scoring script is provided in this repository titled `scoring_script.py`. This script is designed to evaluate the performance of different methods on the challenge data. However, please note that using this script requires a fully analyzed dataset consisting of Ktrans NIfTI files for each synthetic and clinical patient.

To use the scoring script, follow these steps:

1. Ensure that you have a fully analyzed dataset with Ktrans NIfTI files for each synthetic and clinical patient.
2. Modify the script as needed to provide the necessary paths to your Ktrans NIfTI files.
3. Run the script using a compatible Python environment to obtain the evaluation scores for your method.

Please note that the scoring script assumes the availability of the required data files and follows specific guidelines outlined

## Scoring Script Usage

A scoring script is provided in this repository [here](Scoring/challengeScoring.py). This script is designed to evaluate the performance of different Ktrans submissions for the challenge data. However, please note that using this script requires a fully analyzed dataset consisting of Ktrans NIfTI files for each synthetic and clinical patient as seen in uploaded [test folder](Scoring/entryDirectories/constantKtransModel).

1. Clone or pull this directory
2. Analyse the Challenge data set
    - Find the challenge data for download [here](https://osf.io/u7a6f/files)
3. Save Nifti files in specifc name format which is different for synthetic and clinical data sets:
    - Synthetic_P#_Visit#.nii
    - Clinical_P#_Visit#.nii
    - See an example directory in the correct formatting in the constant Ktrans [test case](Scoring/entryDirectories/constantKtransModel)
4. Add your directory to entry_list in the scroring script
5. Run!

