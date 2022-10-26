# DCE-DSC-MRI_Challenges

This repository houses various supplementary data from the OSIPI DCE-MRI Challenge. (https://osf.io/u7a6f/wiki/home/)

OSIPI TF6.2 designed the DCE (OSIPI-DCE) challenge for evaluating DCE-MRI methods that estimate the volume transfer constant (Ktrans) in brain tumors.

## Challenge Overview

The analysis methods submitted by the competing teams were be evaluated and compared with each other for quantification of Ktrans from DCE-MRI in terms of accuracy, repeatability, and reproducibility. The accuracy of Ktrans estimations by the participating teams were scored based on synthetic data, i.e., 2 digital reference objects (DROs) each including 2 visits, specifically designed for this challenge; repeatability was rated based on test-retest DCE-MRI scans of 8 patients from RIDER Neuro MRI Database; reproducibility was assessed based on an independent analysis by a neutral team.

The submissions were be ranked according to a global score reflecting that an ideal method should be accurate AND repeatable AND reproducible.

## Repository Overview
This data and code collection includes:
- The Digital Reference Object parameter maps used the simulate the synthetic patient data
    * Nifti files of the Ktrans maps: DROKtransNifti
    * Compressed numpy array files: pythonArraysDRO
- Challenge DICOM data:
  * Synthetic patinet data provided for the challenge: SyntheticDicom
  - Scoring script used for the challenge:
    * challengeScoring.py: Just scores caluclated with plotting options.
      - Nifti mask files to extract tumor regions-of-interest defined for the challenge: Masks

## Scoring Script Usage
1. Clone of pull this directory
2. Analyse the Challenge data set
3. Save Nifti files in format dependent if synthetic or clinical data sets:
    - Synthetic_P#_Visit#.nii
    - Clinical_P#_Visit#.nii
4. Add your directory to entry_list in the scroring script
5. Run!

