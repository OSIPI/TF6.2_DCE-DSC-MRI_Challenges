# ISMRM-OSIPI DCE-MRI Challenge Repository

This repository accompanies the ISMRM-OSIPI DCE-MRI Challenge and provides supplementary data from the manuscript. OSIPI TF6.2 designed this OSIPI-DCE challenge for evaluating DCE-MRI methods that estimate the volume transfer constant ($K^{trans}$) in brain tumors.

A set of parameter maps and scripts for analyzing the challenge data are provided withi this repository. For access to the challenge guidelines and challenge DICOM data files see the OSF repository (https://osf.io/u7a6f).

# Challenge Overview

 Analysis methods were submitted by competing teams which were evaluated and compared with each other for quantification of $K^{trans}$ from DCE-MRI in terms of accuracy, repeatability, and reproducibility. The accuracy of $K^{trans}$ estimations by the participating teams were scored based on synthetic data, i.e., 2 digital reference objects (DROs) each including 2 visits, specifically designed for this challenge; repeatability was rated based on test-retest DCE-MRI scans of 8 patients from RIDER Neuro MRI Database; reproducibility was assessed based on an independent analysis by a neutral team.

The submissions were then ranked according to a global score reflecting that an ideal method should be accurate AND repeatable AND reproducible.

## Repository Overview
This data and code collection includes:
- The Digital Reference Object parameter maps used the simulate the synthetic patient data:
    * [NIfTI files](SyntheticData/NIfTI)
    * Compressed [numpy array files](SyntheticData/pythonArraysDRO)
- Challenge DICOM data:
    * [Synthetic patient data](SyntheticData/SyntheticDicom) provided for the challenge.
- [Scoring script](Scoring/challengeScoring.py) used for the challenge:
    * This calculates scores with plotting options.
- [NIfTI mask files](Scoring/Masks) to extract tumor regions-of-interest defined for the challenge

### Parameter Maps

All the parameter maps used in this repository can be found in the syntehtic data part of the repository [here](SyntheticData). The are serveral data types available: [NIfTI](SyntheticData/NIfTI) or [numpy array files](SyntheticData/pythonArraysDRO). These maps include:

- **Volume transfer constant ($K^{trans}$)**
- **Plasma Volume Fraction ($v_{p}$)**
- **Leakage Rate (kep)**
- **Arterial Input Function (AIF)**
- **Precontrast Longitudinal Relaxation Rate ($R_{10}$)**
- **Precontrast Transverse Relaxation Rate ($R^{*}_{20}$)**

### Mask NIfTI Files

The mask NIfTI files, which define the regions of interest for the analysis, can be found [here](Scoring/Masks). These mask files should be used to restrict the scoring area to specific regions of interest within the DCE-MRI data.

## Scoring Script

A scoring script is provided in this repository titled [challengeScoring.py](Scoring/challengeScoring.py). This script is designed to evaluate the performance of different methods on the challenge data. However, please note that using this script requires a fully analyzed dataset consisting of Ktrans NIfTI files for each synthetic and clinical patient.

To use the scoring script, follow these steps:

1. Ensure that you have a fully analyzed dataset with $K^{trans}$ NIfTI files for each synthetic and clinical patient.
2. Modify the script as needed to provide the necessary paths to your $K^{trans}$ NIfTI files.
    - Recommended to save as inside a 'TeamName' directory and seperate 'TeamName_neutral' directory for reproduced entry values.
    - See example directories in the correct formatting in the constant $K^{trans}$ [test case](Scoring/entryDirectories).
3. Run the script using a compatible Python environment to obtain the evaluation scores for your method.

Please note that the scoring script assumes the availability of the required data files and follows specific guidelines outlined

## Scoring Script Usage

A scoring script is provided in this repository: [challengeScoring.py](Scoring/challengeScoring.py). This script is designed to evaluate the performance of different $K^{trans}$ submissions for the challenge data. However, please note that using this script requires a fully analyzed dataset consisting of $K^{trans}$ NIfTI files for each synthetic and clinical patient as seen in uploaded [test folder](Scoring/entryDirectories/constantKtransModel).

To use the scoring script, follow these steps:

1. Clone or pull this directory.
2. Analyse the Challenge data set.
    - Find the challenge data for download [here](https://osf.io/u7a6f/files).
    - Follow the Challenge guidelines [here](https://osf.io/qagc3).
3. Save your $K^{trans}$ NIfTI map in specifc naming format which is different for synthetic and clinical data sets:
    - Synthetic_P#_Visit#.nii; Clinical_P#_Visit#.nii
    - See an example directory in the correct formatting in the constant $K^{trans}$ [test case](Scoring/entryDirectories/constantKtransModel).
    - If files are not named correctly the scoring script will not run as expected.
4. Ensure that you have a fully analyzed dataset with $K^{trans}$ NIfTI files for each synthetic and clinical patient.
5. Add your directory to [entryDirectories](Scoring/entryDirectories) in the [Scoring](Scoring) location.
    - Recommended to save as inside a 'TeamName' directory and seperate 'TeamName_neutral' directory for reproduced entry values.
    - See example directories in the correct formatting in the constant $K^{trans}$ [test case](Scoring/entryDirectories).
6. Modify the script as needed to provide the necessary paths to your $K^{trans}$ NIfTI files.
    - Add yout 'TeamNAme' that you have chosen for the directory in [entryDirectories](Scoring/entryDirectories) to `entry_list`.
7. Run the script using a compatible Python environment to obtain the evaluation scores for your method.

### Python Package Requirements
- numpy 1.20.1
- matplotlib 3.3.4
- pydicom 2.1.2
- scipy 1.6.1
- nibabel 3.2.1


