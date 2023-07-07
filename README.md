# ISMRM-OSIPI DCE-MRI Challenge Repository

This repository accompanies the ISMRM-OSIPI DCE-MRI Challenge and provides supplementary data from the manuscript. OSIPI TF6.2 designed this OSIPI-DCE challenge for evaluating DCE-MRI methods that estimate the volume transfer constant ($K^{trans}$) in brain tumors.

A set of parameter maps and scripts for analyzing the challenge data are provided withi this repository. For access to the challenge guidelines and challenge DICOM data files see the OSF repository (https://osf.io/u7a6f).

# Challenge Overview

 Analysis methods were submitted by competing teams which were evaluated and compared with each other for quantification of $K^{trans}$ from DCE-MRI in terms of accuracy, repeatability, and reproducibility. The accuracy of $K^{trans}$ estimations by the participating teams were scored based on synthetic data, i.e., 2 digital reference objects (DROs) each including 2 visits, specifically designed for this challenge; repeatability was rated based on test-retest DCE-MRI scans of 8 patients from [RIDER Neuro MRI Database](https://wiki.cancerimagingarchive.net/display/Public/RIDER+NEURO+MRI); reproducibility was assessed based on an independent analysis by a neutral team.

The submissions were then ranked according to a global score reflecting that an ideal method should be accurate AND repeatable AND reproducible.

# Repository Overview
This data and code collection includes:
- The Digital Reference Object parameter maps used the simulate the synthetic patient data:
    * [NIfTI files](additionalDROData/NIfTI)
    * Compressed [numpy array files](additionalDROData/pythonArraysDRO)
- Challenge DICOM data:
    * [Synthetic patient data](ChallengeDICOMData/Synthetic_Data) provided for the challenge.
    * [Clinical patient data](ChallengeDICOMData/Clinical_Data) provided for the challenge.
    * These data are based on (Synthetic) or directly (Clinical) from the [RIDER NEURO Database from TCIA](https://wiki.cancerimagingarchive.net/display/Public/RIDER+NEURO+MRI)
        * Barboriak, D. (2015). Data From RIDER NEURO MRI. The Cancer Imaging Archive. https://doi.org/10.7937/K9/TCIA.2015.VOSN3HN1
        * Clark, K., Vendt, B., Smith, K., Freymann, J., Kirby, J., Koppel, P., Moore, S., Phillips, S., Maffitt, D., Pringle, M., Tarbox, L., & Prior, F. (2013). The Cancer Imaging Archive (TCIA): Maintaining and Operating a Public Information Repository. In Journal of Digital Imaging (Vol. 26, Issue 6, pp. 1045–1057). Springer Science and Business Media LLC. https://doi.org/10.1007/s10278-013-9622-7 PMCID: PMC3824915
- [Scoring script](Scoring/challengeScoring.py) used for the challenge:
    * This calculates scores with plotting options.
- [NIfTI mask files](Scoring/Masks) to extract tumor regions-of-interest defined for the challenge
- [DRO production code](DRO_Production/main.py) used to produce the synthetic data used in the challenge:
    * All functions used can be found within the [enclosing folder](DRO_Production)
    * To run using challenge ROIs for AIF please unzip the .npy files in the [AIFmasksforDRO](DRO_Production/AIFmasksforDRO) directory

## Parameter Maps

All the parameter maps used in this repository can be found in the synthetic data part of the repository [here](additionalDROData). The are serveral data types available: [NIfTI](additionalDROData/NIfTI) or [numpy array files](additionalDROData/pythonArraysDRO). These maps include $K^{trans}$, $v_{p}$, $k_{ep}$, AIF and $R_{10}$.

## Mask NIfTI Files

The mask NIfTI files, which define the regions of interest for the analysis, can be found [here](Scoring/Masks). These mask files should be used to restrict the scoring area to specific regions of interest within the DCE-MRI data.

# Scoring Script

A scoring script is provided in this repository titled [challengeScoring.py](Scoring/challengeScoring.py). This script is designed to evaluate the performance of different methods on the challenge data. However, please note that using this script requires a fully analyzed dataset consisting of Ktrans NIfTI files for each synthetic and clinical patient.

To use the scoring script, follow these steps:

1. Ensure that you have a fully analyzed dataset with $K^{trans}$ NIfTI files for each synthetic and clinical patient.
2. Modify the script as needed to provide the necessary paths to your $K^{trans}$ NIfTI files.
    - Recommended to save as inside a 'TeamName' directory and seperate 'TeamName_neutral' directory for reproduced entry values.
    - See example directories in the correct formatting in the constant $K^{trans}$ [test case](Scoring/entryDirectories).
3. Run the script using a compatible Python environment to obtain the evaluation scores for your method.

Please note that the scoring script assumes the availability of the required data files and follows specific guidelines outlined

## Python Package Requirements
- numpy 1.20.1
- matplotlib 3.3.4
- pydicom 2.1.2
- scipy 1.6.1
- nibabel 3.2.1

## Scoring Script Usage

A scoring script is provided in this repository: [challengeScoring.py](Scoring/challengeScoring.py). This script is designed to evaluate the performance of different $K^{trans}$ submissions for the challenge data. However, please note that using this script requires a fully analyzed dataset consisting of $K^{trans}$ NIfTI files for each synthetic and clinical patient as seen in uploaded [test folder](Scoring/entryDirectories/constantKtransModel).

To use the scoring script, follow these steps:

1. Clone or pull this directory.
2. Analyse the Challenge data set.
    - Find the challenge data for download in [this repository](ChallengeDICOMData), or at the [OSF page](https://osf.io/u7a6f/files).
    - Follow the [Challenge guidelines](OSIPI_DCE_Challenge_Guidelines.pdf).
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

# DRO Production

The files contained in the [DRO_production](DRO_Production) directory were used to produce the synthetic data used in the challenge.

>**Warning**
>You must unzip the files in [AIFmasksforDRO](DRO_Production/AIFmasksforDRO) otherwise the script will not run correctly, unless you wish to define your own AIF regions.

>**Note**
>To select your own AIF region please adjust the `use_saved` parameter to `0` on `Line 93` of [main.py](DRO_Production/main.py).

To use:
- The overall script file is [`main.py`](DRO_Production/main.py).
- Make sure all relevant files are unzipped (see notes above).
- The script has to be run seperately for synthetic patient 1 and 2.
    * Change `DROnum` to `1` or `2` between runs.
- If you wish to output new DRO data change `save_new_DRO_data` to `1` on `Line 45`.
- If you have selected to chose you own AIF a pop up window will appear.
    * You need to trace a full region outline or it will auto connect the end two points with a stright line.
    * You may select more than one region.
    * Close the pop up after you have selected all regions and the array data will be saved in the script.
- To apply this for other DICOM DCE-MRI data provide the correct paths in `Lines 29-40`.


