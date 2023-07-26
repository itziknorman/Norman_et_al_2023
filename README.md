# Data and code accompanying Norman_et_al_2023

﻿This repository contains intra-cranial EEG recordings and analysis code related to the paper:
Norman et al., "Hippocampal ripples support recognition memory in children".

The iEEG data underwent standard preprocessing, as detailed in the paper, and were stored in EEGLAB datasets within each patient's folder (PP02, PP03, PP04, etc.), under the subfolders "EEGLAB_datasets" and "EEGLAB_datasets_BP" (BP denotes bipolar recordings). 

Additional data structures containing spectrograms, ripples time-stamps, etc., that are utilized in the group-level analyses can be found in the "results/data" folder (stored as .mat files). 

The patients' anatomical data underwent standard preprocessing in Freesurfer [1] and iELVis [2], including the recon-all procedure and hippocampal subfields parcellation. This data is located in the "Freesurfer" folder. 
The main analysis code was written in MATLAB (R2021a). All functions and scripts are organized in subdirectories within the "matlab_scripts" folder, based on their respective subjects of analysis, such as ripple detection, ripple rate and PETH, deconvolution of peri-ripple responses, and more.
The analysis was conducted on a desktop computer equipped with a Intel® Xeon(R) CPU E5-2630 v3 @ 2.40GHz × 16 and 32GB RAM. Matlab's Signal Processing Toolbox is essential.

General instructions:

Download the following zip files -
“Norman_et_al_2023_iEEG_data.zip" 
"Norman_et_al_2023_code.zip"
“Norman_et_al_2023_Freesurfer.zip”

Unzip the content of these files under the same parent directory (e.g. "Norman_et_al_zenodo").
Note that the analysis code assumes that the content of these three zip files is located within a single parent folder.

Additional remarks:

1) The dataset is anonymized. 

2) Before running the analyses, set the correct paths in the "defineParentFolderAndPath.m" located in the folder "matlab_scripts".

3) The "external_matlab_toolboxes" folder encompasses various additional open-source tools developed by others, including EEGLAB, FieldTrip, and others.

Main toolboxes used in the analysis:
•	EEGLAB, version: v2021.1
•	Superbar / Scott Lowe (https://github.com/scottclowe/superbar), version: 1.5.0
•	Export_fig / Oliver J. Woodford & Yair M. Altman (https://github.com/altmany/export_fig)
•	iELVis toolbox / David Groppe (https://github.com/iELVis/iELVis) [2]
•	AFNI Matlab Library / Ziad Saad & Gang Chen (https://sscc.nimh.nih.gov/afni/matlab/)
•	Cbrewer / Charles Robert (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)
•	Mass Univariate ERP Toolbox (https://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox) [4]

5) The majority of the analyses are organized as lengthy scripts divided into distinct sections. The recommended approach for executing the code is to progress through each section one at a time, in sequential order.

6) The following table links the code to the main figures in the paper:

FIGURE #		MATLAB FILE
1f			analysis___ripples_rate_PETH_analysis_final.m
2b			get_ripple_events.m
2c-e 			analysis___draw_ripple_spectrogram_group_result.m
2f-g			analysis___ripples_rate_accross_conditions_and_IRI_distribution
3a-c,e,f		analysis___draw_ripple_spectrogram_group_result.m
3d			analysis___ripple_rate_analysis_old_vs_new_anatomy.m
4			analysis___draw_clip_spectrogram_group_result.m

References
1.	Fischl, B. FreeSurfer. Neuroimage 62, 774–781 (2012).
2.	Groppe, D. M. et al. iELVis: An open source MATLAB toolbox for localizing and visualizing human intracranial electrode data. J. Neurosci. Methods 281, 40–48 (2017).
3.	Delorme, A. & Makeig, S. EEGLAB: An open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. J. Neurosci. Methods 134, 9–21 (2004).
4.	Groppe, D. M., Urbach, T. P. & Kutas, M. Mass univariate analysis of event-related brain potentials/fields I: A critical tutorial review. Psychophysiology 48, 1711–1725 (2011).
