# Code accompanying Norman_et_al_2023

﻿This repository contains code for the main analyses described in Norman et al. (2023), "Hippocampal ripples linked to recognition memory judgments in children".

For the intracranial data and the complete code, please refer to Zenodo DOI (will be made available after publication): 10.5281/zenodo.8187953.

--> Minimal demo: 
A minimal demo for ripple detection can be downloaded from here: 
https://www.dropbox.com/scl/fi/m2nje3e23mrzok81zmfqj/Ripple_detection_demo.zip?rlkey=ra0nudvy9khrj9t3dndz6whbe&dl=0

Unzip the folder, and once it's extracted, run the "get_ripple_events_demo.m" script. Ensure you have EEGLAB installed and set the path as specified in the MATLAB code before running the script.

**General remarks and instructions:**
The analysis code was written in MATLAB (R2021a). All functions and scripts are organized in subdirectories within the "matlab_scripts" folder, based on their respective subjects of analysis, such as ripple detection, ripple rate and PETH, ripple-triggered spectrograms, and more.
The analysis was conducted on a desktop computer equipped with a Intel® Xeon(R) CPU E5-2630 v3 @ 2.40GHz × 16 and 32GB RAM. Matlab's Signal Processing Toolbox is essential.
Most of the analyses are organized as lengthy scripts divided into distinct sections. The recommended approach for executing the code is to progress through each section one at a time, in sequential order.

External open-source toolboxes used in the analysis:
•	EEGLAB, version: v2021.1
•	Superbar / Scott Lowe (https://github.com/scottclowe/superbar), version: 1.5.0
•	Export_fig / Oliver J. Woodford & Yair M. Altman (https://github.com/altmany/export_fig)
•	iELVis toolbox / David Groppe (https://github.com/iELVis/iELVis)
•	AFNI Matlab Library / Ziad Saad & Gang Chen (https://sscc.nimh.nih.gov/afni/matlab/)
•	Cbrewer / Charles Robert (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)
•	Mass Univariate ERP Toolbox (https://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox)

**The following table links the code to the main figures in the paper:**

![image](https://github.com/itziknorman/Norman_et_al_2023/assets/59057794/6ad19145-59c3-40de-88af-db6c425a960b)





