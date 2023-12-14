# Code accompanying Norman_et_al_2023

﻿This repository contains code for the main analyses described in Norman et al. (2023), "Hippocampal ripples linked to recognition memory judgments in children".

For the intracranial data and the complete code, please refer to Zenodo DOI (will be made available after publication): 10.5281/zenodo.8187953.

**General remarks:**
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
Figure #	  Matlab file
1f	  analysis___ripples_rate _PETH_analysis_final.m

2b	  get_ripple_events.m

2c-e	  analysis___draw_ripple_spectrogram_group_results.m

2f-g	  analysis___ripples_rate _across_conditions_and_IRI_distribution.m

3a-c,e,f	  analysis___draw_ripple_spectrogram_group_results.m

3d	  analysis___ripple_rate_analysis_old_vs_new_anatomy.m

4a	  analysis___ripples_rate_PETH_analysis_final.m

4b-f	  analysis___ripple_train_crosscorrelations.m

5	  analysis___draw_ripple_spectrogram_group_results.m

![image](https://github.com/itziknorman/Norman_et_al_2023/assets/59057794/6ad19145-59c3-40de-88af-db6c425a960b)





