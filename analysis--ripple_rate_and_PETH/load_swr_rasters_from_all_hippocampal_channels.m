  
      
%% load data:

DATA = struct;   
% Set data structure:
FN={'Rall','R1','R2','M1','M2','MT'};
for j = 1:length(FN)
    DATA.(FN{j}).subjid = [];
    DATA.(FN{j}).taskid = [];
    DATA.(FN{j}).channelid = [];
    DATA.(FN{j}).stimulus_type = [];
    DATA.(FN{j}).stimulus = [];
    DATA.(FN{j}).RT = [];
    DATA.(FN{j}).response = [];
    DATA.(FN{j}).correct = [];
    DATA.(FN{j}).raster = [];
end
DATA.rippletable = [];
DATA.trialtime = [];

for iSub = 1:length(subjects)
    subjid = subjects{iSub};
    % load data (ripple rasters):
    ripplesdir = fullfile(parentfolder,'results','data','Ripple_times_hamming_2-4std_adjusted_band_20ms_30ms',subjid);
    filename = fullfile(ripplesdir,sprintf('%s_ripple_psth_data_ref_%d_%s_all_hippocampal_channels.mat',subjid,ref_flag,time_locking_event));
    tmp = load(filename);
    DATA = CatStructFields(DATA,tmp.DATA,1);
    fprintf('\nSubject %d - o.k.\n',iSub)
end
DATA.trialtime = DATA.trialtime(1,:);
DATA.Fs = 1000;
Fs = DATA.Fs;
% =========================================================================
