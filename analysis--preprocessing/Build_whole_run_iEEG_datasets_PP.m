%% Build whole run EEGLAB dataset:

clear all
close all
clc;
% set the relevant path:
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;


% add/remove paths:
addpath(fullfile(path_to_toolboxes,'eeglab2021.1'));
rmpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));

[ALLEEG, EEG, CURRENTSET] = eeglab;
warning('off')
subjects={'PP01','PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP14','PP15','PP16','PP17','PP18'};

for subjid = subjects(end)
    subjid = cell2mat(subjid);
    close all;
    clearex('subjects', 'subjid', 'ALLEEG', 'EEG');
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    maindir=fullfile(parentfolder,subjid);
    cd(maindir)
    
    % Select EDF file:
    switch subjid
        case 'PP01'
            edffilename = fullfile(maindir,'iEEG','PP01_19.01.2017_222.edf');
        case 'PP02'
            edffilename = fullfile(maindir,'iEEG','PP02_14122016_1632to1730.edf');            
        case 'PP03'
            edffilename = fullfile(maindir,'iEEG','PP03_06.12.17_secondseeg.edf');
        case 'PP04'
            edffilename = fullfile(maindir,'iEEG','PP04_24.03.2017_2.edf');            
        case 'PP05'
            edffilename = fullfile(maindir,'iEEG','PP05_06.09.17.edf');            
        case 'PP06'
            edffilename = fullfile(maindir,'iEEG','PP06_07122016_AllTasks_PPGoNoGoAuditoryVisual.edf');
        case 'PP07'
            edffilename = fullfile(maindir,'iEEG','PP07_29.11.17.edf');            
        case 'PP08'
            edffilename = fullfile(maindir,'iEEG','PP08.edf');
        case 'PP09'
            edffilename = fullfile(maindir,'iEEG','PP09_01.02.17.edf');
        case 'PP10'
            edffilename = fullfile(maindir,'iEEG','PP10_01.12.2016.edf');
        case 'PP11'
            edffilename = fullfile(maindir,'iEEG','PP11_18.01.18.edf');
        case 'PP13'
            edffilename = fullfile(maindir,'iEEG','PP13_19.03.18.edf');
        case 'PP14'
            edffilename = fullfile(maindir,'iEEG','PP14_04.05.17.edf');
        case 'PP15'
            edffilename = fullfile(maindir,'iEEG','PP15_21.06.18.edf');
        case 'PP16'
            edffilename = fullfile(maindir,'iEEG','PP16_03.01.18.edf');
        case 'PP17'
            edffilename = fullfile(maindir,'iEEG','PP17_15.02.18.edf');
        case 'PP18'
            edffilename = fullfile(maindir,'iEEG','PP18_23.02.18.edf');
                         
    end
    
    outdir = fullfile(maindir,'EEGLAB_datasets','raw'); % ADJUST OUTDIR
    mkdir(outdir);
    setname =[subjid '_pink_panther_entire_run'];
    
    % Load EDF:
    if contains(edffilename,'edf')
        
            EEG = pop_biosig(edffilename,'importevent','on');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', setname);
            eeglab redraw;
            datasetFileName = [EEG.setname '.set'];
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);            
          eeglab redraw;
    elseif contains(edffilename,'set')
        [EEG] = pop_loadset('filename', edffilename);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
        EEG.setname = setname;
        datasetFileName = [EEG.setname '.set'];
        eeglab redraw; 
    end
    trigchannel = zeros(EEG.pnts,1);
    trigchannel([EEG.event.latency]) = 1;
    
    EEG.data(end+1,:) = trigchannel;
    EEG.nbchan = size(EEG.data,1);
    EEG.chanlocs(end+1).labels = 'TRIG';
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = eeg_checkset( EEG );
    % Save set
    EEG = pop_saveset( EEG,  'filename', datasetFileName, 'filepath', outdir);
    disp('data saved')
    
    inFileName = datasetFileName;
end

