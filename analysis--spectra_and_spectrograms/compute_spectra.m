% The following script loads the raw iEEG data and computes multitaper
% spectrum per each channel, individually for each condition.
%
% Author: Yitzhak Norman, 2023

clear all
close all
clc;
% set the relevant path:
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;

% add/remove paths:
addpath(fullfile(path_to_toolboxes,'eeglab2021.1'));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));
[ALLEEG, EEG, CURRENTSET] = eeglab;

subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};

blocks = {'R1','R2','M1','M2','MT'};
DATA=struct;
chINFO={};

for iSub=1:length(subjects)
    
    clearvars -except subjects blocks iSub outdir path_to_toolboxes parentfolder DATA chINFO
    subjid = subjects{iSub};
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    
    % load datasets:
    maindir=fullfile(parentfolder,subjid);
    datasetID='pink_panther';
    
    % Set Manually:
    % 1 = Common Ref;
    % 2 = Bipolar montage;
    
    ref_flag = 1;
    switch ref_flag
        case 1, datadir=fullfile(maindir,'EEGLAB_datasets');
        case 2, datadir=fullfile(maindir,'EEGLAB_datasets_BP');
    end
    
    % LOAD DATASETS: Full time courses
    for iBlock = 1:length(blocks)
        blockID = blocks{iBlock};
        filename = dir(fullfile(datadir,sprintf('%s_%s_preproc*%s.set',subjid,datasetID,blockID)));
        filename = fullfile(filename.folder,filename.name);
        [EEG] = pop_loadset('filename',filename);
        [ALLEEG EEG ~] = eeg_store(ALLEEG, EEG, iBlock);
    end
    
    % Set outdir:
    outdir=fullfile(parentfolder,'results','data','multitaper-spectra',['spectra_data_ref_' num2str(ref_flag)]);
    if ~exist(outdir,'dir')
        mkdir(outdir);
        disp('Creating Output Directory...')
    end
    
    
    %% ========================================
    % add events and epoch and merge:
    
    for iBlock = 1:numel(ALLEEG)
        
        EEG = ALLEEG(iBlock);
        blockID = blocks{iBlock};
        
        if EEG.srate>256, EEG = pop_resample(EEG,256); end
        
        event_list = {};
        e = 1;
        epochlength = 10000; % in ms
        timestamps = (EEG.times(1)+epochlength):epochlength:(EEG.times(end)-epochlength);
        
        % add ripple events to the list:
        for index = 1 : length(timestamps)
            current_event_timing = timestamps(index)./1000;
            event_list{e,1} = sprintf('%s',blockID);
            event_list{e,2} = current_event_timing;
            event_list{e,3} = blockID;
            event_list{e,4} = sprintf('E%d-%s',index,blockID);
            e = e+1;
        end
        
        EEG = pop_importevent( EEG, 'event', event_list, 'fields',{'type', 'latency','blockID','code'}, 'timeunit',1, 'append','yes');
        EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
        [ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG, iBlock);
        fprintf('\n%s - stored OK\n',EEG.setname);
        
        % Epoch data:
        event_ind=find(startsWith({EEG.event.type},{'R1','M1','R2','M2','MT'}));
        epochlim = [-epochlength/2 epochlength/2]./1000; % in sec
        EEG = pop_epoch( EEG,{},epochlim,'eventindices',event_ind);
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, iBlock);
        
    end
    
    if numel(ALLEEG)>1, EEG_mereged = pop_mergeset( ALLEEG, [1:numel(ALLEEG)], 0);    % Merge all blocks together
    else, EEG_mereged = EEG; end
    
    % get the labels of the epochs:
    epoch_list = {};
    mainEventInd = cell2mat(cellfun(@(x)find(cell2mat(x)==0),{EEG_mereged.epoch.eventlatency},'UniformOutput',0));
    for i = 1:EEG_mereged.trials
        main_event_ind = find(cell2mat([EEG_mereged.epoch(i).eventlatency]) == 0);
        epoch_list{i,1} = EEG_mereged.epoch(i).eventtype{main_event_ind};
    end
    
    %% Spectrograms parameters:
    % close all
    addpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')))
    
    params=struct;
    params.Fs=EEG_mereged.srate; % Sampling Rate [Hz]
    % Tapers Configuration:
    T=4; % Window Size [Sec]
    W=1; % Half Band Width (Frequency Resolution) [Hz]
    params.tapers=[(T * W), 7]; % Must be less than (2*T*W)-1
    params.fpass=[0 125]; % [fmin fmax] [HZ]
    params.err=0; %[2 0.05]; % Jackknife error bars [p=0.05]
    params.trialave=0; % Avrage over trials
    params.pad=0;
    
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0);
    ch_ind = find(ismember({EEG_mereged.chanlocs.labels},hippocampus_all_channels));
    
    %% RUN ANALYSIS AND PLOT:
    set_figure_colors;
    close all;
    saveflag = 1;
    counter = 0;
    S0 = []; S1 = []; S2 = []; S3 = []; S4 = [];
    S1_norm = []; S2_norm = []; S3_norm = []; S4_norm = [];
    ch_list = {};
    for cnum = ch_ind
        
        channel = EEG_mereged.chanlocs(cnum).labels;
        counter = counter+1;
        
        fprintf('\n *** Subject: %s, Channel: %s ***',subjid,channel)
        
        data=squeeze(EEG_mereged.data(cnum,:,:));
        data=bsxfun(@minus,data,nanmean(data)); % demean
        [spec0,f_out] = mtspectrumc(data(:,contains(epoch_list,'R1')), params ); % REST 1
        [spec1,f_out] = mtspectrumc(data(:,contains(epoch_list,'M1')), params ); % MOVIE 1
        [spec2,f_out] = mtspectrumc(data(:,contains(epoch_list,'R2')), params ); % REST 2
        [spec3,f_out] = mtspectrumc(data(:,contains(epoch_list,'M2')), params ); % MOVIE 2
        [spec4,f_out] = mtspectrumc(data(:,contains(epoch_list,'MT')), params ); % MEMORY TEST
        
        
        spec0 = spec0.*conj(spec0); % power
        spec1 = spec1.*conj(spec1); % power
        spec2 = spec2.*conj(spec2); % power
        spec3 = spec3.*conj(spec3); % power
        spec4 = spec4.*conj(spec4); % power
        
        
        % concatenate:
        S0(:,counter)=mean(spec0,2);
        S1(:,counter)=mean(spec1,2);
        S2(:,counter)=mean(spec2,2);
        S3(:,counter)=mean(spec3,2);
        S4(:,counter)=mean(spec4,2);
        
        % Normalizes each epoch by the baseline condition:
        S0_norm(:,counter)=mean(bsxfun(@rdivide,spec0,mean(spec0,2)),2);
        S1_norm(:,counter)=mean(bsxfun(@rdivide,spec1,mean(spec0,2)),2);
        S2_norm(:,counter)=mean(bsxfun(@rdivide,spec2,mean(spec0,2)),2);
        S3_norm(:,counter)=mean(bsxfun(@rdivide,spec3,mean(spec0,2)),2);
        S4_norm(:,counter)=mean(bsxfun(@rdivide,spec4,mean(spec0,2)),2);
        ch_list{counter} = strcat(subjid,'_',channel);
        
        
        figure('color','w','name','raw spectra','position',[0 0 200 200]);  hold on;
        plot(f_out,10*log10(S0(:,counter)),'color',COLOR.black);
        plot(f_out,10*log10(S1(:,counter)),'color',COLOR.lightred);
        plot(f_out,10*log10(S2(:,counter)),'color',COLOR.gray);
        plot(f_out,10*log10(S3(:,counter)),'color',COLOR.bordo);
        plot(f_out,10*log10(S4(:,counter)),'color',COLOR.green);
        xlabel('Frequency'); ylabel('Power [10*log10(\muV)]')
        title(channel)
        xlim([0 20])
        %set(gca,'XScale','log')
        figure('color','w','name','normalized spectra','position',[0 0 200 200]); hold on;
        plot(f_out,10*log10(S1_norm(:,counter)),'color',COLOR.orange);
        plot(f_out,10*log10(S2_norm(:,counter)),'color',COLOR.blue);
        plot(f_out,10*log10(S3_norm(:,counter)),'color',COLOR.red);
        plot(f_out,10*log10(S4_norm(:,counter)),'color',COLOR.green);
        xlim([0 20])
        xlabel('Frequency'); ylabel('Relative power change (dB)')
        title(channel)
        %set(gca,'XScale','log')
        drawnow;
    end
    
    
    if saveflag
        save(fullfile(outdir,'spectra_parameters'),'params','f_out')
        save(fullfile(outdir,[subjid '_electrode_details']),'ch_list');
        save(fullfile(outdir,[subjid '_spectra_data']),'S0','S1','S2','S3','S4','S0_norm','S1_norm','S2_norm','S3_norm','S4_norm','-v7.3')
        fprintf('\n ---> DATA SAVED: %s',subjid);
    end
end






