clear all
close all
clc;
% set the relevant path:
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;

% add/remove paths:
addpath(fullfile(path_to_toolboxes,'eeglab2021.1'));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));
rmpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));

[ALLEEG, EEG, CURRENTSET] = eeglab;

subjects={'PP01','PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};

for subjid = subjects(1:end)
    subjid = cell2mat(subjid);
    %%
    close all;
    clearex('subjects', 'subjid', 'ALLEEG', 'EEG', 'parentfolder', 'path_to_toolboxes');
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    
    maindir=fullfile(parentfolder,subjid);
    outdir=fullfile(maindir,'EEGLAB_datasets','preprocessed'); % ADJUST OUTDIR
    inFileName = fullfile(maindir,'EEGLAB_datasets','raw',[subjid '_pink_panther_entire_run.set']);
    mkdir(outdir);
    
    outFileName =[subjid '_pink_panther_preprocessed.set'];
    
    % Load raw EEGLAB dataset:
    [EEG] = pop_loadset('filename', inFileName);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
    eeglab redraw
    
    % Select XLS electrodes file:
    path = fullfile(maindir,'BioImage');
    [channel_location_file]=dir(fullfile(path,'*TAL*'));
    channel_location_file = fullfile(path,channel_location_file.name);

    
    if exist('run_commands','var')
        for k=1:numel(run_commands)
            eval(run_commands{k});
        end
    end
    
    %% set AUX Channels:
    
    auxsearch = {'MKR2','ECG','EKG','EMD','EMG'};
    
    for k = auxsearch
        ind = find(contains({EEG.chanlocs.labels},k,'IgnoreCase',true));
        if isempty(ind), continue; 
        else
            for i = 1:numel(ind)
                if i==1
                    EEG=pop_chanedit(EEG,'changefield',{ind(i) 'labels' cell2mat(k)});   
                elseif i>1
                    EEG=pop_chanedit(EEG,'changefield',{ind(i) 'labels' [cell2mat(k) '-' num2str(i)]});                   
                end
                fprintf('\n Renaming channel %d: %s \n',ind(i),EEG.chanlocs(ind(i)).labels);
            end
        end
    end
        
    if sum(multiStrFind({EEG.chanlocs.labels},'TRIG'))~=1
        fprintf('\nChannels: %d\n', find(contains({EEG.chanlocs.labels},'TRIG')))       
        error('******** There is a problem with the TRIG channel! ********');
    end
    
    aux_channels=[find(contains({EEG.chanlocs.labels},'EOG','IgnoreCase',true)),...      
                  find(contains({EEG.chanlocs.labels},'ECG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'EMG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'TRIG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'EMD','IgnoreCase',true))];
              

    %% Remove non-channels:
      
    if contains(subjid,'ASSUTA')
        load(fullfile(maindir,'BioImage','SUB01ASSUTA_matlab.mat'),'bio_order');
        ch_to_remove=setdiff(find(~conEEGmtains({EEG.chanlocs.labels},bio_order)),aux_channels);
    else
        ch_to_remove=setdiff(find(~multiStrFind({EEG.chanlocs.labels},'EEG')),aux_channels);
    end
    
    for ch={EEG.chanlocs(ch_to_remove).labels}
        fprintf('--- Removing Channel: %s --- \n',cell2mat(ch));
    end
    
    EEG = pop_select(EEG,'nochannel',ch_to_remove);
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;
    
    % Rename EEG channels:
    idx=find(multiStrFind({EEG.chanlocs.labels},'EEG'));
    for ch=idx
        new_label=EEG.chanlocs(ch).labels(5:end);
        fprintf('--- Relabeling Channel: %s --- \n',new_label);
        EEG=pop_chanedit(EEG,'changefield',{ch 'labels' new_label});
    end
    aux_channels=[find(contains({EEG.chanlocs.labels},'EOG','IgnoreCase',true)),...      
                  find(contains({EEG.chanlocs.labels},'ECG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'EMG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'TRIG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'EMD','IgnoreCase',true))];

    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);
    
    % remove extra time after the end of the experiment:
    Define_experiment_end_PP;
    EEG = pop_select(EEG,'nopoint',experiment_end:EEG.pnts);
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);
    % Resample to 500 Hz: (optional)
    % EEG = pop_resample(EEG,500);
    
    % Remove DC from each electrode:
    EEG = pop_rmbase(EEG,[EEG.times(1) EEG.times(end)]);
    EEG.setname=[outFileName(1:end-4) ' - resampled - DC removed'];
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw;
    
    %% Electrodes localization:
    EEG=ReadElectrodeCoord(EEG,channel_location_file,maindir);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;
    
    % HeadPLOT setup:
    if ~exist(fullfile(maindir,[subjid '_spline_file_MNIbrain.spl']),'file')
        mesh=fullfile(parentfolder,'MNI_brain_mesh','MNI_brain_downsampled.mat');
        headplot('setup', EEG.chanlocs, fullfile(maindir,[subjid '_spline_file_MNIbrain.spl']),  'orilocs','on','meshfile',mesh, 'transform',[0 0 0 0 0 -pi/2 1 1 1]);
    end
    % Plot Brain to verify:
%     figure; hold on;
%     mesh=fullfile(parentfolder,'MNI_brain_mesh','MNI_brain_downsampled.mat');
   
%     headplot_itzik(EEG.data,fullfile(maindir,[subjid '_spline_file_MNIbrain.spl']),[],'meshfile',mesh,'electrodes','on', ...
%         'title',subjid,'labels',1,'cbar',0, 'maplimits','absmax','colormap',colormap('Gray'));
%     alpha(0.15)
%
    %% Preprocessing:
    good_channels=setdiff(1:EEG.nbchan,aux_channels);
    figure; spectopo(EEG.data(good_channels,:),0,EEG.srate,'percent',5);
    title('Before Removing Line Noise')
    
    % Remove line noise using the new EEGLAB FIR filter:
    if contains(subjid,'ASSUTA'), notchFreqs=[50:50:450];
    else,  notchFreqs=[50]; end
    filterWidth=1.5; % Hz
    EEG_clean=EEG;
    for f=notchFreqs
        % Adjust the filter order manually! (use the EEGLAB menu to calculate the order)
        [EEG_clean,~,b] = pop_firws(EEG_clean, 'fcutoff', [f-filterWidth f+filterWidth], 'ftype', 'bandstop', 'wtype', 'hamming', 'forder', 1100); % 550 -> 3HZ transition
        figure; freqz(b);
    end
    
    figure; spectopo(EEG_clean.data(good_channels,:),0,EEG.srate,'percent',5);
    title('After Removing Line Noise')
    eegplot(EEG_clean.data(good_channels,:),'color','off','srate',EEG.srate,'winlength',15,'limits',[EEG.times(1) EEG.times(end)])
    
    ch = find(strcmpi({EEG.chanlocs.labels},'RAWTRIG'));
    EEG_clean.data(ch,:) = EEG.data(ch,:); clear ch; % keep the unfiltered trigchannel
    ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
    EEG_clean.data(ch,:) = EEG.data(ch,:); clear ch; % keep the unfiltered trigchannel
    % Store DATA:
    EEG=EEG_clean;
    EEG.setname=[outFileName(1:end-4) ' - Filtered'];
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw;
    
    %% Find noisy channels in 2 steps:
    % (first automatically and then manualy by visual inspection)
    clear noisy_channels
    % automatically exclude any channel with stats above 5 SD:
    excluded_channels=DetectNoisyChannelsForCREF(EEG,aux_channels,1,99,5); 
    good_channels=setdiff(1:EEG.nbchan,excluded_channels);
      
    %% ===================================================================
    % Compute the robust mean reference (PJ Rousseeuw, AM Leroy. Robust Regression and Outlier Detection, 1987)
    signal_channels=setdiff(1:EEG.nbchan,aux_channels);
    mean_good_signal=robustMean(EEG.data(good_channels,:),1,5);
    % ===================================================================
    
    % EMG-from-high-frequency:
    close all;
    disp('Processing EMGHF channel...')
    forder_HP=330;
    HF = eegfilt(mean_good_signal,EEG.srate,100,[],0,forder_HP,0,'fir1',0);
    EMGHF = abs(hilbert(HF));
    ind = find(strcmpi({EEG.chanlocs.labels},'EMGHF'));
    if isempty(ind)
        EEG.data(end+1,:) = EMGHF;
        EEG.nbchan = size(EEG.data,1);
        EEG.chanlocs(end+1).labels = 'EMGHF';
    else
        EEG.data(ind,:) = EMGHF;
    end
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = eeg_checkset( EEG );
    aux_channels=union(aux_channels,find(strcmpi({EEG.chanlocs.labels},'EMGHF')));
    
    
    % Substructing:
    for channel=1:EEG.nbchan
        if  ismember(channel,aux_channels)
            disp(['Skip ch. ' num2str(channel)]);
        else
            disp('*');
            EEG.data(channel,:)= EEG.data(channel,:) - mean_good_signal;
        end
    end
    EEG.setname=[outFileName(1:end-4) ' - Common Ref Subtructed'];
    EEG.data(end+1,:)=mean_good_signal;
    EEG.nbchan = size(EEG.data,1);
    EEG.chanlocs(end+1).labels = 'CREF';
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = eeg_checkset( EEG );
    aux_channels=union(aux_channels,find(strcmpi({EEG.chanlocs.labels},'CREF')));
    
    
    %% Remove DC from each electrode:
    EEG = pop_rmbase(EEG,[EEG.times(1) EEG.times(end)]);
    EEG.setname=[outFileName(1:end-4) ' - DC removed (second time)'];
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw;
    %====================================================================
    eegplot(EEG.data(good_channels,:),'color','off','srate',EEG.srate,'winlength',30,'limits',[EEG.times(1) EEG.times(end)])
    figure; spectopo(EEG.data(good_channels,:),0,EEG.srate,'percent',10);
    title('After common ref.')
    
    %% Re-Check for Noisy Channels:
    
    % Adjust Tresholds Manually:
    % Inspect any channel with stats above 3 SD:
    excluded_channels=DetectNoisyChannelsForCREF(EEG,aux_channels,2,99,3);
    good_channels=setdiff(1:EEG.nbchan,excluded_channels);
    
    
    %% Save Channel List:
    save(fullfile(maindir,[subjid '_channels_list']),'excluded_channels','good_channels');
    disp('channels list was saved')
    
    
    %% Save set:
    EEG.setname=outFileName(1:end-4);
    EEG = pop_saveset( EEG,  'filename', outFileName, 'filepath', outdir);
    disp('data saved')
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    eegplot(EEG.data(good_channels,:),'color','off','srate',EEG.srate,'winlength',30,'limits',[EEG.times(1) EEG.times(end)])
    
end



