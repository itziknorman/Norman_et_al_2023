% The following script loads the raw iEEG data and computes clip-triggered wavelet
% spectrograms.
% Author: Yitzhak Norman, 2022

clear all
close all
clc;
% set the relevant path:
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;

% add/remove paths:
addpath(fullfile(path_to_toolboxes,'eeglab2021.1'));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));
cd(fullfile(parentfolder,'matlab_scripts'));
[ALLEEG, EEG, CURRENTSET] = eeglab;
% =================================
set(0,'DefaultFigureVisible','off')
% =================================
subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};
blocks = {'R1','M1','R2','M2','MT'};


for iSub=1:length(subjects)
    
    clearvars -except subjects blocks iSub outdir path_to_toolboxes parentfolder
    subjid = subjects{iSub};
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    DATA=struct;
    chINFO={};
    % load datasets:    
    maindir=fullfile(parentfolder,subjid);
    datasetID='pink_panther';

    % Set Manually:
    % 1 = Common Ref;
    % 2 = Bipolar montage;
    
    ref_flag = 2;
    switch ref_flag
        case 1, datadir=fullfile(maindir,'EEGLAB_datasets'); refID = '';
        case 2, datadir=fullfile(maindir,'EEGLAB_datasets_BP'); refID = 'BP_montage_';
    end
    
    % Timing of original events in movie: (in sec)
    eventsInMovie = readtable(fullfile(parentfolder,'events.xlsx'));
    [origEventTime,sortind] = sort(eventsInMovie.timeInSec);
   
    
    % LOAD DATASETS: Full time courses
    for iBlock = 1:length(blocks)
        blockID = blocks{iBlock};
        filename = fullfile(datadir,sprintf('%s_%s_preprocessed_%s%s.set',subjid,datasetID,refID,blockID));
        EEG = pop_loadset('filename',filename);
        EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
        if ismember(blockID,{'M1','M2'})
            origEventName = cellstr(string(blockID) + eventsInMovie.name(sortind));
            strind = round(EEG.event(strcmpi({EEG.event.type},'str')).latency); % finds the start point
            margin_correction = EEG.times(strind)./1000;
            event_list = {}; 
            event_list(:,1) = origEventName; 
            event_list(:,2) = num2cell(origEventTime + margin_correction);
            EEG = pop_importevent( EEG, 'event', event_list, 'fields',{ 'type', 'latency'}, 'timeunit',1, 'append','yes');
            EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
            fprintf('\n%s events added sucessfully... O.K\n',blockID)          
        end
       [ALLEEG EEG ~] = eeg_store(ALLEEG, EEG,iBlock);
    end
    
    % Set outdir:
    outdir=fullfile(parentfolder,'results','data','clip_triggered_spectrograms',['ref_' num2str(ref_flag)]);
    if ~exist(outdir,'dir')
        mkdir(outdir);
        disp('Creating Output Directory...')
    end
        
    %% Spectrograms parameters:
    close all
    CM=colormap(jet); close;
    
    % Spectograms parameters:
    params = struct;
    params.freqs=[4 128];
    params.nfreqs=200;
    params.winsize=128;
    params.timesout=1000;
    params.cycles= [3 0.5];
    params.units='abs';
    params.freqscale='linear';
    set_figure_colors;
    
    %% RUN ANALYSIS AND PLOT:
    close all;    
    saveflag=1; % change to 1 if you want to save the figures
    norm_flag=0; % 0 = gain model; 1 = addative model;
    switch norm_flag
        case 0,  powunits='dB'; plim=10; cbar_units='dB'; file_end='wavelet_dB';
        case 1,  powunits='SD from BL';  plim=5; cbar_units='\sigma'; file_end='wavelet_zscore';
    end
    
    % ADJUST FIGDIR and OUTDIR:
    figdir=fullfile(parentfolder,'results','clip_triggered_spectrograms',subjid,['ref_' num2str(ref_flag)]);
    if ~exist(figdir,'dir')
        mkdir(figdir);
        disp('Creating Output Directory...')
    end    
    
    %% MAIN ANALYSIS LOOP:
 
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0);
    
    hipch = find(ismember({EEG.chanlocs.labels},hippocampus_all_channels));
    cref = find(strcmpi({EEG.chanlocs.labels},'CREF'));       
    EEGm = pop_mergeset( ALLEEG, [1:numel(ALLEEG)], 0);

    ET  =  [EEGm.event.latency];   
    event_ind1=find(startsWith({EEGm.event.type},'ClipA')); % old clips
    event_ind2=find(startsWith({EEGm.event.type},'ClipB')); % new clips
    event_ind3=find(startsWith({EEGm.event.type},'M1ClipA')); % orig events video presentation 1
    event_ind4=find(startsWith({EEGm.event.type},'M2ClipA')); % orig events video presentation 2

    %eegplot(EEGm.data,'color','off','srate',EEG.srate,'winlength',30,'limits',[EEG.times(1) EEG.times(end)],'events',EEGm.event)
 
    % Epoch data around memory test clips and orig movie events:
    fullepochlim = [-2 8];
    ix = [event_ind1,event_ind2 event_ind3 event_ind4];
    [EEGe,indices] = pop_epoch( EEGm,{},fullepochlim,'eventindices',ix);
 
    EE = cellfun(@(x,y)x{cell2mat(y)==0},{EEGe.epoch.eventtype},{EEGe.epoch.eventlatency},'UniformOutput',0);
    RESP = cellfun(@(x,y)x(ismember(x,{'old','new'}) & cell2mat(y)>0 & cell2mat(y)==max(cell2mat(y))),{EEGe.epoch.eventtype},{EEGe.epoch.eventlatency},'UniformOutput',0); 
    tmp = cellfun(@(x)isempty(x),RESP,'UniformOutput',0); RESP(cell2mat(tmp)) = {{'n/a'}}; % replace empty values with 0
    RESP = string(RESP); 
    
    old_responses = strcmp(RESP,'old'); 
    new_responses = strcmp(RESP,'new');
    
   
    trial_ind1 = startsWith(EE,'ClipA')&~new_responses;    % OLD, discarding wrong trials          
    trial_ind2 = startsWith(EE,'ClipB')&~old_responses;    % NEw, discarding wrong trials
    trial_ind3 = startsWith(EE,'M1ClipA'); % epochs aligned to original movie events
    trial_ind4 = startsWith(EE,'M2ClipA'); % epochs aligned to original movie events
    

    counter=0;  
    for cnum=[hipch]
        %%
        clc;
        channel=EEG.chanlocs(cnum).labels;
        fprintf('\nAnalizing channel: %s\n',channel)
    
        T = EEGe.times;
        epochlim = [T(1) T(end)];        
        
        [~,~,~,t_out,f_out,~,~,S] = newtimef(EEGe.data(cnum,:,:),numel(T),epochlim, EEGe.srate, 'cycles', params.cycles, 'freqs', params.freqs, 'winsize', params.winsize,...
            'scale',params.units,'plotitc','off','plotersp','off','timesout',params.timesout,'nfreqs',params.nfreqs,'baseline',nan,'freqscale',params.freqscale);
   
        t_out_sec = t_out./1000;
        S = S.*conj(S);
        S1 = S(:,:,trial_ind1); 
        S2 = S(:,:,trial_ind2);               
        S3 = S(:,:,trial_ind3);  
        S4 = S(:,:,trial_ind4);  
        
        S1_norm = nan(size(S1));
        S2_norm = nan(size(S2));
        S3_norm = nan(size(S3));
        S4_norm = nan(size(S4));
        baseline_ind = t_out<0;
        % Normalizes each condition by its baseline power:
        switch norm_flag            
            case 0
                disp('*** Normalization: dB ***')                
                S1_norm = bsxfun(@rdivide,S1, nangeomean(nangeomean(S(:,baseline_ind,trial_ind1|trial_ind2),2),3));
                S2_norm = bsxfun(@rdivide,S2, nangeomean(nangeomean(S(:,baseline_ind,trial_ind1|trial_ind2),2),3));      
                S3_norm = bsxfun(@rdivide,S3, nangeomean(nangeomean(S3(:,baseline_ind,:),2),3));      
                S4_norm = bsxfun(@rdivide,S4, nangeomean(nangeomean(S4(:,baseline_ind,:),2),3));      
                
            case 1
                disp('*** Normalization: zscore ***')       
                % normalize MT clips:
                BL = S(:,baseline_ind,trial_ind1|trial_ind2);
                mstd = sqrt(nanmean(nanvar(BL,0,2),3));  mbase = nanmean(nanmean(BL,2),3);                
                S1_norm = bsxfun(@rdivide,bsxfun(@minus,S1,mbase),mstd); % zscore normalization                
                S2_norm = bsxfun(@rdivide,bsxfun(@minus,S2,mbase),mstd); % zscore normalization
                % normalize M1/M2 orig clips:
                S3_norm = bsxfun(@rdivide,bsxfun(@minus,S3,nanmean(nanmean(S3(:,baseline_ind,:),2),3)),...
                                  sqrt(nanmean(nanvar(S3(:,baseline_ind,:),0,2),3))); % zscore normalization
                S4_norm = bsxfun(@rdivide,bsxfun(@minus,S4,nanmean(nanmean(S4(:,baseline_ind,:),2),3)),...
                                  sqrt(nanmean(nanvar(S4(:,baseline_ind,:),0,2),3))); % zscore normalization
                clear mbase mstd BL  
        end
        %======================================================================
        
        switch norm_flag
            case 0, S1_norm=single(10*log10(S1_norm)); S2_norm=single(10*log10(S2_norm));  S3_norm=single(10*log10(S3_norm)); S4_norm=single(10*log10(S4_norm));
            case 1, S1_norm=single(S1_norm); S2_norm=single(S2_norm); S3_norm=single(S3_norm); S4_norm=single(S4_norm); 
        end
        
        % single trial baseline correction:
        S1_norm = bsxfun(@minus,S1_norm,nanmean(S1_norm(:,baseline_ind,:),2));
        S2_norm = bsxfun(@minus,S2_norm,nanmean(S2_norm(:,baseline_ind,:),2));
        S3_norm = bsxfun(@minus,S3_norm,nanmean(S3_norm(:,baseline_ind,:),2));
        S4_norm = bsxfun(@minus,S4_norm,nanmean(S4_norm(:,baseline_ind,:),2));
        epochlim=[T(1) T(end)];        
       
        % High-frequency broadband power (HFB):
        F=[60 160];
        HFB1=squeeze(nanmean(S1_norm((f_out>F(1) & f_out<F(2)),:,:),1));
        HFB2=squeeze(nanmean(S2_norm((f_out>F(1) & f_out<F(2)),:,:),1));     
        HFB3=squeeze(nanmean(S3_norm((f_out>F(1) & f_out<F(2)),:,:),1));
        HFB4=squeeze(nanmean(S4_norm((f_out>F(1) & f_out<F(2)),:,:),1));    
        
        %======================================================================
        %figure position:
        POS = [0 0 250 200];
        
        % ================
        % draw HFB response:
        H = figure('color','w','name',sprintf('HFB clip response electrode %s',channel),'position',POS); hold on;
        shadedErrorBar(t_out_sec,mean(HFB1,2),std(HFB1,[],2)./sqrt(size(HFB1,2)),{'color',COLOR.blue,'linewidth',1},0.5)
        shadedErrorBar(t_out_sec,mean(HFB2,2),std(HFB2,[],2)./sqrt(size(HFB2,2)),{'color',COLOR.red,'linewidth',1},0.5)
        shadedErrorBar(t_out_sec,mean(HFB3,2),std(HFB3,[],2)./sqrt(size(HFB3,2)),{'color',COLOR.greendark,'linewidth',1},0.5)
        shadedErrorBar(t_out_sec,mean(HFB4,2),std(HFB4,[],2)./sqrt(size(HFB4,2)),{'color',COLOR.turquoise,'linewidth',1},0.5)
        xlabel('Time (s)'); ylabel('Amplitude (dB)');
        if saveflag, save_current_figure(H,figdir,1); end
        % ================
        % set plotting parameters:  
        [f_out_log,f_log_ind] = convert_linfreq_to_logfreq(f_out,params.nfreqs);  
        
        POSCB = [0.875,0.375,0.03,0.3];
        POSAX = [0.15 0.15 0.8 0.8];
        params_plot = struct;       
        params_plot.CM = CM;
        params_plot.clim = [-10 10];
        params_plot.xlim = [-0.5 6];
        params_plot.logscaleflag = 1;
        params_plot.cbflag = 1;
        params_plot.normflag = norm_flag;
        params_plot.position = POSAX;
        params_plot.cb_position = POSCB;
        params_plot.fig_position = POS;
        params_plot.t_out = t_out_sec;
        params_plot.f_out = f_out;
        params_plot.f_out_log = f_out_log;
        params_plot.f_log_ind = f_log_ind;
        params_plot.stimulus_onset = 0;  % in sec
        params_plot.stimulus_dur = 4;  % in sec
        params_plot.ch_label = channel;
        params_plot_grplevel = params_plot; % keep for grp level analysis
        % Figure 1 - OLD:
        params_plot.title = sprintf('OLD [N=%d]',EEGe.trials);
        params_plot.name = sprintf('OLD clip spectrogram electrode %s',params_plot.ch_label);    
        [H1,h1,hcb1] = plot_spectrogram(S1_norm,params_plot);
        if saveflag, save_current_figure(H1,figdir,1); end
        %Figure 2 - NEW:
        params_plot.title = sprintf('NEW [N=%d]',EEGe.trials);
        params_plot.name = sprintf('NEW clip spectrogram electrode %s',params_plot.ch_label);   
        [H2,h2,hcb2] = plot_spectrogram(S2_norm,params_plot);
        if saveflag, save_current_figure(H2,figdir,1); end
        % Figure 3 - OLD-NEW:
        params_plot.title = sprintf('OLD-NEW [N=%d]',EEGe.trials);
        params_plot.name = sprintf('OLD-vs-NEW clip spectrogram electrode %s',params_plot.ch_label);   
        [H3,h3,hcb3] = plot_spectrogram(bsxfun(@minus,mean(S1_norm,3),mean(S2_norm,3)),params_plot);
        if saveflag, save_current_figure(H3,figdir,1); end
                %%      
        counter = counter+1;
        % DATA structure for group analysis:
        % ERP =============================================================
        DATA.ERP{counter}=squeeze(EEGe.data(cnum,:,:));
        DATA.HFB1{counter}=HFB1;
        DATA.HFB2{counter}=HFB2;
        DATA.HFB3{counter}=HFB3;
        DATA.HFB4{counter}=HFB4;
        DATA.normspec1{counter}=S1_norm;
        DATA.normspec2{counter}=S2_norm;        
        chINFO{counter,1}=channel;
        chINFO{counter,2}=subjid;
        chINFO{counter,3}=cnum;
        
        drawnow;        
    end
    
    if saveflag
        electrodes_table = array2table(chINFO,'VariableNames',{'channel','subjid','cnum'});
        save(fullfile(outdir,'clip_spectrogram_parameters'),'params','params_plot','f_out','t_out');
        save(fullfile(outdir,[subjid '_clip_spectrogram_DATA']),'DATA','-v7.3')
        save(fullfile(outdir,[subjid '_electrodes_table']),'electrodes_table')
        fprintf('\n ---> DATA SAVED: %s',subjid);
    end
end






