% The following script loads the raw iEEG data and computes example wavelet
% spectrograms, peri-ripple field potential, and addtional info.
% Author: Yitzhak Norman, 2022, Chang Lab

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
    time_locking_event = 'peak';
    % Set Manually:
    % 1 = Common Ref;
    % 2 = Bipolar montage;
    
    ref_flag = 2;
    switch ref_flag
        case 1, datadir=fullfile(maindir,'EEGLAB_datasets'); refID = '';
        case 2, datadir=fullfile(maindir,'EEGLAB_datasets_BP'); refID = 'BP_montage_';
    end
    
    % LOAD DATASETS: Full time courses
    for iBlock = 1:length(blocks)
        blockID = blocks{iBlock};
        filename = fullfile(datadir,sprintf('%s_%s_preprocessed_%s%s.set',subjid,datasetID,refID,blockID));
        EEG = pop_loadset('filename',filename);
        EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
        eval(sprintf('[ALLEEG EEG ~] = eeg_store(ALLEEG, EEG, %d)',iBlock));
    end
    
    % Set outdir:
    outdir=fullfile(parentfolder,'results','data','ripple_triggered_spectrograms',['ref_' num2str(ref_flag)]);
    if ~exist(outdir,'dir')
        mkdir(outdir);
        disp('Creating Output Directory...')
    end
    
    % load swr timings:
    EEG = ALLEEG(1);
    ripplesdir = fullfile(parentfolder,'results','data','Ripple_times_hamming_2-4std_adjusted_band_20ms_30ms');
    ch_list = {EEG.chanlocs.labels};
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0);
      
    
    % ========================================
    % LOAD and add ripple events:
    
    for iBlock = 1:numel(ALLEEG)
        
        blockID = blocks{iBlock};
        event_list = {};
        e = 1;
        EEG = ALLEEG(iBlock);
        strind = round(EEG.event(1).latency); % finds the start point
        margin_correction = EEG.times(strind)./1000; 
        % the timing of the ripples are relative to the beginning of the
        % experimental run, so you need to add the 5 s margins                           

        for curHipChannel = hippocampus_all_channels
            
            curHipChannel = cell2mat(curHipChannel);
            currentRippleMat = fullfile(ripplesdir,subjid,sprintf('%s %s ripples %s.mat',subjid,blockID,curHipChannel));
            if ~exist(currentRippleMat,'file'), continue; end
            rippleTimeStamps = load(currentRippleMat);
            
            normamp = nanzscore(rippleTimeStamps.ripples.amplitude);
            % add ripple events to the list:
            for index = 1 : length(rippleTimeStamps.ripples.(time_locking_event))
                current_event_timing = rippleTimeStamps.ripples.(time_locking_event)(index) + margin_correction;
                event_list{e,1} = sprintf('Ripple%sE%s-%s',blockID,num2str(index),curHipChannel);
                event_list{e,2} = current_event_timing;
                event_list{e,3} = blockID;
                event_list{e,4} = 'SWR';
                event_list{e,5} = normamp(index);
                e = e+1;
            end
            
        end
        EEG = pop_importevent( EEG, 'event', event_list, 'fields',{'type', 'latency','blockID','code','amplitude'}, 'timeunit',1, 'append','yes');
        EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
        [ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG, iBlock);
        fprintf('\n%s - stored OK\n',EEG.setname);
    end
    
    
    
    %% Spectrograms parameters:
    close all
    CM=colormap(trubo); close;
    
    % Spectograms parameters:
    params = struct;
    params.freqs=[4:220];
    params.nyquistfreq=floor(EEG.srate/2);
    params.nfreqs=200;
    params.winsize=128;
    params.timesout=300;
    params.cycles= [2 15];
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
    figdir=fullfile(parentfolder,'results','ripple_triggered_spectrograms',subjid,['ripple_spectrogram_' file_end '_ref_' num2str(ref_flag) '_' time_locking_event]);
    if ~exist(figdir,'dir')
        mkdir(figdir);
        disp('Creating Output Directory...')
    end
    
    
    %% MAIN ANALYSIS LOOP:
    hipch = find(ismember({EEG.chanlocs.labels},hippocampus_all_channels));
    cref = find(strcmpi({EEG.chanlocs.labels},'CREF'));   
    
    EEGm = pop_mergeset( ALLEEG, [1:numel(ALLEEG)], 0);
    E = {EEGm.event.type};       % event list
    ET  =  [EEGm.event.latency];   
    
    %eegplot(EEGm.data,'color','off','srate',EEG.srate,'winlength',15,'limits',[EEG.times(1) EEG.times(end)],'events',EEGm.event)
    counter = 0;
    for cnum=[hipch]
        %%
        clc;
        channel=EEG.chanlocs(cnum).labels;
        fprintf('\nAnalizing channel: %s\n',channel)
        
        % Ripple-type subdevisions:
        for k = 1:numel(E),if size(E{k},1)>1,E{k} = E{k}(1,:); end; end
        swr_event_ind =find(startsWith(E,'Ripple') & endsWith(E,channel));
        
        % Epoch data around swr events:
        fullepochlim = [-1.5 1.5];
        [EEGe,indices] = pop_epoch( EEGm,{},fullepochlim,'eventindices',swr_event_ind);
        swr_event_ind = swr_event_ind(indices); % updated timelocking-event-indices
        EE = E(swr_event_ind); % updated trial list
        
        % run through all ripple events and get their labels:        
        Rind = struct;
        Rind.blockid = {};
        for i = 1:numel(swr_event_ind)
            str = E{swr_event_ind(i)};
            Rind.blockid{i,1} = str(7:8);
        end
            
        % Ripples:
        subepoch_idx = [abs(EEGe.times)<900];
        T = EEGe.times(subepoch_idx);
        epochlim = [T(1) T(end)];
        
        params.timesout=300;
        [~,~,~,t_out,f_out,~,~,S] = newtimef(EEGe.data(cnum,subepoch_idx,:),numel(T),epochlim, EEGe.srate, 'cycles', params.cycles, 'freqs', params.freqs, ...
            'scale',params.units,'plotitc','off','plotersp','off','timesout',params.timesout,'nfreqs',params.nfreqs,'baseline',nan,'freqscale',params.freqscale);
        t_out_sec = t_out./1000;
        S = S.*conj(S); % power
        
        S(f_out>params.nyquistfreq,:,:) = nan; % remove frequencies above nyquist (in low sample rate patients)
        
        S_norm = nan(size(S)); BL = [];
        % Normalizes each condition by its baseline power:
        switch norm_flag            
            case 0
                disp('*** Normalization: Gain Model ***')
                for k = unique(Rind.blockid)'
                    % Divide by the average of each RUN:
                    ind = contains(Rind.blockid,k);
                    S_norm(:,:,ind) = bsxfun(@rdivide,S(:,:,ind), nangeomean(nangeomean(S(:,:,ind),2),3));
                end
            case 1
                disp('*** Normalization: Addative Model ***')
                % File based:
                for k = unique(Rind.blockid)
                    ind = contains(Rind.blockid,k); % normalized separatly in each RUN:
                    BL = S(:,:,ind);
                    mstd = sqrt(nanmean(nanvar(BL,0,2),3));
                    mbase = nanmean(nanmean(BL,2),3);
                    S_norm(:,:,ind) = bsxfun(@rdivide,bsxfun(@minus,S(:,:,ind),mbase),mstd); % Single trial zscore normalization
                    clear mbase mstd BL
                end
        end
      

        switch norm_flag
            case 0, S_norm=single(10*log10(S_norm));
            case 1, S_norm=single(S_norm); 
        end
        
        % single trial baseline correction:
       % S_norm = bsxfun(@minus,S_norm,nanmean(S_norm(:,baseline_ind,:),2));
     
        epochlim=[T(1) T(end)];        
       
        % High-frequency broadband power (HFB):
        F=[60 160];
        HFB=squeeze(nanmean(S_norm((f_out>F(1) & f_out<F(2)),:,:),1));

        %======================================================================
        %figure position:
        POS = [0 0 200 180];
        
        % ================
        % draw HFB response:
        H = figure('color','w','name',sprintf('HFB clip response electrode %s',channel),'position',POS); hold on;
        shadedErrorBar(t_out_sec,mean(HFB(:,ismember(Rind.blockid,'MT')),2),std(HFB(:,ismember(Rind.blockid,'MT')),[],2)./sqrt(size(HFB,2)),{'color',COLOR.blue,'linewidth',1},0.5)
        shadedErrorBar(t_out_sec,mean(HFB(:,ismember(Rind.blockid,{'M1','M2'})),2),std(HFB(:,ismember(Rind.blockid,{'M1','M2'})),[],2)./sqrt(size(HFB,2)),{'color',COLOR.red,'linewidth',1},0.5)
        shadedErrorBar(t_out_sec,mean(HFB(:,ismember(Rind.blockid,{'R1','R2'})),2),std(HFB(:,ismember(Rind.blockid,{'R1','R2'})),[],2)./sqrt(size(HFB,2)),{'color',COLOR.black,'linewidth',1},0.5)
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
        params_plot.xlim = [-0.5 0.5];
        params_plot.logscaleflag = 0;
        params_plot.cbflag = 1;
        params_plot.verticalLine = 0;
        params_plot.normflag = norm_flag;
        params_plot.position = POSAX;
        params_plot.cb_position = POSCB;
        params_plot.fig_position = POS;
        params_plot.t_out = t_out_sec;
        params_plot.f_out = f_out;
        params_plot.f_out_log = f_out_log;
        params_plot.f_log_ind = f_log_ind;
        params_plot.stimulus_onset = []; % in sec
        params_plot.stimulus_dur = [];   % in sec
        params_plot.ch_label = channel;
        params_plot_grplevel = params_plot; % keep for grp level analysis
       
        %% Figure 1 - All:
        params_plot.title = sprintf('ALL [N=%d]',EEGe.trials);
        params_plot.name = sprintf('ALL ripples spectrogram electrode %s',params_plot.ch_label);    
        [H1,h1,hcb1] = plot_spectrogram(S_norm,params_plot);
        axes(h1); axis square
        % find max freq:
        tmp = squeeze(mean(S_norm,3)); sz = size(tmp);
        tmp = tmp(:); 
        [mx,ind] = max(tmp);
        [fmx, tmx] = ind2sub(sz,ind);
        fprintf('\n max time: %.2f ms, max freq.: %.2f Hz,',t_out(tmx),f_out(fmx))
        hold on; 
        arrow([t_out_sec(tmx)+0.15, f_out(fmx)],[t_out_sec(tmx)+0.05, f_out(fmx)],10,...
                'Width',0.5,'BaseAngle',90,'edgecolor',[0 0 0],'facecolor',[0 0 0]);
        text(t_out_sec(tmx)+0.15, f_out(fmx)+4,sprintf('\n %.1f Hz',f_out(fmx)),...
                    'fontsize',6,'color','k','VerticalAlignment','middle')
        if saveflag, save_current_figure(H1,figdir,1); end
        
        %% Figure 2 - Memory test:
        params_plot.title = sprintf('Memory test [N=%d]',sum(ismember(Rind.blockid,'MT')));
        params_plot.name = sprintf('MT ripples spectrogram electrode %s',params_plot.ch_label);    
        [H2,h2,hcb2] = plot_spectrogram(S_norm(:,:,ismember(Rind.blockid,'MT')),params_plot);
        axes(h2); axis square 
        if saveflag, save_current_figure(H2,figdir,1); end
        % Figure 3 - Movie:
        params_plot.title = sprintf('Movie [N=%d]',sum(ismember(Rind.blockid,{'M1','M2'})));
        params_plot.name = sprintf('Movie ripples spectrogram electrode %s',params_plot.ch_label);    
        [H3,h3,hcb3] = plot_spectrogram(S_norm(:,:,ismember(Rind.blockid,{'M1','M2'})),params_plot);
        axes(h3); axis square 
        if saveflag, save_current_figure(H3,figdir,1); end
        % Figure 4 - Rest:
        params_plot.title = sprintf('Rest [N=%d]',sum(ismember(Rind.blockid,{'R1','R2'})));
        params_plot.name = sprintf('Rest ripples spectrogram electrode %s',params_plot.ch_label);    
        [H4,h4,hcb4] = plot_spectrogram(S_norm(:,:,ismember(Rind.blockid,{'R1','R2'})),params_plot);
        axes(h4); axis square 
        if saveflag, save_current_figure(H4,figdir,1); end
        %%      
        counter = counter+1;
        % DATA structure for group analysis:
        % ERP =============================================================
        DATA.ERP{counter}=squeeze(EEGe.data(cnum,:,:));
        DATA.HFB{counter}=HFB;
        DATA.normspec{counter}=S_norm;
        DATA.Rind{counter}=Rind;
        DATA.T{counter}=EEGe.times;      
        chINFO{counter,1}=channel;
        chINFO{counter,2}=subjid;
        chINFO{counter,3}=cnum;
        
        drawnow;        
    end
    
    if saveflag
        electrodes_table = array2table(chINFO,'VariableNames',{'channel','subjid','cnum'});
        save(fullfile(outdir,'ripple_spectrogram_parameters'),'params','params_plot','f_out','t_out');
        save(fullfile(outdir,[subjid '_ripple_spectrogram_DATA']),'DATA','-v7.3')
        save(fullfile(outdir,[subjid '_electrodes_table']),'electrodes_table')
        fprintf('\n ---> DATA SAVED: %s',subjid);
    end
end






