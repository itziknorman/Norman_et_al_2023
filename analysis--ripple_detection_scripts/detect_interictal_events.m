% The script run IED detection in hippocampal sites for later exclusion of candidate hippocampal ripple events. 
% calling the function LLspikedetector.m (By Kleen-Lab)
%
% Author: Itzik Norman, Chang lab 2022

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

% Load data:
subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};
blocks = {'R1','M1','R2','M2','MT'};

for iSub=1 %:numel(subjects)
    
    
    clearvars -except subjects blocks iSub ALLEEG EEG parentfolder path_to_toolboxes outdir
    subjid = subjects{iSub};
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    PEepoches = {}; PEchlist = {};
    RIPPLES = [];
        
    warning('off','all')    
    
   
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
          
    
    
    %% Defining Hippocampus:
    set_figure_colors;
   
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0);
    saveflag = 1;
    for current_hippocampus = hippocampus_all_channels
                
        cnum1=find(strcmpi({EEG.chanlocs.labels},current_hippocampus)); % CA1 channel, identified anatomically
        channel=EEG.chanlocs(cnum1).labels;
               
        % we use the averaged LFP across all good
        % channels (common ref) to perform control detection of ripples,
        % which allows to discard muscular/electrical artifacts.
        
        if isempty(current_hippocampus)
            error('Missing hippocampal channel!')
        end
        
        for iBlock = 1:length(blocks)
            blockID = blocks{iBlock};
            pathologicalEvents = [];
            EEG = ALLEEG(iBlock);
            outdir = fullfile(parentfolder,'results','data','IEDs',subjid,sprintf('ref%d',ref_flag)); 
            if ~exist(outdir,'dir'), mkdir(outdir); disp('folder created...'); end
            outFileName = fullfile(outdir,sprintf('%s %s IEDs %s.mat',subjid,blockID,channel));
            llw = 0.04; % in Sec
            prc = 99.9; 
            [ets,ech]=LLspikedetector(EEG.data(cnum1,:),EEG.srate,llw,prc,[]);
            onsets = ets(:,1)./EEG.srate; % convert to sec
            offsets = ets(:,2)./EEG.srate; % convert to sec
            
            % save results:
            save(outFileName,'onsets','offsets');
            if isempty(onsets), continue; end
            %% ============================================================
            event_list = {}; e = 1;
            timestamps = onsets;            
            % add events to the dataset:
            for index = 1 : length(timestamps)
                current_event_timing = timestamps(index);                
                event_list{e,1} = sprintf('E%d-%s',index,channel);
                event_list{e,2} = current_event_timing;
                event_list{e,3} = channel;
                e = e+1;
            end
            EEG = pop_importevent( EEG, 'event', event_list, 'fields',{'type', 'latency','channel'}, 'timeunit',1, 'append','yes');
            EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
                       
            epochlim = [-0.5 0.5]; % in sec
            EEGepoched = pop_epoch(EEG,{},epochlim,'eventindices',find(endsWith({EEG.event.type},channel)));
            T = EEGepoched.times;
               
            %% Spectograms parameters:
            params = struct;
            params.freqs=[4 floor(EEG.srate./2)];
            params.nfreqs=range(params.freqs);
            params.winsize=128;
            params.cycles= [1 20];
            params.units='abs';
            params.freqscale='linear';
            params.normflag = 0;
            
            
            epochlim = [T(1) T(end)];
            
            params.timesout=150;
            [~,~,~,t_out,f_out,~,~,S] = newtimef(EEGepoched.data(cnum1,:,:),numel(T),epochlim, EEGepoched.srate, 'cycles', params.cycles, 'freqs', params.freqs, 'winsize', params.winsize,...
                'scale',params.units,'plotitc','off','plotersp','off','timesout',params.timesout,'nfreqs',params.nfreqs,'baseline',nan,'freqscale',params.freqscale);
            
            S = S.*conj(S); % power
            S_norm = nan(size(S)); BL = [];
            % Normalizes each condition by its baseline power:
            switch params.normflag
                case 0
                    disp('*** Normalization: Gain Model ***')                    
                    S_norm = bsxfun(@rdivide,S, nangeomean(S,2));
                case 1
                    disp('*** Normalization: Addative Model ***')
                    mstd = sqrt(nanmean(nanvar(S,0,3),2));
                    mbase = nanmean(nanmean(S,3),2);
                    S_norm = bsxfun(@rdivide,bsxfun(@minus,S,mbase),mstd); % Single trial zscore normalization
                    clear mbase mstd
            end
            % transform into dB:
            switch  params.normflag
                case 0, S_norm=single(10*log10(S_norm));
                case 1, S_norm=single(S_norm);
            end
            epochlim=[T(1) T(end)];
            
            %% Plot figure:
            CM=flipud(cbrewer('div','RdBu',32));
            
            figure('color','w','position',[0 0 800 400],'name',sprintf('%s %s - %s IEDs summary',subjid,blockID,channel));
            H = subplot(2,2,[1,3]); hold on;
            imagesc(T, 1:EEGepoched.trials, squeeze(EEGepoched.data(cnum1,:,:))'); caxis([-500 500])
            axis tight;
            colormap(CM);            
            xlabel('Time (ms)'); ylabel('Events');
            cb = colorbar; 
            cb.Label.String = 'Voltage';
            freezeColors;
            

            subplot(2,2,2); hold on;
            plot(EEGepoched.times,nanmean(EEGepoched.data(cnum1,:,:),3),'k','linesmoothing','on');            
             xlabel('Time (s)'); ylabel('Voltage \muV');
            H = subplot(2,2,4); hold on;
            
            % Plot the AVERAGE SPECTROGRAM:            
            avgspec = nanmean(S_norm,3);
            tlim = [ceil(min(t_out)*10)./10, floor(max(t_out)*10)./10];
            imagesc (t_out, f_out, avgspec); 
            axis tight square
            set(gca,'xtick',[-1000:250:1000],'xlim',tlim,'TickDir','out');
            plot([0 0], get(gca,'Ylim'), ':k','LineWidth', 1);
            tmpy = get(gca, 'ylim');
            tmpx = get(gca, 'xlim');
            tmpc = round(get(gca, 'clim'));
            text(tmpx(1)+0.05*range(tmpx), tmpy(2)-0.075*range(tmpy),sprintf('n=%d',size(S_norm,3)),'FontSize',6,'HorizontalAlignment','left');
           
            [~,globmax] = max(avgspec(:)); 
            [maxFind,maxTind] = ind2sub(size(avgspec),globmax);  
            [~,ii] = max(S_norm(:,maxTind,:),[],1);            
            ii = squeeze(ii);
            mean_pf = mean(f_out(ii));
            sem_pf = std(f_out(ii))./sqrt(length(ii));
            text(tmpx(1)+0.5*range(tmpx), tmpy(2)+0.1*range(tmpy),sprintf('Peak freq.: %.2f Hz',f_out(maxFind)),'FontSize',6,'HorizontalAlignment','center');
            ylabel ('Frequency (Hz)'); xlabel('Time (ms)')
            set(gca,'Xcolor','k','ycolor','k')

            % add color bar;
            clim = [0 10];
            pos = get(H,'position');
            cbpos = [pos(1)+pos(3)*0.8,pos(2)+pos(4)*0.3,pos(3)*0.025,pos(4)*0.4];
            h = axes('position',cbpos);
            add_colorbar_to_spectrogram(h,clim,CM,'dB');
            %suptitle(sprintf('%s %s %s',subjid,blockID,channel)
            
            if saveflag
                set_font_size_and_type;
                save_current_figure(gcf,outdir,0);
            end
    end
    end
end