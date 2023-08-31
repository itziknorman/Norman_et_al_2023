% Hippocampal ripple detection main script.
% detects ripples after exclusion of electrical & muscular artifacts, IED and pathological HFOs (optional).
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
rmpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12/')));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));
[ALLEEG, EEG, CURRENTSET] = eeglab;
% =================================
set(0,'DefaultFigureVisible','off')
set(0,'DefaultAxesFontName', 'Arial')
% =================================
% Load data:
subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};
blocks = {'R1','M1','R2','M2','MT'};

for iSub=1 %:numel(subjects)
    
    
    clearvars -except subjects blocks iSub ALLEEG EEG parentfolder path_to_toolboxes
    subjid = subjects{iSub};
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    RIPPLES = [];
        
    warning('off','all')    
    
   
    % load datasets:
   
    maindir=fullfile(parentfolder,subjid);
    datasetID='pink_panther';
    
    % Set Manually:
    % 1 = Common Ref;
    % 2 = Bipolar montage;
    
    ref_flag = 2;
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
    set_figure_colors_PP;
 
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0);
    saveflag = 0;
    ripples = [];
    for current_hippocampus = hippocampus_all_channels
        
        cnum1=find(strcmpi({EEG.chanlocs.labels},current_hippocampus)); % CA1 channel, identified anatomically
        cnum2=find(strcmpi({EEG.chanlocs.labels},'CREF')); % EMG or common ref. channel
        channel=EEG.chanlocs(cnum1).labels;
        RIPPLES(end+1).channel = channel;
        
        % we use the averaged LFP across all good
        % channels (common ref) to perform control detection of ripples,
        % which allows to discard muscular/electrical artifacts.
        
        if isempty(current_hippocampus)
            error('Missing hippocampal channel!')
        end
        
        %% Data normalization and epoching:
        
        %%%%% Ripple detection PARAMETERS %%%%%
        minDistance=0.030; % in sec
        minRippleDuration=0.020; % in sec [0.015]
        maxRippleDuration=0.200; % in sec
        % =================================================================
        th=[2 4]; % ripple detection thresholds (in std from the mean)
        Fs=EEG.srate; % in Hz
         if Fs >= 500
             ripple_band=[70 180]; % in Hz
             forder_BP = 170; % 10 Hz transition
             ftype = 'bandpass';             
        else
             ripple_band=[70]; % in Hz
             forder_BP = 86; % 10 Hz transition
             ftype = 'highpass';
        end
                  
        % note: the freqency range is a bit wider to account for a transition bandwith of ~10Hz
        % (the band of interest is 80-170 Hz)
        
        IED_band = [20 40]; % see Smith et al (2022) for similar method
        IEDthr = 8;       
        LPcutoff_ripple_band = 40; % in Hz (fixed size LP filter)  
        LPcutoff_IED_band = round(mean(IED_band)/pi); % in Hz 
        wtype = 'hamming';
        warg = 5.653; % beta value, in case you choose Kaiser window
        % =================================================================
        
        % set output folder:
        outdir=fullfile(parentfolder,'results','data',...
            sprintf('Ripple_times_%s_%d-%dstd_adjusted_band_%dms_%dms',...
            wtype,th(1),th(2),minRippleDuration*1000,minDistance*1000),subjid); % ADJUST OUTDIR
        
        figdir=fullfile(outdir,'Detection_Figures');
        if ~exist(outdir,'dir'), mkdir(outdir);  disp('Creating Output Directory...'); end
        if ~exist(figdir,'dir'), mkdir(figdir);  disp('Creating Figures Directory...'); end

     

        for iBlock = 1:numel(ALLEEG)
            blockID = blocks{iBlock};
            clear ripples avgA stdevA avgB stdevB
            close all
            % Compute mean and std of ripple band amplitude across the entire block: 
                 
            % ripple band:
            EEGtmp = ALLEEG(iBlock);    
            EEGtmp = pop_select(EEGtmp,'nochannel',setdiff(1:EEGtmp.nbchan,cnum1));            
            EEGbp = pop_firws(EEGtmp, 'fcutoff', ripple_band, 'ftype', ftype, 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
            rawSignal = double(EEGtmp.data);
            bandpassSignalA = double(EEGbp.data);
                        
            % control detection channel (common ref):
            EEGtmp = ALLEEG(iBlock);    
            EEGtmp = pop_select(EEGtmp,'nochannel',setdiff(1:EEGtmp.nbchan,cnum2));            
            EEGbp = pop_firws(EEGtmp, 'fcutoff', ripple_band, 'ftype', ftype, 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
            bandpassSignalB = double(EEGbp.data);
                   
            % IED band (METHOD II):
            EEGtmp = ALLEEG(iBlock);    
            EEGtmp = pop_select(EEGtmp,'nochannel',setdiff(1:EEGtmp.nbchan,cnum1));            
            EEGbp = pop_firws(EEGtmp, 'fcutoff', IED_band, 'ftype', 'bandpass', 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
            bandpassSignalC = double(EEGbp.data);
            
            absSignalA = bandpassSignalA;
            absSignalA(~isnan(absSignalA)) = abs(hilbert(bandpassSignalA(~isnan(bandpassSignalA)))); % hilbert envelope
            absSignalB = bandpassSignalB;
            absSignalB(~isnan(absSignalB)) = abs(hilbert(bandpassSignalB(~isnan(bandpassSignalB)))); % hilbert envelope
            absSignalC = bandpassSignalC;
            absSignalC(~isnan(absSignalC)) = abs(hilbert(bandpassSignalC(~isnan(bandpassSignalC)))); % hilbert envelope
            
                        
            % Defining thr for clipping (hippocmpus):
            [robustAvgA,robustStdevA] = robustMean(absSignalA,2,th(2));     % robust mean and std used for clipping
            topLimA = robustAvgA+th(2)*robustStdevA;           
            absSignalA(absSignalA>topLimA) = topLimA; % clipping the signal
            squaredSignalA = absSignalA.^2;           % squared amplitude
            
            % Defining thr for clipping (control channel):
            [robustAvgB,robustStdevB] = robustMean(absSignalB,2,th(2));     % robust mean and std used for clipping
            topLimB = robustAvgB+th(2)*robustStdevB;           
            absSignalB(absSignalB>topLimB) = topLimB; % clipping the signal
            squaredSignalB = absSignalB.^2;           % squared amplitude
            
             % Defining thr for clipping (hippocampal IED band):
            [robustAvgC,robustStdevC] = robustMean(absSignalC,2,IEDthr);     % robust mean and std used for clipping
            topLimC = robustAvgC + IEDthr*robustStdevC;           
            absSignalC(absSignalC>topLimC) = topLimC; % clipping the signal
            squaredSignalC = absSignalC.^2;           % squared amplitude
            
            % low pass filter the hippocampal signal:            
            EEGtmp = ALLEEG(iBlock);
            EEGtmp.data(1,:) = squaredSignalA;
            EEGtmp = pop_select(EEGtmp,'nochannel',2:EEGtmp.nbchan);
            assert(size(EEGtmp.data,1)==1)                  
            [EEGtmp] = pop_firws(EEGtmp, 'fcutoff', LPcutoff_ripple_band, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', forder_BP);
            squaredSignalA = double(EEGtmp.data);
            
            % low pass filter the control signal:            
            EEGtmp = ALLEEG(iBlock);
            EEGtmp.data(1,:) = squaredSignalB;
            EEGtmp = pop_select(EEGtmp,'nochannel',2:EEGtmp.nbchan);
            assert(size(EEGtmp.data,1)==1)                  
            [EEGtmp] = pop_firws(EEGtmp, 'fcutoff', LPcutoff_ripple_band, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', forder_BP);
            squaredSignalB = double(EEGtmp.data);            
            
            % low pass filter the IED band signal:            
            EEGtmp = ALLEEG(iBlock);
            EEGtmp.data(1,:) = squaredSignalC;
            EEGtmp = pop_select(EEGtmp,'nochannel',2:EEGtmp.nbchan);
            assert(size(EEGtmp.data,1)==1)                  
            [EEGtmp] = pop_firws(EEGtmp, 'fcutoff', LPcutoff_IED_band, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', forder_BP);
            squaredSignalC = double(EEGtmp.data);     
            
            % Compute mean and std of the clipped & squared time series:
            avgA = nanmean(squaredSignalA);
            stdevA = nanstd(squaredSignalA,[],2);
            avgB = nanmean(squaredSignalB);
            stdevB = nanstd(squaredSignalB,[],2);
            avgC = nanmean(squaredSignalC);
            stdevC = nanstd(squaredSignalC,[],2);
            
            clear bandpassSignalA bandpassSignalB bandpassSignalC squaredSignalA squaredSignalB squaredSignalC
            close all
                     
            
            %% now reprocess the unclipped signal:            
            % ripple band:
            EEGtmp = ALLEEG(iBlock);    
            EEGtmp = pop_select(EEGtmp,'nochannel',setdiff(1:EEGtmp.nbchan,cnum1));            
            EEGbp  = pop_firws(EEGtmp, 'fcutoff', ripple_band, 'ftype', ftype, 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
            rawSignal = double(EEGtmp.data);
            bandpassSignalA = double(EEGbp.data);
            
            % control detection channel (common ref):
            EEGtmp = ALLEEG(iBlock);    
            EEGtmp = pop_select(EEGtmp,'nochannel',setdiff(1:EEGtmp.nbchan,cnum2));            
            EEGbp  = pop_firws(EEGtmp, 'fcutoff', ripple_band, 'ftype', ftype, 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
            bandpassSignalB = double(EEGbp.data);
            
            % IED band:
            EEGtmp = ALLEEG(iBlock);    
            EEGtmp = pop_select(EEGtmp,'nochannel',setdiff(1:EEGtmp.nbchan,cnum2));            
            EEGbp  = pop_firws(EEGtmp, 'fcutoff', IED_band, 'ftype', 'bandpass' , 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
            bandpassSignalC = double(EEGbp.data);
                        
                                    
            % time vector:
            strind = round(ALLEEG(iBlock).event(1).latency); % finds the start point
            disp('adjusting extra time points...'); 
            tail_correction = ALLEEG(iBlock).times(strind)./1000; % convert to sec
            fprintf('\n Extra "margins" of %.2f s', tail_correction);
            T = double(ALLEEG(iBlock).times./1000) - tail_correction; % in Sec
            
            % Hilbert envelope (rectification):
            absSignalA = bandpassSignalA;
            absSignalA(~isnan(absSignalA)) = abs(hilbert(bandpassSignalA(~isnan(bandpassSignalA)))); % hilbert envelope
            absSignalB = bandpassSignalB;
            absSignalB(~isnan(absSignalB)) = abs(hilbert(bandpassSignalB(~isnan(bandpassSignalB)))); % hilbert envelope
            absSignalC = bandpassSignalC;
            absSignalC(~isnan(absSignalC)) = abs(hilbert(bandpassSignalC(~isnan(bandpassSignalC)))); % hilbert envelope
         
            % Squaring the signal:
            squaredSignalA = absSignalA.^2;
            squaredSignalB = absSignalB.^2;
            squaredSignalC = absSignalC.^2;
            
            % FIR filter (after squaring):
            % low pass filter the hippocampal signal:            
            EEGtmp = ALLEEG(iBlock);
            EEGtmp.data(1,:) = squaredSignalA;
            EEGtmp = pop_select(EEGtmp,'nochannel',2:EEGtmp.nbchan);
            assert(size(EEGtmp.data,1)==1)                  
            [EEGtmp] = pop_firws(EEGtmp, 'fcutoff', LPcutoff_ripple_band, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', forder_BP);
            squaredSignalA = double(EEGtmp.data);
            
            % low pass filter the control signal:            
            EEGtmp = ALLEEG(iBlock);
            EEGtmp.data(1,:) = squaredSignalB;
            EEGtmp = pop_select(EEGtmp,'nochannel',2:EEGtmp.nbchan);
            assert(size(EEGtmp.data,1)==1)                  
            [EEGtmp] = pop_firws(EEGtmp, 'fcutoff', LPcutoff_ripple_band, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', forder_BP);
            squaredSignalB = double(EEGtmp.data);                                       
           
            % low pass filter the control signal:            
            EEGtmp = ALLEEG(iBlock);
            EEGtmp.data(1,:) = squaredSignalC;
            EEGtmp = pop_select(EEGtmp,'nochannel',2:EEGtmp.nbchan);
            assert(size(EEGtmp.data,1)==1)                  
            [EEGtmp] = pop_firws(EEGtmp, 'fcutoff', LPcutoff_IED_band, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', forder_BP);
            squaredSignalC = double(EEGtmp.data);   
            
            % ZSCORE:
            squaredSignalNormA = (squaredSignalA-avgA)/stdevA; % z-score hippocampal ripple band amplitude
            squaredSignalNormB = (squaredSignalB-avgB)/stdevB; % z-score the control signal (extracted from common average)
            squaredSignalNormC = (squaredSignalC-avgC)/stdevC; % z-score hippocampal IED band 
            %squaredSignalNormC = (squaredSignalC-nanmean(squaredSignalC))/std(squaredSignalC); % z-score hippocampal IED band 
        
                        
            % load IED onsets (METHOD I): 
            % detected using the line-length algorithm (e.g. Baud et al. 2001), implemented in "detect_interictal_events.m"
            IED_onsets = zeros(size(rawSignal));  % prepare timeseries for IEDs                
            % Try to load external IEDs onsets:
            IEDsdir = fullfile(parentfolder,'results','data','IEDs',subjid,sprintf('ref%d',ref_flag));
            IEDFile = fullfile(IEDsdir,sprintf('%s %s IEDs %s.mat',subjid,blockID,cell2mat(current_hippocampus)));            
            clear IEDs
            loadExternalEventsFlag = 1;  % activate this flag if you want to exclude SWR that coincide with pathological IEDs
            if exist(IEDFile,'file') && loadExternalEventsFlag, IEDs = load(IEDFile); end
            if ~isempty(IEDs.onsets)
                for ii = 1:length(IEDs.onsets)
                    [~,mn] = min(abs((ALLEEG(iBlock).times/ALLEEG(iBlock).srate)-IEDs.onsets(ii)));
                    IED_onsets(mn)=1;
                end
            end
            clear IEDs
                         
            % Now add IED onsets detected via method II: Smith et al 2022 DOI: 10.7554/eLife.73541            
            ind = find(squaredSignalNormC > IEDthr);
            IED_onsets(ind) = 1;
            
            % Exclude sharp transient in the raw EEG:
            dv = [0 diff(rawSignal)];
            [avgD,stdD] = robustMean(dv,[],IEDthr);
            dvnorm = (dv-avgD)./stdD;
            %sharpTransientsInd = dv > 50 * (1000/EEG.srate);  % excluding sharp transients (>50muV/ms)
            sharpTransientsInd = dvnorm > IEDthr;  % excluding sharp transients
            IED_onsets(sharpTransientsInd) = 1;          
            
            % find boundary events (if any):
            boundary_onsets = [];
            boundaryEventInd = find(strcmpi({ALLEEG(iBlock).event.type},'boundary'));
            if ~isempty(boundaryEventInd)
                boundary_latencies = [ALLEEG(iBlock).event(boundaryEventInd).latency]; % in samples
                boundary_latencies = round(boundary_latencies);% make this integers
                boundary_onsets = zeros(size(rawSignal));
                boundary_onsets(boundary_latencies) = 1;
            end
                                               
            % Find ripples:
            [ripples,ripples_stat,rnoise] = ripples_detection_algorithm(squaredSignalNormA,bandpassSignalA,T,Fs,th,minDistance,...
                                            minRippleDuration,maxRippleDuration,squaredSignalNormB,[],IED_onsets,boundary_onsets);
            outfilename = sprintf('%s %s ripples %s',subjid,blocks{iBlock},cell2mat(current_hippocampus));
            eval([sprintf('t_%s',blocks{iBlock}) '=T;'])
            save(fullfile(outdir,outfilename),'ripples','ripples_stat',sprintf('t_%s',blocks{iBlock}))
            
            RIPPLES(end).swramp = nanmean(ripples.amplitude);
            RIPPLES(end).blavg = avgA;
            RIPPLES(end).blstdev = stdevA;
            RIPPLES(end).count = size(ripples,1);
            RIPPLES(end).rnoise = rnoise;
            
            %% Plot
            close all
            dataset_label=sprintf('%s block %d',datasetID,iBlock); dataset_label(dataset_label=='_')=' ';
            tlim = [0 60];
            H1=figure('Name',[subjid ' ' dataset_label ' ripples detection Ch ' channel],'position',[0 100 1000 300],'Color','w','visible','on');
            
            subplot(2,1,2); hold on;
            h2 = plot(T,squaredSignalNormA,'color',COLOR.red,'Linesmoothing','on');
            h1 = plot(T,ones(size(T))*th(2),'k','Linewidth',0.5,'visible','off');
            h0 = plot(T,ones(size(T))*th(2),'b--','Linewidth',0.5);
            if ~isempty(ripples)
                h3=scatter(ripples.peak,ones(size(ripples,1),1)*-1.8,20,'yo','fill','markeredgecolor','k'); hold on
            end
            if length(ripple_band)==2
                ylabel(sprintf('%d-%dHz Amplitude^2 (Z-score)',ripple_band(1),ripple_band(2)))
            else
                ylabel(sprintf('>%d Hz Amplitude^2 (Z-score)',ripple_band(1)))
            end
            xlim([T(1),T(end)])
            ylim([-2,10])
            set(gca,'ytick',[0:2:10])
            xlim(tlim)
         
            if ~isempty(ripples)
                legend([h1 h2 h3],{'Ripple band envelope^2 (zscore)','Ripple band (A.U)','Ripple event'},'location','northeast'); legend boxoff;
            else
                legend([h1 h2],{'Ripple band envelope^2 (zscore)','Ripple band (A.U)'},'location','northeast'); legend boxoff;
            end
    
                
            subplot(2,1,1); hold on;
            title(sprintf('%s Ripples Detection - %s (ripple band envelope)',subjid,dataset_label),'fontweight','normal');
            plot(T,nanzscore(bandpassSignalA)+18,'linewidth',0.25,'color',COLOR.gray); axis tight;
            set(gca,'xtick',[],'ytick',[],'Xcolor','none','Ycolor','none');
            xlim(tlim)
            
            % save figure;
            if saveflag                
                set(0,'DefaultAxesFontName', 'Arial')
                set(gcf,'renderer','painters')
                save_current_figure(gcf,figdir,0)                 
            end
            
            %% plot detection exmple for Fig. 1:
            plotExampleFlag = 1; % change to true if you want to plot a representative ripple detection example
            if ismember(subjid,{'PP08','PP07'}) && strcmpi(current_hippocampus,hippocampus) && iBlock==1 && plotExampleFlag
                figure('name',sprintf('Detection example %s (%s)',subjid,cell2mat(current_hippocampus)),'color','w','position',[0 0 350 250]);
                switch subjid
                    case 'PP08',t1 = 102; t2 = 112; ytop = 10;
                    case 'PP07',t1 = 15; t2 = 21; ytop = 10;
                end
                %
                %
                subplot(2,1,1); hold on;
                h = plot(T,bandpassSignalA,'linewidth',0.5,'color','k'); hold on;
                ylabel(sprintf('%d-%dHz\n(muV)',ripple_band(1),ripple_band(2)))
                xlim([t1,t2]); ylim([-20 20])
                set(gca,'xtick',[],'xcolor','none')
              
                
                subplot(2,1,2); hold on;
                h1 = plot(T,squaredSignalNormA,'k','Linesmoothing','on');
                h0 = plot(T,ones(size(T))*th(2),'--','color',COLOR.red,'Linewidth',0.5);
                if ~isempty(ripples)
                    h2=scatter(ripples.peak,ones(size(ripples,1),1)*ytop,25,'v','fill','markerfacecolor',COLOR.orange,'markeredgecolor',COLOR.black); hold on
                end
                ylabel(sprintf('Power\n(z-score)'))
                xtck = [t1:2:t2]; ytck = [0 4 ytop];
                set(gca,'xlim',[t1,t2],'xtick',xtck,'XtickLabels',{xtck-t1},'ylim',[-2.5,ytop],'ytick',ytck)
                xlabel('Time (s)')
                [h_leg,h_leg_icon] = legend([h1 h0 h2],{'ripple-band power','detection thr.','ripple event'},'location','NorthWest');
             
                legpos = h_leg.Position;
                h_leg.Position = [legpos(1)*0.9 legpos(2)*1.02 legpos(3:4)];
                h = findobj(h_leg_icon,'type','line');
                for ii = 1:length(h)
                    if length(h(ii).XData)==2
                        h(ii).XData(1) = h(ii).XData(1) + 0.25*range(h(ii).XData);
                        h(ii).LineWidth = 1;
                    end
                end
                legend boxoff
                figdir=fullfile(outdir,'Detection_Figures','ripple_detection_example');
                if ~exist(figdir,'dir')
                    mkdir(figdir);
                    disp('Creating Output Directory...')
                end
                set(0,'DefaultAxesFontName', 'Arial')
                set_font_size_and_type;
                if saveflag, save_current_figure(gcf,figdir,0); end            
            end
            
            %% LFP trace of a single event: 
            if ismember(subjid,{'PP07'}) && strcmpi(current_hippocampus,hippocampus) && iBlock==1 && plotExampleFlag
               %%
                for n = 17; %42 %[17,42] %60 %
                     figure('name',sprintf('Detection example %s (%s) ripple # %d',subjid,cell2mat(current_hippocampus),n),'color','w','position',[0 0 150 250]);
                    
                    t1 = ripples.peak(n)-0.08; t2 = ripples.peak(n)+0.2;
                    %t1 = ripples.peak(n)-0.1; t2 = ripples.peak(n)+0.15; 
                    fprintf('\n ripple #%d \n',n)
               
                    subplot(2,1,1); hold on;
                    h = plot(T,bandpassSignalA,'linewidth',0.5,'color','k'); hold on;
                    ytop = 40;%20;
                    Yscalebar = 20;%10;
                    Xscalebar = 0.05;
                    h2a = plot([t2 t2],[ytop-Yscalebar ytop],'k','Linesmoothing','on');                   
                    h2b = plot([t2-Xscalebar t2],[ytop ytop],'k','Linesmoothing','on');
                    text(mean([t2-Xscalebar t2]),ytop,sprintf('%d ms',Xscalebar*1000),...
                        'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',5);
                    text(t2,mean([ytop-Yscalebar ytop]),sprintf(' %d muV',Yscalebar),...
                        'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',5);
                    title(sprintf('%d-%dHz',ripple_band(1),ripple_band(2)),'fontweight','normal')
                    xlim([t1,t2]); ylim([-ytop ytop])
                    set(gca,'xtick',[],'xcolor','none')
                    
                    axis off
                    
                    subplot(2,1,2); hold on;
                    h1 = plot(T,rawSignal-mean(rawSignal(T>t1&T<t2)),'k','Linesmoothing','on');
                    ytop = 250; %150;
                    Yscalebar = 100;
                    Xscalebar = 0.05;
                    h2a = plot([t2 t2],[ytop-Yscalebar ytop],'k','Linesmoothing','on');                   
                    h2b = plot([t2-Xscalebar t2],[ytop ytop],'k','Linesmoothing','on');
                    text(mean([t2-Xscalebar t2]),ytop,sprintf('%d ms',Xscalebar*1000),...
                        'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',5);
                    text(t2,mean([ytop-Yscalebar ytop]),sprintf(' %d muV',Yscalebar),...
                        'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',5);
                    if ~isempty(ripples)
                        h2=scatter(ripples.peak,ones(size(ripples,1),1)*ytop,25,'v','fill','markerfacecolor',COLOR.orange,'markeredgecolor',COLOR.black); hold on
                    end
                    title(sprintf('Raw LFP'),'fontweight','normal')
                    xtck = [t1:2:t2]; ytck = [-200 0 200 400];
                    set(gca,'xlim',[t1,t2],'xtick',xtck,'XtickLabels',{xtck-t1},'ylim',[-ytop,ytop],'ytick',ytck)
                    xlabel('Time (s)')
                    axis off
                    
                    figdir=fullfile(outdir,'Detection_Figures','ripple_detection_example');
                    if ~exist(figdir,'dir')
                        mkdir(figdir);
                        disp('Creating Output Directory...')
                    end
                  
                    
                    if saveflag, save_current_figure(gcf,figdir,0); end
                    drawnow; pause(0.5);
                    
                end
            end
            
        end
    end
    
    %% Find the channel with the strongest SNR: (this part is still in development...)
    % SNR = ripple power divided by the background ripple-band activity
    channels= {RIPPLES.channel};
    swramp  = [RIPPLES.swramp];
    blavg   = [RIPPLES.blavg];
    blstdev   = [RIPPLES.blstdev];
    swrsnr = 10*log10(swramp./blavg); % SNR calculation
    %swrsnr = (swramp-blavg)./blstdev; % Cohen's d alternative
    %swrsnr = swramp./blstdev; %  alternative SNR calculation
    figure('color','w','name',['compare ripple SNR between electrodes ref ' num2str(ref_flag)],'position',[0 0 200 200]);
    hold on; superbar(1:numel(swrsnr),swrsnr);
    set(gca,'xtick',[1:numel(swrsnr)],'xticklabel',{RIPPLES.channel},'fontsize',5,'Xticklabelrotation',45)
    ylabel('Ripple SNR (dB)')
    rr = [RIPPLES.count]./EEG.xmax;
    enlarge_figure_and_move_axes_to_center(gcf,gca,1.2);
    
    [~, indmax] = max(swrsnr);
    fprintf(' \n BEST SNR IN CHANNEL: %s \n ', channels{indmax});
    title(sprintf('BEST SNR IN CHANNEL: %s', channels{indmax}));   

    set_font_size_and_type;
    save_current_figure(gcf,figdir,0);
    
    figure('color','w','name',['compare ripple rate between electrodes ref ' num2str(ref_flag)],'position',[0 0 200 200]);
    hold on; superbar(1:numel(rr),rr);
    set(gca,'xtick',[1:numel(rr)],'xticklabel',{RIPPLES.channel},'fontsize',5,'Xticklabelrotation',45)    
    ylabel('Ripple rate (Hz)')    
    enlarge_figure_and_move_axes_to_center(gcf,gca,1.2);
    [~, indmax] = max(rr);
    fprintf(' \n BEST RATE IN CHANNEL: %s \n ', channels{indmax});
    title(sprintf('BEST RATE IN CHANNEL: %s', channels{indmax}));  
    set_font_size_and_type;
    save_current_figure(gcf,figdir,0);

end
