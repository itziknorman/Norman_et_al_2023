% The following script loads the ripple-triggered wavelets
% spectrograms and draw the panels of Figure 2.
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
addpath(fullfile(path_to_toolboxes,'fieldtrip-20221022/')); ft_defaults;
addpath(genpath(fullfile(path_to_toolboxes,'dmgroppe-Mass_Univariate_ERP_Toolbox-d1e60d4')));
[ALLEEG, EEG, CURRENTSET] = eeglab;

subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};
set(0,'DefaultAxesFontName', 'Arial')
blocks = {'R1','R2','M1','M2','MT'};
% Set Manually:
% 1 = Common Ref;
% 2 = Bipolar montage;
ref_flag = 2;

figdir=fullfile(parentfolder,'results','ripple_triggered_spectrograms','group-level results',['ref_' num2str(ref_flag)]);
if ~exist(figdir,'dir')
    mkdir(figdir);
    disp('Creating Output Directory...')
end

DATA = struct;
% Init the DATA structure:
DATA.ERP=[];
DATA.HFB=[];
DATA.normspec=[];
DATA.Rind=[];
chINFO = {};
for iSub=1:length(subjects)
    
    clearvars -except subjects blocks iSub outdir path_to_toolboxes parentfolder DATAGRP PARAMSGRP chINFO figdir ref_flag
    subjid = subjects{iSub};
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    
    % load datasets:
    maindir=fullfile(parentfolder,subjid);
    datasetID='pink_panther';
    time_locking_event = 'peak';
    
    % Set outdir:
    %indir=fullfile(parentfolder,'results','data','ripple_triggered_spectrograms',['RTS_data_ref_' num2str(ref_flag) '_ripple_' time_locking_event]);
    indir=fullfile(parentfolder,'results','data','ripple_triggered_spectrograms/',['ref_' num2str(ref_flag)]);
    
    
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0);
    load(fullfile(indir,sprintf('%s_electrodes_table.mat',subjid)));
    load(fullfile(indir,sprintf('%s_ripple_spectrogram_DATA.mat',subjid)));
    load(fullfile(indir,'ripple_spectrogram_parameters.mat'));
    
    
    %load(fullfile(indir,sprintf('%s_electrodes_table.mat',subjid)));
    
    if iSub==1, PARAMSGRP = params_plot;
    else, PARAMSGRP = CatStructFields(PARAMSGRP,params_plot,2); end
    
    if iSub==1, DATAGRP = DATA;
    else, DATAGRP = CatStructFields(DATAGRP,DATA,2); end
    
    if iSub==1, chINFO = electrodes_table;
    else,chINFO = cat(1,chINFO,electrodes_table); end
    fprintf('\nloading subject %s ... O.K.',subjid);
    
end


%% RUN ANALYSIS AND PLOT:
set_figure_colors;
close all;
saveflag = 0;
counter = 0;
CM = params_plot.CM;

params_plot.clim = [-8 8];
params_plot.xlim = [-0.3 0.3];
params_plot.ch_label = '';
params_plot.fig_position = [0 0 160 150];
% find peak freq and block-wise ripple rate:
t_out_sec = params_plot.t_out/1000;
[~,peakTimeInd] = min(abs(t_out_sec));
[~,peakFreqInd] = cellfun(@(x)max(x(:,peakTimeInd,:),[],1),DATAGRP.normspec,'UniformOutput',0);
PeakFreq = cellfun(@(x)f_out(squeeze(x)),peakFreqInd,'UniformOutput',0);
PFavg = cellfun(@(x)nanmean(x),PeakFreq,'UniformOutput',0); % peak freq individual elec
PFsem = cellfun(@(x)nanstd(x)./sqrt(length(x)),PeakFreq,'UniformOutput',0); % sem peak freq individual elec

% Figure 1 - Grand Average:
SPavg = cellfun(@(x)nanmean(x,3),DATAGRP.normspec,'UniformOutput',0);
avgspec = nanmean(cat(3,SPavg{:}),3);
params_plot.CM = turbo(64);
params_plot.name = sprintf('Ripple-triggered spectrogram Grand average');
params_plot.title = sprintf('Ripple-triggered spectrogram\n Grand average, n=%d sites',length(SPavg));
params_plot.xlabel = 'Time from ripple peak (s)';

[H1,h1,hcb1] = plot_spectrogram(avgspec,params_plot);
%hcb1.YTickLabel = {'-1','','8'};
axes(h1); axis square
% grand average peak freq:
mean_pf = mean(cell2mat(PFavg));
sd_pf = std(cell2mat(PFavg));
tmpy = get(gca,'ylim'); tmpx = get(gca,'xlim');
text(tmpx(1)+0.5*range(tmpx), tmpy(2)-0.05*range(tmpy),[sprintf('Peak freq.: %.1f ',mean_pf) '(\pm',sprintf('%.1f) Hz',sd_pf)],'FontSize',6,'color','w','HorizontalAlignment','center');
% arrow([t_out_sec(tmx)+0.15, f_out(fmx)],[t_out_sec(tmx)+0.05, f_out(fmx)],10,...
%     'Width',0.5,'BaseAngle',90,'edgecolor',[0 0 0],'facecolor',[0 0 0]);
% text(t_out_sec(tmx)+0.15, f_out(fmx)+4,sprintf('\n %.1f Hz',f_out(fmx)),...
%     'fontsize',6,'color','k','VerticalAlignment','middle')
if saveflag, save_current_figure(H1,figdir,1); end

% Plot ERP:
Tcut = cellfun(@(x)abs(x)<=300,DATAGRP.T,'UniformOutput',0);
ERPavg = cellfun(@(x,y)mean(x(y,:),2),DATAGRP.ERP,Tcut,'UniformOutput',0);
n = cell2mat(cellfun(@(x)size(x,1),DATAGRP.ERP,'UniformOutput',0));
tmp = ERPavg(n==max(n));
erp = mean(cat(2,tmp{:}),2);
t_out_erp = cellfun(@(x,y)x(y),DATAGRP.T,Tcut,'UniformOutput',0);
t_out_erp = t_out_erp{n==max(n)};
Hgranderp = figure('color','w','name',sprintf('ripple erp grand avg'),'position',[0 0 180 180]); hold on;
plot(t_out_erp./1000,erp,'k')
axis square
ylim([-15, 15])
xlabel('Time (s)');
ylabel('Voltage (muV)');
set(gca,'ytick',[-15:5:15])
title({'peri-ripple field potential',''},'fontweight','normal');

if saveflag, save_current_figure(Hgranderp,figdir,1); end



%% Supp Figure - example patients:
subjectToPlot = {'PP04','PP07','PP08'}; 
% load demographic:
demographics = readtable(fullfile(parentfolder,'subj_info_anon.xlsx'));

PrePostExtraTime = 10;
% Conditions duration:
RESTduration = 180.15*2 + 10; % sec
MOVIEduration = 368.5*2 + 10; % sec
MTduration = 324.620 + 10;  % sec
tasks = {'all ripples','quiet rest','watching video','memory test'};
for iSub = 1:length(subjectToPlot)
    subjid = subjectToPlot{iSub};
    Herp(iSub) = figure('color','w','name',sprintf('%s ripple erp',subjid),'position',[0 0 300 160]); hold on;
    [curHipCh,~,~] = define_hippocampal_channels(subjid,ref_flag,0);
    
    for iTask = 1:numel(tasks)
        curTask = tasks{iTask};
        switch curTask
            case 'quiet rest',   blockID = {'R1','R2'}; taskDur = RESTduration; col = COLOR.rest1;
            case 'watching video',  blockID = {'M1','M2'}; taskDur = MOVIEduration; col = COLOR.movie1;
            case 'memory test',  blockID = {'MT'}; taskDur = MTduration; col = COLOR.memorytest;
            case 'all ripples', blockID = {'R1','R2','M1','M2','MT'}; taskDur = RESTduration + MOVIEduration + MTduration; col = COLOR.black;
        end
        col = COLOR.black;
        % Retrieve the data:
        IX = cellfun(@(x)string(x.blockid),DATAGRP.Rind,'UniformOutput',0);
        SPavg = cellfun(@(x,y)nanmean(x(:,:,ismember(y,blockID)),3),DATAGRP.normspec,IX,'UniformOutput',0);
        ERPavg = cellfun(@(x,y)nanmean(x(:,ismember(y,blockID)),2),DATAGRP.ERP,IX,'UniformOutput',0);
        NN = cellfun(@(x,y)size(x(:,:,ismember(y,blockID)),3),DATAGRP.normspec,IX,'UniformOutput',0);
        PFavg = cellfun(@(x,y)nanmean(x(:,ismember(y,blockID))),PeakFreq,IX,'UniformOutput',0); % peak freq individual elec
        PFsem = cellfun(@(x,y)nanstd(x(:,ismember(y,blockID)))./sqrt(length(x(:,ismember(y,blockID)))),PeakFreq,IX,'UniformOutput',0); % sem peak freq individual elec
        channel_ix = strcmpi(chINFO.subjid,subjid) & strcmpi(chINFO.channel,curHipCh);
        
        % plot spectrogram:
        avgspec = SPavg{channel_ix};
        nn = NN{channel_ix};
        ix = IX{channel_ix}; % trial indices
        params_plot.title = sprintf('subject #%d (%d y/o), %s\nn=%d ripples',find(strcmpi(demographics.SubjectID,subjid)),floor(demographics.Age(strcmpi(demographics.SubjectID,subjid))), curTask,nn);
        params_plot.name = sprintf('%s %s %s ripple spectrogram',subjid,curHipCh,curTask);
        [H1,h1,hcb1] = plot_spectrogram(avgspec,params_plot);
        %hcb1.YTickLabel = {'-1','','8'};
        axes(h1); axis square
        
        % find max freq:
        t_out_sec = params_plot.t_out/1000;
        
        mean_pf = PFavg{channel_ix};
        sd_pf = PFsem{channel_ix};
        tmpy = get(gca,'ylim'); tmpx = get(gca,'xlim');
        text(tmpx(1)+0.5*range(tmpx), tmpy(2)-0.05*range(tmpy),sprintf('Peak freq.: %.1f (+%.1f) Hz',mean_pf,sd_pf),'FontSize',6,'color','w','HorizontalAlignment','center');
          
        fprintf('\n peak freq.: %.1f (%.1f) Hz,',mean_pf,sd_pf)

        if saveflag, save_current_figure(H1,figdir,1); end
        
        % plot ERP:       
       
        figure(Herp(iSub)); 
        erp = ERPavg{channel_ix};
        yscale = max([100 max(abs(erp))]);
        erp = (erp-mean(erp)) + iTask*yscale/2;
        
        t_out_erp = DATAGRP.T{channel_ix}./1000;   
        cutind = abs(t_out_erp)<0.5;
        
        hold on; 
        plot(t_out_erp(cutind),erp(cutind),'color',col)
        axis tight square off
        xlim([-0.5 0.5]);    
        % scale bar:
        
        ypos = 0;
        xpos = 0.6;
        Yscalebar = 20;
        Xscalebar = 0.1;
        h2a = plot([xpos xpos],[ypos+Yscalebar ypos],'k','Linesmoothing','on');
        h2b = plot([xpos-Xscalebar xpos],[ypos ypos],'k','Linesmoothing','on');
        if iTask == 1
        text(mean([xpos-Xscalebar xpos]),ypos,sprintf('%d ms',Xscalebar*1000),...
            'HorizontalAlignment','center','VerticalAlignment','top','fontsize',5);
        text(xpos,mean([ypos+Yscalebar ypos]),sprintf(' %d muV',Yscalebar),...
            'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',5);        
        else, set([h2a,h2b],'visible','off')            
        end
        axis tight         
        tmpy = get(gca,'ylim'); tmpx = get(gca,'xlim');
        text(tmpx(1)-0.1*range(tmpx), mean(erp), tasks{iTask},'fontsize',6,'HorizontalAlignment','right')
        title({'peri-ripple field potential';''},'fontweight','normal');

    end
    if saveflag, save_current_figure(Herp(iSub),figdir,1); end
end

