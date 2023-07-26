% The following script loads the clip-triggered spectrograms and plot
% the panels presented in figure 4 (HFB vs ripple rate).
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

blocks = {'R1','R2','M1','M2','MT'};

set(0,'DefaultAxesFontName', 'Arial')
ref_flag = 2; % 1 = Common Ref; 2 = Bipolar montage;

figdir=fullfile(parentfolder,'results','clip_triggered_spectrograms','group-level results',['ref_' num2str(ref_flag)]);
if ~exist(figdir,'dir')
    mkdir(figdir);
    disp('Creating Output Directory...')
end


fprintf('\nloading the time-frequency data\n');
% Set data structure:
DATA = struct;
DATA.HFB1 = [];
DATA.HFB2 = [];
DATA.S1_norm = [];
DATA.S2_norm = [];
chINFO = {};
for iSub=1:length(subjects)
    
    clearvars -except subjects blocks iSub outdir path_to_toolboxes parentfolder DATAGRP PARAMSGRP chINFO figdir ref_flag
    subjid = subjects{iSub};
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    
    % load datasets:
    maindir=fullfile(parentfolder,subjid);
    datasetID='pink_panther';
    
    % Set outdir:
    indir=fullfile(parentfolder,'results','data','clip_triggered_spectrograms',['ref_' num2str(ref_flag)]);

    
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0);
    load(fullfile(indir,sprintf('%s_clip_spectrogram_DATA.mat',subjid)));
    load(fullfile(indir,'clip_spectrogram_parameters.mat'));
    load(fullfile(indir,sprintf('%s_electrodes_table.mat',subjid)));
    
    if iSub==1, PARAMSGRP = params_plot;
    else, PARAMSGRP = CatStructFields(PARAMSGRP,params_plot,2); end
   
    if iSub==1, DATAGRP = DATA;
    else, DATAGRP = CatStructFields(DATAGRP,DATA,2); end
    
    if iSub==1, chINFO = electrodes_table; 
    else,chINFO = cat(1,chINFO,electrodes_table); end
    fprintf('\nloading subject %s ... O.K.',subjid);
    
end

POS = [0 0 200 160];

%% PLOT Group level results:
set_figure_colors;
close all;
saveflag = 1;
counter = 0;

CM = params_plot.CM;
params_plot.clim = [-2 2];
params_plot.ch_label = '';
params_plot.stimulus_onset = 0;
params_plot.xlim = [-1 6];
params_plot.fig_position = POS;
params_plot.logscaleflag = 0;

% Figure 1 - OLD:
SS =cellfun(@(x)nanmean(x,3),DATAGRP.normspec1,'UniformOutput',0);
avgspec = nanmean(cat(3,SS{:}),3);
params_plot.title = sprintf('OLD');
params_plot.name = sprintf('OLD clip spectrogram Group Level');
[H1,h1,hcb1] = plot_spectrogram(avgspec,params_plot);
axes(h1); tmpx=h1.XLim; tmpy=h1.YLim; h1.XLabel.String = 'Time from clip onset (s)';
text(tmpx(2),tmpy(2)-0.15*range(tmpy),sprintf('n=%d',length(SS)),'HorizontalAlignment','right','fontsize',6);
set_font_size_and_type;
if saveflag, save_current_figure(H1,figdir,1); end

% Figure 1 - New:
SS =cellfun(@(x)nanmean(x,3),DATAGRP.normspec2,'UniformOutput',0);
avgspec = nanmean(cat(3,SS{:}),3);
params_plot.title = sprintf('NEW');
params_plot.name = sprintf('NEW clip spectrogram Group Level');
[H2,h2,hcb2] = plot_spectrogram(avgspec,params_plot);
axes(h2); tmpx=h2.XLim; tmpy=h2.YLim; h2.XLabel.String = 'Time from clip onset (s)';
text(tmpx(2),tmpy(2)-0.15*range(tmpy),sprintf('n=%d',length(SS)),'HorizontalAlignment','right','fontsize',6);
set_font_size_and_type;
if saveflag, save_current_figure(H2,figdir,1); end

% Figure 1 - Diff:
SS1 =cellfun(@(x)nanmean(x,3),DATAGRP.normspec1,'UniformOutput',0);
SS2 =cellfun(@(x)nanmean(x,3),DATAGRP.normspec2,'UniformOutput',0);
avgspec = nanmean(cat(3,SS1{:})-cat(3,SS2{:}),3);
params_plot.title = sprintf('OLD-NEW');
params_plot.name = sprintf('OLD-NEW clip spectrogram Group Level');
[H3,h3,hcb3] = plot_spectrogram(avgspec,params_plot);
axes(h3); tmpx=h3.XLim; tmpy=h3.YLim; h3.XLabel.String = 'Time from clip onset (s)';
text(tmpx(2),tmpy(2)-0.15*range(tmpy),sprintf('n=%d',length(SS)),'HorizontalAlignment','right','fontsize',6);
set_font_size_and_type;
if saveflag, save_current_figure(H3,figdir,1); end


%% Stats:
% ================
% Run stats on HFB: 
t_out_sec = t_out./1000;
test_interval = t_out_sec>params_plot.xlim(1) & t_out_sec<params_plot.xlim(2); % -2 to 6 s
HFB1 = cellfun(@(x)nanmean(x,2),DATAGRP.HFB1,'UniformOutput',0);
HFB2 = cellfun(@(x)nanmean(x,2),DATAGRP.HFB2,'UniformOutput',0);
stat1 = cat(2,HFB1{:});
stat2 = cat(2,HFB2{:});
% group by patients:
s = unique(chINFO.subjid);
stat1grp = []; stat2grp = [];
for iSub = 1:length(s)
    stat1grp = cat(2,stat1grp,mean(stat1(test_interval,strcmpi(chINFO.subjid,s{iSub})),2));
    stat2grp = cat(2,stat2grp,mean(stat2(test_interval,strcmpi(chINFO.subjid,s{iSub})),2));
end

% Draw HFB response:
H4 = figure('color','w','name',sprintf('HFB OLD-vs-NEW group-level'),'position',[0 0 220 180]); hold on;
title(sprintf('High-gamma (60-128Hz)'),'fontweight','normal');
h1 = shadedErrorBar(t_out_sec(test_interval),mean(stat1grp,2),std(stat1grp,[],2)./sqrt(size(stat1grp,2)),{'color',COLOR.blue,'linewidth',0.5},0.5)
h2 = shadedErrorBar(t_out_sec(test_interval),mean(stat2grp,2),std(stat2grp,[],2)./sqrt(size(stat2grp,2)),{'color',COLOR.red,'linewidth',0.5},0.5)
xlabel('Time from clip onset (s)'); ylabel('Power (dB)'); 
axis tight; ylim([-0.5 1]); xlim(params_plot.xlim);
tmpx=get(gca,'XLim'); tmpy=get(gca,'YLim'); 
text(tmpx(2),tmpy(2)-0.15*range(tmpy),sprintf('n=%d',length(s)),'HorizontalAlignment','right','fontsize',6);

% compute stats: (cluster-based permutation test)
[sig, pvalHFB, onsets, offsets] = compare_conditions(gca,stat1grp,stat2grp,t_out_sec(test_interval),1);
if any(sig),text(tmpx(2),tmpy(2)-0.25*range(tmpy),sprintf('%.2f s',onsets(1)),'HorizontalAlignment','right','fontsize',6); end

% add legends (optional):
[h_leg,h_leg_icon] = legend([h1.mainLine,h2.mainLine],{'OLD','NEW'},'location','NorthWest');
pos = h_leg.Position;
h_leg.Position = [pos(1)*0.5 0.75 pos(3:4)];
h = findobj(h_leg_icon,'type','line');
for ii = 1:length(h), if length(h(ii).XData)==2, h(ii).XData(1) = h(ii).XData(1) + 0.5*range(h(ii).XData); end; end
legend boxoff
set_font_size_and_type;
%enlarge_figure_and_move_axes_to_center(gcf,gca,1.1);
if saveflag, save_current_figure(H4,figdir,1); end


% =========================================================
%% Draw HFB response: clips vs orig events from the full-length cartoon (suppl. analysis)
t_out_sec = t_out./1000;
test_interval = t_out_sec>params_plot.xlim(1) & t_out_sec<params_plot.xlim(2); % -2 to 6 s

HFB1 = cellfun(@(x)nanmean(x,2),DATAGRP.HFB1,'UniformOutput',0);
HFB3 = cellfun(@(x)nanmean(x,2),DATAGRP.HFB3,'UniformOutput',0);
HFB4 = cellfun(@(x)nanmean(x,2),DATAGRP.HFB4,'UniformOutput',0);
stat1 = cat(2,HFB1{:});
stat3 = cat(2,HFB3{:});
stat4 = cat(2,HFB4{:});
% group by patients:
s = unique(chINFO.subjid);
stat1grp = []; stat3grp = []; stat4grp = [];
for iSub = 1:length(s)
    stat1grp = cat(2,stat1grp,mean(stat1(test_interval,strcmpi(chINFO.subjid,s{iSub})),2));
    stat3grp = cat(2,stat3grp,mean(stat3(test_interval,strcmpi(chINFO.subjid,s{iSub})),2));
    stat4grp = cat(2,stat4grp,mean(stat4(test_interval,strcmpi(chINFO.subjid,s{iSub})),2));
    stat5grp = (stat3grp+stat4grp)./2;
end

H6 = figure('color','w','name',sprintf('HFB clip-vs-orig group-level'),'position',[0 0 220 180]); hold on;
title(sprintf('High-gamma (60-128Hz)'),'fontweight','normal');

h1 = shadedErrorBar(t_out_sec(test_interval),mean(stat1grp,2),std(stat1grp,[],2)./sqrt(size(stat1grp,2)),{'color',COLOR.blue,'linewidth',0.5},0.5);
h2 = shadedErrorBar(t_out_sec(test_interval),mean(stat3grp,2),std(stat3grp,[],2)./sqrt(size(stat3grp,2)),{'color',COLOR.greendark,'linewidth',0.5},0.5);
h3 = shadedErrorBar(t_out_sec(test_interval),mean(stat4grp,2),std(stat4grp,[],2)./sqrt(size(stat4grp,2)),{'color',COLOR.turquoise,'linewidth',0.5},0.5);

xlabel('Time from clip onset (s)'); ylabel('Power (dB)'); 
axis tight; ylim([-0.5 1]); xlim(params_plot.xlim);
tmpx=get(gca,'XLim'); tmpy=get(gca,'YLim'); 
text(tmpx(2),tmpy(2)-0.15*range(tmpy),sprintf('n=%d',length(s)),'HorizontalAlignment','right','fontsize',6);


% compute stats:
[sig, pval, onsets, offsets] = compare_conditions(gca,stat1grp,stat4grp,t_out_sec(test_interval),1);
if any(sig),text(tmpx(2),tmpy(2)-0.25*range(tmpy),sprintf('%.2f s',onsets(1)),'HorizontalAlignment','right','fontsize',6); end

% add legends (optional):
[h_leg,h_leg_icon] = legend([h1.mainLine,h2.mainLine,h2.mainLine],{'recognized OLD','1st presentation','2nd presentation'},'location','NorthWest');
pos = h_leg.Position;
h_leg.Position = [pos(1)*0.5 0.75 pos(3:4)];
h = findobj(h_leg_icon,'type','line');
for ii = 1:length(h), if length(h(ii).XData)==2, h(ii).XData(1) = h(ii).XData(1) + 0.5*range(h(ii).XData); end; end
legend boxoff
set_font_size_and_type;
%enlarge_figure_and_move_axes_to_center(gcf,gca,1.1);
if saveflag, save_current_figure(H6,figdir,0); end


%% load ripple rate:
% Set Manually:
ref_flag=2; % 1 = Common Ref; 2 = Bipolar montage;
time_locking_event = 'stimonset'; % time locking event (stimonset/rt)
% Load rasters and initalize general varaiables:
load_swr_rasters_from_all_hippocampal_channels;


%% Single trial RIPPLE RATE:
% compute overal ripple PSTH to find the latency of overall ripple rate peak (during video clip presentation):
ind = contains(DATA.Rall.stimulus_type,'Clip');

% set PETH parameters:
binsize = 50;
smt = 5;
[r,binscenters]=nanPSTH(logical(DATA.Rall.raster(ind,:))',binsize,smt,DATA.Fs, []);
binscenters = (binscenters/1000)+DATA.trialtime(1);
[~,pk] = max(r);

% plot the response window:
H0 = figure('color','w','name','peak SWR response','position',[0 0 200 200]); hold on;
plot(binscenters,r,'linewidth',2','color',COLOR.red);
set(gca,'ytick',0.3:0.1:1)
tmpy=get(gca,'ylim'); tmpx=get(gca,'xlim');
line([0 4],[tmpy(1) tmpy(1)]+0.01*range(tmpy),'color',COLOR.black,'linewidth',5)
scatter(binscenters(pk),r(pk),50,'ok','linewidth',2);
switch time_locking_event, case 'stimonset', xlim([-2 6]); case 'rt', xlim([-6 2]); end
xlabel('Time from clip onset (s)'); ylabel('Ripple rate (events/s)')
fprintf('\nLatency of peak SWR rate during the video clip: %.2f sec\n',binscenters(pk));

% set the response time window to +/- 1 sec around the peak:
%load(fullfile(parentfolder,'results','stats','sigTimeWindow.mat'));
timeWindow = [binscenters(pk)-1 binscenters(pk)+1];
%timeWindow = [0 2];
tmpy = get(gca,'ylim');

areaHandel = area(timeWindow,[tmpy(2),tmpy(2)],'BaseValue',tmpy(1),'facecolor',COLOR.lightgray,'facealpha',0.5,'edgecolor','none');
uistack(areaHandel,'bottom')
ylim(tmpy);

%% =========================================================================
% Now compute the single trial ripple rate:
[RRtw,trialDur] = compute_single_trial_ripple_rate(DATA,Fs,timeWindow,time_locking_event); % using a fixed time-window centered on ripple rate peak
[RR4s,~] = compute_single_trial_ripple_rate(DATA,Fs,4,time_locking_event);                 % using a fixed time-window of 4s from event onset (for rest epochs)

% normalize within electrode:
norm_flag = 0;
L1 = unique(DATA.Rall.channelid,'stable');
RRnormTW = nan(size(RRtw)); %   single trial ripple rate within the time window of analysis (for the memory test)
RRnorm4s = nan(size(RR4s)); %   single trial ripple rate within a 4 time window (entire epoch length for rest and movie conditions)
for ii=1:length(L1)
    elec = L1(ii);
    k = strcmpi(DATA.Rall.channelid,elec);
    blavg = nanmean(RR4s(k & contains(DATA.Rall.stimulus_type,'R')));
    blstd = nanstd(RR4s(k & contains(DATA.Rall.stimulus_type,'R')));
    switch norm_flag
        case 0 % raw (events/s)
            RRnormTW(k) = RRtw(k);
            RRnorm4s(k) =  RR4s(k);
        case 1 % zscore
            RRnormTW(k) = (RRtw(k)-blavg)./blstd;
            RRnorm4s(k) = (RR4s(k)-blavg)./blstd;
        case 2 % min-max
            blmin = min(RRtw(k)); blmax = max(RRtw(k));
            RRnormTW(k) = (RRtw(k)-blmin)./(blmax-blmin);
            RRnorm4s(k) = (RR4s(k)-blavg)./blstd;
    end
    disp(elec);
end

indClipA = contains(DATA.Rall.stimulus_type,'ClipA') & ~contains(DATA.Rall.response,'new'); % correct old responses
indClipB = contains(DATA.Rall.stimulus_type,'ClipB') & ~contains(DATA.Rall.response,'old'); % correct new responses

[m1,s1,elecName1,n1] = grpstats(RRnorm4s(indClipA),DATA.Rall.channelid(indClipA),{'mean','var','gname','numel'}); % NEW
[m2,s2,elecName2,n2] = grpstats(RRnorm4s(indClipB),DATA.Rall.channelid(indClipB),{'mean','var','gname','numel'}); % OLD
% Compute hedge's G: (OLD-vs-NEW)
df=n1+n2-2;
sP=((n1-1).*s1 + (n2-1).*s2)./(n1+n2-2); % pooled (within-groups) variance
RRes=(m1-m2)./sqrt(sP);

%% Get stats on HFB:
t_out_sec = t_out./1000;
%test_interval = t_out_sec>timeWindow(1) & t_out_sec<=timeWindow(2); % trial response
test_interval = t_out_sec>0 & t_out_sec<=4; % trial response

HFBch = strcat(chINFO.channel,'_',chINFO.subjid);
HFB1 = cellfun(@(x)nanmean(x,2),DATAGRP.HFB1,'UniformOutput',0);
HFB2 = cellfun(@(x)nanmean(x,2),DATAGRP.HFB2,'UniformOutput',0);
HFBs1 = cellfun(@(x)nanvar(x,[],2),DATAGRP.HFB1,'UniformOutput',0);
HFBs2 = cellfun(@(x)nanvar(x,[],2),DATAGRP.HFB2,'UniformOutput',0);
HFBn1 = cellfun(@(x)size(x,2),DATAGRP.HFB1,'UniformOutput',0);
HFBn2 = cellfun(@(x)size(x,2),DATAGRP.HFB2,'UniformOutput',0);

% rearange in matrices:
HFBm1 = cat(2,HFB1{:}); HFBm1 = mean(HFBm1(test_interval,:),1);
HFBm2 = cat(2,HFB2{:}); HFBm2 = mean(HFBm2(test_interval,:),1);
HFBs1 = cat(2,HFBs1{:}); HFBs1 = mean(HFBs1(test_interval,:),1);
HFBs2 = cat(2,HFBs2{:}); HFBs2 = mean(HFBs2(test_interval,:),1);
HFBn1 = cat(2,HFBn1{:});
HFBn2 = cat(2,HFBn2{:}); 

% compute hedge's G: (OLD-vs-NEW) on HFB:
df=HFBn1+HFBn2-2;
sP=((HFBn1-1).*HFBs1 + (HFBn2-1).*HFBs2)./(HFBn1+HFBn2-2); % pooled (within-groups) variance
HFBes=(HFBm1-HFBm2)./sqrt(sP);

% sort the electrodes:
[ElecID,HFBch_sortix] = sort(HFBch);
[~,RRch_sortix] = sort(elecName1);
% make sure the order of electrodes corresponds:
assert(all(strcmpi(elecName1(RRch_sortix),HFBch(HFBch_sortix))))


%% Scatter: HFB-vs-ripple rate
H7 = figure('color','w','name',sprintf('HFB-vs-RR group-level'),'position',[0 0 220 180]); hold on;
title(sprintf('Effect size (OLD>NEW)'),'fontweight','normal');

scatter(HFBes(HFBch_sortix),RRes(RRch_sortix),30,'k.');
axis tight 
%axis([-0.5 1.5 -0.5 1.5])
xlabel('High gamma')
ylabel('Ripple rate')
tmpx=get(gca,'XLim'); tmpy=get(gca,'YLim'); 

[r,pcorr]= corr(HFBes(HFBch_sortix)',RRes(RRch_sortix));

% mixed effet analysis comparing RR and HFB:
[s] = split(ElecID,'_'); s=s(:,2);
tbl = table(double(RRes(RRch_sortix)-HFBes(HFBch_sortix)'),s,ElecID,'VariableNames',{'effectSize','subject','electrode'});
lmeES = fitlme(tbl,'effectSize~1+(1|subject:electrode)');
res = lmeES.Coefficients;
%text(tmpx(2),tmpy(2)-0.15*range(tmpy),sprintf('t(%d)=%.2f, P<%.3f',res.DF,res.tStat,fix(res.pValue*1000)/1000),'HorizontalAlignment','right','fontsize',6);
text(tmpx(2),tmpy(2)-0.05*range(tmpy),sprintf('Rsquare=%.2f, P=%.4f',r^2,pcorr),'HorizontalAlignment','right','fontsize',6);
h_ls = lsline; h_ls.Color=[0.5 0.5 0.5];
%plot([-0.5 1.5], [-0.5,1.5],'k:')
%hl=lsline; set(hl,'color',COLOR.red,'linestyle',':')
axis square
% tmpx = get(gca,'xlim'); tmpy = get(gca,'ylim');
% mn = min([xlim, ylim]); mx = max([xlim,ylim]);
% axis([mn,mx,mn,mx]);
% hold on; plot([mn mx], [mn mx],'k-');
if saveflag, save_current_figure(H7,figdir,0); end


%% histograms:
H7 = figure('color','w','name',sprintf('HFB-vs-RR group-level'),'position',[0 0 220 180]); hold on;
title(sprintf('Effect size histogram (HFB-vs-RR)'),'fontweight','normal');
be = -0.5:0.1:1;
histogram(HFBes(HFBch_sortix),'BinEdges',be,'Normalization','probability'); histogram(RRes(RRch_sortix),'BinEdges',be,'Normalization','probability')


%% Line plot 1: HFB-vs-ripple rate mean

d = [HFBes(HFBch_sortix)',RRes(RRch_sortix)];
ss = split(elecName1(RRch_sortix),'_');

d_mean = grpstats(d, ss(:,2),'mean');
d_max = grpstats(d, ss(:,2),'max');
pval1 = signrank(d_mean(:,1),d_mean(:,2));
pval2 = signrank(d_max(:,1),d_max(:,2));

H8 = figure('color','w','name',sprintf('HFB-vs-RR group-level line plot mean'),'position',[0 0 220 180]); hold on;
title(sprintf('Effect size (OLD>NEW)'),'fontweight','normal');

hold on; 
plot([ones(size(d_mean,1),1),ones(size(d_mean,1),1)*2]',d_mean','ok-','markersize',3,'markerfacecolor',COLOR.lightgray)
title(sprintf('P=%.4f',pval1),'fontweight','normal')
set(gca, 'xlim',[0.5 2.5],'xtick',[1,2],'ylim',[-0.3 1],'ytick',[-2:0.5:2],'xticklabels',{'HFB','R.Rate'})
ylabel('Old>New, effect size (g)');
axis square
if saveflag, save_current_figure(H8,figdir,0); end

H9 = figure('color','w','name',sprintf('HFB-vs-RR group-level line plot max'),'position',[0 0 220 180]); hold on;
title(sprintf('Effect size (OLD>NEW)'),'fontweight','normal');

hold on; 
plot([ones(size(d_max,1),1),ones(size(d_max,1),1)*2]',d_max','ok-','markersize',3,'markerfacecolor',COLOR.lightgray)
title(sprintf('P=%.4f',pval2),'fontweight','normal')
set(gca, 'xlim',[0.5 2.5],'xtick',[1,2],'ylim',[-0.3 1.5],'ytick',[-2:0.5:2],'xticklabels',{'HFB','R.Rate'})
ylabel('Old>New, effect size (g)');
axis square
if saveflag, save_current_figure(H9,figdir,0); end
