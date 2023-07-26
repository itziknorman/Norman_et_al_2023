% The following script loads the multitaper spectra computed in
% "compute_spectra.m" and plot the grand average hippocampal spectra in
% each condition (as a supplementary analysis, currently not included)
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


DATA = struct;
% Set data structure:
DATA.S0 = []; % REST 1
DATA.S1 = []; % MOVIE 1
DATA.S2 = []; % REST 2
DATA.S3 = []; % MOVIE 2
DATA.S4 = []; % MEMORY TEST
DATA.S0_norm = [];
DATA.S1_norm = [];
DATA.S2_norm = [];
DATA.S3_norm = [];
DATA.S4_norm = [];

chINFO = {};
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
    
    % Set outdir:
    indir=fullfile(parentfolder,'results','data','multitaper-spectra',['spectra_data_ref_' num2str(ref_flag)]);
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0);
    load(fullfile(indir,sprintf('%s_electrode_details.mat',subjid)));
    tmp = load(fullfile(indir,sprintf('%s_spectra_data.mat',subjid)));
    DATA = CatStructFields(DATA,tmp,2);
    chINFO = cat(1,chINFO,ch_list');
    
end

outdir=fullfile(parentfolder,'results','multitaper-spectra-group',['spectra_data_ref_' num2str(ref_flag)]);
if ~exist(outdir,'dir')
    mkdir(outdir);
    disp('Creating Output Directory...')
end
set_figure_colors;
close all;
%% load electrodes anatomical info:
L1 = chINFO;
[elecinfo,Hmap] = getHippocampalElectrodeCoord(L1, subjects,0,1);

%% get subfield/hemi info:
ElecID = {};
ElecAnatSF = {};
ElecHemisphere = {};
for ii = 1:length(chINFO)
    tmp = split(chINFO(ii),'_');
    curElecName = tmp(2);
    subjectName = tmp(1);
    curElecName = split(curElecName,'-'); curElecName = curElecName(1);
    subjdir=fullfile(parentfolder, string(subjectName));
    fid = fopen(fullfile(subjdir,'subjectname')); C = textscan(fid,'%s'); fclose(fid);
    FSsubjid = cell2mat(C{1}); % subject FS code
    currentSubElec = find(contains(elecinfo.groupLabels ,FSsubjid) & contains(elecinfo.groupLabels ,curElecName));
    
    [vals,ind] = min(elecinfo.groupDist2sf{currentSubElec,:},[],2);
    if ~elecinfo.groupIsSFparc(currentSubElec) % exclude electrodes that don't have hippocampal parcellation
        ElecAnatSF = [ElecAnatSF; {'N/A'}];
    else
        ElecAnatSF = [ElecAnatSF; elecinfo.groupDist2sf.Properties.VariableNames(ind)];
    end
    ElecHemisphere = [ElecHemisphere; num2str(elecinfo.groupIsLeft(currentSubElec))];
    ElecID = [ElecID; chINFO(ii)];
end


leftInd = contains(ElecHemisphere,'1');
rightInd = contains(ElecHemisphere,'0');


%% plot:
set_figure_colors_PP;
load(fullfile(indir,'spectra_parameters.mat'));
saveflag = 0;

f1 = find(abs(f_out-50)<1.5);

excluded_freq = [f1];
f_out_all = f_out;

f_out(excluded_freq) = nan;
W = 10;
avg0 = nanmean(W*log10(DATA.S0),2); avg0(excluded_freq) = nan;
avg1 = nanmean(W*log10(DATA.S1),2); avg1(excluded_freq) = nan;
avg2 = nanmean(W*log10(DATA.S2),2); avg2(excluded_freq) = nan;
avg3 = nanmean(W*log10(DATA.S3),2); avg3(excluded_freq) = nan;
avg4 = nanmean(W*log10(DATA.S4),2); avg4(excluded_freq) = nan;
sem0 = nanstd(W*log10(DATA.S0),[],2)./sqrt(size(DATA.S0,2)); sem0(excluded_freq) = nan;
sem1 = nanstd(W*log10(DATA.S1),[],2)./sqrt(size(DATA.S1,2)); sem1(excluded_freq) = nan;
sem2 = nanstd(W*log10(DATA.S2),[],2)./sqrt(size(DATA.S2,2)); sem2(excluded_freq) = nan;
sem3 = nanstd(W*log10(DATA.S3),[],2)./sqrt(size(DATA.S3,2)); sem3(excluded_freq) = nan;
sem4 = nanstd(W*log10(DATA.S4),[],2)./sqrt(size(DATA.S4,2)); sem4(excluded_freq) = nan;


figure('color','w','name','raw spectra','position',[0 0 250 200]);  hold on;
ix = find(f_out>0 & f_out<=40);
h0 = plot(log10(f_out(ix)),avg0(ix),'color',COLOR.rest1);
h1 = plot(log10(f_out(ix)),avg1(ix),'color',COLOR.movie1);
h2 = plot(log10(f_out(ix)),avg2(ix),'color',COLOR.rest2);
h3 = plot(log10(f_out(ix)),avg3(ix),'color',COLOR.movie2);
%h4 = plot(log10(f_out(ix)),avg4(ix),'color',COLOR.memorytest);
xlabel('Frequency'); ylabel({'Power sepctral density';'10*log_{10}(\muV^{2}/Hz)'})
set(gca,'xtick',[log10(4) log10(10) log10(20) log10(40) 2],'xticklabel',{'4','10','20','40','100'});
xlim([0 log10(30)])
legend([h1 h3 h0 h2],{'1st screening','2nd screening','Rest1','Rest2'},'location','northeast'); legend boxoff;
set_font_size_and_type;
enlarge_figure_and_move_axes_to_center(gcf,gca,1.1);
if saveflag, save_current_figure(gcf,outdir,0); end


%% Relative power:
% stats:
[s] = split(chINFO,'_'); s=s(:,1);
p1 = [];
p2 = [];
p3 = [];

nbins = 128;
testix = find(f_out>0 & f_out<=128);
[~,f_out_bin,binix] = histcounts(f_out(testix),nbins);
f_out_bin = f_out_bin(2:end);
% fit mixed effects model to each frequency bin (log-spaced):
for ii = unique(binix)
    jj = find(binix==ii);
    pow = log10(nangeomean([DATA.S1_norm(testix(jj),:)'; DATA.S2_norm(testix(jj),:)';  DATA.S3_norm(testix(jj),:)'],2));
    cond = [ones(size(s)); ones(size(s))*2; ones(size(s))*3];
    tbl = table(pow,cond,repmat(leftInd,[3,1]),repmat(s,[3,1]),repmat(chINFO,[3,1]),'VariableNames',{'power','cond','hemi','subject','electrode'});
    lmeHemi = fitlme(tbl,'power ~ 1 + cond + (1|subject:electrode)');
    res1 = lmeHemi.Coefficients;
    p1(ii,1) = lmeHemi.anova.pValue(1);  
    disp(f_out(testix(jj(1))))
end

[~, crit_p1, ~, adj_p1]=fdr_bh(p1,0.05,'pdep','yes');

 
%%
for hemi = 0
    switch hemi
        case 0,  ch_ix = true(size(leftInd));  elecgrp = 'all';
        case 1,  ch_ix = leftInd;  elecgrp = 'LH';
        case 2,  ch_ix = rightInd; elecgrp = 'RH';
    end
    sig = adj_p1<0.05;
    
    % compute the mean speactra:
    avg0 = nanmean(W.*log10(DATA.S0_norm(:,ch_ix)),2); avg0(excluded_freq) = nan; % REST 1
    avg1 = nanmean(W.*log10(DATA.S1_norm(:,ch_ix)),2); avg1(excluded_freq) = nan; % MOVIE 1
    avg2 = nanmean(W.*log10(DATA.S2_norm(:,ch_ix)),2); avg2(excluded_freq) = nan; % REST 2
    avg3 = nanmean(W.*log10(DATA.S3_norm(:,ch_ix)),2); avg3(excluded_freq) = nan; % MOVIE 2
    avg4 = nanmean(W.*log10(DATA.S4_norm(:,ch_ix)),2); avg4(excluded_freq) = nan; % MEMORY TEST
    sem0 = nanstd(W.*log10(DATA.S0_norm(:,ch_ix)),[],2)./sqrt(sum(ch_ix)); sem0(excluded_freq) = nan;
    sem1 = nanstd(W.*log10(DATA.S1_norm(:,ch_ix)),[],2)./sqrt(sum(ch_ix)); sem1(excluded_freq) = nan;
    sem2 = nanstd(W.*log10(DATA.S2_norm(:,ch_ix)),[],2)./sqrt(sum(ch_ix)); sem2(excluded_freq) = nan;
    sem3 = nanstd(W.*log10(DATA.S3_norm(:,ch_ix)),[],2)./sqrt(sum(ch_ix)); sem3(excluded_freq) = nan;
    sem4 = nanstd(W.*log10(DATA.S4_norm(:,ch_ix)),[],2)./sqrt(sum(ch_ix)); sem4(excluded_freq) = nan;
    
    
    %set(gca,'XScale','log')
    ix = find(f_out>0 & f_out<=128);
    figure('color','w','name',sprintf('normalized spectra %s',elecgrp),'position',[300 0 200 200]); hold on;
    h1 = shadedErrorBar(log10(f_out(ix)),avg1(ix),sem1(ix),{'color',COLOR.movie1,'linewidth',0.5},1);
    h2 = shadedErrorBar(log10(f_out(ix)),avg2(ix),sem2(ix),{'color',COLOR.gray,'linewidth',0.5},1);
    h3 = shadedErrorBar(log10(f_out(ix)),avg3(ix),sem3(ix),{'color',COLOR.movie2,'linewidth',0.5},1);
    h4 = shadedErrorBar(log10(f_out(ix)),avg4(ix),sem4(ix),{'color',COLOR.blue,'linewidth',0.5},1);
  
%     h1 = plot(log10(f_out(ix)),avg1(ix),'color',COLOR.movie1);
%     h2 = plot(log10(f_out(ix)),avg2(ix),'color',COLOR.rest2);
%     h3 = plot(log10(f_out(ix)),avg3(ix),'color',COLOR.movie2);

    lowfreqix = find(f_out<20);
    [mnval,mnind] = min(avg1(lowfreqix));
    mntick = log10(f_out(lowfreqix(mnind)));
    mntickl = sprintf('%.1f',f_out(lowfreqix(mnind)));
    plot([mntick,mntick],get(gca,'ylim'),'k--')
    set(h1.patch,'FaceAlpha',0.25);
    set(h2.patch,'FaceAlpha',0.25);
    set(h3.patch,'FaceAlpha',0.25);
    set(h4.patch,'FaceAlpha',0.25);
    set(gca,'xtick',[0 log10(4) mntick log10(20) log10(40) 2],'xticklabel',{'1','4',mntickl,'20','40','100'});
    xlim([0 log10(120)])
    xlabel('Frequency'); ylabel('Relative power change (dB)')
    
    % plot significance:
    if any(sig)
        xx = log10(f_out_bin);
        idx=sprintf('%d',sig);
        onsets = xx(regexp(idx, '1{1,}', 'start'));    offsets = xx(regexp(idx, '1{1,}', 'end'));
        tmp = get(gca,'Ylim');
        X=[onsets;offsets];
        Y=repmat(tmp(1)+range(tmp)*0.02,size(X));
        line(X,Y,'LineWidth',4,'Color',COLOR.orange)
        %text(mean(X),mean(Y)+0.1,sprintf('p<0.05'),'fontsize',6,'horizontalalignment','center','VerticalAlignment','bottom');
    end
    
    
    
     [h_leg,h_leg_icon] = legend([h1.mainLine h3.mainLine h2.mainLine h4.mainLine],{'1st screening','2nd screening','Rest II','Memory Test'},'location','southwest'); legend boxoff;
         
     pos = h_leg.Position;
     h_leg.Position = [pos(1)*0.5 0.7 pos(3:4)];
     h = findobj(h_leg_icon,'type','line');
     for ii = 1:length(h), if length(h(ii).XData)==2, h(ii).XData(1) = h(ii).XData(1) + 0.5*range(h(ii).XData); end; end
     legend boxoff

    set_font_size_and_type;
    
    if saveflag, save_current_figure(gcf,outdir,0); end
    
    
end


