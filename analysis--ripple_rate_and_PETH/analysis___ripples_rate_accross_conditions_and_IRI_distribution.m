% The script computes  basal ripple rate during each experimental condition 
% and runs mixed effects models to compare ripple rate, duration and amp across conditions.

%% 
% =========================================================================
clear all
close all
clc;
% set the relevant path:
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;

% add/remove paths:
addpath(fullfile(path_to_toolboxes,'eeglab2021.1'));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));
addpath(genpath(fullfile(path_to_toolboxes,'dmgroppe-Mass_Univariate_ERP_Toolbox-d1e60d4')));

[ALLEEG, EEG, CURRENTSET] = eeglab;

%==========================================================================
% set path:
outdir=fullfile('parentfolder','results','mean_ripple_rate');
if ~exist(outdir,'dir')
    mkdir(outdir);
    disp('Creating Output Directory...')
end
%==========================================================================

% load eeglab:
[ALLEEG, EEG, CURRENTSET] = eeglab;
saveflag = 0;
ref_flag = 2;
ripplesdir = fullfile(parentfolder,'results','data','Ripple_times_hamming_2-4std_adjusted_band_20ms_30ms');
figdir = fullfile(parentfolder,'results','basal_ripple_rate');
if ~exist(figdir,'dir'), mkdir(figdir); end

%% Load the raw iEEG datasets:
subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};

ALLEEG=[]; EEG=[]; CURRENTSET=1;
allR1=[]; allR2=[]; allM1=[]; allM2=[]; allMT=[];
for i = 1:numel(subjects)
    subjid=subjects{i};
    clearvars -except i subjects subjid ALLEEG EEG CURRENTSET ripplesdir parentfolder path_to_toolboxes ...
        allM1 allM2 allMT allR1 allR2 outdir ref_flag figdir
    maindir = fullfile(parentfolder,subjid);
    
    switch ref_flag
        case 1
            datadir=fullfile(maindir,'EEGLAB_datasets');
        case 2
            datadir=fullfile(maindir,'EEGLAB_datasets_BP');
    end
    figdir
    % Load eeglab datasets:
    filenames=dir(fullfile(datadir,['*pink_panther_preprocessed*MT.set']));
    for SET=1:numel(filenames)
        filename=fullfile(datadir,filenames(SET).name);
        [EEG] = pop_loadset('filename', filename);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
        eeglab redraw;
    end
    
    % load ripple timings:
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjects{i},ref_flag,0);
    
    for hp = 1:length(hippocampus_all_channels)
        currenthippocampus = hippocampus_all_channels{hp};
        M1=load(fullfile(ripplesdir,subjid,[subjid ' M1 ripples ' currenthippocampus '.mat']));
        M2=load(fullfile(ripplesdir,subjid,[subjid ' M2 ripples ' currenthippocampus '.mat']));
        MT=load(fullfile(ripplesdir,subjid,[subjid ' MT ripples ' currenthippocampus '.mat']));
        R1=load(fullfile(ripplesdir,subjid,[subjid ' R1 ripples ' currenthippocampus '.mat']));
        R2=load(fullfile(ripplesdir,subjid,[subjid ' R2 ripples ' currenthippocampus '.mat']));
        % Add field: subject id
        M1.ripples.subjid = repmat({subjid},size(M1.ripples.peak));
        M2.ripples.subjid = repmat({subjid},size(M2.ripples.peak));
        MT.ripples.subjid = repmat({subjid},size(MT.ripples.peak));
        R1.ripples.subjid = repmat({subjid},size(R1.ripples.peak));
        R2.ripples.subjid = repmat({subjid},size(R2.ripples.peak));
        % Add field: electrode id
        M1.ripples.elecid = repmat({strcat(subjid,'_',currenthippocampus)},size(M1.ripples.peak));
        M2.ripples.elecid = repmat({strcat(subjid,'_',currenthippocampus)},size(M2.ripples.peak));
        MT.ripples.elecid = repmat({strcat(subjid,'_',currenthippocampus)},size(MT.ripples.peak));
        R1.ripples.elecid = repmat({strcat(subjid,'_',currenthippocampus)},size(R1.ripples.peak));
        R2.ripples.elecid = repmat({strcat(subjid,'_',currenthippocampus)},size(R2.ripples.peak));
        
        % Concatenate data:
        allM1=cat(1,allM1,M1.ripples);  % movie 1
        allM2=cat(1,allM2,M2.ripples);  % movie 2
        allMT=cat(1,allMT,MT.ripples);  % memory test
        allR1=cat(1,allR1,R1.ripples);  % rest 1
        allR2=cat(1,allR2,R2.ripples);  % rest 2
    end
end

set_figure_colors;

% set experimental duration parametes:
movieDuration = 368.64; restDuration = 180.14;

%% Compute ripple rate across conditions:
binsize = 5; overlap = 0;
Trest = -binsize : (1/EEG.srate) : restDuration + binsize;
Tmovie = -binsize : (1/EEG.srate) : movieDuration + binsize;
% init vars:
n1=[]; n2=[]; n3=[]; n4=[]; n5=[];
rr1={}; rr2={}; rr3={}; rr4={}; 
E = unique(allR1.elecid,'stable');

for ii = 1:numel(E)
    elecid = E{ii};
    subjid = split(elecid,'_'); subjid = subjid{1};
    % compute memory test duration:
    ind = find(contains({ALLEEG.setname}',subjid) & contains({ALLEEG.setname}',{'MT'}));
    EEG = ALLEEG(ind);
    strind = round(EEG.event(strcmpi({EEG.event.type},'str')).latency);
    finind = round(EEG.event(strcmpi({EEG.event.type},'fin')).latency);
    
    tailDuration = strind./EEG.srate;
    memoryTestDuration = round(max(ALLEEG(ind).times)-min(ALLEEG(ind).times))/1000-(2*tailDuration);
    
    ind1 = strcmpi(allR1.elecid,elecid) & (allR1.peak >= 0 & allR1.peak <= restDuration);
    ind2 = strcmpi(allM1.elecid,elecid) & (allM1.peak >= 0 & allM1.peak <= movieDuration);
    ind3 = strcmpi(allR2.elecid,elecid) & (allR2.peak >= 0 & allR2.peak <= restDuration);
    ind4 = strcmpi(allM2.elecid,elecid) & (allM2.peak >= 0 & allM2.peak <= movieDuration);
    ind5 = strcmpi(allMT.elecid,elecid) & (allMT.peak >= 0 & allMT.peak <= memoryTestDuration);
    
    % compute statistics of interests:
    % mean ripple rate across conditions:
    n1(ii,1) = sum(ind1)/ restDuration; % in Hz
    n2(ii,1) = sum(ind2)/ movieDuration;
    n3(ii,1) = sum(ind3)/ restDuration;
    n4(ii,1) = sum(ind4)/ movieDuration;
    n5(ii,1) = sum(ind5)/ memoryTestDuration;
    
    % mean ripple duration across conditions:
    dur1(ii,1) = mean(allR1.fin(ind1)-allR1.str(ind1)) * 1000; % in ms
    dur2(ii,1) = mean(allM1.fin(ind2)-allM1.str(ind2)) * 1000;
    dur3(ii,1) = mean(allR2.fin(ind3)-allR2.str(ind3)) * 1000;
    dur4(ii,1) = mean(allM2.fin(ind4)-allM2.str(ind4)) * 1000;
    dur5(ii,1) = mean(allMT.fin(ind5)-allMT.str(ind5)) * 1000;
    
    % mean ripple amplitude across conditions:
    amp1(ii,1) = mean(sqrt(allR1.amplitude(ind1))); % in microV
    amp2(ii,1) = mean(sqrt(allM1.amplitude(ind2)));
    amp3(ii,1) = mean(sqrt(allR2.amplitude(ind3)));
    amp4(ii,1) = mean(sqrt(allM2.amplitude(ind4)));
    amp5(ii,1) = mean(sqrt(allMT.amplitude(ind5)));    
    
    % Compute sliding window ripple rate:
    swr_time = allR1.peak(strcmpi(allR1.elecid,elecid));
    ripplevec = zeros(size(Trest));
    for k=1:numel(swr_time), current_event = swr_time(k); [~,idx] = min(abs(Trest-current_event));  ripplevec(idx(1)) = 1; end
    %[t1 ,rr1{ii}] = bin_data(ripplevec,Trest,binsize,overlap,'rate');
    rr1{ii} = ripplevec;
    
    swr_time = allM1.peak(strcmpi(allM1.elecid,elecid));
    ripplevec = zeros(size(Tmovie));
    for k=1:numel(swr_time), current_event = swr_time(k); [~,idx] = min(abs(Tmovie-current_event));  ripplevec(idx(1)) = 1; end
    %[t2 ,rr2{ii}] = bin_data(ripplevec,Tmovie,binsize,overlap,'rate');
    rr2{ii} = ripplevec;
   
    swr_time = allR2.peak(strcmpi(allR2.elecid,elecid));
    ripplevec = zeros(size(Trest));
    for k=1:numel(swr_time), current_event = swr_time(k); [~,idx] = min(abs(Trest-current_event));  ripplevec(idx(1)) = 1; end
    %[t3 ,rr3{ii}] = bin_data(ripplevec,Trest,binsize,overlap,'rate');
    rr3{ii} = ripplevec;
   
    swr_time = allM2.peak(strcmpi(allM2.elecid,elecid));
    ripplevec = zeros(size(Tmovie));
    for k=1:numel(swr_time), current_event = swr_time(k); [~,idx] = min(abs(Tmovie-current_event));  ripplevec(idx(1)) = 1; end
    %[t4 ,rr4{ii}] = bin_data(ripplevec,Tmovie,binsize,overlap,'rate');
    rr4{ii} = ripplevec;
    disp(ii);
end

% % average within patient: (optional)
% [n1,s1] = grpstats(n1,str,{'mean','gname'});
% [str] = split(E,'_'); str=str(:,1);
% [n2,s2] = grpstats(n2,str,{'mean','gname'});
% [str] = split(E,'_'); str=str(:,1);
% [n3,s3] = grpstats(n3,str,{'mean','gname'});
% [str] = split(E,'_'); str=str(:,1);
% [n4,s4] = grpstats(n4,str,{'mean','gname'});
% [str] = split(E,'_'); str=str(:,1);
% [n5,s5] = grpstats(n5,str,{'mean','gname'});

%% compute grand average ripple duration:
avg_duration = mean([dur1,dur2,dur3,dur4,dur5],2);
std_duration = std(avg_duration); 
grandavg_duration = mean(avg_duration);
fprintf('\n Mean ripple duration: %.3f (%.3f) \n',grandavg_duration,std_duration)

%% MLE stats: (Figure 2)
saveflag = 0; 
demographics = readtable(fullfile(parentfolder,'subj_info_anon.xlsx'));

[s] = split(E,'_'); s=s(:,1);
[e] = E; %split(E,'_'); e=e(:,2);
[~,ix] = ismember(s,demographics.SubjectID(1:length(unique(s)))); ix(ix==0) = []; assert(all(strcmpi(demographics.SubjectID(ix),s)));
age = demographics.Age(ix); % patient age

cond = [ones(size(n1)); 2*ones(size(n1)); 3*ones(size(n1)); 4*ones(size(n1)); 5*ones(size(n1));]
tbl = table([n1;n2;n3;n4;n5],cond,repmat(age,[5,1]),repmat(s,[5,1]),repmat(e,[5,1]),'VariableNames',{'rippleRate','condition','age','subject','electrode'});
lmeRate = fitlme(tbl,'rippleRate~condition+age+(1|subject:electrode)');
tbl = table([dur1;dur2;dur3;dur4;dur5],cond,repmat(age,[5,1]),repmat(s,[5,1]),repmat(e,[5,1]),'VariableNames',{'rippleDuration','condition','age','subject','electrode'});
lmeDur = fitlme(tbl,'rippleDuration~condition+age+(1|subject)+(1|subject:electrode)');
tbl = table([amp1;amp2;amp3;amp4;amp5],cond,repmat(age,[5,1]),repmat(s,[5,1]),repmat(e,[5,1]),'VariableNames',{'rippleAmplitude','condition','age','subject','electrode'});
lmeAmp = fitlme(tbl,'rippleAmplitude~condition+age+(1|subject)+(1|subject:electrode)');

% prepare outdir:
outdir =fullfile(parentfolder,'results/stats/lme/');
if ~exist(outdir), mkdir(outdir); end
writetable(dataset2table(lmeRate.anova),fullfile(outdir,'lmeRippleRate.csv'))
writetable(dataset2table(lmeDur.anova),fullfile(outdir,'lmeRippleDuration.csv'))
writetable(dataset2table(lmeAmp.anova),fullfile(outdir,'lmeRippleAmplitude.csv'))

% check the model fit:
% F = fitted(lmeRate);
% R = response(lmeRate);
% figure();
% plot(R,F,'rx')
% xlabel('Response')
% ylabel('Fitted')

% BOXPLOT:
H1=figure('Name',['Ripple rates boxplot'],'position',[0 200 200 200],'color','w'); hold on;
res1=lmeRate.anova;
res2=lmeRate.Coefficients;

dataRate = [n1, n2, n3, n4, n5];
boxplot(dataRate,{'R1','M1','R2','M2','TEST'},'plotstyle','traditional','colors',[COLOR.gray; COLOR.red; COLOR.gray; COLOR.red; COLOR.blue],'symbol','.','outliersize',2)
dataRate = {n1, n2, n3, n4, n5};
h=plotSpread(dataRate,'xNames',{'R1','M1','R2','M2','Mem test'},'categoryLabels',{'R1','M1','R2','M2','TEST'},...
    'distributionColors',[COLOR.gray; COLOR.red; COLOR.gray; COLOR.red; COLOR.blue]);
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',2)
title(sprintf('M=%.2f+%.2f\nF(%d,%d)=%.2f, P=%.2f',res2.Estimate(1),res2.SE(1),res1.DF1(2),res1.DF2(2),res1.FStat(2),res1.pValue(2)),'fontweight','normal')
set(gca, 'xlim',[0.5 5.5],'xtick',[1,2,3,4,5],'ylim',[0 0.8],'ytick',[0:0.2:0.8])
ylabel('Ripple rate (Hz)');
axis square
if saveflag, save_current_figure(H1,figdir,0); end


%% IRIs zoom in showing Doublets:
tmaxRest = restDuration; tmaxMovie = movieDuration;

DD1=[];
DD2=[];
DD3=[];
DD4=[];
DD5=[];
DD6=[];

maxt = 2; % in sec
binWidth = 0.05; % in sec
bins_edge=0:binWidth:maxt;
E = unique(allR1.elecid,'stable');

for ii = 1:numel(E)
    
    
    elecid = E{ii};
    subjid = split(elecid,'_'); subjid = subjid{1};
    
    % compute memory test duration:
    ind = find(contains({ALLEEG.setname}',subjid) & contains({ALLEEG.setname}',{'MT'}));
    EEG = ALLEEG(ind);
    strind = round(EEG.event(strcmpi({EEG.event.type},'str')).latency);
    finind = round(EEG.event(strcmpi({EEG.event.type},'fin')).latency);
    
    tailDuration = strind./EEG.srate;
    memoryTestDuration = round(max(ALLEEG(ind).times)-min(ALLEEG(ind).times))/1000-(2*tailDuration);
    
    
    REST1 = diff(allR1.peak(strcmpi(allR1.elecid,elecid) & (allR1.str >= 0 & allR1.str <= tmaxRest)));
    REST2 = diff(allR2.peak(strcmpi(allR2.elecid,elecid) & (allR2.str >= 0 & allR2.str <= tmaxRest)));
    MOVIE1 = diff(allM1.peak(strcmpi(allM1.elecid,elecid) & (allM1.str >= 0 & allM1.str <= tmaxMovie)));
    MOVIE2 = diff(allM2.peak(strcmpi(allM2.elecid,elecid) & (allM2.str >= 0 & allM2.str <= tmaxMovie)));
    MT = diff(allMT.peak(strcmpi(allMT.elecid,elecid) & (allMT.str >= 0 & allMT.str <= memoryTestDuration)));
    
    
    for k=1:6
        switch k
            case 1, Dtmp = REST1;
            case 2, Dtmp = REST2;
            case 3, Dtmp = MOVIE1;
            case 4, Dtmp = MOVIE2;
            case 5, Dtmp = MT;
            case 6, Dtmp = [REST1; REST2; MOVIE1; MOVIE2; MT];
        end        
        n=length(Dtmp);
        prob = histcounts(Dtmp,bins_edge,'Normalization','pdf');
        %prob = counts / (n * binWidth);
        switch k
            case 1, DD1(ii,:) = prob;
            case 2, DD2(ii,:) = prob;
            case 3, DD3(ii,:) = prob;
            case 4, DD4(ii,:) = prob;
            case 5, DD5(ii,:) = prob;
            case 6, DD6(ii,:) = prob;
        end
        
        
    end
end

%% PLOT: distribution of inter-ripple intervals by condition

H2a=figure('Name',['inter-ripple intervals - doublets'],'position',[0 200 1500 180],'color','w'); hold on;
xx = (bins_edge(1:end-1)+0.5*binWidth)*1000 ;
condLabels = {'REST1','REST2','MOVIE1','MOVIE2','MT','ALL'};
for k=1:6
    switch k
        case 1, cc=COLOR.blue;
        case 2, cc=COLOR.blue;
        case 3, cc=COLOR.red;
        case 4, cc=COLOR.red;
        case 5, cc=COLOR.green;
        case 6, cc=COLOR.black;
    end
    subplot(1,6,k); hold on;
    eval(sprintf('Dtmp = DD%d;',k));
    prob = mean(Dtmp);
    probSEM = std(Dtmp)./sqrt(size(Dtmp,1));
    superbar(xx,prob(1:end),'e',probSEM(1:end),'BarFaceColor',cc,'ErrorbarLineWidth',1)
    axis tight square
    xlim([0 maxt*1000]); ylim([0 1.5]);
    title(condLabels{k})
end
%title(sprintf('Mean = %.1f',double(mean(Dall))),'horizontalAlignment','center')
xlabel('Inter-ripple interval (ms)')
ylabel('Probability density')

%% Grand average - IRI distribution (Figure 2):
H2b = figure('Name',['inter-ripple intervals grand average'],'position',[0 200 220 160],'color','w'); hold on;
xx = (bins_edge(1:end-1)+0.5*binWidth)*1000 ;
prob = mean(DD6);
probSEM = std(DD6)./sqrt(size(DD6,1));
superbar(xx./1000,prob,'e',probSEM,'BarFaceColor',cc,'ErrorbarLineWidth',1);%,'BarRelativeGroupWidth',1,'BarEdgeColor','none')
axis tight
xlim([0 maxt]); ylim([0 1.5]);

xlabel('Inter-ripple interval (s)')
ylabel('Probability density')
if saveflag, save_current_figure(H2b,figdir,0); end


%% Re-organize the data for spectral analysis of the ripples time series (supp.):
close all
D1=[]; D2=[]; D3=[]; D4=[];
Crest = []; Cmov=[];
Prest = []; Pmov=[];
for ii = 1:numel(E)
    curData1 = rr1{ii};
    curData2 = rr2{ii};
    curData3 = rr3{ii};
    curData4 = rr4{ii};
   
    D1 = cat(1,D1, curData1);
    D2 = cat(1,D2, curData2);
    D3 = cat(1,D3, curData3);
    D4 = cat(1,D4, curData4);
end

%% ripple time-series spectral analysis:
ws = 15; % window size (sec)
[S1,freq] = spectopo(D1, 0, EEG.srate, 'plot','off','freqfac',8,'winsize',EEG.srate*ws,'overlap',EEG.srate*ws/2);
[S2,freq] = spectopo(D2, 0, EEG.srate, 'plot','off','freqfac',8,'winsize',EEG.srate*ws,'overlap',EEG.srate*ws/2);
[S3,freq] = spectopo(D3, 0, EEG.srate, 'plot','off','freqfac',8,'winsize',EEG.srate*ws,'overlap',EEG.srate*ws/2);
[S4,freq] = spectopo(D4, 0, EEG.srate, 'plot','off','freqfac',8,'winsize',EEG.srate*ws,'overlap',EEG.srate*ws/2);
Srest = (S1+S3)./2;
Smovie = (S2+S4)./2;

logfreq = log10(freq); logfreq(1) = logfreq(2)-(logfreq(3)-logfreq(2));

%% Ripple time series power spectra:
saveflag = 0;
HS = figure('name','ripple rate spectra','color','w','position',[0 0 250 250]);  hold on;
h1 = shadedErrorBar(logfreq, mean(S1),sem(S1).*0,{'color',COLOR.rest1,'linewidth',0.5},0.5); 
h2 = shadedErrorBar(logfreq, mean(S2),sem(S2).*0,{'color',COLOR.movie1,'linewidth',0.5},0.5);
h3 = shadedErrorBar(logfreq, mean(S3),sem(S3).*0,{'color',COLOR.rest2,'linestyle',':','linewidth',0.5},0.5); 
h4 = shadedErrorBar(logfreq, mean(S4),sem(S4).*0,{'color',COLOR.movie2,'linestyle',':','linewidth',0.5},0.5);

axis tight;
xlim([log10(0.01) log10(10)]);
set(gca,'xtick',log10([0.01 0.1 1 10]),'xticklabel',{'0.01','0.1','1','10'})
xlabel('log frequency (Hz)');
ylabel('Power (dB)');

test_ind = freq<10;
[sig, pval, onsets, offsets] = compare_conditions(gca,S2(:,test_ind)',S4(:,test_ind)',logfreq(test_ind),3,0)

if any(sig)
    R_time = logfreq;
    idx=sprintf('%d',sig);
    onsets = R_time(regexp(idx, '1{2,}', 'start'));    offsets = R_time(regexp(idx, '1{2,}', 'end'));
    tmp = get(gca,'Ylim');
    X=[onsets;offsets];
    Y=repmat(tmp(1)+range(tmp)*0.02,size(X));
    line(X,Y,'LineWidth',4,'Color',COLOR.sig)
    %text(mean(X),mean(Y)+0.05,sprintf('p<%.4f',min_cluster_pval),'fontsize',6,'horizontalalignment','center');
end

legend([h1.mainLine; h3.mainLine; h2.mainLine; h4.mainLine;],{'rest-1','rest-2','movie-1','movie-2'}); legend boxoff;
if saveflag, save_current_figure(HS,figdir,0); end


%% save figures;

for F=[H1 H3b] %H1 H2 H3]
    figure(F);
    set_font_size_and_type;
    set(findall(gcf,'-property','FontType'),'FontType','Arial')
    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-pdf','-painters','-nofontswap','-transparent')
    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-jpg','-painters','-r300','-nofontswap')
end
