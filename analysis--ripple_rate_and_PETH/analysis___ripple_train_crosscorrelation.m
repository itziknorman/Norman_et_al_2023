
% The script loads the ripple train in each experimental condition
% and compute crosscorrelation across repeated presentations of the
% animated events. Then the script fits a mixed effects model to examine
% the differentatial involvement of individual subfields and the
% recognition process.

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
outdir=fullfile(parentfolder,'results','crosscorrelation');
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
figdir = outdir;
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

% extract the memory test and tail duration from the first patient:
EEG = ALLEEG(5);
assert(contains(EEG.setname,'MT'))
strind = round(EEG.event(strcmpi({EEG.event.type},'str')).latency);
tailDuration = strind./EEG.srate;
memoryTestDuration = round(max(ALLEEG(1).times)-min(ALLEEG(1).times))/1000-(2*tailDuration);

%% Compute ripple rate across conditions:
Trest =  -tailDuration : (1/EEG.srate) : restDuration + tailDuration;
Tmovie = -tailDuration : (1/EEG.srate) : movieDuration + tailDuration;
Tmemorytest = -tailDuration : (1/EEG.srate) : movieDuration + tailDuration;

% init vars:
n1=[]; n2=[]; n3=[]; n4=[]; n5=[];
rr1={}; rr2={}; rr3={}; rr4={}; rr5={};
E = unique(allR1.elecid,'stable');
oldclips_time = {}; oldclips_label = {};
for ii = 1:numel(E)
    elecid = E{ii};
    subjid = split(elecid,'_'); subjid = subjid{1};
    
    % find the latency of the old clips within the memory test:
    ind = find(contains({ALLEEG.setname}',subjid) & contains({ALLEEG.setname}',{'MT'}));
    EEG = ALLEEG(ind);
    oldclips_time{ii} = round([EEG.event(contains({EEG.event.type},'ClipA')).latency])./EEG.srate;
    oldclips_label{ii} = string({EEG.event(contains({EEG.event.type},'ClipA')).type});
    
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
    rr1{ii} = ripplevec;
    
    swr_time = allM1.peak(strcmpi(allM1.elecid,elecid));
    ripplevec = zeros(size(Tmovie));
    for k=1:numel(swr_time), current_event = swr_time(k); [~,idx] = min(abs(Tmovie-current_event));  ripplevec(idx(1)) = 1; end
    %[t2 ,rrbin2{ii}] = bin_data(ripplevec,Tmovie,binsize,overlap,'rate');
    rr2{ii} = ripplevec;
    
    swr_time = allR2.peak(strcmpi(allR2.elecid,elecid));
    ripplevec = zeros(size(Trest));
    for k=1:numel(swr_time), current_event = swr_time(k); [~,idx] = min(abs(Trest-current_event));  ripplevec(idx(1)) = 1; end
    rr3{ii} = ripplevec;
    
    swr_time = allM2.peak(strcmpi(allM2.elecid,elecid));
    ripplevec = zeros(size(Tmovie));
    for k=1:numel(swr_time), current_event = swr_time(k); [~,idx] = min(abs(Tmovie-current_event));  ripplevec(idx(1)) = 1; end
    rr4{ii} = ripplevec;
    
    swr_time = allMT.peak(strcmpi(allMT.elecid,elecid));
    ripplevec = zeros(size(Tmemorytest));
    for k=1:numel(swr_time), current_event = swr_time(k); [~,idx] = min(abs(Tmemorytest-current_event));  ripplevec(idx(1)) = 1; end
    rr5{ii} = ripplevec;
    disp(ii);
end

%% Cross-correlation between 1st screening, 2nd screening, and memory test (supp.):
close all
D1=[]; D2=[]; D3=[]; D4=[]; D5=[];
D1conv=[]; D2conv=[]; D3conv=[]; D4conv=[];
CCrest=[]; CCmovie=[]; CCrestshuff=[]; CCmovieshuff=[];

Fs = EEG.srate;
nlags = 10;
sigma = 4; % reciprocal of the std (sigma=0.25)
winsize = 4; % in sec
gaussian_kernel = gausswin(winsize*EEG.srate,sigma);
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel);

for iElec = 1:numel(E)
    curData1 = rr1{iElec};
    curData2 = rr2{iElec};
    curData3 = rr3{iElec};
    curData4 = rr4{iElec};
    curData5 = rr5{iElec};
    D1 = cat(1,D1, curData1);
    D2 = cat(1,D2, curData2);
    D3 = cat(1,D3, curData3);
    D4 = cat(1,D4, curData4);
    D5 = cat(1,D5, curData5);
    
    % compute crosscorr accross the entire movie:
    curData1 = conv(curData1,gaussian_kernel,'same');
    curData2 = conv(curData2,gaussian_kernel,'same');
    curData3 = conv(curData3,gaussian_kernel,'same');
    curData4 = conv(curData4,gaussian_kernel,'same');
    D1conv = cat(1,D1conv, curData1);
    D2conv = cat(1,D2conv, curData2);
    D3conv = cat(1,D3conv, curData3);
    D4conv = cat(1,D4conv, curData4);
    
    [CCrest(iElec,:),lag] = xcorr(curData1,curData3,EEG.srate*nlags,'normalized');
    [CCmovie(iElec,:),lag] = xcorr(curData2,curData4,EEG.srate*nlags,'normalized');
    
    
    for iShuff = 1:50
        n = randi(length(curData2));
        shuff1 = circshift(curData1,n);
        shuff2 = circshift(curData2,n);
        % compute crosscorr after convolving the ripple train with a gaussian kernel:
        [CCrestshuff(iElec,:,iShuff),~] = xcorr(shuff1,curData3,EEG.srate*nlags,'normalized');
        [CCmovieshuff(iElec,:,iShuff),~] = xcorr(shuff2,curData4,EEG.srate*nlags,'normalized');
    end
    disp(iElec);
end

% read anatomical and recognition effectSize data:
T = readtable(fullfile(parentfolder,'results/stats/anatomicalModel_old-vs-new.csv'));
str = split(E,'_'); Eanat = strcat(str(:,2),'_',str(:,1));
assert(all(strcmpi(T.electrode,Eanat))==1); % make sure the anat data is aligned


%% plot:
CCrestnorm = bsxfun(@rdivide,CCrest,mean(mean(CCrestshuff,3),2));
CCmovienorm = bsxfun(@rdivide,CCmovie,mean(mean(CCmovieshuff,3),2));
CCrestshuffnorm = bsxfun(@rdivide,CCrestshuff,mean(mean(CCrestshuff,3),2));
CCmovieshuffnorm = bsxfun(@rdivide,CCmovieshuff,mean(mean(CCmovieshuff,3),2));
lagsec = lag./EEG.srate;
[str] = split(E,'_'); str=str(:,1);
ind = true(size(E));
[CCrest_grp,ss] = grpstats(CCrestnorm(ind,:),str(ind),{'mean','gname'});
[CCmovie_grp,ss] = grpstats(CCmovienorm(ind,:),str(ind),{'mean','gname'});
[CCrestshuff_grp,ss] = grpstats(mean(CCrestshuffnorm(ind,:,:),3),str(ind),{'mean','gname'});
[CCmovieshuff_grp,ss] = grpstats(mean(CCmovieshuffnorm(ind,:,:),3),str(ind),{'mean','gname'});

Hccfull = figure('color','w','name','ripples train cross correlations MOVIES','position',[0 0 150 150]); hold on;

shadedErrorBar(lagsec,mean(CCmovie_grp),sem(CCmovie_grp),{'color',COLOR.red},0.5);
SEshuff = sqrt(mean(var(CCmovieshuff_grp,[],2)));   % std across time
SEshuffrest = sqrt(mean(var(CCrestshuff_grp,[],2)));   % std across time
plot([min(lagsec) max(lagsec)],[1 1],'k-')
plot([min(lagsec) max(lagsec)],[1 1]+SEshuff,'k:')
plot([min(lagsec) max(lagsec)],[1 1]-SEshuff,'k:')
%shadedErrorBar(lagsec,mean(CCmovieshuff_grp),sem(CCmovieshuff_grp),{'color',COLOR.gray},0.5);

cond1 = CCmovie_grp';
cond2 = CCmovieshuff_grp';
xlabel('Lag (s)')
axis tight
tmpy = get(gca,'ylim');
set(gca,'ylim',[tmpy(1),ceil(tmpy(2)*10)/10],'xtick',[-10:5:10],'ytick',[1:0.1:2]);
%[sig, pval, onsets, offsets] = compare_conditions(gca,cond1,cond2,lagsec,1,1);
ylabel('Corr. coeff (NU)')
subtitle({'Crosscorrelogram','group level'})

set_font_size_and_type;
if saveflag, save_current_figure(gcf,figdir,0); end

%% epoching around events of interest (those included in  the memory test):

% Timing of original events in movie: (in sec)
eventsInMovie = readtable(fullfile(parentfolder,'events.xlsx'));
[origEventTime,sortind] = sort(eventsInMovie.timeInSec);
eventsInMovie = eventsInMovie(sortind,:);

[eventsInMovie1, epochlim] = epoch(D2,origEventTime+tailDuration,[-1 5],'srate',Fs);
[eventsInMovie2, epochlim] = epoch(D4,origEventTime+tailDuration,[-1 5],'srate',Fs);
epoch_time = epochlim(1):1/Fs:epochlim(2);

% get the corresponding clips from the memory test: (in the same order as
% in the full length screening)
eventsInMemoryTest = [];
for iElec = 1:length(E)
    clipLabel = oldclips_label{iElec};
    clipTime = oldclips_time{iElec};
    [ia,ib] = ismember(eventsInMovie.name,clipLabel);
    ib(ib==0) = []; % remove missing items
    clipLabel = clipLabel(ib)';
    clipTime = clipTime(ib)';
    [tmp, ~] = epoch(D5(iElec,:),clipTime+tailDuration,[-1 5],'srate',Fs);
    
    eventsInMemoryTest(iElec,:,:) = zeros(length(epoch_time),length(eventsInMovie.name));
    eventsInMemoryTest(iElec,:,ia) =  squeeze(tmp);
end

%% cross correlation limited to the events presented in the clips:

xcorr_clip = [];
xcorr_clip_shuff = [];
xcorr_mt1 = [];
xcorr_mt2 = [];
xcorr_mt_shuff = [];

nlags = 10; % in sec
for iElec = 1:size(eventsInMovie1,1)
    
    tmp1 = squeeze(eventsInMovie1(iElec,:,:));
    tmp2 = squeeze(eventsInMovie2(iElec,:,:));
    tmp3 = squeeze(eventsInMemoryTest(iElec,:,:));
    
    % compute cross corr after convolving the ripple train with a gaussian kernel:
    [xcorr_clip(iElec,:),lag_clip] = xcorr(conv(tmp1(:),gaussian_kernel,'same'),...
        conv(tmp2(:),gaussian_kernel,'same'), EEG.srate*nlags,'normalized');
    
    % compute cross corr after convolving the ripple train with a gaussian kernel:
    [xcorr_mt1(iElec,:),lag_clip] = xcorr(conv(tmp1(:),gaussian_kernel,'same'),...
        conv(tmp3(:),gaussian_kernel,'same'), EEG.srate*nlags,'normalized');
    
    % compute cross corr after convolving the ripple train with a gaussian kernel:
    [xcorr_mt2(iElec,:),lag_clip] = xcorr(conv(tmp2(:),gaussian_kernel,'same'),...
        conv(tmp3(:),gaussian_kernel,'same'), EEG.srate*nlags,'normalized');
    
    for iShuff = 1:50
        ix = randperm(size(eventsInMovie1,3));
        tmp1 = squeeze(eventsInMovie1(iElec,:,ix));
        
        % compute crosscorr after convolving the ripple train with a gaussian kernel:
        [xcorr_clip_shuff(iElec,:,iShuff),~] = xcorr(conv(tmp1(:),gaussian_kernel,'same'),...
            conv(tmp2(:),gaussian_kernel,'same'),EEG.srate*nlags,'normalized');
        
        [xcorr_mt_shuff(iElec,:,iShuff),~] = xcorr(conv(tmp1(:),gaussian_kernel,'same'),...
            conv(tmp3(:),gaussian_kernel,'same'),EEG.srate*nlags,'normalized');
        
    end
    disp(iElec)
end


%% Ripple rate correlation between 1st anc 2nd screening:
x_clip = squeeze(sum(eventsInMovie1,2));
y_clip = squeeze(sum(eventsInMovie2,2));
z_clip = squeeze(sum(eventsInMemoryTest,2));

n_clips1 = sum(sum(eventsInMovie1,2),3)./(range(epochlim)*size(eventsInMovie1,3));
n_clips2 = sum(sum(eventsInMovie2,2),3)./(range(epochlim)*size(eventsInMovie2,3));
n_clipsMT = sum(sum(eventsInMemoryTest,2),3)./(range(epochlim)*size(eventsInMemoryTest,3));

[rho1,pval1] = corr(n2,n4,'type','pearson');
[rho2,pval2] = corr(n1,n3,'type','pearson');

figure('color','w','name','ripple rate correlations between screenings','position',[0 0 150 150]); hold on;
scatter(n2,n4,8,'ko')
axis square tight
xtmp = get(gca,'xlim'); ytmp = get(gca,'ylim');
xl = [min([xtmp,ytmp]), max([xtmp,ytmp])];
xlabel({'Screening #1 (ripple/s)'});
ylabel({'Screening #2 (ripple/s)'});
% subtitle(sprintf('r=%.2f, p=10^{%.0f}',rho1,log10(pval1)))
subtitle(sprintf('R^2=%.2f',rho1.^2),'fontsize',6)
xlim(xl); ylim(xl);
plot(xl,xl,':k')
xticks([0.1:0.2:1]);
yticks([0.1:0.2:1]);
set_font_size_and_type;
if saveflag, save_current_figure(gcf,figdir,0); end

figure('color','w','name','ripple rate correlations between rest','position',[0 0 150 150]); hold on;
scatter(n1,n3,8,'ko')
axis square tight
xtmp = get(gca,'xlim'); ytmp = get(gca,'ylim');
%xl = [min([xtmp,ytmp]), max([xtmp,ytmp])];
xlabel({'Rest #1 (ripple/s)'});
ylabel({'Rest #2 (ripple/s)'});
subtitle(sprintf('R^2=%.2f',rho2.^2),'fontsize',6)
xlim(xl); ylim(xl);
plot(xl,xl,':k')
xticks([0.1:0.2:1]);
yticks([0.1:0.2:1]);
set_font_size_and_type;
if saveflag, save_current_figure(gcf,figdir,0); end

%% line plot:
[str] = split(E,'_'); str=str(:,1);
ind = true(size(E));
[rate_diff_movie_grp,ss] = grpstats((n2-n4).^2,str,{'mean','gname'});
[rate_diff_rest_grp,ss] = grpstats((n1-n3).^2,str,{'mean','gname'});

figure('color','w','name','ripple rate variance movie vs rest','position',[0 0 110 135]);
line([ones(size(ss)),2*ones(size(ss))]',[rate_diff_movie_grp,rate_diff_rest_grp]','color',COLOR.gray);
hold on;
scatter(ones(size(ss)),rate_diff_movie_grp,15,'markerfacecolor',COLOR.black,'markeredgecolor','none'); hold on;
scatter(2.*ones(size(ss)),rate_diff_rest_grp,15,'markerfacecolor',COLOR.black,'markeredgecolor','none'); hold on;
axis tight
ytmp = get(gca,'ylim');
line([1,2],[ytmp(2),ytmp(2)].*1.2,'color','k')
xl = [0, ytmp(2).*1.5];
ylabel('ripple rate \Delta');
p_val_diff = signrank(rate_diff_movie_grp,rate_diff_rest_grp);
subtitle(sprintf('P=%.4f',p_val_diff),'fontsize',6)
xlim([0.5 2.5]); ylim(xl);
set(gca,'xtick',[1 2],'XtickLabel',[string('Movie'),string('Rest')]);
set_font_size_and_type;
if saveflag, save_current_figure(gcf,figdir,0); end


%% Exmple ripples train - full-length:

% select the best representative electrode:
% xcorr_clip_norm = bsxfun(@rdivide,xcorr_clip,mean(mean(xcorr_clip_shuff,3),2));
% xcorr_clip_norm(n2<0.2,:)=nan;
% [mx,mxind]  = max(xcorr_clip_norm,[],2);    % only events participating in the memory test
% [~,best_elec] = max(mx);

xcorr_movie_norm = CCmovienorm;
xcorr_movie_norm(n2<0.2,:)=nan;            % exclude sites with a very low ripple rate
[mx,mxind]  = max(xcorr_movie_norm,[],2);  % find the site showing max crosscorr
[~,best_elec] = max(mx);

Htrain = figure('color','w','name','ripple train full','position',[0 0 350 100]); hold on;
tvec = Tmovie;
subsampleind = 1:2:length(tvec);
M = [       zeros(size(subsampleind));
    D2conv(best_elec,subsampleind);
    zeros(size(subsampleind));
    D4conv(best_elec,subsampleind);
    zeros(size(subsampleind));
    zeros(size(subsampleind)); ];

imagesc(tvec(subsampleind),1:6,M);
cl = prctile(abs(M(:)),95);
colormap(cbrewer('seq','Greys',32)); caxis([0 cl]);
xlabel('Time (s)')
ylim([1 5.5]);xlim([min(tvec),max(tvec)]);
set(gca,'ytick',[2,4],'YTickLabel',[string('screening #1'), string('screening #2')])
ch_label = split(E{best_elec},'_');
ch_label = sprintf('%s (%s)',ch_label{1},ch_label{2});
subtitle(sprintf('Ripple train example: %s',string(ch_label)));

if saveflag, save_current_figure(gcf,figdir,0); end

%% example crosscorr - entire cartoon:

Hccfull = figure('color','w','name','ripples train cross correlations MOVIES example electrode','position',[0 0 150 150]); hold on;

plot(lagsec,CCmovienorm(best_elec,:),'color',COLOR.black);
SEshuff = sqrt(mean(var(CCmovieshuff(best_elec,:),[],2)));   % std across time
SEshuffrest = sqrt(mean(var(CCrestshuff(best_elec,:),[],2)));   % std across time
plot([min(lagsec) max(lagsec)],[1 1],'k-')
plot([min(lagsec) max(lagsec)],[1 1]+SEshuff,'k:')
plot([min(lagsec) max(lagsec)],[1 1]-SEshuff,'k:')
xlabel('Lag (s)')
axis tight
tmpy = get(gca,'ylim');
set(gca,'ylim',[floor(tmpy(1)*10)/10,ceil(tmpy(2)*10)/10],'ytick',[1:0.5:2],'xtick',[-10:5:10]);
ylabel('Corr. coeff (NU)')
subtitle({'Crosscorrelogram',''})

set_font_size_and_type;
if saveflag, save_current_figure(gcf,figdir,0); end


%% exmple ripples train - only clips:

[mx,mxind]  = max(xcorr_clip_norm,[],2);
xcorr_clip_norm(n2<0.2,:)=nan;
[~,best_elec] = max(mx);
%[~,best_elec] = min(abs(lagsec(mxind)));

Htrain = figure('color','w','name','ripple train clips','position',[0 0 500 80]); hold on;
tmp1 = squeeze(eventsInMovie1(best_elec,:,:));
tmp2 = squeeze(eventsInMovie2(best_elec,:,:));
tmp3 = squeeze(eventsInMemoryTest(best_elec,:,:));

M = [conv(tmp1(:),gaussian_kernel,'same')'; -conv(tmp2(:),gaussian_kernel,'same')'; conv(tmp3(:),gaussian_kernel,'same')';];
imagesc([1/Fs:1/Fs:size(M,2)]./EEG.srate,2:3,M(1:2,:)); hold on;
cl = prctile(abs(M(:)),95);
cm = cbrewer('div','PiYG',32); cm(15:17,:) = ones(3,3);
colormap(cm); axis tight; caxis([-cl cl]); freezeColors;
imagesc([1/Fs:1/Fs:size(M,2)]./EEG.srate,1,M(3,:)); hold on;
cm = cbrewer('seq','Blues',32); cm(1:3,:) = ones(3,3);
colormap(cm); axis tight; caxis([0 cl]); freezeColors;
xlabel('Time (s)')
set(gca,'ytick',[1:3],'YTickLabel',{'Test',sprintf('first'),sprintf('sec')})
ch_label = split(E{best_elec},'_');
ch_label = join(ch_label',' ');
subtitle(strcat(ch_label,' (20 events)'))
for i = 1:size(eventsInMovie1,3)-1
    tmp = i * size(eventsInMovie1,2)./EEG.srate;
    x = [tmp tmp]';
    y = get(gca,'ylim')';
    hold on; plot(x,y,'k-','linewidth',1);
end
if saveflag, save_current_figure(gcf,figdir,0); end


%% cross correlation:
[str] = split(E,'_'); str=str(:,1);
% normalize by the null distribution:
xcorr_clip_norm = bsxfun(@rdivide,xcorr_clip,mean(mean(xcorr_clip_shuff,3),2));
xcorr_mt1_norm = bsxfun(@rdivide,xcorr_mt1,mean(mean(xcorr_mt_shuff,3),2));
xcorr_mt2_norm = bsxfun(@rdivide,xcorr_mt2,mean(mean(xcorr_mt_shuff,3),2));
xcorr_clip_shuff_norm = bsxfun(@rdivide,xcorr_clip_shuff,mean(mean(xcorr_clip_shuff,3),2));
xcorr_mt_shuff_norm = bsxfun(@rdivide,xcorr_mt_shuff,mean(mean(xcorr_mt_shuff,3),2));
% average within patient:
ind = T.effectSize>0.3; %  strcmpi(T.Hemisphere,'r');
[xcorr_clip_grp,ss] = grpstats(xcorr_clip_norm,str,{'mean','gname'});
[xcorr_mt1_grp,ss] = grpstats(xcorr_mt1_norm,str,{'mean','gname'});
[xcorr_mt2_grp,ss] = grpstats(xcorr_mt2_norm,str,{'mean','gname'});
[xcorr_clip_shuff_grp,ss] = grpstats(mean(xcorr_clip_shuff_norm,3),str,{'mean','gname'});
[xcorr_mt_shuff_grp,ss] = grpstats(mean(xcorr_mt_shuff_norm,3),str,{'mean','gname'});

% ind1 = T.effectSize > -1;
% ind2 = T.effectSize < 0;
Hcc_clips = figure('color','w','name','ripples train cross correlations between clip screenings','position',[0 0 200 200]); hold on;
shadedErrorBar(lagsec,mean(xcorr_clip_norm(ind,:)),sem(xcorr_clip_norm(ind,:)),{'color',COLOR.red},0.5);

% null distribution:
SEshuffclip = sqrt(mean(var(xcorr_clip_shuff_grp,[],2)));   % std across time
plot([min(lagsec) max(lagsec)],[1 1],'k-')
plot([min(lagsec) max(lagsec)],[1 1]+SEshuffclip,'k:')
plot([min(lagsec) max(lagsec)],[1 1]-SEshuffclip,'k:')

cond1 = xcorr_clip_grp';
cond2 = xcorr_clip_shuff_grp';
xlabel('Lag (s)')
axis tight
[sig, pval, onsets, offsets] = compare_conditions(gca,cond1,cond2,lagsec,1,1);
tmpy = get(gca,'ylim');
tmpy = [floor(tmpy(1)*10)/10,ceil(tmpy(2)*10)/10];
set(gca,'ylim',tmpy);
xlabel('Lag (s)')
ylabel('Crosscorr. coeff (NU)')
subtitle({'Ripples crosscorrelogram','1st-vs-2nd screenings'})
if saveflag, save_current_figure(gcf,figdir,0); end

Hcc_clips_mt = figure('color','w','name','ripples train cross correlations with MT','position',[0 0 200 200]); hold on;
col = cbrewer('div','PiYG',5);
shadedErrorBar(lagsec,mean(xcorr_mt1_grp),sem(xcorr_mt1_grp),{'color',col(1,:)},0.5);
shadedErrorBar(lagsec,mean(xcorr_mt2_grp),sem(xcorr_mt2_grp),{'color',col(end,:)},0.5);
SEshuffmt = sqrt(mean(var(xcorr_mt_shuff_grp,[],2)));   % std across time

plot([min(lagsec) max(lagsec)],[1 1],'k-')
plot([min(lagsec) max(lagsec)],[1 1]+SEshuffmt,'k:')
plot([min(lagsec) max(lagsec)],[1 1]-SEshuffmt,'k:')

axis tight
set(gca,'ylim',tmpy);
xlabel('Time (s)')
ylabel('Crosscorr. coeff (NU)')
subtitle({'Ripples crosscorrelogram','screenings vs recognition-test'})

% stats: (NS)
% cond1 = xcorr_mt1_grp';
% cond2 = xcorr_mt_shuff_grp';
% [sig, pval, onsets, offsets] = compare_conditions(gca,cond1,cond2,lagsec,1,1);
% cond1 = xcorr_mt2_grp';
% cond2 = xcorr_mt_shuff_grp';
% [sig, pval, onsets, offsets] = compare_conditions(gca,cond1,cond2,lagsec,1,1);

if saveflag, save_current_figure(gcf,figdir,0); end


%% save figures;

for F=[H1 H3b] %H1 H2 H3]
    figure(F);
    set_font_size_and_type;
    set(findall(gcf,'-property','FontType'),'FontType','Arial')
    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-pdf','-painters','-nofontswap','-transparent')
    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-jpg','-painters','-r300','-nofontswap')
end

