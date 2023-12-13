% The script loads the SWR time stamps across all hippocampal sites
% (saved as rasters) and computes peri-stimulus time histograms and overall
% ripple rates, individually for each condition.
%
% The SWR raster data is constructed beforehand using the script:
% "compute_swr_raster_per_electrode.m"

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

% Set Manually:
ref_flag=2; % 1 = Common Ref; 2 = Bipolar montage;
time_locking_event = 'stimonset'; % time locking event (stimonset/rt)
set_figure_colors;
subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};

% load demographic:
demographics = readtable(fullfile(parentfolder,'subj_info_anon.xlsx'));
set(0,'DefaultAxesFontName', 'Arial')
% Load rasters and initalize general varaiables:
load_swr_rasters_from_all_hippocampal_channels;

%% Load & plot Hippocampal electrodes:
L1 = unique(DATA.Rall.channelid,'stable');
[elecinfo,Hmap] = getHippocampalElectrodeCoord(L1, subjects,0,1);

selectedElecInd = true(size(DATA.Rall.channelid)); % optional variable to limit the analysis to a subset of electrodes
fprintf('number of verified hippocampal electrodes: %d \n \n',length(unique(DATA.Rall.channelid)))

%% Single trial RIPPLE RATE:
close all;

% compute overal ripple PSTH to find the latency of overall ripple rate peak (during video clip presentation):
ind = contains(DATA.Rall.stimulus_type,'Clip') & selectedElecInd;

% set PETH parameters:
binsize = 50;
smt = 5;
[r,binscenters]=nanPSTH(logical(DATA.Rall.raster(ind,:))',binsize,smt,DATA.Fs, []);
binscenters = (binscenters/1000)+DATA.trialtime(1);

% compute bootstrap data:
opt = statset('UseParallel',true);
bs = bootstrp(100,@(x)nanPSTH(x',binsize,smt,Fs,[]),logical(DATA.Rall.raster(ind,:)),'Options',opt);
ix = find(binscenters<4 & binscenters>0);
% peak latency:
[~,pk] = max(r);              % find the actual peak
[~,bspk] = max(bs(:,ix),[],2);    % compute  SE of peak latency
% half height:
halfHeight = (min(bs(:,ix),[],2) + max(bs(:,ix),[],2)) / 2;
bshh = [];
for i = 1:size(bs,1), bshh(i,1) = find(bs(i,ix)>halfHeight(i),1,'first'); end

% plot the response window:
H0 = figure('color','w','name','peak SWR response','position',[0 0 200 200]); hold on;
plot(binscenters,r,'linewidth',2','color',COLOR.red);
set(gca,'ytick',0.3:0.1:1)
tmpy=get(gca,'ylim'); tmpx=get(gca,'xlim');
line([0 4],[tmpy(1) tmpy(1)]+0.01*range(tmpy),'color',COLOR.black,'linewidth',5)
scatter(binscenters(pk),r(pk),50,'ok','linewidth',2);
switch time_locking_event, case 'stimonset', xlim([-2 6]); case 'rt', xlim([-6 2]); end
xlabel('Time from clip onset (s)'); ylabel('Ripple rate (events/s)')
subtitle({sprintf('peak: %.3f (%.3f) s',mean(binscenters(ix(bspk))),std(binscenters(ix(bspk)))),...
          sprintf('halfHeight: %.3f (%.3f) s',mean(binscenters(ix(bshh))),std(binscenters(ix(bshh)))) });

% set the response time window to +/- 1 sec around the peak:
timeWindow = [binscenters(pk)-1 binscenters(pk)+1];
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


%% Group level analysis: bar plot - ripple rate across conditions

%selectedElecInd = CA1elecInd; % use this to restrict the analysis to the CA1 subfield

rr = struct;
ss = struct;
avg = struct;
se = struct;
% Compute grp-level ripple rates - video clips:
indClipA = startsWith(DATA.Rall.stimulus,'ClipA') & ~contains(DATA.Rall.response,'new') &  selectedElecInd;
indClipB = startsWith(DATA.Rall.stimulus,'ClipB') & ~contains(DATA.Rall.response,'old') &  selectedElecInd;
[rr.clipA,ss.clipA] = grpstats(RRnormTW(indClipA),DATA.Rall.channelid(indClipA),{'nanmean','gname'});
[rr.clipB,ss.clipB] = grpstats(RRnormTW(indClipB),DATA.Rall.channelid(indClipB),{'nanmean','gname'});

indM1clipA = startsWith(DATA.Rall.stimulus,'M1ClipA') &  selectedElecInd;
indM2clipA = startsWith(DATA.Rall.stimulus,'M2ClipA') &  selectedElecInd;
[rr.M1clipA,ss.M1clipA] = grpstats(RRnorm4s(indM1clipA),DATA.Rall.channelid(indM1clipA),{'nanmean','gname'});
[rr.M2clipA,ss.M2clipA] = grpstats(RRnorm4s(indM2clipA),DATA.Rall.channelid(indM2clipA),{'nanmean','gname'});

% Compute grp-level ripple rates - full movie/rest:
indM1 = contains(DATA.Rall.stimulus_type,'M1') &  selectedElecInd;
indM2 = contains(DATA.Rall.stimulus_type,'M2') & selectedElecInd;
indR1 = contains(DATA.Rall.stimulus_type,'R1') & selectedElecInd;
indR2 = contains(DATA.Rall.stimulus_type,'R2') &  selectedElecInd;

[rr.M1,ss.M1] = grpstats(RRnorm4s(indM1),DATA.Rall.channelid(indM1),{'nanmean','gname'});
[rr.M2,ss.M2] = grpstats(RRnorm4s(indM2),DATA.Rall.channelid(indM2),{'nanmean','gname'});
[rr.R1,ss.R1] = grpstats(RRnorm4s(indR1),DATA.Rall.channelid(indR1),{'nanmean','gname'});
[rr.R2,ss.R2] = grpstats(RRnorm4s(indR2),DATA.Rall.channelid(indR2),{'nanmean','gname'});

% average within patient:
fn = fieldnames(ss);
for ii = 1:length(fn)
    [str] = split(ss.(fn{ii}),'_'); str=str(:,2);
    [rr.(fn{ii}),ss.(fn{ii})] = grpstats(rr.(fn{ii}),str,{'nanmean','gname'});
    avg.(fn{ii}) = nanmean(rr.(fn{ii}));
    se.(fn{ii}) =  sem(rr.(fn{ii}));
end


% plot figures:

% memory test: old vs new
H1a = figure('Name',sprintf('ripple rate across clips norm flag %d',norm_flag),'position',[0 0 100 180],'color','w'); hold on;
superbar([1 2],[avg.clipB avg.clipA],'E',[se.clipB se.clipA],'BarFaceColor',[COLOR.red; COLOR.blue],'ErrorbarLineWidth',1,'baredgecolor','none')
set(gca, 'xlim',[0.25 2.75],'ylim',[0.2 0.61],'ytick',[0:0.1:1],'xtick',[1 2],'xticklabel',{''})
ylabel('Ripple rate (events/s)')
[p1]=signrank(rr.clipA,rr.clipB);
tmpy=get(gca,'ylim');
line([1;2],[tmpy(2);tmpy(2)].*0.92,'color','k')
text(1.5,tmpy(2),sprintf('P<%.2f\n*',ceil(p1*100)./100),'fontsize',6,'HorizontalAlignment','center')
text([1,2],repmat(tmpy(1)-0.02*range(tmpy),[1,2]),{'NEW', 'OLD'},...
    'fontsize',6,'HorizontalAlignment','right','rotation',45);
fprintf('\n Memory Test p-value < %.4f (n=%d) \n',ceil(p1*100)./100,length(ss.clipA));%,'fontweight','normal')
set(gcf,'units','pixels')
pos = get(gca,'position');
set(gca,'position',pos)
enlarge_figure_and_move_axes_to_center(gcf,gca,1.25);

% memory test vs passive viewing:
H1b = figure('Name',sprintf('ripple rate across clips in movie norm flag %d',norm_flag),'position',[0 0 135 180],'color','w'); hold on;

% stats:
p = [];
p(1) = signrank(rr.M1clipA,rr.M2clipA); % episodic-rest
p(2) = signrank(rr.M1clipA,rr.clipA); % eposodic-math
p(3) = signrank(rr.M2clipA,rr.clipA); % rest-math
[~, crit_p, ~, adj_p]=fdr_bh(p,0.05,'pdep','yes');

P = nan(length(p));
P = squareform(adj_p); P(logical(eye(size(P))))=nan;
P(P>0.05)=nan;

superbar([1 2 3],[avg.M1clipA avg.M2clipA avg.clipA],'E',[se.M1clipA se.M2clipA se.clipA],'BarFaceColor',[COLOR.turquoise; COLOR.turquoise; COLOR.blue],'ErrorbarLineWidth',1,'baredgecolor','none',...
    'P',P,'PStarIcon','*','PStarFontSize',8,'PLineColor','k','PLineWidth',0.1,'PStarIcon','*','PStarShowNS',0,'PStarColor',[0 0 0]);
set(gca, 'xlim',[0.25 3.75],'ylim',[0.2 0.61],'ytick',[0:0.1:1],'xtick',[1 2],'xticklabel',{''})
ylabel('Ripple rate (events/s)')
tmpy=get(gca,'ylim');
text([1,2,3],repmat(tmpy(1)-0.02*range(tmpy),[1,3]),{'1st watch', '2nd watch' 'judged OLD'},...
    'fontsize',6,'HorizontalAlignment','right','rotation',45);
set(gcf,'units','pixels')
pos = get(gca,'position');
set(gca,'position',pos)
enlarge_figure_and_move_axes_to_center(gcf,gca,1.25);

% movie:
H1c = figure('Name',sprintf('ripple rate across movies norm flag %d',norm_flag),'position',[200 0 100 180],'color','w'); hold on;
superbar([1 2],[avg.M1 avg.M2],'E',[se.M1, se.M2],'BarFaceColor',[COLOR.graydark; COLOR.graydark],'ErrorbarLineWidth',1,'baredgecolor','none')
set(gca, 'xlim',[0.25 2.75],'ylim',[0.2 0.6],'ytick',[0:0.1:1],'xtick',[1 2],'xticklabel',{''})
ylabel('Ripple rate (events/s)')
[p2]=signrank(rr.M1,rr.M2);
tmpy=get(gca,'ylim');

text(1.5,tmpy(2),sprintf('P=%.3f\n(n.s)',p2),'fontsize',6,'HorizontalAlignment','center')
text([1,2],repmat(tmpy(1)-0.02*range(tmpy),[1,2]),{'1st full', '2nd full'},...
    'fontsize',6,'HorizontalAlignment','right','rotation',45);
fprintf('\n Movie p-value < %.4f (n=%d) \n',ceil(p2*1000)./1000,length(rr.M1));%,'fontweight','normal')
set(gcf,'units','pixels')
set(gca,'position',pos)
enlarge_figure_and_move_axes_to_center(gcf,gca,1.25);

% rest:
H1d = figure('Name',sprintf('ripple rate across rest norm flag %d',norm_flag),'position',[400 0 100 180],'color','w'); hold on;
superbar([1 2],[avg.R1 avg.R2],'E',[se.R1, se.R2],'BarFaceColor',[COLOR.graydark; COLOR.graydark],'ErrorbarLineWidth',1,'baredgecolor','none')
set(gca, 'xlim',[0.25 2.75],'ylim',[0.2 0.6],'ytick',[0:0.1:1],'xtick',[1 2],'xticklabel',{''})
ylabel('Ripple rate (events/s)')
[p3]=signrank(rr.R1,rr.R2);
tmpy=get(gca,'ylim');
%line([1;2],[tmpy(2);tmpy(2)].*0.92,'color','k')
text(1.5,tmpy(2),sprintf('P=%.2f\n(n.s)',p3),'fontsize',6,'HorizontalAlignment','center')
text([1,2],repmat(tmpy(1)-0.02*range(tmpy),[1,2]),{'REST 1', 'REST 2'},...
    'fontsize',6,'HorizontalAlignment','right','rotation',45);
fprintf('\n REST p-value < %.4f (n=%d) \n',ceil(p3*1000)./1000,length(rr.R1));%,'fontweight','normal')
set(gcf,'units','pixels')
set(gca,'position',pos)
enlarge_figure_and_move_axes_to_center(gcf,gca,1.25);

% memory test scatter:
Hsc = figure('color','w','name','scatter plot ripple rate','position',[600 0 220 220]);
scatter(rr.clipB,rr.clipA,15,'ok','markerfacecolor',COLOR.gray,'Linewidth',1); axis square;
% hold on;
% scatter(rr.clipB(strcmpi(demographics.Coverage,'right')),rr.clipA(strcmpi(demographics.Coverage,'right')),15,'.r','markerfacecolor',COLOR.gray,'Linewidth',1); axis square;
axis([0.1 0.9 0.1 0.9]);
hold on; plot(get(gca,'xlim'),get(gca,'ylim'),'k:','Linewidth',1);
set(gca,'xtick',[0.2:0.2:0.8],'ytick',[0.2:0.2:0.8]);
tmpy=get(gca,'ylim');
text(0.5,tmpy(2)+0.1*range(tmpy),'Memory test','fontweight','normal','HorizontalAlignment','center');
%text(RR1(RR1>RR2),RR2(RR1>RR2),ElecAnatLoc(RR1>RR2))'ErrorbarRelativeWidth',0.25
xlabel({'ripple rate (events/s)','NEW'}); ylabel({'OLD';'ripple rate (events/s)'})
set_font_size_and_type;
enlarge_figure_and_move_axes_to_center(gcf,gca,1.25);

%% Behavioral performance:
singleElecInd = false(size(DATA.Rall.channelid));
for i = 1:length(subjects)
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjects{i},ref_flag,0);
    singleElecInd = singleElecInd | contains(DATA.Rall.channelid,[hippocampus,'_',subjects{i}]);
end

tmpind = contains(DATA.Rall.stimulus_type,'Clip') & ~contains(DATA.Rall.response,'n/a') & singleElecInd;
[grpAcc,gname] = grpstats(DATA.Rall.correct(tmpind),DATA.Rall.subjid(tmpind),{'nanmean','gname'});
[grpHit,gname1] = grpstats(DATA.Rall.correct(tmpind & contains(DATA.Rall.response,'old')),...
    DATA.Rall.subjid(tmpind & contains(DATA.Rall.response,'old')),{'nansum','gname'});
[grpFA,gname2] = grpstats(~DATA.Rall.correct(tmpind & contains(DATA.Rall.response,'old')),...
    DATA.Rall.subjid(tmpind & contains(DATA.Rall.response,'old')),{'nansum','gname'});
[grpCr,~] = grpstats(DATA.Rall.correct(tmpind & contains(DATA.Rall.response,'new')),...
    DATA.Rall.subjid(tmpind & contains(DATA.Rall.response,'new')),{'nansum','gname'});

tmpind = contains(DATA.Rall.stimulus_type,'Clip') & ~isnan(DATA.Rall.RT) & singleElecInd & ismember(DATA.Rall.subjid,[gname1; gname2]);
indA = contains(DATA.Rall.stimulus_type,'ClipA');
indB = contains(DATA.Rall.stimulus_type,'ClipB');
[nOLD] = grpstats(indA(tmpind),DATA.Rall.subjid(tmpind),{'nansum'});
[nNOVEL] = grpstats(indB(tmpind),DATA.Rall.subjid(tmpind),{'nansum'});
grpHitCr = (grpHit+grpCr)./(nOLD+nNOVEL);
grpHit = grpHit./nOLD;
grpFA = grpFA./nNOVEL;
grpCr = grpCr./nNOVEL;
% calculate d-prime:
dp = arrayfun(@(x1,x2,x3,x4)dprime(x1,x2,x3,x4),grpHit,grpFA,nOLD,nNOVEL);

[m1,v1,elecName1,n1] = grpstats(RRnormTW(indClipA),DATA.Rall.channelid(indClipA),{'mean','var','gname','numel'}); % NEW
[m2,v2,elecName2,n2] = grpstats(RRnormTW(indClipB),DATA.Rall.channelid(indClipB),{'mean','var','gname','numel'}); % OLD

% Compute hedge's G effect size: (OLD-vs-NEW)
df=n2+n1-2;
sP=((n2-1).*v2 + (n1-1).*v1)./(n2+n1-2); % pooled (within-groups) variance
es=(m2-m1)./sqrt(sP);

% average effect size within patient:
[str] = split(elecName1,'_'); str=str(:,2);
[ES,esSubj] = grpstats(es,str,{'mean','gname'});

fprintf('\nOverall accuracy in the task: %.3f (+/-%.3f)',mean(grpAcc),std(grpAcc)./sqrt(length(grpAcc)));
fprintf('\nOverall accuracy in the task: %.3f (+/-%.3f)',mean(1-grpAcc),std(1-grpAcc)./sqrt(length(grpAcc)));

% performance bar plot:
meas = grpAcc*100;  % select dp (d-prime) or grpAcc*100 (% correct)
[~,ix] = ismember(gname,demographics.SubjectID); ix(ix==0) = []; assert(all(strcmpi(demographics.SubjectID(ix),gname)));
age = demographics.Age(ix); % patient age
[r,p] = corr(meas,age,'type','pearson');
Hbehavior = figure('Name',sprintf('Performance'),'position',[0 0 350 180],'color','w'); hold on;
superbar([1:length(meas)],meas,'BarFaceColor',COLOR.gray,'ErrorbarLineWidth',1)
%set(gca, 'xlim',[0.25 length(meas)+0.75],'ylim',[0 4],'ytick',[0:1:5],'xtick',[1:length(dp)],'xticklabel',gname)
set(gca, 'xlim',[0.25 length(grpAcc)+0.75],'ylim',[0 100],'ytick',[0:20:100],'xtick',[1:length(grpAcc)],'xticklabel',gname)
ylabel('% correct')

subtitle({sprintf('Overall accuracy: %.2f (+/-%.2f)',mean(meas),std(meas)./sqrt(length(meas))),...
          sprintf('Correlation with age: r=%.2f, p=%.2f',r,p), ''})
rotateXLabels(gca,30)

% scatter plot age vs d-prime:
[r,p] = corr(dp,age,'type','spearman');
Hagecorr = figure('Name',sprintf('Corr with age'),'position',[0 0 180 150],'color','w'); hold on;
scatter(age,dp,15,'ok','fill')
htmp = lsline; htmp.Color = [0 0 0]; htmp.LineStyle = '--';
set(gca, 'xlim',[min(age)-1 max(age)+1],'ylim',[0 4],'ytick',[0:1:4],'xtick',min(age):2:max(age),'xticklabelrotation',0);
xlabel('Age'); ylabel('Discriminabilty (d'')'); 
subtitle(sprintf('r=%.2f, P>%.2f',r,p));
      
% Supplemental analysis: linear regression of performance and age:
[~,ix] = ismember(gname,demographics.SubjectID); ix(ix==0) = []; assert(all(strcmpi(demographics.SubjectID(ix),gname)));
age = demographics.Age(ix);
tbl = table(grpAcc,age,'VariableNames',{'acc','age'});
lm = fitlm(tbl,'acc ~ 1 + age','RobustOpts','off');



%% Memory effect across hemispheres:
indClipA = contains(DATA.Rall.stimulus_type,'ClipB') & selectedElecInd & ~contains(DATA.Rall.response,'old'); % correct old responses
indClipB = contains(DATA.Rall.stimulus_type,'ClipA') & selectedElecInd & ~contains(DATA.Rall.response,'new'); % correct new responses

[RR1] = grpstats(RRnormTW(indClipA),DATA.Rall.channelid(indClipA),{'nanmean'});
[RR2] = grpstats(RRnormTW(indClipB),DATA.Rall.channelid(indClipB),{'nanmean'});
[m1,s1,elecName1,n1] = grpstats(RRnormTW(indClipA),DATA.Rall.channelid(indClipA),{'mean','var','gname','numel'}); % NEW
[m2,s2,elecName2,n2] = grpstats(RRnormTW(indClipB),DATA.Rall.channelid(indClipB),{'mean','var','gname','numel'}); % OLD

ElecID = {};
ElecAnatSF = {};
ElecHemisphere = {};
for ii = 1:length(elecName1)
    tmp = split(elecName1(ii),'_');
    curElecName = tmp(1);
    subjectName = tmp(2);
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
    ElecID = [ElecID; elecName1(ii)];
end


leftInd = contains(ElecHemisphere,'1');
rightInd = contains(ElecHemisphere,'0');

% Compute hedge's G: (OLD-vs-NEW)
df=n2+n1-2;
sP=((n2-1).*s2 + (n1-1).*s1)./(n2+n1-2); % pooled (within-groups) variance
es=(m2-m1)./sqrt(sP);

[s] = split(ElecID,'_'); s=s(:,2);

% loop through all electrodes and get the corresponding patient's age:
age = []; 
for i = 1:length(s), age(i,1) = demographics.Age(strcmpi(s(i),demographics.SubjectID)); end

% build the table for the mixed effects model:
tbl = table(es,leftInd,age,s,ElecID,'VariableNames',{'effectSize','hemi','age','subject','electrode'});
lmeHemi = fitlme(tbl,'effectSize ~ 1 + hemi + age +  (1|subject:electrode)');
res1 = lmeHemi.Coefficients;

Hhemi = figure('color','w','name','bar plot ripple rate across hemisphere','position',[0 0 100 200]);
superbar([1 2],[mean(es(leftInd)), mean(es(rightInd))],'E',[res1.SE(2), res1.SE(2)],...
    'BarFaceColor',[COLOR.lightgray; COLOR.lightgray],'ErrorbarLineWidth',1,'baredgecolor','none','ErrorbarRelativeWidth',0.5)
set(gca, 'xlim',[0.25 2.75],'ylim',[-0.05 0.5],'ytick',[0:0.1:1],'xtick',[1 2],'xticklabel',{''})
ylabel({'Effect size (OLD-vs-NEW)'})
[p_hemi_ranksum]=ranksum(es(leftInd),es(rightInd));
tmpy=get(gca,'ylim');
line([1;2],[tmpy(2);tmpy(2)].*0.92,'color','k','linewidth',1)
subtitle({sprintf('t(%d)=%.2f, P=%.3f',res1.DF(2),res1.tStat(2),res1.pValue(2)); ''},'fontsize',6,'HorizontalAlignment','center')

text(1.5,tmpy(2)-0.05*range(tmpy),'*','fontsize',6,'HorizontalAlignment','center')
%text(0,tmpy(2)+0.05*range(tmpy),sprintf('age: t(%d)=%.2f\nP=%.2f\n*',res1.DF(3),res1.tStat(3),res1.pValue(3)),'fontsize',6,'HorizontalAlignment','center')
%text(1.5,tmpy(2),sprintf('P<10^{%.0f}\n*',fix(log10(p_hemi))),'fontsize',6,'HorizontalAlignment','center')
text([1,2],repmat(tmpy(1)-0.03*range(tmpy),[1,2]),{sprintf('LH (n=%d)',sum(leftInd)), sprintf('RH (n=%d)',sum(rightInd))},...
    'fontsize',6,'HorizontalAlignment','right','rotation',90);

fprintf('\nHEMISPHERE effect p-value < %.4f (n=%d) \n',ceil(res1.pValue(2)*1000)./1000,length(RR1));%,'fontweight','normal')
fprintf('\nAGE effect: t(%d)=%.2f, P=%.2f \n',res1.DF(3),res1.tStat(3),res1.pValue(3));
set(gca,'ylim',tmpy)

outdir =fullfile(parentfolder,'results/stats/lme/');
writetable(dataset2table(lmeHemi.Coefficients),fullfile(outdir,'lmeHemi.csv'))



%% Group analysis - PSTH time-locked to stimulus onset:
warning('off')
set_figure_colors;
%L1 = unique(DATA.Rall.channelid,'stable');

single_elec_flag = 0; % change this to 1 if you want to only include a single CA1 channel from each patient, otherwise keet it 0 (all electrodes)
if single_elec_flag
    
    % find CA1 electrode indices:
    CA1elecInd = false(size(DATA.Rall.channelid));
    selected_subfield = 'CA1';
    for i = 1:length(subjects)
        FSdir = fullfile(parentfolder,'Freesurfer');
        curSubjectDir = fullfile(parentfolder, subjects{i});
        fid = fopen(fullfile(curSubjectDir,'subjectname'));
        C = textscan(fid,'%s'); fclose(fid);
        FSsubjid = cell2mat(C{1}); % subject FS code
        elecInfoDir=fullfile(FSdir,FSsubjid,'hippocampal_channels');
        tmp = load(fullfile(elecInfoDir,'hippocampalElectrodesInfo.mat'));
        [ind] = tmp.elecinfo.groupDist2sf.(selected_subfield)<2;
        %fprintf('\nDistance from subfield: %f',mm);
        if ~any(ind), continue; end
        closest_contacts = split(tmp.elecinfo.groupLabels(ind),'-');
        if size(closest_contacts,2)>1
            closest_contacts = closest_contacts(:,2);
        else
            closest_contacts = closest_contacts(2);
        end
        [~,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjects{i},ref_flag,0);
        ElecInSubfield = hippocampus_all_channels(startsWith(hippocampus_all_channels,closest_contacts));
        CA1elecInd = CA1elecInd | contains(DATA.Rall.channelid,string(ElecInSubfield)'+"_"+string(subjects{i}));
    end
else
    CA1elecInd = selectedElecInd;
end

%% -------------------------------------------------------------------------
% PLOT viewing Ripples PSTHs:
params = struct;
norm_flag = 0;
raster_type = 2;
%CM =  [ones(5,3); cbrewer('seq','Reds',12); COLOR.bordo];
CM =  [ones(5,3); cbrewer('seq','Greys',12); ];

for analysis_flag = 4 %[1,2,3]
    % ---------------------------------------------------------------------
    % General notes:
    % doublet_flag: % -1 = all ripples; 0 = exclude doublets,triplets etc.; 1 = include only doublets
    % doublet_num: 2 = doublet; 3 = triplets; 4 = ...
    % ---------------------------------------------------------------------
    
    switch analysis_flag
        case 1
            % analysis parameteres:
            params.figure_handel = {'H4a','H4b'};
            params.doublet_flag = -1;
            params.doublet_num = 2;
            params.color1 = COLOR.red;
            params.color2 = COLOR.blue;
            params.textlabel1 = '4s video clips';
            params.textlabel2 = 'rest';
            params.sem_flag = 'acrosssubjects';
            %params.binsize = 100;
            params.statistical_test = 1;
            params.xtick = -10:1:10;
           
            params.ytick = [0:0.1:2];
            params.position = [0 200 200 200];
            params.drawshuffleddistribution = 0;
            params.maskExtraTimePoints = 0;
            % select trials: cond1 vs cond2
            params.nconditions = 1;
            params.conditionToCompare = [1,2];
            params.ind1 = CA1elecInd & contains(DATA.Rall.stimulus_type,'Clip');
            params.ind2 = CA1elecInd & contains(DATA.Rall.stimulus_type,'R');
            params.title1 = sprintf('raster_all_clips');
            params.title2 = sprintf('psth_all_clips');
        case 2
            % analysis parameteres:
            params.figure_handel = {'H4c','H4d'};
            params.doublet_flag = -1;
            params.doublet_num = 2;
            params.color1 = COLOR.blue;
            params.color2 = COLOR.red;
            params.textlabel1 = 'OLD';
            params.textlabel2 = 'NEW';
            
            if single_elec_flag, params.sem_flag = 'acrosselectrodes'; raster_type = 1;
            else,  params.sem_flag = 'acrosssubjects'; end
            %params.binsize = 125;
            params.statistical_test = 1;
            params.xtick = -10:1:10;
         
            params.ytick = [0:0.1:2];
            params.position = [0 200 200 200];
            params.drawshuffleddistribution = 0;
            params.maskExtraTimePoints = 0;
            % select trials: cond1 vs cond2rightInd = contains(ElecHemisphere,'0');

            params.nconditions = 2;
            params.conditionToCompare = [1,2];
            params.ind1 = CA1elecInd & contains(DATA.Rall.stimulus_type,'ClipA') & ~contains(DATA.Rall.response,'new');
            params.ind2 = CA1elecInd & contains(DATA.Rall.stimulus_type,'ClipB') & ~contains(DATA.Rall.response,'old');
            if single_elec_flag
                params.title1 = sprintf('raster_OLD-vs-NEW only CA1');
                params.title2 = sprintf('psth_OLD-vs-NEW only CA1');
            else
                params.title1 = sprintf('raster_OLD-vs-NEW');
                params.title2 = sprintf('psth_OLD-vs-NEW');
            end
        case 3
            % analysis parameteres:
            params.figure_handel = {'H4e','H4f'};
            params.doublet_flag = -1;
            params.doublet_num = 2;
            params.color1 = COLOR.blue;
            params.color2 = COLOR.greendark;
            params.color3 = COLOR.turquoise;
            params.textlabel1 = 'memory test';
            params.textlabel2 = '1st presentation';
            params.textlabel3 = '2nd presentation';
            if single_elec_flag, params.sem_flag = 'acrosselectrodes';
            else,  params.sem_flag = 'acrosssubjects'; end
            %params.binsize = 125;
            params.statistical_test = 1;
            params.xtick = -10:1:10;
            
            params.ytick = [0:0.1:2];
            params.position = [0 200 200 200];
            params.drawshuffleddistribution = 0;
            params.maskExtraTimePoints = 0;
            % select trials: cond1 vs cond2
            params.nconditions = 3;
            params.conditionToCompare = [1,2];
            params.ind1 = CA1elecInd & contains(DATA.Rall.stimulus_type,'ClipA') & ~contains(DATA.Rall.response,'new');
            params.ind2 = CA1elecInd & contains(DATA.Rall.stimulus,'M1ClipA');
            params.ind3 = CA1elecInd & contains(DATA.Rall.stimulus,'M2ClipA');
            if single_elec_flag
                params.title1 = sprintf('raster_OLD-vs-ORIG only CA1');
                params.title2 = sprintf('psth_OLD-vs-ORIG only CA1');
            else
                params.title1 = sprintf('raster_OLD-vs-ORIG');
                params.title2 = sprintf('psth_OLD-vs-ORIG');
            end
            
        case 4        
            raster_type = 1;
            % analysis parameteres:
            params.figure_handel = {'H4g','H4h'};
            params.doublet_flag = -1;            
            params.color1 = COLOR.correct;
            params.color2 = COLOR.error;            
            params.textlabel1 = 'correct';
            params.textlabel2 = 'error';            
            if single_elec_flag, params.sem_flag = 'acrosselectrodes';
            else,  params.sem_flag = 'acrosssubjects'; end
            %params.sem_flag = 'acrosstrials';
          
            %params.binsize = 250;
            params.statistical_test = 1;
            params.xtick = -10:1:10;
            params.ylim = [0 0.55];
            params.ytick = [0:0.1:2];
            params.position = [0 200 200 200];
            params.drawshuffleddistribution = 0;
            params.maskExtraTimePoints = 0;
            % select trials: cond1 vs cond2
            params.nconditions = 2;
            params.conditionToCompare = [1,2];
            
            params.ind1 = CA1elecInd & DATA.Rall.correct==1;
            params.ind2 = CA1elecInd & DATA.Rall.correct==0;
            if single_elec_flag
                params.title1 = sprintf('raster_Correct-vs-Error only CA1');
                params.title2 = sprintf('psth_Correct-vs-Error only CA1');
            else
                params.title1 = sprintf('raster_Correct-vs-Error');
                params.title2 = sprintf('psth_Correct-vs-Error');
            end
    end
    
    if ~isfield(params,'ylim')
        switch norm_flag
            case 0, params.ylim = [0.2 0.65];
            case 1, params.ylim = [-0.5 1.5];
            case 2, params.ylim = [0.1 0.5];
        end
    end
    switch time_locking_event
        case 'stimonset'
            if analysis_flag==4, params.xlim = [-4 4.5];
            else params.xlim = [-2 6];
            end
        case 'rt', params.xlim = [-5 2];
    end
    sig = [];
    ind = DATA.trialtime>=params.xlim(1) & DATA.trialtime<=params.xlim(2);
    CHANNELID = DATA.Rall.channelid;
    SUBID = DATA.Rall.subjid;
    RT = DATA.Rall.RT;
    RT(contains(DATA.Rall.response,'n/a')) = nan;
    data = logical(DATA.Rall.raster)'; % all trials
    testCond = params.conditionToCompare;
    clear R1 R2 r1_sem r2_sem
    
    
    % Doublets analysis (keep/remove doublet events):
    if params.doublet_flag>=0
        dt = 300; % in samples
        data = find_doublets(data,dt,params.doublet_flag,params.doublet_num);
    end
    
    % sort trials by RT (optional):
    sortflag = 1;
    if sortflag
        [~,sortind]=sort(RT);
        data = data(:,sortind);
        CHANNELID = CHANNELID(sortind);
        SUBID = SUBID(sortind);
        RT = RT(sortind);
    else
        sortind = 1:size(data,2);
    end
    
    % prepare mask of timepoint outside the trial by NaNs: (optional)
    if params.maskExtraTimePoints
        trialmask = zeros(size(data));
        assert(size(RT,1)==size(data,2));
        for k = 1:size(RT,1)
            switch time_locking_event
                case 'stimonset',  timeind = DATA.trialtime > params.xlim(1)  & DATA.trialtime <= RT(k)+0.5;
                    trialmask(timeind,k)=1;
                case 'rt', timeind = DATA.trialtime < params.xlim(2) & DATA.trialtime >= -(RT(k)+0.5);
                    trialmask(timeind,k)=1;
            end
        end
        data(~trialmask) = 0;
    else, trialmask = ones(size(data));
    end
    
    % GLOBAL vars:
    I = 1000; % number of randomizations
    t = DATA.trialtime(ind);
    % initialize data structures:
    N = struct; N.cond0=0;
    clear condData condChannelid condSubid condRT R r C X Y SEM HL sortedSubjectList
    for k=1:params.nconditions
        condData.(sprintf('cond%d',k)) = [];
        condChannelid.(sprintf('cond%d',k)) = [];
        condSubid.(sprintf('cond%d',k)) = [];
        condRT.(sprintf('cond%d',k)) = [];
        sortedSubjectList.(sprintf('cond%d',k)) = {};
        R.(sprintf('cond%d',k)) = [];
        r.(sprintf('cond%d',k)) = [];
        C.(sprintf('cond%d',k)) = [];
        N.(sprintf('cond%d',k)) = [];
        X.(sprintf('cond%d',k)) = [];
        Y.(sprintf('cond%d',k)) = [];
        SEM.(sprintf('cond%d',k)) = [];
        HL.(sprintf('cond%d',k)) = [];
    end
    
    % Draw ripples raster:
    eval(sprintf('%s=figure;',params.figure_handel{1}));
    set(gcf,'Name',params.title1,'position',params.position+[400 0 0 0],'color','w','renderer','painter');
    hold on;
    NN = [];
    for k=1:params.nconditions
        trialInd = params.(sprintf('ind%d',k))(sortind);
        condData.(sprintf('cond%d',k)) = data(ind,trialInd);
        condChannelid.(sprintf('cond%d',k)) = CHANNELID(trialInd);
        condSubid.(sprintf('cond%d',k)) = SUBID(trialInd);
        condRT.(sprintf('cond%d',k)) = RT(trialInd);
        N.(sprintf('cond%d',k)) = 1:size(condData.(sprintf('cond%d',k)),2);  % trial number
        N.(sprintf('cond%d',k)) = N.(sprintf('cond%d',k))+N.(sprintf('cond%d',k-1))(end);
        TM.(sprintf('cond%d',k)) = trialmask(ind,trialInd);
        % prepare raster plot -
        [X.(sprintf('cond%d',k)),Y.(sprintf('cond%d',k))]=find(condData.(sprintf('cond%d',k)));
        
        % Estimate the optimal bin width (following fieldtrip's ft_spike_psth function):
        % [based on the std of ripples' onset in the current condition]
        Nr = round(sum(condData.(sprintf('cond%d',k))(:))); % ripple count in the current condition
        sd = std(t(X.(sprintf('cond%d',k)))); % std of the current condition
        trial_duration = (max(t)-min(t));
        binsize_method = 'scott';
        bs = optimal_binsize_estimation(Nr,sd,trial_duration,Fs,binsize_method);
        smt = 4; % number of successive time bins to smooth (two runs of moving average: n*smt-n+1 = 6 points filter, where n=number of smooth iterations)
        fprintf('\n*** binsize needs to be approximately %d ms (before smoothing) ***\n',round(bs/smt)/Fs*1000)
        binsize = round(bs/smt); % samples
        
        %%% ploting raster using scatter:
        
        switch raster_type
            case 1 % use scatter
                s = 5; % marker size
                scatter(t(X.(sprintf('cond%d',k))),N.(sprintf('cond%d',k))(Y.(sprintf('cond%d',k))),s,'.',...
                    'markerFaceColor',params.(sprintf('color%d',k)),'markerEdgeColor',params.(sprintf('color%d',k))); hold on;
            case 2 % use densitiy plot
                [NbinX,NbinY] = size(condData.(sprintf('cond%d',k)));
                NtimepointsPerBin = NbinX * 0.025;
                NtrialsPerBin = NbinY * 0.025;
                [nn,xx,yy] = histcounts2(t(X.(sprintf('cond%d',k)))',N.(sprintf('cond%d',k))(Y.(sprintf('cond%d',k)))',...
                    [floor( NbinX / NtimepointsPerBin ),floor( NbinY / NtrialsPerBin )]);
                binsizeX = xx(2)-xx(1);  binsizeY = yy(2)-yy(1);
                % smooth (optional):
                ff = fspecial('gaussian',5,2);
                nn = imfilter(nn, ff, 'replicate');
                imagesc(xx+binsizeX/2, yy+binsizeY/2, nn');
                foo = nn(xx<0,:); % gather the prestimulus counts for adjusting the colorscale later
                NN = cat(1,NN,foo(:));
                colormap(CM);
                selectedElecInd = true(size(DATA.Rall.channelid)); % optional variable to limit the analysis to a subset of electrodes
        end
        
        axis tight
        hold on;
        drawnow;
        
        % add markers at trial onset / resopnse:
        if strcmpi(time_locking_event, 'rt'), scatter(-condRT.(sprintf('cond%d',k)),N.(sprintf('cond%d',k)),1,'o','markerFaceColor',COLOR.black,'markerEdgeColor','none'); hold on;
        else, scatter(condRT.(sprintf('cond%d',k)),N.(sprintf('cond%d',k)),1,'o','markerFaceColor',COLOR.black,'markerEdgeColor','none'); hold on; end
        
        plot([0 0],get(gca,'ylim'),'color',COLOR.black,'linewidth',0.5)
        axis tight;
        set(gca,'xtick',params.xtick,'TickDir','out','xlim',params.xlim,'ytick',[],'yticklabel',{'',''});
        switch time_locking_event
            case 'stimonset', xlabel(sprintf('Time from clip onset (s)'));
            case 'rt', xlabel(sprintf('Time from response onset (s)'));
        end
        
        ylabel('Single trials (pooled across elec.)');
        
        if k~=1
            if raster_type == 1, plot(get(gca,'xlim'),[N.(sprintf('cond%d',k))(1) N.(sprintf('cond%d',k))(1)],'k--','linewidth',1);
            else, plot(get(gca,'xlim'),[yy(1) yy(1)],'k--','linewidth',1);
            end
        end
        
        % add condition labels (optional):
        tmpx=get(gca,'xlim'); tmpy=get(gca,'ylim');
        yText = (N.(sprintf('cond%d',k))(1)+N.(sprintf('cond%d',k))(end))/2;
        text(1.1*tmpx(2),yText,params.(sprintf('textlabel%d',k)),...
            'fontsize',8,'Rotation',90,'horizontalalignment','center','color',params.(sprintf('color%d',k)),'fontweight','bold');
        pos = get(gca, 'Position');
                
    end
    
    if raster_type == 2
   
        clim = [floor(prctile(NN,2.5)/5)*5 ceil(prctile(NN,97.5)/5)*5];
        caxis(clim);
        Hcb = figure('name',sprintf('colorbar_analysis_%d',analysis_flag),'position',[0 0 100 100],'color','w');
        cb=colorbar('location','south'); caxis(clim); colormap(CM);
        cb.Label.String = 'SWR count';
        cb.Limits = clim; cb.Ticks = clim;
        axis off
        
    end

    
    % =====================================================================
    % PSTH:
    eval(sprintf('%s=figure;',params.figure_handel{2}));
    set(gcf,'Name',params.title2,'position',params.position,'color','w','renderer','painter'); hold on;
    hold on;
    
    % compute single trial ripple rate for baseline:
    timeWindow = [0 4];
    [RR,trialDur] = compute_single_trial_ripple_rate(DATA,Fs,timeWindow,time_locking_event);
    
    current_rest_baseline = [];
    for k=1:params.nconditions
        clear L1 L2
        [r.(sprintf('cond%d',k)),binscenters]=nanPSTH(condData.(sprintf('cond%d',k)),binsize,smt,Fs, []);
        binscenters = (binscenters/1000)+t(1);
        
        % compute the mean of each channel:
        [L1,iii,foo] = unique(condChannelid.(sprintf('cond%d',k)));
        [L2] = condSubid.(sprintf('cond%d',k))(iii);
        
        % compute ripple rate separatly in each channel:
        C.(sprintf('cond%d',k)) = [];
        for l = 1:length(L1)
            channelid = L1(l);
            currentdata = condData.(sprintf('cond%d',k))(:,strcmpi(condChannelid.(sprintf('cond%d',k)),channelid));
            currentmask = TM.(sprintf('cond%d',k))(:,strcmpi(condChannelid.(sprintf('cond%d',k)),channelid));
            [currentChRate]=nanPSTH(currentdata,binsize,smt,Fs, currentmask);
            % ripple rate normalization within electrode:
            tmpind = strcmpi(DATA.Rall.channelid,channelid);
            switch norm_flag
                case 0, blavg = nanmean(RR(tmpind & contains(DATA.Rall.stimulus_type,'R')));
                case 1  % zscore
                    blavg = nanmean(RR(tmpind & contains(DATA.Rall.stimulus_type,'R')));
                    blstd = nanstd(RR(tmpind & contains(DATA.Rall.stimulus_type,'R')));
                    currentChRate = (currentChRate-blavg)./blstd;
                    blavg = blavg.*0;
                case 2 % minmax
                    blmax = max(RR(tmpind));
                    blmin = min(RR(tmpind));
                    blavg = nanmean(RR(tmpind & contains(DATA.Rall.stimulus_type,'R')));
                    currentChRate = (currentChRate-blmin)./(blmax-blmin);
                    blavg = (blavg-blmin)./(blmax-blmin);
            end
            C.(sprintf('cond%d',k))(l,:) = currentChRate;
            current_rest_baseline(l,:) = blavg;
        end
        % average within patient:
        [R.(sprintf('cond%d',k)),sortedSubjectList.(sprintf('cond%d',k))] = grpstats(C.(sprintf('cond%d',k)),L2,{'mean','gname'});
        
        % compute rest baseline per each patient:
        if k==1
            rest_trials_ix = contains(DATA.Rall.stimulus_type,'R');
            [rest_baseline,rest_subjid] = grpstats(RR(rest_trials_ix),DATA.Rall.subjid(rest_trials_ix),{'mean','gname'});
            R.('rest_baseline') = repmat(rest_baseline,[1,size(R.(sprintf('cond%d',testCond(1))),2)]);
        end
    end
    
    % For within-subjects analysis, exclude patient with not enough trials / unpaired data:
    if  params.nconditions>1 && (params.statistical_test==1 ||  params.statistical_test==3)
        validSubjects = intersect(sortedSubjectList.(sprintf('cond%d',testCond(1))),...
                                  sortedSubjectList.(sprintf('cond%d',testCond(2))));
        [~,included_ix1] = ismember(validSubjects,sortedSubjectList.(sprintf('cond%d',testCond(1))));
        [~,included_ix2] = ismember(validSubjects,sortedSubjectList.(sprintf('cond%d',testCond(2))));
        [~,included_ix3] = ismember(validSubjects,rest_subjid);
        included_ix1(included_ix1==0)=[];
        included_ix2(included_ix2==0)=[];
        included_ix3(included_ix3==0)=[];
        assert(all(strcmpi(sortedSubjectList.(sprintf('cond%d',testCond(1)))(included_ix1), ...
                           sortedSubjectList.(sprintf('cond%d',testCond(2)))(included_ix2))))
        assert(all(strcmpi(sortedSubjectList.(sprintf('cond%d',testCond(1)))(included_ix1),rest_subjid(included_ix3))));
        R.(sprintf('cond%d',testCond(1))) = R.(sprintf('cond%d',testCond(1)))(included_ix1,:);
        R.(sprintf('cond%d',testCond(2))) = R.(sprintf('cond%d',testCond(2)))(included_ix2,:);
        R.rest_baseline = R.rest_baseline(included_ix3,:);
    end
    
    % Compute SEM and plot the average response:
    AVG = struct;
    for k=1:params.nconditions
        % SEM over trails / subjects / withinsubject differences / electrodes:
        sem_flag = params.sem_flag;
        if strcmpi(params.sem_flag,'withinsubjects') & ~ismember(k,params.conditionToCompare), sem_flag = 'acrosssubjects'; end
        
        switch sem_flag
            case 'acrosssubjects'
                % SEM across subjects
                AVG.(sprintf('cond%d',k)) = nanmean(R.(sprintf('cond%d',k)),1);
                SEM.(sprintf('cond%d',k)) = nanstd(R.(sprintf('cond%d',k)))./sqrt(sum(~all(isnan(R.(sprintf('cond%d',k))),2)));
            case 'withinsubjects'
                % SEM for within-subject comparison:
                AVG.(sprintf('cond%d',k)) = nanmean(R.(sprintf('cond%d',k)),1);
                SEM.(sprintf('cond%d',k)) = [];
                for ii = 1:size(binscenters,2)
                    [~,~,SEM.(sprintf('cond%d',k))(:,ii),nn] = ci_within_subject(R.(sprintf('cond%d',params.conditionToCompare(1)))(:,ii),...
                                                                                 R.(sprintf('cond%d',params.conditionToCompare(2)))(:,ii));
                end
            case 'acrosselectrodes'
                % SEM across subjects
                AVG.(sprintf('cond%d',k)) = nanmean(C.(sprintf('cond%d',k)),1);
                SEM.(sprintf('cond%d',k)) = nanstd(C.(sprintf('cond%d',k)))./sqrt(sum(~all(isnan(C.(sprintf('cond%d',k))),2)));
            case 'acrosstrials'
                AVG.(sprintf('cond%d',k)) = r.(sprintf('cond%d',k));
                % bootstrap SEM across all trials pooled together from all subjects:
                SEM.(sprintf('cond%d',k)) = [];
                opt = statset('UseParallel',true);
                SEM.(sprintf('cond%d',k)) = bootstrp(100,@(x)nanPSTH(x,binsize,smt,Fs,[]),condData.(sprintf('cond%d',k)),'Options',opt);
                SEM.(sprintf('cond%d',k)) = std(SEM.(sprintf('cond%d',k)));
                
        end
        hold all;
        Htmp = [];
        % plot the response:
        if sum(~all(isnan(R.(sprintf('cond%d',k)))))>1
            Htmp = shadedErrorBar(binscenters,AVG.(sprintf('cond%d',k)),SEM.(sprintf('cond%d',k)),...
                {'color',params.(sprintf('color%d',k)),'linewidth',1,'linesmoothing','on'},1);
        end
        HL.(sprintf('cond%d',k)) = Htmp.mainLine;
    end
    
    axis tight
    tmpy = get(gca,'ylim');
    set(gca,'xtick',params.xtick,'ytick',params.ytick,'TickDir','out','xlim',params.xlim);
    if isfield(params,'ylim'), ylim([params.ylim(1) params.ylim(2)]);
    else,  ylim([0 tmpy(2)]);
    end
    
    switch time_locking_event
        case 'stimonset', xlabel(sprintf('Time from clip onset (s)'));
        case 'rt', xlabel(sprintf('Time from response onset (s)'));
    end
    
    set(gca, 'Position',pos)
    
    % =========================================================================
    % STATS:
    hax = gca;
   
       
    if params.statistical_test
        % baseline-vs-response comparison across patients:
        if analysis_flag == 1 || analysis_flag == 4 
            
            for k=1:params.nconditions                
                cond1 = R.(sprintf('cond%d',testCond(k)))'; % video clip resposnse
                cond2 = R.rest_baseline'; % rest baseline
                
                % for comparison with prestim baseline: (optional)
                % baselineWin = binscenters<0;
                % prestimbl = mean(R.(sprintf('cond%d',testCond(1)))(:,baselineWin),2)';
                % cond2 = repmat(bl,[length(binscenters),1]);
                
                [sig, pval, onsets, offsets] = compare_conditions(hax,cond1,cond2,binscenters,params.statistical_test,1,params.(sprintf('color%d',k)));
                if any(sig)
                    sigTimeWindow = [onsets(1); offsets(1)];
                end
                if ~single_elec_flag && analysis_flag == 1, save(fullfile(parentfolder,'results','stats','sigTimeWindow.mat'),'sigTimeWindow'); end
            end
        end
        
        % condition1-vs-condition2 comparison across patients:
        if analysis_flag > 1
            if single_elec_flag
                % comparison across electrodes
                cond1 = C.(sprintf('cond%d',testCond(1)))';
                cond2 = C.(sprintf('cond%d',testCond(2)))';
            else
                
                % comparison across subjects
                cond1 = R.(sprintf('cond%d',testCond(1)))';
                cond2 = R.(sprintf('cond%d',testCond(2)))';
                
            end
            [sig, pval, onsets, offsets] = compare_conditions(hax,cond1,cond2,binscenters,params.statistical_test,1);
                        
        end
    end
    
    handels = [];
    LEGENDS = {};
    for k = 1:params.nconditions
        LEGENDS{k,1} = params.(sprintf('textlabel%d',k));
        handels=cat(1,handels, HL.(sprintf('cond%d',k))(1));
    end
    
    tmpy=get(gca,'ylim');
    set(gca,'ylim',[tmpy(1) tmpy(2)*1.1]);
    
    % draw baseline (optinoal):
    h_baseline = plot(get(gca,'xlim'),[mean(R.rest_baseline(:,1)) mean(R.rest_baseline(:,1))],'--','color',COLOR.black,'linewidth',1);
    if norm_flag==0
        disp('no-normalization, showing ripple rate...');
        ylabel('Ripple rate (events/s)');
    elseif norm_flag==1, ylabel('Ripple rate modulation (fold-change)');
    elseif norm_flag==2, ylabel('Ripple rate (n.u)');
    end
    LEGENDS = cat(1,LEGENDS,{'resting state'});
    handels = cat(1,handels,h_baseline);
    
    % add legends (optional):
    [h_leg,h_leg_icon] = legend(handels,LEGENDS,'location','NorthWest');
    pos = h_leg.Position;
    h_leg.Position = [pos(1)*0.5 0.8 pos(3:4)];
    h = findobj(h_leg_icon,'type','line');
    for ii = 1:length(h), if length(h(ii).XData)==2, h(ii).XData(1) = h(ii).XData(1) + 0.5*range(h(ii).XData); end; end
    legend boxoff
    
    % add stimulus on mark:
    tmpx=get(gca,'xlim');
    tmpy=get(gca,'ylim'); 
%     set(gca,'ylim',[tmpy(1)- 0.05*range(tmpy) tmpy(2)]);
%     tmpy=get(gca,'ylim'); 
    line([0 4],[tmpy(1) tmpy(1)]+0.01*range(tmpy),'color',COLOR.black,'linewidth',5)
    text(tmpx(2),tmpy(2)*0.9,sprintf('n=%d',size(cond1,2)),'fontsize',6,'HorizontalAlignment','right');
end



%% SAVE:
%==========================================================================
% set path:
outdir=fullfile(parentfolder,'results','ripple_rate_PETH_final',sprintf('ref_%d_%s',ref_flag,time_locking_event));
if ~exist(outdir,'dir'), mkdir(outdir);disp('Creating Output Directory...'); end
%=========================================================================
for F = [H4a H4b H4c H4d H4e H4f]; %[Hagecorr Hbehavior Hhemi H4a H4b H4c H4d H0 H1a H1b H1c H1d Hsc Hcb] 
    figure(F);  
    set_font_size_and_type;
    enlarge_figure_and_move_axes_to_center(gcf,gca,1.1);
    save_current_figure(F,outdir,1); 
end

