% The script loads the SWR time stamps across all hippocampal sites
% (saved as rasters) and runs the following analyses:
% (1) compute the mean ripple rate during each video clip
% (2) ripple rate map
% (3) ripple rate old-vs-new effect size map (Figure 3)
% (4) ripple rate across subfields (mixed effects supplemental analysis)
% note: the SWR raster data is constructed beforehand using the script:
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

demographics = readtable(fullfile(parentfolder,'subj_info_anon.xlsx'));

% Load rasters and initalize general varaiables:
load_swr_rasters_from_all_hippocampal_channels;


%% Plot Hippocampal electrodes:
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

%% select effect size measure: 
analysis_flag = 1;
switch analysis_flag
    case 1 % old - new
        ind0 = contains(DATA.Rall.stimulus_type,'R');
        ind1 = contains(DATA.Rall.stimulus_type,'ClipA') & ~contains(DATA.Rall.response,'new');
        ind2 = contains(DATA.Rall.stimulus_type,'ClipB') & ~contains(DATA.Rall.response,'old');
        [m0,s0,elecName0,n0] = grpstats(RRnorm4s(ind0),DATA.Rall.channelid(ind0),{'mean','var','gname','numel'});
        [m1,s1,elecName1,n1] = grpstats(RRnormTW(ind1),DATA.Rall.channelid(ind1),{'mean','var','gname','numel'});
        [m2,s2,elecName2,n2] = grpstats(RRnormTW(ind2),DATA.Rall.channelid(ind2),{'mean','var','gname','numel'});
        % Compute hedge's G:
        df=n1+n2-2;
        sP=((n1-1).*s1 + (n2-1).*s2)./(n1+n2-2); % pooled (within-groups) variance
        es=(m1-m2)./sqrt(sP);
        % p-values:
        pp = []; counter = 1;
        for currentElec = unique(elecName1,'stable')'
            [pp(counter,1)] = ranksum(RRnormTW(ind1&strcmpi(DATA.Rall.channelid,currentElec)), RRnormTW(ind2&strcmpi(DATA.Rall.channelid,currentElec)));
            counter = counter+1;
        end
        [~, crit_p, ~, adj_p]=fdr_bh(pp,0.05,'pdep','yes');
        adj_p = pp;
        fprintf('\nnumber of electrodes showing significant preference to old clips: %d\n',sum(adj_p(es<0)<0.05))
        fprintf('\nnumber of electrodes showing significant preference to old clips: %d\n',sum(adj_p(es>0)<0.05)) 
end

% resort the electrode data:
assert(all(strcmpi(elecName0,elecName1)&strcmpi(elecName1,elecName2)));
includedRows = []; includedSubs = {};
for ii = 1:length(elecName1)
    currentElec = split(elecName1(ii),'_');
    ch=split(currentElec{1},'-');
    subjdir=fullfile(parentfolder, currentElec{2});
    fid = fopen(fullfile(subjdir,'subjectname'));
    C = textscan(fid,'%s'); fclose(fid);
    FSsubjid = cell2mat(C{1}); % subject FS code
    currentElec = strcat(FSsubjid,'-',ch{1});
    includedRows(ii,1) = find(strcmpi(elecinfo.groupLabels,currentElec));
    includedSubs{ii,1} = FSsubjid;   
end

% generate a table with all the anatomical details:
[~,~,includedSubsInd] = unique(includedSubs,'stable');
[predictorLabel,includedCols] = setdiff(elecinfo.groupDist2sf.Properties.VariableNames,{'xxx'},'stable'); % in case you want to exclude one of the subfields
X = elecinfo.groupDist2sf(includedRows,includedCols);
X(elecinfo.groupIsSFparc(includedRows)==0,:)=array2table(nan(sum(elecinfo.groupIsSFparc(includedRows)==0),size(X,2)),'VariableNames',predictorLabel);

%% Plot the ripple rate map:

[elecinfo,Hrate] = getHippocampalElectrodeCoord(L1, subjects,0,1);
set(Hrate,'name',sprintf('ripple-rate map analysis %d',analysis_flag))
set(elecinfo.Helec,'visible','off')
Helec = elecinfo.Helec(includedRows);
set(Helec,'visible','on')

DD = m0; 
clim=[0.2 0.6]; CM = load(fullfile(parentfolder,'Freesurfer','mycmappos.mat'));   
CM = CM.mycmap; 
ncolors = size(CM,1);
thr = max(0,clim(1));
thrColorInd1 = floor(interp1(clim,[1 ncolors],max(-thr,clim(1)),'linear',ncolors));
thrColorInd2 = floor(interp1(clim,[1 ncolors],min(thr,clim(2)),'linear',ncolors));
tmp = repmat(COLOR.lightgray,thrColorInd2-thrColorInd1,1);
CM = [CM(1:thrColorInd1-1,:); tmp; CM(thrColorInd2+1:end,:)];
DD(abs(DD)<thr) = 0; DD(DD>clim(2))=clim(2); DD(DD<clim(1))=clim(1);  

ncolors = size(CM,1);
for ii = 1:size(DD,1)
    colind = floor(interp1(clim,[1 ncolors],DD(ii),'linear',ncolors));
    set(Helec(ii),'facecolor',CM(colind,:));
end
fprintf('\n number of electrodes: %d\n',length(Helec))

% color bar:
Hcb1 = figure('color','w','name',sprintf('colorbar ripple rate analysis %d',analysis_flag),'position',[0 0 100 100]);
cb=colorbar('location','south'); caxis(clim)
colormap(CM);
cb.Label.String = 'SWR rate (events/sec)';
cb.Limits = clim; cb.Ticks = [clim(1),clim(2)];
axis tight
axis off
set(gca,'XColor','none','Ycolor','none')

%% Plot the old-new effect size map:
[elecinfo,Hes] = getHippocampalElectrodeCoord(L1, subjects,0,1);
set(Hes,'name',sprintf('effectsize anlysis %d',analysis_flag))
set(elecinfo.Helec,'visible','off')
Helec = elecinfo.Helec(includedRows);
set(Helec,'visible','on')

DD = es; 
switch analysis_flag
    case 1
        clim=[0 0.7]; thr = 0.1;
        CM = cbrewer('seq','YlOrBr',8); %(cbrewer('seq','OrRd',6));
end

DD(abs(DD)<thr) = thr; DD(DD>clim(2))=clim(2); DD(DD<clim(1))=clim(1);  

% For adding a threshold color:
ncolors = size(CM,1);
thrColorInd1 = floor(interp1(clim,[1 ncolors],max(-thr,clim(1)),'linear',ncolors));
thrColorInd2 = floor(interp1(clim,[1 ncolors],min(thr,clim(2)),'linear',ncolors));
tmp = repmat(COLOR.lightgray,max(thrColorInd2-thrColorInd1,1),1);
CM = [CM(1:thrColorInd1-1,:); tmp; CM(thrColorInd2+1:end,:)];

% clear tmp
ncolors = size(CM,1);
for ii = 1:size(DD,1)
    colind = floor(interp1(clim,[1 ncolors],DD(ii),'linear',ncolors));
    set(Helec(ii),'facecolor',CM(colind,:));
end
fprintf('\n number of electrodes: %d\n',length(es))

% color bar:
Hcb2 = figure('color','w','name',sprintf('colorbar effect size analysis %d',analysis_flag),'position',[0 0 100 100]);
cb=colorbar('location','south'); caxis(clim)
colormap(CM);
cb.Label.String = sprintf('Old > New\n effect size (g)');
cb.Limits = clim; cb.Ticks = [clim(1),clim(2)];
%cb.TickLabels = {'new','old'};

axis tight
axis off
set(gca,'XColor','none','Ycolor','none')

%% Histogram of old-vs-new preferences: (supplemental analysis)
Heshist = figure('color','w','name','memory-type preference histogram','position',[0 0 120 180]);
[counts,bs] = hist(es,[-0.6:0.1:0.6]);
counts1 = counts; counts1(bs<0.05) = nan;
counts2 = counts; counts2(bs>-0.05) = nan;
counts0 = counts; counts0(abs(bs)>0.05) = nan;
hb0=superbar(bs,counts0,'barfacecolor',COLOR.gray,'Orientation','h','BarWidth',0.0985,'BarEdgeColor','none'); hold on;
hb2=superbar(bs,counts2,'barfacecolor',COLOR.blue,'Orientation','h','BarWidth',0.0985,'BarEdgeColor','none');
hb1=superbar(bs,counts1,'barfacecolor',COLOR.red,'Orientation','h','BarWidth',0.0985,'BarEdgeColor','none');
set([hb0; hb1; hb2],'facealpha',0.5)
set(gca,'ylim',[-0.8 0.8],'ytick',[-1:0.2:1])
ylabel({'Effect size (Hedges g)'; '\leftarrow new    old \rightarrow'})
xlabel('Electrode count')


%% multiple linear regression of distances from hippocampal subfields (supplemental analysis):
hemisphere = repmat({''},size(includedRows));
hemisphere(find(elecinfo.groupIsLeft(includedRows)==1)) = {'l'};
hemisphere(find(elecinfo.groupIsLeft(includedRows)==0)) = {'r'};
anatLoc = repmat({''},size(includedRows));
anatLoc(find(elecinfo.groupIsAnt(includedRows)==1)) = {'a'};
anatLoc(find(elecinfo.groupIsAnt(includedRows)==0)) = {'p'};

tmp = split(L1,'_'); s=tmp(:,2);
[~,ageind]  = ismember(s,demographics.SubjectID);
patientAge = demographics.Age(ageind);
% stats data table - recognition effect size:
T = [table(es,'VariableNames',{'effectSize'}),...
     table(m0,'VariableNames',{'basalRate'}),...
     table(m1,'VariableNames',{'oldRate'}),...
     table(m2,'VariableNames',{'newRate'}),...
     table(patientAge,'VariableNames',{'age'}),...
     table(L1,'VariableNames',{'electrode'}),...
     table(num2str(includedSubsInd),'VariableNames',{'subject'}),...
     table(hemisphere,'VariableNames',{'Hemisphere'}), ...
     table(anatLoc,'VariableNames',{'AnatomicalLocation'})];

T = [T X];
lme_anatomy = fitlme(T,'effectSize ~ 1 + CA1 + CA2_3 + subiculum + GCDG + AnatomicalLocation + Hemisphere + AnatomicalLocation*Hemisphere + (1|subject)');

% Export to R (optional):
mkdir(fullfile(parentfolder,'results','stats'));
switch analysis_flag
    case 1, csvFileName = fullfile(parentfolder,'results','stats','anatomicalModel_old-vs-new.csv');
end
writetable(T,csvFileName,'WriteRowNames',false)


%% SAVE FIGURES: 
%set path:
outdir=fullfile(parentfolder,'results','ripple_rate_PETH_final',sprintf('ref_%d_%s',ref_flag,time_locking_event));
if ~exist(outdir,'dir')
    mkdir(outdir);
    disp('Creating Output Directory...')
end
%=========================================================================
for F = [Hes Hcb2, Heshist, Hrate, Hcb1]
    figure(F);
    set_font_size_and_type;      
    save_current_figure(F,outdir,1);
end

