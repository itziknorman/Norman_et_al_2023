%% Set Path:
clear all
close all
clc;


% set the relevant path:
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;

% add/remove paths:
addpath(fullfile(path_to_toolboxes,'eeglab2021.1'));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));
rmpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));

[ALLEEG, EEG, CURRENTSET] = eeglab;


subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};
stats =  table([],[],[],[],[],'VariableNames',{'electrode','subjid','latency','minpval','sumtval'});
SS = [];
FF = [];

for subjid = subjects
    subjid = cell2mat(subjid);
    %======================================================================
    ref_flag=2; % Adjust Manually (1=common; 2=bipolar montage)
    %======================================================================
    close all;
    clearvars -except subjects subjid run ALLEEG EEG ref_flag Mall Rall M1all M2all R1all R2all MTall parentfolder path_to_toolboxes stats SS FF
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    % Load Datasets:
    scenario={'R1','M1','R2','M2','MT'};
    maindir=fullfile(parentfolder,subjid);
    switch ref_flag
        case 1         
            datadir = fullfile(maindir,'EEGLAB_datasets');
            HFBmatdir = fullfile(maindir,'HFB_MAT');     
            
        case 2            
            datadir = fullfile(maindir,'EEGLAB_datasets_BP');
            HFBmatdir= fullfile(maindir,'HFB_MAT_BP');            
    end
    
    % create a new fig folder:
    figdir =  fullfile('D:\ECoG\pink_panther_anonymized\results\spectrograms\',sprintf('ref_%d',ref_flag),subjid); % ADJUST OUTDIR 
    if ~exist(figdir,'dir')
        mkdir(figdir);
        disp('Creating Output Directory...')
    end
    groupfigdir =  fullfile('D:\ECoG\pink_panther_anonymized\results\spectrograms-group-results');
    if ~exist(groupfigdir,'dir')
        mkdir(groupfigdir);
        disp('Creating Output Directory...')
    end
    % define filenames:
    cd(datadir);
    
    % LOAD EEGLAB EEG DATASETS:
    for s=1:numel(scenario)
        filename=dir(['*' scenario{s} '*.set']);
        filename=filename(1).name
        [EEG] = pop_loadset('filename', filename);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
        eeglab redraw;
        
        % load HFB: =======================================================
        %         load(fullfile(HFBmatdir,[subjid '_' scenario{s} '_preprocessed_HFB.mat']));
        %         % dB transformation:
        %         EEG.data=DATA.HFB.amplitude;
        %         EEG.data=bsxfun(@rdivide,EEG.data,nangeomean(nangeomean(EEG.data,2),3));
        %         EEG.data=10.*log10(EEG.data);
        %         % smoothing (triangular window):
        %         for i=1:size(EEG.data,1)
        %         EEG.data(i,:)=nanfastsmooth(EEG.data(i,:),25,3,0.5);
        %         end
        %==================================================================
        % [EEG, ~,b] = pop_firws(EEG, 'fcutoff', 20, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', 826, 'minphase', 0);
        
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end
    
    % epoch memory test:
    idx=find(multiStrFind({EEG.event.type},'Clip'));
    EEG = pop_epoch( EEG,{}, [-1  7],'eventindices',idx,'verbose','no');
    EEG = pop_rmbase( EEG, [-500    0]);
    EEG.setname = ['Movie_clips_epoched'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); eeglab redraw;
    
    % load good channels list:
    switch ref_flag
        case 1
            load(fullfile(maindir,[subjid '_channels_list']))
        case 2
            good_channels=find(multiStrFind({EEG.chanlocs.type},'signal'));
    end
        
    % epoch memory test:
    idx=find(multiStrFind({EEG.event.type},'Clip'));
    EEG = pop_epoch( EEG,{}, [-1  7],'eventindices',idx,'verbose','no');
    EEG = pop_rmbase( EEG, [-500    0]);
    EEG.setname = ['Movie_clips_epoched'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); eeglab redraw;
    
    % load good channels list:
    switch ref_flag
        case 1
            load(fullfile(maindir,[subjid '_channels_list']))
        case 2
            good_channels=find(multiStrFind({EEG.chanlocs.type},'signal'));
    end
    
    %% gather all hippocampal channels:
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0); % Define hippocampus channels
    hpchan = find(contains({EEG.chanlocs.labels},hippocampus_all_channels));
     
    tval = []; pval = []; 
    H =  figure('color','w','position',[0 0 300 200],'name',sprintf('%s hippocampal erps during video clip',subjid)); hold on;
    for i = 1:length(hpchan)
        t = ALLEEG(6).times;
        
        [S1,~] = spectopo(ALLEEG(1).data(hpchan(i),:,:),ALLEEG(1).pnts,ALLEEG(6).srate,'winsize',512,'plot','off','rmdc','on','freqfac',4); 
        [S2,F] = spectopo(ALLEEG(3).data(hpchan(i),:,:),ALLEEG(3).pnts,ALLEEG(3).srate,'winsize',512,'plot','off','rmdc','on','freqfac',4);
        
        ind = ismember(F,1:0.5:256/2);
        S = (S1(ind)+S2(ind))./2;
        F = F(ind);        
        SS = cat(1,SS,S);
        FF = F;  
        tmp = squeeze(ALLEEG(6).data(hpchan(i),:,:)); 
        idx_baseline = t<0;
        idx_stimulus = t>0 & t<4000;       

        % compare variance during the clip with baseline:
        L0 = mean(tmp(idx_baseline,:));
        L1 = tmp(idx_stimulus,:)-L0;
        
        [~,p,~,st] = ttest(L1');
        pval(i,:) = p;
        tval(i,:) = st.tstat;
        lfp(i,:) = mean(tmp(idx_stimulus,:),2);
    end    

    %[~, crit_p, ~, adj_p]=fdr_bh(pval,0.05,'pdep','yes'); % for FDR correctrion
    crit_p = 0.05;     
    clustersize = ceil(0.050 * ALLEEG(6).srate); % set the min cluster size at 50 ms
      
    
    figure(H); hold on;
    for i = 1:size(tval,1)
        elecname = ALLEEG(6).chanlocs(hpchan(i)).labels;
        idx = (pval(i,:) < crit_p);
        idx = sprintf('%d', idx);
        [onset,offset] = regexp(idx, sprintf('1{%d,}',clustersize), 'start'); % find the first time point that is sig for at least 50 ms
        k = find(onset,1,'first');
        onset = onset(k);        
        offset = offset(k);   
        dur = offset-onset;
        latency = nan; minp = nan; sumt = nan;
        if ~isempty(onset) 
            onset=onset(1);
            latency=t(idx_stimulus); latency=latency(onset);
            minp = min(pval(i,onset:offset));
            sumt = sum(tval(i,onset:offset));
            fprintf('\n-- Onset latency channel %s: %.1f ms (duration: %.2f s)', elecname,latency,dur./ALLEEG(6).srate)            
            plot(t(idx_stimulus),lfp(i,:),'k'); axis tight
            hs = scatter(latency,lfp(i,onset),20,'y','fill','markeredgecolor','k');
        else
            plot(t(idx_stimulus),lfp(i,:),'b'); axis tight;
        end      
        xlabel('Time from clip onset (ms)');
        ylabel('Voltage (\muV)')
        title(sprintf('ERP (%s)',subjid));
        drawnow;
        stats = cat(1,stats,array2table({elecname,subjid,latency,minp,sumt},...
            'VariableNames',{'electrode','subjid','latency','minpval','sumtval'}));       
    end
    outdir = fullfile(parentfolder,'results','erp');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-jpg','-painters','-r300','-nofontswap')
end
save(fullfile(outdir,'erpstats.mat'),'stats','SS','FF');
disp(' * * * ');
disp('data saved');

%%    
%     for i = hpchan
%         t=ALLEEG(2).times./1000;
%         indR = find(ALLEEG(2).times >= 0 & ALLEEG(2).times < 185000);
%         indM = find(ALLEEG(2).times >= 0 & ALLEEG(2).times < 370000);
%         R1 = ALLEEG(1).data(i,indR);
%         R2 = ALLEEG(3).data(i,indR);
%         M1 = ALLEEG(2).data(i,indM);
%         M2 = ALLEEG(4).data(i,indM);        
%         tmp = squeeze(ALLEEG(6).data(i,:,:)); 
%         MT=tmp(:)';
%             
%         if length(MT)<size(MTall,2), continue; end
%         R1all = [R1all; R1];
%         R2all = [R2all; R2];
%         M1all = [M1all; M1];
%         M2all = [M2all; M2];
%         MTall = [MTall; MT];        
%     end
% end

% %% PCA:
% close all
% for i = 1: 10
%     [co,sc,~,~,ex]=pca(zscore(MTall,[],2)');
%     tmp = reshape(sc(:,i)',[1,2048,40]);
%     ind1 = find(multiStrFind(clip_order,'ClipA'));
%     ind2 = find(multiStrFind(clip_order,'ClipB'));
%     load(fullfile('D:\ECoG\pink_panther_anonymized','clip_order.mat'));
%     figure('color','w'); hold on;
%     shadedErrorBar(linspace(-1,7,size(tmp,2)),(squeeze(mean(tmp(:,:,ind1),3))'),(squeeze(std(tmp(:,:,ind1),[],3))')./sqrt(length(ind1)),{'color','b'},0.5)
%     shadedErrorBar(linspace(-1,7,size(tmp,2)),(squeeze(mean(tmp(:,:,ind2),3))'),(squeeze(std(tmp(:,:,ind2),[],3))')./sqrt(length(ind2)),{'color','r'},0.5)
%     title(sprintf('PC%d (var explained: %.2f',i,ex(i)));
% end
% 
 
    
%     %% gather all hippocampal channels:
%     [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0); % Define hippocampus channels
%     hpchan = find(contains({EEG.chanlocs.labels},hippocampus_all_channels));
%     smt = ceil(EEG.srate./5); % 200 ms smoothing
%     cm = [];
%     cr = [];
%     counter = 1;
%     for i = hpchan
%         t=ALLEEG(2).times./1000;
%         indR = find(ALLEEG(2).times >= 0 & ALLEEG(2).times < 185000);
%         indM = find(ALLEEG(2).times >= 0 & ALLEEG(2).times < 370000);
%         R1 = ALLEEG(1).data(i,indR);
%         R2 = ALLEEG(3).data(i,indR);
%         M1 = ALLEEG(2).data(i,indM);
%         M2 = ALLEEG(4).data(i,indM);
% %         R1 = tsmovavg(ALLEEG(1).data(i,indR),'t',smt,2);
% %         R2 = tsmovavg(ALLEEG(3).data(i,indR),'t',smt,2);
% %         M1 = tsmovavg(ALLEEG(2).data(i,indM),'t',smt,2);
% %         M2 = tsmovavg(ALLEEG(4).data(i,indM),'t',smt,2);
%         [cr(counter,1),~]=corr(R1',R2','rows','pairwise','type','spearman');
%         [cm(counter,1),~]=corr(M1',M2','rows','pairwise','type','spearman');    
%         counter = counter+1;
%     e params = stuct;
%     Rall = [Rall; cr];
%     Mall = [Mall; cm]; 
% %     R1all = [R1all, R1];
% %     R2all = [R2all, R2];
% %     M1all = [M1all, M1];
% %     M2all = [M2all, M2];
%  %% Spectrograms parameters:
%     close all
%     addpath(genpath('D:\MATLAB_ToolBoxes\chronux_2_12'))
%     
%     params=struct;
%     params.Fs=500; % Sampling Rate [Hz]
%     % Tapers Configuration:
%     T=5; % Window Size [Sec]
%     W=1; % Half Band Width (Frequency Resolution) [Hz]
%     params.tapers=[(T * W), 7]; % Must be less than (2*T*W)-1
%     params.fpass=[0 250]; % [fmin fmax] [HZ]
%     params.err=0;%[2 0.05]; % Jackknife error bars [p=0.05]
%     params.trialave=0; % Avrage over trials
%     params.pad=1;
%     
%    
%     % RUN ANALYSIS AND PLOT:
%     close all;
%     saveflag=0;
%    
%     EEG=ALLEEG(1);
%     for cnum=hpchan
%         
%         channel=EEG.chanlocs(cnum).labels;
%         counter=counter+1;    
%        
%      
%         % RESTING STATE:
%         data=squeeze(ALLEEG(1).data(cnum,:,:)); 
%         data=bsxfun(@minus,data,repmat(nanmean(data),[size(data,1),1])); % demean    
%         [S0_run1,f_out] = mtspectrumc(data, params ); % Rest run1
%         % RESTING STATE:
%         data=squeeze(ALLEEG(3).data(cnum,:,:));
%         data=bsxfun(@minus,data,repmat(nanmean(data),[size(data,1),1])); % demean 
%         [S0_run2] = mtspectrumc( data, params ); % Rest run2
%         
%         S0=cat(2,S0_run1,S0_run2);
%         
%         % Movie:
%         data=squeeze(ALLEEG(2).data(cnum,:,:));
%         data=bsxfun(@minus,data,repmat(nanmean(data),[size(data,1),1])); % demean 
%         [S1_run1] = mtspectrumc( data, params )';
%                 
%         data=squeeze(ALLEEG(4).data(cnum,:,:));
%         data=bsxfun(@minus,data,repmat(nanmean(data),[size(data,1),1])); % demean 
%         [S1_run2] = mtspectrumc( data, params )';
%         S1=cat(2,S1_run1,S1_run2);
%         
%         % MT:
%         data=squeeze(ALLEEG(6).data(cnum,:,:));
%         data=bsxfun(@minus,data,repmat(nanmean(data),[size(data,1),1])); % demean 
%         [S2] = mtspectrumc( data, params )';
%           
%                  
%         
%         % Normalizes each epoch by the baseline of the corresponding run:        
%         fprintf('\n *** Subject: %s, Channel: %s ***',subjid,channel)
%         S0_norm=cat(2,bsxfun(@rdivide,S0_run1,nanmean(S0_run1,2)),bsxfun(@rdivide,S0_run2,nanmean(S0_run2,2)));
%         S1_norm=cat(2,bsxfun(@rdivide,S1_run1,nanmean(S0_run1,2)),bsxfun(@rdivide,S1_run2,nanmean(S0_run2,2)));
%         S2_norm=bsxfun(@rdivide,S2,nanmean(S0,2));      
% % PLOT:
%          figure; hold on;
%         plot(f_out,10.*log10(nanmean(S0,2)),'b','linewidth',2)
%         plot(f_out,10.*log10(nanmean(S1,2)),'r','linewidth',2)
%         plot(f_out,10.*log10(nanmean(S2,2)),'m','linewidth',2)
%           xlim([0 max(f_out)]) 
%                     
%         figure; hold on;
%         plot(f_out,10.*log10(nanmean(S0_norm,2)),'b','linewidth',2)
%         plot(f_out,10.*log10(nanmean(S1_norm,2)),'r','linewidth',2)
%         plot(f_out,10.*log10(nanmean(S2_norm,2)),'m','linewidth',2)
% 
%         
%         xlim([0 max(f_out)])        
%     end
%        
% 
% end
% 
% 
% %% Correlation coeef. dist.:
% 
% set_figure_colors;
% figure('color','w','name','cross-movie-correlation EEG','position',[0 0 200 200]); hold on;
% 
% histogram(Rall,12,'facecolor',COLOR.gray);
% histogram(Mall,12,'facecolor',COLOR.orange);
% 
% tmpx = get(gca,'xlim');
% tmpy = get(gca,'ylim');
% h1 =  scatter([mean(Mall) mean(Mall)],[tmpy(2) tmpy(2)],30,'v','markerfacecolor',COLOR.orange,'markeredgecolor','none','linewidth',2);
% h2 =  scatter([mean(Rall) mean(Rall)],[tmpy(2) tmpy(2)],30,'v','markerfacecolor',COLOR.black,'markeredgecolor','none','linewidth',2);
% set(gca,'ylim',tmpy*1.2)
% tmpx = get(gca,'xlim');
% tmpy = get(gca,'ylim');
% text(tmpx(1)+0.05*range(tmpx),tmpy(2),sprintf('n=%d, P<10^{%d}',length(Mall),fix(log10(signrank(Mall,Rall)))),'HorizontalAlignment','left','fontsize',8)
% xlabel('Correlation coef.')
% ylabel('Electrode count')
% set_font_size_and_type;
% export_fig(fullfile(groupfigdir,get(gcf,'name')),'-nocrop','-pdf','-painters','-nofontswap','-transparent');
% export_fig(fullfile(groupfigdir,get(gcf,'name')),'-nocrop','-jpg','-painters','-r300','-nofontswap')
% 
% %% PLOT spectrograms:
% % ADJUST OUTDIR:
% outdir=fullfile(figdir,'clips');
% if ~exist(outdir,'dir')
%     mkdir(outdir);
%     disp('Creating Output Directory...')
% end
% 
% EEG = ALLEEG(6); 
% nepochs=length(EEG.epoch);
% A_epochs=[];
% B_epochs=[];
% for i=1:nepochs
%     event=EEG.epoch(i).eventtype;
%     if any(contains(event,'A'))
%         A_epochs=[A_epochs i];
%     elseif any(contains(event,'B'))
%         B_epochs=[B_epochs i];   
%     end
% end
%    
% % Spectograms parameters:
% params = struct;
% if EEG.srate>=500
%     params.freqs=[3 250];
%     params.nfreqs=240;
% else
%     params.freqs=[3 125];
%     params.nfreqs=120;
% end
% params.winsize=64;
% params.cycles=[3 20];
% params.units='abs';
% params.freqscale='log';
% params.timesout=1200;
% epochlim=[EEG.times(1) EEG.times(end)];
% [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,2); % Define hippocampus channels
% hpchan = find(contains({EEG.chanlocs.labels},hippocampus_all_channels));
% clim=[-3 3];
% 
% for cnum = hpchan
%     
%     channel=EEG.chanlocs(cnum).labels;
%      
%     [~,~,~,t_out,f_out,~,~,S] = newtimef(EEG.data(cnum,:,:),EEG.pnts,epochlim, EEG.srate, 'cycles', params.cycles, 'freqs', params.freqs, 'winsize', params.winsize,...
%             'scale',params.units,'plotitc','off','plotersp','off','timesout',params.timesout,'nfreqs',params.nfreqs,'baseline',nan,'freqscale',params.freqscale);
%         
%     S=S.*conj(S); % power
%   
%     % Normalize:
%     BL_time=t_out<0 & t_out>-1000;
%     %S_norm=bsxfun(@rdivide,bsxfun(@minus,S,nanmean(nanmean(S(:,BL_time,:),2),3)),sqrt(nanmean(nanvar(S(:,BL_time,:),[],2),3)));  
%     S_norm=bsxfun(@rdivide,bsxfun(@minus,S,nanmean(S(:,BL_time,:),2)),nanstd(S(:,BL_time,:),[],2));  
%     %%   
%     H=figure('Name',['Clip Response'],'position',[0 0 800 1200],'color','w');
%     set(0,'DefaultAxesFontName', 'Arial');    
%  
%     l='Familiar Video';
%     subplot(4,1,1); hold on;
%     imagesc (t_out, f_out, nanmean(S_norm(:,:,A_epochs),3));
%     plot([0 0], get(gca,'Ylim'), '-k','LineWidth', 2,'LineSmoothing','on'); 
%     plot([4000 4000], get(gca,'Ylim'), '--k','LineWidth', 1,'LineSmoothing','on');  
%     axis xy; axis ([min(t_out) max(t_out) min(f_out) max(f_out)]); 
%     caxis (clim); 
%     title(sprintf('%s (N=%d pictures)',l,numel(A_epochs)),'FontWeight','bold','FontSize',8);
%     xlabel ('Time [ms]', 'fontsize', 6);   ylabel ('Freq [Hz]', 'fontsize', 6);
%     cbar; title ('zscore', 'fontsize', 6);
%     
%     l='Novel Video';
%     subplot(4,1,2); hold on;
%     imagesc (t_out, f_out, nanmean(S_norm(:,:,B_epochs),3));
%     plot([0 0], get(gca,'Ylim'), '-k','LineWidth', 2,'LineSmoothing','on'); 
%     plot([4000 4000], get(gca,'Ylim'), '--k','LineWidth', 1,'LineSmoothing','on');  
%     axis xy; axis ([min(t_out) max(t_out) min(f_out) max(f_out)]); 
%     caxis (clim); 
%     title(sprintf('%s (N=%d pictures)',l,numel(B_epochs)),'FontWeight','bold','FontSize',8);
%     xlabel ('Time [ms]', 'fontsize', 6);   ylabel ('Freq [Hz]', 'fontsize', 6);
%     cbar; title ('zscore', 'fontsize', 6);   
%     
%        
%     f_idx=f_out>50&f_out<150;
%     BLPA=squeeze(nanmean(S_norm(f_idx,:,A_epochs)));
%     BLPB=squeeze(nanmean(S_norm(f_idx,:,B_epochs)));
%     
%     tmp=[BLPA(:); BLPB(:)];
%     [AVG,STD]=robustMean(tmp,1,5,0);     
%     outlier_th=AVG+10*STD;
%     BLPA(BLPA>outlier_th)=nan;
%     BLPB(BLPB>outlier_th)=nan;
%     
%     for k=1:size(BLPA,2)
%         BLPA(:,k)=nanfastsmooth(BLPA(:,k),25,2,0);
%     end
%     for k=1:size(BLPB,2)
%         BLPB(:,k)=nanfastsmooth(BLPB(:,k),25,2,0);
%     end
%     subplot(4,1,3); hold on;    
% 
%    % h1=shadedErrorBar(t_out,mean(BLPA,2),nanstd(BLPA,[],2)./sqrt(size(BLPA,2)),{'k','linesmoothing','on','linewidth',1},0);    
%    % h2=shadedErrorBar(t_out,mean(BLPB,2),nanstd(BLPB,[],2)./sqrt(size(BLPB,2)),{'r','linesmoothing','on','linewidth',1},0);
%     
%     h1=shadedErrorBar(t_out,nanmean(BLPA,2),nanstd(BLPA,[],2).*0,{'k','linesmoothing','on','linewidth',2},0);    
%     h2=shadedErrorBar(t_out,nanmean(BLPB,2),nanstd(BLPB,[],2).*0,{'r','linesmoothing','on','linewidth',2},0);
%  
%     plot([0 0], get(gca,'Ylim'), '-k','LineWidth', 1,'LineSmoothing','on');  
%     plot([4000 4000], get(gca,'Ylim'), '--k','LineWidth', 1,'LineSmoothing','on');  
%     xlabel('Time (ms)')
%     ylabel('Power (zscore)')
%     xlim([min(t_out) max(t_out)])
%     title(sprintf('HFB'),'FontWeight','bold','FontSize',8);
%    
%     h_legend=legend([h1.mainLine h2.mainLine],{'Familiar','Novel'});
%     set(h_legend,'fontsize',8,'Location','NorthEast','LineWidth',2,'Box','off'); 
%     
%     subplot(4,1,4); hold on;    
%     t_idx=EEG.times>min(t_out)&EEG.times<max(t_out); 
%     plot(EEG.times(t_idx),mean(EEG.data(cnum,t_idx,A_epochs),3),'k')
%     plot(EEG.times(t_idx),mean(EEG.data(cnum,t_idx,B_epochs),3),'r')
%     plot([0 0], get(gca,'Ylim'), '-k','LineWidth', 1,'LineSmoothing','on');  
%     plot([4000 4000], get(gca,'Ylim'), '--k','LineWidth', 1,'LineSmoothing','on');  
%     xlabel('Time (ms)')
%     ylabel('Voltage (/mV)')
%     xlim([min(t_out) max(t_out)])
%     title(sprintf('ERP'),'FontWeight','bold','FontSize',8);
%     
%     
%     suptitle(sprintf('Video-clip Response (%s)',channel));
%     figname = sprintf('%s elec %s(%s) clip response',subjid,sprintf('%.3d', cnum),channel); 
% 
%     
%     
%     export_fig(fullfile(outdir,figname),'-jpg','-r150','-painters')      
% %     hgexport(gcf, fullfile(outdir,figname), hgexport('factorystyle'), 'Format', 'png')
%     disp(['Channel # ' num2str(cnum)])
% %    close;
% end
% 
% 
% 
% %% Correlation:
% smt=50;
% % X=smoothts(ALLEEG(2).data(good_channels,:),'b',smt);
% % Y=smoothts(ALLEEG(4).data(good_channels,:),'b',smt);
% close all
% Xall=[];
% Yall=[];
% for i = hpchan
% t=ALLEEG(2).times./1000;
% ind = find(ALLEEG(2).times >= 0 & ALLEEG(2).times < 365000);
% X = tsmovavg(ALLEEG(2).data(i,ind),'t',smt,2);
% Y = tsmovavg(ALLEEG(4).data(i,ind),'t',smt,2);
% [r,p]=corr(X',Y','rows','pairwise','type','pearson'); 
% figure; hold on; plot(t(ind),X,'r');  plot(t(ind),Y,'b');  %plot(t,((D-nanmean(D))/nanstd(D)),'m')
% title(sprintf('Channel: %s, r=%.2f',EEG.chanlocs(i).labels,r))
% Xall=[Xall; X];
% Yall=[Yall; Y];
% end
% C=corr(Xall',Yall','rows','pairwise','type','pearson');
% figure; imagesc(C); caxis([-0.1 0.1])
% 
% 
% %%
% smt=250;
% X=[]; Y=[];
% for i=1:numel(good_channels)
%     ch=good_channels(i);
%     X(i,:)=nanfastsmooth(ALLEEG(2).data(ch,:),smt,3,0.5);
%     Y(i,:)=nanfastsmooth(ALLEEG(4).data(ch,:),smt,3,0.5);
% end
% %X=circshift(fliplr(X),20000,2);
% N=min(size(X,2),size(Y,2));
% COR=corr(X(:,2500:N-2500)',Y(:,2500:N-2500)');
% 
% figure('name','Cross-movie HFB Correlation','color','w','position',[0 500 750 200]);
% 
% 
% ch=1:numel(good_channels);
% 
% hold on; bar(ch,diag(COR),1); axis tight;
% xlabel('Electrode #')
% ylabel('Correlation coefficient (r)');
% ylim([-0.2 0.8])
% 
% E={};
% for i=1:numel(good_channels)
%     tmp=EEG.chanlocs(good_channels(i)).labels;
%     tmp=tmp(1:2);
%     E{i}=tmp;
% end
% [E,idx]=unique(E,'stable');
% y=get(gca,'ylim');
% for k=1:numel(E)
%     text(idx(k)+0.5,y(2),E{k},'fontsize',8);
%     plot([idx(k) idx(k)]-0.5, get(gca,'ylim'),':k')
% end
% 
% %subplot(2,1,2); hold on; imagesc(COR); axis square tight; caxis([-0.5 0.5])
% % save figure:
% set(findall(gcf,'-property','FontSize'),'FontSize',8)
% set(findall(gcf,'-property','FontType'),'FontType','Arial')
% figname=sprintf('%s_%s',subjid,get(gcf,'name')); rgb2cm;
% export_fig(fullfile(figdir,figname),'-nocrop','-pdf','-painters','-transparent')
% export_fig(fullfile(figdir,figname),'-nocrop','-jpg','-opengl','-r150','-nofontswap')
% 
% 
% 
% %%
% close all
% smt=101;
% condition='movie';
%  X=[]; Y=[];
% switch condition
%     case 'rest'
%        
%         for i=1:size(ALLEEG(1).data,1)
%             X(i,:)=nanfastsmooth(ALLEEG(1).data(i,:),smt,2,0.5);
%             Y(i,:)=nanfastsmooth(ALLEEG(3).data(i,:),smt,2,0.5);
%         end
%         N=min(size(X,2),size(Y,2));
%         t=(ALLEEG(1).times(2500:N-2500)/1000)-5;
%         
%     case 'movie'
% 
%         for i=1:size(ALLEEG(2).data,1)
%             X(i,:)=nanfastsmooth(ALLEEG(2).data(i,:),smt,2,0.5);
%             Y(i,:)=nanfastsmooth(ALLEEG(4).data(i,:),smt,2,0.5);
%         end
%         N=min(size(X,2),size(Y,2));
%         t=(ALLEEG(2).times(2500:N-2500)/1000)-5;
% end
% 
% 
% %cnum=find(strcmpi({EEG.chanlocs.labels},'TP1-TP2'));
% 
% for cnum=find(contains({EEG.chanlocs.labels},hippocampus_all_channels))
%     clc;
%     channel=EEG.chanlocs(cnum).labels;
%     TC1=(X(cnum,:));
%     TC2=(Y(cnum,:));
%     
%     TC1=TC1(2500:N-2500);
%     TC2=TC2(2500:N-2500);
%     
%     figure('name',[subjid ' Cross Correlation ' channel],'color','w','position',[0 0 1000 600]);
%     subplot(2,2,1:2); hold on;
%     h1=plot(t,TC1,'r-','linesmoothing','off')
%     h2=plot(t,TC2,'k-','linesmoothing','off')
%     axis tight; xlabel('Time(s)'); ylabel('HFB Amplitude (\muV)')
%     [r,p]=corrcoef(TC1,TC2,'rows','complete');
%     r=r(2);
%     xlim([10 200])
%     title(sprintf('Channel: %s (r=%.3f)',channel,r));
%     legend([h1,h2],{sprintf('%s 1',condition),sprintf('%s 2',condition)}); legend boxoff
%     
%     subplot(2,2,3); hold on;
%     [xcf,lags,~]=crosscorr(TC1,TC2,2500);
%     lags=lags*0.002;
%     plot(lags,xcf,'k','linesmoothing','on','linewidth',1); hold on;
%     [~,m]=max(xcf);
%     strtemp=sprintf('%d ms',lags(m)*1000);
%     text(double(lags(m)+0.5),double(xcf(m)),strtemp,'units','data')
%     stem(lags(m),xcf(m),'r')    
%     xlabel('Lag (s)'); ylabel('Correlation (r)')
%     
%         
%     subplot(2,2,4); hold on;
%     p=1:1:99;
%     prc1=prctile(TC1,p);
%     prc2=prctile(TC2,p);
%     scatter(prc1,prc2,'k+'); hold on;  
%     axis tight
%     axis([min([prc1 prc2])-0.1 max([prc1 prc2])+0.1 min([prc1 prc2])-0.1 max([prc1 prc2])+0.1])
%     h_ls=lsline;
%     set(h_ls,'color',[0 0 1])
%     axis square 
%     plot(get(gca,'xlim'),get(gca,'ylim'),':')
%     xlabel('Movie1'); ylabel('Movie2')
%     %set(gca,'xtick',0:0.2:10,'ytick',0:0.2:10)
%     
%     
%     
%     
%     % save figure:
%     path=fullfile(figdir,'cross-corr');
%     if ~exist(path,'dir')
%         mkdir(path);
%         disp('Creating Output Directory...')
%     end
%     set(findall(gcf,'-property','FontSize'),'FontSize',8)
%     set(findall(gcf,'-property','FontType'),'FontType','Arial')
%     figname=sprintf('%s_%s_rest',subjid,get(gcf,'name')); rgb2cm;
%      %export_fig(fullfile(path,figname),'-nocrop','-pdf','-painters','-transparent')
%    % export_fig(fullfile(path,figname),'-nocrop','-jpg','-opengl','-r150','-nofontswap')
% end
% 
% %% Peak triggered average:
% % Follow-up analysis:
% for cnum=find(multiStrFind({EEG.chanlocs.labels},hippocampus))
%     clc;
%     channel=EEG.chanlocs(cnum).labels;
% %     TC1=10*log10(double(bsxfun(@rdivide,X(cnum,:),nanmedian(X(cnum,:)))));
% %     TC2=10*log10(double(bsxfun(@rdivide,Y(cnum,:),nanmedian(Y(cnum,:)))));
%     TC1=double(X(cnum,:));
%     TC2=double(Y(cnum,:));
%     
%     TC1=TC1(2500:N-2500);
%     TC2=TC2(2500:N-2500);
% 
%     figure('name',[subjid ' Cross Correlation ' channel],'color','w','position',[0 0 1000 600]);
%     subplot(2,2,1:2); hold on;
%     h1=plot(t,TC1,'r-','linesmoothing','off');
%     h2=plot(t,TC2,'k-','linesmoothing','off');
%     [pks1,locs1]=findpeaks(zscore(TC1),'MINPEAKHEIGHT',1.5,'MINPEAKDISTANCE',500);
%     [pks2,locs2]=findpeaks(zscore(TC2),'MINPEAKHEIGHT',1.5,'MINPEAKDISTANCE',500);
%     scatter(t(locs1),TC1(locs1),'bo')
%     scatter(t(locs2),TC2(locs2),'go')
%     axis tight
%     title(sprintf('Channel: %s (r=%.3f)',channel,r));
%     legend([h1,h2],{sprintf('%s 1',condition),sprintf('%s 2',condition)}); legend boxoff
%     
%     epochlim=[-2 2];
%     [E1,newtime1]=epoch(TC1,t(locs1),epochlim,'srate',500);
%     [E2,newtime2]=epoch(TC2,t(locs1),epochlim,'srate',500);
%     subplot(2,2,3); hold on;
%     tepoch=newtime1(1):0.002:newtime1(2);
%     shadedErrorBar(tepoch,squeeze(nanmean(E1,3)),squeeze(nanstd(E1,[],3)./sqrt(size(E1,3))),'b',0.5)
%     shadedErrorBar(tepoch,squeeze(nanmean(E2,3)),squeeze(nanstd(E2,[],3)./sqrt(size(E2,3))),'k',0.5)
%  
%     
%     % save figure:
%     path=fullfile(figdir,'cross-corr');
%     if ~exist(path,'dir')
%         mkdir(path);
%         disp('Creating Output Directory...')
%     end
%     set(findall(gcf,'-property','FontSize'),'FontSize',8)
%     set(findall(gcf,'-property','FontType'),'FontType','Arial')
%     figname=sprintf('%s_%s_rest',subjid,get(gcf,'name')); rgb2cm;
%      %export_fig(fullfile(path,figname),'-nocrop','-pdf','-painters','-transparent')
%    % export_fig(fullfile(path,figname),'-nocrop','-jpg','-opengl','-r150','-nofontswap')
% end
% 
% 
% 
% 
% 
% %%
% close all
% 
% %cnum=find(strcmpi({EEG.chanlocs.labels},'HA8'));
% EEG=ALLEEG(6);
% %[EEG, ~,b] = pop_firws(EEG, 'fcutoff', 10, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', 826, 'minphase', 0);
% smt=101;
% for cnum=find(multiStrFind({EEG.chanlocs.labels},'HP8'))
%     channel=EEG.chanlocs(cnum).labels;
%     A=smoothts(squeeze(EEG.data(cnum,:,multiStrFind({EEG.epoch.eventtype},'A')))','b',smt)';
%     B=smoothts(squeeze(EEG.data(cnum,:,multiStrFind({EEG.epoch.eventtype},'B')))','b',smt)';
%     
%     BLavg=mean(mean(EEG.data(cnum,EEG.times>-1000 & EEG.times<0,:),2),3);
%     BLstd=sqrt(mean(var(EEG.data(cnum,EEG.times>-1000 & EEG.times<0,:),[],2),3));
%     % zscore:
%     A=bsxfun(@rdivide,bsxfun(@minus,A,BLavg),BLstd);
%     B=bsxfun(@rdivide,bsxfun(@minus,B,BLavg),BLstd);
%     
%     
%     figure('color','w','position',[0 0 750 750],'name',[subjid ' Recogntition Test Channel ' channel]);
%     
%     subplot(3,1,1); hold on;
%     pos=get(gca,'position');
%     set(gca,'position',[pos(1),pos(2),pos(3)*0.91,pos(4)])
%     
% %     [avgA,stdA]=robustMean(A,2,3,0);
% %     [avgB,stdB]=robustMean(B,2,3,0);
%     avgA=mean(A,2); avgB=mean(B,2);
%     stdA=std(A,[],2); stdB=std(B,[],2);
%     h1=shadedErrorBar(EEG.times/1000,avgA',stdA./sqrt(size(A,2))',{'color','r','linewidth',1,'linesmoothing','on'},1);
%     h2=shadedErrorBar(EEG.times/1000,avgB',stdB./sqrt(size(B,2))',{'color','k','linewidth',1,'linesmoothing','on'},1);
%     
%     xlabel('Time from clip onset(s)'); ylabel('HFB Amplitude (zscore)'); xlim([-1 6]) 
%     plot([0 0],get(gca,'ylim'),'k-','linewidth',2)
%     plot([4 4],get(gca,'ylim'),'k-','linewidth',1)
%     legend([h1.mainLine; h2.mainLine],{'same movie','different movie'}); legend boxoff
%     title(sprintf('Channel %s',channel),'fontsize',14)
%     
%     subplot(3,1,2); hold on;
% 
%     imagesc(EEG.times/1000,1:20,[A']);
%     ylabel('Trial #');
%     axis tight; caxis([-5 5]); xlim([-1 6])
%     plot([0 0],get(gca,'ylim'),'k-','linewidth',2)
%     plot([4 4],get(gca,'ylim'),'k-','linewidth',1)
%     title('Clips - SAME movie');
%     hc=colorbar('location','eastoutside');
%     
%     subplot(3,1,3); hold on;
%     pos=get(gca,'position');
%     set(gca,'position',[pos(1),pos(2)+0.2*pos(2),pos(3),pos(4)])
%     imagesc(EEG.times/1000,1:20,[B']);
%     ylabel('Trial #'); xlabel('Time from clip onset(s)');
%     axis tight; caxis([-5 5]); xlim([-1 6])
%     plot([0 0],get(gca,'ylim'),'k-','linewidth',2)
%     plot([4 4],get(gca,'ylim'),'k-','linewidth',1)    
%     title('Clips - DIFFERENT movie');
%     hc=colorbar('location','eastoutside'); 
%     
%     % save figure:
%     path=fullfile(figdir,'Recognition test');
%     if ~exist(path,'dir')
%         mkdir(path);
%         disp('Creating Output Directory...')
%     end
%     set(findall(gcf,'-property','FontSize'),'FontSize',8)
%     set(findall(gcf,'-property','FontType'),'FontType','Arial')
%     figname=get(gcf,'name'); rgb2cm;
% %     export_fig(fullfile(path,figname),'-nocrop','-pdf','-painters','-transparent')
% %     export_fig(fullfile(path,figname),'-nocrop','-jpg','-painters','-r150','-nofontswap')
% %     close;
% end
% 
% 
% %% PLOT:
% % Get indices of A and B epochs
% close all
% EEG=ALLEEG(6);
% 
% outdir = fullfile(figdir,'spectrograms');
% if ~exist(outdir,'dir')
%     mkdir(outdir);
%     disp('Creating Output Directory...')
% end
% 
% %==========================================================================
% % figure;
% % [PP3,~,~,t,f_out] = newtimef({EEG.data(cnum,:,A_epochs), EEG.data(cnum,:,B_epochs)},size(EEG.data,2), [EEG.times(1) EEG.times(end)], EEG.srate, 'cycles', cycles, 'freqs', f, 'winsize', winsize,'baseline',[-1000 0],...
% %  'scale',p_units,'plotitc','on','plotersp','on','timesout',timesout,'trialbase','full','padratio',padratio,'alpha',0.01,'pcontour','on');
% %==========================================================================
% 
% 
% % Spectogram Parameters:
% f=[3 200];
% winsize=64;
% timesout=1000;
% padratio=4;
% p_units='log';
% cycles=[1 20];
% ersplim=[-10 10];
% categories={'A','B'};
% 
% for cnum=find(multiStrFind({EEG.chanlocs.labels},'TP1')) %1:EEG.nbchan
%     
%     channel=EEG.chanlocs(cnum).labels;
%     H=figure('Name',['Spectrum'],'position',[0 0 1200 800],'color','w');
%     set(0,'DefaultAxesFontName', 'Arial');
% 
%     plotidx=0;
%     [PP1,~,~,t,f_out] = newtimef(EEG.data(cnum,:,A_epochs),size(EEG.data,2), [EEG.times(1) EEG.times(end)], EEG.srate, 'cycles', cycles, 'freqs', f, 'winsize', winsize,'baseline',[-1500 0],...
%                                   'scale',p_units,'plotitc','off','plotersp','off','timesout',timesout,'trialbase','full','padratio',padratio);
%     [PP2,~,~,t,f_out] = newtimef(EEG.data(cnum,:,B_epochs),size(EEG.data,2), [EEG.times(1) EEG.times(end)], EEG.srate, 'cycles', cycles, 'freqs', f, 'winsize', winsize,'baseline',[-1500 0],...
%                          'scale',p_units,'plotitc','off','plotersp','off','timesout',timesout,'trialbase','full','padratio',padratio);
%                      
%                      
%     subplot(3,1,1); hold on;
%     imagesc (t, f_out, PP1);
%     plot([0 0], get(gca,'Ylim'), '-k','LineWidth', 2,'LineSmoothing','on'); 
%     plot([4000 4000], get(gca,'Ylim'), '--k','LineWidth', 1,'LineSmoothing','on');  
%     axis xy; axis ([min(t) max(t) min(f_out) max(f_out)]);   
%     caxis (ersplim); 
%     title(sprintf('Same movie (N=%d clips)',numel(A_epochs)),'FontWeight','bold','FontSize',9);
%     xlabel ('Time [ms]', 'fontsize', 7);   ylabel ('Freq [Hz]', 'fontsize', 7);
%     cbar; title ('Db', 'fontsize', 8);
%    
%     subplot(3,1,2); hold on;
%     imagesc (t, f_out, PP2);
%     plot([0 0], get(gca,'Ylim'), '-k','LineWidth', 2,'LineSmoothing','on'); 
%     plot([4000 4000], get(gca,'Ylim'), '--k','LineWidth', 1,'LineSmoothing','on');  
%     axis xy; axis ([min(t) max(t) min(f_out) max(f_out)]); 
%     caxis (ersplim); 
%     title(sprintf('Different movie (N=%d clips)',numel(B_epochs)),'FontWeight','bold','FontSize',9);
%     xlabel ('Time [ms]', 'fontsize', 7);   ylabel ('Freq [Hz]', 'fontsize', 7);
%     cbar; title ('Db', 'fontsize', 8);
%    
%     
%     
%     t_idx=EEG.times>min(t)&EEG.times<max(t);
%     
%     subplot(3,1,3); hold on;
%     plot(EEG.times(t_idx),mean(EEG.data(cnum,t_idx,A_epochs),3),'-k','LineWidth', 1,'LineSmoothing','off')
%     plot(EEG.times(t_idx),mean(EEG.data(cnum,t_idx,B_epochs),3),'-r','LineWidth', 1,'LineSmoothing','off')
%     plot([0 0], get(gca,'Ylim'), '-k','LineWidth', 1,'LineSmoothing','on');  
%     plot([4000 4000], get(gca,'Ylim'), '--k','LineWidth', 1,'LineSmoothing','on');  
%     xlabel('Time (ms)')
%     ylabel('Voltage (\muV)')
%     xlim([min(t) max(t)])
%     title(sprintf('Grand ERP (N=%d Clips)  -  %s',numel(A_epochs),channel),'FontWeight','bold','FontSize',9);
%     suptitle(sprintf('Video-clip Response (%s)',channel));
%     figname = sprintf('%s elec %s(%s) clip response',subjid,sprintf('%.3d', cnum),channel); 
%     
%     %export_fig(fullfile(outdir,figname),'-jpg','-r150','-painters')      
% 
%     disp(['Channel # ' num2str(cnum)])
%   
% end
% 
%     