close all;
clear all;
clc; warning('off');

% set the relevant path:
path_to_toolboxes = 'D:\MATLAB_ToolBoxes\';
path_to_project = 'D:\ECoG\pink_panther_anonymized\';
addpath(fullfile(path_to_toolboxes,'eeglab14_1_2b'));
addpath(genpath(fullfile(path_to_project,'matlab_scripts')));
rmpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));
rmpath(genpath(fullfile('D:\ECoG\MMR','matlab_scripts')));

[ALLEEG, EEG, CURRENTSET, ~] = eeglab;

subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18','SUB01ASSUTA'};

for subjid = {'SUB01ASSUTA'} %subjects(7:end-3)
    subjid = cell2mat(subjid);
    close all;
    clearex('subjects', 'subjid', 'ALLEEG', 'EEG', 'path_to_project','path_to_toolboxes');
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    
    maindir=fullfile(path_to_project,subjid);
    outdir=fullfile(maindir,'EEGLAB_datasets_BP','preprocessed'); % ADJUST OUTDIR
    inFileName = fullfile(maindir,'EEGLAB_datasets','raw',[subjid '_pink_panther_entire_run.set']);
    mkdir(outdir);
   
    outFileName =[subjid '_pink_panther_preprocessed_BP_montage.set'];
    
    % Select XLS electrodes file:
    path = fullfile(maindir,'BioImage');
    [channel_location_file]=dir(fullfile(path,'*TAL*'));
    channel_location_file = fullfile(path,channel_location_file.name);
    
    % load individual brain:
    S_brain = struct;
    S_brain.plotsurf = 'pial';
    S_brain.layout = 'compact';
    S_brain.surfacealpha = 1;
    S_brain.meshdir = fullfile(path_to_project,'SUMA_meshData');
    
    % Load raw EEGLAB dataset:
    [EEG] = pop_loadset('filename', inFileName);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
    eeglab redraw
    
    if exist('run_commands','var')
        for k=1:numel(run_commands)
            eval(run_commands{k});
        end
    end
    
   
   %% set AUX Channels:
    
    auxsearch = {'MKR2','ECG','EKG','EMD','EMG'};
    
    for k = auxsearch
        ind = find(contains({EEG.chanlocs.labels},k,'IgnoreCase',true));
        if isempty(ind), continue; 
        else
            for i = 1:numel(ind)
                if i==1
                    EEG=pop_chanedit(EEG,'changefield',{ind(i) 'labels' cell2mat(k)});   
                elseif i>1
                    EEG=pop_chanedit(EEG,'changefield',{ind(i) 'labels' [cell2mat(k) '-' num2str(i)]});                   
                end
                fprintf('\n Renaming channel %d: %s \n',ind(i),EEG.chanlocs(ind(i)).labels);
            end
        end
    end
        
    if sum(multiStrFind({EEG.chanlocs.labels},'TRIG'))~=1
        fprintf('\nChannels: %d\n', find(contains({EEG.chanlocs.labels},'TRIG')))       
        error('******** There is a problem with the TRIG channel! ********');
    end
    
    aux_channels=[find(contains({EEG.chanlocs.labels},'EOG','IgnoreCase',true)),...      
                  find(contains({EEG.chanlocs.labels},'ECG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'EMG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'TRIG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'EMD','IgnoreCase',true))];
              

    %% Remove non-channels:
    if contains(subjid,'ASSUTA')
        load(fullfile(maindir,'BioImage','SUB01ASSUTA_matlab.mat'),'bio_order');
        ch_to_remove=setdiff(find(~contains({EEG.chanlocs.labels},bio_order)),aux_channels);
    else
        ch_to_remove=setdiff(find(~multiStrFind({EEG.chanlocs.labels},'EEG')),aux_channels);
    end
    
    for ch={EEG.chanlocs(ch_to_remove).labels}
        fprintf('--- Removing Channel: %s --- \n',cell2mat(ch));
    end
    
    EEG = pop_select(EEG,'nochannel',ch_to_remove);
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;
    
    % Rename EEG channels:
    idx=find(multiStrFind({EEG.chanlocs.labels},'EEG'));
    for ch=idx
        new_label=EEG.chanlocs(ch).labels(5:end);
        fprintf('--- Relabeling Channel: %s --- \n',new_label);
        EEG=pop_chanedit(EEG,'changefield',{ch 'labels' new_label});
    end
    aux_channels=[find(contains({EEG.chanlocs.labels},'EOG','IgnoreCase',true)),...      
                  find(contains({EEG.chanlocs.labels},'ECG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'EMG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'TRIG','IgnoreCase',true)),...
                  find(contains({EEG.chanlocs.labels},'EMD','IgnoreCase',true))];
    
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);
    
    % remove extra time after the end of the experiment:
    Define_experiment_end_PP;
    EEG = pop_select(EEG,'nopoint',experiment_end:EEG.pnts);
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);
    
    % Resample to 500 Hz: (optional; try to avoid)
    % EEG = pop_resample(EEG,500);
    % EEG.setname=[outFileName(1:end-4) ' - resampled'];
    % [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    
    eeglab redraw;
    
    %% Anatomical location (native space):
    load(fullfile(maindir,[subjid '_channels_list'])); % Load good_channels list
   
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,1,0); % Define hippocampus channels
    FSflag = 0;
    
    % Check if there is a FreeSurfer reconstructed brain:     
    elocDir = fullfile(S_brain.meshdir,subjid);
    if exist(fullfile(elocDir,'electrodes.mat'),'file'), FSflag = 1; end

    % Load electrode locations: (MNI)
    EEG = ReadElectrodeCoord(EEG,channel_location_file,maindir);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;
    
    % Plot FS Brain (optional):
    if FSflag
        [S_brain,H_brain,SUMAsrf] = plot_FS_brain_master(subjid,S_brain);
        elocDir=fullfile(S_brain.meshdir,subjid);
        load(fullfile(elocDir,'SUMAprojectedElectrodes.mat'))
        ch_labels=[SUMAprojectedElectrodes.elecNames];
        ePlot1=ch_labels;
        ePlot2=ch_labels(multiStrFind(ch_labels,{hippocampus(isstrprop(hippocampus, 'alpha'))}));
        eSize=1.5; % 1.5 radius in mm
        eColor1=[0 0 0];
        eColor2=[1 0 0];
        textFlag=1;
        S_brain=plot_Fs_electrode_master(subjid,H_brain,S_brain,SUMAsrf,elocDir,ePlot1,eColor1,eSize,textFlag,'gouraud','top');
        S_brain=plot_Fs_electrode_master(subjid,H_brain,S_brain,SUMAsrf,elocDir,ePlot2,eColor2,eSize,textFlag,'gouraud','top');
        set(findall(gcf,'type','patch'),'facealpha',0.75)
    end
    
    % Plot MNI Brain to verify:
    figure; hold on;
    mesh=fullfile(path_to_project,'MNI_brain_mesh','MNI_brain_downsampled.mat');
    headplot_itzik(EEG.data,fullfile(maindir,[subjid '_spline_file_MNIbrain.spl']),[],'meshfile',mesh,'electrodes','on', ...
        'title',subjid,'labels',1,'cbar',0, 'maplimits','absmax','colormap',colormap('Gray'));
    alpha(0.15)
      
    %==================================================================
    
    % keep only channels that are connected and are in the brain (verified visually):
    if FSflag
        load(fullfile(elocDir,'electrodes.mat'))
        inBrain_channels = find(ismember({EEG.chanlocs.labels},electrodes.elecNames));
        good_channels = intersect(good_channels, inBrain_channels);
        % get some extra electrode info (optional):
        %         channelFSlabel = {};
        %         numOfvoxels = [];
        %         for k = good_channels
        %             ch_label = EEG.chanlocs(k).labels;
        %             electrodeInd = find(strcmpi(electrodes.elecNames,ch_label));
        %             if isempty(electrodeInd),  channelFSlabel{k}='N/A'; continue; end
        %             currentROI = electrodes.aparcaseg.bestLabel.labels(electrodeInd);
        %             channelFSlabel(k,1) = currentROI;
        %             numOfvoxels(k) = electrodes.aparcaseg.bestLabel.NumOfVoxel(electrodeInd);
        %         end
    else
        inBrain_channels = find([EEG.chanlocs.X] ~= 0 & [EEG.chanlocs.Y] ~= 0 & [EEG.chanlocs.Z] ~= 0);
        good_channels = intersect(good_channels, inBrain_channels);
    end
    
    % Compute common cortical average after 50Hz removal:
    if contains(subjid,'ASSUTA'), notchFreqs=[50:50:450];
    else,  notchFreqs=[50]; end
    filterWidth=1.5; % Hz
    EEG_clean=EEG;
    for f=notchFreqs
        % Adjust the filter order manually! (use the EEGLAB menu to calculate the order)
        [EEG_clean,~,b] = pop_firws(EEG_clean, 'fcutoff', [f-filterWidth f+filterWidth], 'ftype', 'bandstop', 'wtype', 'hamming', 'forder', 1100);
        figure; freqz(b);
    end
    
    %% common average:
    CREF = robustMean(EEG_clean.data(good_channels,:),1,5);
    
    % Compute cortical contacts average for control ripple detection (optional):
    %   cortical_channels = find(multiStrFind(channelFSlabel,'ctx')&~multiStrFind(channelFSlabel,'unknown'));
    %   mean_good_signal=robustMean(EEG.data(cortical_channels,:),1,5);
    % alternative:
    mean_good_signal = CREF;
    
    % Compute EMG-from-high-frequency:  (Schomburg et al, 2014; Watson et al. 2016);
    close all;
    disp('Processing EMG channel...')
    forder_HP=330;
    HF = eegfilt(mean_good_signal,EEG.srate,100,[],0,forder_HP,0,'fir1',0);
    EMGHF = abs(hilbert(HF));
    ind = find(strcmpi({EEG.chanlocs.labels},'EMGHF'));
    if isempty(ind)
        EEG.data(end+1,:) = EMGHF;
        EEG.nbchan = size(EEG.data,1);
        EEG.chanlocs(end+1).labels = 'EMGHF';
    else
        EEG.data(ind,:) = EMGHF;
    end
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = eeg_checkset( EEG );
    aux_channels=union(aux_channels,find(strcmpi({EEG.chanlocs.labels},'EMGHF')));
    
    
    % Save common average as CREF channel:
    EEG.data(end+1,:) = CREF;
    EEG.nbchan = size(EEG.data,1);
    EEG.chanlocs(end+1).labels = 'CREF';
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = eeg_checkset( EEG );
    aux_channels=union(aux_channels,find(strcmpi({EEG.chanlocs.labels},'CREF')));
    clear EEG_clean
    
    %% Bipolar montage:
    %==========================================================================
    
    XYZ = [];
    good_ch_labels={EEG.chanlocs(good_channels).labels};
    % For Freesurfer native space coord:
    if exist('electrodes','var')
        for i=1:numel(good_ch_labels)
            XYZ(i,:) = electrodes.coord.afniXYZ(strcmpi(electrodes.elecNames,good_ch_labels{i}),:);
        end
    else
        % For MNI coords:
         for i=1:numel(good_ch_labels)
            ch_idx = strcmpi({EEG.chanlocs.labels},good_ch_labels{i});
            XYZ(i,:)=[cell2mat({EEG.chanlocs(ch_idx).X}),...
                cell2mat({EEG.chanlocs(ch_idx).Y}),...
                cell2mat({EEG.chanlocs(ch_idx).Z})];
        end
    end
    
      
    counter=1;
    sig=[]; ref=[];
    new_labels={};
    clc;
    % STEP 1: all contacts in the hippocampal depth electrode are paired with a nearby WM contact -
    if isempty(hippocampus_all_channels), fprintf('\n --> No hippocampal channels, moving on... \n');
    else
        
        hippocampal_shank = unique(cellfun(@(x)x(isstrprop(x, 'alpha')),hippocampus_all_channels,'UniformOutput',false));
        if isempty(WM_ref), error('Missing a WM contact!'); end
        
        for k = 1:numel(hippocampal_shank)
            
            % hippocampal contacts:
            hipp_array = good_ch_labels(multiStrFind(good_ch_labels,hippocampal_shank(k))); % get all good contacts on the hippocampal electrode shank
            
            %% Choose the closest white matter reference site:
            tmp1 = find(multiStrFind(good_ch_labels,hippocampal_shank(k)));  % index within good_ch_labels and XYZ
            tmp2 = find(multiStrFind(good_ch_labels,WM_ref));                % index within good_ch_labels and XYZ
            d = [];
            for ii = 1:numel(tmp2)
                d(1,ii) = mean(sqrt(sum((bsxfun(@minus,XYZ(tmp1,:),XYZ(tmp2(ii),:)).^2),2)));
            end
            [~,ind] = min(d);
            
            channel2 = good_ch_labels{tmp2(ind)};                   % channel label
            cnum2 = find(strcmpi({EEG.chanlocs.labels},channel2));  % real index (within the EEGLAB dataset)
            
            hipp_array = setdiff(hipp_array,channel2);
            
            for i = 1:numel(hipp_array)
                channel1 = hipp_array{i};
                cnum1 = find(strcmpi({EEG.chanlocs.labels},channel1)); % index within the EEG dataset
                if str2num(channel1(isstrprop(channel1, 'digit'))) >= 8, continue; end % skip the superficial hippocampus-electrode contacts (contact #8 and above)
                % Usually, Da/Dh/Dp/Hb/Hh contacts 1-5 are located within the
                % hippocampus (since they are the deepest contacts)
                
                % Calculate distance in mm (for sanity check):
                tmp1 = find(strcmpi(good_ch_labels,channel1));
                tmp2 = find(strcmpi(good_ch_labels,channel2));
                d = sqrt(sum((bsxfun(@minus,XYZ(tmp1,:),XYZ(tmp2,:)).^2),2));
                
                % Sanity Check:
                if ~ismember(cnum2,good_channels)
                    error(sprintf('\n --> Wrong Channel: %s \n Please Verify... \n',channel2));
                end
                fprintf('\n *** Channel: %s [%3d] - %s [%3d]  (%.2f mm) *** \n',channel1,cnum1,channel2,cnum2,d);
                new_labels{cnum1}=sprintf('%s-%s',channel1,channel2);
                sig(counter)=cnum1;
                ref(counter)=cnum2;
                counter=counter+1;
            end
        end
    end
    
    % STEP 2: process all remaining channels -
    
    for i=1:numel(good_ch_labels)
        
        channel1=good_ch_labels{i};
        cnum1=find(strcmpi({EEG.chanlocs.labels},channel1)); % index within the EEG dataset
        
        current_array=good_ch_labels(multiStrFind(good_ch_labels,channel1(isstrprop(channel1, 'alpha'))));
        if ismember(cnum1,sig), continue; end
        if ismember(cnum1,ref), continue; end % to use only unique channels
        % Choose:
        current_array=setdiff(current_array,{EEG.chanlocs([sig, ref]).labels}); % only unique pairs
        % current_array=setdiff(current_array,{EEG.chanlocs(sig).labels});        % allow duplicates
        current_array=setdiff(current_array,channel1);
        
        if isempty(current_array)
            fprintf('\n --> Skipping Channel: %s    (last contact in the strip)\n',channel1);
            continue;
        end
        
        % Find the cloest channel on the strip to serve as reference:
        dist=[];
        for k=1:numel(current_array)
            dist(k)=sqrt(sum((bsxfun(@minus,XYZ(i,:),XYZ(strcmpi(good_ch_labels,current_array{k}),:)).^2),2));
        end
        [d,idx]=min(dist);
        channel2=current_array{idx};
        cnum2=find(strcmpi({EEG.chanlocs.labels},channel2));
        
        % exclude pairs that are >20 mm apart from each other
        if  d>20
            fprintf('\n --> Skipping Channel: %s [%3d] - %s [%3d]  (%.2f mm) \n',channel1,cnum1,channel2,cnum2,d);
            continue;
        end
        
        % Sanity Check:
        if ~ismember(cnum2,good_channels)
            error(sprintf('\n --> Wrong Channel: %s \n Please Verify... \n',channel2));
        end
        fprintf('\n *** Channel: %s [%3d] - %s [%3d]  (%.2f mm) *** \n',channel1,cnum1,channel2,cnum2,d);
        new_labels{cnum1}=sprintf('%s-%s',channel1,channel2);
        sig(counter)=cnum1;
        ref(counter)=cnum2;
        counter=counter+1;
    end
    
    figure;
    scatter(sig,ref,30,'.k')
    xlabel('Sig ch.'); ylabel('Ref ch.');
    title('Indices of all electrode pairs')
    
    %% Re-referenceing:
    reref_data=zeros(numel(sig),EEG.pnts);
    for i=1:numel(sig)
        reref_data(i,:)=EEG.data(sig(i),:)-EEG.data(ref(i),:);
        fprintf('\n Subtracting Electrodes: %3d - %3d \n',sig(i),ref(i));
    end
    
    for i=1:EEG.nbchan
        if ismember(i,sig)
            EEG=pop_chanedit(EEG,'changefield',{i 'labels' new_labels{i}});
            EEG=pop_chanedit(EEG,'changefield',{i 'type' 'signal'});
            EEG.data(i,:)=reref_data(sig==i,:);
        end
        if ismember(i,aux_channels)
            EEG=pop_chanedit(EEG,'changefield',{i 'type' 'aux'});
        end
    end
    
    EEG = pop_select(EEG,'channel',[aux_channels sig]);
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw;
    
    %==========================================================================
    % Remove DC from each electrode:
    EEG = pop_rmbase(EEG,[EEG.times(1) EEG.times(end)]);
    EEG.setname=[outFileName(1:end-4) ' - DC removed'];
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw;
    %==========================================================================
    
    % Remove line noise using the new EEGLAB FIR filter (this time from the bipolar derivations):
    good_channels=find(strcmpi({EEG.chanlocs.type},'signal'));
    %figure; spectopo(EEG.data(good_channels,:),0,EEG.srate,'percent',10,'title','Before Removing Line Noise');
    notchFreqs=[50];
    filterWidth=1.5; % Hz
    EEG_clean=EEG;
    for f=notchFreqs
        % Adjust the filter order manually! (use the EEGLAB menu to calculate the order)
        [EEG_clean,~,b] = pop_firws(EEG_clean, 'fcutoff', [f-filterWidth f+filterWidth], 'ftype', 'bandstop', 'wtype', 'hamming', 'forder', 1100);
        figure; freqz(b);
    end
%    figure; spectopo(EEG_clean.data(good_channels,:),0,EEG.srate,'percent',10,'title','After Removing Line Noise');
    
    % Store DATA:
    EEG=EEG_clean;
    EEG.setname=[outFileName(1:end-4) ' - Filtered'];
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw;
    eegplot(EEG.data(good_channels,:),'color','off','srate',EEG.srate,'winlength',15,'limits',[EEG.times(1) EEG.times(end)])
    
    
    %% Save set:
    
    EEG.setname=outFileName(1:end-4);
    EEG = pop_saveset( EEG,  'filename', outFileName, 'filepath', outdir);
    disp('data saved')
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    eegplot(EEG.data(strcmpi({EEG.chanlocs.type},'signal'),:),'color','off','srate',EEG.srate,'winlength',15,'limits',[EEG.times(1) EEG.times(end)])
    
end



