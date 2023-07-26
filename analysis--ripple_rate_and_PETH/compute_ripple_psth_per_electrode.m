clear all
close all
clc;
% set the relevant path:
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;

% add/remove paths:
addpath(fullfile(path_to_toolboxes,'eeglab2021.1'));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));

% load eeglab:
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Parameters to set manually:
ref_flag=2; % 1 = Common Ref; 2 = Bipolar montage;
time_locking_event = 'stimonset'; % time locking event (stimonset/rt)

%==========================================================================
% set path to ripples information:
ripplesdir = fullfile(parentfolder,'results','data','Ripple_times_hamming_2-4std_adjusted_band_20ms_30ms');

% =================================================================
subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};

        
for i = 1:numel(subjects)
    subjid=subjects{i};
    clearex('i','subjects', 'subjid',  'ALLEEG' , 'EEG', 'ref_flag', 'path_to_toolboxes','parentfolder',...
            'ripplesdir','scenario','time_locking_event','EXC')
    maindir=fullfile(parentfolder,subjid);
    datasetID='pink_panther';
    
    switch ref_flag
        case 1
            datadir=fullfile(maindir,'EEGLAB_datasets');
        case 2
            datadir=fullfile(maindir,'EEGLAB_datasets_BP');
    end    
    ALLEEG=[]; EEG=[]; CURRENTSET=1;    
    
    % Load eeglab datasets:
    filenames=dir(fullfile(datadir,sprintf('*%s*.set',datasetID)));
    for SET=1:numel(filenames)
        filename=fullfile(datadir,filenames(SET).name);
        [EEG] = pop_loadset('filename', filename);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
        eeglab redraw;
    end
    
    DATA = struct;   
    % Set data structure:
    FN={'Rall','R1','R2','M1','M2','MT'};
    for j = 1:length(FN)
        DATA.(FN{j}).subjid = [];
        DATA.(FN{j}).taskid = [];
        DATA.(FN{j}).channelid = [];        
        DATA.(FN{j}).stimulus_type = [];
        DATA.(FN{j}).stimulus = [];
        DATA.(FN{j}).RT = [];
        DATA.(FN{j}).response = [];
        DATA.(FN{j}).correct = [];
        DATA.(FN{j}).raster = [];
    end
    DATA.rippletable = [];
    DATA.Fs = 1000;
    eventsInMovie = readtable(fullfile(parentfolder,'events.xlsx'));

    % -----------------------------------------------------------------
    for SET = 1:numel(ALLEEG)
        EEG = ALLEEG(SET);
        task_label = EEG.setname(end-5:end-4); % task label: M1 M2 R1 R2 MT         
        strind = round(EEG.event(strcmpi({EEG.event.type},'str')).latency); % finds the start point
        finind = round(EEG.event(strcmpi({EEG.event.type},'fin')).latency); % finds the end point
        margin_correction = EEG.times(strind)./1000;       
        blockDuration = EEG.times(finind)./1000;
        
        if strcmpi(task_label,'MT')
            stimulus = {EEG.event(multiStrFind({EEG.event.type},{'ClipA','ClipB'})).type}';
            stimulus_type = cellfun(@(x)x(1:5),stimulus,'UniformOutput',0);      
            stimulus_onset = EEG.times(round([EEG.event(multiStrFind({EEG.event.type},{'ClipA','ClipB'})).latency]))./1000;
            
            rt = EEG.times(round([EEG.event(multiStrFind({EEG.event.type},{'old','new'})).latency]))./1000;
            response = {EEG.event(multiStrFind({EEG.event.type},{'old','new'})).type}';
            valid_responses = {}; % find valid response
            valid_rt = []; 
            index = [];
            for k = 1:length(stimulus_onset)
                ii = find( (rt-stimulus_onset(k) < 8) & (rt-stimulus_onset(k) > 0), 1, 'last');  
                if isempty(ii), valid_responses{k,1} = 'n/a';  valid_rt(k,1) = stimulus_onset(k)+4; % if there was no response, take the trial offset
                else, valid_rt(k,1) = rt(ii); valid_responses{k,1} = response{ii};  end
            end
            
            RT = valid_rt - margin_correction;
            stimulus_onset = stimulus_onset' - margin_correction;
            RT = bsxfun(@minus, RT, stimulus_onset);
            response = valid_responses;            
            correctness = (contains(response,'old') & contains(stimulus_type,'ClipA')) | (contains(response,'new') & contains(stimulus_type,'ClipB')) ;
                       
        elseif any(ismember(task_label,{'M1','M2'}))
            % Timing of event boundaries: (in sec)
%             stimulus = {'b0';'b1';'b2';'b3';'b4';'b5';'b6'}; 
%             stimulus_type = {'B';'B';'B';'B';'B';'B';'B'};
%             stimulus_onset = [0; 39.534; 121.087; 203.495; 242.895; 333.215; 364.815];
          
            epochlength = 4;
            % stimulus_onset = [epochlength:epochlength:(blockDuration-epochlength)]'; % non-overlaping X-sec-separated rest markers for epoching           stimulus = cellstr(string('movepoch')+(1:length(stimulus_onset)))';  
                     
            % Timing of original events in movie: (in sec)
            [origEventTime,sortind] = sort(eventsInMovie.timeInSec);
            origEventName = cellstr(string(task_label) + eventsInMovie.name(sortind));             
            % add the inter-events intervals:
            IEIsTime = [];
            for ii = 1:length(origEventName)-1
                IEIsTime = [IEIsTime, origEventTime(ii)+epochlength:epochlength:(origEventTime(ii+1)-epochlength)];
            end
            IEIsName = cellstr(string('IEI')+(1:length(IEIsTime)))';            
            allEventsTime = [origEventTime; IEIsTime'];
            allEventsName = [origEventName; IEIsName]; 
            
            stimulus_onset = allEventsTime;
            stimulus = allEventsName;
            stimulus_type = repmat({task_label},[length(stimulus_onset),1]);
            response = repmat({'n/a'},[length(stimulus_onset),1]);
            correctness = nan(length(stimulus_onset),1);              
            RT = ones(length(stimulus_onset),1) .* epochlength;
        
        elseif any(ismember(task_label,{'R1','R2'}))
            % Extract ten equally spaced resting-state epochs for reference: (in sec)
            epochlength = 4;
            stimulus_onset = [epochlength:epochlength:(blockDuration-epochlength)]'; % non-overlaping X-sec-separated rest markers for epoching
            stimulus = cellstr(string('restepoch')+(1:length(stimulus_onset)))';
            stimulus_type = repmat({task_label},[length(stimulus_onset),1]);
            response = repmat({'n/a'},[length(stimulus_onset),1]);
            correctness = nan(length(stimulus_onset),1);              
            RT = ones(length(stimulus_onset),1) .* epochlength;
       
        else
            fprintf('\nSkipping dataset: %s\n',EEG.setname); 
            continue;
        end
        
        taskid = repmat({task_label},[size(stimulus,1),1]);
           
        % load detected ripples:
        [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,0);
        swr_time = [];
        for ii = 1:numel(hippocampus_all_channels)
            channelid = hippocampus_all_channels{ii};
            filename = dir(fullfile(ripplesdir,subjid,sprintf('%s %s ripples %s.mat',subjid,task_label,channelid)));
            fprintf('\nProcessing channel %s... \n',hippocampus_all_channels{ii})
            if isempty(filename), error('Could not find ripples data for this channel'); end
            tmp = load(fullfile(ripplesdir,subjid,filename.name));
            swr_time = tmp.ripples;
            swr_time.subjid = repmat({subjid},size(swr_time.peak));
            clear tmp;            
            
            T = -10 : (1/DATA.Fs) : EEG.times(finind)/1000 + 10; % time vector (+ 10 sec margins)
            epochlim = [-8 8];
            DATA.rippletable = cat(1,DATA.rippletable,swr_time); % concatenate ripples table
            ripplevec = zeros(size(T));
            if ~isempty(swr_time)
                for k=1:numel(swr_time.peak)
                    current_event = swr_time.peak(k);
                    [~,idx]=min(abs(T-current_event));
                    ripplevec(idx(1))=swr_time.amplitude(k);
                end
            end

            switch time_locking_event
                case 'stimonset', [epocheddata, newtime, included_ind] = epoch( ripplevec, stimulus_onset-T(1), epochlim,'srate',DATA.Fs);
                case 'rt',    [epocheddata, newtime, included_ind] = epoch( ripplevec, stimulus_onset+RT-T(1), epochlim,'srate',DATA.Fs);
            end
            epocheddata = squeeze(epocheddata)';
            if size(epocheddata,2)==1, epocheddata=epocheddata'; end
            t = newtime(1):(1/DATA.Fs):newtime(2); % in sec

            DATA.trialtime = t;
            fieldtxt = task_label;                
                     
            % continue here:
            if ~isempty(epocheddata)
                DATA.Rall.raster = cat(1,DATA.Rall.raster,epocheddata);
                DATA.Rall.taskid = cat(1,DATA.Rall.taskid,taskid(included_ind));
                DATA.Rall.subjid = cat(1,DATA.Rall.subjid,repmat({subjid},size(epocheddata,1),1));
                DATA.Rall.channelid = cat(1,DATA.Rall.channelid,repmat({[channelid '_' subjid]},size(epocheddata,1),1));
                DATA.Rall.stimulus_type = cat(1,DATA.Rall.stimulus_type,stimulus_type(included_ind));
                DATA.Rall.stimulus = cat(1,DATA.Rall.stimulus,stimulus(included_ind));
                DATA.Rall.RT = cat(1,DATA.Rall.RT,RT(included_ind));
                DATA.Rall.response = cat(1,DATA.Rall.response,response(included_ind));
                DATA.Rall.correct = cat(1,DATA.Rall.correct,correctness(included_ind));
            end
        end
    end
    %% save the data:
    filename = fullfile(ripplesdir,subjid,sprintf('%s_ripple_psth_data_ref_%d_%s_all_hippocampal_channels.mat',subjid,ref_flag,time_locking_event));
    save(filename,'DATA') % save the data
    fprintf('\ndata saved - \n%s\n',filename)
end

