%% This script cut the experiment into the various conditions: REST x2, MOVIE x2, and MEMORY TEST 


clear all
close all
clc;

% set the relevant path:
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;
clc; warning('off');

addpath(fullfile(path_to_toolboxes,'eeglab14_1_2b'));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));
rmpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));

[ALLEEG, EEG, CURRENTSET, ~] = eeglab;

subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};

for subjid = subjects(1:end) 
    subjid = cell2mat(subjid);
    close all;
    clearvars -except subjects subjid ALLEEG EEG parentfolder path_to_toolboxes
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    SCENARIOS={'R1', 'M1', 'R2', 'M2','MT'};
    
    % Define Directories:
    maindir=fullfile(parentfolder,subjid);
    outdir=fullfile(maindir,subjid,'EEGLAB_datasets');
    ref_flag=2; % 1 = Common Ref; 2 = Bipolar montage;
    switch ref_flag
        case 1
            indir=fullfile(maindir,'EEGLAB_datasets','preprocessed');
            outdir=fullfile(maindir,'EEGLAB_datasets');
            file_ending='preprocessed';
        case 2
            indir=fullfile(maindir,'EEGLAB_datasets_BP','preprocessed');
            outdir=fullfile(maindir,'EEGLAB_datasets_BP');
            file_ending='preprocessed_BP_montage';
    end
    
    for scenario = SCENARIOS(1:5)
        
        scenario=cell2mat(scenario);
        clearvars -except scenario subjid maindir indir outdir ALLEEG EEG file_ending parentfolder
        CURRENTSET=1;
        % define in/out filenames:
        inFileName=fullfile(indir,[subjid '_pink_panther_' file_ending '.set']);
        outFileName=[subjid '_pink_panther_' file_ending '_' scenario '.set'];
        
        % LOAD preprocessed dataset:
        
        [EEG] = pop_loadset('filename', inFileName);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
        eeglab redraw
        
        
        %% Find the Triggers:
        
        trigChannel=EEG.data(strcmpi({EEG.chanlocs.labels},'TRIG'),:);
        % Get scenario times:
        [relevantTrig,logFileName] = subjects_exp_times_PP(EEG,subjid,scenario,maindir);
     
        figure('units','normalized','outerposition',[0 0 0.5 0.5]); hold on; 
        plot(trigChannel,'b','linewidth',3); plot(relevantTrig,'r-','linewidth',1);
        hold on; plot(get(gca,'xlim'),[0.8 0.8],'r--'); axis tight
        fprintf('\n PRESS ANY KEY TO CONTINUE... \n')

        %%
              
        % Parameteres:       
        mph=0.8;  % minimal trigger peak height in standard deviations (e.g. 5 or 10)
        minDistance=20;
        [~,triggers] = findpeaks(double(relevantTrig), 'MINPEAKHEIGHT', mph,'MINPEAKDISTANCE',minDistance);
        figure; plot(EEG.times, relevantTrig,'k'); hold on; scatter(EEG.times(triggers),relevantTrig(triggers),'ro')
        
        
        %% Event-trigger coupling:        
             
        % PP08 fix:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     [~,keep]=unique(event_labels,'stable');
        %     event_labels=event_labels(keep);
        %     triggers=triggers(keep);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EEG = pop_editeventvals(EEG,'delete',1:numel(EEG.event));
        trigger_diff = diff(triggers);
        figure; bar(trigger_diff,'b'); 
        if ~strcmpi(scenario,'MT')               
            [~,ind1]=max(trigger_diff);
            ind2 = ind1+1;              
        else                  
            load(fullfile(parentfolder,'clip_order.mat'));     
            % exceptions:
            ind = trigger_diff > 7.6*EEG.srate;           
            if contains(subjid,'ASSUTA')
                clipind = find(ind);
                ind1 = clipind(1) - 1;
                ind2 = clipind(end) + 1;
            else
                clipind = sprintf('%d',ind);
                clipind = regexp(clipind, sprintf('1{%d,}',sum(ind)), 'start'):regexp(clipind, sprintf('1{%d,}',sum(ind)), 'end');  
                if strcmpi(subjid, 'PP08')                
                    ind1 = clipind(1) - 2;
                else
                    ind1 = clipind(1) - 1;
                end
                ind2 = clipind(end) + 1;
            end
            
            clipTriggers = triggers(clipind);
            
            % add response from log file (if existing):
            if ~exist(logFileName,'file')
                warning('no log file - compute from the events in the EDF');
                e=1;
                event_list_clips = {};
                for ii = 1:numel(clipind)
                    event_list_clips{e,1} = clip_order{ii};
                    event_list_clips{e,2} = EEG.times(triggers(clipind(ii)))./1000;
                    e=e+1;
                end
               
                % Responses recovered from video (xlsx file):
                path = fullfile(maindir, 'log');
                [xlsFileName]=dir(fullfile(path,['Memory_test_responses_from_video*']));
                xlsFileName=fullfile(path,xlsFileName.name);
                
                [~,~,D] = xlsread(xlsFileName);
                for jj = 1:numel(D)
                    if ~ischar(D{jj}), D{jj} = num2str(D{jj}); disp(D{jj}); end 
                end
                
                T = cell2table(D(2:end,:),'VariableNames',{'TriggerNum','EventType','Code','Time','RelTime','Rec_onset'});               
                [event_list_responses]=EventTriggerCouplingMemoryTest_without_logfile({'00:00:00'},T,EEG);
       
            else
                [event_list_responses,event_list_clips]=EventTriggerCouplingMemoryTest(clipTriggers,logFileName,subjid,EEG.times);
            end     
               
            
        end
            
        strind = triggers(ind1);
        finind = triggers(ind2);

        warning(['Sanity Check - Size of TC: ' num2str((finind-strind)/EEG.srate) ' seconds']);
        beep;  
        fprintf('\n PRESS ANY KEY TO CONTINUE... \n')
            
        % Add str/fin events: 
        event_list = {};       
        event_list{1,1} = 'str';
        event_list{1,2} = EEG.times(strind)./1000;
        event_list{2,1} = 'fin';
        event_list{2,2} = EEG.times(finind)./1000;
        EEG = pop_importevent( EEG, 'event', event_list, 'fields',{ 'type', 'latency'}, 'timeunit',1, 'append','no');     
        [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
          
        
        if strcmpi(scenario,'MT')
            EEG = pop_importevent( EEG, 'event', [event_list_clips; event_list_responses], 'fields',{ 'type', 'latency'}, 'timeunit',1, 'append','yes');
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency          
        end
      
        % Erase extra timepoints:
        tail=5; % in seconds
        ii = find(strcmpi({EEG.event.type},'str'));
        cut_pre=[0, EEG.event(ii).latency/EEG.srate-tail]; clear ii        
        ii = find(strcmpi({EEG.event.type},'fin'));
        cut_post=[EEG.event(ii).latency/EEG.srate+tail, EEG.xmax]; clear ii
        
        EEG = pop_select(EEG,'notime',[cut_pre; cut_post]);
        EEG.setname=[subjid '_' scenario '_preprocessed'];
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
       
        eegplot(EEG.data,'color','off','srate',EEG.srate,'winlength',30,'limits',[EEG.times(1) EEG.times(end)],'events',EEG.event)

        %% Save:
        EEG.setname=outFileName;
        EEG = pop_saveset( EEG, 'filename', outFileName, 'filepath', outdir,'version','7.3');
        disp('data saved')
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        eeglab redraw;
    end
end
    

