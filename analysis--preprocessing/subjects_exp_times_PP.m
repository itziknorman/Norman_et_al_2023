
function [relevantTrig, logFileName] = subjects_exp_times_PP(EEG,subjid,scenario,maindir)

logFileName = [];

% Conditions duration:
RESTduration = 180.15; % sec
MOVIEduration = 368.5; % sec
MTduration = 324.620;  % sec

switch subjid
    
    case 'PP08'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'16:08:42'},{'16:13:18'},245,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
                
            case 'M1'
                [onsetind,offsetind] = time2indices({'16:08:42'},{'16:18:18'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'16:08:42'},{'16:25:24'},245,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                [onsetind,offsetind] = time2indices({'16:08:42'},{'16:30:40'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'16:08:42'},{'16:42:54'},340,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        % =========================================================================
        % =========================================================================
        
    case 'PP06'
        
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'10:37:40'},{'10:42:40'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
                
            case 'M1'
                [onsetind,offsetind] = time2indices({'10:37:40'},{'10:46:10'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'10:37:40'},{'10:52:38'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                [onsetind,offsetind] = time2indices({'10:37:40'},{'10:55:50'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'10:37:40'},{'11:03:10'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        % =========================================================================
        % =========================================================================
        
    case 'PP04'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'08:47:12'},{'08:58:40'},180,EEG);
                EEG.data(trig_ch,onsetind:offsetind) = 0;
                EEG.data(trig_ch,onsetind) = 1;
                EEG.data(trig_ch,offsetind) = 1;
                EraseWrongTriggers=[1 onsetind-1; offsetind+1 EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'08:47:12'},{'09:07:40'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'08:47:12'},{'09:14:25'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'08:47:12'},{'09:17:42'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_2*']));
                logFileName=fullfile(path,logFileName.name);
                [onsetind,offsetind] = time2indices({'08:47:12'},{'09:28:32'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        % =========================================================================
        % =========================================================================
        
    case 'PP03'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'13:42:13'},{'13:48:32'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'13:42:13'},{'13:51:44'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'13:42:13'},{'13:58:07'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'13:42:13'},{'14:01:24'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'13:42:13'},{'14:08:22'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        
        % =========================================================================
        % =========================================================================
        
    case 'PP01'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'14:14:39'},{'14:22:58'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'14:14:39'},{'14:26:14'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'14:14:39'},{'14:32:30'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'14:14:39'},{'14:35:38'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'14:14:39'},{'14:42:47'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        % =========================================================================
        % =========================================================================
        
    case 'PP02'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'16:32:35'},{'16:36:29'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'16:32:35'},{'16:39:52'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'16:32:35'},{'16:46:16'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'16:32:35'},{'16:49:36'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+1500 offsetind-1500;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'16:32:35'},{'16:57:00'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        % =========================================================================
        % =========================================================================
        
    case 'PP05'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'14:28:37.000'},{'14:29:41.440'},RESTduration,EEG);
                EEG.data(trig_ch,onsetind) = 1;
                EEG.data(trig_ch,offsetind) = 1;
                EraseWrongTriggers=[1 onsetind-1; offsetind+1 EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'14:28:37.000'},{'14:32:50.380'},MOVIEduration,EEG);
                EEG.data(trig_ch,onsetind) = 1;
                EEG.data(trig_ch,offsetind) = 1;
                EraseWrongTriggers=[1 onsetind-1; offsetind+1 EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'14:28:37.000'},{'14:39:06.470'},RESTduration,EEG);
                EEG.data(trig_ch,onsetind) = 1;
                EEG.data(trig_ch,offsetind) = 1;
                EraseWrongTriggers=[1 onsetind-1; offsetind+1 EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                [onsetind,offsetind] = time2indices({'14:28:37.000'},{'14:42:15.380'},MOVIEduration,EEG);
                EEG.data(trig_ch,onsetind) = 1;
                EEG.data(trig_ch,offsetind) = 1;
                EraseWrongTriggers=[1 onsetind-1; offsetind+1 EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                logFileName=fullfile(path,'n/a');
                
                [onsetind,offsetind] = time2indices({'14:28:37.000'},{'14:49:03.690'},MTduration,EEG);
                EEG.data(trig_ch,onsetind) = 1;
                EEG.data(trig_ch,offsetind) = 1;
                EraseWrongTriggers=[1 onsetind-1; offsetind+1 EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        % ==================================================================
        
    case 'PP07'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'15:14:23'},{'15:17:53'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'15:14:23'},{'15:21:04'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'15:14:23'},{'15:27:23'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'15:14:23'},{'15:30:38'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'15:14:23'},{'15:38:09'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        % =================================================================
        % =================================================================
        
    case 'PP09'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'08:23:07'},{'08:37:34'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'08:23:07'},{'08:40:51'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'08:23:07'},{'08:47:08'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'08:23:07'},{'08:50:21'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'08:23:07'},{'08:57:35'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        % =================================================================
        % =================================================================
        
    case 'PP10'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'15:46:00'},{'15:54:40'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'15:46:00'},{'15:59:07'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'15:46:00'},{'16:05:32'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'15:46:00'},{'16:08:41'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'15:46:00'},{'16:15:40'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        % =================================================================
        % =================================================================
        
    case 'PP11'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'08:46:16'},{'08:49:50'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'08:46:16'},{'08:53:01'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'08:46:16'},{'08:59:39'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'08:46:16'},{'09:03:14'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+2500 offsetind-2500;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'08:46:16'},{'09:10:58'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
    case 'PP11'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'08:46:16'},{'08:49:50'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'08:46:16'},{'08:53:01'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'08:46:16'},{'08:59:39'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'08:46:16'},{'09:03:14'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+2500 offsetind-2500;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'08:46:16'},{'09:10:58'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
    case 'PP15'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'14:43:18'},{'14:47:19'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'14:43:18'},{'14:51:09'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'14:43:18'},{'14:57:35'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                [onsetind,offsetind] = time2indices({'14:43:18'},{'15:01:12'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+2500 offsetind-2500;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'14:43:18'},{'15:08:55'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
    case 'PP16'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'16:10:56'},{'16:17:43'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'16:10:56'},{'16:21:58'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'16:10:56'},{'16:28:24'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                [onsetind,offsetind] = time2indices({'16:10:56'},{'16:31:41'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+2500 offsetind-2500;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'16:10:56'},{'16:40:31'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
    case 'PP17'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'15:30:18'},{'15:40:00'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'15:30:18'},{'15:43:27'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'15:30:18'},{'15:50:07'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                [onsetind,offsetind] = time2indices({'15:30:18'},{'15:53:17'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+2500 offsetind-2500;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'15:30:18'},{'16:00:57'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
    case 'PP18'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        switch scenario
            case 'R1'
                [onsetind,offsetind] = time2indices({'14:06:37'},{'14:09:19'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                [onsetind,offsetind] = time2indices({'14:06:37'},{'14:12:37'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                [onsetind,offsetind] = time2indices({'14:06:37'},{'14:18:59'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                [onsetind,offsetind] = time2indices({'14:06:37'},{'14:22:17'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+2500 offsetind-2500;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                
                [onsetind,offsetind] = time2indices({'14:06:37'},{'14:29:26'},335,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
      
        case 'SUB01ASSUTA'
        trig_ch = find(strcmpi({EEG.chanlocs.labels},'TRIG'));
        cond = {'prerf-1.1','prerf-1.2','postrf-1.1','postrf-1.2','postrf-2.1','postrf-2.2','postrf-2.3'};

        switch scenario
            
             case cond(1:7)          
                onsetind = 1; 
                offsetind = min(EEG.pnts,900*1000); % take only the first 15 minutes (14.9 min in the first dataset)
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+1500 offsetind-1500;];
                relevantTrig = EEG.data(trig_ch,:);
%                 relevantTrig(onsetind)=1; relevantTrig(offsetind)=1;
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
                           
            
            case 'R1'            
                onsetind = 1633287; 
                offsetind = 1633287 + 185*1000; 
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+1500 offsetind-1500;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M1'
                onsetind = 2170264; 
                offsetind = 2170264 + 375*1000; 
                %[onsetind,offsetind] = time2indices({'15:34:11'},{'16:08:50'},375,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+5000 offsetind-5000;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'R2'
                onsetind = 2572578; 
                offsetind = 2572578 + 185*1000; 
                %[onsetind,offsetind] = time2indices({'15:32:30'},{'16:16:00'},185,EEG);
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+40000 offsetind-10000;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'M2'
                
                onsetind = 2791206; 
                offsetind = 2791206 + 375*1000; 
%                 [onsetind,offsetind] = time2indices({'15:32:30'},{'16:19:00'},375,EEG);       
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts; onsetind+25000 offsetind-25000;];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                relevantTrig(EraseWrongTriggers(3,1):EraseWrongTriggers(3,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
            case 'MT'
                % log file for responses:
                path = fullfile(maindir, 'log');
                [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
                logFileName=fullfile(path,logFileName.name);
                onsetind = 3258637; 
                offsetind = 3258637 + 335*1000;    
                %[onsetind,offsetind] = time2indices({'15:32:30'},{'16:27:00'},335,EEG);    
                EraseWrongTriggers=[1 onsetind; offsetind EEG.pnts];
                relevantTrig = EEG.data(trig_ch,:);
                relevantTrig(EraseWrongTriggers(1,1):EraseWrongTriggers(1,2))=0;
                relevantTrig(EraseWrongTriggers(2,1):EraseWrongTriggers(2,2))=0;
                figure; hold on; plot(EEG.data(trig_ch,:),'b'); plot(relevantTrig,'r','linewidth',3);
        end
        
        
        
end

end



