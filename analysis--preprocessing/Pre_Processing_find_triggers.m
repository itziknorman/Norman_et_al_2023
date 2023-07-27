%% Find triggers:

% This script runs peak detection on the trigger channel to find the triggers 
% that were sent during the task:

trigChannel = false(EEG.pnts,1);
trigChannel([EEG.event.latency]) = true;

% Alternative (in some patients):
%  trigChannel=EEG.data(strcmpi({EEG.chanlocs.labels},'EVENT'),:);

figure; plot(trigChannel)

% Get scenario times:
scenario = 'MT';
[EraseWrongTriggers,cfg] = subjects_exp_times_PP(EEG,initials,scenario);

%  Clear unrelevant time intervals:
if ~isempty(EraseWrongTriggers)
    for i=1:size(EraseWrongTriggers,1)
        trigChannel(EraseWrongTriggers(i,1):EraseWrongTriggers(i,2))=0;
    end
end

% Parameteres:

path=fullfile(maindir,'log');
switch scenario
    case 'R1'
        [logFileName]=dir(fullfile(path,['*-rest_1*']));
    case 'M1'
        [logFileName]=dir(fullfile(path,['*-movie_1*']));
    case 'R2'
        [logFileName]=dir(fullfile(path,['*-rest_2*']));
    case 'M2'
        [logFileName]=dir(fullfile(path,['*-movie_2*']));
    case 'MT'
        [logFileName]=dir(fullfile(path,['*-memory_test_1*']));
end


if isempty(logFileName)
    error('no file');
else
    logFileName=logFileName.name;
end
[event_labels,event_presentation_onsets]=EventTriggerCoupling(triggers,fullfile(path,logFileName),initials,EEG.times);

