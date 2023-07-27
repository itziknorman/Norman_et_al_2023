function [event_list_responses]=EventTriggerCouplingMemoryTest_without_logfile(file_onset,T,EEG)

% this function reads the triggeres' onset from a table (only relevant when the Presentation log
% file is missing) and produces event lists that can be imported to eeglab


event_list_responses={};
e = 1;
for k = 1:size(T,1)
    tmp=T.EventType{k,1};
    if ismember(tmp,{'Response'})
        tmp = T.Code{k};
        switch tmp
            case '2'
                l = 'old';
            case '3'
                l = 'new';
            otherwise
                l = [];
        end
        if isempty(l),continue; end
        event_list_responses{e,1} = l;
        tt = datestr(str2num(T.RelTime{k}),'HH:MM:SS');
        [ind,~] = time2indices(file_onset,{tt},0,EEG);
        tt=EEG.times(ind)/1000;
        event_list_responses{e,2} = tt;  % convert to EEGLAB time;
        e=e+1;
    end
end

