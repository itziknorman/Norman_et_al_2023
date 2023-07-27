function [event_list_responses,event_list_clips]=EventTriggerCouplingMemoryTest(triggers,logfilePath,subjid,t)

% this function reads the triggers' onset from the Presentation log file
% and produces event lists that can be imported to eeglab

fileID = fopen(logfilePath);

% read the first 8 columns:
C = textscan( fileID , '%s %s %s %s %s %s %s %s' );

idx=multiStrFind(C{3},{'Video','Response'})&multiStrFind(C{4},{'2','3','Clip'});

tmp={};
col=0;
for i=1:numel(C)
    col=col+1;
    row=0;
    for j=find(idx)'
        row=row+1;
        tmp{row,col}=C{i}{j};
    end
end
T = cell2table(tmp,'VariableNames',{'Subject','Trial','EventType','Code','Time','TTime','Uncertainty','Duration'});

event_name=[];
event_index=[];
event_time=[];
for i=1:size(T,1)
    tmp=T.EventType{i,1};
    if ismember(tmp,{'Video'})
        event_time(end+1,1)=double(str2num(T.Time{i})/10);
        event_index(end+1,1)=i;
        event_name{end+1,1}=T.Code{i};
        %  if strcmp(T.Code{i},'DING'), break; end
    end
end

% Conversion to presentation time (based on the first 10 triggers):

x = t(triggers(1:10));
y = event_time(1:10)';
[r,m,b] = regression(x,y);
event_time=event_time-b;


event_list_responses={};
e = 1;
for k=1:size(T,1)
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
        event_list_responses{e,2} = (((double(str2num(T.Time{k}))/10)) - b)./1000;  % convert to EEGLAB time;
        e=e+1;
    end
end

event_list_clips={};
e = 1;
for k=1:size(T,1)
    tmp=T.EventType{k,1};
    if ismember(tmp,{'Video'})
        l = T.Code{k};
        if isempty(l),continue; end
        event_list_clips{e,1} = l;
        event_list_clips{e,2} = (((double(str2num(T.Time{k}))/10)) - b)./1000;  % convert to EEGLAB time;
        e=e+1;
    end
end

