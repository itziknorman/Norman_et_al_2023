%% Find triggers:
triggers = [EEG.event.latency];
y = zeros(size(EEG.times));
y(triggers) = 1;
index = find(y); index_diff = [diff(index)];
figure('position',[0 0 2000 500]); hold on; bar(index_diff,'k');

% Find relevant event based on inter-trigger intervals:
Fs = EEG.srate;
movie_duration = 368; % in seconds
clip_duration = 8; % in seconds
rest_duration = 180; % in seconds
load(fullfile('D:\ECoG\Paris\Pink_Panther_Data','clip_order.mat'));

eventIndex = struct;
clipind = sprintf('%d',abs(index_diff - (8 * Fs)) < 50);
clipind = regexp(clipind, '1{39,}', 'start'):regexp(clipind, '1{39,}', 'end')+1;
for ii = 1:numel(clipind)
    eventIndex.(clip_order{ii}) = clipind(ii);
end
eventIndex.memorytest_str = eventIndex.(clip_order{1}) - 1;
eventIndex.memorytest_fin = eventIndex.(clip_order{end}) + 1;
eventIndex.movie1_str = find(abs(index_diff - (368 * Fs)) < 500,1,'first'); % find first movie
eventIndex.movie1_fin = eventIndex.movie1_str + 1;

eventIndex.movie2_str = find(abs(index_diff - (368 * Fs)) < 500,1,'last'); % find second movie
eventIndex.movie2_fin = eventIndex.movie2_str + 1;
  
% resting state:
rest_periods = find(abs(index_diff - (180 * Fs)) < 500,3,'first');

% =============================================
assert(rest_periods(1) < eventIndex.movie1_str)
% =============================================

eventIndex.rest1_str = rest_periods(1); % find 1st rest
eventIndex.rest1_fin = eventIndex.rest1_str + 1; 

% =============================================
assert(rest_periods(2) < eventIndex.movie2_str)
% =============================================

eventIndex.rest2_str = rest_periods(2); % find 2nd rest
eventIndex.rest2_fin = eventIndex.rest2_str + 1; 

if numel(rest_periods)>2
    % =============================================
    assert(rest_periods(3) > eventIndex.memorytest_fin)
    % =============================================
    eventIndex.rest3_str = rest_periods(3); % find 3rd rest
    eventIndex.rest3_fin = eventIndex.rest3_str + 1; 
end 

