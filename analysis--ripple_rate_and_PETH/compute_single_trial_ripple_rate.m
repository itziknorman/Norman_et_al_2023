function [RR,trialDur] = compute_single_trial_ripple_rate(DATA,Fs,timeWindow,time_locking_event)
% compute single trial ripple rate:

% "timeWindow" selects the time window for analysis relative to "time_locking_event":
% [] => the entire trial duration
% [>1] => interval from [rt-timeWindow] to response
% [<-1] => interval from onset to [rt-timeWindow]
% [t1 t2] => a specific interval


trialDur = [DATA.Rall.RT];

assert(all(~isnan(trialDur))); % make sure no NaNs
RR = [];
for k = 1:length(trialDur)
    if isempty(timeWindow)
        switch time_locking_event
            case 'stimonset',  timeind = DATA.trialtime > 0 & DATA.trialtime <= trialDur(k);
            case 'rt', timeind = DATA.trialtime <= 0 & DATA.trialtime > -trialDur(k);
        end
    elseif length(timeWindow)==2
        assert(timeWindow(2)>timeWindow(1))      
        % use a fixed, predefined time window:
        timeind = DATA.trialtime <= timeWindow(2) & DATA.trialtime > timeWindow(1);       
        
    elseif length(timeWindow)==1 && abs(timeWindow)>=1
        % trial duration, if negative: [0.5 to rt-timeWindow]; if positive: [-timeWindow to RT])
        win = abs(timeWindow);
        trialdur = trialDur(k);
        switch time_locking_event
            case 'stimonset'
                if timeWindow>0, timeind = DATA.trialtime > 0 &  DATA.trialtime <= trialdur;
                else, timeind = DATA.trialtime > min(0.5,trialdur) &  DATA.trialtime <= 0; end
            case 'rt'
                if timeWindow>0, timeind = DATA.trialtime <= 0 & DATA.trialtime > -win;
                else, timeind = DATA.trialtime > min(0.5,-trialdur) & DATA.trialtime <= -win; end
        end
        
    else
        error('wrong input, check the value of "timeWindow"')
    end
    % compute ripple rate:
    winlength = sum(timeind)./Fs;
    RR(k,1) = sum(logical(DATA.Rall.raster(k,timeind)),2)./winlength;  % compute ripple rate (events/sec)
    
end

fprintf('\n single-trial SWR rate was computed (O.K.) \n');

