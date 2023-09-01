
function [ripples,ripples_stat,r]=ripples_detection_algorithm(signal,BP,t,Fs,th,minDistance,minRippleDuration,maxRippleDuration,noise_ch,fastRippleSignal,IED_onsets,boundary_onsets)
% detecting ripples based on Strack et al. 2014 method (doi: 10.1016/j.neuron.2014.06.023)  
% signal: normalized squared bandpassed LFP signal from the hippocampus         
% BP: raw bandpass signal
% t: corresponding time vector in sec
% th: thresholds for event expansion and detection in stdev [onset/offset peak]
% minDistance: in sec
% minRippleDuration; min ripple duration in Sec
% maxRippleDuration: max ripple duration in Sec
% noise_ch: channel to exclude muscular/electrical artifacts identifed as ripples
% fastRippleSignal: (optional) bandpass signal to exclude fast-ripple-band transients that coincide with rippples
% IED_onsets: timestamps of IEDs, detected using the linelength algorithm & IED-band ampltiude transients detector
% boundary_onsets: recording boundaries to avoid edge artifacts
%
% Author: Yitzhak Norman 10/03/20

if size(signal,2)<size(signal,1),signal=signal'; end
if size(noise_ch,2)<size(noise_ch,1),noise_ch=noise_ch'; end
if size(fastRippleSignal,2)<size(fastRippleSignal,1),fastRippleSignal=fastRippleSignal'; end

[pks,locs] = findpeaks(signal, 'MINPEAKHEIGHT', th(2));
ENV = BP;
ENV(~isnan(ENV)) = abs(hilbert(BP(~isnan(BP)))).^2;   % to calculate ripple amplitude (squared)

% Rejection based on a control reference channel:
if ~isempty(noise_ch)
    [r,p]=corr(signal',noise_ch','rows','pairwise');
    fprintf('\n --- correlation with noise ch.: %.2f \n',r)
    rej_count=0;
    pre_rej_count=size(locs,2);    
    fprintf('\n => rejecting noise artefacts... \n')
    [~,noise_locs] = findpeaks(noise_ch, 'MINPEAKHEIGHT', th(2));
    %plot_event_triggered_spectrogram(noise_ch,Fs,noise_locs,[-1 1],[-10 10],'noise channel events')
    
    % Ignore global electrical artifacts:
    for i=1:numel(noise_locs)
        tmp=find(abs(locs-noise_locs(i))< (0.25 * Fs));
        if ~isempty(tmp)
            locs(tmp)=[];
            pks(tmp)=[];
            rej_count=rej_count+1;
        end
    end
    fprintf('\n *** rejected %d / %d events based on noise channel correlation \n',rej_count,pre_rej_count)
    ripples_stat.noise_rejection = rej_count/pre_rej_count;
end

% fast ripple rejection:
if ~isempty(fastRippleSignal)
    rej_count=0;
    pre_rej_count=size(locs,2);    
    fprintf('\n => rejecting fast ripple events... \n')
    [~,fast_ripples_locs] = findpeaks(fastRippleSignal, 'MINPEAKHEIGHT', th(2)); % use the same detection creteria as ripples 
    %plot_event_triggered_spectrogram(fastRippleSignal,Fs,setdiff(fast_ripples_locs,noise_locs),[-1 1],[-10 10],'fast ripple events')
    % Ignore ripples coninciding (<50ms) with pathological fast ripples:
    for i=1:numel(fast_ripples_locs)
        tmp=find(abs(locs-fast_ripples_locs(i))< (0.05 * Fs));
        if ~isempty(tmp)
            locs(tmp)=[];
            pks(tmp)=[];
            rej_count=rej_count+1;
        end
    end
    fprintf('\n *** rejected %d / %d events that coincide with fast ripples \n',rej_count,pre_rej_count)
    ripples_stat.fast_ripples_rejection = rej_count/pre_rej_count;
end

% IED rejection: 
% METHOD I: IEDs detected based on Estellar et al 2001 (DOI:10.1109/IEMBS.2001.1020545)
% METHOD II: IEDs detected based on Smith et al 2022 (DOI: 10.7554/eLife.73541)
if any(IED_onsets)
    rej_count=0;
    pre_rej_count=size(locs,2); 
    % Ignore IED-ripples events (candidate ripples that coincide within 1s window centered on the IED):
    IED_locs = find(IED_onsets);    
    for i=1:numel(IED_locs)
        tmp=find(abs(locs-IED_locs(i))< (0.5 * Fs));
        if ~isempty(tmp)
            locs(tmp)=[];
            pks(tmp)=[];
            rej_count=rej_count+1;
        end
    end
    fprintf('\n *** rejected %d / %d events that coincide with IEDs \n',rej_count,pre_rej_count)
    ripples_stat.IED_rejection = rej_count/pre_rej_count;
end

% Recording boundaries rejection (to avoid edge-related artiacts): 
if any(boundary_onsets)
    rej_count=0;
    pre_rej_count=size(locs,2); 
    % Ignore IED-ripples events (candidate ripples that coincide within 2s window centered on the IED):
    boundary_locs = find(boundary_onsets);    
    for i=1:numel(boundary_locs)
        tmp=find(abs(locs-boundary_locs(i))< (2 * Fs));
        if ~isempty(tmp)
            locs(tmp)=[];
            pks(tmp)=[];
            rej_count=rej_count+1;
        end
    end
    fprintf('\n *** rejected %d / %d events that were close to file boundaries \n',rej_count,pre_rej_count)
    ripples_stat.boundary_rejection = rej_count/pre_rej_count;
end

% Exclusion based on duration creteria:
counter=1;
ripples=nan(1,4);
ripples=array2table(ripples,'VariableNames',{'str','peak','fin','amplitude'});
ripples(1,:)=[];
for k=locs
  
    % find the starting point of the peak:
    stop=0;
    str=k;
    while ~stop && ~str==0
        if str==1
            break
        end
        str=str-1;
        if signal(str)<th(1), stop=1; end
    end
    
    % find the ending point of the peak:
    stop=0;
    fin=k;
    while ~stop && ~(fin==numel(signal))
        fin=fin+1;
        if signal(fin)<th(1), stop=1; end
    end
    
    % Detect negative peak position for each ripple (closest to ripple's power peak)
    minIndex = [];
    [~,minpos] = findpeaks(-double(BP(str:fin)));
    if isempty(minpos), [~,minpos] = min(BP(str:fin)); minpos = minpos(1); end
    [~,maxamp] = max(double(ENV(str:fin)));
    minpos=minpos-1; 
    maxamp=maxamp-1;
    [~,tmp] = min(abs(minpos-maxamp));
    minIndex=minpos(tmp);
    peakPosition = min((str + minIndex),numel(signal));  
    
    try
        M = [t(str), t(peakPosition), t(fin), ENV(peakPosition)];
        M = round(M.*1000)./1000;
        ripples(counter,:)=array2table(M);
    catch
        disp(ripples);
        fprintf('\n Error has occured in event # %d \n',counter);
    end
    counter=counter+1;
end
disp(['After detection by thresholding: ' num2str(size(ripples,1)) ' events.']);
if isempty(ripples),return; end


% Merge ripples if inter-ripple period is less than minDistance:
ripples_edit=ripples;
rej=zeros(size(ripples,1),1);
for k = 2:size(ripples,1)
    if (ripples.peak(k)-ripples.peak(k-1)) < minDistance        
        ripples_edit.fin(k-1) = ripples.fin(k); % Merge the events
        rej(k)=1;
    end
end

if any(rej), ripples_edit(find(rej),:)=[]; end
ripples=ripples_edit;
disp(['After ripple merge: ' num2str(size(ripples,1)) ' events.']);
if isempty(ripples),return; end

% duration test:
duration = ripples.fin-ripples.str;
ripples(duration<minRippleDuration,:) = [];
disp(['After min duration test: ' num2str(size(ripples,1)) ' events.']);
duration = ripples.fin-ripples.str;
ripples(duration>maxRippleDuration,:) = [];
disp(['After max duration test: ' num2str(size(ripples,1)) ' events.']);

end
