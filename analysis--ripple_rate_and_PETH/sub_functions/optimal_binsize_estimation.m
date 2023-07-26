function [bs] = optimal_binsize_estimation(Nr,sd,trial_duration,Fs,binsize_method)
% Compute the optimal binsize. 
% input:
% Nr - overall number of ripples/spikes
% sd - standard dev.
% trial_duration - in sec.
% Fs - sampling freq.
% binsize_method - 'scott', 'sqrt', 'shimazaki'
% return optimal binsize in samples

if strcmpi(binsize_method,'scott')
    bs = round((3.49*sd/(Nr^(1/3)))*Fs); % in samples
    % Original paper: David W. Scott 1979, On optimal and Data-Based Histograms
elseif strcmpi(binsize_method,'sqrt')
    k = ceil(sqrt(Nr));
    bs = round((trial_duration/k)*Fs); % in samples
elseif strcmpi(binsize_method,'shimazaki')
    minN = 10; % require at least 10 bins
    [~,bs,~,~,~]=sshist(t(x0),minN:100);  % convert to samples
    % Original paper: Shimazaki and Shinomoto (2007) http://dx.doi.org/10.1162/neco.2007.19.6.1503
end

end