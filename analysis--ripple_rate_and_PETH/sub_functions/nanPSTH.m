function  [r,binscenters] = nanPSTH(data, binsize, wsize, fs, trialmask)
% The function computes the peri-stimulus time histogram (PSTH)from SWRs events.
% The function outputs the trial averaged ripple rate as a function of time.
% ===============================================================
% [r,binscenters]  = nanPSTH(data, binsize, wsize, fs, trialmask)
% ===============================================================
% data - events matrix (samples X trials); 1=event, 0=no-event, NaN=ignored timepoint
% binsize - binwidth (samples)
% wsize - size of triangular smoothing window (bins).
% fs - sampling rate (hz)
% trialmask - matrix, same size as 'data', specifying the timepoints that should be taken into account.
% NOTE: trialmask is optional. If you don't need replace it by empty vector [].
% Author: Itzik Norman
% normanik@gmail.com
% @ Weizmann Institute, Israel

narginchk(5,5);
nargoutchk(1,3);
trialmask = []; % (!)
if isempty(trialmask), trialmask = ones(size(data)); end
triallen = sum(trialmask,2);
% Compute PSTH:
binsizems = floor(binsize/fs*1000); % in msec
lastBinms = floor(binsizems * floor(size(data,1)./binsize));
% t = 1/fs*1000 : 1/fs*1000 : size(data,1)/fs*1000;
t = linspace(0, size(data,1)/fs*1000,  size(data,1));
t = ceil(t);
edges = 0 : binsizems : lastBinms;
edges = ceil(edges);
binscenters = edges(1:end-1)+(binsizems/2);
r = nan(size(binscenters));
for k = 1:length(binscenters)
    timeind = t>=edges(k) & t<edges(k+1);
    assert((sum(timeind)-binsize)<ceil(1/fs*1000));
    nSamples = nansum(nansum(trialmask(timeind,:),2));
    if nSamples == 0, continue; end % require at least 1 trial to compute PSTH
    r(k) = nansum(nansum(data(timeind,:),2)) ./ (nSamples/fs);
end
if nargout < 2, clear binscenters; end
%fprintf('\n PSTH succesfully computed \n');
r = nanfastsmooth(r,wsize,2,0.5);



