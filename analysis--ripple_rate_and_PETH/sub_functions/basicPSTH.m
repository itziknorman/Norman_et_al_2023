function  [binsizems,binscenters,r,edges] = basicPSTH(eventind, binsize, wsize, fs, ntrials, triallen)
% PSTH Computes the peri-stimulus time histogram from spike times.
% The routine plots the trial averaged spike rate as a function of time.
% H = PSTH(TIMES, BINSIZE, FS,NTRIALS,TRIALLEN)
% TIMES - spike times (samples)
% BINSIZE - binwidth (samples)
% wsize - size of triangular smoothing window (bins).
% FS - sampling rate (hz)
% NTRIALS - number of trials
% TRIALLEN - length of a trial (samples)
%
% An example:
% spike times specified in continuous time -
% here we have 2 trials and a trial length of 1000 samples
% t = [10, 250, 900, 1300, 1600, 2405, 2900];
% r = psth(t,10,5,1000,3,1000) ;
%
% Author: Itzik Norman
% normanik@gmail.com
% @ Weizmann Institute, Israel

narginchk(6,6);
nargoutchk(1,4);

% Compute PSTH:
% lastBin = binsize * ceil((triallen)*(1000/(fs*binsize)));
triallenms = (triallen*1000/fs);
binsizems = binsize*1000/fs;
lastBinms = binsizems * ceil(triallenms/binsizems);

edges = 0 : binsizems : lastBinms;
binscenters = edges(1:end-1)+(binsizems/2);

x = (mod(eventind-1,triallen)+1)*(1000/fs);
count = histc(x,edges);
r =  count / (ntrials*binsizems/1000);
r = r(1:end-1);
r = nanfastsmooth(r,wsize,2,0.5); 




