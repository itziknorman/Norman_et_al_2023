function  [r] = basicPSTHforBS(eventind, binsize, wsize, fs, ntrials, triallen)
% PSTH Computes the peri-stimulus time histogram from spike times.
% The routine plots the trial averaged spike rate as a function of time.
% r = PSTH(TIMES, BINSIZE, FS,NTRIALS,TRIALLEN)
% TIMES - spike times (samples)
% BINSIZE - binwidth (samples)
% wsize - size of a triangular smoothing window (bins).
% FS - sampling rate (hz)
% NTRIALS - number of trials
% TRIALLEN - length of a trial (samples)
%
% Author: Itzik Norman
% normanik@gmail.com
% @ Weizmann Institute, Israel
narginchk(6,6);

%Compute PSTH
triallenms = (triallen*1000/fs);
binsizems = binsize*1000/fs;
lastBinms = binsizems * ceil(triallenms/binsizems);

edges = 0 : binsizems : lastBinms;
binscenters = edges(1:end-1)+binsizems/2;

x = (mod(eventind-1,triallen)+1)*(1000/fs);
count = histc(x,edges);
r =  count / (ntrials*(binsizems/1000));
r = r(1:end-1);
r = nanfastsmooth(r,wsize,2,0.5); 

