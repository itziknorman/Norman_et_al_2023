function [maxClusterSize,meanSuffledRateDiff,maxClusterActualdata,sig,MSD,clusterPval]=PSTH_permutation_stats_shuffled_labels(data1,data2,I,binsize,smt,Fs,clusterTh)
% Computes the shuffled peri-ripple time histogram by randomly shifting ripple times.
% The function compute the baseline (mean) difference when randomly dividing the data into two groups, and compute for each
% iteration the max cluster size that significantly differ from 0 (e.g. at
% p<0.01). Then the clusters observed in the actual data are compared to
% the distribution of random cluster. 
% the output is a vecotr of max cluster sizes, meanShuffledRateDiff,realRateDiff, 
% significant temporal clusters, Mean Squared Diff of the shuffeld data and clusters P-values 
%
% [maxClusterSize,meanSuffledRateDiff,maxClusterActualdata,sig,MSD,clusterPval]=PSTH_permutation_stats_shuffled_labels(data1,data2,I,binsize,smt,Fs,clusterTh)
%
% data1 - spike matrix (samples X trials). 1=spike; 0=no spike;
% data2 - spike matrix for permutation (samples X trials). 1=spike; 0=no spike;
% I - num of iterations
% binsize - size of bins (ms).
% smt - smoothing window (number of bins)
% Fs - sampling rate (Hz)
% clusterTh - stdev
%
% Author: Itzik Norman
% normanik@gmail.com
% Weizmann Institute of Science, Israel

if isempty(clusterTh), clusterTh=1.96; end % one sided ztest at p<0.05: clusterTh=1.65;

fprintf('\n PERMUTING DATA (%d iterations) \n',I);
assert(size(data1,1)==size(data2,1));

[~,binscenter,realRate1,~]=basicPSTH(find(data1),binsize,smt,Fs,size(data1,2),size(data1,1));
[~,binscenter,realRate2,~]=basicPSTH(find(data2),binsize,smt,Fs,size(data2,2),size(data2,1));
ActualRateDiff = realRate1-realRate2;

nbins=length(ActualRateDiff);
shuffledRate1=zeros(nbins,I);
shuffledRate2=zeros(nbins,I);
% Shuffled PSTHs produced by circularly jittering ripple timing in each
% trial by a random amount: (see a similar method in Rothschild et al, NN 2016)
data=[data1 data2];
N1=size(data1,2); N2=size(data2,2);
rng('shuffle');
for i=1:I
    permdata=data(:,randperm(size(data,2)));
    [~,~,shuffledRate1(:,i),~]=basicPSTH(find(permdata(:,1:N1)),binsize,smt,Fs,N1,size(permdata,1));
    [~,~,shuffledRate2(:,i),~]=basicPSTH(find(permdata(:,N1+1:end)),binsize,smt,Fs,N2,size(permdata,1));
    if mod(i,100)==0, fprintf(' %d ',i); end
end

% preventing edge effects:
shuffledRate1(1,:)=nan;
shuffledRate1(end,:)=nan;
shuffledRate2(1,:)=nan;
shuffledRate2(end,:)=nan;
% max squraed differences under H0 should be zero (no difference between conditions):
suffledRateDiff = bsxfun(@minus,shuffledRate1,shuffledRate2);
meanSuffledRateDiff=nanmean(suffledRateDiff(:)); % should be zero
assert(abs(meanSuffledRateDiff)<0.01,'Mean shuffled differences should be closer to zero... add more permutations')
MSD=sqrt(nanmean(nanmean(bsxfun(@minus,suffledRateDiff,meanSuffledRateDiff).^2)));
TH=[meanSuffledRateDiff-(clusterTh*MSD) meanSuffledRateDiff+(clusterTh*MSD)];
TEMP=suffledRateDiff<TH(1)|suffledRateDiff>TH(2);
for i=1:size(TEMP,2)
    [str,fin] = regexp(num2str(TEMP(:,i))', '1{1,}', 'start','end');
    if ~isempty(str)
        l=(fin-str)+1;
        maxClusterSize(i)=max(l);
    else
        maxClusterSize(i)=0;
    end
end

TEMP=ActualRateDiff<TH(1)|ActualRateDiff>TH(2);
[str,fin] = regexp(num2str(TEMP)', '1{1,}', 'start','end');
l=(fin-str)+1; % cluster size
sig=zeros(size(TEMP));
fprintf('\n %d cluster were found in the actual data \n',numel(str))

clusterPval = [];
for k=1:length(l)
    clusterPval(k) = (sum(maxClusterSize>l(k))+1) / (I+1); % compute p-value from the permutations
    if clusterPval(k)<0.05, sig(str(k):fin(k))=1; end
end
maxClusterActualdata=max(l);

end
