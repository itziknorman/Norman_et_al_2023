function [maxClusterSize,meanSuffledRate,maxClusterActualdata,sig,MSD,clusterPval]=PSTH_permutation_stats(data1,data2,I,binsize,smt,Fs,clusterTh)
% Computes the shuffled peri-stimulus time histogram by randomly shifting ripple times.
% The function computes the baseline (mean) shuffled ripple rate and computes for each
% iteration the max cluster size that significantly differs from the baseline rate.
% the output is a vecotr of max cluster sizes, meanSuffledRate, realRate
% and the significant time bins (clusters)
%
% [maxClusterSize,meanSuffledRate,maxClusterActualdata,sig,MSD,clusterPval] = PSTH_permutation_stats(data1,data2,I,binsize,smt,Fs,clusterTh)
%
% data1 - ripple matrix (samples X trials). 1=ripple; 0=no ripple;
% data2 - ripple matrix for permutation (samples X trials). 1=ripple; 0=no ripple;
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
nbins=ceil(size(data2,1)/binsize);
%B=zeros(nbins,I); 
shuffledRate=zeros(nbins,I);
[~,~,realRate,~]=basicPSTH(find(data1),binsize,smt,Fs,size(data1,2),size(data1,1));
rng('shuffle');
nshift=randi([binsize,size(data2,1)],[I,size(data2,2)]);
for i=1:I
    permdata=nan(size(data2,1),size(data2,2));
    for k=1:size(data2,2), permdata(:,k)=circshift(data2(:,k),nshift(i,k)); end    
    [~,~,shuffledRate(:,i),~]=basicPSTH(find(permdata),binsize,smt,Fs,size(permdata,2),size(permdata,1));
    if mod(i,100)==0, fprintf(' %d ',i); end
end
% preventing edge effects:
shuffledRate(1,:)=nan;
shuffledRate(end,:)=nan;
% max squraed differences from baseline:

meanSuffledRate=nanmean(shuffledRate(:));
%MSD=sqrt(nanmean(nanmean(bsxfun(@minus,shuffledRate,meanSuffledRate).^2)));
MSD=sqrt(nanmean(nanvar(shuffledRate)));
TH=[meanSuffledRate-(clusterTh*MSD) meanSuffledRate+(clusterTh*MSD)];
TEMP=shuffledRate<TH(1)|shuffledRate>TH(2);
for i=1:size(TEMP,2)
    [str,fin] = regexp(num2str(TEMP(:,i))', '1{1,}', 'start','end');
    if ~isempty(str)
        l=(fin-str)+1;
        maxClusterSize(i)=max(l);
    else
        maxClusterSize(i)=0;
    end
end

TEMP=realRate<TH(1)|realRate>TH(2);
[str,fin] = regexp(num2str(TEMP)', '1{1,}', 'start','end');
l=(fin-str)+1; % cluster size
sig=zeros(size(TEMP));

clusterPval = [];
for k=1:length(l)
    clusterPval(k) = (sum(maxClusterSize>l(k))+1) / (I+1); % compute p-value from the permutations
    if clusterPval(k)<0.05, sig(str(k):fin(k))=1; end
end
maxClusterActualdata=max(l);

end
