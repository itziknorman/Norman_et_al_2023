
function [sig, pval, onsets, offsets] = compare_conditions(h,cond1,cond2,T,stat_method,plot_sig_mark)
% statistical test beween a pair of two conditions.
% h = axes handle
% cond1
% cond2
% T = timestamps
% stat_method:
% # 1=cluster-based paired;
% # 2=cluster-based unpaired
% # 3=wilcoxon signed-rank
% # 4=wilcoxon ranksum
if nargin<6, plot_sig_mark = 1; end
axes(h);
set_figure_colors;
sig = [];
pval = [];
onsets = [];
offsets = [];

if stat_method == 1
    % Two-sided paired-samples cluster-based permutation test using Groppe et al. 2011 implementation:
    tmp = cond1-cond2;
    tmp = reshape(tmp, 1, size(tmp,1), size(tmp,2));
    [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1(tmp,0,5000,0.05,0,0.05,2,[],0);
    sig = pval < 0.05; min_cluster_pval = min(pval);
    
elseif stat_method == 2
    % Two-sided unpaired-samples cluster-based permutation test using Groppe et al. 2011 implementation:
    tmp1 = reshape(cond1, 1, size(cond1,1), size(cond1,2));
    tmp2 = reshape(cond2, 1, size(cond2,1), size(cond2,2));
    [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm2(tmp1,tmp2,0,5000,0.05,0,0.05,2,[],0);
    
    sig = pval < 0.05; min_cluster_pval = min(pval);
    
elseif stat_method == 3
    % Wilcoxon signed rank-
    p = []; z = [];
    for i = 1:size(cond1,1)
        [p(i),~,temp] = signrank(cond1(i,:),cond2(i,:),'tail','both');
    end
    [~, crit_p, ~, pval]=fdr_bh(p,0.05,'pdep','yes');
    min_cluster_pval = 0.05;
    sig = pval < min_cluster_pval;
    
    
    
elseif stat_method == 4
    % Wilcoxon rank sum -
    p = []; z = [];
    for i = 1:size(cond1,1)
        [p(i),~,temp] = ranksum(cond1(i,:)',cond2(i,:)','tail','both');
    end
    [~, crit_p, ~, pval]=fdr_bh(p,0.05,'pdep','yes');
    sig = pval < 0.05;
end

if any(sig) & plot_sig_mark
    R_time = T;
    idx=sprintf('%d',sig);
    onsets = R_time(regexp(idx, '1{2,}', 'start'));    offsets = R_time(regexp(idx, '1{2,}', 'end'));
    tmp = get(gca,'Ylim');
    X=[onsets;offsets];
    Y=repmat(tmp(1)+range(tmp)*0.01,size(X));
    line(X,Y,'LineWidth',4,'Color',COLOR.sig)
    text(mean(X),mean(Y)+0.05,sprintf('p<%.4f',min_cluster_pval),'fontsize',6,'horizontalalignment','center');
end
end