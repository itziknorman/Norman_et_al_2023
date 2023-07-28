function  [H,h,hcb] = plot_spectrogram(S,params_plot)
% plotting spectrogram based on eeglab
% S = spectrogram matrix (freq X time X trials)
% params_plot = structure with various plotting parameters;
%
% Author: Itzik Norman (2022)

H = figure('color','w','name',params_plot.name,'position',params_plot.fig_position); 
hold on;
if isempty(params_plot.clim), params_plot.clim = [min(S(:)), max(S(:))]; end
if ~isfield(params_plot,'stimulus_onset'), params_plot.stimulus_onset = []; end
if ~isfield(params_plot,'stimulus_dur'), params_plot.stimulus_dur = []; end
if ~isfield(params_plot,'xlabel'), params_plot.xlabel = 'Time (s)'; end
if ~isfield(params_plot,'ylabel'), params_plot.ylabel = 'Frequency (Hz)'; end

if params_plot.logscaleflag
    imagesclogy (params_plot.t_out, params_plot.f_out_log, mean(S(params_plot.f_log_ind,:,:),3));
    set(gca,'ytick',[6 10 30 120],'ylim',[min(params_plot.f_out), 50])
    caxis (params_plot.clim); colormap(params_plot.CM);
else
    imagesc(params_plot.t_out, params_plot.f_out, mean(S,3)); caxis (params_plot.clim); colormap(params_plot.CM);
end

h = gca;
axis tight; set(gca,'position',params_plot.position);
set(gca,'position',params_plot.position);
xlim(params_plot.xlim)

%plot(get(gca,'Xlim'),[6 6], '-k','LineWidth', 0.5);
%plot(get(gca,'Xlim'),[10 10], '-k','LineWidth', 0.5);

if ~isempty(params_plot.stimulus_onset), plot([params_plot.stimulus_onset params_plot.stimulus_onset], get(gca,'Ylim'), '-k','LineWidth', 2); end
if ~isempty(params_plot.stimulus_dur), plot([params_plot.stimulus_dur params_plot.stimulus_dur], get(gca,'Ylim'), '-k','LineWidth', 2); end

title(params_plot.title,'FontSize',8,'FontWeight','normal');
ylabel (params_plot.ylabel, 'fontsize', 8); xlabel(params_plot.xlabel);

tmpx=get(gca,'xlim'); 
tmpy=get(gca,'ylim');
text(tmpx(2)-range(tmpx)*0.05,tmpy(2)-range(tmpy)*0.1, ...
     sprintf('%s',params_plot.ch_label),'fontsize',6,'HorizontalAlignment','right')
enlarge_figure_and_move_axes_to_center(gcf,gca,1.25);

% add colorbar:
if params_plot.cbflag
    cb = axes('position', params_plot.cb_position);
    set(cb,'box','off','color',[1 1 1],'Xcolor','k','Ycolor','k')
    hcb = cbar(cb,1:size(params_plot.CM,1),params_plot.clim,3); 
    hcb.Title.FontWeight = 'normal';
    switch params_plot.normflag
        case 0,  hcb.Title.String = 'db';
        case 1,  hcb.Title.String = 'zscore';
    end
  
    hcb.Title.FontSize = 6;    
    set(hcb,'Xcolor','k','Ycolor','k');
    colormap(params_plot.CM);
end

