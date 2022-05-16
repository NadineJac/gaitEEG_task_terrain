function plot_allChanERSP(ERSP, plotID, chanlocs, f_axis, t_axis, clim, stringTitle, unit, visible)% EV)
% plot_allChanERSP in layout, similar to BS
% ERSP = pnts, chan, freq

% lab = EV{1}; % lab = params.allGaitEV; names of all gait events
% lab = EV{2}; % lat = params.newLatms; latency of all warped to gait events in ms, handed over in one cell

% make a few more variables accessible?
freqs = f_axis;
times = t_axis;

f_ticks = 10:10:f_axis(end);
t_ticks = linspace(times(1), times(end), 11);
t_tickLab = 0:10:100; %labels in percent

fig = figure('visible',visible); set(fig,'units', 'centimeters','Position', [0 0 40 25])
for ch = 1:height(plotID) % plot each chan seperatly
    
    idxChan = find(strcmp(plotID.ID, chanlocs(ch).labels));
    subplot(10,11,plotID.num(idxChan))
    
    data = squeeze(ERSP(:,ch,:));
    contourf(times, freqs, data', 50,'linecolor','none'); hold on % plot ERSP
%     if nargin(8) %plot lines for events
%         plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
%     end
        
    caxis([-clim clim]); %        colorbar
    colormap jet, box off
    
    title(chanlocs(ch).labels)
    xticks(t_ticks); yticks(f_ticks); xticklabels(t_tickLab);
end

subplot(10,11,110)
data = squeeze(mean(ERSP,2));
contourf(times, freqs,data', 50,'linecolor','none'); hold on % plot ERSP

%  if nargin(8) %plot lines for events
%         plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
%     end

caxis([-clim clim]);
colormap  jet, box off

title('Channel mean')
xticks(t_ticks); yticks(f_ticks); xticklabels(t_tickLab);

h = colorbar('Position',[0.92 0.54 0.01 0.05]);
ylabel(h, unit)

han=axes(fig,'visible','off');
han.XLabel.Visible='on';xlabel(han,'Gait cycle (%)');
han.YLabel.Visible='on';ylabel(han,'Frequency (Hz)');

set(findall(gcf,'-property','FontSize'),'FontSize',6)

sgtitle(stringTitle,'interp', 'none')
end

