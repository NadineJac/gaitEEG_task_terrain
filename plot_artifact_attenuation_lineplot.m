function plot_artifact_attenuation_lineplot(PATH, cfg)
% plot group average ERSP and GPM across all conditions at Cz, also add
% mean spectra and topographies
%
% load file_in from PATHIN and store file_out in PATHOUT. Reports saved as
% reportNames
% 
% INPUT
% - PATH:   structure with all paths leading to diff processing steps
% generated in naj_neurCorGait_paths
% - cfg:    structure with all processing variables generated in naj_neurCorGait_config
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-12-2022

% directories
PATHIN = PATH.sPCAsensor;
PATHOUT = PATH.plotCleaning;

reportNames = {'artifactAttenuation'};

%% load group averaged data (all conditions together)
% without artifact attenuation
load(fullfile(PATHIN,['sub-all_rawERSP_sensor.mat']), 'gaitERSP');
ERSP_uncor = gaitERSP.ERSP_uncor;
GPM_uncor = gaitERSP.GPM_uncor;
% with artfact attenuation
load(fullfile(PATHIN,['sub-all_sPCA_sensor.mat']), 'gaitERSP');
ERSP = gaitERSP.ERSP;
GPM = gaitERSP.GPM;

%% plot
% parameters
f_axis  = cfg.f_axis; %frequency axis
chanIdx = 24;       % channel to plot --> here: 24 = Cz
clim1   = 2.25;     % colorbar limits for ERSP
clim2   = 0.75;     % colorbar limits for GPM
colorbr = 0;        % flag for plotting colorbar -> here: 0 = no
beta    = 20:30;    % frequecies for which to plot topographies -> here: beta from 20-30Hz
faceAlpha = .3;     % transparency to indicate frequencies for topography in spectrum

figure, set(gcf, 'units', 'centimeters', 'position', [0 0 27.5 15]);
row= 4; col = 7;

%% ERSP
p=1; subplot(row, col,[p,p+1, p+col, p+col+1])% ERSP without artifact attenuation
plot_gaitERSP(mean(ERSP_uncor,4), chanIdx, f_axis,clim1, colorbr, cfg);
line([1 1; 100 100], [beta(1), beta(end); beta(1), beta(end)],'color', [.5 .5 .5]);

p=3; subplot(row, col,[p,p+1, p+col, p+col+1]) % ERSP with artifact attenuation
plot_gaitERSP(mean(ERSP,4), chanIdx, f_axis, clim1, colorbr, cfg);
line([1 1; 100 100], [beta(1), beta(end); beta(1), beta(end)],'color', [.5 .5 .5]);
c = colorbar(subplot(row, col,[p,p+1, p+col, p+col+1]),'Position', [0.93 0.78 0.014 0.13]);
ylabel(c, 'change to standing BL                           ')
title(c, 'dB')

%% GPM
p=2*col+1; subplot(row, col,[p,p+1, p+col, p+col+1])% GPM without artifact attenuation
plot_gaitERSP(mean(GPM_uncor,4), chanIdx, f_axis,clim2, colorbr, cfg);
line([1 1; 100 100], [beta(1), beta(end); beta(1), beta(end)],'color', [.5 .5 .5]);
xlabel('Gait cycle (%)'), ylabel('Frequency (Hz)')

p=2*col+3; subplot(row, col,[p,p+1, p+col, p+col+1])% GPM with artifact attenuation
plot_gaitERSP(mean(GPM,4), chanIdx, f_axis, clim2, colorbr, cfg);
line([1 1; 100 100], [beta(1), beta(end); beta(1), beta(end)],'color', [.5 .5 .5]);
c = colorbar(subplot(row, col,[p,p+1, p+col, p+col+1]),'Position', [0.93 0.35 0.014 0.13]);
ylabel(c, 'change to mean gait cycle                      ')
title(c, 'dB')

%% spectra
p=5; subplot(row, col,[p,p+col]), hold on
% highlight frequencies of beta topo
patch([-2 -2 6 6], [beta(1), beta(end), beta(end), beta(1)], ...
    [.5 .5 .5], 'FaceAlpha', faceAlpha, 'EdgeColor', 'none')
% mean spectra of ERSP without artifact attenuation
plot(squeeze(mean(ERSP_uncor(:,chanIdx,:,:),[1,2,4])),cfg.f_axis, ...
    'color', cfg.orange, 'linewidth', 1.5);
% mean spectra ERSP with artifact attenuation
plot(squeeze(mean(ERSP(:,chanIdx,:,:),[1,2,4])),cfg.f_axis,...
    'color', cfg.blue, 'linewidth', 1.5);
line([0 0], [f_axis(1), f_axis(end)], 'color', [.5 .5 .5])
ylim([f_axis(1), f_axis(end)]);
grid on, box off
% legend({'without', 'with'}), legend box off
% xlabel('power')

p=2*col+5; subplot(row, col,[p,p+col]); hold on
% highlight frequencies of beta topo
patch([0 0 1 1], [beta(1), beta(end), beta(end), beta(1)], ...
    [.5 .5 .5], 'FaceAlpha', faceAlpha, 'EdgeColor', 'none')
% mean absolute spectra of GPM without artifact attenuation
plot(squeeze(mean(abs(GPM_uncor(:,chanIdx,:,:)),[1,2,4])),cfg.f_axis, ...
    'color', cfg.orange, 'linewidth', 1.5);
% mean absolute spectra of GPM with artifact attenuation
plot(squeeze(mean(abs(GPM(:,chanIdx,:,:)),[1,2,4])),cfg.f_axis,...
    'color', cfg.blue, 'linewidth', 1.5);
line([0 0], [f_axis(1), f_axis(end)], 'color', [.5 .5 .5])
ylim([f_axis(1), f_axis(end)]);
grid on, box off
xlabel('power (dB)')

%% topoplots: pnts x chan x freq x subj
% ERSP
p=col-1; subplot(row, col,p) % without artifact attenuation
tmp = squeeze(mean(ERSP_uncor(:,:,ismember(cfg.f_axis,beta),:),[1,3,4]));
topoplot(tmp, gaitERSP.chanlocs, 'maplimits', [-clim1 clim1],...
    'emarker2', {chanIdx, '*', 'k'});

p=col; subplot(row, col,p)% with artifact attenuation
tmp = squeeze(mean(ERSP(:,:,ismember(cfg.f_axis,beta),:),[1,3,4]));
topoplot(tmp, gaitERSP.chanlocs, 'maplimits', [-clim1 clim1],...
    'emarker2', {chanIdx, '*', 'k'});

% GPM
p=3*col-1; subplot(row, col,p)% without artifact attenuation
tmp = squeeze(mean(abs(GPM_uncor(:,:,ismember(cfg.f_axis,beta),:)),[1,3,4]));
topoplot(tmp, gaitERSP.chanlocs, 'maplimits', [-clim2 clim2],...
    'emarker2', {chanIdx, '*', 'k'});

p= 3*col; subplot(row, col,p)% with artifact attenuation
tmp = squeeze(mean(abs(GPM(:,:,ismember(cfg.f_axis,beta),:)),[1,3,4]));
topoplot(tmp, gaitERSP.chanlocs, 'maplimits', [-clim2 clim2],...
    'emarker2', {chanIdx, '*', 'k'});

try colormap(cmap_mandrill); catch  colormap jet; end %if somewhere in path: set colormap to mandill (from brainstorm)

set(findall(gcf,'-property','FontSize'),'FontSize',10)

print(fullfile(PATHOUT, reportNames{1}), '-dpng');
print(fullfile(PATHOUT, reportNames{1}), '-deps');%close;
end
