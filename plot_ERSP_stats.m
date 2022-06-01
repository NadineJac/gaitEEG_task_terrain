function plot_ERSP_stats(sii, PATH, cfg)
% plot condition specific group average GPMs at Cz from 4 to 60 Hz
% highlight FDR(ANOVA) and cluster corrected (t-test) differences 
% load file_in from PATHIN and store file_out in PATHOUT. Reports saved as
% reportNames
% 
% INPUT
% - sii: vector of subject IDs to be processed
% - PATH:   structure with all paths leading to diff processing steps
% generated in naj_neurCorGait_paths
% - cfg:    structure with all processing variables generated in naj_neurCorGait_config
%
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-12-2022

% directories
PATHIN = PATH.sPCAsensor;
PATHOUT = PATH.plotConds;

% stats
file_in = {'2x2rmANOVA_Cz.mat',...
    'depT_ME_surface.mat',...
    'depT_ME_task.mat'};

% reports
reportNames = {'2x2rmANOVA_Cz.png',...
    'depT_ME_surface.png',...
    'depT_ME_task.png'};

%% ROI
chanIdx = cfg.COI; % channel to test: Cz
% look into all time points;

% plot params
clim = 0.75; % color axis for GPM
colorbr = 0;

%% load and prepare data
for ci = 1:length(cfg.cond_names)
    load(fullfile(PATHIN,['sub-all_sPCA_sensor_', cfg.cond_names{ci},'.mat']), 'gaitERSP');
    GPM{ci} = gaitERSP.GPM(:,:,:,sii);
    %1 DT easy, 2 ST easy, 3 DTdifficult, 4 ST difficult
end
DTeasy = 1; STeasy = 2; DTdifficult = 3; STdifficult = 4;

%%% main effects
% surface:
% even = mean(STeasy+ DTeasy)
GPM_even = mean(cat(5,GPM{STeasy},GPM{DTeasy}),5);
GPM_uneven = mean(cat(5,GPM{STdifficult},GPM{DTdifficult}),5);

contrasts.ME_surface = GPM_even-GPM_uneven;
contrasts.data.even = GPM_even;
contrasts.data.uneven = GPM_uneven;

% task
GPM_ST = mean(cat(5,GPM{STeasy},GPM{STdifficult}),5);
GPM_DT = mean(cat(5,GPM{DTeasy},GPM{DTdifficult}),5);

contrasts.ME_task = GPM_ST-GPM_DT;
contrasts.data.ST = GPM_ST;
contrasts.data.DT = GPM_DT;

% save contrasts
save(fullfile(PATHIN,['sub-all_contrasts.mat']), 'contrasts');

%% load stats
for fi = 1:length(file_in)
    load(fullfile(PATHIN,file_in{fi}));
end

%% plot GPMs
figure, set(gcf, 'units', 'centimeters', 'position', [0 0 20 20]);
c=6; %number of colums of subplot
% STeven = 2
p=1; subplot(6,c,[p,p+1,p+c,p+c+1]); ci = STeasy; hold on;
plot_gaitERSP(mean(GPM{ci},4), chanIdx, freqs, clim, colorbr, cfg); axis square
%title(cfg.cond_names{ci})

% DTeven = 1
p=3; subplot(6,c,[p,p+1,p+c,p+c+1]); ci = DTeasy; hold on;
plot_gaitERSP(mean(GPM{ci},4), chanIdx, freqs, clim, colorbr, cfg);axis square
%title(cfg.cond_names{ci})

% STuneven = 4
p=13; subplot(6,c,[p,p+1,p+c,p+c+1]); ci = STdifficult; hold on;
plot_gaitERSP(mean(GPM{ci},4), chanIdx, freqs, clim, colorbr, cfg); hold on;axis square
%title(cfg.cond_names{ci})

% DTuneven = 3
p=15; subplot(6,c,[p,p+1,p+c,p+c+1]); ci = DTdifficult; hold on
plot_gaitERSP(mean(GPM{ci},4), chanIdx, freqs, clim, colorbr, cfg);axis square
%title(cfg.cond_names{ci})

%% outline effects in GPMs
%terrain
p=11; subplot(6,c,[p,p+1,p+c,p+c+1]); hold on
plot_gaitERSP(mean(contrasts.ME_surface,4), chanIdx, freqs, clim, colorbr, cfg);axis square
try l{1} = contour(cfg.t_axis, freqs, ANOVA.results.pinter{1}<.05,1,...
    ':k','linewidth',1);end % ANOVA
try l{2}= contour(cfg.t_axis, freqs, ME_surface.results.pinter{1}<.05,1,...
    'k','linewidth',.5);end %t test
%title('even - uneven')

% task
p=26; subplot(6,c,[p,p+1,p+c,p+c+1]); hold on
plot_gaitERSP(mean(contrasts.ME_task,4), chanIdx, freqs, clim, colorbr, cfg);axis square
try contour(cfg.t_axis, freqs, ANOVA.results.pinter{2}<.05,1,...
    ':k','linewidth',1); end % ANOVA
try contour(cfg.t_axis, freqs, ME_task.results.pinter{2}<.05,1,...
    'k','linewidth',.5);end%t test
%title('ST - DT')

% interaction
p=29; subplot(6,c,[p,p+1,p+c,p+c+1]); hold on
plot_gaitERSP(mean(contrasts.ME_surface-contrasts.ME_task,4), chanIdx, freqs, clim, colorbr, cfg);axis square
try contour(cfg.t_axis, freqs, ANOVA.results.pinter{3}<.05,1,...
    'k','linewidth',1);end % ANOVA
%title([{'(even - uneven) -'};{'(ST - DT)'}])
xlabel('Gait cycle (%)'); ylabel('Frequency (Hz)')

c = colorbar(subplot(6,c,[p,p+1,p+c,p+c+1]),'Position', [0.88 0.38 0.018 0.10]);
ylabel(c, 'change to mean gait cycle')
set(findall(gcf,'-property','FontSize'),'FontSize',8)
%sgtitle('Grand average gait-phase related power modulations at Cz across conditions')
print(fullfile(PATHOUT,reportNames{1}), '-dpng'); %close

%% cluster topography
figure; set(gcf, 'units', 'centimeters', 'position', [0 0 4 4]); 
%title('topo cluster');
data = mean(contrasts.ME_surface(:,:,cfg.FOI,:),4);% avg over subjects --> pnts x chan x FOI
idxPixel = ME_surface.results.pinter{1}<.05; % idx of pixels that are part of cluster
idx = permute(repmat(idxPixel, [1,1,size(data,2)]), [2,3,1]);
data(~idx) = nan; % only keep cluster data
data = squeeze(mean(data,[1,3], 'omitnan')); % average over time and freqs --> channel vector
topoplot(data, gaitERSP.chanlocs,...
    'maplimits', [-clim clim], 'emarker2', {chanIdx, '*', 'k'});
try colormap(cmap_mandrill); catch  colormap jet; end
end
