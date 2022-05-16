function task_footprint(sii, PATH, cfg)
% Calculate gait-related artifact footprint with and without artifact
% attenuation. Footprint calculated with
% https://github.com/NadineJac/gaitEEGfootprint/tree/master/gaitEEGfootprint_lite
% and described in: Jacobsen NSJ, Blum S, Witt K, Debener S.
% A walk in the park? Characterizing gait-related artifacts in mobile EEG recordings.
% Eur J Neurosci. 2020;
%
% load file_in from PATHIN and store file_out in PATHOUT. Reports saved as
% reportNames
%
% INPUT
% - sii: vector of subject IDs to be processed
% - PATH:   structure with all paths leading to diff processing steps
% generated in naj_neurCorGait_paths
% - cfg:    structure with all processing variables generated in naj_neurCorGait_config
%
% overview of performed tasks:
% - do for data with and without artifact attenuation
%   - load data
%   - calculate gait-artifact realted footprint features B to F 
% - radar plot features with and without artifact attenuation
% - calculate euclidean distance between the footprint feature vectors
% with and without artifact attenuation
% - raincloud plot of distances
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-12-2022

% directories
PATHIN = PATH.sPCAsensor;
PATHOUT = PATH.footprint;

% load processing info
load(cfg.info_file);

% find all files
FILES = dir([PATHIN, filesep,'*sensor.mat']);
for ci = 1:length(FILES) % loop through all processing stages
    
    % filenames:
    % input file
    file_in = fullfile(PATHIN, FILES(ci).name);
    % output file
    file_out = fullfile(PATHOUT, [FILES(ci).name(1:end-4), '_footprint.mat']); % file
    reportName = {'sub-all_footprint_pre_post_cleaning';'sub-all_distances_pre_post_cleaning'};
    
    % check whether output file exists, if so, only run script if it has
    % been modified in the meantime
    runScript = naj_scriptModified([mfilename('fullpath'),'.m'], file_out);
    
    %% calculate footprint features
    if runScript
        
        % load data: use group averaged data already
        load(file_in);
        
        % calculate footprint features
        % gait ERSP.TF = times x chans x freqs --> TFdata = chans x freqs x pnts (HS to HS)
        if contains(file_in, 'raw') % use uncorrected data
            ERSP = permute(gaitERSP.ERSP_uncor(:,cfg.idx_EEG_chan,:,:), [2 3 1 4]);
        else % otherwise data after sPCA
            ERSP = permute(gaitERSP.ERSP(:,cfg.idx_EEG_chan,:,:), [2 3 1 4]);
        end
        
        for si = sii % subject loop
            ERSPsi = ERSP(:,:,:,si); % get subjects data
            
            % calculate features
            % skip feature A as no ERP after sPCA available!
            % B) correlation across frequencies --------------------------
            feature(si,1) = B_Rfreq(ERSPsi);
            
            % C) power ratio lateral/medial channels -------------------------
            feature(si,2) = C_lateralPowRatio(ERSPsi, cfg.lateralChanIdx);
            
            % D) power at neck electrodes contralateral to HS/ipsi --------------
            feature(si,3) = D_neckChanRatio(ERSPsi,cfg.neckChanL, cfg.neckChanR, cfg.pntsLHS, cfg.pntsRHS);
            
            % E) power double support/single supp gait cycle power -------------
            feature(si,4) = E_doubleSuppRatio(ERSPsi, cfg.pntsDouble);
            
            % F) S/W power ratio --------------------------------------
            feature(si,5) = F_swRatio(ERSPsi);
        end % subject loop
        
        ID = {procInfo(sii).ID};%get IDs from all subjects
        FOOTPRINT = table(feature(:,1),feature(:,2), feature(:,3), feature(:,4), feature(:,5),...
            'VariableNames', { 'B', 'C', 'D', 'E', 'F'},...
            'RowNames', ID');
        save(file_out, 'FOOTPRINT'); % save footprint features
        FEATURES{ci} = feature;
    end %run script
end % processing stage loop

%% radar plot
%  get spider plot from https://de.mathworks.com/matlabcentral/fileexchange/59561-spider_plot
allFootprint = cat(3,FEATURES{:}); % concatenate data of all subj
feature = squeeze(mean(allFootprint, 1));% average across subjects -> plot mean
axes_limits = repmat([floor(min(feature,[],'all'));ceil(max(feature,[],'all'))],1,5);
spider_plot(feature',...
    'AxesLabels', cfg.axes_labels,...
    'AxesInterval', cfg.axes_interval,...
    'AxesPrecision', cfg.axes_precision,...
    'AxesDisplay', cfg.axes_display,...
    'AxesLimits', axes_limits,...
    'FillOption', cfg.fill_option,...%
    'Color', cfg.edgeColors,...
    'LineWidth', cfg.line_width,...
    'Marker', cfg.marker_type,...
    'AxesFontSize', cfg.axes_font_size,...
    'LabelFontSize', cfg.label_font_size,...
    'AxesLabelsEdge', 'none');

% legend
lgd = legend({'pre', 'post'}, 'Fontsize', cfg.label_font_size);
title(lgd,'Artifact attenuation', 'Fontsize', cfg.label_font_size);
legend('boxoff');
set(lgd, 'Position',[0.45 0.059 0.18 0.06]);

set(gcf, 'units', 'centimeters', 'Position',[0 0 10 12])
print(fullfile(PATHOUT, reportName{1}),'-dpng'); %close;

%% euclidean diatances between footprint feature vectors
% for statistical evaluation
distFootprint = nan(length(sii),1);
for si = 1:length(distFootprint)
    distFootprint(si) = norm(allFootprint(si,:,2)-allFootprint(si,:,1)); % calculate distances
end
footprint = table(distFootprint, 'VariableNames', {'dist'}); 
writetable(footprint, [PATHOUT filesep 'distances']); % save to evaluate in JASP

%% Rain cloud plot of distances
% with RainCloudPlots from https://github.com/RainCloudPlots/RainCloudPlots
grey = [.6 .6 .6];
axes_font_size = 10;

figure; set(gcf, 'units', 'centimeters', 'Position',[0 0 10 5])
h = raincloud_plot(distFootprint, 'box_on',1,'box_dodge', 1, 'box_dodge_amount',...
    .3, 'dot_dodge_amount', .3, 'color', grey, 'cloud_edge_col', grey,'line_width', 1);
ylabel('density (a.u.)'), xlabel('distance (a.u.)');
xlim([0 10]);
box off;
set(gca, 'FontSize', axes_font_size);
print(fullfile(PATHOUT, reportName{2}),'-dpng'); %close;

disp(['Done with ', mfilename])

end %end function