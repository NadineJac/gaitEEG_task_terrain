%% naj_neurCorGait_master
% This is the master script calling all other scripts for the publication
% [final title and link here]
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-11-2022

%% Preperation

% data
% download data from doi:10.18112/openneuro.ds003039.v1.0.2
% - EEG from 19 young pb performing a gait task
%   - terrain: lawn + paved
%   - task: only walk (ST) walk + perform self-paced button -presses (DT)
% (additional data not analyzed: gait initiation with both legs, self-paced
% button pressing while standing)

% required toolboxes (used version in parenthesis, not tested for other
% versions)
% - EEGLAB (version 2020.0)
%	- "AMICA" v1.5.2
%	- "ICLabel" v1.3 
%	- "clean_rawdata" v2.3 
%	- "dipfit" v4.1
%	- "fullRankAveRef" v0.10 
%	- "postAmicaUtility" v2.1 
%   - zapline-plus from https://github.com/MariusKlug/zapline-plus (accessed Oct 2021)
%   - noise tools from http://audition.ens.fr/adc/NoiseTools/ (Accessed Sep
%   2021)
% - Brainstorm (version Jan 2020)
% - Rain cloud plots from https://github.com/RainCloudPlots/RainCloudPlots
% (accessed Feb-2020)
% - spider plot v17.8 from https://www.mathworks.com/matlabcentral/fileexchange/59561-spider_plot
% - gait artifact related footprint from https://github.com/NadineJac/gaitEEGfootprint/tree/master/gaitEEGfootprint_lite

%% config files
% set paths for input an all derivates
PATH = naj_neurCorGait_paths; % set for author's machines --> define for further machines
% set and store all (pre)processing variables in one .mat file
cfg = naj_neurCorGait_config(fullfile(PATH.config, 'naj_neurCorGait_cfg.mat'));
sii = cfg.subjects; % subjects to process

%% tasks
% task_processing_step are functions w/ inputs:
% - sii: vector of subject IDs to be processed
% - PATH:   structure with all paths leading to diff processing steps
% generated in naj_neurCorGait_paths
% - cfg:    structure with all processing variables generated in naj_neurCorGait_config
%
% each task
% - checks whether the script was modiefied after the files were generated
% - loops through subjects or just processes one
% - saves figures illustrating the processing
% naj_neurCorGait_task_template(sii, PATH, cfg) 

%% EEG preprocessing
% - downsample
% - filter
% - reject channels
% - ASR w/ custom BL
% - clean line noise (zaplineplus)
% - interpolate channels
% - CAR
task_prep_ICA(sii, PATH, cfg); 

% - w/interim data
%   - cut 1s epochs and reject (threshold, probability)
%    - run AMICA (rank reduced)
% - backproject to original data
% - fit dipoles (no thresh)
% - IC label
task_ICA_decomp(sii, PATH, cfg);

% - reject eye components
task_ICA_rej(sii, PATH, cfg);

% - standing BL
%   - extract, mark epochs exeeding threshold, compute average power
% - gait ERSP
%   - extract data of all gait cond, mark epochs exeeding threshold,
%   compute ERSP
%   - run sPCA on averaged ERSP, keep eigenvalues
%   - extract data of EACH gait cond, mark epochs exceeding threshold,
%   compute ERSP, correct EACH cond w/ previously stored eigenvalues
%   - plot (un)corrected ERSPs and GPMs
% - save timeseries data for standing BL + walk conds for source sPCA
task_sPCA_sensor(sii, PATH, cfg);
task_GA_ERSP(sii, PATH, cfg); % accumulate data across subjects

%% EEG analysis
%%% data quality %%%
% calculate gait ERSP and GPM without artifact attenuation (as comparison)
task_rawERSP_sensor(sii, PATH, cfg); %check!

% footprint %
% calculate gait-artifact realted footprint features B to F before and
% after artifact attenuation and calculate, plot them in a radar plot
% calculate the euclidean distance between the footprint feature vectors
% with and without artifact attenuation and plot them as a raincloud plot
task_footprint(sii, PATH, cfg);

% artifact attenuation %
% plot group average ERSP and GPM across all conditions at Cz, also add
% mean spectra and beta band topographies
plot_artifact_attenuation_lineplot(PATH, cfg);

%%% condition differences of GPM at Cz? %%%
% compare GPMs from 6 to 40 Hz at Cz with a 2 (terrain: even, uneven) x 
% 2(task: ST, DT) rmANOVA and follow up with post-hoc dependent samples
% t-tests
task_ERSP_stats(sii, PATH, cfg);

% plot condition specific group average GPMs at Cz from 4 to 60 Hz
% highlight FDR(ANOVA) and cluster corrected (t-test) differences 
% plot topographies of differences
plot_ERSP_stats(sii, PATH, cfg); 

% prepare Brainstorm protocol and conpute inverse kernel for projection of 
% GA ERSPs and GPMs
task_prepKernel_BS(PATH, cfg);

% project group average ERSPs and GPMs onto cortex
% using brainstorm and previously computed inversion kernel
task_GA_ERSP_projectSources(PATH);

%% behavioral analysis
%%% gait performance %%%
% calculate stride time and stride time variability, save as .csv
% check whether steps are "plausible" i.e. not too long, correct oder of events?
naj_neurCorGait_gait_perf(sii, PATH, cfg)