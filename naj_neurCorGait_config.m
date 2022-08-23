function cfg = naj_neurCorGait_config(file_out)%% Config file
% Configuration parameters for the study.
% These are all the relevant parameters for the analysis organized
% according to the script they are used in
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-11-2022

subjects = 1:19;
tmp = dir(file_out);
info_file = fullfile(tmp.folder, 'naj_neurCorGait_info.mat');
naj_neurCorGait_info(subjects,info_file);

%% task_prepICA ____________________________
% prune data
startExperiment = {'start_restEEG'};

% sampling rate: downsample EEG data
sampling_rate = 250;

% channel indices of diff sensors
idx_EEG_chan = 1:64;
idx_headAcc_chan = 65:67;

% Band-pass filter limits.
bandpass_fmin = 0.2;  % Hz passband
bandpass_fmax = 60;  % Hz passband --> cut-off 67.5 Hz

% ASR parameters
ASR_baseline_events = {'start_standing', 'end_standing'};
% otherwise default, still specify here!
ASR_cutoff = 20;
ASR_FlatlineCriterion = 5;
ASR_ChannelCriterion = 0.8;
ASR_LineNoiseCriterion = 4;

% zap line plus
zapLinePlus.noisefreqs = 50; % line noise freq
zapLinePlus.maxfreq = 61; %highest freq consinderes as noise

%% task_ICA_decompose ________________________________
% Since we are performing ICA on the continuous data,
% it is important that the lower bound is at least 1Hz.
% keep in mind that eeglab does not take the cut-off but the passband as
% input!
% params for processing of interim dataset (EEGTMP)
ICA_bandpass_fmin = 2;  % 1 Hz cut-off

% rejection params
ICA_REJsd = 3; % strd for rejection of improbale epoch
ICA_REJthresh = 350; %muV amplitude for epoch rejection

% AMICA params
% ICA parameters
AMICA_numprocs    = 1;    % # of nodes (default = 1)
AMICA_max_threads = 2;    % # of threads per node
AMICA_num_models  = 1;    % # of models of mixture ICA
AMICA_max_iter    = 2000; % max number of learning steps

%% task_ICArej
%IC label thresholds
IC_eye = .9;
IC_muscle = .9;

%% spectral denoising
V_tsh = 300; % threshold for marking "bad" data sections to be excluded, as Seeber in his scripts
gait_event = 'RightHS'; % event for extracting gait epochs and warping
gait_event_order = {'RightHS', 'LeftTO', 'LeftHS', 'RightTO', 'RightHS'};% order of gait events within each gait cycle
gait_timeNextHs = [.5 1.5]; % time of next RHS in s

% wavelet parameters for time-frequency decomposition
f_axis = 2:2:60; % data filtered at 60 so going up until 80 does not make sense!
t_axis = 1:100;
FWHM = log2(f_axis);
N_freq = length(f_axis);
N_ds = 8;

% events for extracting data
cond_stand ={'start_standing', 'end_standing'};
cond_walk = {'start_easy_button', 'end_easy_button';...
    'start_easy', 'end_easy';...
    'start_difficult_button', 'end_difficult_button'
    'start_difficult', 'end_difficult'};
f_surface = {'even', 'even', 'uneven', 'uneven'};
f_task = {'DT', 'ST', 'DT', 'ST'};
cond_names = {'DTeasy', 'STeasy', 'DTdifficult','STdifficult'};



% gait ERPs
% gait_epoch = [0 1.5]; % length of epoch in s
% gait_timeNextEV = [25 1500]; %delay (ms) ind which next gait events can occur

%% footprint ________________________________________
gait_event_newLat   = [1 18 50 68 100]; % % new latencies in pnts
pntsRHS             = gait_event_newLat(1):gait_event_newLat(2);          % double support following right-heel strike
pntsLHS             = gait_event_newLat(3):gait_event_newLat(4);         % double support following left-heel strike
pntsDouble          = [pntsRHS,pntsLHS];

% channel indices (dim 1 of the time-frequency decomposed data)
% [ADAPT] channel selection basel on you layout (here: custom 64ch layout)
lateralChanIdx = [1,4,5,6,9,10,15,16,17,18,20,21,26,27,30,31,32,33,35,37,41,44,45,46,48,49,51,52,53,57,61,64]; % index of channels labelled as lateral
neckChanR      = [49 52 18 51];  % index of channels located over the right side of the neck
neckChanL      = [45 48 46 16];  % index of channels located over the left side of the neck

% task_head_acc
idx_head_acc_chan = 65:67;

% radar plot parameters
axes_interval = 4;% Axes properties
axes_precision = 2;
axes_display = 'one';
marker_type = 'none';
axes_font_size = 10;
label_font_size = 12;
axes_labels = {'B','C','D','E','F'}; % Axes labels
fill_option = 'on';
fill_transparency = 1;
orange = [.9 .6 0];
lightorange = [251 240 217]/255;
blue = [0 .45 .7];
lightblue = [217 234 244]/255;
edgeColors = [orange;blue];
colors = [lightorange;lightblue];
line_width = 1;
line_style ={'--','-'};
axes_limits = repmat([-1; 2],1,6);

%% source imaging
BSprotocolName = 'naj_neurCorrGait';%BS file
load('cfg_BSsensorTF.mat');% load template
TFsensor.Comment = ['Avg. Power', num2str(f_axis(1)), '-',  num2str(f_axis(end)),'Hz'];
TFsensor.Freqs = f_axis;
TFsensor.Options.MorletFwhmTc = FWHM;
TFsensor.Time = linspace(0,100,100);
TFsensor.TFmask = ones(length(f_axis), length(TFsensor.Time));
save('cfg_BSsensorTF.mat', 'TFsensor');

% source
load('cfg_BSsourceTF.mat');
TFsource.Time = linspace(0, 100, 100);
TFsource.Freqs = f_axis;
TFsource.TFmask = ones(length(TFsource.RowNames),length(TFsource.Freqs));
save('cfg_BSsourceTF.mat', 'TFsource');

%% stats
freqs = 6:2:40; % frequencies to test in Hz
FOI = ismember(f_axis, freqs);% their indices
COI = 24; % channel to test: Cz

stats.depT_surface = load('cfg_depT_surface.mat');
stats.depT_task = load('cfg_depT_task.mat');
stats.depT_interaction = load('cfg_depT_interaction.mat');

%% save
save(file_out);
cfg = load(file_out);

end
