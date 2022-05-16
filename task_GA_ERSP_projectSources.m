function task_GA_ERSP_projectSources(PATH, cfg)
% project group average ERSPs and GPMs onto cortex
% using brainstorm and previously computed inversion kernel
% load file_in from PATHIN and store file_out in PATHOUT. Reports saved as
% reportNames
%
% INPUT
% - PATH:   structure with all paths leading to diff processing steps
% generated in naj_neurCorGait_paths
% - cfg:    structure with all processing variables generated in naj_neurCorGait_config
%
% overview of performed tasks:
% - load templates to load sensor data into brainstorm
% - load templates to load source data into brainstorm
% - load precomputed inversion kernel
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-12-2022

PATHIN = PATH.sPCAsensor;

% start brainstorm to check results
run(fullfile(PATH.brainstorm, 'brainstorm.m')) %start brainstorm

% check whether protocol exists
ProtocolName = cfg.BSprotocolName;
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Unknown protocol: ' ProtocolName]);
    % create cfg.BSprotocolName
    % settings here: default anatomy, one channel file per subj
end

% Select the current procotol
gui_brainstorm('SetCurrentProtocol', iProtocol);

% load templates for adding data to BS database
load('cfg_BSsensorTF.mat');
load('cfg_BSsourceTF.mat');

% sub ID
ID = 'sub-all'; 

file_Inv = dir(fullfile(PATH.brainstorm_db, cfg.BSprotocolName, 'data', ID, '@default_study','*KERNEL*'));
Inv_kernel = load([file_Inv.folder filesep file_Inv.name]); % load

% update links in template
TFsource.DataFile = fullfile(ID, '@default_study', file_Inv.name);
TFsource.HeadModelFile = Inv_kernel.HeadModelFile;
TFsource_avg = TFsource;
TFsource_avg.Time = 0;

% load data 
load(fullfile(PATHIN, 'sub-all_sPCA_sensor.mat'));
for di = {'ERSP_uncor', 'GPM_uncor','ERSP', 'GPM'}%data loop
    % import sensor data, to navigate recordings
    TFsensor.TF = permute(mean(gaitERSP.(di{:}),4), [2, 1, 3]);% chan x pnts x freqs
    
    % add data to BS, reload vis BS GUI if doesn't show up
    TFsensor.Comment = di{:};
    add_bst_timefreq(ProtocolName, ID,di{:},TFsensor, 'data')
    
    % project cleaned TF from sensor to source space
    for fi = 1:size(TFsensor.TF, 3) % frequencies
        TFsource.TF(:,:,fi) = [mean(TFsensor.TF(:,:,fi),4)'*Inv_kernel.ImagingKernel.']';
    end
    % add data to BS, reload BS via GUI if doesn't show up
    TFsource.Comment = [di{:},'_cortex'];
    add_bst_timefreq(ProtocolName, ID,di{:},TFsource, 'results')
    
    % average over time
    if contains(di{:}, 'ERSP')
        TFsource_avg.TF = mean(TFsource.TF,2);
        % add data to BS, reload vis BS GUI if doesn't show up
        TFsource_avg.Comment = [di{:},'_cortex_avgTime'];
        add_bst_timefreq(ProtocolName, ID,di{:},TFsource_avg, 'results')
    end
end

disp(['Done with ', mfilename])
end %end function