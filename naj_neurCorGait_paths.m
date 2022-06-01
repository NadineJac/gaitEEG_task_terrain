function PATH = naj_neurCorGait_paths()
% set all relevant folders, create them if they do not exist yet
% add toolbox paths, start eeeeglab w/o GUI
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-11-2022

if strcmp(getenv('COMPUTERNAME'), 'LAPTOP-796NEANK')
    PATH.main = 'E:\nadine\otto_projects\naj_neuralCorrGait';
    PATH.brainstorm = 'E:\nadine\MATLAB\brainstorm3';
    PATH.eeglab = 'E:\nadine\MATLAB\eeglab2020_0';
     PATH.brainstorm_db = 'E:\nadine\brainstorm_db'; % path to brainstorm database where protocol is prepared
    
elseif strcmp(getenv('COMPUTERNAME'), 'DESKTOP-5V16SQU')
    PATH.main = 'D:\PhD\otto_projects\naj_neuralCorrGait';
    PATH.brainstorm = 'D:\PhD\MATLAB\brainstorm3\';
    PATH.eeglab = 'D:\PhD\MATLAB\eeglab2020_0';
    PATH.brainstorm_db = 'D:\PhD\brainstorm_db'; % path to brainstorm database where protocol is prepared
else
    error(['Define study directories for ',getenv('COMPUTERNAME'), ' in ', mfilename()]);
end

PATH.code       = fullfile(PATH.main, 'code');
PATH.codebase   = fullfile(PATH.code, 'codebase');
PATH.sourcedata = fullfile(PATH.main, 'sourcedata');
PATH.rawdata    = fullfile(PATH.main); 

PATH.config    = fullfile(PATH.code); 

% neurophsysiological derivates
PATH.derivates  = fullfile(PATH.main,'derivates');
PATH.prepICA    = fullfile(PATH.derivates, 'prep_ICA');
PATH.ICAdecomp  = fullfile(PATH.derivates, 'ICA_decomp');
PATH.ICAclean   = fullfile(PATH.derivates, 'ICA_rej');
PATH.sPCAsensor = fullfile(PATH.derivates, 'sPCAsensor');
PATH.source     = fullfile(PATH.derivates, 'source');
PATH.footprint  = fullfile(PATH.derivates, 'footprint');
PATH.study      = fullfile(PATH.derivates, 'study');

% behavioral derivates
PATH.behave     = fullfile(PATH.derivates, 'behave');

% results & figures
PATH.plotCleaning = fullfile(PATH.derivates, 'figures', 'artifactAttenuation');
PATH.plotConds = fullfile(PATH.derivates, 'figures', 'ERSP_conds');

addpath(PATH.eeglab); addpath(PATH.brainstorm)
addpath(genpath(PATH.code)); % add w/ subfolders
cd(PATH.code);
createFolders(PATH);

eeglab nogui