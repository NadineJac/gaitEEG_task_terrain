function naj_neurCorGait_task_template(sii, PATH, cfg)
% QUICK description
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
% - 
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-12-2022

% directories
PATHIN = PATH.rawdata;
PATHOUT = PATH.prepICA;

% load processing info
load(cfg.info_file);

for si = sii % loop through all subjects
    
    % sub ID
    ID = sprintf('sub-%03d',si);
    
    PATHOUTsi = fullfile(PATHOUT, ID);
    if ~exist(PATHOUTsi), mkdir(PATHOUTsi);end
    
    % filenames:
    % input file
    file_in = fullfile(PATHIN, ID,'eeg',[ID,'_task-neurCorrYoung_eeg.set']);
    % output file
    file_out = fullfile(PATHOUTsi,[ID,'_prepICA.set']); % file
    % report files
    reportNames = {'reportName1.png'};
    
    % check whether output file exists, if so, only run script if it has
    % been modified in the meantime
    %also check scripts that are being called!
    runScript = naj_scriptModified([mfilename('fullpath'),'.m'], file_out);
    if runScript
    disp('entering loop')
    % do tasks
    % save report figs
    
    % save data
    end
    fprintf('Done with sub-%03d!\n', si);
end % subject loop

save(cfg.info_file, 'procInfo');
disp(['Done with ', mfilename])
end %end function