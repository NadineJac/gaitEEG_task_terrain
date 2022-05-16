function task_GA_ERSP(sii, PATH, cfg)
% Calculate and save group average gait ERSPs and GPMs from subjecsts (sii)
% across all gait conditions together and for each condition (cfg.cond_walk) alone
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
% - for all conditions together and condition-specific
%   - Load and average subject-specific ERSPs and GPMs
%   - save group average
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-12-2022

% directories
PATHIN = PATH.sPCAsensor;
PATHOUT = PATH.plotCleaning;

% load processing info
load(cfg.info_file);

for si = sii % loop through all subjects
    
    % sub ID
    ID = sprintf('sub-%03d',si);
    
    % filenames:
    % input file
    file_in = fullfile(PATHIN,ID, [ID,'_sPCA_sensor.mat']);
    % output file
    file_out = fullfile(PATHIN,['sub-all_sPCA_sensor.mat']); % file
    % report files
    reportNames = {'reportName1.png'};
    
    % check whether output file exists, if so, only run script if it has
    % been modified in the meantime
    runScript = naj_scriptModified([mfilename('fullpath'),'.m'], file_out);
    if runScript
        load(file_in)
        TMP.ERSP(:,:,:,si) = gaitERSP.ERSP;
        TMP.GPM(:,:,:,si) = gaitERSP.GPM;
        TMP.ERSP_uncor(:,:,:,si) = gaitERSP.ERSP_uncor;
        TMP.GPM_uncor(:,:,:,si) = gaitERSP.GPM_uncor;
        TMP.chanlocs = gaitERSP.chanlocs;
        
        for ci = 1:length(cfg.cond_walk) % condition loop
            load([file_in(1:end-4), '_', cfg.cond_names{ci}, '.mat'], 'gaitERSP');
            TMP_cond{ci}.ERSP(:,:,:,si) = gaitERSP.ERSP;
            TMP_cond{ci}.GPM(:,:,:,si)  = gaitERSP.GPM;
            TMP_cond{ci}.chanlocs       = gaitERSP.chanlocs;
            
            if si == sii(end) % condition-specific group averages in last iteration
             gaitERSP = TMP_cond{ci};   
             save([file_out(1:end-4), '_', cfg.cond_names{ci}, '.mat'], 'gaitERSP');
            end
        end
    end
    fprintf('Done with sub-%03d!\n', si);
end % subject loop

% save group averages (all condition together) after last subject is
% processed
gaitERSP = TMP;
save(file_out, 'gaitERSP');
save(cfg.info_file, 'procInfo');

disp(['Done with ', mfilename])
end % function
