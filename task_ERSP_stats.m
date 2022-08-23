function task_ERSP_stats(sii, PATH, cfg)
% compare GPMs from 6 to 40 Hz at Cz with a 2 (terrain: even, uneven) x 
% 2(task: ST, DT) rmANOVA and follow up with post-hoc dependent samples
% t-tests
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
% v1.1 last changed Aug-23-2022: deleted ANOVA, included dep T-test for
% interaction

% directories
PATHIN = PATH.sPCAsensor;
PATHOUT = PATH.plotConds;

% directories
PATHIN = PATH.sPCAsensor;
PATHOUT = PATH.sPCAsensor;

% reports
file_out = {'depT_ME_surface.mat',...
    'depT_ME_task.mat',...
    'depT_interaction.mat'};

%% channel and frequency ROI
chanIdx = cfg.COI; % channel to test: Cz
% look into all time points;

%% load and prepare data
for ci = 1:length(cfg.cond_names)
    load(fullfile(PATHIN,['sub-all_sPCA_sensor_', cfg.cond_names{ci},'.mat']), 'gaitERSP');
    GPM{ci} = gaitERSP.GPM(:,:,:,sii);
    %1 DT easy, 2 ST easy, 3 DTdifficult, 4 ST difficult
end
DTeasy = 1; STeasy = 2; DTdifficult = 3; STdifficult = 4;

% assemble all ersp struct of size:freqs x pnts x (chan )x sub
% 2x2 cell
allersp{1,1} = squeeze(permute(GPM{STeasy}(:,chanIdx, cfg.FOI,:),[3 1 2 4]));%ST even
allersp{1,2} = squeeze(permute(GPM{DTeasy}(:,chanIdx, cfg.FOI,:),[3 1 2 4]));%DT even
allersp{2,1} = squeeze(permute(GPM{STdifficult}(:,chanIdx, cfg.FOI,:),[3 1 2 4]));%ST uneven
allersp{2,2} = squeeze(permute(GPM{DTdifficult}(:,chanIdx, cfg.FOI,:),[3 1 2 4]));%DT ueven

%%% main effects
% surface:
% even = mean(STeasy+ DTeasy)
GPM_even = mean(cat(4,allersp{1,1}, allersp{1,2}),4);
GPM_uneven = mean(cat(4,allersp{2,1}, allersp{2,2}),4);
ME_surface = GPM_even-GPM_uneven;

% task
GPM_ST = mean(cat(4,allersp{1,1}, allersp{2,1}),4);
GPM_DT = mean(cat(4,allersp{1,2}, allersp{2,2}),4);
ME_task = GPM_ST-GPM_DT;

data_interaction = {ME_surface, ME_task};
%% main effect terrain --> dependent samples t-test
ME_surface.cfg = cfg.stats.depT_surface.stats;  % load respective configuration
ME_surface.data = allersp;                      % store data
[results.pcond, results.pgroup, results.pinter, results.statcond, results.statgroup, results.statinter] = ...
    std_stat(allersp, ME_surface.cfg);          % perfrom stat test
ME_surface.results = results;                   % add results to structure
save(fullfile(PATHOUT,file_out{1}), 'ME_surface');% save results

%% main effect task --> dependent samples t-test
ME_task.cfg = cfg.stats.depT_task.stats;    	% load respective configuration
ME_task.data = allersp;                         % store data
[results.pcond, results.pgroup, results.pinter, results.statcond, results.statgroup, results.statinter] = ...
    std_stat(allersp, ME_task.cfg);             % perfrom stat test
ME_task.results = results;                      % add results to structure
save(fullfile(PATHOUT,file_out{2}), 'ME_task'); % save results

%% interaction terrain and task --> dependent samples t-test of differences
interaction.cfg = cfg.stats.depT_interaction.stats; % load respective configuration
interaction.data = data_interaction;                         % store data
[results.pcond, results.pgroup, results.pinter, results.statcond, results.statgroup, results.statinter] = ...
    std_stat(data_interaction , interaction.cfg);         % perfrom stat test
interaction.results = results;                      % add results to structure
save(fullfile(PATHOUT,file_out{3}), 'interaction'); % save results

end
