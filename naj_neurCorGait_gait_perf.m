function naj_neurCorGait_gait_perf(sii, PATH, cfg)
% calculate stride time and stride time variability, save as .csv
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
% - load data
% - extract walking conditions
% - get valid strides
% - claculate papameters
% - save in long and short format
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-12-2022

% directories
PATHIN = PATH.rawdata;
PATHOUT = PATH.behave;
if ~exist(PATHOUT), mkdir(PATHOUT);end

% output files
file_out =  {fullfile(PATHOUT,['sub-all_gait_performance.csv']),...
    fullfile(PATHOUT,['sub-all_strideTime.csv']),...
    fullfile(PATHOUT,['sub-all_strideTimeVar.csv'])}; % file


% set up tables
gait_perf_long = table('Size', [length(sii)*4, 5], ...
    'VariableTypes', {'string', 'string','string','double','double'},...
    'VariableNames', {'ID', 'f_surface','f_task', 'AvgStrideTime', 'strideTimeVar'});
strideTime_short = table('Size', [length(sii), 5], ...
    'VariableTypes', {'string', 'double', 'double','double','double'},...
    'VariableNames', [{'ID'}, cfg.cond_names]);
strideTimeVar_short = table('Size', [length(sii), 5], ...
    'VariableTypes', {'string', 'double', 'double','double','double'},...
    'VariableNames', [{'ID'}, cfg.cond_names]);

for si = sii % loop through all subjects
    
    % sub ID
    ID = sprintf('sub-%03d',si);
    
    % filenames:
    % input file
    file_in = fullfile(PATHIN, ID,'eeg',[ID,'_task-neurCorrYoung_eeg.set']);
    
    disp('entering loop')
    
    % load EEG
    EEG = pop_loadset(file_in);
    
    count = si*4-4; %counter for long format;
    for ci = 1:length(cfg.cond_walk) %condition loop
        count = count+1; % update counter
        
        % store participant ID and respective condition
        gait_perf_long.ID(count)         = ID;
        gait_perf_long.f_surface(count)  = cfg.f_surface(ci);
        gait_perf_long.f_task(count)     = cfg.f_task(ci);
        
        % select data from condition
        FROM    = [EEG.event(ismember({EEG.event.type}, cfg.cond_walk(ci,1))).latency];
        TO      =  [EEG.event(ismember({EEG.event.type}, cfg.cond_walk(ci,2))).latency];
        EEG_block  = pop_select(EEG, 'point', [FROM; TO]' );
        
        % calculate parameters
        % check whether steps are "plausible" i.e. not too long (check), correct
        % oder of events? (could epoch data and use only "correct
        idxHS = find(strcmp({EEG_block.event.type}, cfg.gait_event));
        strideTime = [];
        for cycle_cnt = 1:length(idxHS)-1 % stride loop
            cycle_edge = round([EEG_block.event(idxHS(cycle_cnt)).latency,...
                EEG_block.event(idxHS(cycle_cnt+1)).latency-1]); % first and last frame of gait cycle
            cycle_event = {EEG_block.event([idxHS(cycle_cnt):idxHS(cycle_cnt+1)]).type}; % labels of all events within this cycle
            % only keep labels of gait events to check their order:
            cycle_gaitEvent = cycle_event(contains(cycle_event,cfg.gait_event_order));
            strideDur = (cycle_edge(2)-cycle_edge(1))/EEG.srate;
            if cfg.gait_timeNextHs(1) <= strideDur &&... % check time until next HS
                    strideDur <= cfg.gait_timeNextHs(2) &&...
                    all(ismember(cfg.gait_event_order,cfg.gait_event_order))% oder of gait events correct
                strideTime(end+1) = strideDur;
            end
        end
        
        %             w/o checks
        %             strideTime   = [diff([EEG_block.event(strcmp({EEG_block.event.type},'RightHS')).latency]/EEG.srate)]; % stride time
        
        gait_perf_long.AvgStrideTime(count)  = mean(strideTime);               % average stride time
        gait_perf_long.strideTimeVar(count)  = ...
            (std(strideTime)/gait_perf_long.AvgStrideTime(count))*100; % stride time variability
        
        %  also store in short format (for JASP)
        strideTime_short.(ci+1)(si)     = gait_perf_long.AvgStrideTime(count);
        strideTimeVar_short.(ci+1)(si)  = gait_perf_long.strideTimeVar(count);
    end % end condition loop
fprintf('Done with sub-%03d!\n', si);
end % subject loop

% save data
writetable(gait_perf_long, file_out{1});
writetable(strideTime_short, file_out{2});
writetable(strideTimeVar_short, file_out{3});

disp(['Done with ', mfilename])
end %end function
