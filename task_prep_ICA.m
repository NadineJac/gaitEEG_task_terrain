function task_prep_ICA(sii, PATH, cfg) % also call w/paths?
% prepare data for running ICA on it
% load file_in from PATHIN and store file_out in PATHOUT. Reports saved as
% reportNames
% 
% overview of performed tasks:
% - prune data
% - downsample
% - filter
% - reject channels
% - ASR w/ custom BL
% - clean line noise
% - interpolate channels
% - CAR
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-11-2022

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
    reportNames = {...
        '_chanRej.png',...
        '_ASR.png',...
        '_zapLinePlus.png'};
      
    % check whether output file exists, if so, only run script if it has
    % been modified in the meantime
    % ToDo: also check scripts that are being called!
    scriptName = [mfilename('fullpath'),'.m'];
    runScript = naj_scriptModified(scriptName, file_out);
    if runScript
        
        %% load data
        EEG = pop_loadset(file_in);
        
        %% prune data
        % discard gait initiation
        FROM = EEG.event(strcmp({EEG.event.type}, cfg.startExperiment)).latency-60*EEG.srate;
        TO = EEG.event(end).latency;
        EEG = pop_select(EEG, 'point', [FROM TO]);
        
        %% only keep EEG chans
        EEG = pop_select(EEG, 'channel', cfg.idx_EEG_chan);
        
        %% downsample
        EEG = pop_resample( EEG, cfg.sampling_rate);
        
        %% bandpass filter
        data_unfiltered = EEG.data;
        EEG = pop_eegfiltnew(EEG, 'locutoff',cfg.bandpass_fmin);
        
        % 135 Hz LPF: performing 57 point lowpass filtering, transition band width: 30 Hz
        % passband edge(s): 120 Hz, cutoff frequency(ies) (-6 dB): 135 Hz, (zero-phase, non-causal)
        EEG = pop_eegfiltnew(EEG, 'hicutoff',cfg.bandpass_fmax);
        
        % eegplot(EEG.data, 'data2', data_unfiltered);% visualize effect of
        % filtering
        clearvars data_unfiltered % free up some memory
        
        %% channel rejection w/ ASR
        % store original chanlocs in EEG.etc
        EEG.etc.orgChanlocs = EEG.chanlocs;
        EEGorg = EEG; % keep for vis only
        
        EEG = pop_clean_rawdata(EEG, ...
            'FlatlineCriterion',cfg.ASR_FlatlineCriterion,...
            'ChannelCriterion',cfg.ASR_ChannelCriterion,...
            'LineNoiseCriterion',cfg.ASR_LineNoiseCriterion,...
            'Highpass','off',...
            'BurstCriterion','off',...
            'WindowCriterion','off',...
            'BurstRejection','off',...
            'Distance','Euclidian');
        
        if ~isfield(EEG.etc, 'clean_channel_mask')
            EEG.etc.clean_channel_mask = ones(1,EEG.nbchan);
        else
            EEG.etc.rmChanLab = {EEG.etc.orgChanlocs(EEG.etc.clean_channel_mask == 0).labels};
            procInfo(si).rmChanLab =EEG.etc.rmChanLab;
        end
        EEG.etc.rmChanNum = sum(EEG.etc.clean_channel_mask == 0); % store number of removed channels
        procInfo(si).rmChanNum =EEG.etc.rmChanNum;
        
        visChanRej(EEGorg, EEG.chanlocs);
        set(gcf, 'units', 'centimeters', 'position', [0 0 30 20]);
        print(fullfile(PATHOUTsi, [ID, reportNames{1}]), '-dpng'); close;
        clearvars EEGorg
        
        %% clean w/ ASR
        % ASR params
        EEG = ASR(EEG, cfg.ASR_baseline_events, cfg.ASR_cutoff);
        print(fullfile(PATHOUTsi, [ID, reportNames{2}]), '-dpng'); close;
        
        %% clean line noise
        % from https://github.com/MariusKlug/zapline-plus as of 26.10.21
       clean_data_with_zapline_plus_eeglab_wrapper(EEG, cfg.zapLinePlus);
       print(fullfile(PATHOUTsi, [ID, reportNames{3}]), '-dpng'); close;
        
        %% interpolate channels w EEG.etc.orgChanlocs
        EEG = pop_interp(EEG, EEG.etc.orgChanlocs, 'spherical');
        
        %% full rank average Ref
        EEG = fullRankAveRef(EEG); % instead of EEG = pop_reref(EEG, []);
        
        %% save
        pop_saveset(EEG, file_out);
    end
    fprintf('Done with sub-%03d!', si');
end % subject loop

save(cfg.info_file, 'procInfo');
disp(['Done with ', mfilename])
end %end function