function task_sPCA_sensor(sii, PATH, cfg)
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
% - standing BL
%   - extract, mark epochs exeeding threshold, compute average power
% - gait ERSP
%   - extract data of all gait cond, mark epochs exceeding threshold,
%   compute ERSP
%   - run sPCA on averaged ERSP, keep eigenvalues
%   - extract data of EACH gait cond, mark epochs exeeding threshold,
%   compute ERSP, correct EACH cond w/ previously stored eigenvalues
%   - plot (un)corrected ERSPs and GPMs
% - save timeseries data for standing BL + walk conds for source sPCA
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-12-2022

% directories
PATHIN = PATH.ICAclean;
PATHOUT = PATH.sPCAsensor;

% load processing info
load(cfg.info_file);

for si = sii % loop through all subjects
    
    % sub ID
    ID = sprintf('sub-%03d',si);
    
    PATHOUTsi = fullfile(PATHOUT, ID);
    if ~exist(PATHOUTsi), mkdir(PATHOUTsi);end
    
    % filenames:
    % input file
    file_in = fullfile(PATHIN, ID,[ID,'_ICArej.set']);
    % output file
    file_out = fullfile(PATHOUTsi,[ID,'_sPCA_sensor.mat']); % file
    
    % report files
    reportNames = {'_standBL_power.png',...
        '_allChan_ERSP.png',...
        '_allChan_GPM.png'};
    
    % check whether output file exists, if so, only run script if it has
    % been modified in the meantime
    runScript = naj_scriptModified([mfilename('fullpath'),'.m'], file_out);
    if runScript
        EEG = pop_loadset(file_in);
        
        %% standing BL TF
        % extract data
        EEG_stand = pop_select(EEG, 'point', [EEG.event(ismember({EEG.event.type}, cfg.cond_stand)).latency]);
        
        % only keep data not exceeding threshold
        EEG_stand = threshContinous(EEG_stand, cfg.V_tsh);
        procInfo(si).standBL_valid = sum(EEG_stand.etc.valid_eeg)/ EEG_stand.pnts;
        EEG_stand.data = EEG_stand.data(:,EEG_stand.etc.valid_eeg); 
        
        [F_Rest, Noise_cov] = baselineF(EEG_stand, cfg.f_axis, cfg.FWHM); % chan x freq
        print(fullfile(PATHOUTsi, [ID, reportNames{1}]),'-dpng');close;
        
        %% gait cycle ERSPs
        % extract data (all gait conditions together)
        FROM = [EEG.event(ismember({EEG.event.type}, cfg.cond_walk(:,1))).latency];
        TO =  [EEG.event(ismember({EEG.event.type}, cfg.cond_walk(:,2))).latency];
        EEG_block = pop_select(EEG, 'point', [FROM; TO]' );
        EEG_block = threshContinous(EEG_block, cfg.V_tsh);
        
        % time-frequency transformation w/ morlet-wavelets
        [Gait_avg, ERSP, GPM, cycle_cnt, valid_cycle_cnt] = gait_ersp(EEG_block, F_Rest,...
            cfg.N_freq, cfg.f_axis, cfg.FWHM,...
            cfg.gait_event, cfg.gait_timeNextHs, cfg.gait_event_order);
        
        %% denoising: spectral PCA as suggested by Seeber et al., 2015
        % spectral PCA denoising, removing broadband component
        % please check function for futher details
        [ERSP_corr, GPM_corr, PSC1, ~,V] = specPCAdenoising(ERSP);
        
        % save all info together
        gaitERSP.ID        = ID;
        gaitERSP.cond      = cfg.cond_names;
        gaitERSP.Noise_cov = Noise_cov;% noise cov for kernel computation
        gaitERSP.F_Rest    = F_Rest;
        gaitERSP.TF        = Gait_avg;
        gaitERSP.ERSP_uncor = ERSP;
        gaitERSP.GPM_uncor  = GPM;
        gaitERSP.ERSP       = ERSP_corr;
        gaitERSP.GPM        = GPM_corr;
        gaitERSP.PSC1       = PSC1;
        gaitERSP.numStrides= cycle_cnt;
        gaitERSP.numValidStrides = valid_cycle_cnt;
        gaitERSP.chanlocs   = EEG.chanlocs;
        save(file_out, 'gaitERSP');
        
        procInfo(si).numStrides = gaitERSP.numStrides;
        procInfo(si).numValidStrides = gaitERSP.numValidStrides;
        
        %% condition specific spectral PCA denoising
        for ci = 1:length(cfg.cond_walk) % condition loop
            
                    % extract data for each gait condition
            FROM = [EEG.event(ismember({EEG.event.type}, cfg.cond_walk(ci,1))).latency];
            TO =  [EEG.event(ismember({EEG.event.type}, cfg.cond_walk(ci,2))).latency];
            EEG_block = pop_select(EEG, 'point', [FROM; TO]' );
            
            EEG_block = threshContinous(EEG_block, cfg.V_tsh); 
            
            % calculate condition specific gait ERSP
            [Gait_avg, ERSP, GPM, cycle_cnt, valid_cycle_cnt] = gait_ersp(EEG_block, F_Rest,...
                cfg.N_freq, cfg.f_axis, cfg.FWHM,...
                cfg.gait_event, cfg.gait_timeNextHs, cfg.gait_event_order);
            
            % correct with PCA results from all Conds together (V)
            [ERSP_corr, GPM_corr, PSC1] = specPCAdenoising(ERSP, V);
            
            % save all info together
            gaitERSP.ID         = ID;
            gaitERSP.cond       = cfg.cond_names{ci};
            gaitERSP.Noise_cov  = Noise_cov;% noise cov for kernel computation
            gaitERSP.F_Rest     = F_Rest;
            gaitERSP.TF         = Gait_avg;
            gaitERSP.ERSP       = ERSP;
            gaitERSP.GPM        = GPM;
            gaitERSP.ERSP       = ERSP_corr;
            gaitERSP.GPM        = GPM_corr;
            gaitERSP.PSC1       = PSC1;
            gaitERSP.numStrides = cycle_cnt;
            gaitERSP.numValidStrides = valid_cycle_cnt;
            gaitERSP.chanlocs   = EEG.chanlocs;
            save([file_out(1:end-4), '_', cfg.cond_names{ci}, '.mat'], 'gaitERSP');
            
            % save info on number of strides
            procInfo(si).(cfg.cond_names{ci}).numStrides = gaitERSP.numStrides;
            procInfo(si).(cfg.cond_names{ci}).numValidStrides = gaitERSP.numValidStrides;
            
            
            %% plot ERSPs
            % plotID: table with info about chanNames and plottingGrid
            load(fullfile(PATH.code, 'chanSubplotID_LA64.mat'));
            
            % ERSP uncorrected
            plot_allChanERSP(ERSP, plotID, gaitERSP.chanlocs,cfg.f_axis, cfg.t_axis, 9, ...
                [cfg.cond_names{ci}, ' uncorrected ERSP of ', ID], 'Power change to standing BL (dB)',...
                'off');
            print(fullfile(PATHOUTsi, [ID, '_', cfg.cond_names{ci}, '_uncorr', reportNames{2}]),'-dpng');
            close;
            
            % ERSP corrected
            plot_allChanERSP(ERSP_corr, plotID, gaitERSP.chanlocs,cfg.f_axis, cfg.t_axis,...
                5, [cfg.cond_names{ci}, ' ERSP of ', ID], 'Power change to standing BL (dB)',...
                'off');
            print(fullfile(PATHOUTsi, [ID, '_', cfg.cond_names{ci}, reportNames{2}]),'-dpng');
            close;
            
            % GPM uncorrected
            plot_allChanERSP(GPM, plotID, gaitERSP.chanlocs,cfg.f_axis,cfg.t_axis,...
                2, [cfg.cond_names{ci}, ' uncorrected GPM of ', ID], 'Power change to mean gait cycle BL (dB)',...
                'off')
            print(fullfile(PATHOUTsi, [ID, '_', cfg.cond_names{ci}, '_uncorr', reportNames{3}]),'-dpng');
            close;
            
            % GPM corrected
            plot_allChanERSP(GPM_corr, plotID, gaitERSP.chanlocs,cfg.f_axis,cfg.t_axis,...
                1, [cfg.cond_names{ci}, ' GPM of ', ID], 'Power change to mean gait cycle BL (dB)',...
                'off')
            print(fullfile(PATHOUTsi, [ID, '_', cfg.cond_names{ci}, reportNames{3}]),'-dpng');
            close;
        end % condition loop
    end
    fprintf('Done with sub-%03d!\n', si);
end % subject loop

save(cfg.info_file, 'procInfo');
disp(['Done with ', mfilename])
end % end function