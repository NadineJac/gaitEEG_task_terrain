function task_ICA_decomp(sii, PATH, cfg)
% perform AMICA on interim data and back-project obtained weights
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
% - w/interim data
%    - cut 1s epochs and reject (threshold, probability)
%    - run AMICA (rank reduced)
% - backproject to original data
% - fit dipoles (no thresh)
% - IC label
%
% Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-11-2022

% directories
PATHIN = PATH.prepICA;
PATHOUT = PATH.ICAdecomp;

% load processing info
load(cfg.info_file);

for si = sii % loop through all subjects
    
    % sub ID
    ID = sprintf('sub-%03d',si);
    
    PATHOUTsi = fullfile(PATHOUT, ID);
    if ~exist(PATHOUTsi), mkdir(PATHOUTsi);end
    
    % filenames:
    % input file
    file_in = fullfile(PATHIN, ID,[ID,'_prepICA.set']);
    % output file
    file_out = fullfile(PATHOUTsi,[ID,'_ICAdecomp.set']); % file
    % report files
    reportNames = {'_allComps.png'};
    
    % check whether output file exists, if so, only run script if it has
    % been modified in the meantime
    %also check scripts that are being called!
    runScript = naj_scriptModified([mfilename('fullpath'),'.m'], file_out);
    if runScript
        disp('entering loop')
        %% load data
        EEG = pop_loadset(file_in);
        
        %% Interim data
        % save interim results in EEGTEMP
        EEGTEMP = EEG;
        
        % filter before ICA
        disp(fprintf('Filter interim dataset: %3i Hz...', cfg.ICA_bandpass_fmin))
        EEGTEMP = pop_eegfiltnew(EEG, 'locutoff',cfg.ICA_bandpass_fmin);
        
        disp('Reject bad segments...')
        % pseudo-epochs, 1 second long, unrelated to task structure
        EEGTEMP = eeg_regepochs(EEGTEMP,1);
        
        % remove epochs
        EEGTEMP = pop_eegthresh(EEGTEMP,1,[1:EEGTEMP.nbchan],...
            -1*cfg.ICA_REJthresh,cfg.ICA_REJthresh,EEGTEMP.xmin,EEGTEMP.xmax,0,1); % threshold
        EEGTEMP = pop_jointprob(EEGTEMP,1,1:size(EEGTEMP.data,1),cfg.ICA_REJsd,cfg.ICA_REJsd,0,1); % probability
        
        %% run AMICA
        disp('Run AMICA...')
        % attention: AMICA data has to be chans x frames, i.e. not epoched!
        outdir  = fullfile(PATHOUTsi, 'eegamicaout');
        runamica15(EEGTEMP.data(:,:),...
            'num_models',cfg.AMICA_num_models,...
            'outdir',outdir, ...
            'max_iter', cfg.AMICA_max_iter,...
            'numprocs', cfg.AMICA_numprocs,...
            'max_threads', cfg.AMICA_max_threads,...
            'pcakeep', EEG.nbchan-(EEG.etc.rmChanNum+1));
        
        % load computed and stored weights
        disp('Done! Load weights into original data.')
        EEG = pop_loadmodout(EEG,outdir);
        
        %% fit dipoles 
        % --> output only used for comparison of post ICA processing (supplemental), otherwise may be discarded
        
        % set up head model
        disp('Dipole fitting...')
        EEG = pop_dipfit_settings( EEG,...
            'hdmfile',[PATH.eeglab, '\\plugins\\dipfit4.1\\standard_BEM\\standard_vol.mat'],...
            'coordformat','MNI',...
            'mrifile',[PATH.eeglab, '\\plugins\\dipfit4.1\\standard_BEM\\standard_mri.mat'],...
            'chanfile',[PATH.eeglab, '\\plugins\\dipfit4.1\\standard_BEM\\elec\\standard_1005.ced'],...
            'chansel',[1:EEG.nbchan]);
        
        % estimate dipoles (may be excloded based on position and residual variance whne setting up a study as well)
        EEG = pop_multifit(EEG, 1:EEG.nbchan,'threshold',100);
        
        %% IC label
        disp(['Labelling components of ', ID])
        EEG = iclabel(EEG);
        
        %% update info
        EEG.setname = file_out;
        
        %% visualization
        disp('Done. Plotting ICs now')
        
        % plot topographies w/ dipole locations
        pop_topoplot(EEG,0, [1:size(EEG.icawinv,2)],EEG.setname,[] ,1,'electrodes','on','iclabel','on');
        print(fullfile(PATHOUTsi, [ID, reportNames{1}]), '-dpng'); % save
        close;
        
        %% save
        disp(['Saving dataset ', EEG.setname])
        pop_saveset(EEG, file_out);
    end
    fprintf('Done with sub-%03d!', si');
end % end subject loop

save(cfg.info_file, 'procInfo');
disp(['Done with ', mfilename])
end % end function