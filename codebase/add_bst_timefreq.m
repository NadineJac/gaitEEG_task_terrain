% Martin Seeber, BCI Lab 2012
% Institute of Knowledge Discovery, TU Graz
% adapted by Nadine Jacobsen, Jan 2022: added DataType input, commented out
% ll24-25

function add_bst_timefreq(Protocol,Subject,Condition,sMat, DataType)

    % %% ===== start BRAINSTORM =====
    % Add brainstorm.m path to the path
    addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    % If brainstorm is not running yet: start_clac brainstorm without the GUI
    if ~brainstorm('status')
        brainstorm nogui
    end
    
    [sStudy, iStudy] = bst_get('StudyWithCondition', [Subject '\' Condition]);
    
    if isempty(iStudy)
        iStudy = db_add_condition(Subject, Condition);
        [sStudy, iStudy] = bst_get('StudyWithCondition', [Subject '\' Condition]);
    end
    
    sMat.DataType=DataType; %changed by naj recordings';
%     sMat.Device = 'Unknown';
%     sMat.ChannelFlag = ones(length(sMat.F(:,1)),1);    
    sMat.nAvg=1;
    
    OutputFile=fullfile(bst_get('BrainstormDbDir'),Protocol,'data', Subject,Condition, ['timefreq_' sMat.Comment]);
    
    save(OutputFile, '-struct', 'sMat');

    iStudyToRedraw = iStudy;

    db_add_data(iStudyToRedraw, OutputFile, sMat);
    
    iStudyToRedraw = unique(iStudyToRedraw);
    % Unload all datasets
    bst_memory('UnloadAll', 'Forced', 'KeepScouts');
    %Update results links in target study
    db_links('Study', iStudyToRedraw);
    % Update tree 
    panel_protocols('UpdateNode', 'Study', iStudyToRedraw);
    % Select target study as current node
    panel_protocols('SelectStudyNode', iStudyToRedraw(1) );
    % Save database
    db_save();