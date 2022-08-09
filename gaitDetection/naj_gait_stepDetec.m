function struct = naj_gait_prelimStepDetec(PATH, currentSubj,varargin)
% preliminary detection of gait events (initial contact and toe off)
% using peak detection of detrended and low-pass filtered acceleration data from both feet
% 1. find steps: peaks of vertical acceleration higher than set thershold
% 2. IC: within 2 steps, look for minima of anterior-posterior acceleration
% 3. TO:  get AP acceleration and find two highest peaks (toe off is inbetween, but closer to
% first peak)
%
% Caveats:
% - detection is thershold dependent and might therefore crash (individual
% acceleration patterns)
% - detection has no plausability checking, hence some misdetections might
% occur
%
% To-do:
%- add counter for not properly detected steps instead of letting
% code crash?

load([PATH,'\rawdata\SUBJ.mat'])
PATHTMP = [PATH, '\data\'];

% experimental conditions during which walking occurs
CONDS = {'initiation','difficult', 'difficult_button','easy','easy_button'};
% placement of Faros sensors and belonging channels
SIDE = {'Left','Right'};
CHAN = [1:3; 4:6];
gait_event_order = {'RightHS', 'LeftTO', 'LeftHS', 'RightTO', 'RightHS'};% order of gait events within each gait cycle
gait_timeNextHs = [.5 1.5]; % time of next RHS in s

% LPF filter design
LPF1 = 6; % [in Hz] high-frequency cut-off
LPF2 = 30; % [in Hz] high-frequency cut-off
Order = 2;
srate  = 500;
[b_low1,a_low1]=butter(Order,LPF1/(srate/2),'low');
[b_low2,a_low2]=butter(Order,LPF2/(srate/2),'low');

% thresholds for peakdetection
Step_thresh     = 600; %default, overwritten by optional input
minPeakDist = 0.5; %min distance between steps in s

for i = 1:2:length(varargin) % work for a list of name-value pairs
    if strcmp(varargin{i}, 'Step_thresh') % check if is character
        Step_thresh = varargin{i+1}; % override or add parameters to structure.
    end
end

% loop through subjects
for s = currentSubj
    FILENAME = ['naj_gait_',SUBJ.ID{s}];
    disp(['Findings steps of ', FILENAME,'...']);
    
    % load data
    load([PATH,'\rawdata\', FILENAME], 'EEG');
    
    % only keep data from experiment
    evalc('EEG = pop_select(EEG, ''point'', [EEG.event(find(strcmp({EEG.event.type}, ''start_experiment''),1,''last'')).latency-1 EEG.event(strcmp({EEG.event.type}, ''end_experiment'')).latency]);');
    
    % seperate acceleration data
    Acc = pop_select(EEG, 'nochannel', [1:67]);
    Acc2 = Acc;
    
    % detrend
    data = detrend(double(Acc.data'));
    
    % LPF
    data2 =filtfilt(b_low2,a_low2, data); % lowpass 30 Hz
    Acc.data = data2';
    
    data =filtfilt(b_low1,a_low1, data); % lowpass 6 Hz (only to detect number of steps)
    Acc2.data = data';
    
    % extract only walking sections
    for c = 1:length(CONDS)
        START =EEG.event(strcmp({Acc.event.type}, ['start_', CONDS{c}])).latency;
        TO =EEG.event(strcmp({Acc.event.type}, ['end_', CONDS{c}])).latency;
        evalc('TMP = pop_select(Acc, ''point'', [START TO]);');
        evalc('TMP2 = pop_select(Acc2, ''point'', [START TO]);');
        
        % detect gait events and add them to structure
        % process left and right foot seperatly
        for f = 1:length(SIDE)
            if strcmp(CONDS{c}, 'initiation')
                [~, idxMS] =findpeaks(TMP2.data(CHAN(f,3),:), 'MinPeakHeight', Step_thresh,'minPeakDist', minPeakDist*Acc.srate); %adjust?
                for ev = 1:length(idxMS) %ignore first and last 2 peaks?
                    EEG.event(end+1).type = [SIDE{f},'Step'];
                    EEG.event(end).latency = START+idxMS(ev);
                    EEG.event(end).duration = 1;
                end
            else
                % get all peaks of the 6 HZ LPF filtered signal vertical acceleraton that are higher than
                % threshold and further apart than the min distance
                [~, idxMS] = findpeaks(TMP2.data(CHAN(f,3),:), 'MinPeakHeight', Step_thresh, 'minPeakDist', minPeakDist*Acc.srate);
                
                %                  %del later
                %                     dat = TMP.data(CHAN(f,3),:);
                %                     figure; plot(dat), hold on;
                %                     plot(idxMS, dat(idxMS),'*');
                
                for ev = 2:length(idxMS)-1 %loop through all peak, to-do: ignore first and last 2 peaks(acceleration and decceleration) or later in pipeline?
                    
                    FROM = idxMS(ev)-.1*Acc.srate;
                    
                    %%% HEAL STRIKE %%%
                    % find peak in the 30 Hz LPF filtered vertical
                    % acceleration signal in the surronding 200 ms
                    [~, idxHS] = findpeaks(TMP.data(CHAN(f,3),FROM:FROM+.2*Acc.srate), 'MinPeakHeight', Step_thresh, 'NPeaks', 1);
                    if ~isempty(idxHS) % if present add as event
                        EEG.event(end+1).type = [SIDE{f},'HS'];
                        EEG.event(end).latency = START+FROM+idxHS;
                        EEG.event(end).duration = 1;
                    end
                    
                    %%% TOE OFF %%%
                    % get AP acceleration of half a second before
                    shift =.5*Acc.srate;
                    tmp = TMP.data(CHAN(f,1),FROM-shift:FROM-.1*Acc.srate);
                    % find two highest peaks (toe off is inbetween, but closer to
                    % first peak)
                    [~, idxTO] = findpeaks(tmp, 'NPeaks', 2, 'minPeakHeight', 200, 'SortStr' ,'descend');
                    
                    if ~isempty(idxTO)
                        EEG.event(end+1).type = [SIDE{f},'TO'];
                        EEG.event(end).latency = START+FROM-shift+mean(idxTO);% add TO in the middle of two peaks --> improve
                        EEG.event(end).duration = 1;
                    end
                    clearvars idxIC idxTO
                end
            end
        end
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');
    
    % check percaentage of valid steps (per cond)
    for c = 2:length(CONDS) %skip initiation
        disp(CONDS{c})
        START =EEG.event(strcmp({EEG.event.type}, ['start_', CONDS{c}])).latency;
        TO =EEG.event(strcmp({EEG.event.type}, ['end_', CONDS{c}])).latency;
        evalc('TMP = pop_select(EEG, ''point'', [START TO]);');
        EEG.etc.strides(c).cond = CONDS{c};
        [EEG.etc.strides(c).number,EEG.etc.strides(c).percValid] = ...
            validStrides(TMP, gait_event_order,gait_timeNextHs);
    end
    struct = EEG.etc.strides;
    % check results
    %     pop_eegplot( EEG, 1, 0, 1);
    
    % save
    EEG.setname = [FILENAME,'_steps'];
    disp('...done. Saving dataset now.')
    save([PATHTMP, EEG.setname], 'EEG');
end
end
