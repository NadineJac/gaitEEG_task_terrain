function EEG = ASR(EEG, BLevent, sd)
% run ASR (org, euclid dist) with custom BL
FROM = EEG.event(strcmp({EEG.event.type}, BLevent{1})).latency;
% TO = FROM+EEG.srate*60; %for pilot, late subtitute by:
TO = EEG.event(strcmp({EEG.event.type}, BLevent{2})).latency;
mybaseline = pop_select(EEG, 'point', [FROM:TO]); 

% extract data to be cleaned: remaining data of experiment after snipped extracted for ASR
% (has to remain continous)
data_old = EEG.data; % store data for visualization

% clean bursts, no extra hpf!!!!
EEG = clean_asr(EEG,sd,[],[],[],mybaseline);

% visualize
figure(); eegplot(EEG.data, 'data2', data_old);
end
