function cohensd = task_clusterEffectSize(sii, PATH, cfg)
% load contrast
PATHIN = PATH.sPCAsensor;
load(fullfile(PATHIN,['sub-all_contrasts.mat']), 'contrasts');

% stats
file_in = {'2x2rmANOVA_Cz.mat',...
    'depT_ME_surface.mat',...
    'depT_ME_task.mat'};
for fi = 1:length(file_in)
    load(fullfile(PATHIN,file_in{fi}));
end

% ROI
freqs = 6:2:40; % frequencies in Hz
FOI = ismember(cfg.f_axis, freqs);%their indices
chanIdx = 24;

% cluster effect size --> does not generalize for other tests or several
% clusters yet
data = contrasts.ME_surface(:,chanIdx,FOI,:);% pnts x chan x FOI x subj
idxPixel = ME_surface.results.pinter{1}<.05; % idx of pixels that are part of cluster
idx = permute(repmat(idxPixel, [1,1,size(data,2), size(data,4)]), [2,3,1,4]);%pnts x chan x freq x subj
data(~idx) = nan; % only keep cluster data
data = squeeze(mean(data,[1,2,3], 'omitnan')); % average over chan,time and freqs --> subj vector

%calculate cohen's d
cohensd = mean(data)/std(data);

end