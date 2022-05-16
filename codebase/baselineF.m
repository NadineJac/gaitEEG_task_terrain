function [F_Rest, Noise_cov] = baselineF(EEG, f_axis, FWHM)

data = permute(EEG.data, [2,1]); % pnts x chans

% CAR (common average refrence), make sure you have 'clean' data before
data = bsxfun(@minus, data, mean(data,2));

% compute covariance matrix for inverse kernel computations (brainstorm)
Noise_cov = cov(data);

% Time-frequency analysis (function adapted by Seeber from brainstorm)
TF_Rest = morlet_transform_fast(data,[0,1/EEG.srate],f_axis,1,FWHM,'n');

% average over time, keep magnitude (not power -> would amplify outliers)
F_Rest(1,:,:) = squeeze(mean(abs(TF_Rest(10+1:end-10,:,:)),1)); 
% average over time, keep magnitude (not power -> would amplify outliers)

% visualize
figure(); set(gcf, 'position', [0 0 600 500]);
plot(f_axis, squeeze(F_Rest)', 'k-');
ylabel('Amplitude (\muV)'), ylabel('Frequency (Hz)');
grid on; box off
title('Mean standing baseline power of each channel');
end