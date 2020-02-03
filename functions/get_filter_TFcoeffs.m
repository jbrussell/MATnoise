function [ b_mat, a_mat ] = get_filter_TFcoeffs( frange, dt )
% Get filter coefficients for use with filtfilt

frange = sort(frange,'ascend');

% Build frequency vector
df = frange(1)/4; % 1/4 of lowest frequency
freqs = frange(1):df:frange(2);

% Loop over all frequencies
fN = 1/2/dt;
b_mat = [];
a_mat = [];
for ifreq = 1:length(freqs)-1
    [b,a] = butter(2,[freqs(ifreq)/fN, freqs(ifreq+1)/fN]);
    b_mat = [b_mat; b];
    a_mat = [a_mat; a];
end

end

