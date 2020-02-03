function [ data_FTN ] = FTN( datan, b, a )
% Frequency-Time Normalization
% Shen et al. 2012 (doi:10.1785/0120120023)
%
% This is to be used instead of spectral whitening and one-bit
% normalization. No need to segments containing earthquakes.
%
% b, a : Precalculated transfer function coefficients of shape [Nfreq,Ncoeff];
%

% Loop over all frequencies and filter
dataf = zeros(length(datan),size(a,1)-1);
datan = detrend(datan);
datan = cos_taper(datan);
for ifreq = 1:size(a,1)
    % MATLAB built in zero-phase filter
%     dataf(:,ifreq) = filtfilt(b(ifreq,:),a(ifreq,:),datan);

    % Faster FilterM function (https://www.mathworks.com/matlabcentral/fileexchange/32261-filterm)
    dataf(:,ifreq) = FiltFiltM(b(ifreq,:),a(ifreq,:),datan);
    
    dataf(:,ifreq) = detrend(dataf(:,ifreq));
    dataf(:,ifreq) = cos_taper(dataf(:,ifreq));
end
data_FTN = sum( dataf ./ abs(hilbert(dataf)) ,2);
data_FTN = detrend(data_FTN);
data_FTN = cos_taper(data_FTN);

% Plots
if 0
    dt = 1;
    figure(222);
    t = [0:dt:length(datan)-1];
    NFFT = length(datan);
    f = 1/dt/2*linspace(0,1,NFFT/2+1);
    % raw
    fftdatan = abs(fft(datan));
    subplot(3,2,1);
    plot(t,datan,'-k'); axis tight;
    subplot(3,2,2);
    plot(f,fftdatan(1:NFFT/2+1),'-b'); axis tight;
    % one-bit normalization
    data_OBN = runwin_norm(datan);
    data_OBN = real(ifft( spectrumwhiten_smooth(fft(data_OBN),0.001) ));
    fftdata_OBN = abs(fft(data_OBN));
    subplot(3,2,3);
    plot(t,data_OBN,'-k'); axis tight;
    subplot(3,2,4);
    plot(f,fftdata_OBN(1:NFFT/2+1),'-b'); axis tight;
    % Frequency-time normalization
    fftdata_FTN = abs(fft(data_FTN));
    subplot(3,2,5);
    plot(t,data_FTN,'-k'); axis tight;
    subplot(3,2,6);
    plot(f,fftdata_FTN(1:NFFT/2+1),'-b'); axis tight;
    
    pause;
end

end

