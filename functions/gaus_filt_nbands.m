function [ data_gaus, faxis, gaus_filters ] = gaus_filt_nbands( data,periods,dt,min_width,max_width )
% JBR -- 7/25/16
% Apply gaussian filter
Nt = length(data);
[gaus_filters,faxis] = build_gaus_filter(1./periods,dt,Nt,min_width,max_width);
for ip = 1:length(periods)
    nband = data(:) .* gaus_filters(:,ip);
    data_gaus(:,ip) = nband;
end

% for ip = 1:length(periods)
% 		nband = fft_win_xcor .* [gaus_filters(:,ip); zeros(Nt-length(gaus_filters(:,ip)),1)];
% 		nband = ifft(nband);
% 		nband = 2*real(nband);
% 		nband_win_xcors(:,ip) = nband;
% end % end of periods loop
end

