function [ ccf_ ] = ftan_win( ccf,max_grv,min_grv )
% Does ftan windowing in time domain and calculates spectra for each group
% velocity. Outputs spectra for each velocity window for a station pair.
%
% grv_f_out(igrv,ifreq)
%
% JBR 3/5/18
%
grvs = min_grv:0.01:max_grv;

for igrv = 1:length(grvs)
    % Define window limits
    wmin_grv = grvs(igrv)-7*
    wmax_grv = 
    % Window CCF
    % positive side first
    time = [0:N-1];
    t1_pos = sta1sta2_dist{ista1}(nstapair)/max_grv;
    t2_pos = sta1sta2_dist{ista1}(nstapair)/min_grv;
    Iset0_pos = (time >= t1_pos) & (time <= t2_pos);
    ccf_ifft_full_pos = ccf_ifft_full;
    ccf_ifft_full_pos(Iset0_pos)= cos_taper(ccf_ifft_full_pos(Iset0_pos));
    ccf_ifft_full_pos(~Iset0_pos)= 0;
    % negative side
    time = [-N:-1];
    t1_neg = sta1sta2_dist{ista1}(nstapair)/-max_grv;
    t2_neg = sta1sta2_dist{ista1}(nstapair)/-min_grv;
    Iset0_neg = (time <= t1_neg) & (time >= t2_neg);
    ccf_ifft_full_neg = ccf_ifft_full;
    ccf_ifft_full_neg(Iset0_neg)= cos_taper(ccf_ifft_full_neg(Iset0_neg));
    ccf_ifft_full_neg(~Iset0_neg)= 0;

    ccf_ifft_full = ccf_ifft_full_pos + ccf_ifft_full_neg;

    ccf_win = fft(ccf_ifft_full);

    %ccf_ifft = real(ifft(2*ccf([1:N/2+1]),N)); % inverse FFT to get time domain
    ccf_ifft = real(ifft(2*ccf_win([1:N/2+1]),N)); % inverse FFT to get time domain

    if 1
        %rearrange and keep values corresponding to lags: -(len-1):+(len-1)
        ccf_ifft_plot = [ccf_ifft(end-N+2:end) ; ccf_ifft(1:N)];
        plot(ccf_ifft_plot
        pause;
    end
end


end

