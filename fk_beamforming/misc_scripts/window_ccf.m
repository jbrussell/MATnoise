function [ccf_win] = window_ccf(ccf, grv_min, grv_max, r, dt, isfigure)
% Window cross-correlation function in the time-domain
%
% Input:
% - ccf: cross-correlation function in the frequency domain
% - grv_min: [km/s] minimum group velocity
% - grv_max: [km/s] maximum group velocity
% - r: [km] station separation
% - dt: [s/sample] sample rate
% - isfigure: Plot waveforms? 1 or 0

ccf_ifft = ifft(ccf);
ccf_ifft = fftshift(ccf_ifft); % Rearrange for windowing
NN= length(ccf_ifft);
lag = ([0:NN-1]-floor(NN/2))*dt;  % build lagtime vector for plotting
lag = [lag(lag<0), lag(lag>=0)];
win_min = r./grv_max;
win_max = r./grv_min;
if win_min < 15; win_min = 0; end
if win_max < 50; win_max = 50; end
% signal_ind = ((lag>-win_max & lag<-win_min) | (lag>win_min & lag<win_max));
signal_ind_neg = lag>-win_max & lag<-win_min;
signal_ind_pos = lag>win_min & lag<win_max;

ccf_ifft_win = ccf_ifft;
% Taper edges of group velocity window
ccf_ifft_win(signal_ind_neg) = cos_taper(ccf_ifft_win(signal_ind_neg));
ccf_ifft_win(signal_ind_pos) = cos_taper(ccf_ifft_win(signal_ind_pos));
% Zero outside of window
ccf_ifft_win(~(signal_ind_neg | signal_ind_pos)) = 0;

if isfigure
    figure(101); clf;
    subplot(2,1,1);
    plot(lag,ccf_ifft,'-','color',[0 0 0],'linewidth',2); hold on;
    plot(lag,ccf_ifft_win,'--','color',[1 0 0],'linewidth',2); hold on;
    plot([win_min win_min],[-1 1]*max(abs(ccf_ifft)),'-b','linewidth',2);
    plot(-[win_min win_min],[-1 1]*max(abs(ccf_ifft)),'-b','linewidth',2);
    plot([win_max win_max],[-1 1]*max(abs(ccf_ifft)),'-b','linewidth',2);
    plot(-[win_max win_max],[-1 1]*max(abs(ccf_ifft)),'-b','linewidth',2);
    xlabel('Time (s)');
    xlim((win_max+100)*[-1 1]);
    set(gca,'linewidth',1.5,'fontsize',16);
    title('Unfiltered CCF');
end

% Transform back to frequency domain
ccf_ifft_win = ifftshift(ccf_ifft_win); % Shift window back
ccf_win = fft(ccf_ifft_win);

end

