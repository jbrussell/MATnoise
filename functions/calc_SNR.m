function [ snr, signal_ind, snr_pos, snr_neg ] = calc_SNR( ccf, grv_min, grv_max, r, dt, isfigure_snr )
    % Calc SNR within the specified group velocity window
    % JBR 7/8/18
    
    snrdata = real(ifft((ccf)));
    snrdata = fftshift(snrdata);
    NN= length(snrdata);
    lag = ([0:NN-1]-floor(NN/2))*dt;  % build lagtime vector for plotting
    lag = [lag(lag<0), lag(lag>=0)];
    win_min = r./grv_max;
    win_max = r./grv_min;
    if win_min < 15; win_min = 0; end
    if win_max < 50; win_max = 50; end
    signal_ind = ((lag>-win_max & lag<-win_min) | (lag>win_min & lag<win_max));
    signal_amp = sum(snrdata(signal_ind).^2)/length(snrdata(signal_ind));
    noise_amp = sum(snrdata(~signal_ind).^2)/length(snrdata(~signal_ind));
    snr = signal_amp/noise_amp;
    
    % Positive SNR
    signal_ind_pos = (lag>win_min & lag<win_max);
    signal_amp_pos = sum(snrdata(signal_ind_pos).^2)/length(snrdata(signal_ind_pos));
    noise_amp_pos = sum(snrdata(~signal_ind_pos & lag>=0).^2)/length(snrdata(~signal_ind_pos & lag>=0));
    snr_pos = signal_amp_pos/noise_amp_pos;
    
    % Negative SNR
    signal_ind_neg = (lag>-win_max & lag<-win_min);
    signal_amp_neg = sum(snrdata(signal_ind_neg).^2)/length(snrdata(signal_ind_neg));
    noise_amp_neg = sum(snrdata(~signal_ind_neg & lag<=0).^2)/length(snrdata(~signal_ind_neg & lag<=0));
    snr_neg = signal_amp_neg/noise_amp_neg;
    
    
    if isfigure_snr
        figure(101); clf;
        subplot(2,1,1);
        plot(lag,snrdata,'-','color',[0.5 0.5 0.5],'linewidth',2); hold on;
        plot([win_min win_min],[-1 1]*max(abs(snrdata)),'-r','linewidth',2);
        plot(-[win_min win_min],[-1 1]*max(abs(snrdata)),'-r','linewidth',2);
        plot([win_max win_max],[-1 1]*max(abs(snrdata)),'-r','linewidth',2);
        plot(-[win_max win_max],[-1 1]*max(abs(snrdata)),'-r','linewidth',2);
        plot(lag(signal_ind),snrdata(signal_ind),'-k','linewidth',2);
        xlabel('Lag (s)');
        title(sprintf('SNR %.2f',snr));
        xlim([-300 300]);
        set(gca,'linewidth',1.5,'fontsize',16);
    
        subplot(2,1,2);
        plot(lag,snrdata,'-','color',[0.5 0.5 0.5],'linewidth',2); hold on;
        plot([win_min win_min],[-1 1]*max(abs(snrdata)),'-r','linewidth',2);
        plot(-[win_min win_min],[-1 1]*max(abs(snrdata)),'-r','linewidth',2);
        plot([win_max win_max],[-1 1]*max(abs(snrdata)),'-r','linewidth',2);
        plot(-[win_max win_max],[-1 1]*max(abs(snrdata)),'-r','linewidth',2);
        plot(lag(signal_ind),snrdata(signal_ind),'-k','linewidth',2);
        xlabel('Lag (s)');
        xlim([-2000 2000]);
        set(gca,'linewidth',1.5,'fontsize',16);
    end
    
    end
    
    