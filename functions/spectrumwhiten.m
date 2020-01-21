function fout=spectrumwhiten(fin, pctwl)
    % pctwl : waterlevel as percent of max (default was 0.001)
    
    % Calculate smooth spectra
    fsmooth = smooth(abs(fin),10);
    
    % Apply waterlevel to smoothed spectra
    waterlevel = pctwl*max(fsmooth);
    fsmooth(fsmooth<waterlevel) = waterlevel;
     
    % Whiten magnitudes
    mag = abs(fin)./fsmooth;
    
    % Phase vector
    phase = unwrap(angle(fin));
    
    % Construct white spectrum
    fout = mag.*exp(1i*phase);

    
    if 0
        figure(5); clf;
        plot(abs(fin),'-k'); hold on;
        plot(abs(fout),'-r');
    end

end
