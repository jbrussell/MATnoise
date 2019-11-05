function fout=spectrumwhiten(fin)

	eps1=max(abs(fin)).*1.e-12;   % numerical water level to deal w/ 0 frequencies
	ta1=abs(fin+eps1);
	taper=1.0;                %normalize by 1/the amplitude spectrum
	fout=fin.*taper./ta1; 

end
