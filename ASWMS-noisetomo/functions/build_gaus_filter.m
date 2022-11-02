function [gausf,faxis] = build_gaus_filter(centf,dt,N,minwidth,maxwidth,isfigure)
% Function to build up gaussian filters within low frequency bands
% the output gausf has a structure as: 
% gausf(:,ifreq)=exp(-(faxis-centf(ifreq)).^2./2./(width(ifreq)*centf(ifreq)).^2);

faxis = [0:floor(N/2)]/dt/N;

if ~exist('isfigure','var')
	isfigure = 0;
end

% Linear varify the width of bands
if max(centf)-min(centf)==0
	width = maxwidth
else
	width = maxwidth - (maxwidth-minwidth)./(max(centf)-min(centf)).*(centf-min(centf));
end

for ifreq = 1:length(centf)
		gausf(:,ifreq)=exp(-(faxis-centf(ifreq)).^2./2./(width(ifreq)*centf(ifreq)).^2);
end 

if isfigure
	figure(isfigure)
	clf
	[xi yi] = ndgrid(faxis,centf);
	surface(xi,yi,gausf);
	shading flat;
	xlabel('Frequency axis');
	ylabel('center frequency');
	xlim([0 max(centf)*1.5])
end

end
