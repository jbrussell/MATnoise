function err = besselerr(tw,xsp,varargin)
% Error function for fitting two nearby station pair's frequency domain xspectrum
%
% Modified on 2/2020 by JBR
% Normalize the bessel function by its analytic envelope when is_normbessel=1 to 
% avoid fitting amplitude variations. This was a problem particularly at the 
% endpoints where signal dies out. The caveat is that low SNR parts of
% the cross-spectrum will now be amplified, so want to make sure you're only
% fitting where there's really signal!
% 
% jbrussell 8/20/2020
% is_normbessel=1 SHOULD BE USED WITH CAUTION. I have found that it tends to
% produce rougher dispersion curves if there is any hint of noise in the Bessel.
% Probably best to set is_normbessel=0.
%
% jbrussell 8/28/2020
% Removed flatness constraint (1st derivative). After some testing, found that
% it actually produces rougher dispersion curves. Still not sure why that is.
%

global tN
global waxis
global twloc
global weight 

Isfigure=0;
interpmethod = 'linear';

% damp = [1; 1; 1; 1]; % deprecated flatness constraint
damp = [1; 1; 1];
is_normbessel = 0; % no bessel normalization by default
if nargin>2 % if option is provided
     damp = varargin{1};
     if length(varargin)==2
        is_normbessel = varargin{2};
     end
end

tw1 = interp1(twloc,tw(1:tN),waxis,interpmethod);

x1 = waxis.*tw1;

% First part of error: fit the Bessel Function
%F1 = normalise(xsp);
F1 = xsp;
if is_normbessel
    F1 = F1./abs(hilbert(F1));
end
be = besselj(0,x1);
be = be./mean(abs(be)).*mean([abs(F1)]);
if is_normbessel
%     be = be./abs(hilbert(be));
    be = be./SmoothAnalyticEnv(waxis/2/pi,be);
end
F1z =  be(:) - F1(:); 
F1z = F1z.*weight(:);

% Second part of error: dispersion curve smoothness
sm = del2(tw1);
sm = sm./mean(abs(sm)).*mean(abs(F1z));

% % dispersion curve flatness (deprecated... for some reason produces rougher dispersion curves)
% fl = diff(tw1);
% fl = ( fl./mean(abs(fl)).*mean(abs(F1z)) );

% Third part of error: x1 has to be always increasing
dx = diff(tw1);
dx(find(dx>0)) = 0;

% % ORIGINAL
% err = [F1z(:); sm(:)*0.2; dx(:)*10];
% % ADD FLATNESS
% err = [F1z(:)*damp(1); sm(:)*0.2*100*damp(2); fl(:)*0.2*100*damp(3); dx(:)*10*1*damp(4)];
err = [F1z(:)*damp(1); sm(:)*0.2*100*damp(2); dx(:)*10*1*damp(3)];

% err = err./mean(abs(err))*1;

if Isfigure>0
    figure(Isfigure)
    clf
    subplot(2,1,1)
    hold on
    plot(waxis/2/pi,F1,'b','linewidth',2);
    plot(waxis/2/pi,be,'r','linewidth',2);
%     ylim([-0.02 0.02]);
%	title(real part of stacked ccf) 
% 	title( num2str(sum(err.^2)./length(err)./sum(F1.^2)*length(F1)));
%    subplot(3,1,2)
%   plot(waxis/2/pi,(F1z(:)/mean(abs(F1.*weight(:)))).^2,'k')
%        subplot(3,1,3)
%    plot(x1)
    pause;
end
end
