function err = besselerr_dist(c,xsp,omega,r,varargin)
% Error function for fitting average phase velocity as a function of
% distance for a single frequency
%
% Modified on 10/2021 by JBR
%

Isfigure=0;

% damp = [1; 1; 1; 1]; % deprecated flatness constraint
is_normbessel = 0; % no bessel normalization by default
if nargin>4 % if option is provided
     is_normbessel = varargin{1};
end

% x1 = omega.*tw;
x1 = omega .* r ./ c;

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


% % ORIGINAL
err = F1z(:);

if Isfigure>0
    figure(Isfigure)
    clf
    subplot(2,1,1)
    hold on
    plot(r,F1,'.','linewidth',2);
    plot(r,be,'.','linewidth',2);
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
