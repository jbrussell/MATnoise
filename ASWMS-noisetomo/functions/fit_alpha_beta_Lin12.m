% This function is used for fitting the azimuthal dependence of the
% corrected amplitude decay term to extract the attenuation coefficent
% (alpha) and the local site amplification (beta)
%
% Equation 9 in Bao et al. (2016) GJI; doi:10.1093/gji/ggw151
%
% github.com/jbrussell
% 2021-05

function [fitstr, alpha, beta_tau, azi_maxamp]=fit_alpha_beta_Lin12(azi,amp,varargin)
	if ~isempty(varargin)
		ampstd=varargin{1};
    else
        ampstd=amp;
		ampstd(:)=1;
    end
    
    n=0;
    for i=1:length(azi)
        if ~isnan(amp(i))
            n=n+1;
            fazi(n)=azi(i);
            famp(n)=amp(i);
            fampstd(n)=ampstd(i);
        end
    end
	% Initial condition
	para0(1)=0;
	para0(2)=0;
	para0(3)=nanmean(amp);
	% Lower Boundary
	paraL(1)=abs(nanmean(amp))*-100;
	paraL(2)=-180;
	paraL(3)=abs(nanmean(amp))*-100;
	% Upper Boundary
	paraH(1)=abs(nanmean(amp))*100;
	paraH(2)=180;
	paraH(3)=abs(nanmean(amp))*100;
% 	para=lsqnonlin(@errfun(para,azi,amp,ampstd),para0,paraL,paraH,azi,amp,ampstd);
%     ampfun=fittype('a*(1+d*cosd(2*(x-e)))');
%     ampfun=fittype('a*sind(x) + b*cosd(x) - c');
    ampfun=fittype('a*cosd(x-b) - c');
    s = fitoptions(ampfun);
    s.StartPoint=para0;
    s.Lower=paraL;
    s.Upper=paraH;
    s.Weights=1./fampstd;
%     s.Method='NonlinearLeastSquares';
    if size(fazi,1)==1
        fazi=fazi';
    end
    if size(famp,1)==1
        famp=famp';
    end
    
    fit1=fit(fazi,famp,ampfun,s);
        
    fitstr=fit1;
	beta_tau=fitstr.a;
	azi_maxamp=fitstr.b;
	alpha=fitstr.c;
    
    if 0
        figure(100); clf;
        plot(fazi,famp,'xb'); hold on;
        x = [0:360];
%         plot(x,dlnbeta_dx*sind(x) + dlnbeta_dy*cosd(x) - alpha,'-r');
        plot(x,beta_tau*cosd(x-azi_maxamp) - alpha,'-r');
        pause;
    end

end


% function err=errfun(para,azi,amp,ampstd)
% 	mphv=aziphv(para,azi);
% 	errs=(mphv-amp)./ampstd;
% 	err=nansum(errs.^2);
% end
