% This function is used for fitting the azimuthal anisotropy for the certain grid on 
% a Rayleigh wave tomography map.
% written by Ge Jin
% jinwar@gmail.com
% Apr, 2012
%
% JBR 10/1/16 - Fit both 4 theta and 2 theta
%

function [fitstr, isophv, A2_2, A4_2, phi2_2, phi4_2]=fit_azi_anisotropy2theta4theta_2(azi,phV,comp,varargin)
	if ~isempty(varargin{1})
		phvstd=varargin{1};
    else
        phvstd=phV;
		phvstd(:)=1;
    end
    phvstd(phvstd==0) = nanmean(phvstd);
    
    n=0;
    for i=1:length(azi)
        if ~isnan(phV(i))
            n=n+1;
            fazi(n)=azi(i);
            fphv(n)=phV(i);
            fphvstd(n)=phvstd(i);
        end
    end
    
	% Initial condition
	para0(1)=nanmean(phV);
	para0(2)=0.01;
	para0(3)=0.01;
    para0(4)=90;
    para0(5)=90;
	% Lower Boundary
	paraL(1)=nanmean(phV).*0.8;
	paraL(2)=0.0;
	%paraL(3)=-20;
    paraL(3)=0;
    paraL(4)=0;
    paraL(5)=0;
	% Upper Boundary
	paraH(1)=nanmean(phV).*1.2;
	paraH(2)=0.1;
	%paraH(3)=200;
    paraH(3)=0.1;
    paraH(4)=180;
    paraH(5)=180;
    
% 	para=lsqnonlin(@errfun(para,azi,phV,phvstd),para0,paraL,paraH,azi,phV,phvstd);
%     if comp == 'Z' || comp == 'R'
%         phvfun=fittype('a*(1+d*cosd(2*(x-e)))');
%     elseif comp == 'T'
%         phvfun=fittype('a*(1+d*cosd(4*(x-e)))');
%     end
%     phvfun = fittype('a*(1+d2*cosd(2*(x-e))+d4*cosd(4*(x-e)))','coeff',{'a','d2','d4','e'});
    phvfun = fittype('a*(1+d2*cosd(2*(x-e2))+d4*cosd(4*(x-e4)))','coeff',{'a','d2','d4','e2','e4'});
    s = fitoptions(phvfun);
    s.StartPoint=para0;
    s.Lower=paraL;
    s.Upper=paraH;
    s.Weights=1./fphvstd.^2;
    
%     s.Method='NonlinearLeastSquares';
    if size(fazi,1)==1
        fazi=fazi';
    end
    if size(fphv,1)==1
        fphv=fphv';
    end
    
    fit1=fit(fazi,fphv,phvfun,s);
    
    
    fitstr=fit1;
	isophv=fitstr.a;
	A2_2=fitstr.d2;
    A4_2=fitstr.d4;
	phi2_2=fitstr.e2;
    phi4_2=fitstr.e4;

end


function err=errfun(para,azi,phV,phvstd)
	mphv=aziphv(para,azi);
	errs=(mphv-phV)./phvstd;
	err=nansum(errs.^2);
end
