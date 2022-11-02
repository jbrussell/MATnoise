% This function is used for fitting the azimuthal anisotropy for the certain grid on 
% a Rayleigh wave tomography map.
% written by Ge Jin
% jinwar@gmail.com
% Apr, 2012

function [fitstr, isophv, A_2, phi_2]=fit_azi_anisotropy(azi,phV,varargin)
	if ~isempty(varargin)
		phvstd=varargin{1};
    else
        phvstd=phV;
		phvstd(:)=1;
    end
    
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
	para0(3)=90;
	% Lower Boundary
	paraL(1)=nanmean(phV).*0.8;
	paraL(2)=0.0;
	paraL(3)=-20;
	% Upper Boundary
	paraH(1)=nanmean(phV).*1.2;
	paraH(2)=0.1;
	paraH(3)=200;
% 	para=lsqnonlin(@errfun(para,azi,phV,phvstd),para0,paraL,paraH,azi,phV,phvstd);
    phvfun=fittype('a*(1+d*cosd(2*(x-e)))');
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
	A_2=fitstr.d;
	phi_2=fitstr.e;

end


function err=errfun(para,azi,phV,phvstd)
	mphv=aziphv(para,azi);
	errs=(mphv-phV)./phvstd;
	err=nansum(errs.^2);
end
