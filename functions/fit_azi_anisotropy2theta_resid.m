% This function is used for fitting the azimuthal anisotropy for the certain grid on 
% a Rayleigh wave tomography map.
% written by Ge Jin
% jinwar@gmail.com
% Apr, 2012
%
% JBR 2/16/19 - GeoPRISMS TEI
%             - This version inverts phase velocity residuals (in units of
%             percent/100) for the 2-theta terms. Assumes the istropic
%             part has already been subtracted out!
%

function [fitstr, A2_2, phi2_2]=fit_azi_anisotropy2theta_resid(azi,resid,varargin)
	if ~isempty(varargin{1})
		residstd=varargin{1};
        residstd(residstd==0) = 1;
    else
        residstd=resid;
		residstd(:)=1;
    end
    
    n=0;
    for i=1:length(azi)
        if ~isnan(resid(i))
            n=n+1;
            fazi(n)=azi(i);
            fresid(n)=resid(i);
            fresidstd(n)=residstd(i);
        end
    end
    
	% Initial condition
	para0(1)=0.01;
    para0(2)=90;
	% Lower Boundary
	paraL(1)=0.0;
	%paraL(3)=-20;
    paraL(2)=0;
	% Upper Boundary
	paraH(1)=0.1;
	%paraH(3)=200;
    paraH(2)=190;
    
% 	para=lsqnonlin(@errfun(para,azi,resid,residstd),para0,paraL,paraH,azi,resid,residstd);
%     if comp == 'Z' || comp == 'R'
%         residfun=fittype('a*(1+d*cosd(2*(x-e)))');
%     elseif comp == 'T'
%         residfun=fittype('a*(1+d*cosd(4*(x-e)))');
%     end
%     residfun = fittype('a*(1+d2*cosd(2*(x-e))+d4*cosd(4*(x-e)))','coeff',{'a','d2','d4','e'});
    residfun = fittype('d2*cosd(2*(x-e2))','coeff',{'d2','e2'});
    s = fitoptions(residfun);
    s.StartPoint=para0;
    s.Lower=paraL;
    s.Upper=paraH;
    s.Weights=1./fresidstd.^2;
    
%     s.Method='NonlinearLeastSquares';
    if size(fazi,1)==1
        fazi=fazi';
    end
    if size(fresid,1)==1
        fresid=fresid';
    end
    
    fit1=fit(fazi,fresid,residfun,s);
    
    
    fitstr=fit1;
	A2_2=fitstr.d2;
	phi2_2=fitstr.e2;
    if phi2_2 > 180
        phi2_2 = phi2_2 - 180;
    end

end

