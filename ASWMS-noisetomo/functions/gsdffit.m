function [para,resnorm,residual, exitflag] = gsdffit(xcor,lagtime,freq,nfit)
% Matlab function that does similar stuff as gsdf_fit.c. Fitting the narrow-banded cross-correlation 
% function with a five - parameter wavelet 
%	A=para(1);  is the amplitude
%	w=para(2);	is the oemga, 2*pi/T
%	sigma=para(3);	is the frequency band width
%	tp=para(4);	is the phase delay
%	tg=para(5);	is the group delay
% Written by Ge Jin
% Lamont, Colubmia University
% jinwar@gmail.com 
% April, 2012
% Modified by Ge Jin for industry data
% Conocophillips
% June, 2012
    isfigure = 0;
	xcor = xcor(:);
	lagtime = lagtime(:);

    T=1/freq;
	N=length(xcor);
	N=ceil(N/2);
    maxi=find(xcor==max(xcor),1);
    delta=lagtime(2)-lagtime(1);
    fitN=round(nfit*T/delta);
	if maxi-fitN <= 0 || maxi+fitN > length(xcor)
		para = zeros(5,1);
		resnorm = 0;
		residual = 0;
		exitflag = -8;
		return;
	end
	cutxcor = xcor(maxi-fitN:maxi+fitN);
	taxis = lagtime(maxi-fitN:maxi+fitN);

	w = 2*pi/T;
    
    maxi=find(cutxcor==max(cutxcor),1);
    t0=taxis(maxi);
	amp0=cutxcor(maxi);
	cutxcor_norm=cutxcor./amp0;

	para0=[1,w,w*.1,t0,t0];
	paraL=[0.5,w*.5,0,t0-T,taxis(1)];
	paraH=[2,w*2,w*.1*10,t0+T,taxis(end)];

	options=optimset('Display','off');

	[para,resnorm,residual, exitflag] = lsqcurvefit(@gsdfwavelet,para0, ...
								taxis,cutxcor_norm,paraL,paraH,options);
%	[para,resnorm,residual, exitflag] = lsqcurvefit(@gsdfwavelet,para0,taxis,cutxcor_norm);

	para(1)=para(1)*amp0;
	resnorm = resnorm*amp0^2;
	residual = residual*amp0;
%     
    if isfigure && T==41
        figure(99)
        fiterr = resnorm./para(1)^2./nfit./T;
        clf
        hold on
        plot(lagtime,xcor);
        plot(taxis,gsdfwavelet(para,taxis),'r--')
        title([num2str(T),' s (fiterr = ',num2str(fiterr),')']);
        xlim([-8 8]*T/delta);
%         pause;
    end

end

