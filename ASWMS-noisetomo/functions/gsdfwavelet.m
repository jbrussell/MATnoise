function wavelet = gsdfwavelet(para,t)
% function with a five - parameter wavelet
%   A=para(1);  is the amplitude
%   w=para(2);  is the oemga, 2*pi/T
%   sigma=para(3);  is the frequency band width
%   tp=para(4); is the phase delay
%   tg=para(5); is the group delay
% Written by Ge Jin
% jinwar@gmail.com
% April, 2012
	A=para(1);
	w=para(2);
	sigma=para(3);
	tp=para(4);
	tg=para(5);
	wavelet = A.*exp(-0.5*((t-tg).*sigma).^2).*cos(w.*(t-tp));
end
