function [wavefield] = ricker_wavelet(t,r,vel,f)
%RICKER_WAVELET
%
% t: time axis [s]
% r: distance from source [km]
% vel: velocity of medium [km/s]
% f: frequency [Hz]
%
% jbrussell - 7/2023

% Ricker wavefield 
% https://wiki.seg.org/wiki/Dictionary:Ricker_wavelet
wavefield = (1-2.*pi^2*f.^2.*(t - r./vel).^2) .* exp(-pi^2.*f.^2.*(t - r./vel).^2);


end

