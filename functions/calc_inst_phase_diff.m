function [dphi] = calc_inst_phase_diff(ccf1,ccf2)
% Calculate differential phase between two time series
% 
% -ive = signal 2 leads signal 1 (retrograde if 1=Z and 2=R)
% +ive = signal 1 leads signal 2 (prograde if 1=Z and 2=R)
%
% jbrussell 10/2025

ccf1_analytic = hilbert(ccf1);
ccf2_analytic = hilbert(ccf2);
phase_1 = angle(ccf1_analytic);
phase_2 = angle(ccf2_analytic); 
dphi = 180/pi*wrapToPi(phase_1 - phase_2); % Element-wise phase difference

end

