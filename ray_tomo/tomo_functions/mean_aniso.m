function [ A_mean, phi_mean ] = mean_aniso(A, phi)
% Calculate the mean anisotropy given vectors of strength and orientation
%

if length(A) ~= length(phi)
    disp('A and PHI must be same length');
    return
end

% Vectors along length
Ac = A.*cosd(2*phi);
As = A.*sind(2*phi);

Asm = nanmean(As);
Acm = nanmean(Ac);

phi_mean = 0.5*atan2d(Asm,Acm);
A_mean = sqrt(Acm.^2 + Asm.^2);

end

