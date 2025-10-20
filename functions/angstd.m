function [angstd_rad] = angstd(angledata_rad)
% Calculate circular std
x_sum = sum(cos(angledata_rad));
y_sum = sum(sin(angledata_rad));
N = size(angledata_rad,1);
R = sqrt(x_sum.^2 + y_sum.^2) ./ N;
angstd_rad = sqrt(-2 * log(R));
end

