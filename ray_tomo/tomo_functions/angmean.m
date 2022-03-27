function [meanangle_rad] = angmean(angledata_rad)
meanangle_rad = atan2(mean(sin(angledata_rad)),mean(cos(angledata_rad)));
end
