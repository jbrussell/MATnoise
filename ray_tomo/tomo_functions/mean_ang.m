function [ azi_mean ] = mean_ang(lat, lon)
% Calculate the mean azimuth for a path (in lat,lon) by averaging the
% vectors
%

if length(lat) ~= length(lon)
    disp('LAT and LON must be same length');
    return
end

[arclen, az] = distance(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end),wgs84Ellipsoid);

% Vectors along length
Vx = arclen.*cosd(az);
Vy = arclen.*sind(az);

azi_mean = atan2d(mean(Vy),mean(Vx));

end

