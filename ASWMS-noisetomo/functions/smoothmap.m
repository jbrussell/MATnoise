function smap=smoothmap(xi,yi,map,D)
% This function is used to make a gaussian smoothing on a grd map

% D is average distance

smap = map;
% smap2 = map;
for i=1:length(xi(:))
	dist = distance(xi(i),yi(i),xi,yi);
	dist = deg2km(dist);
    
%     % Old boxcar smoothing function
% 	ind = find(dist<D);
% 	smap(i) = nanmean(map(ind));
    
    % New gaussian smoothing function
    gaus2d = 1./(2*pi*D.^2) .* exp(-0.5*dist.^2 ./ D.^2);
    gaus2d(isnan(map)) = nan;
    smap(i) = nansum(map(:).*gaus2d(:)) ./ nansum(gaus2d(:));
    
end
ind = find(isnan(map(:)));
smap(ind) = NaN;
% smap2(ind) = NaN;

if 0 
    figure(1); clf;
    subplot(2,2,1);
    imagesc(map);
    cb = colorbar;
    subplot(2,2,3);
    imagesc(smap);
    colorbar;
    caxis(cb.Limits);
    subplot(2,2,4);
    imagesc(smap2);
    colorbar;
    caxis(cb.Limits);
end

