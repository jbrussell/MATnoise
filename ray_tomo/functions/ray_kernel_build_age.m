function [G,count_all] = ray_kernel_build_age(ray, xnode, ynode, agenode, agebins)
%input: ray: nray*4 matrix; each row: [x1,y1,x2,y2]
%       xnode and ynode are grid axis for lat and long
% Output: G is the kernal for Vx and Vy
% g(i,l) = r(i,l) is length of i th ray in lth pixel;
% build data kernel
% written by Yang Zha, modified for ray theory for ambient noise analysis
% by Ge Jin, jinwar@gmail.com
%
% JBR 10/6/17
% Modified to use age model parameterization
%

[nrow,ncol]=size(ray);
nray = nrow;
Nx=length(xnode);
Ny=length(ynode);
Nm = Nx*Ny;
xmin = min(xnode);
ymin = min(ynode);
xmax = max(xnode);
ymax = max(ynode);
dr = deg2km(mean(diff(xnode)))/10;

Nx_age = 1; % 
Nm_age = length(agebins);
xmin_age = 0;
age_min = min(agebins);
xmax_age = 1;
age_max = max(agebins);
Dx_age = xmax_age-xmin_age;
D_age = age_max-age_min;
bins_age = [1:Nm_age];
G=spalloc(nray,Nm_age,nray*Nm_age);
count_all = zeros(1,Nm_age);

Dx = xmax - xmin;
Dy = ymax - ymin;

bins=[1:Nm];

for i = 1:nray
    lat1 = ray(i,1);
    lon1 = ray(i,2);
    lat2 = ray(i,3);
    lon2 = ray(i,4);
    %r = distance(lat1,lon1,lat2,lon2)*d2r;
    [r azi] = distance(lat1,lon1,lat2,lon2);
	r = deg2km(r);

	% set segment length, 1km
%     if r<dr
%         continue;
%     end    
	Nr = floor(r/dr);
	% AGORITHEM BY W.MENKE, MATLAB BOOK CODE 12-5    

	% I use a sloppy way of computing the length of the ra
	% in each pixel.  I subdivide the ray into Nr pieces, and
	% assign each piece to exactly one pixel, the one its
	% closest to

	[lat_way,lon_way] = gcwaypts(lat1,lon1,lat2,lon2,Nr);
    [LAT_way, LON_way] = meshgrid(lat_way,lon_way);
    [XNODE, YNODE] = meshgrid(xnode,ynode);
    AGE_way = interp2(XNODE,YNODE,agenode',LAT_way,LON_way);
    
	% mid point location of segment
	xv = 0.5*(lat_way(1:Nr)+lat_way(2:Nr+1));
	yv = 0.5*(lon_way(1:Nr)+lon_way(2:Nr+1));
    [XV,YV] = meshgrid(xv,yv);
    agev = interp2(LAT_way,LON_way,AGE_way,XV,YV);
	% way-point of each ray, for small area they can be approximated 
	% linear intevals both in lat and lon;
	% faster way, or so we hope
    % calculate the array indices of all the ray pieces
    %xv = x1 + (x2-x1)*[1:Nr]'/Nr;
    %yv = y1 + (y2-y1)*[1:Nr]'/Nr;
    
    iage = 1+floor( (Nm_age-1)*(agev-age_min)/D_age);
    iage = iage(:,1);
    qvage = iage;
    
    count = hist(qvage,bins_age);
    icount = find( count~=0 );
    G(i,icount) = G(i,icount) + count(icount)*dr;
    count_all(icount) = count_all(icount)+count(icount)*dr;
    
    % plot ages and paths
    if 0
        figure(40); clf;
        subplot(1,2,1);
        hold on;
        contour(ynode,xnode,agenode,'ShowText','on');
        plot(lon1,lat1,'or');
        plot(lon2,lat2,'or');
        plot(lon_way,lat_way,'-k')
        subplot(1,2,2);
        hist(qvage,bins_age)
        pause
    end
    
end


return
end
