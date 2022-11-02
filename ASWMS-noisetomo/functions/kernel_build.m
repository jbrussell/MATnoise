function G = kernel_build(ray, xnode, ynode)
%input: ray: nray*4 matrix; each row: [x1,y1,x2,y2]
%       xnode and ynode are grid axis for lat and long
% Output: G is the kernal for Vx and Vy
% g(i,l) = r(i,l) is length of i th ray in lth pixel;
% build data kernel
% written by Yang Zha, modified for fitting phase gradient surface 
% by Ge Jin, jinwar@gmail.com
%
% jbrussell 11/20: Use more accurate method for counting portion of ray within 
% each grid cell

[nrow,ncol]=size(ray);
nray = nrow;
Nx=length(xnode);
Ny=length(ynode);
Nm = Nx*Ny;
xmin = min(xnode);
ymin = min(ynode);
xmax = max(xnode);
ymax = max(ynode);

Dx = xmax - xmin;
Dy = ymax - ymin;
G=spalloc(nray,Nm*2,2*nray*Nx); % for each ray, maximum number of pixels to be sampled is 2*Nx
%G=zeros(nray,Nm);
bins=[1:Nm];

for i = 1:nray
    dr = deg2km(mean(diff(xnode)))/1e3;
    lat1 = ray(i,1);
    lon1 = ray(i,2);
    lat2 = ray(i,3);
    lon2 = ray(i,4);
    %r = distance(lat1,lon1,lat2,lon2)*d2r;
    [r, azi] = distance(lat1,lon1,lat2,lon2,referenceEllipsoid('GRS80'));
    r = r/1000;

	% set segment length
    if r<dr
        continue;
    end    
	Nr = round(r/dr);
    % AGORITHEM BY W.MENKE, MATLAB BOOK CODE 12-5    

    % I use a sloppy way of computing the length of the ra
    % in each pixel.  I subdivide the ray into Nr pieces, and
    % assign each piece to exactly one pixel, the one its
    % closest to

    [lat_way,lon_way] = gcwaypts(lat1,lon1,lat2,lon2,Nr);
    [dr, azi_ray] = gc_raydr_km(lat_way,lon_way);
    % mid point location of segment
    xv = 0.5*(lat_way(1:Nr)+lat_way(2:Nr+1));
    yv = 0.5*(lon_way(1:Nr)+lon_way(2:Nr+1));
    % way-point of each ray, for small area they can be approximated by 
    % linear intevals both in lat and lon;
    if( 0 ) % slow but sure way
        for ir = 1:Nr
            x = x1 + (x2-x1)*i/Nr;
            y = y1 + (y2-y1)*i/Nr;
            ix = 1+floor( Nx*(x-xmin)/Dx );
            iy = 1+floor( Ny*(y-ymin)/Dy );
            q = (ix-1)*Ny + iy;
            G(k,q) = G(k,q) + dr;
        end
    else % faster way, or so we hope
        % calculate the array indices of all the ray pieces
        %xv = x1 + (x2-x1)*[1:Nr]'/Nr;
        %yv = y1 + (y2-y1)*[1:Nr]'/Nr;
        
        ixv = 1+floor( (Nx-1)*(xv-xmin)/Dx );
        iyv = 1+floor( (Ny-1)*(yv-ymin)/Dy );
        qv = (ixv-1)*Ny + iyv;
        % sum binned dr values
        qv = sort(qv);    
        drq = zeros(size(unique(qv)'));
        azq = zeros(size(unique(qv)'));
        ii = 0;
        for iq = unique(qv)'
            ii = ii+1;
            drq(ii) = sum(dr(qv==iq));
            azq(ii) = angmean(azi_ray(qv==iq)*pi/180)*180/pi;
        end
        % now count of the ray segments in each pixel of the
        % image, and use the count to increment the appropriate
        % element of G.  The use of the hist() function to do
        % the counting is a bit weird, but it seems to work
        count=hist(qv,bins); 
        icount = find( count~=0 );
        if length(drq) ~= length(icount)
            % Something is wrong...
            continue
        end
        % G(i,2*icount-1) = G(i,2*icount-1) + count(icount)*dr*cosd(azi);
        % G(i,2*icount) = G(i,2*icount) + count(icount)*dr*sind(azi);
        G(i,2*icount-1) = G(i,2*icount-1) + drq.*cosd(azq);
        G(i,2*icount) = G(i,2*icount) + drq.*sind(azq);
    end
end


return

%%
    function [dr_ray,azi_ray] = gc_raydr_km(lat_way,lon_way)
        % Calculate dr vector in units of km for lat-lon waypoints using great circle
        % approximations along each segment. (If assume straight rays, can
        % accumulate errors of ~20% !)
        % JBR 5/8/2020
        %
        [dr_ray,azi_ray] = distance(lat_way(1:end-1),lon_way(1:end-1),...
                                 lat_way(2:end),lon_way(2:end),referenceEllipsoid('GRS80'));
        dr_ray = dr_ray / 1000;
    end
    
end
