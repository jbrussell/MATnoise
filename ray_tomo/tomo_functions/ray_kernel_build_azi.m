function [G, G_hits] = ray_kernel_build_azi(ray, xnode, ynode)
%input: ray: nray*4 matrix; each row: [x1,y1,x2,y2]
%       xnode and ynode are grid axis for lat and long
% Output: G is the kernal for Vx and Vy
% g(i,l) = r(i,l) is length of i th ray in lth pixel;
% build data kernel
% written by Yang Zha, modified for ray theory for ambient noise analysis
% by Ge Jin, jinwar@gmail.com
%
% JBR 7/22/18 - modified to include 2D azimuthal anisotropy
% g(i,l) = r(i,l)*[cosd(2*azi), sind(2*azi)]
%
% When plotting using surface(), the color of each square corresponds to
% the value of the lower left corner. Therefore, in order to plot fast
% directions in map view, need to add dy/2 to y-coordinate and dx/2 to
% x-coordinate
%
% o: Where grid places fast direction
% *: Where we want to plot fast direction
%
%             -----------
%            |           |
%            |           |
%            |     *     | ^
%            |           | |
%            |           | | dy/2
%            o-----------  v
%            <---->
%             dx/2
%
isplot = 0;

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

Dx = xmax - xmin;
Dy = ymax - ymin;
Gc=spalloc(nray,Nm,nray*Nx); % for each ray, maximum number of pixels to be sampled is 2*Nx
Gs=spalloc(nray,Nm,nray*Nx); % for each ray, maximum number of pixels to be sampled is 2*Nx
G_hits=spalloc(nray,Nm,nray*Nx); % for each ray, maximum number of pixels to be sampled is 2*Nx
%G=zeros(nray,Nm);
bins=[1:Nm];

for i = 1:nray
    lat1 = ray(i,1);
    lon1 = ray(i,2);
    lat2 = ray(i,3);
    lon2 = ray(i,4);
    %r = distance(lat1,lon1,lat2,lon2)*d2r;
    [r, azi] = distance(lat1,lon1,lat2,lon2,referenceEllipsoid('GRS80'));
    r = r/1000;

	% set segment length, 1km
    if r<dr
        continue;
    end    
	Nr = floor(r/dr);
	% AGORITHEM BY W.MENKE, MATLAB BOOK CODE 12-5    

	% I use a sloppy way of computing the length of the ra
	% in each pixel.  I subdivide the ray into Nr pieces, and
	% assign each piece to exactly one pixel, the one its
	% closest to

	[lat_way,lon_way] = gcwaypts(lat1,lon1,lat2,lon2,Nr);
	% mid point location of segment
	xv = 0.5*(lat_way(1:Nr)+lat_way(2:Nr+1));
	yv = 0.5*(lon_way(1:Nr)+lon_way(2:Nr+1));
	% way-point of each ray, for small area they can be approximated by 
	% linear intevals both in lat and lon;
	% faster way, or so we hope
    % calculate the array indices of all the ray pieces
    %xv = x1 + (x2-x1)*[1:Nr]'/Nr;
    %yv = y1 + (y2-y1)*[1:Nr]'/Nr;
    
    ixv = 1+floor( (Nx-1)*(xv-xmin)/Dx );
    iyv = 1+floor( (Ny-1)*(yv-ymin)/Dy );
    qv = (ixv-1)*Ny + iyv;
    % now count of the ray segments in each pixel of the
    % image, and use the count to increment the appropriate
    % element of G.  The use of the hist() function to do
    % the counting is a bit weird, but it seems to work
    count=hist(qv,bins); 
    icount = find( count~=0 );
    Gc(i,icount) = Gc(i,icount) + count(icount)*dr*cosd(2*azi);
    Gs(i,icount) = Gs(i,icount) + count(icount)*dr*sind(2*azi);
    G_hits(i, icount) = G_hits(i,icount) + count(icount)*dr;
    
    % plot paths
    if isplot
%         figure(40);
%         if i == 1
%             clf;
%         end
%         subplot(1,2,1);
%         hold on;
% %         contour(ynode,xnode,agenode,'ShowText','on');
%         plot(lon1,lat1,'or');
%         plot(lon2,lat2,'or');
%         plot(lon_way,lat_way,'-k')
%         subplot(1,2,2);
%         hist(qv,bins)
        
        % PLOT RAYS AND RAY DENSITY
        figure(40);
        set(gcf,'position',[55   285   919   420]);
        
        %%%%%% PLOT OLD VERSION WITH VOXELS DEFINED BY CORNERS
        subplot(1,2,1);
        if i == 1
            clf;
            ax1 = subplot(1,2,1);
        end
        hold on;
%         contour(ynode,xnode,agenode,'ShowText','on');
        plot(ax1,lon1,lat1,'or');
        plot(ax1,lon2,lat2,'or');
        hl(i) = plot(ax1,lon_way,lat_way,'-k');
        if i == nray
            for ii=1:Nx
                for jj=1:Ny 
                    n=Ny*(ii-1)+jj;
                    raydense(ii,jj) = sum(G_hits(:,n));
                end
            end
%             hi = imagesc(ynode,xnode,raydense); % ynode=lon, xnode=lat
%             set(hi,'alphadata',.5);
            hi1 = surface(ax1,ynode,xnode,raydense,'FaceAlpha',0.5,'Linestyle','none'); % ynode=lon, xnode=lat
        end
%         pause
        
        %%%%%% PLOT NEW VERSION WITH VOXELS DEFINED BY MIDPOINT
        ax2 = subplot(1,2,2);
        hold on;
%         contour(ynode,xnode,agenode,'ShowText','on');
        plot(ax2,lon1,lat1,'or');
        plot(ax2,lon2,lat2,'or');
        hl(i) = plot(ax2,lon_way,lat_way,'-k');
        if i == nray
            for ii=1:Nx
                for jj=1:Ny 
                    n=Ny*(ii-1)+jj;
                    raydense(ii,jj) = sum(G_hits(:,n));
                end
            end
            hi2 = imagesc(ax2,ynode,xnode,raydense); % ynode=lon, xnode=lat
            set(hi2,'alphadata',.5);
            set(ax2,'XLim',ax1.XLim,'YLim',ax1.YLim);
%             ax2.XLim = ax1.XLim;
%             ax2.YLim = ax1.YLim;
        end
    end
end
G = [Gc, Gs];

return
end
