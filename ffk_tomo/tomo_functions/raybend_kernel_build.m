function G = raybend_kernel_build(ray, xnode, ynode,dtds,xnodek,ynodek,freq)
%input: ray: nray*4 matrix; each row: [x1,y1,x2,y2]
%       xnode and ynode are grid axis for lat and long
% Output: G is the kernal for Vx and Vy
% g(i,l) = r(i,l) is length of i th ray in lth pixel;
% build data kernel
% written by Yang Zha, modified for ray theory for ambient noise analysis
% by Ge Jin, jinwar@gmail.com
%
% JBR 5/2020: Trace center of banana-donut finite frequency kernel to extract the
% bent ray path and build kernel for ray-based inversion.
% global windir
global kernel_path
issave = 1; % save kernels?
% kernel_path = ['./SEM2D_FFK_save/',windir,'/']; % path to saved kernels
type = 'raybend';

[nrow,ncol]=size(ray);
nray = nrow;
Nx=length(xnode);
Ny=length(ynode);
Nm = Nx*Ny;
xmin = min(xnode);
ymin = min(ynode);
xmax = max(xnode);
ymax = max(ynode);

% Check if kernels have already been calculated
filename = [kernel_path,type,'_',num2str(1./freq),'s_',num2str(Nx),'x',num2str(Ny),'.mat'];
if exist(filename)
    temp = load(filename);
    G = temp.G;
    return
end

Dx = xmax - xmin;
Dy = ymax - ymin;
G=spalloc(nray,Nm,nray*Nx); % for each ray, maximum number of pixels to be sampled is 2*Nx
%G=zeros(nray,Nm);
bins=[1:Nm];

for i = 1:nray
    dr = deg2km(mean(diff(xnode)))/10;
    lat1 = ray(i,1);
    lon1 = ray(i,2);
    lat2 = ray(i,3);
    lon2 = ray(i,4);
    
    [y_mesh,x_mesh,vg] = vec2mesh(ynodek,xnodek,dtds(i,:));
    
    %r = distance(lat1,lon1,lat2,lon2)*d2r;
    [r, azi] = distance(lat1,lon1,lat2,lon2,referenceEllipsoid('GRS80'));
    r = r/1000;

	% set segment length, 1km
    if r<dr
        continue;
    end    
	Nr = floor(r/dr);
	% AGORITHEM BY W.MENKE, MATLAB BOOK CODE 12-5   
    % JBR 5/2020: Modify to use bent rays

	% I use a sloppy way of computing the length of the ra
	% in each pixel.  I subdivide the ray into Nr pieces, and
	% assign each piece to exactly one pixel, the one its
	% closest to

% 	[lat_way,lon_way] = gcwaypts(lat1,lon1,lat2,lon2,Nr);
    
    % Trace curved rays
    ray_waypts = trace_kernel(ray(i,:),Nr,y_mesh,x_mesh,vg);
    lat_way = ray_waypts(:,1);
    lon_way = ray_waypts(:,2);
    % recalculate Nr and dr
    Nr = length(lon_way)-1;
    dr = gc_raydr_km(lat_way,lon_way);
%     dr = deg2km(sqrt(diff(lat_way).^2+diff(lon_way).^2));
%     dr = mean(dr);
    
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
    
    % sum binned dr values
    qv = sort(qv);    
    drq = zeros(size(unique(qv)'));
    ii = 0;
    for iq = unique(qv)'
        ii = ii+1;
        drq(ii) = sum(dr(qv==iq));
    end
    
    % now count of the ray segments in each pixel of the
    % image, and use the count to increment the appropriate
    % element of G.  The use of the hist() function to do
    % the counting is a bit weird, but it seems to work
    count=hist(qv,bins); 
    icount = find( count~=0 );
%     G(i,icount) = G(i,icount) + count(icount)*dr;
    G(i,icount) = G(i,icount) + drq;
    
    if 0      
        % PLOT RAYS AND RAY DENSITY
        figure(39);
        if i == 1
            clf;
        end
        hold on;
%         contour(ynode,xnode,agenode,'ShowText','on');
        plot(lon1,lat1,'or');
        plot(lon2,lat2,'or');
        plot(lon_way,lat_way,'-k')
%         if i == nray
            for ii=1:Nx
                for jj=1:Ny 
                    n=Ny*(ii-1)+jj;
                    raydense(ii,jj) = sum(G(:,n));
                end
            end
%             hi = imagesc(ynode,xnode,raydense); % ynode=lon, xnode=lat
%             set(hi,'alphadata',.5);
            hi = surface(ynode,xnode,raydense,'FaceAlpha',0.5,'Linestyle','none'); % ynode=lon, xnode=lat
%         end
        drawnow;
        pause
    end
end

if issave
    if ~exist(kernel_path)
        mkdir(kernel_path)
    end
    save(filename,'G','xnode','ynode');
end

return

%%
    function [ray_waypts] = trace_kernel(ray,Nr,x_mesh,y_mesh,vg)
        % Trace central well of finite frequency kernel to extract
        % equivalent ray path. This is done by taking transections through
        % finite frequency kernel that are perpendicular to the
        % great-circle path connecting station pairs. The findpeaks()
        % function is used to identify the ray path. Directly near the
        % stations the finite frequency kernels are numerically unstable to
        % can't trust findpeaks(), so I've just discarded 20% of the points
        % on the edges of the ray while keeping the 2 endpoints.
        %
        % vg: finite frequency kernel with amplitude flipped
        %
        % JBR 5/7/2020
        
        la1 = ray(1,1);
        lo1 = ray(1,2);
        la2 = ray(1,3);
        lo2 = ray(1,4);
        ray_lat = [];
        ray_lon = [];
        Nr = floor(Nr*1.2);
        
        gcpath = gcwaypts(la1,lo1,la2,lo2,Nr);
        [~,az] = distance(la1,lo1,la2,lo2);
        
        Nskip = floor(Nr*0.2); % point to skip near stations due to numerical instabilities in the kernels...
        for ir=(1+Nskip):(size(gcpath,1)-Nskip)
            ypt = gcpath(ir,1);
            xpt = gcpath(ir,2);
            
            % Define transect perpendicular to the great circle path
            % connecting the two stations
            R = 3; % Length of transect in degrees
            x1 = xpt+R*sind(az+90); y1 = ypt+R*cosd(az+90);
            x2 = xpt+R*sind(az-90); y2 = ypt+R*cosd(az-90);
            x = linspace(x1,x2,250);
            y = linspace(y1,y2,250);
            
            igood = (y>=min(y_mesh(:,1)) & y<=max(y_mesh(:,1))) & ...
                    (x>=min(x_mesh(1,:)) & x<=max(x_mesh(1,:)));
            x = x(igood);
            y = y(igood);

            % Interpolate kernel along transect
            kern = interp2(x_mesh,y_mesh,vg,x,y);
            
%             figure(41);
%             Nlvls = 80;
%             lvls = linspace(-prctile(abs(vg(:)),99.9),prctile(abs(vg(:)),99.9),Nlvls);
%             contourf(x_mesh,y_mesh,vg,lvls,'linestyle','none'); hold on; %reshape(dtds,mmx,mmy)' 
%             plot(x,y,'-r');
            
            [pks,locs] = findpeaks(kern);
            if ~isempty(pks)
                [~,Imax] = max(pks);
                ray_lat = [ray_lat; y(locs(Imax))];
                ray_lon = [ray_lon; x(locs(Imax))];
            end
            
%             figure(41); clf;
%             findpeaks(kern)

        end
        ray_lat = [la1; ray_lat; la2];
        ray_lon = [lo1; ray_lon; lo2];
%         ray_lat = smooth(ray_lat,Nr/10);
%         ray_lon = smooth(ray_lon,Nr/10);
        ray_lat = smooth(ray_lat,Nr/3);
        ray_lon = smooth(ray_lon,Nr/3);
        ray_waypts = [ray_lat, ray_lon];
        
        if 0
            figure(40); clf;
            [law,low] = gcwaypts(la1,lo1,la2,lo2,Nr);
            dray = sum(gc_raydr_km(law,low));
            draybend = sum(gc_raydr_km(ray_lat,ray_lon));
            dgc = distance(la1,lo1,la2,lo2,referenceEllipsoid('GRS80'))/1000;
            Nlvls = 80;
            lvls = linspace(-prctile(abs(vg(:)),99.9),prctile(abs(vg(:)),99.9),Nlvls);
            contourf(x_mesh,y_mesh,vg,lvls,'linestyle','none'); hold on; %reshape(dtds,mmx,mmy)' 
            plot(low,law,'-k');
            plot(ray_lon,ray_lat,'-r');
%             title([num2str((draybend-dray)./dray*100),'%']);
            title([num2str((draybend-dgc)./dgc*100),'%']);
            % plot(
            caxis([min(lvls) max(lvls)]);
            pause;
        end
    end
%%
    function [dr_ray] = gc_raydr_km(lat_way,lon_way)
        % Calculate dr vector in units of km for lat-lon waypoints using great circle
        % approximations along each segment. (If assume straight rays, can
        % accumulate errors of ~20% !)
        % JBR 5/8/2020
        %
        dr_ray = distance(lat_way(1:end-1),lon_way(1:end-1),...
                                 lat_way(2:end),lon_way(2:end),referenceEllipsoid('GRS80'))/1000;
    end
    
end
