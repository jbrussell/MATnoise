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

function [dr_ray] = gc_raydr_km(lat_way,lon_way)
    % Calculate dr vector in units of km for lat-lon waypoints using great circle
    % approximations along each segment. (If assume straight rays, can
    % accumulate errors of ~20% !)
    % JBR 5/8/2020
    %
    dr_ray = distance(lat_way(1:end-1),lon_way(1:end-1),...
                             lat_way(2:end),lon_way(2:end),referenceEllipsoid('GRS80'))/1000;
end
