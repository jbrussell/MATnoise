% Extract dispersion estimates at the pixel nearest to each station.
%
% jbrussell - 11/2025

clear all;

% Parameters
Npixel_dist = 2; % Station must be within this many pixels of the nearest good pixel

% Path to raytomo dispersion .mat file
matname = 'raytomo_9_35s_1wl_ZZ_0S_LRT.mat';
outpath = './phvmaps/';
path2mat = [outpath,'/',matname];
temp = load(path2mat);
parameters = temp.parameters;
gridsize = parameters.gridsize;
raytomo = temp.raytomo;
xi = temp.xi;
yi = temp.yi;
xnode = temp.xnode;
ynode = temp.ynode;
lalim = parameters.lalim;
lolim = parameters.lolim;

% Periods to plot;
per_ind = [5 20 30];
r = 0.03; % for plotting only

% Load station info
station_list = parameters.station_list;
[sta.net, sta.name, sta.lat, sta.lon, sta.dep] = textread(station_list,'%s %s %f %f %f');
netname = {};
for ista = 1:length(sta.name)
    netname{ista,1} = [sta.net{ista},'.',sta.name{ista}];
end
sta.name = netname;

%% Loop through stations and extract dispersion curves at nearest pixel

periods = [raytomo.period];
stadisp = [];
for ista = 1:length(sta.name)
    rmat = distance(sta.lat(ista),sta.lon(ista),xi,yi,referenceEllipsoid('GRS80'))/1000;
    
    phv = zeros(size(periods));
    phv_std = zeros(size(periods));
    for iper = 1:length(periods)
        phv_map = raytomo(iper).GV;
        phv_map_std = raytomo(iper).GV_std;
        
        % Get distance matrix and set bad pixel locations to nan
        rmat_nan = rmat;
        rmat_nan(isnan(phv_map)) = nan;
        % Find nearest good pixel to array
        [rmin,imin] = min(rmat_nan(:));
        imin = imin(1);
        rmin = rmin(1);

        % If too far away from a good pixel, skip
        if rmin/deg2km(gridsize) > Npixel_dist
            display([sta.name{ista},' ',num2str(periods(iper)),'s is too far from a good pixel'])
            phv(iper) = nan;
            phv_std(iper) = nan;
            continue
        end
        
        phv(iper) = phv_map(imin);
        phv_std(iper) = phv_map_std(imin);

    end

    stadisp(ista).sta = sta.name{ista};
    stadisp(ista).lat = sta.lat(ista);
    stadisp(ista).lon = sta.lon(ista);
    stadisp(ista).z = sta.dep(ista);
    stadisp(ista).periods = periods;
    stadisp(ista).phv = phv;
    stadisp(ista).phv_std = phv_std;
end

save([outpath,'/',strrep(matname,'raytomo','stadisp')],'stadisp')

%% Plot map of station dispersion

Mp = 2; Np = 3;
fig33 = figure(33);
% set(gcf,'position',[94     1   599   704]);
set(gcf,'position',[1    1   1244   704],'color','w');
clf
ii = 0;
cmap = tomo_cmap(200);
sta_lats = [stadisp.lat];
sta_lons = [stadisp.lon];
for ip=per_ind
    ii = ii + 1;
    subplot(Mp,Np,ii)
    ax = worldmap(lalim, lolim);
    setm(gca,'MapProjection','mercator','FLineWidth',1.5,'FontSize',13)
    tightmap
%     set(ax, 'Visible', 'on')
    set(gcf,'color','w')
    setm(gca,'FFaceColor',[0.9 0.9 0.9])
%     set(gca,'Color',[0.7 0.7 0.7])
%     surfacem(xi,yi,raytomo(ip).GV);
    avgv = nanmean(raytomo(ip).GV(:));
    levels = linspace(avgv*(1-r), avgv*(1+r),100);
    contourfm(xi,yi,raytomo(ip).GV,levels,'edgecolor','none');
%     drawlocal
    title([num2str((periods(ip))),' s'],'fontsize',15)
    % text(0.05,0.85,[num2str(round(Tperiods(ip))),' s'],'fontsize',15,'fontweight','bold','Units','normalized','HorizontalAlignment','left');
    caxis([avgv*(1-r) avgv*(1+r)])
    cb = colorbar;
    ylabel(cb,'Phase Velocity (km/s)','fontsize',15);
    posax = get(ax,'Position');
    pos=get(cb,'Position');
    set(cb,'Position',[pos(1)+0.03 pos(2) pos(3)*0.8 pos(4)],'linewidth',1.5,'fontsize',15);
    set(gca,'Position',[posax(1) posax(2:4)],'fontsize',15);
%     colormap(seiscmap)
%     rbc = flip(redbluecmap);
%     rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
%     colormap(rbc);
    colormap(cmap);    
    hold on;
    plotm(sta_lats,sta_lons,'ok','markerfacecolor',[0 0 0]);
%     [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);

    
    % gather values at each period
    phv_stas = [];
    for ista = 1:length(stadisp)
        phv_stas(ista) = stadisp(ista).phv(ip);
%         phv_std_stas(ista) = stadisp(ista).phv_std(ip);
    end

    subplot(Mp,Np,ii+Np)
    ax = worldmap(lalim, lolim);
    setm(gca,'MapProjection','mercator','FLineWidth',1.5,'FontSize',13)
    tightmap
%     set(ax, 'Visible', 'on')
    set(gcf,'color','w')
    setm(gca,'FFaceColor',[0.9 0.9 0.9])
%     set(gca,'Color',[0.7 0.7 0.7])
%     surfacem(xi,yi,raytomo(ip).GV);
    avgv = nanmean(raytomo(ip).GV(:));
%     drawlocal
    title([num2str((periods(ip))),' s'],'fontsize',15)
    % text(0.05,0.85,[num2str(round(Tperiods(ip))),' s'],'fontsize',15,'fontweight','bold','Units','normalized','HorizontalAlignment','left');
    caxis([avgv*(1-r) avgv*(1+r)])
    cb = colorbar;
    ylabel(cb,'Phase Velocity (km/s)','fontsize',15);
    posax = get(ax,'Position');
    pos=get(cb,'Position');
    set(cb,'Position',[pos(1)+0.03 pos(2) pos(3)*0.8 pos(4)],'linewidth',1.5,'fontsize',15);
    set(gca,'Position',[posax(1) posax(2:4)],'fontsize',15);
%     colormap(seiscmap)
%     rbc = flip(redbluecmap);
%     rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
%     colormap(rbc);
    colormap(cmap);    
    hold on;
    scatterm(sta_lats,sta_lons,80,phv_stas,'o','filled','markeredgecolor',[0 0 0]);
%     [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
end
