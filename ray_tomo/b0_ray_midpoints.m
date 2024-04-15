% Color interstation rays by phase velocity
%
% https://github.com/jbrussell
clear; close all;

%======================= PARAMETERS =======================%
setup_parameters_tomo;
comp = parameters.comp; % = {'ZZ'};
xspdir = parameters.xspdir; % = 'phv_dir';
windir = parameters.windir; % = 'window3hr'; 
N_wl = parameters.N_wl;
frange = parameters.frange; % = [1/10 1/5]; % [Hz]

% QC parameters
snr_tol = parameters.snr_tol; % = 3; % minimum signal-to-noise
r_tol_min = parameters.r_tol_min; % = 90; % [km] minimum station separation
r_tol_max = parameters.r_tol_max; % = 600; % [km] maximum station separation
err_tol = parameters.err_tol; % = 0.5; % maximum misfit of bessel fit between observed and synthetic

%==========================================================%
%%
% figure output path
phv_fig_path = ['./figs/',windir,'/fullStack/raytomo/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/'];
if ~exist(phv_fig_path)    
    mkdir(phv_fig_path);
end

% Set up geometry parameters
setup_parameters;
station_list = parameters.station_list;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx = length(xnode);
Ny = length(ynode);

% Load station info
[sta.name, sta.lat, sta.lon, sta.dep] = textread(station_list,'%s %f %f %f');

fiterrtol = parameters.fiterrtol;
maxerrweight = parameters.maxerrweight;
polyfit_dt_err = parameters.polyfit_dt_err;
smweight0 = parameters.smweight0;


dterrtol = parameters.dterrtol;

raydensetol = parameters.raydensetol;
r = parameters.r;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
Nx = length(xnode);
Ny = length(ynode);

% read in bad station list, if existed
if exist('badsta.lst')
    badstnms = textread('badsta.lst','%s');
    badstaids = find(ismember({stainfo.staname},badstnms));
    disp('Found Bad stations:')
    disp(badstnms)
end

% Initialize the xsp structure
% Xsp_path = ['../Xsp/',windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/'];
Xsp_path = ['../Xsp/',windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',num2str(N_wl),'wl_',xspdir,'/'];
xspfiles = dir([Xsp_path,'*_xsp.mat']);

disp('Looking at Xsp Files')
for ixsp = 1:length(xspfiles)
    
    temp = load([Xsp_path,xspfiles(ixsp).name]);
    xspinfo = temp.xspinfo;
    
    if ixsp ==1
%         Tperiods = (2*pi)./temp.twloc;
        Tperiods = xspinfo.per_start;
        waxis = temp.waxis;
        twloc = temp.twloc;
        xspinfo.isgood = 0;
        xspsum = xspinfo;
    else
        xspinfo.isgood = 0;
        xspsum = [xspsum;xspinfo];
    end
    clear temp

    if xspinfo.snr >= snr_tol && xspinfo.r >= r_tol_min && xspinfo.r <= r_tol_max && xspinfo.sumerr <= err_tol
        xspsum(ixsp).isgood = 1;
    end
    
    if rem(ixsp,500)==0
        disp(['Looking at #',num2str(ixsp),' of ',num2str(length(xspfiles))])
    end
end % end of loop ixsp'
%%
% Loop through periods
for ip=1:length(Tperiods)
    disp(' ');
    disp(['Inverting Period: ',num2str(Tperiods(ip))]);
    clear rays dt fiterr mat phaseg err raydense dist snrs phv phv_cor
    raynum = 0;
    
    for ixsp = 1:length(xspsum)
        if xspsum(ixsp).isgood ==0
            continue;
        end
        
        raynum = raynum+1;
        rays(raynum,1) = xspsum(ixsp).lat1;
        rays(raynum,2) = xspsum(ixsp).lon1;
        rays(raynum,3) = xspsum(ixsp).lat2;
        rays(raynum,4) = xspsum(ixsp).lon2;
        
        dist(raynum) = deg2km(distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4)));
        snrs(raynum) = xspsum(ixsp).snr;
        
        % Midpoint of rays
        midlat(raynum) = (rays(raynum,1)+rays(raynum,3))./2;
        midlon(raynum) = (rays(raynum,2)+rays(raynum,4))./2;
        
        % Phase velocity of ray
        phv(raynum) = xspsum(ixsp).c(ip); % km/s
        
        % Azimuth of rays
        [~,azi(raynum)]=distance(xspsum(ixsp).lat1,xspsum(ixsp).lon1,xspsum(ixsp).lat2,xspsum(ixsp).lon2);
        if azi(raynum) > 180
            azi(raynum) = azi(raynum) - 360;
        end
        
    end
    dat(ip).rays = rays;
    dat(ip).dist = dist;
%     dat(ip).dt = dt; 
    dat(ip).midlat = midlat;
    dat(ip).midlon = midlon;
    dat(ip).phv = phv;
%     dat(ip).phv_cor = phv_cor;
    dat(ip).period = Tperiods(ip);
    dat(ip).snrs = snrs;
end

%% Plot interstation rays colored by phase velocity

Mp = 4; Np = 5;

fig20 = figure(20);
set(gcf,'position',[1    1   1244   704]);
clf
rbc = flip(redbluecmap);
% rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
for ip=1:length(Tperiods)
subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')

    if isempty(find(~isnan(dat(ip).phv)))
        disp(['All nan values for ',num2str(Tperiods(ip)),'... skipping'])
        continue
    end
    
%     avgv = nanmean(dat(ip).phv_cor);
    avgv = nanmean(dat(ip).phv);
    vels = linspace(avgv*(1-r),avgv*(1+r),size(rbc,1));
    clrs = [];
    for ixsp = 1:length(dat(ip).phv)
        lat1 = dat(ip).rays(ixsp,1);
        lon1 = dat(ip).rays(ixsp,2);
        lat2 = dat(ip).rays(ixsp,3);
        lon2 = dat(ip).rays(ixsp,4);
%         [~,iclr] = min(abs(vels - dat(ip).phv_cor(ixsp)));
        [~,iclr] = min(abs(vels - dat(ip).phv(ixsp)));
%         plotm([lat1 lat2],[lon1 lon2],dat(ip).phv(ixsp),'color',rbc(iclr,:),'linewidth',1.5);
        clrs(ixsp,:) = rbc(iclr,:);
        hold on;
    %     drawlocal
    end
    h = plotm([dat(ip).rays(:,1) dat(ip).rays(:,3)]',[dat(ip).rays(:,2) dat(ip).rays(:,4)]','linewidth',1.5);
    set(h,{'color'},num2cell(clrs,2));
    title([num2str(Tperiods(ip))],'fontsize',15)
    colorbar
%     colormap(seiscmap)
    colormap(rbc);
    caxis([vels(1) vels(end)]);
    
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
%     [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
    drawnow;
end
save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_rays.pdf'],fig20,1000);

%% Plot interstation mid points colored by phase velocity

fig22 = figure(22); % without azi corr
set(gcf,'position',[1    1   1244   704]);
clf
for ip=1:length(Tperiods)
    subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')

    if isempty(find(~isnan(dat(ip).phv)))
        disp(['All nan values for ',num2str(Tperiods(ip)),'... skipping'])
        continue
    end

%     surfacem(xi,yi,raytomo(ip).err);
%     drawpng
% scatterm(dat(ip).midlat,dat(ip).midlon,30,dat(ip).phv,'filled'); 
scatterm(dat(ip).midlat,dat(ip).midlon,60,dat(ip).phv,'filled'); 
% plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);

hold on;
% drawlocal
title([num2str(Tperiods(ip))],'fontsize',15)

% avgv = nanmean(dat(ip).phv_cor);
avgv = nanmean(dat(ip).phv);
caxis([avgv*(1-r) avgv*(1+r)])
colorbar
% colormap(seiscmap)
% colormap(flip(jet))
rbc = flip(redbluecmap);
% rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
colormap(rbc);

% [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
end

% save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_ray_midpts_azicorr.pdf'],fig21,1000);
save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_ray_midpts.pdf'],fig22,1000);