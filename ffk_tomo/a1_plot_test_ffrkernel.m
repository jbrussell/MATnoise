%  Plot 2-D finite-frequency kernels for surface waves using equations (11) and
% (12) from Lin & Ritzwoller (2010) doi: 10.1111/j.1365-246X.2010.04643.x
%
%
setup_parameters_tomo

% Index of period
per_ind = [2:2:20];
Mp=4; Np=4; % for plotting

% station pairs to generate kernels for
sta1 = 'WW04'; sta2 = 'EE04';

% Settings for kernels
kmode = -1; % <0 : bandwidth; =0 instantaneous; >0 Fresnel zones
nfreqs = 30; % number of frequencies to include in bandwidth average
bw_frac = 0.25; % fraction of frequency to include in bandwidth
tphase_in = []; 
% type = 'empirical'; % Requires tphase_in path
type = 'analytical'; % theoretical kernels from Lin & Ritzwoller (2010)

%%
% Set up geometry parameters
% setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
gridsize_azi = parameters.gridsize_azi;
station_list = parameters.station_list;

% Load station info
[sta.name, sta.lat, sta.lon, sta.dep] = textread(station_list,'%s %f %f %f');

fiterrtol = parameters.fiterrtol;
maxerrweight = parameters.maxerrweight;
polyfit_dt_err = parameters.polyfit_dt_err;
damp0 = parameters.damp0;
smweight0 = parameters.smweight0;
smweight0_azi = parameters.smweight0_azi;
flweight0_azi = parameters.flweight0_azi;
damp0_azi = parameters.damp0_azi;
is_wlsmooth = parameters.is_wlsmooth;
dterrtol = parameters.dterrtol;
stderrfac = parameters.stderrfac; % remove measurements with err(i) > std(err)*stderrfac
raydensetol = parameters.raydensetol;
raydensetol_azi = parameters.raydensetol_azi;
if ~is_raydensity_thresh 
    raydensetol = 1;
    raydensetol_azi = 1;
end
r = parameters.r;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
Nx = length(xnode);
Ny = length(ynode);

%%
% Initialize the xsp structure
% Xsp_path = './Xsp/';
% Xsp_path = ['../Xsp/',windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/'];
Xsp_path = ['../Xsp/',windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',num2str(N_wl),'wl_',xspdir,'/'];
xspfiles = dir([Xsp_path,'*_xsp.mat']);

disp('Looking at Xsp Files')
for ixsp = 1:length(xspfiles)
    
    temp = load([Xsp_path,xspfiles(ixsp).name]);
    xspinfo = temp.xspinfo;
    
    if ixsp ==1
        Tperiods = xspinfo.per_start;
        waxis = temp.waxis;
        twloc = temp.twloc;
        xspinfo.isgood = zeros(size(Tperiods));
        xspinfo.waxis = waxis;
        xspsum = xspinfo;
        wavelength = xspinfo.c .* xspinfo.per;
    else
        waxis = temp.waxis;
        twloc = temp.twloc;
        xspinfo.isgood = zeros(size(Tperiods));
        xspinfo.waxis = waxis;
        xspsum = [xspsum;xspinfo];
        wavelength = xspinfo.c .* xspinfo.per;
    end
    clear temp
    dep1 = sta.dep(strcmp(xspsum(ixsp).sta1,sta.name));
    dep2 = sta.dep(strcmp(xspsum(ixsp).sta2,sta.name));

    
    % 	xspinfo(ixsp).isgood = 0;
%     if xspsum(ixsp).sumerr < errlevel ...
%             && xspsum(ixsp).snr > snrtol && xspsum(ixsp).coherenum > mincoherenum
%         xspsum(ixsp).isgood = 1;
%     end
    for ip = 1:length(xspinfo.per)
        if ~is_rtolmin_wavelength && xspinfo.snr >= snr_tol && xspinfo.r >= r_tol_min && xspinfo.r <= r_tol_max && xspinfo.sumerr <= err_tol
            xspsum(ixsp).isgood(ip) = 1;
        elseif  is_rtolmin_wavelength && xspinfo.snr >= snr_tol && xspinfo.r >= wavelength(ip)*wl_fac && xspinfo.r <= r_tol_max && xspinfo.sumerr <= err_tol
            xspsum(ixsp).isgood(ip) = 1;
        end

        if isempty(cell2mat(strfind(sta.name,xspsum(ixsp).sta1))) || isempty(cell2mat(strfind(sta.name,xspsum(ixsp).sta2))) || dep1>min_dep || dep2>min_dep
            xspsum(ixsp).isgood(ip) = 0;
        end
    end
    xspsum(ixsp).isgood = logical(xspsum(ixsp).isgood .* xspsum(ixsp).isgood_wl);
    
    if rem(ixsp,500)==0
        disp(['Looking at #',num2str(ixsp),' of ',num2str(length(xspfiles))])
    end
end % end of loop ixsp'

%%

ii = 1;
% Np = length(per_ind); % columns
% Mp = 2; % rows
cmap = tomo_cmap(200);

for iper = 1:length(per_ind)
    
    period = Tperiods(per_ind(iper));
    freq=1/period;
    
    Ixsp = find(strcmp({xspsum.sta1},sta1) & strcmp({xspsum.sta2},sta2));
%     ista1 = strcmp(sta.name,sta1);
%     ista2 = strcmp(sta.name,sta2);
    
    lat1 = xspsum(Ixsp).lat1;
    lon1 = xspsum(Ixsp).lon1;
    lat2 = xspsum(Ixsp).lat2;
    lon2 = xspsum(Ixsp).lon2;
    X1 = [lon1 lat1];
    X2 = [lon2 lat2];
    
    phv = xspsum(Ixsp).c(per_ind(iper));
    xg = ynode;
    yg = xnode;
    
    switch type
    case 'analytical'
        [dtds, ~]=ffrq_2dkernel_lalo_fb(kmode,nfreqs,bw_frac,freq,phv,X1,X2, ynode, xnode);
    case 'empirical'
        sta1 = xspsum(Ixsp).sta1;
        sta2 = xspsum(Ixsp).sta2;
        fil1 = dir([tphase_in,num2str(1/freq),'/phasedelay_',sta1,'*',num2str(1/freq),'s.mat']);
        fil2 = dir([tphase_in,num2str(1/freq),'/phasedelay_',sta2,'*',num2str(1/freq),'s.mat']);
        temp = load([fil1.folder,'/',fil1.name]);
        s1 = temp.run.lalo;
        temp = load([fil2.folder,'/',fil2.name]);
        s2 = temp.run.lalo;
        [dtds, ~]=ffrq_2dkernel_lalo_fb(kmode,nfreqs,bw_frac,freq,phv,X1,X2, xg, yg,s1,s2);
    end

    % % Get the kernel (instantaneous kernel)
    %  [dtds, ier]=ffrq_2dkernel_lalo(kmode,freq, vref,X1,X2, xg, yg);
    % 
    % % Get the kernel (empirical)
    %  [dtds_emp, ier]=ffrq_2dkernel_lalo(kmode,freq, vref,X1,X2, xg, yg,s1,s2);

    [x_mesh,y_mesh,vg] = vec2mesh(xg,yg,dtds);

    %% Get bent rays extracted from kernel
    ynodek = min(ynode):0.05:max(ynode);
    xnodek = min(xnode):0.05:max(xnode);
    dr = deg2km(mean(diff(xnode)))/10;
    
    switch type
    case 'analytical'
        [dtds_ray, ~]=ffrq_2dkernel_lalo_fb(1,[],bw_frac,freq,phv,X1,X2, ynodek, xnodek);
    case 'empirical'
        [dtds_ray, ~]=ffrq_2dkernel_lalo_fb(1,[],bw_frac,freq,phv,X1,X2, ynodek, xnodek,s1,s2);
    end
    [y_meshk,x_meshk,vg_ray] = vec2mesh(ynodek,xnodek,dtds_ray*-1);

    [r azi] = distance(X1(2), X1(1), X2(2), X2(1));
    r = deg2km(r);

    Nr = floor(r/dr);
    ray = [X1(2) X1(1) X2(2) X2(1)];
    ray_waypts = trace_kernel(ray,Nr,y_meshk,x_meshk,vg_ray);
    lat_way = ray_waypts(:,1);
    lon_way = ray_waypts(:,2);

    %%
    figure(11)
    if iper == 1
        clf;
        set(gcf,'color','w','Position',[197    69   898   631]);
    end

    subplot(Mp,Np,iper);
    Nlvls = 80;
    lvls = linspace(-prctile(abs(vg(:)),99.9),prctile(abs(vg(:)),99.9),Nlvls);
    [~,h]=contourf(x_mesh,y_mesh,vg,lvls); hold on; %reshape(dtds,mmx,mmy)' 
    % contour(x_mesh,y_mesh,vg,[0.001 0.001],'-k')
    set(gca,'fontsize',16,'linewidth',1.5)
    set(h,'linestyle','none')
    % colormap('jet')
    colormap(redblue(Nlvls))
    plot(sta.lon,sta.lat,'ok','markerfacecolor','k','markersize',3); hold on;
%     plot(run.sta.lon,run.sta.lat,'k^')
    plot([X1(1) X2(1)],[X1(2) X2(2)],'k^','markerfacecolor','y')
    [lat,lon] = gcwaypts(X1(2),X1(1),X2(2),X2(1));
    plot(lon,lat,'k-','linewidth',1.5)
    % colorbar;
    title(type);
    caxis([min(lvls) max(lvls)]);
    axis equal 
    % xlim([75 225])
    title([num2str(period),' s']);


    %%
    figure(9);
    if iper == 1
        clf;
        set(gcf,'color','w','Position',[197    69   898   631]);
    end
    subplot(Mp,Np,iper);
    
    [~,I] = min(abs(yg - (X1(2)+X2(2))/2));
    plot(xg,vg(I,:),'-k','linewidth',2); hold on;
%     plot(xg,vg_emp(I,:),'-r','linewidth',2); hold on;
    xlabel('X');
%     legend({'Analytical','Empirical'},'box','off');
    set(gca,'LineWidth',1.5,'FontSize',12);
    % ylim([-0.04 0.04]);
    title([num2str(period),' s']);

    %%
%     figure(10); clf;
%     [~,h]=contourf(run.out.xpts_mesh,run.out.ypts_mesh,vs_mesh,100); %reshape(dtds,mmx,mmy)' 
%     set(gca,'fontsize',12)
%     set(h,'linestyle','none')
%     colormap(flip(jet));
%     hold on
%     plot(run.out.xsta,run.out.ysta,'k^','markerfacecolor','w','markersize',10)
%     % colorbar;
%     colorbar;
%     axis equal
    
    
end
