% Script to do the ray theory tomography based on the ambient noise measurement
% Written by Ge Jin, jinwar@gmail.com
% Nov 2012
%
% Modified by NJA, April 2016
%
% Modified by JBR, 10/3/17
% Added azimuthally anisotropic terms to the inversion. Now inverts for
% 2theta and 4theta coefficients (averaged over the entire array) jointly
% with the 2D isotropic part.
%
% phv(theta,freq) = phv_iso(freq) + Ac2(freq)*cos(2*theta) + As2(freq)*sin(2*theta)
%                                 + Ac4(freq)*cos(4*theta) + As4(freq)*sin(4*theta)
%

clear; close all;

comp = {'ZZ'};
xspdir = 'test_1.6win_avg'; %'Nomelt3inttaper_iso.s0to333_br1avg'; %'4.0_S1_10pers_avg'; %'Nomelt3inttaper_iso.s0to333_br1avg'; %'4.0_S0_waverage';

% comp = {'TT'};
% xspdir = 'Nomelt3inttaper_iso.t0to500_br0avg'; %'Nomelt3inttaper_iso_waverage'; %'4.0_S0_waverage'; %'Nomelt3inttaper_iso_T1'; %'Nomelt3inttaper_iso_waverage'; %'Nomelt3inttaper_iso'; %'4.2kmsstart';

stafile = 'stations_all'; %'stations_shallow'; %'stations_deep'; %'stations_all'

% aniso_data = 'phv_2theta4theta_wRMS_SNRtol0_disttol200_errtol0.7.mat';
frange = [1/40 1/15]; %[1/10 1/4]; %[1/30 1/12]; %[0.1 0.25];

average_vel = 3.8; % (km/s) For calculating wavelength for determining r_tol_min

% QC parameters
snr_tol = 3; %0; % 5
is_rtolmin_wavelength = 0; wl_fac = 1.0; % determine distance tolerance by wavelength?
r_tol_min = 90; %90; %150; %200; %200; % 100 km
r_tol_max = 600;
% err_tol = 0.5; % FOR AGU17
err_tol = 0.3; %100; %0.7; %100;

% % 1D
% snr_tol = 3; %0; % 5
% r_tol_min = 150; %90; %150; %200; %200; % 100 km
% r_tol_max = 600;
% % err_tol = 0.5; % FOR AGU17
% err_tol = 0.5; %100; %0.7; %100;

% Norm damping for azimuthal anisotropy
% damp_azi = [1 1 1e10 1e10]; % [2c 2s 4c 4s] % Damping individual parameters
% aziweight = 1; % global weight

fastdir = 95; % Fast direction for azimuthal anisotropy (only for plotting purposes);

iscompare_aniso = 0; % compare to old anisotropic measurements

windir = 'window3hr_LH_Zcorr_tiltonly'; %'window3hr_LH_Zcorr'; %'window3hr_LH_Zcorr'; %'window0.2hr'; %'window24hr_specwhite';

% Save results?
isoutput = 1;

%%

% Load color scale
load seiscmap.mat

% Load station info
[sta.name, sta.lat, sta.lon, sta.dep] = textread(stafile,'%s %f %f %f');

% Load anisotropy data (from old inversion)
if iscompare_aniso
    load(['./aniso_DATA/',xspdir,'/',aniso_data]);
end

% Set up geometry parameters
setup_parameters_tomo;
setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
gridsize_azi = parameters.gridsize_azi;

% Set up error parameters
% errlevel = parameters.errlevel;
% snrtol = parameters.snrtol;
% mincoherenum = parameters.mincoherenum;
fiterrtol = parameters.fiterrtol;

% refv = parameters.refv;
% distrange = parameters.distrange;

maxerrweight = parameters.maxerrweight;
polyfit_dt_err = parameters.polyfit_dt_err;
smweight0 = parameters.smweight0;
smweight0_azi = parameters.smweight0_azi;
flweight0_azi = parameters.flweight0_azi;


dterrtol = parameters.dterrtol;

raydensetol = parameters.raydensetol;
raydensetol_azi = parameters.raydensetol_azi;
r = parameters.r;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
Nx = length(xnode);
Ny = length(ynode);

xnode_azi=lalim(1):gridsize_azi:lalim(2);
ynode_azi=lolim(1):gridsize_azi:lolim(2);
[xi_azi yi_azi] = ndgrid(xnode_azi,ynode_azi);
Nx_azi = length(xnode_azi);
Ny_azi = length(ynode_azi);

% savefile = parameters.savefile;
savefile = ['test'];

% figure output path
phv_fig_path = ['./figs/',windir,'/fullStack/raytomo_azi2theta_2D/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/',stafile,'/'];
if ~exist(phv_fig_path)    
    mkdir(phv_fig_path);
end

% read in bad station list, if existed
if exist('badsta.lst')
    badstnms = textread('badsta.lst','%s');
    badstaids = find(ismember({stainfo.staname},badstnms));
    disp('Found Bad stations:')
    disp(badstnms)
end

%% Set up initial smoothing kernels (second derivative)
% Isotropic smoothing kernels
F_iso = smooth_kernel_build(xnode, ynode, Nx*Ny);
F = sparse(Nx*Ny*2,Nx*Ny+Nx_azi*Ny_azi*2);
F(1:end,1:Nx*Ny) = F_iso;
% Azimuthal smoothing kernels
F_azi = smooth_kernel_build(xnode_azi, ynode_azi, Nx_azi*Ny_azi);
F_azic = sparse(Nx_azi*Ny_azi*2,Nx*Ny+Nx_azi*Ny_azi*2);
F_azic(1:end,Nx*Ny+[1:Nx_azi*Ny_azi]) = F_azi;
F_azis = sparse(Nx_azi*Ny_azi*2,Nx*Ny+Nx_azi*Ny_azi*2);
F_azis(1:end,Nx*Ny+Nx_azi*Ny_azi+[1:Nx_azi*Ny_azi]) = F_azi;

%% Set up initial flattening kernels (first derivative)
% Azimuthal flattening kernels
J_azi = flat_kernel_build(xnode_azi, ynode_azi, Nx_azi*Ny_azi);
J_azic = sparse(Nx_azi*Ny_azi*2,Nx*Ny+Nx_azi*Ny_azi*2);
J_azic(1:end,Nx*Ny+[1:Nx_azi*Ny_azi]) = J_azi;
J_azis = sparse(Nx_azi*Ny_azi*2,Nx*Ny+Nx_azi*Ny_azi*2);
J_azis(1:end,Nx*Ny+Nx_azi*Ny_azi+[1:Nx_azi*Ny_azi]) = J_azi;

%%
% Initialize the xsp structure
% Xsp_path = './Xsp/';
Xsp_path = ['../Xsp/',windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/'];
xspfiles = dir([Xsp_path,'*_xsp.mat']);

disp('Looking at Xsp Files')
for ixsp = 1:length(xspfiles)
    
    temp = load([Xsp_path,xspfiles(ixsp).name]);
    xspinfo = temp.xspinfo;
    
    if ixsp ==1
        Tperiods = (2*pi)./temp.twloc;
        waxis = temp.waxis;
        twloc = temp.twloc;
        xspinfo.isgood = zeros(size(Tperiods));
        xspsum = xspinfo;
        wavelength = 3.8*Tperiods;
    else
        xspinfo.isgood = zeros(size(Tperiods));
        xspsum = [xspsum;xspinfo];
    end
    clear temp

    
    % 	xspinfo(ixsp).isgood = 0;
%     if xspsum(ixsp).sumerr < errlevel ...
%             && xspsum(ixsp).snr > snrtol && xspsum(ixsp).coherenum > mincoherenum
%         xspsum(ixsp).isgood = 1;
%     end

    for ip = 1:length(Tperiods)
        if ~is_rtolmin_wavelength && xspinfo.snr >= snr_tol && xspinfo.r >= r_tol_min && xspinfo.r <= r_tol_max && xspinfo.sumerr <= err_tol
            xspsum(ixsp).isgood(ip) = 1;
        elseif  is_rtolmin_wavelength && xspinfo.snr >= snr_tol && xspinfo.r >= wavelength(ip)*wl_fac && xspinfo.r <= r_tol_max && xspinfo.sumerr <= err_tol
            xspsum(ixsp).isgood(ip) = 1;
        end

        if isempty(cell2mat(strfind(sta.name,xspsum(ixsp).sta1))) || isempty(cell2mat(strfind(sta.name,xspsum(ixsp).sta2)))
            xspsum(ixsp).isgood(ip) = 0;
        end
    end
    
    if rem(ixsp,500)==0
        disp(['Looking at #',num2str(ixsp),' of ',num2str(length(xspfiles))])
    end
end % end of loop ixsp'


% Loop through periods
for ip=1:length(Tperiods)
    disp(' ');
    disp(['Inversing Period: ',num2str(Tperiods(ip))]);
    clear rays dt fiterr mat phaseg err raydense dist azi mat_azi phv
    raynum = 0;

    for ixsp = 1:length(xspsum)
        if xspsum(ixsp).isgood(ip) ==0
            continue;
        end
%         if xspsum(ixsp).r > refv*Tperiods(ip)*distrange(2)...
%                 || xspsum(ixsp).r < refv*Tperiods(ip)*distrange(1)
%             continue;
%         end
        
        raynum = raynum+1;
        rays(raynum,1) = xspsum(ixsp).lat1;
        rays(raynum,2) = xspsum(ixsp).lon1;
        rays(raynum,3) = xspsum(ixsp).lat2;
        rays(raynum,4) = xspsum(ixsp).lon2;
        
        dist(raynum) = deg2km(distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4)));
        dt(raynum) = xspsum(ixsp).tw(ip);
        phv(raynum) = dist(raynum)./dt(raynum);
        
        dep1 = sta.dep(strcmp(xspsum(raynum).sta1,sta.name));
        dep2 = sta.dep(strcmp(xspsum(raynum).sta2,sta.name));
        dep(raynum) = mean([dep1 dep2]);
        
        %JRB - load azimuthal anisotropy
        [~,azi(raynum)]=distance(xspsum(ixsp).lat1,xspsum(ixsp).lon1,xspsum(ixsp).lat2,xspsum(ixsp).lon2);
        if azi(raynum) > 180
            azi(raynum) = azi(raynum) - 360;
        end
        if iscompare_aniso
            A2 = aniso.A2(ip);
            A4 = aniso.A4(ip);
            phi2 = aniso.phi2(ip);
            phi4 = aniso.phi4(ip);
            c_iso = aniso.c_iso(ip);       
        end
%         if comp{1}(1) == 'Z'
%             phv_cor = dist(raynum)./dt(raynum) - A2*c_iso*cosd(2*(azi - phi2));
%         elseif comp{1}(1) == 'T'
%             phv_cor = dist(raynum)./dt(raynum) - A2*c_iso*cosd(2*(azi-phi2)) - A4*c_iso*cosd(4*(azi-phi4)); 
%         end
%         dt(raynum) = dist(raynum)./phv_cor;
        
        err = smooth((abs(xspsum(ixsp).err)./mean(abs(xspsum(ixsp).xsp))).^2,round(length(waxis)/length(twloc)));
        fiterr(raynum) = interp1(waxis(:),err(:),twloc(ip)); 
        % Fix the fact that last period always breaks (JBR 9/29/17)
        if isnan(fiterr(raynum))
            [~,I] = min(abs(twloc(ip)-waxis(:)));
            fiterr(raynum) = err(I);
        end
        csnum(raynum) = xspsum(ixsp).coherenum;
        snr(raynum) = xspsum(ixsp).snr;
        errays(raynum,1) = xspsum(ixsp).lat1;
        errays(raynum,2) = xspsum(ixsp).lon1;
        errays(raynum,3) = xspsum(ixsp).lat2;
        errays(raynum,4) = xspsum(ixsp).lon2; 
        errays(raynum,5) = fiterr(raynum);
        
        % JBR - Build azimuthal part of data kernel
%         mat_azi(raynum,:) = dist(raynum) * [cosd(2*azi(raynum)), sind(2*azi(raynum)), cosd(4*azi(raynum)), sind(4*azi(raynum)) ];
   
    end
    if size(dt,1) ~=raynum
        dt = dt';
    end
    
    % Building the isotropic data kernel
    disp('Start building the kernel');
    tic
    mat_iso=ray_kernel_build(rays,xnode,ynode);   
    toc
    [mat_azi, mat_azi_hits] = ray_kernel_build_azi(rays,xnode_azi,ynode_azi);
    
    % JBR - Combine isotropic and anisotropic
    mat = [mat_iso, mat_azi];
    
    % Calculate the weighting matrix
    W = sparse(length(dt),length(dt));
    for i=1:length(dt)
        W(i,i)=1./fiterr(i);
    end
    ind = find(W > maxerrweight);
    W(ind) = maxerrweight;
    ind = find(W < 1/fiterrtol);
    W(ind) = 0;
    for i=1:length(dt)
        W(i,i)=W(i,i).*(csnum(i).^0.5);
    end
    para = polyfit(dist(:),dt,1);
    polyerr = polyval(para,dist(:)) - dt;
    errind = find(abs(polyerr) > polyfit_dt_err);
    for i = errind
        W(i,i) = 0;
    end
    
    % calculate the smoothing weight
    NR=norm(F,1);
    NA=norm(W*mat,1);
    smweight = smweight0*NA/NR;
    
    NR=norm(F_azic,1);
    NA=norm(W*mat,1);
    smweight_azi = smweight0_azi*NA/NR;
    
    NR=norm(J_azic,1);
    NA=norm(W*mat,1);
    flweight_azi = flweight0_azi*NA/NR;
    
    disp('start inverse');
%     A=[W*mat; smweight*F; aziweight*F_azi_damp; aziweight*F_azi_smooth];
%     rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1); zeros(size(F_azi_smooth,1),1)];
%     A=[W*mat;smweight*F;aziweight*F_azi_damp];
%     rhs=[W*dt;zeros(size(F,1),1);zeros(size(F_azi_damp,1),1)];
% %     A=[W*mat; smweight*F; aziweight*F_azi_smooth];
% %     rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_smooth,1),1)];
%     A=[W*mat; smweight*F; aziweight*F_azi_damp];
%     rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1)];
%     A=[W*mat; smweight*F; smweight_azi*F_azic; smweight_azi*F_azis];
%     rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1)];
    A=[W*mat; smweight*F; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];
    rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1); zeros(size(J_azic,1),1); zeros(size(J_azis,1),1)];

    phaseg=(A'*A)\(A'*rhs);
    %        toc
    %        disp('Done');
    
    
    % Iteratively down weight the measurement with high error
    niter=1;
    
    while niter < 2
        niter=niter+1;
        err = mat*phaseg - dt;

        stderr=std(err);
        if stderr > dterrtol
            stderr = dterrtol;
        end
        ind = find(diag(W)==0);
        disp('Before iter:');
        disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
        disp(['Bad Measurement Number: ', num2str(length(ind))]);
        for i=1:length(err)
            if abs(err(i)) > 2*stderr
                W(i,i)=0;
            end
        end
        ind = find(diag(W)==0);
        disp('After iter:');
        disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
        disp(['Bad Measurement Number: ', num2str(length(ind))]);
        
        % Rescale the smooth kernel
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0*NA/NR;
        
        NR=norm(F_azic,1);
        NA=norm(W*mat,1);
        smweight_azi = smweight0_azi*NA/NR;
        
        NR=norm(J_azic,1);
        NA=norm(W*mat,1);
        flweight_azi = flweight0_azi*NA/NR;
        
        % Invert
%         A=[W*mat;smweight*F];
%         rhs=[W*dt;zeros(size(F,1),1)];        
%         A=[W*mat;smweight*F;aziweight*F_azi];
%         rhs=[W*dt;zeros(size(F,1),1);zeros(size(F_azi,1),1)];
%         A=[W*mat; smweight*F; aziweight*F_azi_damp; aziweight*F_azi_smooth];
%         rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1); zeros(size(F_azi_smooth,1),1)];
%         A=[W*mat;smweight*F;aziweight*F_azi_damp];
%         rhs=[W*dt;zeros(size(F,1),1);zeros(size(F_azi_damp,1),1)];
% %         A=[W*mat; smweight*F; aziweight*F_azi_smooth];
% %         rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_smooth,1),1)];
%         A=[W*mat; smweight*F; aziweight*F_azi_damp];
%         rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1)];
%         A=[W*mat; smweight*F; smweight_azi*F_azic; smweight_azi*F_azis];
%         rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1)];
        A=[W*mat; smweight*F; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];
        rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1); zeros(size(J_azic,1),1); zeros(size(J_azis,1),1)];

        phaseg=(A'*A)\(A'*rhs);

        
    end
    
    % Anisotropic terms from model vector
    phaseg_azic = phaseg(Nx*Ny+1 : Nx*Ny+Nx_azi*Ny_azi);
    phaseg_azis = phaseg(Nx*Ny+Nx_azi*Ny_azi+1 : Nx*Ny+Nx_azi*Ny_azi*2);
    
    % Isotropic phase velocity
    phv_iso = dist'./(mat_iso*phaseg(1:Nx*Ny));
    
    %        disp(' Get rid of uncertainty area');
    
    Igood = find(diag(W)~=0);
    mat_good = mat(Igood,:);
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            %raydense(i,j) = sum(mat(:,n));
            raydense(i,j) = sum(mat_good(:,n));
            if raydense(i,j) < raydensetol
                phaseg(n)=NaN;
            end
        end
    end
    for i=1:Nx_azi
        for j=1:Ny_azi 
            n=Ny_azi*(i-1)+j;
            raydense_azi(i,j) = sum(mat_azi_hits(:,n));
            if raydense_azi(i,j) < raydensetol_azi
                phaseg_azic(n)=NaN;
                phaseg_azis(n)=NaN;
            end
        end
    end
    
    % Convert into phase velocity
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            GV(i,j)= 1./phaseg(n);
        end
    end
    for i=1:Nx_azi
        for j=1:Ny_azi
            n=Ny_azi*(i-1)+j;
            Ac(i,j) = -phaseg_azic(n);
            As(i,j) = -phaseg_azis(n);
        end
    end
    
    % JBR - Get Azimuthal coefficients from phaseg (s/km)
    phv_av = nanmean(GV(:));
%     phv_av = nanmedian(GV(:));
    phv_avstd = nanstd(GV(:));
    slow_av = 1/phv_av;
%     slow_av = 1/mean(phv_iso);
    Ac2 = Ac./slow_av;
    As2 = As./slow_av;
    A2 = sqrt(Ac2.^2+As2.^2);
    phi2 = 1/2*atan2d(As2,Ac2);

    raytomo(ip).GV = GV;
    raytomo(ip).mat = mat;
    raytomo(ip).raydense = raydense;
    raytomo(ip).period = Tperiods(ip);
    raytomo(ip).w = diag(W);
    raytomo(ip).err = err;
    raytomo(ip).rays = rays;
    raytomo(ip).fiterr = fiterr;
    raytomo(ip).dt = dt;
    raytomo(ip).smweight0 = smweight0;
    raytomo(ip).smweight0_azi = smweight0_azi;
    raytomo(ip).flweight0_azi = flweight0_azi;
    %JBR    
    raytomo(ip).phv_iso = phv_iso;    
    raytomo(ip).phv_av = phv_av;
    raytomo(ip).phv_avstd = phv_avstd;
    raytomo(ip).Ac2 = Ac2;
    raytomo(ip).As2 = As2;
%     raytomo(ip).Ac4 = Ac4;
%     raytomo(ip).As4 = As4;
    raytomo(ip).A2 = A2;
%     raytomo(ip).A4 = A4;
    raytomo(ip).phi2 = phi2;
%     raytomo(ip).phi4 = phi4;
    raytomo(ip).phv = phv;
    raytomo(ip).azi = azi;
    
    
    if 0
        figure(1)
        clf
        ax = worldmap(lalim, lolim);
        set(ax, 'Visible', 'off')
        surfacem(xi,yi,raytomo(ip).GV);
%         drawlocal
        title([num2str(Tperiods(ip))],'fontsize',15)
        avgv = nanmean(raytomo(ip).GV(:));
        caxis([avgv*(1-r) avgv*(1+r)])
        colorbar
        colormap(seiscmap)
        
%         pause;
    end
    
end % end of period loop

lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
% isoutput = 1;
if isoutput
    save(savefile,'raytomo','xnode','ynode');
    save('coor.mat','xi','yi','xnode','ynode','gridsize','lalim','lolim');
end

%% Azimuthal anisotropy (%)

Mp = 3; Np = 4;
fig16 = figure(16);
set(gcf,'position',[1    1   1244   704]);
clf
vperc = [-r r];
for ip=1:length(Tperiods)
    subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    set(gcf,'color',[0.9 0.9 0.9])
    A2 = raytomo(ip).A2;
%     surfacem(xi,yi,resid);
    levels = linspace(0,0.03,10)*100;
    surfacem(xi_azi,yi_azi,A2*100,'Linestyle','none');
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    caxis([min(levels) max(levels)])
    colorbar
%     colormap(seiscmap)
    rbc = flip(redbluecmap);
%     rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
%     colormap(rbc);
    colormap('parula');
    
    u=raytomo(ip).A2 .* cosd(raytomo(ip).phi2)*20;
	v=raytomo(ip).A2 .* sind(raytomo(ip).phi2)*20./cosd(mean(lalim));
	[m n]=size(xi_azi);
    hold on;
    xpts = [];
    ypts = [];
	for ix=1:m
		for iy=1:n
% 			if avgphv_aniso(ip).aniso_azi_std(ix,iy) < 40 && avgphv_aniso(ip).aniso_strength(ix,iy)>0.02
%                 [xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2]
%                 [yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2]
% % 			geoshow([xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2]+gridsize_azi/2,...
% % 					[yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2]+gridsize_azi/2,'Color','k','linewidth',2);
%             plotm([yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2],...
% 					[xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2],'k-','linewidth',2);
% 			end
            xpts = [xpts, [xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2]+gridsize_azi/2, nan];
            ypts = [ypts, [yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2]+gridsize_azi/2, nan];
        end
	end
    plotm(xpts,ypts,'Color','k','linewidth',2);
    hold on;
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
end
% save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_raytomo_2Dazimuthal.pdf'],fig16,1000);
% stop
%% Phase Velocity Maps (km/s)
% % Load seafloor age
% load('age_grid.mat');
% 
% Mp = 3; Np = 4;
% fig17 = figure(17);
% set(gcf,'position',[1    1   1244   704]);
% clf
% for ip=1:length(Tperiods)
%     subplot(Mp,Np,ip)
%     ax = worldmap(lalim, lolim);
%     set(ax, 'Visible', 'off')
%     set(gcf,'color',[0.9 0.9 0.9])
% %     surfacem(xi,yi,raytomo(ip).GV);
%     avgv = nanmean(raytomo(ip).GV(:));
%     levels = linspace(avgv*(1-r), avgv*(1+r),10);
%     contourfm(xi,yi,raytomo(ip).GV,levels);
% %     drawlocal
%     title([num2str(Tperiods(ip))],'fontsize',15)
%     caxis([avgv*(1-r) avgv*(1+r)])
%     colorbar
% %     colormap(seiscmap)
%     rbc = flip(redbluecmap);
% %     rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
%     colormap(rbc);
%     
%     hold on;
%     plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
%     [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
% end
% save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_raytomo.pdf'],fig17,1000);

%% Phase Velocity Maps (%)
% Load seafloor age
load('age_grid.mat');

Mp = 3; Np = 4;
fig19 = figure(19);
set(gcf,'position',[1    1   1244   704]);
clf
vperc = [-r r];
for ip=1:length(Tperiods)
    subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    set(gcf,'color',[0.9 0.9 0.9])
    avgv = nanmean(raytomo(ip).GV(:));
    resid = (raytomo(ip).GV-avgv)./avgv;
    levels = linspace(vperc(1),vperc(2),10)*100;
    contourfm(xi,yi,resid*100,levels);
    title([num2str(Tperiods(ip))],'fontsize',15)
    caxis([min(levels) max(levels)])
    colorbar
%     colormap(seiscmap)
    rbc = flip(redbluecmap);
%     rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
    colormap(rbc);
    
    hold on;
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
    [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
end
% save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_raytomo_perc.pdf'],fig19,1000);

% stop
%%
% RAY DENSITY
fig18 = figure(18);
set(gcf,'position',[1    1   1244   704]);
clf

for ip=1:length(Tperiods)
subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).raydense);
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    colorbar
    colormap(flip(hot));
    caxis([0 500])
end
save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_raydense.pdf'],fig18,1000);

% % ERRORS ON XSP
% fig19 = figure(19)
% set(gcf,'position',[1    1   1244   704]);
% clf
% 
% for ip=1:length(Tperiods)
%     subplot(Mp,Np,ip)
%     ax = worldmap(lalim, lolim);
%     set(ax, 'Visible', 'off')
%     clear rays
%     rays = raytomo(ip).rays;
% %     surfacem(xi,yi,raytomo(ip).err);
% %     drawpng
% scatterm((rays(:,1)+rays(:,3))./2,(rays(:,2)+rays(:,4))./2,30,raytomo(ip).fiterr,'filled')
% % drawlocal
% title([num2str(Tperiods(ip))],'fontsize',15)
% colorbar
% caxis([ 0 5])
% 
% 
% end

%% PLOT AZIMUTHAL MODEL

fig4 = figure(4); clf;
% set(gcf,'position',[4         325        1239         380],'color','w');
set(gcf,'position',[6   220   490   485]);
periods = Tperiods;
% Old Values
if iscompare_aniso
    wRMS_2A = aniso.wRMS_2A;
    wRMS_4A = aniso.wRMS_4A;
    err_phi2 = aniso.err_phi2;
    err_phi4 = aniso.err_phi4;
    A2_2 = aniso.A2;
    A4_2 = aniso.A4;
    phi2_2 = aniso.phi2;
    phi4_2 = aniso.phi4;
end
for iper = 1:length(periods)
    % New values
    A2_rt(iper) = nanmean(raytomo(iper).A2(:));
%     A4_rt(iper) = raytomo(iper).A4;
    phi2_rt(iper) = nanmean(raytomo(iper).phi2(:));
%     phi4_rt(iper) = raytomo(iper).phi4;
    phv_av_rt(iper) = raytomo(iper).phv_av;
    phv_avstd_rt(iper) = raytomo(iper).phv_avstd;
end
% A2_rt = sqrt(Ac2.^2 + As2.^2)./phv_av_rt;
% A4_rt = sqrt(Ac4.^2 + As4.^2)./phv_av_rt;
% phi2_rt = 1/2*atand(As2./Ac2);
% phi4_rt = 1/4*atand(As4./Ac4);


% peak-to-peak
subplot(2,1,1); hold on;
if iscompare_aniso
    h3(2) = errorbar(periods,A4_2*2*100,wRMS_4A*100,'--ob','linewidth',2);
    h3(1) = errorbar(periods,A2_2*2*100,wRMS_2A*100,'--o','color',[0 0.7 0],'linewidth',2);
end
% h3(2) = plot(periods,A4_rt*2*100,'-ob','linewidth',2);
h3(1) = plot(periods,A2_rt*2*100,'-o','color',[0 0.7 0],'linewidth',2);
xlim(flip(1./frange));
ylim([0 5]);
set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',18);
xlabel('Period (s)','fontsize',18);
ylabel('Peak-to-peak amp (%)','fontsize',18);
legend(h3,{'2\theta','4\theta'},'location','northwest');

% Azimuth
subplot(2,1,2); hold on;
for iper = 1:length(periods)
    % OLD AZIMUTHAL PARAMTERS
    if iscompare_aniso
        phi2_vec(1) = phi2_2(iper);
        phi2_vec(2) = phi2_2(iper)+180;
        phi2_vec(3) = phi2_2(iper)-180;
        if comp{1}(1) == 'Z'
            [dif, I] = min(abs(phi2_vec-fastdir));
            phi2_2(iper) = phi2_vec(I);
            if phi2_2(iper) < fastdir && dif > 10
                phi2_2(iper) = phi2_2(iper)+180;
            end
        elseif comp{1}(1) == 'T'
            [dif, I] = min(abs(phi2_vec-fastdir+90));
            phi2_2(iper) = phi2_vec(I);
            if phi2_2(iper) < fastdir && dif
                phi2_2(iper) = phi2_2(iper)+180;
            end
        end


        phi4_vec(1) = phi4_2(iper);
        phi4_vec(2) = phi4_2(iper)+90;
        phi4_vec(3) = phi4_2(iper)+180;
        phi4_vec(4) = phi4_2(iper)+270;
        phi4_vec(5) = phi4_2(iper)-90;
        phi4_vec(6) = phi4_2(iper)-180;
        phi4_vec(7) = phi4_2(iper)-270;
        [~, I] = min(abs(phi4_vec-fastdir+45));
        phi4_2(iper) = phi4_vec(I);
        if phi4_2(iper) < fastdir
            phi4_2(iper) = phi4_2(iper)+90;
        end
    end
    
    % NEW AZIMUTHAL PARAMTERS
    phi2_vec(1) = phi2_rt(iper);
    phi2_vec(2) = phi2_rt(iper)+180;
    phi2_vec(3) = phi2_rt(iper)-180;
    if comp{1}(1) == 'Z'
        [dif, I] = min(abs(phi2_vec-fastdir));
        phi2_rt(iper) = phi2_vec(I);
        if phi2_rt(iper) < fastdir && dif > 10
            phi2_rt(iper) = phi2_rt(iper)+180;
        end
    elseif comp{1}(1) == 'T'
        [dif, I] = min(abs(phi2_vec-fastdir+90));
        phi2_rt(iper) = phi2_vec(I);
        if phi2_rt(iper) < fastdir
            phi2_rt(iper) = phi2_rt(iper)+180;
        end
    end
    phi2_rt(phi2_rt>160) = phi2_rt(phi2_rt>160)-180;


%     phi4_vec(1) = phi4_rt(iper);
%     phi4_vec(2) = phi4_rt(iper)+90;
%     phi4_vec(3) = phi4_rt(iper)+180;
%     phi4_vec(4) = phi4_rt(iper)+270;
%     phi4_vec(5) = phi4_rt(iper)-90;
%     phi4_vec(6) = phi4_rt(iper)-180;
%     phi4_vec(7) = phi4_rt(iper)-270;
%     [~, I] = min(abs(phi4_vec-fastdir+45));
%     phi4_rt(iper) = phi4_vec(I);
%     if phi4_rt(iper) < fastdir
%         phi4_rt(iper) = phi4_rt(iper)+90;
%     end
end
plot(periods,ones(size(periods))*fastdir,'--k','linewidth',2);
plot(periods,ones(size(periods))*fastdir-90,'--','color',[0.5 0.5 0.5],'linewidth',2);
plot(periods,ones(size(periods))*fastdir+90,'--','color',[0.5 0.5 0.5],'linewidth',2);
if iscompare_aniso
    errorbar(periods,phi4_2,err_phi4,'--ob','linewidth',2);
    errorbar(periods,phi2_2,err_phi2,'--o','color',[0 0.7 0],'linewidth',2);
end
% plot(periods,phi4_rt,'-ob','linewidth',2);
plot(periods,phi2_rt,'-o','color',[0 0.7 0],'linewidth',2);
ylabel('Fast Direction (%)','fontsize',18);
ylim([fastdir-130 fastdir+130]);
xlim(flip(1./frange));
set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',18);
xlabel('Period (s)','fontsize',18);


% % PLOT COSINE AND SINE
% fig5 = figure(5);
% clf
% plot(periods,Ac2./phv_av_rt*100,'-k','linewidth',2); hold on;
% plot(periods,As2./phv_av_rt*100,'-r','linewidth',2);
% plot(periods,Ac4./phv_av_rt*100,'--k','linewidth',2);
% plot(periods,As4./phv_av_rt*100,'--r','linewidth',2);
% ylabel('A (%)','fontsize',18);
% title('Azimuthal Coefficients','fontsize',18);
% xlim(flip(1./frange));
% set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',18);
% xlabel('Period (s)','fontsize',18);
% legend({'A_{c2}','A_{s2}','A_{c4}','A_{s4}'},'fontsize',13,'box','off');

% save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_A_phi_plots.pdf'],fig4,1000);
stop
%% Plot phase velocities
for ip = 1:length(Tperiods)
    avgv(ip) = nanmean(raytomo(ip).GV(:));
    avgv_std(ip) = nanstd(raytomo(ip).GV(:));
end

fig2 = figure(2);clf;
set(gcf,'position',[320     2   508   703]);

subplot(2,1,1);
hold on; box on;
if iscompare_aniso
    errorbar(Tperiods,aniso.c_iso,aniso.err_c_iso*2,'-k','linewidth',2);
end
try 
    plot(Tperiods,xspinfo.c_start,'ok','linewidth',2);
catch
    display('No starting data')
end
errorbar(Tperiods,phv_av_rt,phv_avstd_rt*2,'-r','linewidth',2);
% errorbar(Tperiods,mean([vertcat(raytomo(:).phv)],2),std([vertcat(raytomo(:).phv)],0,2)*2,'-b','linewidth',2);
title('Isotropic Phase Velocity');
xlabel('Period (s)','fontsize',16);
ylabel('Phase Velocity (km/s)','fontsize',16);
set(gca,'fontsize',16,'linewidth',1.5);
legend({'Starting','Raytomo Avg.'},'location','southeast','fontsize',12,'box','off');
xlim(flip(1./frange));
if comp{1}(1) == 'Z'
    ylim([3.4 4.3]);
elseif comp{1}(1) == 'T'
    ylim([3.8 4.7]);
end

if iscompare_aniso
    subplot(2,1,2);
    box on; hold on;
    resids = (phv_av_rt-aniso.c_iso)./aniso.c_iso*100;
    resids2 = (mean([vertcat(raytomo(:).phv)],2)'-aniso.c_iso)./aniso.c_iso*100;
    patch([aniso.periods fliplr(aniso.periods)],[aniso.err_c_iso./aniso.c_iso*2*100 fliplr(-aniso.err_c_iso./aniso.c_iso*2*100)],[0.8 0.8 0.8],'linestyle','none');
    errorbar(Tperiods,resids,phv_avstd_rt./phv_av_rt*100*2,'-r','linewidth',2);
    % errorbar(Tperiods,resids2,std([vertcat(raytomo(:).phv)],0,2)./mean([vertcat(raytomo(:).phv)],2)*2*100,'-b','linewidth',2)
    title('Residual (%)');
    xlabel('Period (s)','fontsize',16);
    ylabel('\delta c/c','fontsize',16);
    set(gca,'fontsize',16,'linewidth',1.5);
    xlim([3.5 10.5]);
    ylim([-2 2]);
end

save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_compareisophv.pdf'],fig2,1000);

%% Plot Azimuthal Data
fig6 = figure(6); clf;
% set(gcf,'position',[10         248        1203         457]);
set(gcf,'position',[10          11        1203         695]);
for iper = 1:length(periods)
    azi = raytomo(iper).azi;
%     dphv = (raytomo(iper).phv - phv_av_rt(iper))./phv_av_rt(iper);
     dphv = (raytomo(iper).phv' - raytomo(iper).phv_iso) ./ raytomo(iper).phv_iso;
    subplot(3,4,iper); hold on;
%     if comp{1}(1) == 'Z' || comp{1}(1) == 'R'
%         c = 2; % 2 theta
%         e_patty = 78;
%     elseif comp{1}(1) == 'T'
%         c = 4; % 4 theta
%         e_patty = 78-45;
%     end
    x = [-180:180];
    % PHV FIT = a*(1+d*cosd(c*(x-e)))
    if iscompare_aniso
        h2(1) = plot(x,A2_2(iper)*cosd(2*(x-phi2_2(iper)))*100+A4_2(iper)*cosd(4*(x-phi4_2(iper)))*100,'--','color',[0.5 0.5 0.5],'linewidth',3);
    end
    h2(2) = plot(x,A2_rt(iper)*cosd(2*(x-phi2_rt(iper)))*100+A4_rt(iper)*cosd(4*(x-phi4_rt(iper)))*100,'-','color',[0.5 0.5 0.5],'linewidth',3);
%     h2(1) = plot(x,d4*cosd(4*(x-e4))*100,'-b','linewidth',3);
%     h2(2) = plot(x,d2*cosd(2*(x-e2))*100,'-','color',[0 0.7 0],'linewidth',3);
%     plot(azi,dphv*100,'xr','linewidth',1); hold on;
%     scatter(azi,dphv*100,30,dep,'filled'); hold on;
    scatter(azi,dphv*100,30,dist,'filled'); hold on;
    
    if iper == 4
        ax = get(gca);
        pos = ax.Position;
        colorbar;
        set(gca,'Position',pos);
    end
    if iper == 1
%         legend(h2,{'2\theta + 4\theta'},'location','northwest');
    end
    title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
    if iper == length(periods)
        xlabel('Azimuth (degrees)','fontsize',15);
    end
    if iper == 1
        ylabel('\delta{c}/c (%)','fontsize',15);
    end
    set(gca,'fontsize',18);
    xlim([-180 180]);
    ylim([-10 10]);
    %ylim([3.8 4.8]);
    %box on;
    
end

save2pdf([phv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_sinplots_24theta.pdf'],fig6,1000);
