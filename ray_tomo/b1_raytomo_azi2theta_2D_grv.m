% Script to do ray theory tomography using ambient-noise group velocity measurements. 
% Solves for azimuthal anisotropy and isotropic velocity on a 2D grid.
%
% grv(theta,freq) = grv_iso(freq) + Ac2(freq)*cos(2*theta) + As2(freq)*sin(2*theta)
%                                 + Ac4(freq)*cos(4*theta) + As4(freq)*sin(4*theta)
%
% Written by Ge Jin, jinwar@gmail.com
% Nov 2012
% JBR: Modified in 2018 to include azimuthal anisotropy
% https://github.com/jbrussell
% JBR: 2025 - Modified for group velocity!
clear; close all;

%%
%======================= PARAMETERS =======================%
% Save results?
isoutput = 1;
savefile = ['test'];
fastdir = 78; % Fast direction for azimuthal anisotropy (only for plotting purposes);
iscompare_aniso = 0; % compare to old anisotropic measurements
Mp = 4; % # of subplot rows
Np = 4; % # of subplot columns

setup_parameters_tomo;
comp = parameters.grv.comp; % {'ZZ'};
% ftan_grv_path = parameters.grv.ftan_grv_path; % FTAN grv project
windir = parameters.grv.windir; %'window3hr'; 
% N_wl = parameters.N_wl;
frange = parameters.grv.frange; %[1/10 1/5]; % [Hz]
per_ind = parameters.grv.per_ind; % [1:12]; % index of periods to consider
nBranches = parameters.grv.nBranches; % number of mode branches to consider
branch2analyze = parameters.grv.branch2analyze; % which branch do you want to analyze?

% QC parameters
snr_tol = parameters.grv.snr_tol; % 3; % minimum signal-to-noise
corr_tol = parameters.grv.corr_tol; % Minimum coherence required between grv estimates from -lag, +lag, stacked
is_rtolmin_wavelength = parameters.grv.is_rtolmin_wavelength; % 0; 
wl_fac = parameters.grv.wl_fac; % 1.0; % determine distance tolerance by wavelength?
r_tol_min = parameters.grv.r_tol_min; %90; % [km] minimum station separation
r_tol_max = parameters.grv.r_tol_max; %600; % [km] maximum station separation
% err_tol = parameters.err_tol; %0.5; % maximum misfit of bessel fit between observed and synthetic
is_raydensity_thresh = parameters.grv.is_raydensity_thresh; % Apply raydensity threshold to wipe out poorly constrained grid cells?
min_dep = parameters.grv.min_dep; %= 9999; %-3500 for min station depth to use

%==========================================================%
%%
%%
% Set up geometry parameters
setup_parameters;
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
smweight0 = parameters.smweight0;
smweight0_azi = parameters.smweight0_azi;
flweight0_azi = parameters.flweight0_azi;
damp0_azi = parameters.damp0_azi;
is_wl_smooth = parameters.is_wl_smooth; % wavelength dependent smoothing
dterrtol = parameters.dterrtol;
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

xnode_azi=lalim(1):gridsize_azi:lalim(2);
ynode_azi=lolim(1):gridsize_azi:lolim(2);
[xi_azi yi_azi] = ndgrid(xnode_azi,ynode_azi);
Nx_azi = length(xnode_azi);
Ny_azi = length(ynode_azi);

% figure output path
grv_fig_path = ['./figs/',windir,'/fullStack/raytomo_azi2theta_2D/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s/'];
if ~exist(grv_fig_path)    
    mkdir(grv_fig_path);
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

%% Anisotropy norm damping
Nxyazi = Nx_azi*Ny_azi*2;
Nxy = Nx*Ny;
Areg_azi = zeros(Nxyazi,Nxy+Nxyazi);
azipart = eye(Nxyazi,Nxyazi); %.* diag(damp_azi);
Areg_azi(1:Nxyazi,Nxy+1:Nxy+Nxyazi) = azipart;
F_azi_damp = Areg_azi;

%%
% Initialize the ftan structure
ftan_path = ['../ftan/',windir,'/fullStack/ftan',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_br',num2str(nBranches),'/'];
ftanfiles = dir([ftan_path,'*_ftan.mat']);

disp('Looking at FTAN Files')
ftansum = [];
for iftan = 1:length(ftanfiles)
    
    temp = load([ftan_path,ftanfiles(iftan).name]);
    ftan = temp.ftan;
    
    if iftan ==1
        Tperiods = ftan.periods;
        wavelength = ftan.stack.grv(:,branch2analyze) .* ftan.periods(:);
        ftan.isgood = zeros(size(Tperiods));
        ftansum = ftan;
    else
        wavelength = ftan.stack.grv(:,branch2analyze) .* ftan.periods(:);
        ftan.isgood = zeros(size(Tperiods));
        ftansum = [ftansum;ftan];
    end
    
    clear temp
    dep1 = sta.dep(strcmp(ftansum(iftan).sta1,sta.name));
    dep2 = sta.dep(strcmp(ftansum(iftan).sta2,sta.name));

    
    % 	ftan(iftan).isgood = 0;
%     if ftansum(iftan).sumerr < errlevel ...
%             && ftansum(iftan).snr > snrtol && ftansum(iftan).coherenum > mincoherenum
%         ftansum(iftan).isgood = 1;
%     end
    for ip = 1:length(ftan.periods)
        if ~is_rtolmin_wavelength && ftan.snr(ip) >= snr_tol && ftan.r >= r_tol_min && ftan.r <= r_tol_max && ftan.corr_avg(branch2analyze)>corr_tol %&& ftan.sumerr <= err_tol
            ftansum(iftan).isgood(ip) = 1;
        elseif  is_rtolmin_wavelength && ftan.snr(ip) >= snr_tol && ftan.r >= wavelength(ip)*wl_fac && ftan.r <= r_tol_max && ftan.corr_avg(branch2analyze)>corr_tol %&& ftan.sumerr <= err_tol
            ftansum(iftan).isgood(ip) = 1;
        end

        if isempty(cell2mat(strfind(sta.name,ftansum(iftan).sta1))) || isempty(cell2mat(strfind(sta.name,ftansum(iftan).sta2))) || dep1>min_dep || dep2>min_dep
            ftansum(iftan).isgood(ip) = 0;
        end
    end
    % Remove all nan values and make logical
    ftan.isgood(isnan(ftan.stack.grv(:,branch2analyze))) = 0;
    ftansum(iftan).isgood = logical(ftansum(iftan).isgood); %.* ftansum(iftan).isgood_wl);
    
    if rem(iftan,500)==0
        disp(['Looking at #',num2str(iftan),' of ',num2str(length(ftanfiles))])
    end
end % end of loop iftan'


% Loop through periods
for ip=1:length(Tperiods)
    disp(' ');
    disp(['Inverting Period: ',num2str(Tperiods(ip))]);
    clear rays dt fiterr mat groupg err raydense dist azi mat_azi grv grv_std dt_std
    raynum = 0;
    
    for iftan = 1:length(ftansum)
        if ftansum(iftan).isgood(ip) ==0
            continue;
        end
%         if ftansum(iftan).r > refv*Tperiods(ip)*distrange(2)...
%                 || ftansum(iftan).r < refv*Tperiods(ip)*distrange(1)
%             continue;
%         end
        
        raynum = raynum+1;
        rays(raynum,1) = ftansum(iftan).stapairsinfo.lats(1);
        rays(raynum,2) = ftansum(iftan).stapairsinfo.lons(1);
        rays(raynum,3) = ftansum(iftan).stapairsinfo.lats(2);
        rays(raynum,4) = ftansum(iftan).stapairsinfo.lons(2);
        
        % dist(raynum) = deg2km(distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4)));
        dist(raynum) = distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4),referenceEllipsoid('GRS80'))/1000;
        dt(raynum) = ftansum(iftan).stack.tg(ip);
        grv(raynum) = dist(raynum)./dt(raynum);
        
        % convert uncertainty in velocity to uncertainty in time
        % dt = |r / v^2 * dv| = t^2 / r * dv
        grv_std(raynum,1) = 0.01; %ftansum(iftan).grv_std(ip);
        dt_std(raynum,1) = abs( dt(raynum).^2 / dist(raynum) * grv_std(raynum) );
        
        dep1 = sta.dep(strcmp(ftansum(iftan).sta1,sta.name));
        dep2 = sta.dep(strcmp(ftansum(iftan).sta2,sta.name));
        dep(raynum) = mean([dep1 dep2]);
        
        %JRB - load azimuthal anisotropy
        Nr = 100;
        [lat_way,lon_way] = gcwaypts(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4),Nr);
        [dr_ray,azi_ray] = distance(lat_way(1:end-1),lon_way(1:end-1),...
                             lat_way(2:end),lon_way(2:end),referenceEllipsoid('GRS80'));
        rayazi_mean = angmean(azi_ray(:)*pi/180)*180/pi;
        azi(raynum) = rayazi_mean;
        if azi(raynum) > 180
            azi(raynum) = azi(raynum) - 360;
        end
%         if comp{1}(1) == 'Z'
%             grv_cor = dist(raynum)./dt(raynum) - A2*c_iso*cosd(2*(azi - phi2));
%         elseif comp{1}(1) == 'T'
%             grv_cor = dist(raynum)./dt(raynum) - A2*c_iso*cosd(2*(azi-phi2)) - A4*c_iso*cosd(4*(azi-phi4)); 
%         end
%         dt(raynum) = dist(raynum)./grv_cor;
        
%         csnum(raynum) = ftansum(iftan).corr_avg(branch2analyze);
        snr(raynum) = ftansum(iftan).snr(ip);
%         errays(raynum,1) = ftansum(iftan).lat1;
%         errays(raynum,2) = ftansum(iftan).lon1;
%         errays(raynum,3) = ftansum(iftan).lat2;
%         errays(raynum,4) = ftansum(iftan).lon2; 
%         errays(raynum,5) = fiterr(raynum);
        
        % JBR - Build azimuthal part of data kernel
%         mat_azi(raynum,:) = dist(raynum) * [cosd(2*azi(raynum)), sind(2*azi(raynum)), cosd(4*azi(raynum)), sind(4*azi(raynum)) ];
   
    end
    if raynum == 0
        disp(['No good data for ',num2str(Tperiods(ip))]);
        continue
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
    dt_err_frac = dt_std./dt; % error fraction of total travel time
    for i=1:length(dt)
        % W(i,i)=1./fiterr(i);
%         W(i,i)=1./dt_std(i);
        W(i,i)=1./dt_err_frac(i);
    end
    dt_err_frac_thresh = prctile(dt_err_frac,25); % 25th percentile cutoff
    ind = find(W > 1./dt_err_frac_thresh); % avoid having a few measurements dominate the inversion
    W(ind) = 1./dt_err_frac_thresh;
    % ind = find(W > maxerrweight);
    % W(ind) = maxerrweight;
    % ind = find(W < 1/fiterrtol);
    % W(ind) = 0;
    % for i=1:length(dt)
    %     W(i,i)=W(i,i).*(csnum(i).^0.5);
    % end
    para = polyfit(dist(:),dt,1);
    polyerr = polyval(para,dist(:)) - dt;
    errind = find(abs(polyerr) > polyfit_dt_err);
    for i = errind
        W(i,i) = 0;
    end
    
    if 0 % Set [1] to plot
        % Plot data weighting threshold
        figure(99); clf;
        w = diag(W);
        plot(dist,1./dt_err_frac,'.r'); hold on;
        plot(dist,w,'ob'); hold on;
        plot([min(dist),max(dist)],1./dt_err_frac_thresh*[1 1],'--g','linewidth',2);
        xlabel('Distance (km)');
        ylabel('Weighting');
        legend({'before threshold','after threshold','threshold'},'location','eastoutside')
        pause;
    end
    
    if is_wl_smooth
        wavelength = nanmean(grv) .* Tperiods(ip);
        wl_smooth_factor = wavelength ./ deg2km(gridsize);
    else
        wl_smooth_factor = 1;
    end

    % calculate the smoothing weight
    NR=norm(F,1);
    NA=norm(W*mat,1);
    smweight = smweight0*NA/NR * wl_smooth_factor;
    
    NR=norm(F_azi_damp,1);
    NA=norm(W*mat,1);
    damp_azi = damp0_azi*NA/NR;
    
    NR=norm(F_azic,1);
    NA=norm(W*mat,1);
    smweight_azi = smweight0_azi*NA/NR * wl_smooth_factor;
    
    NR=norm(J_azic,1);
    NA=norm(W*mat,1);
    flweight_azi = flweight0_azi*NA/NR * wl_smooth_factor;
    
    disp('start inverse');
    A=[W*mat; smweight*F; damp_azi*F_azi_damp; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];
    rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1); zeros(size(J_azic,1),1); zeros(size(J_azis,1),1)];

    groupg=(A'*A)\(A'*rhs);
    %        toc
    %        disp('Done');
    
    
    % Iteratively down weight the measurement with high error
    niter=1;
    
    while niter < 2
        niter=niter+1;
        err = mat*groupg - dt;

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

        if is_wl_smooth
            wavelength = 1./nanmean(groupg(1:Nx*Ny)) .* Tperiods(ip);
            wl_smooth_factor = wavelength ./ deg2km(gridsize);
        else
            wl_smooth_factor = 1;
        end
        
        % Rescale the smooth kernel
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0*NA/NR * wl_smooth_factor;
        
        NR=norm(F_azi_damp,1);
        NA=norm(W*mat,1);
        damp_azi = damp0_azi*NA/NR;
        
        NR=norm(F_azic,1);
        NA=norm(W*mat,1);
        smweight_azi = smweight0_azi*NA/NR * wl_smooth_factor;
        
        NR=norm(J_azic,1);
        NA=norm(W*mat,1);
        flweight_azi = flweight0_azi*NA/NR * wl_smooth_factor;
        
        % Invert
        A=[W*mat; smweight*F; damp_azi*F_azi_damp; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];
        rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1); zeros(size(J_azic,1),1); zeros(size(J_azis,1),1)];

        groupg=(A'*A)\(A'*rhs);

        
    end

    % Calculate model resolution and chi2
    Ginv = (A'*A)\mat'*W;
    R = Ginv * W*mat; % model resolution
    Rdiag = diag(R);
    [~,~,resol] = vec2mesh(ynode,xnode,Rdiag(1:Nx*Ny));
    % degrees of freedom
    v = length(dt) - trace(R);
    % normalized chi2 uncertainties
    res = (mat*groupg - dt);
    res(diag(W)==0) = nan;
    rms_res = sqrt(nanmean(res.^2));
    dt_std(dt_std<0.10*rms_res) = 0.10*rms_res;
    chi2 = nansum(res.^2./dt_std.^2)/v;
    
    % Calculate model uncertainties
    % slo_std = diag(Ginv*diag(rms_res.^2)*Ginv').^(1/2);
    % slo_std = diag(inv(A'*A)).^(1/2);
    A_dt = A;
    A_dt(1:length(dt),:) = diag(1./dt) * A(1:length(dt),:); % scale by dt to get units correct
    slo_std = diag(inv(A_dt'*A_dt)).^(1/2);
    % convert from dslow to dv
    grv_std = groupg(1:Nx*Ny).^(-2) .* slo_std(1:Nx*Ny);
    [~,~,GV_std] = vec2mesh(ynode,xnode,grv_std(1:Nx*Ny));
    
    % Model Rougness
    R2 = nanmean((F*groupg).^2);
    
    % Anisotropic terms from model vector
    groupg_azic = groupg(Nx*Ny+1 : Nx*Ny+Nx_azi*Ny_azi);
    groupg_azis = groupg(Nx*Ny+Nx_azi*Ny_azi+1 : Nx*Ny+Nx_azi*Ny_azi*2);
    
    % Isotropic group velocity
    grv_iso = dist'./(mat_iso*groupg(1:Nx*Ny));
    
    %        disp(' Get rid of uncertainty area');
    
    Igood = find(diag(W)~=0);
    mat_good = mat(Igood,:);
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            %raydense(i,j) = sum(mat(:,n));
            raydense(i,j) = sum(mat_good(:,n));
            if raydense(i,j) < raydensetol
                groupg(n)=NaN;
            end
        end
    end
    Igood = find(diag(W)~=0);
    mat_azi_hits_good = mat_azi_hits(Igood,:);
    for i=1:Nx_azi
        for j=1:Ny_azi 
            n=Ny_azi*(i-1)+j;
%             raydense_azi(i,j) = sum(mat_azi_hits(:,n));
            raydense_azi(i,j) = sum(mat_azi_hits_good(:,n));
            if raydense_azi(i,j) < raydensetol_azi
                groupg_azic(n)=NaN;
                groupg_azis(n)=NaN;
            end
        end
    end
    
    % Convert into group velocity
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            GV(i,j)= 1./groupg(n);
        end
    end
    for i=1:Nx_azi
        for j=1:Ny_azi
            n=Ny_azi*(i-1)+j;
            Ac(i,j) = -groupg_azic(n);
            As(i,j) = -groupg_azis(n);
        end
    end
    
    % JBR - Get Azimuthal coefficients from groupg (s/km)
    grv_av = nanmean(GV(:));
%     grv_av = nanmedian(GV(:));
    grv_avstd = nanstd(GV(:));
    slow_av = 1/grv_av;
%     slow_av = 1/mean(grv_iso);
    if sum(size(GV)==size(Ac)) == 2
        GV_aznode = GV;
    else
        GV_aznode = interp2(yi,xi,GV,yi_azi,xi_azi);
    end
    % get azimuthal anisotropy
    Ac2 = Ac.*GV_aznode; % s/km -> %
    As2 = As.*GV_aznode; % s/km -> %
    % Ac2 = Ac./slow_av;
    % As2 = As./slow_av;
    A2 = sqrt(Ac2.^2+As2.^2);
    phi2 = 1/2*atan2d(As2,Ac2);

    raytomo(ip).GV = GV;
    raytomo(ip).GV_std = GV_std;
    raytomo(ip).resol = resol;
    raytomo(ip).chi2 = chi2;
    raytomo(ip).res = res;
    raytomo(ip).mat = mat;
    raytomo(ip).raydense = raydense;
    raytomo(ip).period = Tperiods(ip);
    raytomo(ip).w = diag(W);
    raytomo(ip).err = err;
    raytomo(ip).rays = rays;
%     raytomo(ip).fiterr = fiterr;
    raytomo(ip).dt = dt;
    raytomo(ip).smweight0 = smweight0 * wl_smooth_factor;
    raytomo(ip).smweight0_azi = smweight0_azi * wl_smooth_factor;
    raytomo(ip).flweight0_azi = flweight0_azi * wl_smooth_factor;
    raytomo(ip).damp0_azi = damp0_azi;
    %JBR    
    raytomo(ip).grv_iso = grv_iso;    
    raytomo(ip).grv_av = grv_av;
    raytomo(ip).grv_avstd = grv_avstd;
    raytomo(ip).Ac2 = Ac2;
    raytomo(ip).As2 = As2;
%     raytomo(ip).Ac4 = Ac4;
%     raytomo(ip).As4 = As4;
    raytomo(ip).A2 = A2;
%     raytomo(ip).A4 = A4;
    raytomo(ip).phi2 = phi2;
%     raytomo(ip).phi4 = phi4;
    raytomo(ip).R2 = R2;
    raytomo(ip).grv = grv;
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

fig16 = figure(16);
set(gcf,'position',[1    1   1244   704],'color','w');
clf
vperc = [-r r];
ii = 0;
for ip=per_ind
    if sum(~isnan(raytomo(ip).GV(:))) == 0
        continue
    end
    ii = ii + 1;
    subplot(Mp,Np,ii)
    ax = worldmap(lalim, lolim);
    setm(gca,'MapProjection','mercator','FLineWidth',1.5,'FontSize',13)
    tightmap
    set(ax, 'Visible', 'on')
    set(gcf,'color','w')
    setm(gca,'FFaceColor',[0.9 0.9 0.9])
    A2 = raytomo(ip).A2;
%     surfacem(xi,yi,resid);
    levels = linspace(0,0.02,10)*100;
%     surfacem(xi_azi,yi_azi,A2*100,'Linestyle','none');
    pcolorm(xi_azi,yi_azi,A2*100,'Linestyle','none');
%     contourfm(xi_azi,yi_azi,A2*100,levels,'edgecolor','none');
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    caxis([min(levels) max(levels)])
    cb = colorbar;
    ylabel(cb,'Azimuthal Anisotropy (%)','fontsize',15);
    colormap('parula');
    
    scale = 20 * 4;
    u=raytomo(ip).A2 .* cosd(raytomo(ip).phi2)*scale;
	v=raytomo(ip).A2 .* sind(raytomo(ip).phi2)*scale;%./cosd(mean(lalim));
    u_bg=(raytomo(ip).A2+0.3/100) .* cosd(raytomo(ip).phi2)*scale;
	v_bg=(raytomo(ip).A2+0.3/100) .* sind(raytomo(ip).phi2)*scale;%./cosd(mean(lalim));
	[m n]=size(xi_azi);
    hold on;
    xpts = []; xpts_bg = [];
    ypts = []; ypts_bg = [];
	for ix=1:m
		for iy=1:n
            xpts = [xpts, [xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2]+gridsize_azi/2, nan];
            ypts = [ypts, [yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2]+gridsize_azi/2, nan];
            xpts_bg = [xpts_bg, [xi_azi(ix,iy)-u_bg(ix,iy)/2 xi_azi(ix,iy)+u_bg(ix,iy)/2]+gridsize_azi/2, nan];
            ypts_bg = [ypts_bg, [yi_azi(ix,iy)-v_bg(ix,iy)/2 yi_azi(ix,iy)+v_bg(ix,iy)/2]+gridsize_azi/2, nan];
        end
    end
%     plotm(xpts_bg,ypts_bg,'-','Color',[0 0 0],'linewidth',4);
    plotm(xpts,ypts,'-','Color',[0.9 0 0],'linewidth',2);
    hold on;
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
    % Plot reference
    refstick = scale*0.02;
    plotm([min(lalim) min(lalim)]+abs(diff(lalim))*0.15,[max(lolim)-refstick/2 max(lolim)+refstick/2]-abs(diff(lolim))*0.15,'-','Color',[0.9 0 0],'linewidth',2);
    textm(min(lalim)+abs(diff(lalim))*0.09,max(lolim)-abs(diff(lolim))*0.15,'2%','fontsize',12,'HorizontalAlignment', 'center');
end
save2pdf([grv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_raytomo_2Dazimuthal.pdf'],fig16,1000);


%% Group Velocity Maps (km/s)
% Load seafloor age
% load('age_grid.mat');

% Mp = 3; Np = 2;
fig17 = figure(17);
% set(gcf,'position',[94     1   599   704]);
set(gcf,'position',[1    1   1244   704],'color','w');
clf
ii = 0;
cmap = tomo_cmap(200);
for ip=per_ind
    if sum(~isnan(raytomo(ip).GV(:))) == 0
        continue
    end
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
    title([num2str((Tperiods(ip))),' s'],'fontsize',15)
    % text(0.05,0.85,[num2str(round(Tperiods(ip))),' s'],'fontsize',15,'fontweight','bold','Units','normalized','HorizontalAlignment','left');
    caxis([avgv*(1-r) avgv*(1+r)])
    cb = colorbar;
    ylabel(cb,'Group Velocity (km/s)','fontsize',15);
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
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
%     [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
end
save2pdf([grv_fig_path,'',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_raytomo.pdf'],fig17,1000);

%% Group Velocity Maps (%)
% Load seafloor age
% load('age_grid.mat');

fig19 = figure(19);
set(gcf,'position',[1    1   1244   704]);
clf
vperc = [-r r];
ii = 0;
for ip=per_ind
    if sum(~isnan(raytomo(ip).GV(:))) == 0
        continue
    end
    ii = ii + 1;
    subplot(Mp,Np,ii)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    set(gcf,'color',[0.9 0.9 0.9])
    avgv = nanmean(raytomo(ip).GV(:));
    resid = (raytomo(ip).GV-avgv)./avgv;
%     surfacem(xi,yi,resid);
    % levels = linspace(vperc(1),vperc(2),10)*100;
    levels = [-r:0.0025:r]*100;
    contourfm(xi,yi,resid*100,levels);
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
%     caxis([-r r])
    caxis([min(levels) max(levels)])
    cb = colorbar;
    ylabel(cb,'{\delta}c/c (%)','fontsize',15);
%     colormap(seiscmap)
    rbc = flip(redblue(length(levels)-1));
%     rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
    colormap(rbc);
    
    hold on;
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
%     [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
end
save2pdf([grv_fig_path,'',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_raytomo_perc.pdf'],fig19,1000);

%%
% MODEL RESOLUTION
fig21 = figure(21);
set(gcf,'position',[1    1   1244   704],'color','w');
clf
ii = 0;
for ip=per_ind
    if sum(~isnan(raytomo(ip).GV(:))) == 0
        continue
    end
    ii = ii + 1;
subplot(Mp,Np,ii); hold on;
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
%     surfacem(xi,yi,(raytomo(ip).resol));
    surfacem(xi,yi,(raytomo(ip).resol)./prctile(abs(raytomo(ip).resol(:)),99));
%     contourfm(xi,yi,(raytomo(ip).resol)./prctile(abs(raytomo(ip).resol(:)),99),[0:0.05:1],'edgecolor','none');
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    colorbar
    colormap(flip(hot));
%     caxis([0 prctile(abs(raytomo(ip).resol(:)),99)])
%     caxis([0 0.02]);
    caxis([0 0.5]);
%     RES = (raytomo(ip).resol)./max(abs(raytomo(ip).resol(:)));
%     contourm(xi,yi,RES,[0.1],'-g');
    plotm(sta.lat,sta.lon,'ob','markerfacecolor',[0 0 1]);
end
% save2pdf([grv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_ModelResolution.pdf'],fig21,1000);


% MODEL UNCERTAINTY
fig22 = figure(22);
set(gcf,'position',[1    1   1244   704],'color','w');
clf
ii = 0;
cmap = tomo_cmap(200);
for ip=per_ind
    if sum(~isnan(raytomo(ip).GV(:))) == 0
        continue
    end
    ii = ii + 1;
subplot(Mp,Np,ii); hold on;
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
%     surfacem(xi,yi,raytomo(ip).GV_std);
    surfacem(xi,yi,raytomo(ip).GV_std./raytomo(ip).GV*100);
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    colorbar
    colormap(cmap);
%     caxis([0 500])
    caxis([0 1]);
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
end
% save2pdf([grv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_ModelUncertainty.pdf'],fig22,1000);


%%
% RAY DENSITY
fig18 = figure(18);
set(gcf,'position',[94     1   599   704]);
clf
ii = 0;
for ip=per_ind
    ii = ii + 1;
    subplot(Mp,Np,ii)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).raydense ./ deg2km(gridsize));
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    cb = colorbar;
    ylabel(cb,'Hits','fontsize',15);
    colormap(flip(hot));
    caxis([0 20])
end
save2pdf([grv_fig_path,'',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_raydense.pdf'],fig18,1000);

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
    grv_av_rt(iper) = nanmean(raytomo(iper).grv_av);
    grv_avstd_rt(iper) = nanmean(raytomo(iper).grv_avstd);
end
% A2_rt = sqrt(Ac2.^2 + As2.^2)./grv_av_rt;
% A4_rt = sqrt(Ac4.^2 + As4.^2)./grv_av_rt;
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
% plot(periods,Ac2./grv_av_rt*100,'-k','linewidth',2); hold on;
% plot(periods,As2./grv_av_rt*100,'-r','linewidth',2);
% plot(periods,Ac4./grv_av_rt*100,'--k','linewidth',2);
% plot(periods,As4./grv_av_rt*100,'--r','linewidth',2);
% ylabel('A (%)','fontsize',18);
% title('Azimuthal Coefficients','fontsize',18);
% xlim(flip(1./frange));
% set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',18);
% xlabel('Period (s)','fontsize',18);
% legend({'A_{c2}','A_{s2}','A_{c4}','A_{s4}'},'fontsize',13,'box','off');

% save2pdf([grv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_A_phi_plots.pdf'],fig4,1000);

%% Plot group velocities
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
    plot(Tperiods,ftan.c_start,'ok','linewidth',2);
catch
    display('No starting data')
end
errorbar(Tperiods,grv_av_rt,grv_avstd_rt*2,'-r','linewidth',2);
% errorbar(Tperiods,mean([vertcat(raytomo(:).grv)],2),std([vertcat(raytomo(:).grv)],0,2)*2,'-b','linewidth',2);
title('Isotropic Group Velocity');
xlabel('Period (s)','fontsize',16);
ylabel('Group Velocity (km/s)','fontsize',16);
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
    resids = (grv_av_rt-aniso.c_iso)./aniso.c_iso*100;
    resids2 = (mean([vertcat(raytomo(:).grv)],2)'-aniso.c_iso)./aniso.c_iso*100;
    patch([aniso.periods fliplr(aniso.periods)],[aniso.err_c_iso./aniso.c_iso*2*100 fliplr(-aniso.err_c_iso./aniso.c_iso*2*100)],[0.8 0.8 0.8],'linestyle','none');
    errorbar(Tperiods,resids,grv_avstd_rt./grv_av_rt*100*2,'-r','linewidth',2);
    % errorbar(Tperiods,resids2,std([vertcat(raytomo(:).grv)],0,2)./mean([vertcat(raytomo(:).grv)],2)*2*100,'-b','linewidth',2)
    title('Residual (%)');
    xlabel('Period (s)','fontsize',16);
    ylabel('\delta c/c','fontsize',16);
    set(gca,'fontsize',16,'linewidth',1.5);
    xlim([3.5 10.5]);
    ylim([-2 2]);
end

save2pdf([grv_fig_path,comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_compareisogrv.pdf'],fig2,1000);