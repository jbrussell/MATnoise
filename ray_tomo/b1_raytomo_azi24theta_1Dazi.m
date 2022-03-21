% Script to do the ray theory tomography based on the ambient noise measurement. 
% Solves for azimuthal anisotropy on a 1D grid and isotropic velocity on a 2D grid.
%
% phv(theta,freq) = phv_iso(freq) + Ac2(freq)*cos(2*theta) + As2(freq)*sin(2*theta)
%                                 + Ac4(freq)*cos(4*theta) + As4(freq)*sin(4*theta)
%
% Written by Ge Jin, jinwar@gmail.com
% Nov 2012
% JBR: Modified in 2018 to include azimuthal anisotropy
% https://github.com/jbrussell
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
comp = parameters.comp; % {'ZZ'};
xspdir = parameters.xspdir; % 'phv_dir'; %'Nomelt3inttaper_iso.s0to333_br1avg'; %'4.0_S1_10pers_avg'; %'Nomelt3inttaper_iso.s0to333_br1avg'; %'4.0_S0_waverage';
windir = parameters.windir; %'window3hr'; 
N_wl = parameters.N_wl;
frange = parameters.frange; %[1/10 1/5]; % [Hz]
per_ind = parameters.per_ind; % [1:12]; % index of periods to consider

% QC parameters
snr_tol = parameters.snr_tol; % 3; % minimum signal-to-noise
is_rtolmin_wavelength = parameters.is_rtolmin_wavelength; % 0; 
wl_fac = parameters.wl_fac; % 1.0; % determine distance tolerance by wavelength?
r_tol_min = parameters.r_tol_min; %90; % [km] minimum station separation
r_tol_max = parameters.r_tol_max; %600; % [km] maximum station separation
err_tol = parameters.err_tol; %0.5; % maximum misfit of bessel fit between observed and synthetic
is_raydensity_thresh = parameters.is_raydensity_thresh; % Apply raydensity threshold to wipe out poorly constrained grid cells?
min_dep = parameters.min_dep; %= 9999; %-3500 for min station depth to use

% Damping parameters
azi_bin_deg = parameters.azi_bin_deg; % 20; % (degrees) size of azimuthal data bin
% Norm damping for azimuthal anisotropy
damp_azi = parameters.damp_azi; % [1 1 1e10 1e10]; % [2c 2s 4c 4s] % Damping individual parameters
aziweight = parameters.aziweight; %1; % global weight


%% Figure output
% figure output path
phv_fig_path = ['./figs/',windir,'/fullStack/raytomo_azi24theta_TEI19_1D/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/'];
if ~exist(phv_fig_path)    
    mkdir(phv_fig_path);
end

%%


% Load anisotropy data (from old inversion)
if iscompare_aniso
    load(['./aniso_DATA/',xspdir,'/',aniso_data]);
end

% Set up geometry parameters
setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
station_list = parameters.station_list;

% Load station info
[sta.name, sta.lat, sta.lon, sta.dep] = textread(station_list,'%s %f %f %f');

fiterrtol = parameters.fiterrtol;
maxerrweight = parameters.maxerrweight;
polyfit_dt_err = parameters.polyfit_dt_err;
smweight0 = parameters.smweight0;
dterrtol = parameters.dterrtol;
raydensetol = parameters.raydensetol;
if ~is_raydensity_thresh 
    raydensetol = 1;
end
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


%% Set up kernels, load data, and do the inversion
% Set up initial smoothing kernel
[i,j] = ndgrid(1:Nx,2:(Ny-1));
ind = j(:) + Ny*(i(:)-1);
dy = diff(ynode)*cosd(mean(xnode));  % correct smoothing for latitude;
dy1 = dy(j(:)-1);
dy2 = dy(j(:));
Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
    [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],Nx*Ny+4,Nx*Ny+4);
[i,j] = ndgrid(2:(Nx-1),1:Ny);
ind = j(:) + Ny*(i(:)-1);
dx = diff(xnode);
dx1 = dx(i(:)-1);
dx2 = dx(i(:));
Areg = [Areg;sparse(repmat(ind,1,3),[ind-Ny,ind,ind+Ny], ...
    [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],Nx*Ny+4,Nx*Ny+4)];
F=Areg;

%JBR - add anisotropic terms to smoothing matrix (add 4 columns and 4 rows)
% damping matrix
Areg_azi = zeros(4,Nx*Ny+4);
azipart = eye(4,4) .* diag(damp_azi);
Areg_azi(1:4,Nx*Ny+1:Nx*Ny+4) = azipart;
F_azi_damp = Areg_azi;

% smoothing matrix
F00 = zeros(4,Nx*Ny+4);
F00_part = 2*eye(4,4);
Fup = -1*[zeros(4-1,1) eye(4-1) ;zeros(1,4) ];
Fdown = -1*[zeros(1,4); eye(4-1) zeros(4-1,1) ];
F00_part = F00_part+Fup+Fdown;
% F00_part(1,:) = 0; F00_part(end,:) = 0;
F00(1:4,Nx*Ny+1:Nx*Ny+4) = F00_part;
F_azi_smooth = F00;


% Initialize the xsp structure
% Xsp_path = './Xsp/';
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

% cnt = 0;
% for ii = 1:length(xspsum)
%     cnt = cnt + xspsum(ii).isgood(1);
% end

% Loop through periods
for ip=1:length(Tperiods)
    disp(' ');
    disp(['Inversing Period: ',num2str(Tperiods(ip))]);
    clear rays dt fiterr mat phaseg err raydense dist azi mat_azi phv phv_std dt_std
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
        
        stas{raynum,1} = xspsum(ixsp).sta1;
        stas{raynum,2} = xspsum(ixsp).sta2;
        
        % dist(raynum) = deg2km(distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4)));
        dist(raynum) = distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4),referenceEllipsoid('GRS80'))/1000;
        dt(raynum) = xspsum(ixsp).tw(ip);
        phv(raynum) = dist(raynum)./dt(raynum);
        
        % convert uncertainty in velocity to uncertainty in time
        % dt = |r / v^2 * dv| = t^2 / r * dv
        phv_std(raynum,1) = xspsum(ixsp).c_std(ip);
        dt_std(raynum,1) = abs( dt(raynum).^2 / dist(raynum) * phv_std(raynum) );
        
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
        
        twloc = xspsum(ixsp).twloc;
        waxis = xspsum(ixsp).waxis;
        err = smooth((abs(xspsum(ixsp).err)./mean(abs(xspsum(ixsp).xsp_norm))).^2,round(length(waxis)/length(twloc)));
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
        mat_azi(raynum,:) = dist(raynum) * [cosd(2*azi(raynum)), sind(2*azi(raynum)), cosd(4*azi(raynum)), sind(4*azi(raynum)) ];
   
    end
    if size(dt,1) ~=raynum
        dt = dt';
    end
    
    % Building the isotropic data kernel
    disp('Start building the kernel');
    tic
    mat_iso=ray_kernel_build(rays,xnode,ynode);   
    toc
    
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
    smweight = smweight0;
    NR=norm(F,1);
    NA=norm(W*mat,1);
    smweight = smweight0*NA/NR;
    
    disp('start inverse');
%     A=[W*mat; smweight*F; aziweight*F_azi_damp; aziweight*F_azi_smooth];
%     rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1); zeros(size(F_azi_smooth,1),1)];
%     A=[W*mat;smweight*F;aziweight*F_azi_damp];
%     rhs=[W*dt;zeros(size(F,1),1);zeros(size(F_azi_damp,1),1)];
% %     A=[W*mat; smweight*F; aziweight*F_azi_smooth];
% %     rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_smooth,1),1)];
    A=[W*mat; smweight*F; aziweight*F_azi_damp];
    rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1)];

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
        A=[W*mat; smweight*F; aziweight*F_azi_damp];
        rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1)];

        phaseg=(A'*A)\(A'*rhs);

        
    end
    
    % Calculate model resolution and chi2
    Ginv = (A'*A)\mat'*W.^2;
    R = Ginv * mat; % model resolution
    Rdiag = diag(R);
    [~,~,resol] = vec2mesh(ynode,xnode,Rdiag(1:Nx*Ny));
    % degrees of freedom
    v = length(dt) - trace(R);
    % normalized chi2 uncertainties
    res = (mat*phaseg - dt);
    res(diag(W)==0) = nan;
    rms_res = sqrt(nanmean(res.^2));
    dt_std(dt_std<0.10*rms_res) = 0.10*rms_res;
    chi2 = nansum(res.^2./dt_std.^2)/v;
    
    % Calculate model uncertainties
    slo_std = diag(Ginv*diag(rms_res.^2)*Ginv').^(1/2);
    % convert from dslow to dv
    phv_std = phaseg(1:Nx*Ny).^(-2) .* slo_std(1:Nx*Ny);
    [~,~,GV_std] = vec2mesh(ynode,xnode,phv_std(1:Nx*Ny));
    
    % Model Rougness
    R2 = nanmean((F*phaseg).^2);
    
    % Isotropic phase velocity
    phv_iso = dist'./(mat_iso*phaseg(1:Nx*Ny));
    
    %        disp(' Get rid of uncertainty area');
    
    % Calculate ray density from "good" measurements (weight ~= 0)
    Igood = find(diag(W)~=0);
    mat_good = mat(Igood,:);
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
%             raydense(i,j) = sum(mat(:,n));
            raydense(i,j) = sum(mat_good(:,n));
            if raydense(i,j) < raydensetol
                phaseg(n)=NaN;
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
    
    % JBR - Get Azimuthal coefficients from phaseg (s/km)
    phv_av = 1./nanmean(1./GV(:)); %nanmean(GV(:));
%     phv_av = nanmedian(GV(:));
    phv_avstd = nanstd(GV(:));
    slow_av = 1/phv_av;
%     slow_av = 1/mean(phv_iso);
    Ac2 = -phaseg(Nx*Ny+1)./slow_av;
    As2 = -phaseg(Nx*Ny+2)./slow_av;
    Ac4 = -phaseg(Nx*Ny+3)./slow_av;
    As4 = -phaseg(Nx*Ny+4)./slow_av;
    A2 = sqrt(Ac2^2+As2^2);
    A4 = sqrt(Ac4^2+As4^2);
    phi2 = 1/2*atan2d(As2,Ac2);
    phi4 = 1/4*atan2d(As4,Ac4);

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
    raytomo(ip).stas = stas;
    raytomo(ip).sta_dep = dep;
    raytomo(ip).fiterr = fiterr;
    raytomo(ip).dt = dt;
    raytomo(ip).smweight0 = smweight0;
    %JBR    
    raytomo(ip).phv_iso = phv_iso;    
    raytomo(ip).phv_av = phv_av;
    raytomo(ip).phv_avstd = phv_avstd;
    raytomo(ip).Ac2 = Ac2;
    raytomo(ip).As2 = As2;
    raytomo(ip).Ac4 = Ac4;
    raytomo(ip).As4 = As4;
    raytomo(ip).A2 = A2;
    raytomo(ip).A4 = A4;
    raytomo(ip).phi2 = phi2;
    raytomo(ip).phi4 = phi4;
    raytomo(ip).R2 = R2;
    raytomo(ip).phv = phv;
    raytomo(ip).azi = azi;
    raytomo(ip).isgood = raytomo(ip).w~=0;
    
    % Do 1-D sinusoidal fitting on residuals
    dphv = (phv' - phv_iso) ./ phv_iso;
    varargin = sqrt(diag(1./W));
    isgood = raytomo(ip).isgood;
    [fitstr{ip}, A2_2(ip), phi2_2(ip)] = fit_azi_anisotropy2theta_resid(azi(isgood)',dphv(isgood),varargin(isgood));
    parastd{ip}=confint(fitstr{ip});    
    err_2p2p(ip) = parastd{ip}(2,1) - fitstr{ip}.d2;
    err_phi2(ip) = parastd{ip}(2,2) - fitstr{ip}.e2;
    raytomo(ip).dphv = dphv;
    raytomo(ip).A2_pts = A2_2(ip);
    raytomo(ip).phi2_pts = phi2_2(ip);
    raytomo(ip).A2_err = err_2p2p(ip);
    raytomo(ip).phi2_err = err_phi2(ip);
    
    % WEIGHTED RMS ERROR
    w = diag(W);
    w = w(isgood);
    dobs_2 = dphv(isgood); 
    dpre_2 = A2_2(ip)*cosd(2*(azi(isgood)'-phi2_2(ip)));
    wRMS_2A(ip) = sqrt(sum(w.*(dobs_2-dpre_2).^2)/sum(w));
    raytomo(ip).A2_wRMS = wRMS_2A(ip);
    
    % Average measurements in azimuthal bins to reduce sampling bias
    azi_good = azi(isgood);
    dphv_good = dphv(isgood);
    [~,Isort] = sort(azi_good);
    azi_srt = azi_good(Isort);
    dphv_srt = dphv_good(Isort);
    w_srt = w(Isort);
    bins = [-180:azi_bin_deg:180];
    dphv_bin = 0; azi_bin = 0; dphv_bin_err = 0; azi_bin_err = 0;
    for ibin = 1:length(bins)-1
        I_bin = azi_srt>=bins(ibin) & azi_srt<bins(ibin+1);
        wbin = w_srt(I_bin);
        wbinnorm = wbin/sum(wbin);
        % Weighted means and standard deviations
        dphv_bin(ibin) = nansum(wbin .* dphv_srt(I_bin)) / nansum(wbin);
        dphv_bin_err(ibin) = sqrt( var(dphv_srt(I_bin) , wbinnorm) );
        azi_bin(ibin) = (bins(ibin)+bins(ibin+1))/2;
%         azi_bin(ibin) = nansum(wbin' .* azi_srt(I_bin)) / nansum(wbin);
%         azi_bin_err(ibin) = sqrt( var(azi_srt(I_bin) , wbinnorm') );
    end
    raytomo(ip).dphv_bin = dphv_bin;
    raytomo(ip).dphv_bin_err = dphv_bin_err;
    raytomo(ip).azi_bin = azi_bin;
    raytomo(ip).azi_bin_err = azi_bin_err;
    varargin = dphv_bin_err;
    [fitstr{ip}, A2_2_bin(ip), phi2_2_bin(ip)] = fit_azi_anisotropy2theta_resid(azi_bin',dphv_bin',varargin');
    parastd{ip}=confint(fitstr{ip});    
    err_2p2p_bin(ip) = parastd{ip}(2,1) - fitstr{ip}.d2;
    err_phi2_bin(ip) = parastd{ip}(2,2) - fitstr{ip}.e2;
    raytomo(ip).A2_bin = A2_2_bin(ip);
    raytomo(ip).phi2_bin = phi2_2_bin(ip);
    raytomo(ip).A2_bin_err = err_2p2p_bin(ip);
    raytomo(ip).phi2_bin_err = err_phi2_bin(ip);
    % WEIGHTED RMS ERROR
    dobs_2 = dphv(isgood); 
    dpre_2 = A2_2_bin(ip)*cosd(2*(azi(isgood)'-phi2_2_bin(ip)));
    wRMS_2A_bin(ip) = sqrt(sum(w.*(dobs_2-dpre_2).^2)/sum(w));
    raytomo(ip).A2_bin_wRMS = wRMS_2A_bin(ip);
    
    
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

%% Phase Velocity Maps (km/s)
% Load seafloor age
% load('age_grid.mat');

% Mp = 3; Np = 2;
fig17 = figure(17);
% set(gcf,'position',[94     1   599   704]);
set(gcf,'Position',[94     1   719   704]);
clf
ii = 0;
cmap = tomo_cmap(200);
for ip=per_ind
    ii = ii + 1;
    subplot(Mp,Np,ii)
    ax = worldmap(lalim, lolim);
    setm(gca,'MapProjection','mercator','FLineWidth',1.5,'FontSize',13)
    tightmap
    set(ax, 'Visible', 'on')
    set(gcf,'color','w')
    setm(gca,'FFaceColor',[0.9 0.9 0.9])
%     set(gca,'Color',[0.7 0.7 0.7])
%     surfacem(xi,yi,raytomo(ip).GV);
    avgv = nanmean(raytomo(ip).GV(:));
    levels = linspace(avgv*(1-r), avgv*(1+r),100);
    contourfm(xi,yi,raytomo(ip).GV,levels,'edgecolor','none');
%     drawlocal
%     title([num2str(round(Tperiods(ip))),' s'],'fontsize',15)
    text(0.05,0.85,[num2str(round(Tperiods(ip))),' s'],'fontsize',15,'fontweight','bold','Units','normalized','HorizontalAlignment','left');
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
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
%     [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
end
save2pdf([phv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_raytomo_TEI19.pdf'],fig17,1000);

%% Phase Velocity Maps (%)
% Load seafloor age

% Mp = 3; Np = 4;
fig19 = figure(19);
set(gcf,'position',[94     1   599   704]);
clf
vperc = [-r r];
ii = 0;
for ip=per_ind
    ii = ii + 1;
    subplot(Mp,Np,ii)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    set(gcf,'color',[0.9 0.9 0.9])
    avgv = nanmean(raytomo(ip).GV(:));
    resid = (raytomo(ip).GV-avgv)./avgv;
%     surfacem(xi,yi,resid);
    levels = linspace(vperc(1),vperc(2),10)*100;
    contourfm(xi,yi,resid*100,levels);
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
%     caxis([-r r])
    caxis([min(levels) max(levels)])
    colorbar
%     colormap(seiscmap)
    rbc = flip(redbluecmap);
%     rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
    colormap(rbc);
    
    hold on;
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
%     [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
end
save2pdf([phv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_raytomo_TEI19_perc.pdf'],fig19,1000);

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
% save2pdf([phv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_ModelResolution.pdf'],fig21,1000);


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
% save2pdf([phv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_ModelUncertainty.pdf'],fig22,1000);


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
    surfacem(xi,yi,raytomo(ip).raydense);
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    colorbar
    colormap(flip(hot));
    caxis([0 500])
end
save2pdf([phv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_raydense_TEI19.pdf'],fig18,1000);

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

%% Ray Paths
fig20 = figure(20);
set(gcf,'position',[94     1   599   704]);
clf
% rbc = flip(redbluecmap);
% rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
ii = 0;
for ip=per_ind
    ii = ii + 1;
    cmap = repmat([0,0,0],length(raytomo(ip).rays(raytomo(ip).isgood,1)),1);
%     cmap = flipud(cmocean('balance',length(raytomo(ip).rays(raytomo(ip).isgood,1))));
%     cmap = flipud(cptcmap('GMT_no_green','ncol',length(raytomo(ip).rays(raytomo(ip).isgood,1))));
    subplot(Mp,Np,ii)
    ax = worldmap(lalim, lolim);
    setm(gca,'MapProjection','mercator','FLineWidth',1.5,'FontSize',13)
    tightmap
    set(ax, 'Visible', 'on')
    set(gcf,'color','w')
    
%     avgv = nanmean(dat(ip).phv_cor);
%     avgv = nanmean(dat(ip).phv);
%     vels = linspace(avgv*(1-r),avgv*(1+r),size(rbc,1));
    residvec = linspace(-2,2,size(cmap,1));
    clrs = [];
    lat1 = [];
    lon1 = [];
    lat2 = [];
    lon2 = [];
    resid = [];
    jj = 0;
    for ixsp = 1:length(raytomo(ip).phv)
        if ~raytomo(ip).isgood(ixsp) %|| raytomo(ip).sta_dep(ixsp) > -4500
            continue
        end
        jj = jj + 1;
        lat1(jj,1) = raytomo(ip).rays(ixsp,1);
        lon1(jj,1) = raytomo(ip).rays(ixsp,2);
        lat2(jj,1) = raytomo(ip).rays(ixsp,3);
        lon2(jj,1) = raytomo(ip).rays(ixsp,4);
        resid(jj,1) = (raytomo(ip).phv(ixsp)' - raytomo(ip).phv_iso(ixsp)) ./ raytomo(ip).phv_iso(ixsp) * 100;
%         [~,iclr] = min(abs(vels - dat(ip).phv_cor(ixsp)));
        [~,iclr] = min(abs(resid(jj) - residvec));
%         plotm([lat1 lat2],[lon1 lon2],dat(ip).phv(ixsp),'color',rbc(iclr,:),'linewidth',1.5);
%         clrs(jj,:) = cmap(iclr,:);
        clrs(jj,:) = cmap(iclr,:);
    %     drawlocal
    end
    h = plotm([lat1 lat2]',[lon1 lon2]','linewidth',1.5); hold on;
    set(h,{'color'},num2cell(clrs,2));
%     title([num2str(Tperiods(ip),'%.1f')],'fontsize',15)
    textm(36.5,-76.7,[num2str(round(Tperiods(ip))),' s'],'fontsize',15,'fontweight','bold');
%     colorbar
%     colormap(cmap);
%     caxis([residvec(1) residvec(end)]);
    
%     plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0.7 0.7 0.7]);
%     [c,h] = contourm(age_grid.LAT,age_grid.LON,age_grid.AGE,'k','LevelStep',5);
    drawnow;
end
save2pdf([phv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_raypaths.pdf'],fig20,1000);


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
    A2_rt(iper) = raytomo(iper).A2;
    A4_rt(iper) = raytomo(iper).A4;
    phi2_rt(iper) = raytomo(iper).phi2;
    phi4_rt(iper) = raytomo(iper).phi4;
    phv_av_rt(iper) = raytomo(iper).phv_av;
    phv_avstd_rt(iper) = raytomo(iper).phv_avstd;
    A2_pts(iper) = raytomo(iper).A2_pts;
    phi2_pts(iper) = raytomo(iper).phi2_pts;
    A2_pts_err(iper) = raytomo(iper).A2_err;
    phi2_pts_err(iper) = raytomo(iper).phi2_err;
    A2_pts_wRMS(iper) = raytomo(iper).A2_wRMS;
    
    % Binned estimates
    A2_bin(iper) = raytomo(iper).A2_bin;
    phi2_bin(iper) = raytomo(iper).phi2_bin;
    A2_bin_err(iper) = raytomo(iper).A2_bin_err;
    phi2_bin_err(iper) = raytomo(iper).phi2_bin_err;
    A2_bin_wRMS(iper) = raytomo(iper).A2_bin_wRMS;
   
end
phi2_bin(phi2_bin>150) = phi2_bin(phi2_bin>150)-180;
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
h3(2) = plot(periods,A4_rt*2*100,'-ob','linewidth',2);
h3(1) = plot(periods,A2_rt*2*100,'-o','color',[1 0 0],'linewidth',2); hold on;
% h3(1) = errorbar(periods,A2_pts*2*100,A2_pts_wRMS*100,'o','color',[1 0 0],'linewidth',2,'markerfacecolor',[1 0 0]);
% h3(1) = errorbar(periods,A2_bin*2*100,A2_bin_wRMS*100,'d','color',[0 0 0],'linewidth',1.5,'markerfacecolor',[0 0 0]);
% xlim(flip(1./frange));
% xlim([17 35]);
xlim([1./frange(2)-1 1./frange(1)+1]);
ylim([0 8]);
set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',18,'box','on');
xlabel('Period (s)','fontsize',18);
ylabel('Peak-to-peak amp (%)','fontsize',18);
% legend(h3,{'2\theta','4\theta'},'location','northwest');

% Azimuth
subplot(2,1,2); hold on;
ii = 0;
for iper = 1:length(phi2_rt)
    ii = ii + 1;
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
%         if phi2_rt(iper) < fastdir && dif > 10
%             phi2_rt(iper) = phi2_rt(iper)+180;
%         end
    elseif comp{1}(1) == 'T'
        [dif, I] = min(abs(phi2_vec-fastdir+90));
        phi2_rt(iper) = phi2_vec(I);
%         if phi2_rt(iper) < fastdir
%             phi2_rt(iper) = phi2_rt(iper)+180;
%         end
    end
%     phi2_rt(phi2_rt>160) = phi2_rt(phi2_rt>160)-180;


    phi4_vec(1) = phi4_rt(iper);
    phi4_vec(2) = phi4_rt(iper)+90;
    phi4_vec(3) = phi4_rt(iper)+180;
    phi4_vec(4) = phi4_rt(iper)+270;
    phi4_vec(5) = phi4_rt(iper)-90;
    phi4_vec(6) = phi4_rt(iper)-180;
    phi4_vec(7) = phi4_rt(iper)-270;
    [~, I] = min(abs(phi4_vec-fastdir+45));
    phi4_rt(iper) = phi4_vec(I);
    if phi4_rt(iper) < fastdir
        phi4_rt(iper) = phi4_rt(iper)+90;
    end
end
plot(periods,ones(size(periods))*fastdir,'--k','linewidth',2);
plot(periods,ones(size(periods))*fastdir+45,'--','color',[0.5 0.5 0.5]*0,'linewidth',2);
% plot(periods,ones(size(periods))*fastdir-90,'--','color',[0.5 0.5 0.5],'linewidth',2);
plot(periods,ones(size(periods))*fastdir+90,'--','color',[0.5 0.5 0.5]*0,'linewidth',2);
if iscompare_aniso
    errorbar(periods,phi4_2,err_phi4,'--ob','linewidth',2);
    errorbar(periods,phi2_2,err_phi2,'--o','color',[0 0.7 0],'linewidth',2);
end
plot(periods,phi4_rt,'-ob','linewidth',2);
plot(periods,phi2_rt,'-o','color',[1 0 0],'linewidth',2); hold on;
% errorbar(periods,phi2_pts,phi2_pts_err,'o','color',[1 0 0],'linewidth',2,'markerfacecolor',[1 0 0]);
% errorbar(periods,phi2_bin,phi2_bin_err,'d','color',[0 0 0],'linewidth',1.5,'markerfacecolor',[0 0 0]);
% errorbar(periods,phi2_bin-180,phi2_bin_err,'d','color',[0 0 0],'linewidth',1.5,'markerfacecolor',[0 0 0]);
% errorbar(periods,phi2_bin+180,phi2_bin_err,'d','color',[0 0 0],'linewidth',1.5,'markerfacecolor',[0 0 0]);
ylabel('Fast Direction (\circ)','fontsize',18);
ylim([fastdir-78 fastdir+130]);
% xlim(flip(1./frange));
% xlim([17 35]);
xlim([1./frange(2)-1 1./frange(1)+1]);
% ylim([-30 150]);
set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',18,'box','on');
xlabel('Period (s)','fontsize',18);

xs = get(gca,'xlim');
text(xs(2)+0.3,fastdir,'FSD','fontsize',18,'fontweight','bold');
text(xs(2)+0.3,fastdir-90,'Margin','color',[0.5 0.5 0.5],'fontsize',18,'fontweight','bold');


% PLOT COSINE AND SINE
fig5 = figure(5);
clf
plot(periods,Ac2./phv_av_rt*100,'-k','linewidth',2); hold on;
plot(periods,As2./phv_av_rt*100,'-r','linewidth',2);
plot(periods,Ac4./phv_av_rt*100,'--k','linewidth',2);
plot(periods,As4./phv_av_rt*100,'--r','linewidth',2);
ylabel('A (%)','fontsize',18);
title('Azimuthal Coefficients','fontsize',18);
xlim(flip(1./frange));
set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',18);
xlabel('Period (s)','fontsize',18);
legend({'A_{c2}','A_{s2}','A_{c4}','A_{s4}'},'fontsize',13,'box','off');

save2pdf([phv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_A_phi_plots_TEI19.pdf'],fig4,1000);

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

save2pdf([phv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_compareisophv_TEI19.pdf'],fig2,1000);

%% Plot Azimuthal Data
fig6 = figure(6); clf;
% set(gcf,'position',[10         248        1203         457]);
% set(gcf,'position',[94     1   599   704]);
set(gcf,'Position',[157    67   580   638]);
ii = 0;
for iper = per_ind
    ii = ii + 1;
    azi = raytomo(iper).azi;
    dep = raytomo(ip).sta_dep;
%     dphv = (raytomo(iper).phv - phv_av_rt(iper))./phv_av_rt(iper);
     dphv = (raytomo(iper).phv' - raytomo(iper).phv_iso) ./ raytomo(iper).phv_iso;
     
     % Binned values
    dphv_bin = raytomo(iper).dphv_bin;
    dphv_bin_err = raytomo(iper).dphv_bin_err;
    azi_bin = raytomo(iper).azi_bin;
    azi_bin_err = raytomo(iper).azi_bin_err;
     
    subplot(Mp,Np,ii); hold on;
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
    h2(2) = plot(x,A2_rt(iper)*cosd(2*(x-phi2_rt(iper)))*100+A4_rt(iper)*cosd(4*(x-phi4_rt(iper)))*100,'-','color',[1 0 0],'linewidth',5); hold on;
%     h2(2) = plot(x,A2_bin(iper)*cosd(2*(x-phi2_bin(iper)))*100,'-','color',[0 0 1],'linewidth',3); hold on;
%     h2(3) = plot(x,A2_pts(iper)*cosd(2*(x-phi2_pts(iper)))*100,'-','color',[0 1 0],'linewidth',3);
%     plot(azi,dphv*100,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
    
    isgood = raytomo(iper).isgood; %& dep' > -3000;
    plot(azi(isgood),dphv(isgood)*100,'o','color',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7],'markersize',5); hold on;
    errorbar(azi_bin,dphv_bin*100,dphv_bin_err*100,'dk','markerfacecolor',[0 0 0],'markersize',6,'linewidth',1.5);
%     plot(azi_bin,dphv_bin*100,'db','markersize',5,'linewidth',1);
    text(-170,7,[num2str(round(periods(iper))),' s'],'fontsize',16,'fontweight','bold','color',[0, 0, 0]);
    
    if iper == 4
        ax = get(gca);
        pos = ax.Position;
%         colorbar;
        set(gca,'Position',pos);
    end
    if sum(ii == [1 2 3 4]) > 0
        set(gca,'XTickLabel',[]);
    end
    if sum(ii == [2 4 6]) > 0
        set(gca,'YTickLabel',[]);
    end
    if iper == 1
%         legend(h2,{'2\theta + 4\theta'},'location','northwest');
    end
%     title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
    if ii == 5
        xlabel('Azimuth (\circ)','fontweight','bold','fontsize',15);
        vec_pos = get(get(gca, 'XLabel'), 'Position');
        set(get(gca, 'XLabel'), 'Position', vec_pos + [+230 -6 0]);
    end
    if ii == 3
        ylabel('\delta{c}/c (%)','fontweight','bold','fontsize',15);
    end
    set(gca,'fontsize',18,'linewidth',1.5,'box','on','TickDir','out','TickLength',[0.03, 0.005]);
    xlim([-180 180]);
    ylim([-8 8]);
    %ylim([3.8 4.8]);
    %box on;
    
    axpos = get(gca,'Position');
    set(gca,'Position',[axpos(1) axpos(2) axpos(3)*1.1 axpos(4)*1.2])
    
end

save2pdf([phv_fig_path,'TEI19_',comp{1}(1),'_','r',num2str(r_tol_min),'_',num2str(r_tol_max),'_snr',num2str(snr_tol),'_err',num2str(err_tol),'_sinplots_24theta_TEI19.pdf'],fig6,1000);
