% Read in the eventcs structures and apply eikonal tomography.
% Written by Ge Jin, jinwar@gmail.com
% 2013.1.16
%
% jbrussell 11/20/2020 : Invert all events simultaneously for 2-D isotropic 
% velocity and 2-D azimuthal anisotropy of the form:
% c(theta,freq) = c_iso(freq) + Ac2(freq)*cos(2*theta) + As2(freq)*sin(2*theta)
%
clear
%plot native

% setup parameters
setup_parameters
setup_ErrorCode

% JBR
% cohere_tol = parameters.cohere_tol;

is_offgc_propagation = parameters.is_offgc_propagation; % Account for off-great-circle propagation using eikonal tomography maps using eikonal results (a6_eikonal_eq)? Otherwise will assume great-circle propagation.
off_azi_tol = parameters.off_azi_tol; % [deg] tolerance for propagation off great-circle

% Smoothing parameters
flweight_array = 0*ones(length(parameters.periods)); %100*ones(length(parameters.periods)); %parameters.flweight_array
dterrtol = 2;    % largest variance of the inversion error allowed
inverse_err_tol = 2; %2  % count be number of standard devition
% Main QC parameters
fiterr_tol = 1e-2; % wavelet fitting error, throw out measurements greater than this
maxstadist = 600; % maximum epicentral distance difference
minstadist = 50; % minimum epicentral distance difference
cohere_tol = 0.80; % 0.65 % minimum coherence allowed
min_stadist_wavelength = 0.33; %0.5; % minimum number of wavelengths required between stations
max_stadist_wavelength = 999; % maximum number of wavelengths allowed between stations
ref_phv = 4*ones(1,length(parameters.periods)); % reference velocity for calculating wavelength
APM = 114; % absolute plate motion (GSRM 2.1; NNR) https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html
FSD = 75; 

% Norm damping for azimuthal anisotropy
% damp_azi = [1 1 1e10 1e10]; % [2c 2s 4c 4s] % Damping individual parameters
aziweight = 1*0; % global weight
smweight0_azi = 5 * 10;
flweight0_azi = 0; %0.1 * 1;

% debug setting
isfigure = 1;
isdisp = 1;
is_overwrite = 1;

% % input path
% eventcs_path = './CSmeasure/';
% % output path
% eikonl_ani_output_path = './eikonal/';

workingdir = parameters.workingdir;
% input path
eventcs_path = [workingdir,'CSmeasure/'];
% output path
eikonl_propazi_output_path = [workingdir,'eikonal_propazi/'];
eikonl_ani_output_path = [workingdir];


comp = parameters.component;
lalim=parameters.lalim;
lolim=parameters.lolim;
gridsize=parameters.gridsize;
gridsize_azi = parameters.gridsize_azi;
periods = parameters.periods;
raydensetol=parameters.raydensetol;
smweight_array = parameters.smweight_array;
% flweight_array = parameters.flweight_array; % JBR
Tdumpweight0 = parameters.Tdumpweight;
Rdumpweight0 = parameters.Rdumpweight;
fiterrtol = parameters.fiterrtol;
% dterrtol = parameters.dterrtol;
isRsmooth = parameters.isRsmooth;
% inverse_err_tol = parameters.inverse_err_tol;
min_amp_tol  = parameters.min_amp_tol;

% setup useful variables
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi, yi]=ndgrid(xnode,ynode);

xnode_azi=lalim(1):gridsize_azi:lalim(2);
ynode_azi=lolim(1):gridsize_azi:lolim(2);
Nx_azi=length(xnode_azi);
Ny_azi=length(ynode_azi);
[xi_azi, yi_azi]=ndgrid(xnode_azi,ynode_azi);

%% Setup Smoothing and damping kernels
% Isotropic velocity second derivative smoothing
disp('initial the smoothing kernel')
Nxyazi = Nx_azi*Ny_azi*2;
Nxy = Nx*Ny*2;
Fsm = smooth_kernel_build(xnode, ynode, Nx*Ny);
F=sparse(Nxy*2,Nxy+Nxyazi);
for n=1:size(Fsm,1)
    ind=find(Fsm(n,:)~=0);
    F(2*n-1,2*ind-1)=Fsm(n,ind);
    F(2*n,2*ind)=Fsm(n,ind);
end
% Azimuthal second derivative smoothing kernels
F_azi = smooth_kernel_build(xnode_azi, ynode_azi, Nx_azi*Ny_azi);
F_azic = sparse(Nxyazi,Nxy+Nxyazi);
F_azic(1:end,Nxy+[1:Nx_azi*Ny_azi]) = F_azi;
F_azis = sparse(Nxyazi,Nxy+Nxyazi);
F_azis(1:end,Nxy+Nx_azi*Ny_azi+[1:Nx_azi*Ny_azi]) = F_azi;

% Isotropic velocity first derivative "flatness" kernel
F1 = flat_kernel_build2(xnode, ynode, Nx*Ny);
Ftemp = sparse(Nxy*2,Nxy+Nxyazi);
for n=1:size(F1,1)
    ind=find(F1(n,:)~=0);
    Ftemp(2*n-1,2*ind-1)=F1(n,ind);
    Ftemp(2*n,2*ind)=F1(n,ind);
end
F1 = Ftemp;

% Azimuthal first derivative flattening kernels
J_azi = flat_kernel_build2(xnode_azi, ynode_azi, Nx_azi*Ny_azi);
J_azic = sparse(Nxyazi,Nxy+Nxyazi);
J_azic(1:end,Nxy+[1:Nx_azi*Ny_azi]) = J_azi;
J_azis = sparse(Nxyazi,Nxy+Nxyazi);
J_azis(1:end,Nxy+Nx_azi*Ny_azi+[1:Nx_azi*Ny_azi]) = J_azi;

% Norm damping of anisotropy
Areg_azi = zeros(Nxyazi,Nxy+Nxyazi);
azipart = eye(Nxyazi,Nxyazi); %.* diag(damp_azi);
Areg_azi(1:Nxyazi,Nxy+1:Nxy+Nxyazi) = azipart;
F_azi_damp = Areg_azi;

%% read in bad station list, if existed
if exist('badsta.lst')
	badstnms = textread('badsta.lst','%s');
	disp('Found Bad stations:')
	disp(badstnms)
end

csmatfiles = dir([eventcs_path,'/*cs_',comp,'.mat']);
% Loop through the periods
clear eventphv_ani;
for ip = 1:length(periods)
	clear mat_azi azi rays ddist dt w azi_vec
    disp(periods(ip));
	smweight0 = smweight_array(ip);
	flweight0 = flweight_array(ip); % JBR
 	raynum = 0;
	for ie = 1:length(csmatfiles)
	%for ie = 30
		% read in data and set up useful variables
		temp = load([eventcs_path,csmatfiles(ie).name]);
		eventcs =  temp.eventcs;
% 		disp(eventcs.id)
		evla = eventcs.evla;
		evlo = eventcs.evlo;

		matfilename = [eikonl_ani_output_path,'/eikonal_ani_',comp,'.mat'];
		if exist(matfilename,'file') && ~is_overwrite
            disp(['Exist ',matfilename,', skip']);
			continue;
        end
        if is_offgc_propagation==1
            eikonal_in = [eikonl_propazi_output_path,'/',eventcs.id,'_eikonal_',comp,'.mat'];
            if ~exist(eikonal_in,'file')
                error('No propagation azimuth found. Need to first run a6_a0_eikonal_eq_GetPropAzi.m');
            end
            temp = load(eikonal_in);
            phase_lat = temp.eventphv(ip).GVx; % phase slowness in x-direction
            phase_lon = temp.eventphv(ip).GVy; % phase slowness in y-direction
                        
            % Use event eikonal tomography results to get propagation azimuth
            azimat = 90 - atan2d(phase_lat,phase_lon);
            azimat(azimat<0) = azimat(azimat<0) + 360;
            [~, azimat_ev] = distance(xi,yi,evla,evlo,referenceEllipsoid('GRS80'));
            azimat(isnan(azimat)) = azimat_ev(isnan(azimat));
            % Ensure that propagation azimuth is not too far from great-circle
            diff_az = angdiff(azimat*pi/180,azimat_ev*pi/180)*180/pi;
            if mean(abs(diff_az(:))) > off_azi_tol
                % disp('off_azi_tol exceeded... skipping');
                continue
            end
        else
            [~, azimat] = distance(xi,yi,evla,evlo,referenceEllipsoid('GRS80'));               
        end

		if exist('badstnms','var')
			badstaids = find(ismember(eventcs.stnms,badstnms));
		else
			badstaids = [];
		end

		% Build the rotation matrix
		razi = azimuth(xi+gridsize/2,yi+gridsize/2,evla,evlo,referenceEllipsoid('GRS80'))+180;
		% R = sparse(Nxy+Nxyazi,Nxy+Nxyazi);
		% for i=1:Nx
		% 	for j=1:Ny
		% 		n=Ny*(i-1)+j;
		% 		theta = razi(i,j);
		% 		R(2*n-1,2*n-1) = cosd(theta);
		% 		R(2*n-1,2*n) = sind(theta);
		% 		R(2*n,2*n-1) = -sind(theta);
		% 		R(2*n,2*n) = cosd(theta);
		% 	end
		% end

		% Calculate the relative travel time compare to one reference station
% 		travel_time = Cal_Relative_dtp(eventcs);

		% Build the ray locations
		for ics = 1:length(eventcs.CS)
			raynum = raynum+1;
			
			if eventcs.CS(ics).cohere(ip)<cohere_tol && eventcs.CS(ics).isgood(ip)>0
                eventcs.CS(ics).isgood(ip) = ErrorCode.low_cohere;
            end
			if (abs(eventcs.CS(ics).ddist) < ref_phv(ip)*periods(ip)*min_stadist_wavelength || ...
                abs(eventcs.CS(ics).ddist) > ref_phv(ip)*periods(ip)*max_stadist_wavelength) && eventcs.CS(ics).isgood(ip)>0
				eventcs.CS(ics).isgood(ip) = ErrorCode.min_stadist_wavelength;
            end
            if (abs(eventcs.CS(ics).ddist) > maxstadist || ...
                abs(eventcs.CS(ics).ddist) < minstadist) && eventcs.CS(ics).isgood(ip)>0
                eventcs.CS(ics).isgood(ip) = -13;
            end
            if (eventcs.CS(ics).fiterr(ip) > fiterr_tol)
                eventcs.CS(ics).isgood(ip) = -14;
            end
            if eventcs.CS(ics).cohere(ip)>=cohere_tol && eventcs.CS(ics).isgood(ip)==ErrorCode.low_cohere
                eventcs.CS(ics).isgood(ip) = 1;
            end
			if eventcs.CS(ics).isgood(ip) > 0 
				dt(raynum,:) = eventcs.CS(ics).dtp(ip);
				w(raynum,:) = 1;
%                 w(raynum,:) = eventcs.CS(ics).fiterr(ip)^(-1/2);
%                 w(raynum,:) = eventcs.CS(ics).fiterr(ip)^(-1/2) * eventcs.CS(ics).ddist^(1/2);
			else
				dt(raynum,:) = eventcs.CS(ics).dtp(ip);
				w(raynum,:) = 0;
			end
			if sum(ismember([eventcs.CS(ics).sta1 eventcs.CS(ics).sta2],badstaids)) > 0
				w(raynum,:) = 0;
            end
            
            evinfo.evid{raynum,:} = eventcs.id;
            evinfo.evla(raynum,:) = eventcs.evla;
            evinfo.evlo(raynum,:) = eventcs.evlo;
            evinfo.evdp(raynum,:) = eventcs.evdp;
            evinfo.Mw(raynum,:) = eventcs.Mw;
			
			rays(raynum,1) = eventcs.stlas(eventcs.CS(ics).sta1);
			rays(raynum,2) = eventcs.stlos(eventcs.CS(ics).sta1);
			rays(raynum,3) = eventcs.stlas(eventcs.CS(ics).sta2);
			rays(raynum,4) = eventcs.stlos(eventcs.CS(ics).sta2);
		
			% JBR - Build azimuthal part of data kernel
			ddist(raynum,:) = eventcs.CS(ics).ddist;
	        % dt(ics) = eventcs.dtp(ics);
	        % phv(ics) = ddist(ics)./dt(ics);
            mean_stala = mean([eventcs.stlas(eventcs.CS(ics).sta1), eventcs.stlas(eventcs.CS(ics).sta2)]);
            mean_stalo = mean([eventcs.stlos(eventcs.CS(ics).sta1), eventcs.stlos(eventcs.CS(ics).sta2)]);
            
            for i=1:Nx
                for j=1:Ny
                    n=Ny*(i-1)+j;
                    azi(raynum,n)=azimat(i,j);
                end
            end
        end
    end
    azi(azi<0) = azi(azi<0) + 360;
    W = sparse(length(w),length(w));
    for id = 1:length(w)
        if w(id) > 0
            W(id,id) = w(id);
        end
    end

    % Build the kernel
    disp('Building up ray path kernel')
    tic
        % mat_iso=kernel_build(rays,xnode,ynode);
        mat_iso=kernel_build_RT(rays,xnode,ynode,azi);

        % Build 2-D azimuthal part of kernel
        [mat_azi, mat_azi_hits] = kernel_build_azi(rays,xnode_azi,ynode_azi,azi);
    toc

    mat = [mat_iso, mat_azi];

    % build dumping matrix for ST
    % dumpmatT = R(2:2:2*Nx*Ny,:);

    % build dumping matrix for SR
    % dumpmatR = R(1:2:2*Nx*Ny-1,:);


    % Normalize smoothing kernel
    NR=norm(F,1);
    NA=norm(W*mat,1);
    smweight = smweight0*NA/NR;

    % JBR - Normalize flatness kernel
    NR=norm(F1,1);
    NA=norm(W*mat,1);
    flweight = flweight0*NA/NR;

    % % Normalize dumping matrix for ST
    % NR=norm(dumpmatT,1);
    % NA=norm(W*mat,1);
    % dumpweightT = Tdumpweight0*NA/NR;
    % 
    % % Normalize dumping matrix for SR
    % NR=norm(dumpmatR,1);
    % NA=norm(W*mat,1);
    % dumpweightR = Rdumpweight0*NA/NR;

    % Rescale azimuthal anisotropy damping
    NR=norm(F_azi_damp,1);
    NA=norm(W*mat,1);
    aziweight0 = aziweight*NA/NR;

    %Rescale azimuthal anisotropy smoothness
    NR=norm(F_azic,1);
    NA=norm(W*mat,1);
    smweight_azi = smweight0_azi*NA/NR;

    %Rescale azimuthal anisotropy flatness
    NR=norm(J_azic,1);
    NA=norm(W*mat,1);
    flweight_azi = flweight0_azi*NA/NR;

    % Set up matrix on both side
    % A=[W*mat;smweight*F;flweight*F1;dumpweightT*dumpmatT;dumpweightR*dumpmatR;aziweight0*F_azi_damp; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];
    A=[W*mat;smweight*F;flweight*F1;aziweight0*F_azi_damp; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];

    avgv = eventcs.avgphv(ip);
    % rhs=[W*dt;zeros(size(F,1),1);zeros(size(F1,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
    % rhs=[W*dt;zeros(size(F,1),1);zeros(size(F1,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv;zeros(size(F_azi_damp,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1); zeros(size(J_azic,1),1); zeros(size(J_azis,1),1)];
    rhs=[W*dt;zeros(size(F,1),1);zeros(size(F1,1),1);zeros(size(F_azi_damp,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1); zeros(size(J_azic,1),1); zeros(size(J_azis,1),1)];


    % Least square inversion
    phaseg=(A'*A)\(A'*rhs);

    % Iteratively down weight the measurement with high error
    niter=0;
    ind = find(diag(W)==0);
    if isdisp
        disp(['Before iteration'])
        disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
        disp(['Bad Measurement Number: ', num2str(length(ind))]);
    end
    niter=1;
    while niter < 2
        niter=niter+1;
        isgood_mat = diag(diag(W)>0);
        err = mat*phaseg - dt;
% 			err = W*err;
        err = isgood_mat*err;
%             err = W*err;
        % stderr=std(err);
        stderr=std(err(err~=0));
        if stderr > dterrtol
            stderr = dterrtol;
        end
        for i=1:length(err)
            if abs(err(i)) > inverse_err_tol*stderr  || abs(err(i))==0
                W(i,i)=0;
            else
                W(i,i)=1./stderr;
            end
        end        
        ind = find(diag(W)==0);
        if isdisp
            disp('After:')
            disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
            disp(['Bad Measurement Number: ', num2str(length(ind))]);
        end

        % Rescale the smooth kernel
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0*NA/NR;

        % JBR - Normalize flatness kernel
        NR=norm(F1,1);
        NA=norm(W*mat,1);
        flweight = flweight0*NA/NR;

        % % rescale dumping matrix for St
        % NR=norm(dumpmatT,1);
        % NA=norm(W*mat,1);
        % dumpweightT = Tdumpweight0*NA/NR;
        % 
        % % rescale dumping matrix for SR
        % NR=norm(dumpmatR,1);
        % NA=norm(W*mat,1);
        % dumpweightR = Rdumpweight0*NA/NR;

        % Rescale azimuthal anisotropy damping
        NR=norm(F_azi_damp,1);
        NA=norm(W*mat,1);
        aziweight0 = aziweight*NA/NR;

        %Rescale azimuthal anisotropy smoothness
        NR=norm(F_azic,1);
        NA=norm(W*mat,1);
        smweight_azi = smweight0_azi*NA/NR;

        %Rescale azimuthal anisotropy flatness
        NR=norm(J_azic,1);
        NA=norm(W*mat,1);
        flweight_azi = flweight0_azi*NA/NR;

        % A=[W*mat;smweight*F;flweight*F1;dumpweightT*dumpmatT;dumpweightR*dumpmatR;aziweight0*F_azi_damp; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];
        A=[W*mat;smweight*F;flweight*F1;aziweight0*F_azi_damp; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];

        % rhs=[W*dt;zeros(size(F,1),1);zeros(size(F1,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
        rhs=[W*dt;zeros(size(F,1),1);zeros(size(F1,1),1);zeros(size(F_azi_damp,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1); zeros(size(J_azic,1),1); zeros(size(J_azis,1),1)];
        phaseg=(A'*A)\(A'*rhs);
    end	
    
    % Estimate travel-time residuals
    dt_res = dt - mat*phaseg;

    % Calculate the kernel density
    %sumG=sum(abs(mat),1);
    ind=1:Nx*Ny;
    rayW = W;
    rayW(find(rayW>1))=1;
    raymat = rayW*mat;
    sumG(ind)=sum((raymat(:,2*ind).^2+raymat(:,2*ind-1).^2).^.5,1);
    clear raydense
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            raydense(i,j)=sumG(n);
        end
    end

    %        disp(' Get rid of uncertainty area');
    fullphaseg = phaseg;
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            if raydense(i,j) < raydensetol %&& ~issyntest
                phaseg(2*n-1)=NaN;
                phaseg(2*n)=NaN;
            end
        end
    end

    % Isotropic phase velocity
    phv_iso = ddist./(mat_iso*phaseg(1:Nx*Ny*2));

    % Change phaseg into phase velocity
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            GVx(i,j)= phaseg(2*n-1);
            GVy(i,j)= phaseg(2*n);
        end
    end

    % Anisotropic terms from model vector
    phaseg_azic = phaseg(Nx*Ny*2+[1:Nx_azi*Ny_azi]);
    phaseg_azis = phaseg(Nx*Ny*2+Nx_azi*Ny_azi+[1:Nx_azi*Ny_azi]);
    for i=1:Nx_azi
        for j=1:Ny_azi 
            n=Ny_azi*(i-1)+j;
            raydense_azi(i,j) = sum(mat_azi_hits(:,n));
            if raydense_azi(i,j) < raydensetol
                phaseg_azic(n)=NaN;
                phaseg_azis(n)=NaN;
            end
        end
    end
    [~,~,Ac] = vec2mesh(ynode_azi,xnode_azi,-phaseg_azic);
    [~,~,As] = vec2mesh(ynode_azi,xnode_azi,-phaseg_azis);

    GV=(GVx.^2+GVy.^2).^-.5;
    GV_azi = angle(GVx + GVy.*sqrt(-1));
    GV_azi = rad2deg(GV_azi);
    GV_av = 1./nanmean(1./GV(:));
    slow_av = 1./GV_av;
    if sum(size(GV)==size(Ac)) == 2
        GV_aznode = GV;
    else
        GV_aznode = interp2(yi,xi,GV,yi_azi,xi_azi);
    end
    % get azimuthal anisotropy
    Ac2 = Ac.*GV_aznode; % s/km -> %
    As2 = As.*GV_aznode; % s/km -> %
    A2 = sqrt(Ac2.^2+As2.^2);
    phi2 = 1/2*atan2d(As2,Ac2);
    phi2(phi2<0) = phi2(phi2<0)+180;

    % save the result in the structure
    eventphv_ani(ip).rays = rays;
    eventphv_ani(ip).w = diag(W);
    eventphv_ani(ip).goodnum = length(find(eventphv_ani(ip).w>0));
    eventphv_ani(ip).badnum = length(find(eventphv_ani(ip).w==0));
    eventphv_ani(ip).dt = dt;
    eventphv_ani(ip).dt_res = dt_res; % data residuals
    eventphv_ani(ip).GV = GV;
    eventphv_ani(ip).GV_azi = GV_azi;
    eventphv_ani(ip).GV_av = GV_av;
    eventphv_ani(ip).GVx = GVx;
    eventphv_ani(ip).GVy = GVy;
    eventphv_ani(ip).phv_iso = phv_iso;
    eventphv_ani(ip).raydense = raydense;
    eventphv_ani(ip).lalim = lalim;
    eventphv_ani(ip).lolim = lolim;
    eventphv_ani(ip).gridsize = gridsize;
    eventphv_ani(ip).gridsize_azi = gridsize_azi;
    eventphv_ani(ip).azi = azi;
    eventphv_ani(ip).rays = rays; 
    eventphv_ani(ip).ddist = ddist;
    eventphv_ani(ip).id = evinfo.evid;
    eventphv_ani(ip).evla = evinfo.evla;
    eventphv_ani(ip).evlo = evinfo.evlo;
    eventphv_ani(ip).evdp = evinfo.evdp;
    eventphv_ani(ip).Mw = evinfo.Mw;
    eventphv_ani(ip).period = periods(ip);
    % eventphv_ani(ip).traveltime = travel_time(ip).tp;
    % eventphv_ani(ip).stlas = eventcs.stlas;
    % eventphv_ani(ip).stlos = eventcs.stlos;
    % eventphv_ani(ip).stnms = eventcs.stnms;
    eventphv_ani(ip).isgood = eventphv_ani(ip).w>0;
    eventphv_ani(ip).phv = ddist./dt;
    eventphv_ani(ip).A2 = A2;
    eventphv_ani(ip).phi2 = phi2;

    disp(['Period:',num2str(periods(ip)),', Goodnum:',num2str(eventphv_ani(ip).goodnum),...
            'Badnum:',num2str(eventphv_ani(ip).badnum)]);
end % end of periods loop

matfilename = [eikonl_ani_output_path,'/eikonal_ani2D_',comp,'.mat'];
save(matfilename,'eventphv_ani');
disp(['Save the result to: ',matfilename])
    
%% Make event list
fid = fopen([workingdir,'azi2D_evlist.txt'],'w');

for ip = 7 % 50 s
    isgood = eventphv_ani(ip).isgood;
    evids_all = eventphv_ani(ip).id(isgood);
    [evids,I] = unique(evids_all);
    evlas = eventphv_ani(ip).evla(isgood); evlas = evlas(I);
    evlos = eventphv_ani(ip).evlo(isgood); evlos = evlos(I);
    evdps = eventphv_ani(ip).evdp(isgood); evdps = evdps(I);
    evMws = eventphv_ani(ip).Mw(isgood);   evMws = evMws(I);
    for iev = 1:length(evids)
        evid = evids{iev};
        evla = evlas(iev);
        evlo = evlos(iev);
        evdp = evdps(iev);
        evMw = evMws(iev);
        fprintf(fid,'%13s %12f %12f %5.1f %5.2f\n',evid,evla,evlo,evdp,evMw);
    end
%     display(length(evids))
end
fclose(fid);
    
 
 %% Phase velocity
 
 N=3; M = floor(length(periods)/N) +1;
figure(88)
clf
for ip = 1:length(periods)
    subplot(M,N,ip)
    ax = worldmap(lalim+[-gridsize_azi +gridsize_azi], lolim+[-gridsize_azi +gridsize_azi]);
    set(ax, 'Visible', 'off')
    h1=surfacem(xi,yi,eventphv_ani(ip).GV);
    % set(h1,'facecolor','interp');
%			load pngcoastline
%			geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',2)
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
    avgv = nanmean(eventphv_ani(ip).GV(:));
    if isnan(avgv)
        continue;
    end
    r = 0.1;
    caxis([avgv*(1-r) avgv*(1+r)])
    colorbar
    load seiscmap
    colormap(seiscmap)

%             % Plot ray paths
%             jj=0;
%             for ii = 1:length(eventphv_ani(ip).w)
%                 if ~eventphv_ani(ip).isgood(ii) %|| raytomo(ip).sta_dep(ixsp) > -4500
%                     continue
%                 end
%                 jj = jj + 1;
%                 lat1(jj,1) = eventphv_ani(ip).rays(ii,1);
%                 lon1(jj,1) = eventphv_ani(ip).rays(ii,2);
%                 lat2(jj,1) = eventphv_ani(ip).rays(ii,3);
%                 lon2(jj,1) = eventphv_ani(ip).rays(ii,4);
%             end
%             hold on;
%             h = plotm([lat1 lat2]',[lon1 lon2]','-k','linewidth',0.5); hold on;

end
drawnow;
 
%% Anisotropy maps
figure(57)
clf
for ip = 1:length(periods)
	subplot(M,N,ip);
	ax = worldmap(lalim+[-gridsize_azi +gridsize_azi], lolim+[-gridsize_azi +gridsize_azi]);
	set(ax, 'Visible', 'off')
    title([num2str(periods(ip)),' s'])
% 	h1=surfacem(xi,yi,avgphv_aniso(ip).isophv);
% 	h1=surfacem(xi,yi,avgphv_aniso(ip).aniso_strength);
    % Plot zero-to-peak anisotropy strength
    h1=pcolorm(xi_azi-gridsize_azi/2,yi_azi-gridsize_azi/2,eventphv_ani(ip).A2*100,'Linestyle','none');
	colorbar
	colormap(parula)
	drawnow
%     caxis([0 0.05]);

    scale = 50;
    u=eventphv_ani(ip).A2.*cosd(eventphv_ani(ip).phi2)*scale;
    v=eventphv_ani(ip).A2.*sind(eventphv_ani(ip).phi2)*scale;%./cosd(mean(lalim));
    [m n]=size(xi_azi);
    hold on;
    xpts = [];
    ypts = [];
    for ix=1:m
        for iy=1:n
            xpts = [xpts, [xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2]+gridsize_azi/2, nan];
            ypts = [ypts, [yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2]+gridsize_azi/2, nan];
        end
    end
    %     plotm(xpts_bg,ypts_bg,'-','Color',[0 0 0],'linewidth',4);
    plotm(xpts,ypts,'-','Color',[0.9 0 0],'linewidth',2);
    hold on;
    % Plot reference
    %     refstick = scale*0.02;
    %     plotm([min(lalim) min(lalim)]+abs(diff(lalim))*0.15,[max(lolim)-refstick/2 max(lolim)+refstick/2]-abs(diff(lolim))*0.15,'-','Color',[0.9 0 0],'linewidth',2);
    %     textm(min(lalim)+abs(diff(lalim))*0.09,max(lolim)-abs(diff(lolim))*0.15,'2%','fontsize',12,'HorizontalAlignment', 'center');
end

 %%
figure(58);
set(gcf,'position',[351   677   560   668]);
clf
clear avgv avgv_std aniso_str aniso_str_std aniso_azi aniso_azi_std
for ip = 1:length(periods)
    avgv(ip) = eventphv_ani(ip).GV_av;
    avgv_std(ip) = nanstd(eventphv_ani(ip).GV(:));
    aniso_str(ip) = nanmean(eventphv_ani(ip).A2(:));
    aniso_str_std(ip) = nanstd(eventphv_ani(ip).A2(:));
    phi2 = eventphv_ani(ip).phi2(:);
    phi2(phi2<0) = phi2(phi2<0) + 180;
    aniso_azi(ip) = nanmean(phi2(:));
    aniso_azi_std(ip) = nanstd(phi2(:));
end
%plot native
subplot(3,1,1); hold on;
errorbar(periods,avgv,avgv_std*2,'-or');
ylim([3.85 4.4]);
xlim([20 150]);
ylabel('c (km/s)');
%plot native
subplot(3,1,2); hold on;
errorbar(periods,aniso_str*100*2,aniso_str_std*100*2,'-or');
ylim([0 5]);
xlim([20 150]);
ylabel('2A');
%plot native
subplot(3,1,3);
plot([periods(1),periods(end)],FSD*[1 1],'--k','linewidth',1.5); hold on;
plot([periods(1),periods(end)],APM*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1.5);
errorbar(periods,aniso_azi,aniso_azi_std*2,'-or');
errorbar(periods,aniso_azi+180,aniso_azi_std*2,'-or');
errorbar(periods,aniso_azi-180,aniso_azi_std*2,'-or');
ylim([50 180]);
xlim([20 150]);
ylabel('\phi');
xlabel('Periods (s)');

%% Plot residuals
clear residuals
for ip = 1:length(eventphv_ani)
    isgood = eventphv_ani(ip).isgood;
    dt_res = eventphv_ani(ip).dt_res(isgood);
    residuals(ip).rms_dt_res = rms(dt_res(:));
    residuals(ip).mean_dt_res = mean(dt_res(:));
    residuals(ip).dt_res = dt_res(:);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  WARNING! THIS PLOTS ALL DATA AT ONCE AND MAY REQUIRE A LOT OF MEMORY!
%     figure(87); clf; set(gcf,'color','w','position',[1035         155         560         781]);
%     for ip = 1:length(periods)
%         subplot(2,1,1);
%         plot(periods(ip),residuals(ip).mean_dt_res,'o','color',[0.7 0.7 0.7]); hold on;
%         plot(periods(ip),nanmean(residuals(ip).mean_dt_res),'rs','linewidth',2,'markersize',10);
%         ylabel('mean (dt_{obs}-dt_{pre})')
%         set(gca,'linewidth',1.5,'fontsize',15);
%         subplot(2,1,2);
%         plot(periods(ip),residuals(ip).rms_dt_res,'o','color',[0.7 0.7 0.7]); hold on;
%         plot(periods(ip),nanmean(residuals(ip).rms_dt_res),'rs','linewidth',2,'markersize',10);
%         xlabel('Period (s)');
%         ylabel('RMS (dt_{obs}-dt_{pre})')
%         set(gca,'linewidth',1.5,'fontsize',15);
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(86); clf;
for ip = 1:length(periods)
    subplot(M,N,ip)
    histogram(residuals(ip).dt_res);
    title([num2str(periods(ip)),' s'],'fontsize',15)
    xlabel('Residual (s)');
end
