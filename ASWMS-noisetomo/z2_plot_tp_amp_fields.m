% Estimate travel time surface and gradient as well as amplitude gradient
% for use in eq (4) of Bao et al. (2016) GJI
% github.com/jbrussell
% 2021-05

clear;

isoverwrite = 1;
isfigure = 1;
is_save_amp_fig = 1;

min_Mw = 5.5; % minimum magnitude
min_Ngrcells = 20; % minimum numbe of grid cells required in order to use event
azi_bin_deg = 30; % [deg] size of azimuthal bins
min_nbin = 10; % minimum number of measurements in order to include bin

% setup parameters
setup_parameters

is_eikonal_ampgrad_norm = parameters.is_eikonal_ampgrad_norm;

r = 0.05;

workingdir = parameters.workingdir;
eventcs_path = [workingdir,'CSmeasure/'];
phase_v_path = [workingdir,'eikonal/'];
helmholtz_path = [workingdir,'helmholtz/'];
helmholtz_stack_file = [workingdir,'helmholtz_stack_',parameters.component];
traveltime_path = [workingdir,'traveltime/'];

if ~exist(traveltime_path,'dir')
	mkdir(traveltime_path);
end

% load stacked phase velocity map
load(helmholtz_stack_file);

% set up useful variables
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
xnode = lalim(1):gridsize:lalim(2);
ynode = lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
alpha_range = parameters.alpha_range;
alpha_search_grid = parameters.alpha_search_grid;
periods = parameters.periods;

eventfiles = dir([traveltime_path,'/*_traveltime_',parameters.component,'.mat']);

load seiscmap

if exist('badampsta.lst','file')
	badstnms = textread('badampsta.lst','%s');
	disp('Found Bad amplitude stations:')
	for ista = 1:length(badstnms)
	disp(badstnms(ista))
	end
end

clear attenuation
evcnt = 0;
for ie = 1:length(eventfiles)
    clear amp_term azi
    %for ie = 59
    % read in data for this event
    clear eventphv eventcs traveltime
    load(fullfile(traveltime_path,eventfiles(ie).name));
    eventid = traveltime(1).id;
    disp(eventid);    
    eventcsfile = [eventcs_path,'/',eventid,'_cs_',parameters.component,'.mat'];
    if exist(eventcsfile,'file')
        load(eventcsfile);
    else
        disp(['Cannot find CS file for ',eventid,', Skipped']);
        continue;
    end
    eventeikonalfile = [phase_v_path,'/',eventid,'_eikonal_',parameters.component,'.mat'];
    if exist(eventeikonalfile,'file')
        load(eventeikonalfile);
    else
        disp(['Cannot find eikonal file for ',eventid,', Skipped']);
        continue;
    end
    helmholtzfile = [helmholtz_path,'/',eventid,'_helmholtz_',parameters.component,'.mat'];
    if exist(helmholtzfile,'file')
        temp = load(helmholtzfile);
        helmholtz = temp.helmholtz;
    else
        disp(['Cannot find Helmholtz file for ',eventid,', Skipped']);
        continue;
    end

    if traveltime(1).Mw < min_Mw
        continue
    end
        
%         % Remove grid cells with outlier propagation azimuth
%         max_degrelmean = 20;
%         Inan = isnan(traveltime(ip).GV_cor);
%         azi_prop = 90 - atan2d(traveltime(ip).tp_gradlat',traveltime(ip).tp_gradlon');
%         azi_prop(azi_prop<0) = azi_prop(azi_prop<0)+360;
%         azi_prop(Inan) = nan;
% %         Ibad_prop = abs(azi_prop-nanmean(azi_prop(:))) > max_degrelmean;
% %         azi_prop(Ibad_prop) = nan;      
%         azi_GV = angle(eventphv(ip).GVx + eventphv(ip).GVy.*sqrt(-1));
% 		azi_GV = rad2deg(azi_GV) + 180;
%         absDiffDeg = @(a,b) abs(diff(unwrap([a,b]/180*pi)*180/pi));
%         for ii=1:length(xnode)
%             for jj=1:length(ynode)
%                 diffdeg(ii,jj) = absDiffDeg(azi_GV(ii,jj),azi_prop(ii,jj));
%             end
%         end
%         Ibad_prop = diffdeg > max_degrelmean;
%         avgphv(ip).GV_cor(Ibad_prop) = nan;
%         traveltime(ip).GV_cor(Ibad_prop) = nan;
        
        evcnt = evcnt+1;
      
    for ip = 1:length(avgphv)
        Inotnan = find(~isnan(traveltime(ip).GV_cor));
%         Inotnan = find(~isnan(azi_prop));
        if length(Inotnan) < min_Ngrcells
            mat(evcnt).tp(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).tp_ang(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).tp_grad(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).tp_gradlat(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).tp_gradlon(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).tp_laplat(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).tp_laplon(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).tp_lap(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).amp(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).amp_grad(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).amp_gradlat(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).amp_gradlon(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            if is_eikonal_ampgrad_norm
                mat(evcnt).amp_grad_ampnorm(:,:,ip) = nan(size(traveltime(ip).GV_cor));
                mat(evcnt).amp_gradlat_ampnorm(:,:,ip) = nan(size(traveltime(ip).GV_cor));
                mat(evcnt).amp_gradlon_ampnorm(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            end
            mat(evcnt).amp_laplat(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).amp_laplon(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).amp_lap(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).tp_focus(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).amp_decay(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).corr_amp_decay(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).amp_term(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).azi(:,:,ip) = nan(size(traveltime(ip).GV_cor));
            mat(evcnt).id = helmholtz(1).id;
            mat(evcnt).Mw = helmholtz(1).Mw;
            continue
        end
        
		%% fit the amplitude surface
		% reset the arrays
		clear stlas stlos tp
		stlas = eventcs.stlas;
		stlos = eventcs.stlos;
		stnms = eventcs.stnms;
        
        % Load amplitude fields
        amp = helmholtz(ip).ampmap;
        amp_grad = helmholtz(ip).amp_grad;
        amp_gradlat = helmholtz(ip).amp_gradlat;
        amp_gradlon = helmholtz(ip).amp_gradlon;
        [~,amp_laplat,~]=delm(xi,yi,amp_gradlat);
        [~,~,amp_laplon]=delm(xi,yi,amp_gradlon);
        amp_grad_ampnorm = helmholtz(ip).amp_grad_ampnorm;
        amp_gradlat_ampnorm = helmholtz(ip).amp_gradlat_ampnorm;
        amp_gradlon_ampnorm = helmholtz(ip).amp_gradlon_ampnorm;
        
        % Load travel-time fields
        tp_grad = traveltime(ip).tp_grad;
        tp_gradlat = traveltime(ip).tp_gradlat;
        tp_gradlon = traveltime(ip).tp_gradlon;
        tp_lap = traveltime(ip).tp_lap;
        [~,tp_laplat,~]=delm(xi,yi,tp_gradlat);
        [~,~,tp_laplon]=delm(xi,yi,tp_gradlon);
        
        % Get propagation azimuth at each grid cell
        azi_prop = traveltime(ip).tp_ang;
        azi_prop(azi_prop<0) = azi_prop(azi_prop<0)+360;
        
        % Get structural phase velocity
        phv = avgphv(ip).GV_cor ;
        phv(isnan(traveltime(ip).GV_cor)) = nan;

%         amp(Inan)=nan; amp_grad(Inan)=nan; amp_gradlat(Inan)=nan; amp_gradlon(Inan)=nan;
%         tp_lap(Inan)=nan; tp_grad(Inan)=nan; tp_gradlat(Inan)=nan; tp_gradlon(Inan)=nan;
        
%         figure(1);
%         subplot(1,2,1);
%         tp_grad(isnan(traveltime(ip).GV)) = nan;
%         imagesc(tp_grad);
%         cb = colorbar;
%         subplot(1,2,2);
%         imagesc(1./traveltime(ip).GV);
%         colorbar;
%         caxis(cb.Limits);

        % Calculate terms from Bao et al. (2016) equation 4
        % Amplitude decay term
        if is_eikonal_ampgrad_norm
            amp_decay = 2*(amp_gradlat_ampnorm.*tp_gradlat + amp_gradlon_ampnorm.*tp_gradlon);
        else
            amp_decay = 2*(amp_gradlat.*tp_gradlat + amp_gradlon.*tp_gradlon) ./ amp;
        end
        % Focusing correction term
        tp_focus = tp_lap;
        
        % smooth the terms
        if length(find(~isnan(amp_decay(:)))) < 5 || length(find(~isnan(tp_focus(:)))) < 5
            amp_decay = nan(size(amp_decay));
            tp_focus = nan(size(tp_focus));
        else            
            for ii=1
                smD=max([300 periods(ip).*parameters.refv]);
                amp_decay = gridfit_jg_geo(xi(:),yi(:),amp_decay(:),xnode,ynode,...
                                    'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal')';
                tp_focus = gridfit_jg_geo(xi(:),yi(:),tp_focus(:),xnode,ynode,...
                                    'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal')';
            end
        end
        Inan = isnan(phv);
        amp_decay(Inan) = nan;
        tp_focus(Inan) = nan;
        
        % Corrected amplitude decay
        corr_amp_decay = amp_decay + tp_focus;
        
        amp_term(:,:,ip) = (phv/2) .* corr_amp_decay;        
        azi(:,:,ip) = azi_prop;
        
        % Remove large outliers
        Ibad = amp_term>abs(nanmedian(amp_term(:)))*100 | amp_term<abs(nanmedian(amp_term(:)))*-100;
        amp_term(Ibad) = nan;
        
        nanmat = nan(size(phv));
        nanmat(~isnan(phv)) = 1;
        
        mat(evcnt).tp(:,:,ip) = traveltime(ip).tpmap .* nanmat;
        mat(evcnt).tp_ang(:,:,ip) = traveltime(ip).tp_ang .* nanmat;
        mat(evcnt).tp_grad(:,:,ip) = tp_grad .* nanmat;
        mat(evcnt).tp_gradlat(:,:,ip) = tp_gradlat .* nanmat;
        mat(evcnt).tp_gradlon(:,:,ip) = tp_gradlon .* nanmat;
        mat(evcnt).tp_laplat(:,:,ip) = tp_laplat .* nanmat;
        mat(evcnt).tp_laplon(:,:,ip) = tp_laplon .* nanmat;
        mat(evcnt).tp_lap(:,:,ip) = tp_lap .* nanmat;
        mat(evcnt).amp(:,:,ip) = amp .* nanmat;
        mat(evcnt).amp_grad(:,:,ip) = amp_grad .* nanmat;
        mat(evcnt).amp_gradlat(:,:,ip) = amp_gradlat .* nanmat;
        mat(evcnt).amp_gradlon(:,:,ip) = amp_gradlon .* nanmat;
        if is_eikonal_ampgrad_norm
            mat(evcnt).amp_grad_ampnorm(:,:,ip) = amp_grad_ampnorm .* nanmat;
            mat(evcnt).amp_gradlat_ampnorm(:,:,ip) = amp_gradlat_ampnorm .* nanmat;
            mat(evcnt).amp_gradlon_ampnorm(:,:,ip) = amp_gradlon_ampnorm .* nanmat;
        end
        mat(evcnt).amp_laplat(:,:,ip) = amp_laplat .* nanmat;
        mat(evcnt).amp_laplon(:,:,ip) = amp_laplon .* nanmat;
        mat(evcnt).amp_lap(:,:,ip) = helmholtz(ip).amp_lap .* nanmat;
        mat(evcnt).tp_focus(:,:,ip) = tp_focus .* nanmat;
        mat(evcnt).amp_decay(:,:,ip) = amp_decay .* nanmat;
        mat(evcnt).corr_amp_decay(:,:,ip) = corr_amp_decay .* nanmat;
        mat(evcnt).amp_term(:,:,ip) = amp_term(:,:,ip) .* nanmat;
        mat(evcnt).azi(:,:,ip) = azi(:,:,ip) .* nanmat;
        mat(evcnt).id = helmholtz(1).id;
        mat(evcnt).Mw = helmholtz(1).Mw;
        
%         if strcmp(eventid,'201112161254')
%             eventid
%             figure(63); clf; set(gcf,'position',[146+75          1         726        1024],'color','w');
%             N=3; M = floor(length(periods)/N)+1;
%             sgtitle(['Travel Time',' : ',mat(evcnt).id,' Mw',num2str(mat(evcnt).Mw)],'fontsize',18,'fontweight','bold');
%             for iper = 1:length(periods)    
%                 subplot(M,N,iper)
%                 ax = worldmap(lalim, lolim);
%                 surfacem(xi,yi,traveltime(iper).tpmap); hold on;
%                 quiverm(xi,yi,traveltime(iper).tp_gradlat,traveltime(ip).tp_gradlon,'-k')
%                 title([num2str(periods(iper)),' s'],'fontsize',15)
%                 cb = colorbar;
%         %         caxis([-1e-4 4e-4]);
%                 colormap(flip(seiscmap))
%             end
%         end
        
    end
    
    %% Do curve fitting (eq 9 in Bao et al. 2016)
%     
%     % Unbinned 2-D sinusoidal fit
%     for ix = 1:length(xnode)
%         for iy = 1:length(ynode)
%             amps = squeeze(amp_term(ix,iy,:));
%             azis = squeeze(azi(ix,iy,:));
%             if length(find(~isnan(amps)))>10
%                 [para, alpha, dlnbeta_dx, dlnbeta_dy]=fit_alpha_beta(azis,amps);
%                 parastd=confint(para,.95);
%                 dlnbeta_dx_err = (parastd(2,1)-parastd(1,1))/2;
%                 dlnbeta_dy_err = (parastd(2,2)-parastd(1,2))/2;
%                 alpha_err = (parastd(2,3)-parastd(1,3))/2;
%             else
%                 alpha = nan;
%                 dlnbeta_dx = nan;
%                 dlnbeta_dy = nan;
%                 alpha_err = nan;
%                 dlnbeta_dx_err = nan;
%                 dlnbeta_dy_err = nan;
%                 para = [];
%             end            
%             attenuation(ip).alpha_2d(ix,iy) = alpha;
%             attenuation(ip).dlnbeta_dx_2d(ix,iy) = dlnbeta_dx;
%             attenuation(ip).dlnbeta_dy_2d(ix,iy) = dlnbeta_dy;
%             attenuation(ip).alpha_2d_err(ix,iy) = alpha_err;
%             attenuation(ip).dlnbeta_dx_2d_err(ix,iy) = dlnbeta_dx_err;
%             attenuation(ip).dlnbeta_dy_2d_err(ix,iy) = dlnbeta_dy_err;
%             attenuation(ip).para_2d{ix,iy} = para;
%             attenuation(ip).amp_term_2d = amp_term;
%             attenuation(ip).azi = azi;
%             attenuation(ip).period = periods(ip);
%             attenuation(ip).evcnt = evcnt;
%         end
%     end
    
end

%% Plot 1D average alpha

for iev = 1:length(mat)
    
%     if strcmp(mat(iev).id,'201112161254')
%         mat(iev).id
%     else
%         continue
%     end

    %% Plot maps of Amplitude
    figure(60); clf; set(gcf,'position',[146           1         726        1024],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle(['Amplitude',' : ',mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold');
    for ip = 1:length(periods)    
        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp(:,:,ip)); hold on;
        quiverm(xi,yi,mat(iev).amp_gradlat(:,:,ip),mat(iev).amp_gradlon(:,:,ip),'-k')
        title([num2str(periods(ip)),' s'],'fontsize',15)
        cb = colorbar;
%         caxis([-1e-4 4e-4]);
        colormap(flip(seiscmap))
    end
    
    %% Plot maps of grad(Amplitude)
    figure(61); clf; set(gcf,'position',[146+25           1         726        1024],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle(['\nabla Amplitude',' : ',mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold');
    for ip = 1:length(periods)    
        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_grad(:,:,ip));
        title([num2str(periods(ip)),' s'],'fontsize',15)
        cb = colorbar;
%         caxis([-1e-4 4e-4]);
        colormap(flip(seiscmap))
    end
    
    %% Plot maps of lap(Amplitude)
    figure(62); clf; set(gcf,'position',[146+50           1         726        1024],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle(['\nabla^2 Amplitude',' : ',mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold');
    for ip = 1:length(periods)    
        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_lap(:,:,ip));
        title([num2str(periods(ip)),' s'],'fontsize',15)
        cb = colorbar;
%         caxis([-1e-4 4e-4]);
        colormap(flip(seiscmap))
    end
    
    %% Plot maps of traveltime
    figure(63); clf; set(gcf,'position',[146+75          1         726        1024],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle(['Travel Time',' : ',mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold');
    for ip = 1:length(periods)    
        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).tp(:,:,ip)); hold on;
        quiverm(xi,yi,mat(iev).tp_gradlat(:,:,ip),mat(iev).tp_gradlon(:,:,ip),'-k')
        title([num2str(periods(ip)),' s'],'fontsize',15)
        cb = colorbar;
%         caxis([-1e-4 4e-4]);
        colormap(flip(seiscmap))
    end
    
    %% Plot maps of grad(traveltime)
    figure(64); clf; set(gcf,'position',[146+100           1         726        1024],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle(['\nabla Travel Time',' : ',mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold');
    for ip = 1:length(periods)    
        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).tp_grad(:,:,ip));
        title([num2str(periods(ip)),' s'],'fontsize',15)
        cb = colorbar;
%         caxis([-1e-4 4e-4]);
        colormap(flip(seiscmap))
    end
    
    %% Plot maps of lap(traveltime)
    figure(65); clf; set(gcf,'position',[146+125           1         726        1024],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle(['\nabla^2 Travel Time',' : ',mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold');
    for ip = 1:length(periods)    
        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).tp_lap(:,:,ip));
        title([num2str(periods(ip)),' s'],'fontsize',15)
        cb = colorbar;
%         caxis([-1e-4 4e-4]);
        colormap(flip(seiscmap))
    end
    
    %% Plot maps of Amp Decay
    figure(66); clf; set(gcf,'position',[146+150           1         726        1024],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle(['$(2\nabla{A} \cdot \nabla{\tau})/A$',' : ',mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold','interpreter','latex');
    for ip = 1:length(periods)    
        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_decay(:,:,ip));
        title([num2str(periods(ip)),' s'],'fontsize',15)
        cb = colorbar;
%         caxis([-1e-4 4e-4]);
        colormap(flip(seiscmap))
    end
    
    %% Plot maps of Amp Term
    figure(67); clf; set(gcf,'position',[146+175           1         726        1024],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle(['r.h.s. of eq (9)',' : ',mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold');
    for ip = 1:length(periods)    
        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_term(:,:,ip));
        title([num2str(periods(ip)),' s'],'fontsize',15)
        cb = colorbar;
%         caxis([-1e-4 4e-4]);
        colormap(flip(seiscmap))
    end
    %%
    figdir = [workingdir,'/figs/fields/'];
    if ~exist(figdir)
        mkdir(figdir);
    end
    save2pdf([figdir,mat(iev).id,'_amp0_',parameters.component,'.pdf'],60,100);
    save2pdf([figdir,mat(iev).id,'_amp1_',parameters.component,'.pdf'],61,100);
    save2pdf([figdir,mat(iev).id,'_amp2_',parameters.component,'.pdf'],62,100);
    save2pdf([figdir,mat(iev).id,'_tt0_',parameters.component,'.pdf'],63,100);
    save2pdf([figdir,mat(iev).id,'_tt1_',parameters.component,'.pdf'],64,100);
    save2pdf([figdir,mat(iev).id,'_tt2_',parameters.component,'.pdf'],65,100);
    save2pdf([figdir,mat(iev).id,'_ampdecay_',parameters.component,'.pdf'],66,100);
    save2pdf([figdir,mat(iev).id,'_ampterm_',parameters.component,'.pdf'],67,100);
    
    %% Plot maps of grad(traveltime)
    figure(68); clf; set(gcf,'position',[246         551        1013         474],'color','w');
    sgtitle([mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold');
    for ip = 1 %1:length(periods)    
        subplot(2,3,1)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).tp_gradlon(:,:,ip));
        title('\nabla\tau_x','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,2)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).tp_gradlat(:,:,ip));
        title('\nabla\tau_y','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,3)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).tp_grad(:,:,ip));
        quiverm(xi,yi,mat(iev).tp_gradlat(:,:,ip),mat(iev).tp_gradlon(:,:,ip),'-k')
        title('|\nabla\tau|','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,4)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).tp_laplon(:,:,ip));
        title('\partial\nabla\tau_x / \partial{x}','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,5)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).tp_laplat(:,:,ip));
        title('\partial\nabla\tau_y / \partial{y}','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,6)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).tp_lap(:,:,ip));
        title('\nabla^2\tau','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
    end
    %% Plot maps of grad(amp)
    figure(69); clf; set(gcf,'position',[246         551        1013         474],'color','w');
    sgtitle([mat(iev).id,' Mw',num2str(mat(iev).Mw)],'fontsize',18,'fontweight','bold');
    for ip = 1 %1:length(periods)    
        subplot(2,3,1)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_gradlon(:,:,ip));
        title('\nabla{A}_x','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,2)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_gradlat(:,:,ip));
        title('\nabla{A}_y','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,3)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_grad(:,:,ip));
        quiverm(xi,yi,mat(iev).amp_gradlat(:,:,ip),mat(iev).amp_gradlon(:,:,ip),'-k')
        title('|\nabla{A}|','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,4)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_laplon(:,:,ip));
        title('\partial\nabla{A}_x / \partial{x}','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,5)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_laplat(:,:,ip));
        title('\partial\nabla{A}_y / \partial{y}','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
        
        subplot(2,3,6)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,mat(iev).amp_lap(:,:,ip));
        title('\nabla^2A','fontsize',15)
        cb = colorbar;
        colormap(flip(seiscmap))
    end
    
    %%
    save2pdf([figdir,mat(iev).id,'_tt3example_',parameters.component,'.pdf'],68,100);
    save2pdf([figdir,mat(iev).id,'_amp3example_',parameters.component,'.pdf'],69,100);
end

%% Save
% matfilename = fullfile(traveltime_path,[eventphv(1).id,'_attenuation_',parameters.component,'.mat']);
% save(matfilename,'attenuation');
% fprintf('\n');
% disp(['Saved to ',matfilename]);
% if is_save_amp_fig
%     figdir = [workingdir,'/figs/attenuation/'];
%     if ~exist(figdir)
%         mkdir(figdir);
%     end
%     save2pdf([figdir,eventphv(1).id,'_attenuation_',parameters.component,'.pdf'],39,100);
%     save2pdf([figdir,eventphv(1).id,'_attenuation_',parameters.component,'_StaAmps.pdf'],40,100);
% end