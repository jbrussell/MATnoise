%% This script is similar to stack_phv.m, but adding azimuthal anisotropy inversion into the result.
%
% written by Ge Jin, LDEO
% jinwar@gmail.com, ge.jin@ldeo.columbia.edu
%


clear;
setup_parameters

workingdir = parameters.workingdir;
phase_v_path = [workingdir,'eikonal/'];
r = 0.05;
isfigure = 0;

APM = 114; % absolute plate motion (GSRM 2.1; NNR) https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html
FSD = 75; 

comp = parameters.component;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
min_phv_tol = parameters.min_phv_tol;
max_phv_tol = parameters.max_phv_tol;
is_raydense_weight = parameters.is_raydense_weight;
min_event_num = parameters.min_event_num;
err_std_tol = parameters.err_std_tol;
smsize = parameters.smsize;
off_azi_tol = parameters.off_azi_tol;
is_one_phi = parameters.is_one_phi;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

for ip=1:length(periods)
	avgphv(ip).sumV = zeros(Nx,Ny);
	avgphv(ip).sumweight = zeros(Nx,Ny);
	avgphv(ip).GV_std = zeros(Nx,Ny);
	avgphv(ip).eventnum = zeros(Nx,Ny);
	avgphv(ip).xi = xi;
	avgphv(ip).yi = yi;
	avgphv(ip).xnode = xnode;
	avgphv(ip).ynode = ynode;
	avgphv(ip).period = periods(ip);
end

phvmatfiles = dir([phase_v_path,'/*_',comp,'.mat']);
eventnum =  length(phvmatfiles);

GV_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));
azi_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));
raydense_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));

% Gather information
for ie = 1:length(phvmatfiles)
	temp = load([phase_v_path,phvmatfiles(ie).name]);
	eventphv = temp.eventphv;
	disp(eventphv(1).id);
	evla(ie) = eventphv(ip).evla;
	evlo(ie) = eventphv(ip).evlo;
	for ip=1:length(periods)
        GV_mat(:,:,ie,ip) = eventphv(ip).GV;
        raydense_mat(:,:,ie,ip) = eventphv(ip).raydense;
		azi = angle(eventphv(ip).GVx + eventphv(ip).GVy.*sqrt(-1));
		azi = rad2deg(azi);
		azi_mat(:,:,ie,ip) = azi;
	end
	gcazi_mat(:,:,ie) = azimuth(xi,yi,evla(ie),evlo(ie));
end

for ip = 1:length(periods)
    isophv=zeros(Nx,Ny);
    isophv_std=zeros(Nx,Ny);
    aniso_strength=zeros(Nx,Ny);
    aniso_strength_std=zeros(Nx,Ny);
    aniso_azi=zeros(Nx,Ny);
    aniso_azi_std=zeros(Nx,Ny);
    aniso_1phi_strength=zeros(Nx,Ny);
    aniso_1phi_azi=zeros(Nx,Ny);	
	% start to inverse azimuthal anisotropy grid by grid
	for mi=1:Nx
		disp(['ip:',num2str(ip),' process: ',num2str(mi/Nx)]);
        for mj=1:Ny
            n=0;
            clear phV_best azi phV dist gcazi
			for ie = 1:eventnum
				avgV=GV_mat(mi,mj,ie,ip);
                lowi=max(1,mi-smsize);
                upi=min(Nx,mi+smsize);
                lowj=max(1,mj-smsize);
                upj=min(Ny,mj+smsize);
                for ii=lowi:upi
                    for jj=lowj:upj
                        if ~isnan(GV_mat(ii,jj,ie,ip))
                            n=n+1;
                            azi(n)=azi_mat(ii,jj,ie,ip);
							gcazi(n) = gcazi_mat(ii,jj,ie);
                            phV(n)=GV_mat(ii,jj,ie,ip);
                        end
                    end
                end
			end  % enent loop

            if n < min_event_num*((2*smsize).^2)
                isophv(mi,mj)=NaN;
                isophv_std(mi,mj)=NaN;
                aniso_strength(mi,mj)=NaN;
                aniso_strength_std(mi,mj)=NaN;
                aniso_azi(mi,mj)=NaN;
                aniso_azi_std(mi,mj)=NaN;
                continue;
            end

			% Get rid of significant off circle events
			diffazi = azi - gcazi;
			ind = find(diffazi>180);
			if ~isempty(ind)
				diffazi(ind) = diffazi(ind) - 360;
			end
			ind = find(diffazi<-180);
			if ~isempty(ind)
				diffazi(ind) = diffazi(ind) + 360;
			end
			ind = find(abs(diffazi)> off_azi_tol);
			if ~isempty(ind)
				azi(ind) = [];
				phV(ind) = [];
			end

            if n < min_event_num*((2*smsize).^2)
                isophv(mi,mj)=NaN;
                isophv_std(mi,mj)=NaN;
                aniso_strength(mi,mj)=NaN;
                aniso_strength_std(mi,mj)=NaN;
                aniso_azi(mi,mj)=NaN;
                aniso_azi_std(mi,mj)=NaN;
                continue;
            end

            if is_one_phi
                [para fiterr]=fit_azi_anisotropy_1phi(azi,phV);
            else
                [para fiterr]=fit_azi_anisotropy(azi,phV);
            end
			parastd=confint(para);
            isophv(mi,mj)=para.a;
            isophv_std(mi,mj)=(parastd(2,1)-parastd(1,1))/2;
            aniso_strength(mi,mj)=para.d;
            aniso_azi(mi,mj)=para.e;
            if para.e > 180
                aniso_azi(mi,mj)=para.e-180;
            elseif para.e < 0
                aniso_azi(mi,mj)=para.e+180;
            end
            if is_one_phi
                aniso_1phi_strength(mi,mj)=para.b;
                aniso_1phi_azi(mi,mj)=para.c;
                aniso_strength_std(mi,mj)=(parastd(2,4)-parastd(1,4))/2;
                aniso_azi_std(mi,mj)=(parastd(2,5)-parastd(1,5))/2;
            else
                aniso_strength_std(mi,mj)=(parastd(2,2)-parastd(1,2))/2;
                aniso_azi_std(mi,mj)=(parastd(2,3)-parastd(1,3))/2;
            end          
            if is_one_phi && isfigure
                figure(11)
                clf
                hold on
                plot(azi,phV,'x');
                allazi = -200:200;
                plot(allazi,para.a*(1+para.b*cosd(allazi-para.c)+para.d*cosd(2*(allazi-para.e))),'r')
            end
		end  % mj loop
	end % mi loop
	avgphv_aniso(ip).isophv = isophv;
	avgphv_aniso(ip).isophv_std = isophv_std;
	avgphv_aniso(ip).aniso_strength = aniso_strength;
	avgphv_aniso(ip).aniso_azi = aniso_azi;
	avgphv_aniso(ip).parameters = parameters;
	avgphv_aniso(ip).xi = xi;
	avgphv_aniso(ip).yi = yi;
	if is_one_phi
		avgphv_aniso(ip).aniso_1phi_strength=aniso_1phi_strength;
		avgphv_aniso(ip).aniso_1phi_azi=aniso_1phi_azi;
		avgphv_aniso(ip).aniso_strength_std=aniso_strength_std;
		avgphv_aniso(ip).aniso_azi_std=aniso_azi_std;
	else
		avgphv_aniso(ip).aniso_strength_std=aniso_strength_std;
		avgphv_aniso(ip).aniso_azi_std=aniso_azi_std;
	end          
    
    % Get average values from 2-D maps
    avgphv_aniso(ip).isophv_2d_mean = nanmean(avgphv_aniso(ip).isophv(:));
    % Report which ever is largest, std of map or mean of stds
    avgphv_aniso(ip).isophv_2d_mean_err = max([nanstd(avgphv_aniso(ip).isophv(:)), nanmean(avgphv_aniso(ip).isophv_std(:))]);
    % Get values from center of array
    latc = mean(xi(~isnan(avgphv_aniso(ip).isophv))); % center latitude
    lonc = mean(yi(~isnan(avgphv_aniso(ip).isophv))); % center longitude
    [~,ilat] = min(abs(xnode-latc));
    [~,ilon] = min(abs(ynode-lonc));
    isophv_2d_block = avgphv_aniso(ip).isophv(ilat+[-1:1],ilon+[-1:1]);
    isophv_2d_block_std = avgphv_aniso(ip).isophv_std(ilat+[-1:1],ilon+[-1:1]);
    avgphv_aniso(ip).isophv_2d_center = nanmean(isophv_2d_block(:));
    avgphv_aniso(ip).isophv_2d_center_err = max([nanstd(isophv_2d_block(:)), nanmean(isophv_2d_block_std(:))]);
    avgphv_aniso(ip).latc = latc;
    avgphv_aniso(ip).lonc = lonc;
    
    [ A_mean, phi_mean ] = mean_aniso(avgphv_aniso(ip).aniso_strength(:), avgphv_aniso(ip).aniso_azi(:));
    avgphv_aniso(ip).aniso_strength_2d_mean = A_mean;
    avgphv_aniso(ip).aniso_azi_2d_mean = phi_mean;
    avgphv_aniso(ip).aniso_strength_2d_mean_err = max([nanstd(avgphv_aniso(ip).aniso_strength(:)), nanmean(avgphv_aniso(ip).aniso_strength_std(:))]);
    avgphv_aniso(ip).aniso_azi_2d_mean_err = nanmean(avgphv_aniso(ip).aniso_azi_std(:));
    aniso_strength_2d_block = avgphv_aniso(ip).aniso_strength(ilat+[-1:1],ilon+[-1:1]);
    aniso_strength_2d_block_std = avgphv_aniso(ip).aniso_strength_std(ilat+[-1:1],ilon+[-1:1]);
    aniso_azi_2d_block = avgphv_aniso(ip).aniso_azi(ilat+[-1:1],ilon+[-1:1]);
    aniso_azi_2d_block_std = avgphv_aniso(ip).aniso_azi_std(ilat+[-1:1],ilon+[-1:1]);
    [ A_mean, phi_mean ] = mean_aniso(aniso_strength_2d_block(:), aniso_azi_2d_block(:));
    avgphv_aniso(ip).aniso_strength_2d_center = A_mean;
    avgphv_aniso(ip).aniso_strength_2d_center_err = max([nanstd(aniso_strength_2d_block(:)), nanmean(aniso_strength_2d_block_std(:))]);
    avgphv_aniso(ip).aniso_azi_2d_center = phi_mean;
    avgphv_aniso(ip).aniso_azi_2d_center_err = nanmean(aniso_azi_2d_block_std(:));
    
end % end of period loop

filename = [workingdir,'eikonal_stack_aniso_',comp,'.mat'];
save(filename,'avgphv_aniso');

%%
N=3; M = floor(length(periods)/N)+1;
figure(56)
clf
for ip = 1:length(periods)
	subplot(M,N,ip);
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv_aniso(ip).isophv);
	colorbar
	load seiscmap
	colormap(seiscmap)
	drawnow
	avgv = nanmean(avgphv_aniso(ip).isophv(:));
	caxis([avgv*(1-r) avgv*(1+r)])
    title([num2str(periods(ip)),' s'])
end

%%

figure(57)
clf
for ip = 1:length(periods)
	subplot(M,N,ip);
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
% 	h1=surfacem(xi,yi,avgphv_aniso(ip).isophv);
	h1=surfacem(xi,yi,avgphv_aniso(ip).aniso_strength*100);
	colorbar
	load seiscmap
	colormap(parula)
	drawnow
 	avgv = nanmean(avgphv_aniso(ip).isophv(:));
% 	caxis([avgv*(1-r) avgv*(1+r)])
    caxis([0 5]);
    
    scale = 50;
%     u=eventphv_ani(ip).A2.*cosd(eventphv_ani(ip).phi2)*scale;
% 	v=eventphv_ani(ip).A2.*sind(eventphv_ani(ip).phi2)*scale;%./cosd(mean(lalim));
    u=avgphv_aniso(ip).aniso_strength.*cosd(avgphv_aniso(ip).aniso_azi)*scale;
	v=avgphv_aniso(ip).aniso_strength.*sind(avgphv_aniso(ip).aniso_azi)*scale;%./cosd(mean(lalim));
	[m n]=size(xi);
    hold on;
    xpts = [];
    ypts = [];
    for ix=1:m
        for iy=1:n
            xpts = [xpts, [xi(ix,iy)-u(ix,iy)/2 xi(ix,iy)+u(ix,iy)/2]+gridsize/2, nan];
            ypts = [ypts, [yi(ix,iy)-v(ix,iy)/2 yi(ix,iy)+v(ix,iy)/2]+gridsize/2, nan];
        end
    end
%     plotm(xpts_bg,ypts_bg,'-','Color',[0 0 0],'linewidth',4);
    plotm(xpts,ypts,'-','Color',[0.9 0 0],'linewidth',2);
    hold on;
    title([num2str(periods(ip)),' s'])
    % Plot reference
%     refstick = scale*0.02;
%     plotm([min(lalim) min(lalim)]+abs(diff(lalim))*0.15,[max(lolim)-refstick/2 max(lolim)+refstick/2]-abs(diff(lolim))*0.15,'-','Color',[0.9 0 0],'linewidth',2);
%     textm(min(lalim)+abs(diff(lalim))*0.09,max(lolim)-abs(diff(lolim))*0.15,'2%','fontsize',12,'HorizontalAlignment', 'center');
end

%% Plot 1d estimates
figure(59);
set(gcf,'position',[351   677   560   668]);
clf
clear avgv avgv_std aniso_str aniso_str_std aniso_azi aniso_azi_std aniso_azi_bin_2std aniso_azi_bin aniso_str_bin_2std aniso_str_bin avgv_bin_2std avgv_bin
for ip = 1:length(periods)
    avgv(ip) = avgphv_aniso(ip).isophv_2d_mean;
    avgv_std(ip) = avgphv_aniso(ip).isophv_2d_mean_err;
    aniso_str(ip) = avgphv_aniso(ip).aniso_strength_2d_mean;
    aniso_str_std(ip) = avgphv_aniso(ip).aniso_strength_2d_mean_err;
    aniso_azi(ip) = avgphv_aniso(ip).aniso_azi_2d_mean;
    aniso_azi_std(ip) = avgphv_aniso(ip).aniso_azi_2d_mean_err;
    
    avgv_center(ip) = avgphv_aniso(ip).isophv_2d_center;
    avgv_center_std(ip) = avgphv_aniso(ip).isophv_2d_center_err;
    aniso_str_center(ip) = avgphv_aniso(ip).aniso_strength_2d_center;
    aniso_str_center_std(ip) = avgphv_aniso(ip).aniso_strength_2d_center_err;
    aniso_azi_center(ip) = avgphv_aniso(ip).aniso_azi_2d_center;
    aniso_azi_center_std(ip) = avgphv_aniso(ip).aniso_azi_2d_center_err;
    
end
%plot native
subplot(3,1,1); hold on;
errorbar(periods,avgv,avgv_std,'-or');
errorbar(periods,avgv_center,avgv_center_std,'-ob');
ylim([3.85 4.4]);
% xlim([20 150]);
ylabel('c (km/s)');
legend({'2-D Average','2-D Center'},'location','northwest');
%plot native

subplot(3,1,2); hold on;
errorbar(periods,aniso_str*100*2,aniso_str_std*100,'-or');
errorbar(periods,aniso_str_center*100*2,aniso_str_center_std*100,'-ob');

ylim([0 5]);
% xlim([20 150]);
ylabel('2A');
%plot native
subplot(3,1,3); hold on;
plot([periods(1),periods(end)],FSD*[1 1],'--k','linewidth',1.5); hold on;
plot([periods(1),periods(end)],APM*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1.5);
errorbar(periods,aniso_azi,aniso_azi_std,'-or');
errorbar(periods,aniso_azi+180,aniso_azi_std,'-or');
errorbar(periods,aniso_azi-180,aniso_azi_std,'-or');
errorbar(periods,aniso_azi_center,aniso_azi_center_std,'-ob');
errorbar(periods,aniso_azi_center+180,aniso_azi_center_std,'-ob');
errorbar(periods,aniso_azi_center-180,aniso_azi_center_std,'-ob');
ylim([50 140]);
% xlim([20 150]);
ylabel('\phi');
xlabel('Periods (s)');
