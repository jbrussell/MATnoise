%% This script is similar to stack_phv.m, but adding azimuthal anisotropy inversion into the result.
%
% written by Ge Jin, LDEO
% jinwar@gmail.com, ge.jin@ldeo.columbia.edu
%
% This version assumes a 1-D velocity structure (i.e., inputs measurements
% from a6_eikonal_eq_flat) and includes weighting based on number of good
% GSDF measurements that go into each phv estimate.
% JBR - 11/19



clear;
%plot native
% close 11
setup_parameters

% phase_v_path = './eikonal/'
workingdir = parameters.workingdir;
phase_v_path = [workingdir,'eikonal/'];
eikonl_propazi_output_path = [workingdir,'eikonal_propazi/'];

is_offgc_propagation = parameters.is_offgc_propagation; % Account for off-great-circle propagation using eikonal tomography maps? Otherwise will assume great-circle propagation.

dc_thresh = 5; % [%] remove velocity perturbations larger than this
min_goodnum = 5; % minimum number of GSDF measurements
min_Mw = 5.0; % minimum magnitude
% max_evdp = 20; % [km] maximum event depth
min_nodes_resolved = 0; % minimum number of resolved nodes in final model
min_nbin = 5; % minimum number of measurements to include bin
azi_bin_deg = parameters.azi_bin_deg_ani;

r = 0.05;
isfigure = 1;
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
% smsize = parameters.smsize;
off_azi_tol = parameters.off_azi_tol;
is_one_phi = parameters.is_one_phi;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);
smsize = Nx+Ny;

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

GV_mat = nan(Nx,Ny,length(phvmatfiles),length(periods));
azi_mat = nan(Nx,Ny,length(phvmatfiles),length(periods));
raydense_mat = nan(Nx,Ny,length(phvmatfiles),length(periods));
weights = nan(length(phvmatfiles),length(periods));
% Gather information
for ie = 1:length(phvmatfiles)
	temp = load([phase_v_path,phvmatfiles(ie).name]);
	eventphv = temp.eventphv;
	if is_offgc_propagation
        eikonl_propazi_output_path = [workingdir,'eikonal_propazi/'];
        temp = load([eikonl_propazi_output_path,'/',eventphv(1).id,'_eikonal_',comp,'.mat']);
        eventphv_propazi = temp.eventphv;
    end
	disp(eventphv(1).id);
	evla(ie) = eventphv(ip).evla;
	evlo(ie) = eventphv(ip).evlo;
	for ip=1:length(periods)
        if eventphv(ip).goodnum<min_goodnum || ...
            eventphv(ip).Mw<min_Mw 
                    % eventphv(ip).evdp>max_evdp
            continue;
        end
        GV_mat(:,:,ie,ip) = eventphv(ip).GV;
        raydense_mat(:,:,ie,ip) = eventphv(ip).raydense;
		if is_offgc_propagation
            azi = angle(eventphv_propazi(ip).GVx + eventphv_propazi(ip).GVy.*sqrt(-1));
        else
            azi = angle(eventphv(ip).GVx + eventphv(ip).GVy.*sqrt(-1));
        end
		azi = rad2deg(azi);
		azi_mat(:,:,ie,ip) = azi;
        weights(ie,ip) = eventphv(ip).goodnum.^(-1/2); % jbr;
	end
	gcazi_mat(:,:,ie) = azimuth(xi,yi,evla(ie),evlo(ie));
end
weights(isinf(weights)) = 0;

% %%
figure(11); clf;
N=4; M = floor(length(periods)/N)+1;
for ip = 1:length(periods)
    isophv=nan(Nx,Ny);
    isophv_std=nan(Nx,Ny);
    aniso_strength=nan(Nx,Ny);
    aniso_strength_std=nan(Nx,Ny);
    aniso_azi=nan(Nx,Ny);
    aniso_azi_std=nan(Nx,Ny);
    aniso_1phi_strength=nan(Nx,Ny);
    aniso_1phi_azi=nan(Nx,Ny);	
    phV = [];
    azi = [];
    gcazi = [];
	% start to inverse azimuthal anisotropy grid by grid
	for mi=1
		disp(['ip:',num2str(ip),' process: ',num2str(mi/Nx)]);
        for mj=1
            n=0;
            clear phV_best azi phV dist gcazi
			for ie = 1:eventnum
				diff_az = angdiff(azi_mat(:,:,ie,ip)*pi/180,gcazi_mat(:,:,ie)*pi/180)*180/pi;
                if nanmean(abs(diff_az(:))) > off_azi_tol
                    continue
                end
				avgV=GV_mat(mi,mj,ie,ip);
                % Make sure enough grid nodes are resolved
                I_resolv = ~isnan(GV_mat(:,:,ie,ip));
                N_resolved = sum(I_resolv(:));
                if N_resolved < min_nodes_resolved
                    continue;
                end
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
                            weight(n) = weights(ie,ip);
                        end
                    end
                end
			end  % enent loop

%             if n < min_event_num*((2*smsize).^2)
            % if n < min_event_num*((0.5*smsize).^2)
            %     isophv(mi,mj)=NaN;
            %     isophv_std(mi,mj)=NaN;
            %     aniso_strength(mi,mj)=NaN;
            %     aniso_strength_std(mi,mj)=NaN;
            %     aniso_azi(mi,mj)=NaN;
            %     aniso_azi_std(mi,mj)=NaN;
            %     continue;
            % end

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

%             if n < min_event_num*((2*smsize).^2)
            % if n < min_event_num*((0.5*smsize).^2)
            %     isophv(mi,mj)=NaN;
            %     isophv_std(mi,mj)=NaN;
            %     aniso_strength(mi,mj)=NaN;
            %     aniso_strength_std(mi,mj)=NaN;
            %     aniso_azi(mi,mj)=NaN;
            %     aniso_azi_std(mi,mj)=NaN;
            %     continue;
            % end
			
			% Get rid of large outliers
			phV(abs((phV-nanmedian(phV))/nanmedian(phV))*100>dc_thresh) = nan;

            if is_one_phi
                [para fiterr]=fit_azi_anisotropy_1phi(azi,phV);
            else
                [para fiterr]=fit_azi_anisotropy(azi,phV,weight);
            end
            
			parastd=confint(para,.95);
            isophv(mi,mj)=para.a;
            % isophv_std(mi,mj)=parastd(2,1)-parastd(1,1);
			isophv_std(mi,mj)=parastd(2,1)-para.a;
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
                aniso_strength_std(mi,mj)=parastd(2,4)-parastd(1,4);
                aniso_azi_std(mi,mj)=parastd(2,5)-parastd(1,5);
            else
                % aniso_strength_std(mi,mj)=parastd(2,2)-parastd(1,2);
                % aniso_azi_std(mi,mj)=parastd(2,3)-parastd(1,3);				
				aniso_strength_std(mi,mj)=parastd(2,2)-para.d;
                aniso_azi_std(mi,mj)=parastd(2,3)-para.e;
            end         
            
            if ~is_one_phi
                % Fit binned measurements
                azi(azi<0) = azi(azi<0)+360;
                w = ones(size(azi));
                [~,Isort] = sort(azi);
                azi_srt = azi(Isort);
                phv_srt = phV(Isort);
                w_srt = w(Isort); %phv_std(Isort); % weight by std
                bins = [0:azi_bin_deg:360];
                phv_bin = 0; azi_bin = 0; phv_bin_err = 0; azi_bin_err = 0;
                for ibin = 1:length(bins)-1
                    I_bin = azi_srt>=bins(ibin) & azi_srt<bins(ibin+1);
                    if sum(I_bin)<min_nbin
                        I_bin = false(size(I_bin));
                    end
                    wbin = w_srt(I_bin);
                    wbinnorm = wbin/sum(wbin);
                    % Weighted means and standard deviations
%                     phv_bin(ibin) = nansum(wbin .* phv_srt(I_bin)) / nansum(wbin);
%                     phv_bin_err(ibin) = sqrt( var(phv_srt(I_bin) , wbinnorm) );
                    phv_bin(ibin) = nanmean(phv_srt(I_bin));
                    phv_bin_err(ibin) = nanstd(phv_srt(I_bin));
                    azi_bin(ibin) = (bins(ibin)+bins(ibin+1))/2;
                    weight(ibin) = sum(I_bin);
                end

            %     [para fiterr]=fit_azi_anisotropy(fit_azi_bin.meas(ip).azi, fit_azi_bin.meas(ip).phv, fit_azi_bin.meas(ip).phv_std.^(1/2) .* sqrt(weight));
            %     [para fiterr]=fit_azi_anisotropy(fit_azi_bin.meas(ip).azi, fit_azi_bin.meas(ip).phv, sqrt(weight));
                if length(find(~isnan(phv_bin))) > 3
                    [para_bin fiterr]=fit_azi_anisotropy(azi_bin, phv_bin);
                    parastd=confint(para_bin,.95);
                    avgphv_aniso(ip).bin.phv = phv_bin;
                    avgphv_aniso(ip).bin.phv_std = phv_bin_err;
                    avgphv_aniso(ip).bin.azi = azi_bin;
                    avgphv_aniso(ip).bin.periods=periods(ip);
                    avgphv_aniso(ip).bin.c_iso=para_bin.a;
                    avgphv_aniso(ip).bin.c_iso_95=parastd(2,1)-para_bin.a;
                    avgphv_aniso(ip).bin.A2=para_bin.d;
                    avgphv_aniso(ip).bin.A2_95=parastd(2,2)-para_bin.d;
                    avgphv_aniso(ip).bin.phi2=para_bin.e;
                    avgphv_aniso(ip).bin.phi2_95=parastd(2,3)-para_bin.e;
                else
                    avgphv_aniso(ip).bin.periods=periods(ip);
                    avgphv_aniso(ip).bin.c_iso=nan;
                    avgphv_aniso(ip).bin.c_iso_95=nan;
                    avgphv_aniso(ip).bin.A2=nan;
                    avgphv_aniso(ip).bin.A2_95=nan;
                    avgphv_aniso(ip).bin.phi2=nan;
                    avgphv_aniso(ip).bin.phi2_95=nan;
                end
            end
            
            
            if is_one_phi && isfigure
                figure(11)
                clf
                hold on
                plot(azi,phV,'x');
                allazi = -200:200;
                plot(allazi,para.a*(1+para.b*cosd(allazi-para.c)+para.d*cosd(2*(allazi-para.e))),'r')
            elseif ~is_one_phi && isfigure
				%plot native 
                figure(11)
                subplot(M,N,ip);
%                 clf
                hold on
                azi(azi<0) = azi(azi<0)+360;
                dv = (phV-isophv(mi,mj))./isophv(mi,mj)*100;
                dv_bin = (avgphv_aniso(ip).bin.phv -isophv(mi,mj))./isophv(mi,mj)*100;
                dv_bin_std = (avgphv_aniso(ip).bin.phv_std)./isophv(mi,mj)*100;
%                 plot(azi,phV,'.');
                plot(azi,dv,'.');
                allazi = 0:360; %-200:200;
                plot(allazi,para.d*cosd(2*(allazi-para.e))*100,'r');
                errorbar(avgphv_aniso(ip).bin.azi,dv_bin,dv_bin_std,'ok');
                title([num2str(periods(ip)),' s']);
%                 plot(allazi,para.a*(1+para.d*cosd(2*(allazi-fastdir_plot))),'--k');
                ylim([-5 5]);
                xticks([0 90 180 270 360]);
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
    avgphv_aniso(ip).period = periods(ip);
	if is_one_phi
		avgphv_aniso(ip).aniso_1phi_strength=aniso_1phi_strength;
		avgphv_aniso(ip).aniso_1phi_azi=aniso_1phi_azi;
		avgphv_aniso(ip).aniso_strength_std=aniso_strength_std;
		avgphv_aniso(ip).aniso_azi_std=aniso_azi_std;
	else
		avgphv_aniso(ip).aniso_strength_std=aniso_strength_std;
		avgphv_aniso(ip).aniso_azi_std=aniso_azi_std;
	end          
end % end of period loop

filename = [workingdir,'eikonal_stack_aniso1D_',comp,'.mat'];
save(filename,'avgphv_aniso');

%%
% %%
%plot native
figure(58);
set(gcf,'position',[351   677   560   668]);
clf
clear avgv avgv_std aniso_str aniso_str_std aniso_azi aniso_azi_std aniso_azi_bin_2std aniso_azi_bin aniso_str_bin_2std aniso_str_bin avgv_bin_2std avgv_bin
for ip = 1:length(periods)
    avgv(ip) = nanmean(avgphv_aniso(ip).isophv(:));
    avgv_std(ip) = nanmean(avgphv_aniso(ip).isophv_std(:));
    aniso_str(ip) = nanmean(avgphv_aniso(ip).aniso_strength(:));
    aniso_str_std(ip) = nanmean(avgphv_aniso(ip).aniso_strength_std(:));
    aniso_azi(ip) = nanmean(avgphv_aniso(ip).aniso_azi(:));
    aniso_azi_std(ip) = nanmean(avgphv_aniso(ip).aniso_azi_std(:));
    
    avgv_bin(ip) = avgphv_aniso(ip).bin.c_iso;
    avgv_bin_std(ip) = avgphv_aniso(ip).bin.c_iso_95;
    aniso_str_bin(ip) = avgphv_aniso(ip).bin.A2;
    aniso_str_bin_std(ip) = avgphv_aniso(ip).bin.A2_95;
    aniso_azi_bin(ip) = avgphv_aniso(ip).bin.phi2;
    aniso_azi_bin_std(ip) = avgphv_aniso(ip).bin.phi2_95;
end
%plot native
subplot(3,1,1); hold on;
errorbar(periods,avgv,avgv_std,'-or');
errorbar(periods,avgv_bin,avgv_bin_std,'-ob');
ylim([3.85 4.4]);
% xlim([20 150]);
ylabel('c (km/s)');
%plot native

subplot(3,1,2); hold on;
errorbar(periods,aniso_str*100*2,aniso_str_std*100,'-or');
errorbar(periods,aniso_str_bin*100*2,aniso_str_bin_std*100,'-ob');

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
errorbar(periods,aniso_azi_bin,aniso_azi_bin_std,'-ob');
errorbar(periods,aniso_azi_bin+180,aniso_azi_bin_std,'-ob');
errorbar(periods,aniso_azi_bin-180,aniso_azi_bin_std,'-ob');
ylim([50 140]);
% xlim([20 150]);
ylabel('\phi');
xlabel('Periods (s)');


