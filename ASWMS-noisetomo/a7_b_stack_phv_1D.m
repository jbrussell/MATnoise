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

dc_thresh = 5; % [%] remove velocity perturbations larger than this
min_goodnum = 5; % minimum number of GSDF measurements
min_Mw = 5.0; % minimum magnitude
% max_evdp = 20; % [km] maximum event depth
min_nodes_resolved = 0; % minimum number of resolved nodes in final model

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
		azi = angle(eventphv(ip).GVx + eventphv(ip).GVy.*sqrt(-1));
		azi = rad2deg(azi);
		azi_mat(:,:,ie,ip) = azi;
        weights(ie,ip) = eventphv(ip).goodnum.^(-1/2); % jbr;
	end
	gcazi_mat(:,:,ie) = azimuth(xi,yi,evla(ie),evlo(ie));
end
weights(isinf(weights)) = 0;

% %%
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
			
			% Get rid of large outliers
			phV(abs((phV-nanmedian(phV))/nanmedian(phV))*100>dc_thresh) = nan;


            isophv(mi,mj)=nanmean(phV);
            % isophv_std(mi,mj)=parastd(2,1)-parastd(1,1);
			isophv_std(mi,mj)=nanstd(phV);
            aniso_strength(mi,mj)=0;
            aniso_azi(mi,mj)=0;			
            aniso_strength_std(mi,mj)=0;
            aniso_azi_std(mi,mj)=0;
            
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
%                 plot(azi,phV,'.');
                plot(azi,dv,'.');
                allazi = 0:360; %-200:200;
                plot(allazi,aniso_strength(mi,mj)*cosd(2*(allazi-aniso_azi(mi,mj)))*100,'r');
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

filename = [workingdir,'eikonal_stack_1D_',comp,'.mat'];
save(filename,'avgphv_aniso');

%%
% %%
%plot native
% qfile = './qfiles/pa5_5km.s0to66.q';
% mode0 = readMINEOS_qfile(qfile,0);
% mode1 = readMINEOS_qfile(qfile,1);
% mode2 = readMINEOS_qfile(qfile,2);
% mode3 = readMINEOS_qfile(qfile,3);

figure(58);
set(gcf,'position',[351   677   560   668]);
clf
clear avgv avgv_std aniso_str aniso_str_std aniso_azi aniso_azi_std
for ip = 1:length(periods)
    avgv(ip) = nanmean(avgphv_aniso(ip).isophv(:));
    avgv_std(ip) = nanmean(avgphv_aniso(ip).isophv_std(:));
    aniso_str(ip) = nanmean(avgphv_aniso(ip).aniso_strength(:));
    aniso_str_std(ip) = nanmean(avgphv_aniso(ip).aniso_strength_std(:));
    aniso_azi(ip) = nanmean(avgphv_aniso(ip).aniso_azi(:));
    aniso_azi_std(ip) = nanmean(avgphv_aniso(ip).aniso_azi_std(:));
end
%plot native
subplot(3,1,1); hold on;
% plot(mode0.T,mode0.phv,'-','color',[0.8,0.8,0.8],'linewidth',3);
% plot(mode1.T,mode1.phv,'-','color',[0.8,0.8,0.8],'linewidth',3);
% plot(mode2.T,mode2.phv,'-','color',[0.8,0.8,0.8],'linewidth',3);
% plot(mode3.T,mode3.phv,'-','color',[0.8,0.8,0.8],'linewidth',3);
errorbar(periods,avgv,avgv_std,'-or');
ylim([3.85 4.4]);
% xlim([20 150]);
ylabel('c (km/s)');
%plot native
subplot(3,1,2); hold on;
errorbar(periods,aniso_str*100*2,aniso_str_std*100,'-or');
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
ylim([50 140]);
% xlim([20 150]);
ylabel('\phi');
xlabel('Periods (s)');


