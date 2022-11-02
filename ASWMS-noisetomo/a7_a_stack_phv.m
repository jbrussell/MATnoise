% Program to stack phase velocity maps from each event

clear;

isfigure = 1;

setup_parameters

workingdir = parameters.workingdir;
% phase_v_path = './eikonal/'
phase_v_path = [workingdir,'eikonal/'];

r = 0.10;
load seiscmap

comp = parameters.component;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
min_csgoodratio = parameters.min_csgoodratio;
min_phv_tol = parameters.min_phv_tol;
max_phv_tol = parameters.max_phv_tol;
is_raydense_weight = parameters.is_raydense_weight;
err_std_tol = parameters.err_std_tol;
min_event_num = parameters.min_event_num;
issmoothmap = parameters.issmoothmap;
smooth_wavelength = parameters.smooth_wavelength;
event_bias_tol = parameters.event_bias_tol;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

phvmatfiles = dir([phase_v_path,'/*_eikonal_',comp,'.mat']);

GV_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));
raydense_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));

for ie = 1:length(phvmatfiles)
	temp = load([phase_v_path,phvmatfiles(ie).name]);
	eventphv = temp.eventphv;
	event_ids(ie) = {eventphv(1).id};
			
	disp(eventphv(1).id);
	for ip=1:length(periods)
        ind = find(eventphv(ip).GV < min_phv_tol);
        eventphv(ip).GV(ind) = NaN;
        ind = find(eventphv(ip).GV > max_phv_tol);
        eventphv(ip).GV(ind) = NaN;
		if eventphv(ip).goodnum./eventphv(ip).badnum < min_csgoodratio(ip)
			disp('not enough good cs measurement');
			eventphv(ip).GV(:) = NaN;
        end
		
        GV_mat(:,:,ie,ip) = eventphv(ip).GV;

        raydense_mat(:,:,ie,ip) = eventphv(ip).raydense;
	end
end

avgphv = average_GV_mat(GV_mat, raydense_mat, parameters);

% Calculate std, remove the outliers
GV_mat = 1./GV_mat;
for ip=1:length(periods)
	for i = 1:Nx
		for j=1:Ny
			avgphv(ip).GV_std(i,j) = nanstd(GV_mat(i,j,:,ip));
			ind = find( abs(GV_mat(i,j,:,ip) - 1./avgphv(ip).GV(i,j)) > err_std_tol*avgphv(ip).GV_std(i,j));
			GV_mat(i,j,ind,ip) = NaN;
		end
	end
end
GV_mat = 1./GV_mat;

% calculate the averaged phase velocity again
avgphv = average_GV_mat(GV_mat, raydense_mat, parameters);

% remove bias events
for ip=1:length(periods)
avg_GV = avgphv(ip).GV;
mean_phv = nanmean(avg_GV(:));
badnum = 0;
for ie=1:length(event_ids)
	GV = GV_mat(:,:,ie,ip);
	diff_phv = GV-avg_GV;
	diff_percent = nanmean(diff_phv(:))/mean_phv*100;
	if abs(diff_percent) > event_bias_tol;
% 		matfile = dir(fullfile('eikonal',[char(event_ids(ie)),'*.mat']));
% 		load(fullfile('eikonal',matfile(1).name));
        matfile = dir(fullfile(workingdir,'eikonal',[char(event_ids(ie)),'*.mat']));
		load(fullfile(workingdir,'eikonal',matfile(1).name));
		evla = eventphv(1).evla;
		evlo = eventphv(1).evlo;
		epi_dist = distance(evla,evlo,mean(lalim),mean(lolim));
		badnum = badnum+1;
		ind = find(~isnan(GV(:)));
		stemp = sprintf('remove %s: id %d, ip %d, dist %f, bias: %f percent, good pixels: %d', char(event_ids(ie)),ie,ip,epi_dist, diff_percent,length(ind));
		disp(stemp)
		GV_mat(:,:,ie,ip) = NaN;
	end
end
end

% calculate the averaged phase velocity again
avgphv = average_GV_mat(GV_mat, raydense_mat, parameters);

% re-Calculate std
for ip=1:length(periods)
	for i = 1:Nx
		for j=1:Ny
			avgphv(ip).GV_std(i,j) = nanstd(GV_mat(i,j,:,ip));
		end
	end
end

% fill in information
for ip=1:length(periods)
	avgphv(ip).xi = xi;
	avgphv(ip).yi = yi;
	avgphv(ip).xnode = xnode;
	avgphv(ip).ynode = ynode;
	avgphv(ip).period = periods(ip);
end

if issmoothmap
	disp(['Smoothing map based on wavelength']);
	for ip=1:length(periods)
		disp(ip);
		D = smooth_wavelength*nanmean(avgphv(ip).GV(:))*periods(ip);
		GV = smoothmap(xi,yi,avgphv(ip).GV,D);
		GV(find(isnan(avgphv(ip).GV))) = NaN;
		avgphv(ip).GV = GV;
	end	
end


save([workingdir,'eikonal_stack_',comp,'.mat'],'avgphv','GV_mat','GV_mat','raydense_mat','event_ids');


% plot section

if isfigure

N=3; M = floor(length(periods)/N)+1;
figure(89)
clf
title('stack for dynamics phv')
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	avgv = nanmean(avgphv(ip).GV(:));
	if isnan(avgv)
		continue;
	end
	caxis([avgv*(1-r) avgv*(1+r)])
	colorbar
	load seiscmap
	colormap(seiscmap)
end
drawnow;

figure(90)
clf
title('Std for dynamics phv')
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_std);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	load seiscmap
	colormap(seiscmap)
	meanstd = nanmean(avgphv(ip).GV_std(:));
	if ~isnan(meanstd)
		caxis([0 2*meanstd])
	end
end
drawnow;


figure(95)
clf
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).sumweight);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	load seiscmap
	colormap(seiscmap)
end
drawnow;
end
