% Estimate travel time surface and gradient as well as amplitude gradient
% for use in eq (4) of Bao et al. (2016) GJI
% github.com/jbrussell
% 2021-05

clear;
setup_parameters

isoverwrite = 1;
isfigure = 0;
is_save_amp_fig = 1;

% Damping and smoothing parameters for gradient map inversion (only used if is_eikonal_ampgrad = 1)
dampweight0 = 0.1; % damping weight towards reference amplitude (as fraction of G matrix norm)
smweight0 = 0.1; % smoothing weight (as fraction of G matrix norm)

is_eikonal_phasegrad = parameters.is_eikonal_phasegrad; % 1: use eikonal tomography values for phase gradient; 0: use travel-time field estimates

r = 0.05;

% input path and files
workingdir = parameters.workingdir;
eventcs_path = [workingdir,'CSmeasure/'];
eikonal_data_path = [workingdir,'eikonal/'];
eikonal_stack_file = [workingdir,'eikonal_stack_',parameters.component];
helmholtz_path = [workingdir,'helmholtz/'];
traveltime_path = [workingdir,'traveltime/'];

if ~exist(traveltime_path,'dir')
	mkdir(traveltime_path);
end

% load stacked phase velocity map
load(eikonal_stack_file);

% set up useful variables
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
xnode = lalim(1):gridsize:lalim(2);
ynode = lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
% tp_var_tol = parameters.tp_var_tol;
alpha_range = parameters.alpha_range;
alpha_search_grid = parameters.alpha_search_grid;
periods = parameters.periods;

eventfiles = dir([eikonal_data_path,'/*_eikonal_',parameters.component,'.mat']);

load seiscmap

if exist('badampsta.lst','file')
	badstnms = textread('badampsta.lst','%s');
	disp('Found Bad amplitude stations:')
	for ista = 1:length(badstnms)
	disp(badstnms(ista))
	end
end

for ie = 1:length(eventfiles)
%for ie = 59
	% read in data for this event
	clear eventphv eventcs traveltime;
	load(fullfile(eikonal_data_path,eventfiles(ie).name));
	eventid = eventphv(1).id;
	matfilename = fullfile(traveltime_path,[eventphv(1).id,'_traveltime_',parameters.component,'.mat']);
	if exist(matfilename,'file') && ~isoverwrite
		disp(['exist: ',matfilename,', skip!'])
		continue;
	end
	disp(eventid);
	eventcsfile = [eventcs_path,'/',eventid,'_cs_',parameters.component,'.mat'];
	if exist(eventcsfile,'file')
		load(eventcsfile);
	else
		disp(['Cannot find CS file for ',eventid,', Skipped']);
		continue;
    end
%     eventhelmholtzfile = [helmholtz_path,'/',eventid,'_helmholtz_',parameters.component,'.mat'];
% 	if exist(eventhelmholtzfile,'file')
% 		load(eventhelmholtzfile);
% 	else
% 		disp(['Cannot find Helmholtz file for ',eventid,', Skipped']);
% 		continue;
% 	end
	if length(eventphv) ~= length(eventcs.avgphv)
		disp('Inconsist of period number for CS file and eikonal file');
		continue;
	end

	for ip = 1:length(eventphv)
		%% fit the amplitude surface
		% reset the arrays
		clear stlas stlos tp
		stlas = eventcs.stlas;
		stlos = eventcs.stlos;
		stnms = eventcs.stnms;
		if exist('badstnms','var')
			list_badstaids = find(ismember(eventcs.stnms,badstnms));
		else
			list_badstaids = [];
		end
		tp = zeros(1,length(stlas));
        if length(eventcs.stlas) ~= length(eventphv(ip).traveltime)
            error('something is wrong... check number of stations');
        end
		for ista = 1:length(eventcs.stlas)
% 			if eventcs.autocor(ista).exitflag(ip)>0
                tp(ista) = eventphv(ip).traveltime(ista);
% 			else
% 				tp(ista) = NaN;
% 			end
        end

		% get rid of bad stations
		badstaids = find(isnan(tp));
		stlas(badstaids) = [];
		stlos(badstaids) = [];
		tp(badstaids) = [];
		badstanum = 0; badstaids = [];
		for ista = 1:length(tp)
			if stlas(ista) < lalim(1) || stlas(ista) > lalim(2) || ...
					stlos(ista) < lolim(1) || stlos(ista) > lolim(2) || ismember(ista,list_badstaids);
				badstanum = badstanum+1;
				badstaids(badstanum) = ista;
				continue;
			end
			dist = distance(stlas(ista),stlos(ista),stlas,stlos);
			dist = deg2km(dist);
			nearstaids = find(dist > parameters.minstadist & dist < parameters.maxstadist );
			nearstaids(find(ismember(nearstaids,badstaids))) = [];
			if isempty(nearstaids)
				badstanum = badstanum+1;
				badstaids(badstanum) = ista;
				continue;
			end
			meantp = median(tp(nearstaids));
% 			if tp(ista) < meantp./tp_var_tol | tp(ista) > meantp.*tp_var_tol
% 				badstanum = badstanum+1;
% 				badstaids(badstanum) = ista;
% 			end
		end
		stlas(badstaids) = [];
		stlos(badstaids) = [];
		tp(badstaids) = [];
        if length(tp(~isnan(tp)))<3
            tpmap = nan(size(xi'));
            mesh_xi = xi';
            mesh_yi = yi';
        else
%             [tpmap,mesh_xi,mesh_yi]=gridfit_jg_geo(stlas,stlos,tp,xnode,ynode,...
%                                 'smooth',2,'regularizer','del4','solver','normal');
                            
            % Use traveltime gradient maps from inversion to solve for
            % maps of traveltime. The inversion includes damping towards 
            % the reference traveltime map and enforcing second derivative smoothing)
            [tpmap_ref,mesh_xi,mesh_yi]=gridfit_jg_geo(stlas,stlos,tp,xnode,ynode,...
                                'smooth',2,'regularizer','del4','solver','normal');
            tp_gradlat = -eventphv(ip).GVx; % phase slowness in x-direction
            tp_gradlon = -eventphv(ip).GVy; % phase slowness in y-direction
            [tpmap] = inv_delm(xi,yi,tp_gradlat,tp_gradlon,tpmap_ref',dampweight0,smweight0);
            tpmap(isnan(tp_gradlat)) = nan;
            tpmap_ref(isnan(tp_gradlat')) = nan;
            % add mean back in such that average of amplitude map equals average of station amplitudes
%                 ampmap = ampmap + nanmean(amps);
            tpmap = tpmap';
            mesh_xi = xi';
            mesh_yi = yi';
        end

		%% Calculate the traveltime and amplitude fields
        
        if is_eikonal_phasegrad == 1
            tp_grad = 1./eventphv(ip).GV'; % phase slowness magnitude
            tp_gradlat = -eventphv(ip).GVx; % phase slowness in x-direction
            tp_gradlon = -eventphv(ip).GVy; % phase slowness in y-direction
            % [~,tp_laplat,~]=delm(xi,yi,tp_gradlat);
            % [~,~,tp_laplon]=delm(xi,yi,tp_gradlon);
            [~,tp_laplat,tp_laplon]=del2m_grad_sph(xi,yi,tp_gradlat,tp_gradlon);
            tp_lap = tp_laplat + tp_laplon;
            tp_ang = 90 - atan2d(tp_gradlat,tp_gradlon);
            tp_gradlat = tp_gradlat';
            tp_gradlon = tp_gradlon';
            tp_lap = tp_lap';
            tp_ang = tp_ang';
            
            tp_grad_err = eventphv(ip).dtau_err';
            tp_gradlat_err = eventphv(ip).dtaux_err';
            tp_gradlon_err = eventphv(ip).dtauy_err';
            tp_lap_err = ( (tp_laplat'.*tp_gradlat_err).^2 + (tp_laplon'.*tp_gradlon_err).^2 ).^0.5; % propagate errors to Laplacian
            phv_err = eventphv(ip).phv_err';
        else
            [tp_grad,tp_gradlat,tp_gradlon]=delm_sph(mesh_xi',mesh_yi',tpmap');
            tp_ang = 90 - atan2d(tp_gradlat,tp_gradlon);  
            
            [tp_lap,tp_laplat,tp_laplon]=del2m_sph(mesh_xi',mesh_yi',tpmap');

            tp_lap = tp_lap';
            % tp_laplat = tp_laplat';
            % tp_laplon = tp_laplon';
            tp_grad = tp_grad';
            tp_gradlat = tp_gradlat';
            tp_gradlon = tp_gradlon';
            tp_ang = tp_ang';
            tp_grad_err = nan(size(tp_grad));
            tp_gradlat_err = nan(size(tp_grad));
            tp_gradlon_err = nan(size(tp_grad));
            tp_lap_err = nan(size(tp_grad));
            phv_err = nan(size(tp_grad));
        end
        
        tp_lap(isnan(eventphv(ip).GV)') = nan;
        tp_grad(isnan(eventphv(ip).GV)') = nan;
        tp_gradlat(isnan(eventphv(ip).GV)') = nan;
        tp_gradlon(isnan(eventphv(ip).GV)') = nan;
        tp_ang(isnan(eventphv(ip).GV)') = nan;
        
		% prepare the avg phase velocity and event phase velocity
		avgGV = avgphv(ip).GV;
		if sum(size(avgGV)==size(xi)) < 2
			avgGV = interp2(avgphv(ip).xi,avgphv(ip).yi,avgphv(ip).GV,xi,yi,'linear',NaN);
		end
		eventGV = eventphv(ip).GV;

		% fill in informations
		traveltime(ip).evla = eventphv(ip).evla;
		traveltime(ip).evlo = eventphv(ip).evlo;
        traveltime(ip).Mw = eventphv(ip).Mw;
		traveltime(ip).raydense = eventphv(ip).raydense;
		traveltime(ip).goodnum = eventphv(ip).goodnum;
		traveltime(ip).badnum = eventphv(ip).badnum;
		traveltime(ip).id = eventphv(ip).id;
		traveltime(ip).xi = xi;
		traveltime(ip).yi = yi;
% 		traveltime(ip).GV_cor = helmholtz(ip).GV_cor;
		traveltime(ip).GV = eventGV;
        traveltime(ip).phv_err = phv_err';
		traveltime(ip).tpmap = tpmap';
        traveltime(ip).tp_lap = tp_lap';
        traveltime(ip).tp_grad = tp_grad';
        traveltime(ip).tp_gradlat = tp_gradlat';
        traveltime(ip).tp_gradlon = tp_gradlon';
        traveltime(ip).tp_ang = tp_ang';
        traveltime(ip).tp = tp;
        traveltime(ip).tp_grad_err = tp_grad_err';
        traveltime(ip).tp_gradlat_err = tp_gradlat_err';
        traveltime(ip).tp_gradlon_err = tp_gradlon_err';
        traveltime(ip).tp_lap_err = tp_lap_err';
		traveltime(ip).period = periods(ip);
        traveltime(ip).stainfo.stlas = stlas;
        traveltime(ip).stainfo.stlos = stlos;
% 		bestalphas(ip,ie) = bestalpha;

		% plot to check
        if isfigure
			figure(37)
			clf
                        set(gcf,'renderer','zbuffer');
% 			subplot(2,2,1)
% 			ax = worldmap(lalim, lolim);
% 			surfacem(xi,yi,eventGV);
%             if ~isnan(nanmean(eventGV(:)))
%                 caxis([nanmean(eventGV(:))*(1-r) nanmean(eventGV(:))*(1+r)])
%             end
% 			colorbar
% 			title('before cor');
% 			subplot(2,2,2)
% 			ax = worldmap(lalim, lolim);
% 			surfacem(xi,yi,GV_cor);
%             if ~isnan(nanmean(GV_cor(:)))
%                 caxis([nanmean(GV_cor(:))*(1-r) nanmean(GV_cor(:))*(1+r)])
%             end
% 			colorbar
% 			title('after cor');
			nanind = find(isnan(eventGV(:)));
			tpmap = tpmap';
			tpmap(nanind) = NaN;
            tp_grad = tp_grad';
			tp_grad(nanind) = NaN;
            tp_gradlat = tp_gradlat';
			tp_gradlat(nanind) = NaN;
            tp_gradlon = tp_gradlon';
			tp_gradlon(nanind) = NaN;
			tp_lap = tp_lap';
			tp_lap(nanind) = NaN;
			subplot(2,2,1)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,tpmap);
			title('travel time map')
            if ~isempty(stlas) 
                plotm(stlas,stlos,'v')
                la_gc = [];
                lo_gc = [];
                for ista = 1:length(stlas)
                    [la,lo]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,stlas(ista),stlos(ista));
                    la_gc = [la_gc; la; nan];
                    lo_gc = [lo_gc; lo; nan];
                end
                plotm(la_gc,lo_gc,'-k');
            end
            colormap(seiscmap)
			colorbar
            subplot(2,2,3)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,tp_grad);
			colorbar
% 			[temp bestalphai] = min(alpha_errs);
			title('\nabla \tau_p')
            drawnow;
			subplot(2,2,4)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,tp_lap);
			colorbar
% 			[temp bestalphai] = min(alpha_errs);
			title('\nabla^2 \tau_p')
            drawnow;
        end % end of isfigure
    end  % loop of period
           
            
    if is_save_amp_fig
        figure(39);
        for ip = 1:length(traveltime)
            nanind = find(isnan(traveltime(ip).GV(:)));
            tpmap = traveltime(ip).tpmap;
			tpmap(nanind) = NaN;
            if ip == 1
                clf;
                set(gcf,'Position',[84           3         744        1022]);
                sgtitle([eventphv(ip).id,' M',num2str(eventphv(ip).Mw)],'fontweight','bold','fontsize',18)

                axes('Position',[.4 .005 .35*.6 .4*.6])
                landareas = shaperead('landareas.shp','UseGeoCoords',true);
                ax = axesm('eqdazim', 'Frame', 'on', 'Grid', 'off');
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';
                box off;
                % setm(ax,'Origin',[mean(lalim),mean(lolim)])
                setm(ax,'Origin',[mean(lalim),mean(lolim)],'FLatLimit',[-125 125]+mean(lalim),'FLonLimit',[],'MapLonLimit',[-125 125]+mean(lolim))
                geoshow(ax, landareas,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
                for ii = [30 60 90 120]
                    [latc,longc] = scircle1(mean(lalim),mean(lolim),ii);
                    plotm(latc,longc,'-','color',[0.6 0.6 0.6],'linewidth',1)
                end
                [la_gcev,lo_gcev]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,mean(lalim),mean(lolim));
                plotm(la_gcev,lo_gcev,'-k','linewidth',2);
                plotm(mean(lalim),mean(lolim),'p','color',[0 0.2 0.4],'MarkerFaceColor',[0 0.5 1],'MarkerSize',24,'linewidth',1);
                plotm(eventphv(ip).evla,eventphv(ip).evlo,'o','color',[0.4 0 0],'MarkerFaceColor',[0.85 0 0],'MarkerSize',10,'linewidth',1);
            end
            la_gc = [];
            lo_gc = [];
            for ista = 1:length(traveltime(ip).stainfo.stlas)
                [la,lo]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,traveltime(ip).stainfo.stlas(ista),traveltime(ip).stainfo.stlos(ista));
                la_gc = [la_gc; la; nan];
                lo_gc = [lo_gc; lo; nan];
            end
            N=3; M = floor(length(periods)/N)+1;
            subplot(M,N,ip)
            ax = worldmap(lalim, lolim);
            surfacem(xi,yi,tpmap);
            if ~isempty(traveltime(ip).stainfo.stlas) 
                plotm(traveltime(ip).stainfo.stlas,traveltime(ip).stainfo.stlos,'v');
                plotm(la_gc,lo_gc,'-k');
            end
            quiverm(xi,yi,traveltime(ip).tp_gradlat,traveltime(ip).tp_gradlon,'-k')
            title([num2str(periods(ip)),' s'],'fontsize',15)
            cb = colorbar;
            clim = cb.Limits;
            colormap(seiscmap)
        end

        figure(40);
        for ip = 1:length(traveltime)
            if ip == 1
                clf;
                set(gcf,'Position',[84           3         744        1022]);
                sgtitle([eventphv(ip).id,' M',num2str(eventphv(ip).Mw)],'fontweight','bold','fontsize',18)

                axes('Position',[.4 .005 .35*.6 .4*.6])
                landareas = shaperead('landareas.shp','UseGeoCoords',true);
                ax = axesm('eqdazim', 'Frame', 'on', 'Grid', 'off');
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';
                box off;
                % setm(ax,'Origin',[mean(lalim),mean(lolim)])
                setm(ax,'Origin',[mean(lalim),mean(lolim)],'FLatLimit',[-125 125]+mean(lalim),'FLonLimit',[],'MapLonLimit',[-125 125]+mean(lolim))
                geoshow(ax, landareas,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
                for ii = [30 60 90 120]
                    [latc,longc] = scircle1(mean(lalim),mean(lolim),ii);
                    plotm(latc,longc,'-','color',[0.6 0.6 0.6],'linewidth',1)
                end
                [la_gcev,lo_gcev]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,mean(lalim),mean(lolim));
                plotm(la_gcev,lo_gcev,'-k','linewidth',2);
                plotm(mean(lalim),mean(lolim),'p','color',[0 0.2 0.4],'MarkerFaceColor',[0 0.5 1],'MarkerSize',24,'linewidth',1);
                plotm(eventphv(ip).evla,eventphv(ip).evlo,'o','color',[0.4 0 0],'MarkerFaceColor',[0.85 0 0],'MarkerSize',10,'linewidth',1);
            end
            la_gc = [];
            lo_gc = [];
            for ista = 1:length(traveltime(ip).stainfo.stlas)
                [la,lo]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,traveltime(ip).stainfo.stlas(ista),traveltime(ip).stainfo.stlos(ista));
                la_gc = [la_gc; la; nan];
                lo_gc = [lo_gc; lo; nan];
            end
            subplot(M,N,ip)
            ax = worldmap(lalim, lolim);
            if ~isempty(traveltime(ip).stainfo.stlas) 
                scatterm(traveltime(ip).stainfo.stlas,traveltime(ip).stainfo.stlos,100,traveltime(ip).tp,'v','filled','markeredgecolor',[0 0 0]);
                plotm(la_gc,lo_gc,'-k');
            end
            title([num2str(periods(ip)),' s'],'fontsize',15)
            colorbar
            caxis(clim);
            colormap(seiscmap)
        end
    end
        
            
	matfilename = fullfile(traveltime_path,[eventphv(1).id,'_traveltime_',parameters.component,'.mat']);
	save(matfilename,'traveltime');
	fprintf('\n');
	disp(['Saved to ',matfilename]);
    if is_save_amp_fig
        figdir = [workingdir,'/figs/traveltime/'];
        if ~exist(figdir)
            mkdir(figdir);
        end
        save2pdf([figdir,eventphv(1).id,'_traveltime_',parameters.component,'.pdf'],39,100);
        save2pdf([figdir,eventphv(1).id,'_traveltime_',parameters.component,'_StaAmps.pdf'],40,100);
    end
end  % loop of events
