% Calculate array average PSD from deployment averages
%
% Must first run a00a_psd_staavg_Z_net.m to calculate station averages
%
% jbrussell - 11/2025

clear;
setup_parameters;
IsFigure = 0;
IsFigure_GAUS = 0; % Plot frequency domain filtered and unfiltered

%======================= PARAMETERS =======================%
comp = 'ZZ'; %'ZZ'; %'RR'; %'TT';

stalist = parameters.stalist;
nsta = parameters.nsta;
nsta = length(stalist);
winlength = parameters.winlength;
figpath = [parameters.figpath,'/psd/'];
%fig_winlength_path = [figpath,'window',num2str(winlength),'hr/fullStack/'];
% custom directory names
    fig_winlength_path = [figpath,windir,'/fullStack/'];

%------------ PATH INFORMATION -------------%
% Build File Structure: cross-correlations
% psd_path = parameters.ccfpath;
psd_path = ['./psd/'];
psd_winlength_path = [psd_path,'window',num2str(winlength),'hr/'];
psd_singlestack_path = [psd_winlength_path,'single/'];
psd_singlestack_path = [psd_winlength_path,'dayStack/'];
psd_monthstack_path = [psd_winlength_path,'monthStack/'];
psd_fullstack_path = [psd_winlength_path,'fullStack/'];

psd_stack_path = psd_fullstack_path;

%% Load Depths
STAS = stalist;
LATS = stalat;
LONS = stalon;
DEPTHS = staz;



%%
psd_path = [psd_stack_path,'psd',comp,'/',];
iista = 0;
psd_all = [];
stas_all = {};
%------------ LOAD DATA AND PLOT IN TIME DOMAIN -------------%
for ista1=1:nsta % loop over all stations
    sta1=char(stalist(ista1,:));
    sta1dir=[psd_path,sta1]; % dir to have all cross terms about this central station
    nstapair = 0;
    files = dir([sta1dir,'/*.mat']);
    disp(['Working on sta1: ',sta1])
        
                      
    filename = sprintf('%s/%s_f.mat',sta1dir,sta1);
    
    if ~exist(filename,'file') % check that ccf file exists
        disp(['not exist ',filename])
        continue;
    end
    
    %----------- LOAD DATA -------------%
    data = load(filename);
    dt = data.stapairsinfo.dt;

    if max(abs(data.fftS)) > 1 || min(abs(data.fftS)) < eps
        display(['Bad station, skipping... ',sta1])
        continue
    end

    T = size(data.fftS,2);
    faxis = [0:(T-mod(T-1,2))/2 , -(T-mod(T,2))/2:-1]/dt/T;
    I_fpos = find(faxis>=0);

    Noct_sm = 10; % 1/N octave smoothing (smaller number means more smoothing)
    psdS = smoothSpectrum(abs(data.fftS(I_fpos)),faxis(I_fpos),Noct_sm);
%     fftS_all(nsta,:) = smooth(data.fftS,100);

    if max(psdS) > 1 || min(psdS) < eps
        display(['Bad station after smoothing, skipping... ',sta1])
        continue
    end

    iista = iista + 1;
    
    % Convert back to double-sided (0,+,-)
    psd_all(iista,:) = [psdS(:); flip(psdS(2:end))']';

    stas_all{iista} = sta1;

%     T = length(fftS1Z_full);
%     faxis = [0:(T-mod(T-1,2))/2 , -(T-mod(T,2))/2:-1]/dt/T;
%     ind = find(faxis>0);
%     plot(faxis(ind),abs(fftS1Z_full(ind)),'-','color',[0.7 0.7 0.7],'linewidth',2);
%     plot(faxis(ind),smooth(abs(fftS1Z_full(ind)),100),'-k','linewidth',2);
    
end % ista1

psd_avg = median(psd_all,1);

save([psd_path,'/psd_network_average.mat'],'psd_avg','psd_all','faxis');

%% %----------- PLOT ALL PSDs -------------%

f102 = figure(102); clf;
set(gcf, 'Color', 'w');
box on; hold on;
plot(faxis(I_fpos),psd_all(:,I_fpos),'-','color',[0.7 0.7 0.7]);
plot(faxis(I_fpos),psd_avg(I_fpos),'-r','linewidth',2);
xlabel('Frequency (Hz)','fontsize',15);
ylabel('PSD (displacement)','fontsize',15);
set(gca,'fontsize',15,'XScale','log','YScale','log');

save2pdf(fig[figpath,'psd',comp,'_array_average.pdf']name,f102,500);
% export_fig(figname,'-pdf','-q100','-p0.02','-painters',f102)
% print(f102,'-dpdf',figname);
