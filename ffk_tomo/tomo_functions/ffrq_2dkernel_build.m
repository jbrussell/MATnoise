function G = ffrq_2dkernel_build(kmode, ray, nfreqs,bw_frac,freq,c0, xg, yg, stas, tphase_in)
%input: ray: nray*4 matrix; each row: [x1,y1,x2,y2]
%       xnode and ynode are grid axis for lat and long
%
% kmode: =0 Single frequency, all Fresnel zones
%      : >0 Single frequency calculated up to the kth Fresnel zone
%      : <0 Averge over bandwidth
% 
% Output: G is the kernel for Vx and Vy
% 
% jbrussell 4/2020
global kernel_path
issave = 1; % save kernels?
% kernel_path = ['./SEM2D_FFK_save/',windir,'/']; % path to saved kernels

if ~isempty(tphase_in)
    type = 'empirical';
else
    type = 'analytical';
end

[nrow,ncol]=size(ray);
nray = nrow;
Nx=length(xg);
Ny=length(yg);
Nm = Nx*Ny;
G=spalloc(nray,Nm,nray*Nx); % for each ray, maximum number of pixels to be sampled is 2*Nx

% Check if kernels have already been calculated
if kmode >= 0
    filename = [kernel_path,'2Dkernels_',type,'_',num2str(1./freq),'s_',num2str(kmode),'k_',num2str(Nx),'x',num2str(Ny),'.mat'];
else
    filename = [kernel_path,'2Dkernels_',type,'_',num2str(1./freq),'s_',num2str(kmode),'k_',num2str(bw_frac),'bwfrac_',num2str(Nx),'x',num2str(Ny),'.mat'];
end
if exist(filename)
    temp = load(filename);
    G = temp.G;
    return
end

for ir = 1:nray
    lat1 = ray(ir,1);
    lon1 = ray(ir,2);
    lat2 = ray(ir,3);
    lon2 = ray(ir,4);
    X1 = [lon1 lat1];
    X2 = [lon2 lat2];
    
    switch type
    case 'analytical'
        [dtds, ~]=ffrq_2dkernel_lalo_fb(kmode,nfreqs,bw_frac,freq,c0(ir),X1,X2, xg, yg);
    case 'empirical'
        sta1 = stas{ir,1};
        sta2 = stas{ir,2};
        fil1 = dir([tphase_in,num2str(1/freq),'/phasedelay_',sta1,'*',num2str(1/freq),'s.mat']);
        fil2 = dir([tphase_in,num2str(1/freq),'/phasedelay_',sta2,'*',num2str(1/freq),'s.mat']);
        temp = load([fil1.folder,'/',fil1.name]);
        s1 = temp.run.lalo;
        temp = load([fil2.folder,'/',fil2.name]);
        s2 = temp.run.lalo;
        [dtds, ~]=ffrq_2dkernel_lalo_fb(kmode,nfreqs,bw_frac,freq,c0(ir),X1,X2, xg, yg,s1,s2);
    end
    
	G(ir,:) = dtds;
    
    if 0      
        [x_mesh,y_mesh,vg] = vec2mesh(xg,yg,dtds);
        % PLOT RAYS AND RAY DENSITY
        figure(39); clf;
        Nlvls = 80;
%         lvls = linspace(-0.05,0.05,Nlvls);
        lvls = linspace(-prctile(abs(vg(:)),99.9),prctile(abs(vg(:)),99.9),Nlvls);
        [~,h]=contourf(x_mesh,y_mesh,vg,lvls); hold on; %reshape(dtds,mmx,mmy)' 
        % contour(x_mesh,y_mesh,vg,[0.001 0.001],'-k')
        set(gca,'fontsize',16,'linewidth',1.5)
        set(h,'linestyle','none')
        % colormap('jet')
        colormap(redblue(Nlvls))
%         plot(run.out.xsta,run.out.ysta,'k^')
        plot([X1(1) X2(1)],[X1(2) X2(2)],'k^','markerfacecolor','y')
        plot([X1(1) X2(1)],[X1(2) X2(2)],'k-','linewidth',1.5)
        % colorbar;
        title('Analytical');
        caxis([min(lvls) max(lvls)]);
        axis equal 
        % xlim([75 225])
        pause;
    end
end

if issave
    if ~exist(kernel_path)
        mkdir(kernel_path)
    end
    save(filename,'G','xg','yg');
end


return
end
