function [ tw0_int, N_zc ] = find_zeros( tw, xsp, r )
% Find zeros of bessel function to determine phase velocity dispersion
%
% J. Russell
% github.com/jbrussell
global tN
global waxis
global twloc

interpmethod = 'linear';
tw1 = interp1(twloc,tw(1:tN),waxis,interpmethod);
x1 = waxis.*tw1;

% tw0 = interp1(b,tw,0);

% Function to index zero crossings of cross-spectra
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
I0 = zci(xsp); % zeros index
w0 = waxis(I0);
% tw0 = tw1(I0);
if isempty(w0)
    tw0_int = [];
    N_zc = 0;
    disp('No zero crossings in cross spectra... skipping');
    return;
end
N_zc = length(I0);

% Get zeros of bessel function
N0 = length(I0);
k = 10;
Zn = besselzero(0, N0+k, 1);

% Calculate dispersion from zero crossings
tw0_int_mat = nan(k+1,tN);
for ii = 1:k+1
    c = r.*w0./Zn(ii:N0+ii-1);
    tw0 = r./c;
    tw0_int = interp1(w0,tw0,twloc);
    tw0_int_mat(ii,:) = tw0_int;
end

% Index value closest to reference velocity model averaging over periods >= 80% of longest period
I_good = find(~isnan(tw0_int_mat(1,:)));
if isempty(I_good)
    tw0_int = [];
    disp('Not enough data to interpolate... skipping');
    return
end
pers = 1./(twloc/2/pi);
I_avg = find(pers >= max(pers(I_good))*0.8 & ~isnan(tw0_int_mat(1,:)));
c_k = mean(r./tw0_int_mat(:,I_avg), 2);
c_ref = mean(r./tw(I_avg));
[~,Ik] = min(abs(c_k-c_ref));
tw0_int = tw0_int_mat(Ik,:);

if 1
    figure(59)
    clf
    subplot(2,1,1)
    hold on
    plot(waxis/2/pi,xsp,'b','linewidth',2);
    plot(w0/2/pi,xsp(I0),'xr','linewidth',2);
    
    subplot(2,1,2);
    hold on;
    plot(twloc/2/pi,r./tw,'-o','color',[0 0 0],'linewidth',2);
%     plot(w0/2/pi,r./tw0,'-or');
    plot(twloc/2/pi,r./tw0_int_mat,'-o','color',[0.7 0.7 0.7]);
    plot(twloc/2/pi,r./tw0_int,'-or','linewidth',2);
    ylim([1.5 4.2]);
end
end

