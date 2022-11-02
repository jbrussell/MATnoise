% Invert gradient field to obtain the best fitting scalar field.
% An additional constraint is required. We assume a mean value equal to
% zero.
%
% % INPUT:
% xi = meshgrid (latitude)
% yi = mesghrid (longitude)
% dAx = gradient in A along x direction
% dAy = gradient of A along y direction
% Aref = reference amplitude map to damp towards
% dampweight0 = = damping weight towards reference amplitude (as fraction of G matrix norm)
% smweight0 = = smoothing weight (as fraction of G matrix norm)
%
% OUTPUT: 
% A = desired scalar field values
%
% jbrussell 2021
function [Amap]=inv_delm(xi,yi,dAlat_map,dAlon_map,Aref_map,dampweight0,smweight0)

[Nla,Nlo]=size(xi);
A=zeros(Nla,Nlo);
N = Nla*Nlo;

% Do latitudes first (x)
dAlat = zeros(N,1);
Aref = zeros(N,1);
Glat = zeros(N,N);
for ila = 1:Nla
    for ilo = 1:Nlo
        ii=Nlo*(ila-1)+ilo;
        
        % Build data kernel
        if ila-1 == 0
            % endpoint, left derivative
            dla=vdist(xi(ila,ilo),yi(ila,ilo),xi(ila+1,ilo),yi(ila+1,ilo))/1e3;
            Glat(ii,ii) = -1 ./ dla;
            Glat(ii,ii+Nlo) = 1 ./ dla;
        elseif ila+1 > Nla
            % endpoint, right derivative
            dla=vdist(xi(ila-1,ilo),yi(ila-1,ilo),xi(ila,ilo),yi(ila,ilo))/1e3;
            Glat(ii,ii-Nlo) = -1 ./ dla;
            Glat(ii,ii) = 1 ./ dla;
        else
            dla=vdist(xi(ila-1,ilo),yi(ila-1,ilo),xi(ila+1,ilo),yi(ila+1,ilo))/1e3;
            Glat(ii,ii-Nlo) = -1 ./ (dla);
            Glat(ii,ii+Nlo) = 1 ./ (dla);
        end
        
        % Data vector
        dAlat(ii,1) = dAlat_map(ila,ilo);
        Aref(ii,1) = Aref_map(ila,ilo);
    end
end

% Now longitudes (y)
dAlon = zeros(N,1);
Glon = zeros(N,N);
for ila = 1:Nla
    for ilo = 1:Nlo
        ii=Nlo*(ila-1)+ilo;
        
        % Build data kernel
        if ilo-1 == 0
            % endpoint, left derivative
            dlo=vdist(xi(ila,ilo),yi(ila,ilo),xi(ila,ilo+1),yi(ila,ilo+1))/1e3;
            Glon(ii,ii) = -1 ./ dlo;
            Glon(ii,ii+1) = 1 ./ dlo;
        elseif ilo+1 > Nlo
            % endpoint, right derivative
            dlo=vdist(xi(ila,ilo-1),yi(ila,ilo-1),xi(ila,ilo),yi(ila,ilo))/1e3;
            Glon(ii,ii-1) = -1 ./ dlo;
            Glon(ii,ii) = 1 ./ dlo;
        else
            dlo=vdist(xi(ila,ilo-1),yi(ila,ilo-1),xi(ila,ilo+1),yi(ila,ilo+1))/1e3;
            Glon(ii,ii-1) = -1 ./ (dlo);
            Glon(ii,ii+1) = 1 ./ (dlo);
        end
        
        % Data vector
        dAlon(ii,1) = dAlon_map(ila,ilo);
    end
end


% Combine lat and lon
G = [Glat; Glon];
dA = [dAlat; dAlon];
% G = [Glat; Glon; Gconstraint];
% dA = [dAlat; dAlon; dconstraint];

% remove nan values;
inan = find(isnan(dA));
% igood = find(~isnan(dA));
G(inan,:) = [];
dA(inan,:) = [];
% % Remove unsampled grid points
% inode_nodata = find(sum(G,1)==0);
% inode_good = find(sum(G,1)~=0);
% G(:,inode_nodata) = [];

% Damp towards reference amplitude map
Gconstraint = eye(N);
dconstraint = Aref;
J = Gconstraint;
j = dconstraint;
% Rescale weight
NR=norm(J,1);
NA=norm(G,1);
dampweight = dampweight0*NA/NR;

% Add smoothing kernel
xnode = xi(:,1)';
ynode = yi(1,:);
F = smooth_kernel_build(xnode, ynode, N);
f = zeros(size(F,1),1);
% Rescale the smooth kernel
NR=norm(F,1);
NA=norm(G,1);
smweight = smweight0*NA/NR;

% Construct Full G matrix
H = [G; dampweight*J; smweight*F];
h = [dA; dampweight*j; smweight*f];

% Invert for receiver amplitude terms
% std_err = std(sqrt(dAlat_map(:).^2 + dAlon_map(:).^2)/5);
% W = diag(1./std_err).^2;
% F = W.^(0.5)*G;
% f = W.^(0.5)*dA;
% A = (F'*F)\F'*f;
% A = (G'*G)\G'*dA;
A = (H'*H)\H'*h;

% % Estimate chi2 misfit
% dA_pre = G * A;
% e = (dA - dA_pre) ./ dA_err;
% chi2 = (e'*e)/length(dA);

% Form matrix
Amap = nan(size(xi));
for ila = 1:Nla
    for ilo = 1:Nlo
        ii=Nlo*(ila-1)+ilo;
        Amap(ila,ilo) = A(ii);
    end
end

if 0
    figure(999); clf
    ax = worldmap([min(xi(:)) max(xi(:))], [min(yi(:)) max(yi(:))]);
    surfacem(xi,yi,Amap); colorbar;
    
    figure(1000); clf
    subplot(1,3,1);
    worldmap([min(xi(:)) max(xi(:))], [min(yi(:)) max(yi(:))]);
    surfacem(xi,yi,dAlat_map); colorbar;
    
    subplot(1,3,2);
    worldmap([min(xi(:)) max(xi(:))], [min(yi(:)) max(yi(:))]);
    surfacem(xi,yi,dAlon_map); colorbar;
    
    subplot(1,3,3);
    worldmap([min(xi(:)) max(xi(:))], [min(yi(:)) max(yi(:))]);
    surfacem(xi,yi,sqrt(dAlat_map.^2+dAlon_map.^2)); colorbar;
end

end

