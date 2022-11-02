function [zmap_mesh,mesh_xi,mesh_yi] = gridfit_jg_geo(x,y,z,xnode,ynode,varargin)
% Update to gridfitting function that correctly accounts for
% Earth curvature by converting to ENU coordinates before applying
% the gridfit function. Otherwise, at high latitudes an anomalous 1-theta
% azimtuhal sinudoid is introduced in the Laplacian maps.
%
% jbrussell - 5/26/2022
%
% EXAMPLE:
% [ampmap,mesh_xi,mesh_yi]=gridfit_jg_geo(stlas,stlos,amps,xnode,ynode,...
%                     'smooth',2,'regularizer','del4','solver','normal');

[xi, yi] = ndgrid(xnode,ynode);
mesh_xi = xi';
mesh_yi = yi';
% Convert from geographic to ENU for surface fitting
olat = mean(xnode);
olon = mean(ynode);
[yE, xN, ~] = geodetic2enu(x, y, zeros(size(y)), olat, olon, 0, referenceEllipsoid('GRS80'));
[yim_E, xim_N, ~] = geodetic2enu(xi, yi, zeros(size(xi)), olat, olon, 0, referenceEllipsoid('GRS80'));
dm = min( [min(min(abs(diff(xim_N,1)))), min(min(abs(diff(yim_E',1))))]);
xnodem_N = min(xim_N(:))-dm : dm : max(xim_N(:))+dm;
ynodem_E = min(yim_E(:))-dm : dm : max(yim_E(:))+dm;
[zmap,mesh_xim,mesh_yim]=gridfit_jg(xN/1000,yE/1000,z,xnodem_N/1000,ynodem_E/1000,...
                                    varargin{:});
% Convert ENU back to geographic and sample at even grid spacing
[mesh_xig, mesh_yig, ~] = enu2geodetic(mesh_yim*1000, mesh_xim*1000, zeros(size(mesh_xim)), olat, olon, 0, referenceEllipsoid('GRS80'));
F = scatteredInterpolant(mesh_xig(:),mesh_yig(:),zmap(:));
zmap_mesh = F(mesh_xi,mesh_yi);

end



