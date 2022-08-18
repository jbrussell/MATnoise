function [dtds, ier]=ffrq_2dkernel_lalo(kmode,freq,c0,X1,X2,xg,yg,varargin)
%function [dtds, ier]=ffrq_2dkernel(kmode,freq,X1,X2, mmx,mmy,dxgrd,dygrd)
%  Calculate finite frequency kernel for 2D surface wave/ANT, for stations
%  on a grid.  Based on Lin & Ritzwoller 2010 kernel (eq. 11 and 12), 
%  optionally truncated to kth Fresnel zone.
%
%INPUT
%  kmode = 0 for full single-frequency theoretical kernel,
%        = 1 for just 1st Fresnel zone, probably more "appropriate"
%  freq  = circular frequency (Hz)
%  c0    = reference velocity at freq (km/s)
%  X1(2), X2(2) =  (x,y) for the two ray endpoints
%  mmx,mmy   = # of nodes in X,Y for grid
%  dxgrd,dygrd = node spacing in X,Y
%  (optional input:)
%  s1, s2  = structures containing forward (s1) and reverse (s2) calculation 
%            of travel-time surfaces. If not included, will compute
%            analytical kernels
%
%OUTPUT
%  dtds(mmx*mmy,1) = all kernel entries
%  ier = 0 if all OK, <0 if bad
% Some coordinate notes:
%  Grid is cartesian, evenly spaced nodes. (x,y)=(0,0) is 1/2 node from 1st
%    grid point.
%  For node j in x, i in y,  array offset is joff=j+(i-1)*mmx (incr. x,
%  then y)
% GAA 3/17
%
% Josh Russell 4/14: Added ability to calculate empirical kernels if
% forward and reverse phase-delay fields s1 and s2 are included.
%

% A more efficient code for kmode==1 would only do rest of calculations for
% a subset of nodes that are actually going to be nonzero
% Build grid nodes
mmx = length(xg);
mmy = length(yg);
dx = range(xg)/(mmx-1);
dy = range(yg)/(mmy-1);
[ypt_mesh,xpt_mesh]=meshgrid(yg,xg);
xpt=reshape(xpt_mesh,mmx*mmy,1);
ypt=reshape(ypt_mesh,mmx*mmy,1);

% Set some variables
ier=0; 
om=2*pi*freq;
qpi=pi/4;
hpi=pi/2;

%  kernel geometric elements
L=distance(X1(2),X1(1),X2(2),X2(1),referenceEllipsoid('GRS80'))/1000;     % Ray length
r1 = distance(X1(2),X1(1),ypt,xpt,referenceEllipsoid('GRS80'))/1000;
r2 = distance(X2(2),X2(1),ypt,xpt,referenceEllipsoid('GRS80'))/1000;
% L=sqrt(sum((X1-X2).^2));     % Ray length
% r1 = sqrt( (xpt-X1(1)).^2 + (ypt-X1(2)).^2 );
% r2 = sqrt( (xpt-X2(1)).^2 + (ypt-X2(2)).^2 );

type = 'analytical';
if ~isempty(varargin)
    type = 'empirical';
    s1 = varargin{1};
    s2 = varargin{2};
    % Interpolate traveltime surfaces to desired grid
    tt1 = interp2(s1.lat_mesh,s1.lon_mesh,s1.tt_mesh,ypt_mesh,xpt_mesh);
    tt2 = interp2(s2.lat_mesh,s2.lon_mesh,s2.tt_mesh,ypt_mesh,xpt_mesh);
    tt1_2 = interp2(ypt_mesh,xpt_mesh,tt1,X2(2),X2(1)); % index tt1 at receiver location
    tt1=reshape(tt1,mmx*mmy,1);
    tt2=reshape(tt2,mmx*mmy,1);
%     c0 = L./(tt1_2-qpi/om);
end

% Full kernel
kay=om/c0;        % wavenumber
prefac=2*kay;      % Kernel for d(time)/d(slowness)
kp8=kay*pi*8;
afac=sqrt(kp8.*r1.*r2./L);

if strcmp(type,'analytical') 
    % Eq. 11 of Lin & Ritzwoller (2010)
    phs=kay.*(L-r1-r2) + qpi;
elseif strcmp(type,'empirical')
    % Eq. 12 of Lin & Ritzwoller (2010)
    phs=om.*(tt1_2-tt2-tt1) + hpi;
end
dtds = prefac.*cos(phs)./afac;

% truncate higher order Fresnel zones higher than kmode
if kmode>0   
     dtds = dtds.*(abs(phs)<kmode*pi/2);
end


[x_mesh,y_mesh,dtds_mesh] = vec2mesh(xg,yg,dtds);

% Set normalization to that of ray path 
%renorm=(L/c0)/sum(dtds);  %THIS for inverting for fractional slowness varn
% renorm=L/sum(dtds);      % THIS would invert for slowness perturbation GAA 5/17 
renorm = L/trapz(x_mesh(1,:),trapz(y_mesh(:,1),dtds_mesh)); % THIS would invert for slowness perturbation, using 2D integration JBR 4/20
% renorm = L/c0/trapz(x_mesh(1,:),trapz(y_mesh(:,1),dtds_mesh)); % THIS would invert for slowness perturbation, using 2D integration (actually this is probably the correct normalization...) JBR 4/20
dtds=dtds.*renorm*dx*dy; % need to multiply by dx*dy for integration

return
