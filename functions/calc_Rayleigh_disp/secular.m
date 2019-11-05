function d = secular(k,om,thk,dns,cvp,cvs)

% This function calculates the absolute value of the secular function for
% a particular frequency and wavenumber.

% Copyright 1999 by Glenn J. Rix and Carlo G. Lai

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Check to see if the trial phase velocity is equal to the shear wave velocity
% or compression wave velocity of one of the layers
epsilon = 0.0001;
while any(abs(om/k-cvs)<epsilon) | any(abs(om/k-cvp)<epsilon)
   k = k * (1+epsilon);
end   

[e11,e12,e21,e22,du,mu,nus,nup] = psv(thk,dns,cvp,cvs,om,k);
[td,tu,rd,ru] = modrt(e11,e12,e21,e22,du);
[Td,Rd] = genrt(td,tu,rd,ru);

% Note that the absolute value of the secular function is calculated
d = abs(det(e21(:,:,1) + e22(:,:,1)*du(:,:,1)*Rd(:,:,1))/(nus(1)*nup(1)*mu(1)^2));
