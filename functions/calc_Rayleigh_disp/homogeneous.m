function cvr = homogeneous(cvp,cvs)

% This function calculates the Rayleigh phase velocity in a homogeneous half space

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

% Define Poisson's Ratio  
nu = 0.5*((cvp*cvp-2*cvs*cvs)/(cvp*cvp-cvs*cvs));

% Define Coefficients of Rayleigh's Equation
a =  1;
b = -8;
c =  8*(3-2*(cvs*cvs)/(cvp*cvp));
d = 16*((cvs*cvs)/(cvp*cvp)-1);

% Solution of Rayleigh Equation
p   = [a b c d];
x   = roots(p);
cr  = cvs*sqrt(x);

% Determine which of the roots is correct using the estimated velocity (Achenbach, 1973)
crest = cvs*((0.862+1.14*nu)/(1+nu));
index = find(abs(cr-crest) == min(abs(cr-crest)));
cvr = cr(index);
if isempty(cvr)
   error('No root found for homogeneous half space')
end
