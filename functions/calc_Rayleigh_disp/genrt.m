function [Td,Rd] = genrt(td,tu,rd,ru)

% This function calculates the generalized R/T coefficients

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

[m n N] = size(td);

% Initialize 2x2xN matrices
Td = zeros(2,2,N);
Rd = zeros(2,2,N);

% Calculate the Td and Rd matrices for the Nth layer
Td(:,:,N) = td(:,:,N);
Rd(:,:,N) = rd(:,:,N);

% Loop through the first N-1 layers in reverse order
for j = N-1:-1:1
   Td(:,:,j) = (eye(2) - ru(:,:,j)*Rd(:,:,j+1))\td(:,:,j);
   Rd(:,:,j) = rd(:,:,j) + tu(:,:,j)*Rd(:,:,j+1)*Td(:,:,j);
end
