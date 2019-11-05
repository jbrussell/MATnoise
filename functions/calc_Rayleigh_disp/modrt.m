function [td,tu,rd,ru] = modrt(e11,e12,e21,e22,du)

% This function calculates the modified R/T coefficients

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

% Determine the number of layers, N, not including the half space
[m n N] = size(du);

% Initialize a 4x4xN matrix
X = zeros(4,4,N);

% Loop through the first N-1 layers
for j = 1:N-1
   A = [e11(:,:,j+1) -e12(:,:,j); e21(:,:,j+1) -e22(:,:,j)];
   B = [e11(:,:,j) -e12(:,:,j+1); e21(:,:,j) -e22(:,:,j+1)];
   L = [du(:,:,j) zeros(2); zeros(2) du(:,:,j+1)];
   X(:,:,j) = A\(B*L);
end

% Calculate the Nth layer
A = [e11(:,:,N+1) -e12(:,:,N); e21(:,:,N+1) -e22(:,:,N)];
B = [e11(:,:,N)*du(:,:,N); e21(:,:,N)*du(:,:,N)];
X(:,1:2,N) = A\B;

% Extract R/T submatrices
td = X(1:2,1:2,:);
ru = X(1:2,3:4,:);
rd = X(3:4,1:2,:);
tu = X(3:4,3:4,:);
