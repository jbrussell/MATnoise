function cr = modal(freq,thk,dns,cvp,cvs,crmin,crmax)

% This function calculates the modal phase velocities in an elastic,
% vertically heterogeneous medium using search techniques.

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

% Establish global parameters
global TOL MAXROOT NUMINC

% Initialize a matrix to store modal phase velocities
cr = zeros(length(freq),MAXROOT);

% Loop through the frequencies
for j = 1:length(freq)
   if length(cvs) == 1
      cr(j,1) = homogeneous(cvp,cvs);
   else
      numroot = 0;
      om = 2*pi*freq(j);
      
      % Establish the search parameters
      kmax = om/crmin;
      kmin = om/crmax;
      dk = (kmax - kmin)/NUMINC;
      
      % Establish the first and second points
      k1 = kmax;
      f1 = secular(k1,om,thk,dns,cvp,cvs);
      k2 = kmax - dk;
      f2 = secular(k2,om,thk,dns,cvp,cvs);
      
      % Establish an arbitrary high value for kold
      kold = 1.1*kmax;
                  
      % Loop through the remaining points
      for m = 2:NUMINC-1
         k3 = kmax - m*dk;
         f3 = secular(k3,om,thk,dns,cvp,cvs);
         
         % Determine if a minimum is bracketed
         if (f2 < f1) & (f2 < f3)
            
            % Use golden search/parabolic interpolation to refine minimun
            [ktrial,ftrial] = fminbnd('secular',k3,k1,optimset('TolX',1e-12,'Display','off'),om,thk,dns,cvp,cvs);
                           
         	% Check to see if ktrial is a zero and different from the previous zero
         	if (ftrial < TOL & abs((ktrial-kold)/kold) > TOL)
            	numroot = numroot + 1;
            	cr(j,numroot) = om/ktrial;
            	kold = ktrial;
            end
         end   
         
         % Break out of loop of maxroots is reached
         if numroot == MAXROOT
            break;
         end
         
         k1 = k2; f1 = f2;
         k2 = k3; f2 = f3;
         
      end   
   end
end   
