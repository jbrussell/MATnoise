function [e11,e12,e21,e22,du,mu,nus,nup] = psv(thk,dns,cvp,cvs,om,k)

% This function calculates the E and Lambda matrices (up-going and 
% down-going matrices) for the P-SV case. Note that a separate function,
% updown, is provided for calculating the Lambda matrices for use in
% determining the displacement-stress vectors.

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

cvs2 = cvs.^2; cvp2 = cvp.^2;
mu = dns.*cvs2;

e11 = zeros(2,2,length(cvs));
e12 = zeros(2,2,length(cvs));
e21 = zeros(2,2,length(cvs));
e22 = zeros(2,2,length(cvs));
du = zeros(2,2,length(thk));

if om == 0

   kappa = (1.0 + cvs2./cvp2)./(1.0 - cvs2./cvp2);
   kmu = k*mu;
   
   e11(1,1,:) = ones(1,length(cvs));
   e11(1,2,:) = e11(1,1,:);
   e12(1,1,:) = e11(1,1,:);
   e12(1,2,:) = e11(1,1,:);
   e11(2,1,:) = -(kappa - 1.0);
   e11(2,2,:) = e11(1,1,:);
   e12(2,1,:) = -e11(2,1,:);
   e12(2,2,:) = -e11(1,1,:);
   e21(1,1,:) = (kappa - 3.0).*kmu;
   e21(1,2,:) = -2*kmu;
   e22(1,1,:) = -e21(1,1,:);
   e22(1,2,:) = -e21(1,2,:);
   e21(2,1,:) = (kappa - 1.0).*kmu;
   e21(2,2,:) = -2*kmu;
   e22(2,1,:) = e21(2,1,:);
   e22(2,2,:) = e21(2,2,:);
   
   du(1,1,:) = exp(-k*thk);
   du(2,2,:) = exp(-k*thk);
   du(2,1,:) = -k*thk.*exp(-k*thk);

else
   
   k2 = k^2; om2 = om^2;
   
   ks2 = om2./cvs2;
   nus = sqrt(k2-ks2);
   index = find(imag(-i*nus) > 0);
   nus(index) = -nus(index);
   gammas = nus/k;
   
   kp2 = om2./cvp2;
   nup = sqrt(k2-kp2);
   index = find(imag(-i*nup) > 0);
   nup(index) = -nup(index);
   gammap= nup/k;

   chi = 2.0*k - ks2/k;

   e11(1,1,:) = -ones(1,length(cvs));
   e11(1,2,:) = gammas;
   e12(1,1,:) = e11(1,1,:);
   e12(1,2,:) = gammas;
   e11(2,1,:) = -gammap;
   e11(2,2,:) = -e11(1,1,:);
   e12(2,1,:) = gammap;
   e12(2,2,:) = e11(1,1,:);
   e21(1,1,:) = 2*mu.*nup;
   e21(1,2,:) = -mu.*chi;
   e22(1,1,:) = -e21(1,1,:);
   e22(1,2,:) = -e21(1,2,:);
   e21(2,1,:) = -e21(1,2,:);
   e21(2,2,:) = -2*mu.*nus;
   e22(2,1,:) = -e21(1,2,:);
   e22(2,2,:) = e21(2,2,:);
   
   du(1,1,:) = exp(-nup(1:length(thk)).*thk);
   du(2,2,:) = exp(-nus(1:length(thk)).*thk);
   
end
