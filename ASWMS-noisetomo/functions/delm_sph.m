% Gradient of field in spherical coordinates
% jbrussell 2022
function [dA,dAlat,dAlon]=delm_sph(xi,yi,A)
% xi = latitude grid (deg)
% yi = longitude grid (deg)
% A = Field of interest
%
% dA = Magintude of gradient of A
% dAlat = Gradient along x direction
% dAlon = Gradient along y direction

[m,n]=size(xi);
dAlat=nan(m,n);
dAlon=nan(m,n);
R = 6371;
for i=2:m-1
	for j=2:n-1
% 		h1=km2deg(vdist(xi(i-1,j),yi(i-1,j),xi(i+1,j),yi(i+1,j))/1e3)*pi/180;
%         xh1=km2deg(vdist(xi(i-1,j),yi(i-1,j),xi(i,j),yi(i,j))/1e3)*pi/180;
%         xh2=km2deg(vdist(xi(i+1,j),yi(i+1,j),xi(i,j),yi(i,j))/1e3)*pi/180;
        xh1=sqrt((xi(i-1,j)-xi(i,j))^2 + (yi(i-1,j)-yi(i,j))^2)*pi/180;
        xh2=sqrt((xi(i+1,j)-xi(i,j))^2 + (yi(i+1,j)-yi(i,j))^2)*pi/180;
        dAlat(i,j)=(-A(i-1,j)*xh2/(xh1*(xh1+xh2)) - A(i,j)*(xh1-xh2)/(xh1*xh2) + A(i+1,j)*xh1/(xh2*(xh1+xh2)))/R;
        
% 		h2=km2deg(vdist(xi(i,j-1),yi(i,j-1),xi(i,j+1),yi(i,j+1))/1e3)*pi/180;
%         yh1=km2deg(vdist(xi(i,j-1),yi(i,j-1),xi(i,j),yi(i,j))/1e3)*pi/180;
%         yh2=km2deg(vdist(xi(i,j+1),yi(i,j+1),xi(i,j),yi(i,j))/1e3)*pi/180;
        yh1=sqrt((xi(i,j-1)-xi(i,j))^2 + (yi(i,j-1)-yi(i,j))^2)*pi/180;
        yh2=sqrt((xi(i,j+1)-xi(i,j))^2 + (yi(i,j+1)-yi(i,j))^2)*pi/180;
        dAlon(i,j)=(-A(i,j-1)*yh2/(yh1*(yh1+yh2)) - A(i,j)*(yh1-yh2)/(yh1*yh2) + A(i,j+1)*yh1/(yh2*(yh1+yh2)))/R/cosd(xi(i,j));
	end
end
dA = (dAlon.^2 + dAlat.^2).^0.5;

% now fill in the edges
for i=1:m
	dA(i,1)=dA(i,2);
	dA(i,n)=dA(i,n-1);
end
for j=1:n
	dA(1,j)=dA(2,j);
	dA(m,j)=dA(m-1,j);
end

end

