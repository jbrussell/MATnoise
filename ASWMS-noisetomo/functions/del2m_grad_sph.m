% jbrussell 2021
function [d2A,d2Alat,d2Alon]=del2m_grad_sph(xi,yi,dAlat,dAlon)
% dA = magnitude of gradient of A
% dAx = gradient in A along x direction
% dAy = gradient of A along y direction

[m,n]=size(xi);
d2Alat=nan(m,n);
d2Alon=nan(m,n);
R = 6371;
for i=2:m-1
	for j=2:n-1
        hkm = vdist(xi(i-1,j),yi(i-1,j),xi(i+1,j),yi(i+1,j))/1e3;
%         hdeg = km2deg(hkm);
        hdeg = sqrt((xi(i-1,j)-xi(i+1,j)).^2 + (yi(i-1,j)-yi(i+1,j)).^2);
        h1=hdeg*pi/180;
% 		h1=km2deg(vdist(xi(i-1,j),yi(i-1,j),xi(i+1,j),yi(i+1,j))/1e3)*pi/180;
		d2Alat(i,j)=(dAlat(i+1,j)*cosd(xi(i+1,j))-dAlat(i-1,j)*cosd(xi(i-1,j)))*(hkm/hdeg*180/pi)/h1/R^2/cosd(xi(i,j));
        
        hkm = vdist(xi(i,j-1),yi(i,j-1),xi(i,j+1),yi(i,j+1))/1e3;
%         hdeg = km2deg(hkm);
        hdeg = sqrt((xi(i,j-1)-xi(i,j+1)).^2 + (yi(i,j-1)-yi(i,j+1)).^2);
        h2=hdeg*pi/180;
% 		h2=km2deg(vdist(xi(i,j-1),yi(i,j-1),xi(i,j+1),yi(i,j+1))/1e3)*pi/180;
		d2Alon(i,j)=(dAlon(i,j+1)-dAlon(i,j-1))*(hkm/hdeg*180/pi)/h2/R^2/cosd(xi(i,j))^2;
	end
end
d2A = d2Alat + d2Alon;

% now fill in the edges
for i=1:m
	d2A(i,1)=d2A(i,2);
	d2A(i,n)=d2A(i,n-1);
end
for j=1:n
	d2A(1,j)=d2A(2,j);
	d2A(m,j)=d2A(m-1,j);
end

end

