% Laplacian of field in spherical coordinates
% jbrussell 2022
function [d2A,d2Alat,d2Alon]=del2m_sph(xi,yi,A)
% xi = latitude grid (deg)
% yi = longitude grid (deg)
% A = Field of interest
%
% d2A = Laplacian
% d2Alat = Second derivative along x direction
% d2Alon = Second derivative along y direction

[m,n]=size(xi);
d2Alat=nan(m,n);
d2Alon=nan(m,n);
dAlat=nan(m,n);
dAlon=nan(m,n);
R = 6371;
for i=2:m-1
	for j=2:n-1
        xh1=sqrt((xi(i-1,j)-xi(i,j))^2 + (yi(i-1,j)-yi(i,j))^2)*pi/180;
        xh2=sqrt((xi(i+1,j)-xi(i,j))^2 + (yi(i+1,j)-yi(i,j))^2)*pi/180;
        dAlat(i,j)=(-A(i-1,j)*xh2/(xh1*(xh1+xh2)) - A(i,j)*(xh1-xh2)/(xh1*xh2) + A(i+1,j)*xh1/(xh2*(xh1+xh2)));
        
        yh1=sqrt((xi(i,j-1)-xi(i,j))^2 + (yi(i,j-1)-yi(i,j))^2)*pi/180;
        yh2=sqrt((xi(i,j+1)-xi(i,j))^2 + (yi(i,j+1)-yi(i,j))^2)*pi/180;
        dAlon(i,j)=(-A(i,j-1)*yh2/(yh1*(yh1+yh2)) - A(i,j)*(yh1-yh2)/(yh1*yh2) + A(i,j+1)*yh1/(yh2*(yh1+yh2)));
    end
end
for i=2:m-1
	for j=2:n-1
        %% Latitude
        xh1=sqrt((xi(i-1,j)-xi(i,j))^2 + (yi(i-1,j)-yi(i,j))^2)*pi/180;
        xh2=sqrt((xi(i+1,j)-xi(i,j))^2 + (yi(i+1,j)-yi(i,j))^2)*pi/180;
		d2Alat(i,j)=(-dAlat(i-1,j)*cosd(xi(i-1,j))*xh2/(xh1*(xh1+xh2)) - dAlat(i,j)*cosd(xi(i,j))*(xh1-xh2)/(xh1*xh2) + dAlat(i+1,j)*cosd(xi(i+1,j))*xh1/(xh2*(xh1+xh2)))/R^2/cosd(xi(i,j));
        
        %% Longitude
        yh1=sqrt((xi(i,j-1)-xi(i,j))^2 + (yi(i,j-1)-yi(i,j))^2)*pi/180;
        yh2=sqrt((xi(i,j+1)-xi(i,j))^2 + (yi(i,j+1)-yi(i,j))^2)*pi/180;
        d2Alon(i,j)=(-dAlon(i,j-1)*yh2/(yh1*(yh1+yh2)) - dAlon(i,j)*(yh1-yh2)/(yh1*yh2) + dAlon(i,j+1)*yh1/(yh2*(yh1+yh2)))/R^2/cosd(xi(i,j))^2;
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

