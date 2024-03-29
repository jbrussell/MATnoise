% jbrussell 2021
function [dA,dAlat,dAlon]=delm(xi,yi,A)
% dA = magnitude of gradient of A
% dAx = gradient in A along x direction
% dAy = gradient of A along y direction

[m,n]=size(xi);
dAlat=nan(m,n);
dAlon=nan(m,n);

for i=2:m-1
	for j=2:n-1
		h1=vdist(xi(i-1,j),yi(i-1,j),xi(i+1,j),yi(i+1,j))/1e3/2;
		dAlat(i,j)=(A(i+1,j)-A(i-1,j))/h1/2;
		h2=vdist(xi(i,j-1),yi(i,j-1),xi(i,j+1),yi(i,j+1))/1e3/2;
		dAlon(i,j)=(A(i,j+1)-A(i,j-1))/h2/2;
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

