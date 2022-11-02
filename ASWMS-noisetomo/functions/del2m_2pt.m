function dA=del2m_2pt(xi,yi,A)

[m,n]=size(xi);
dA=nan(m,n);


for i=2:m-1
	for j=2:n-1
		h1=vdist(xi(i-1,j),yi(i-1,j),xi(i+1,j),yi(i+1,j))/1e3/2;
		dA(i,j)=(A(i+1,j)-2*A(i,j)+A(i-1,j))/h1.^2;
		h2=vdist(xi(i,j-1),yi(i,j-1),xi(i,j+1),yi(i,j+1))/1e3/2;
		dA(i,j)=dA(i,j)+(A(i,j+1)-2*A(i,j)+A(i,j-1))/h2.^2;
	end
end

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

