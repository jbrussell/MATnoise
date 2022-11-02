function dA=del2m_cart(xi,yi,A)

[m,n]=size(xi);
dAx=nan(m,n);
dAy=nan(m,n);

for i=2:m-1
	for j=2:n-1
		xh1=sqrt((xi(i-1,j)-xi(i,j)).^2 + (yi(i-1,j)-yi(i,j)).^2);
	    xh2=sqrt((xi(i+1,j)-xi(i,j)).^2 + (yi(i+1,j)-yi(i,j)).^2);
		dAx(i,j)=2.*(xh2*A(i-1,j)-(xh1+xh2)*A(i,j)+xh1*A(i+1,j))/(xh1*xh2*(xh1+xh2));
	    
	    yh1=sqrt((xi(i,j-1)-xi(i,j)).^2 + (yi(i,j-1)-yi(i,j)).^2);
	    yh2=sqrt((xi(i,j+1)-xi(i,j)).^2 + (yi(i,j+1)-yi(i,j)).^2);
		dAy(i,j)=2.*(yh2*A(i,j-1)-(yh1+yh2)*A(i,j)+yh1*A(i,j+1))/(yh1*yh2*(yh1+yh2));
	end
end
dA = dAx + dAy;

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

