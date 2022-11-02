function para = wpolyfit(x,y,w,N)
% function to apply weighted least square polynomial fitting to the data
%
x=x(:);
y=y(:);
w = w(:);
if length(x)~=length(y) || length(x)~=length(w)
	disp('The input vectors should have same length!');
	return 
end
mat = zeros(length(x),N+1);
for i=1:length(x)
	for j=1:N+1
		mat(i,j) = x(i)^(N+1-j);
	end
end
W = diag(w);

para = (mat'*W*mat)\(mat'*W*y);
end

