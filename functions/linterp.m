function [ YI ] = linterp(X,Y,XI)
% YI = LINTERP(X,Y,XI) interpolates to find YI, the values of the
%     underlying function Y at the points in the array XI. X must be a
%     vector of length N.
% this function differs from the simple interp1 matlab function in that it
% can accept a vector x with multiple values of y (e.g. at the top of one
% layer and the bottom of another)
%
% Z. Eilon   May 2015
X = X(:); Y = Y(:); XI = XI(:);

YI = zeros(size(XI));

% use find_ilay to find where things fit - will not work properly for when
% a member of X matches a member of XI
[ilay] = find_ilay(XI,X);

% linearly interpolate
YI = (XI - X(ilay)).*(Y(ilay+1)-Y(ilay))./(X(ilay+1)-X(ilay)) + Y(ilay);
% now sort out coincident elements
olap = intersect(X,XI);
for i = 1:length(olap)
YI(XI==olap(i)) = mean(Y(X==olap(i)));
end

end

