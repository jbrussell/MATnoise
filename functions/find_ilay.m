function [ilay] = find_ilay(r,Rb)
% [ilay] = find_ilay(r,Rb)
%
% Function to find the indices that each element of vector r would slot
% into vector Rb - originally conceived as a solution to the problem of
% having a series points at different radii and wanting to know which
% layers each of them were in, where the boundaries of the layers
% (including the top and bottom) are given by Rb. 
% 
% For example, if 3 layers were given by boundaries: [0;100;200;300] the
% point 50 would be in layer 1, and 217 would be in layer 3.
% thus, find_ilay([50;217],[0;100;200;300]) = [1;3]
%
% If a point in r is on a boundary, it is put into the upper layer, unless
% it is right at the max, then it is included in outermost layer

if any(r < min(Rb)) || any(r > max(Rb))
    error('r must be within extremes of Rb')
end

r = r(:);
Rb = Rb(:);

N = length(r);
Nlay = length(Rb);

[~,ilay] = max((ones(N,1)*Rb' - r*ones(1,Nlay))>0,[],2);
ilay = ilay-1;
ilay(ilay==0) = Nlay-1; % if any are on outer edge, say in outermost layer

end