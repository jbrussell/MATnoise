function F = flat_kernel_build2(xnode, ynode, N)
% Build flattening kernel for tomography problem using 2 points
% (i.e., smooth first derivative)
% 
% N is the total number of model parameters in the final G matrix
% JBR 7/22/18
%
% xnode : latitude nodes (yes, it appears backwards...)
% ynode : longitude nodes (yes, it appears backwards...)

Nx = length(xnode);
Ny = length(ynode);

% Smoothing in y
[i,j] = ndgrid(1:Nx,1:(Ny-1));
ind = j(:) + Ny*(i(:)-1);
% dy = diff(ynode);
% dy2 = dy(j(:));
dy2 = km2deg(distance(xnode(i),ynode(j),xnode(i),ynode(j+1),referenceEllipsoid('GRS80'))/1000);
Areg = sparse(repmat(ind,1,2),[ind,ind+1], ...
    [-1./dy2, 1./dy2],N,N);
% Smoothing in x
[i,j] = ndgrid(1:(Nx-1),1:Ny);
ind = j(:) + Ny*(i(:)-1);
% dx = diff(xnode);
% dx2 = dx(i(:));
dx2 = km2deg(distance(xnode(i),ynode(j),xnode(i+1),ynode(j),referenceEllipsoid('GRS80'))/1000);
Areg = [Areg;sparse(repmat(ind,1,2),[ind,ind+Ny], ...
    [-1./dx2, 1./dx2],N,N)];
F=Areg;

if 0
    for ieq = 1:size(F,1)
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                Fgrid(i,j)= F(ieq,n);
            end
        end
        figure(99); clf;
        imagesc(xnode,ynode,Fgrid);
        pause;
    end
end

return
end
