function F = flat_kernel_build(xnode, ynode, N)
% Build flattening kernel for tomography problem using 3 points
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
[i,j] = ndgrid(1:Nx,2:(Ny-1));
ind = j(:) + Ny*(i(:)-1);
% dy = diff(ynode);
% dy2 = dy(j(:));
dy1 = km2deg(distance(xnode(i(:)),ynode(j(:)),xnode(i(:)),ynode(j(:)-1),referenceEllipsoid('GRS80'))/1000);
dy2 = km2deg(distance(xnode(i(:)),ynode(j(:)),xnode(i(:)),ynode(j(:)+1),referenceEllipsoid('GRS80'))/1000);
Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
    [-dy2./(dy1.*(dy1+dy2)), -(dy1-dy2)./(dy1.*dy2), dy1./(dy2.*(dy1+dy2))],N,N);
% Smoothing in x
[i,j] = ndgrid(2:(Nx-1),1:Ny);
ind = j(:) + Ny*(i(:)-1);
% dx = diff(xnode);
% dx2 = dx(i(:));
dx1 = km2deg(distance(xnode(i(:)),ynode(j(:)),xnode(i(:)-1),ynode(j(:)),referenceEllipsoid('GRS80'))/1000);
dx2 = km2deg(distance(xnode(i(:)),ynode(j(:)),xnode(i(:)+1),ynode(j(:)),referenceEllipsoid('GRS80'))/1000);
Areg = [Areg;sparse(repmat(ind,1,3),[ind-Ny,ind,ind+Ny], ...
    [-dx2./(dx1.*(dx1+dx2)), -(dx1-dx2)./(dx1.*dx2), dx1./(dx2.*(dx1+dx2))],N,N)];
% F=Areg;

F=sparse(Nx*Ny*2*2,Nx*Ny*2);
for n=1:size(Areg,1)
    ind=find(Areg(n,:)~=0);
    F(2*n-1,2*ind-1)=Areg(n,ind);
    F(2*n,2*ind)=Areg(n,ind);
end

if 0
    for ieq = 1:size(F,1)
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                Fgridx(i,j)= F(ieq,2*n-1);
                Fgridy(i,j)= F(ieq,2*n);
            end
        end
        figure(99); clf;
        subplot(1,2,1);
        imagesc(xnode,ynode,Fgridx);
        subplot(1,2,2);
        imagesc(xnode,ynode,Fgridy);
        pause;
    end
end

return
end
