function [ x_mesh, y_mesh, z_mesh ] = vec2mesh( xg, yg, zg )
% Convert vector to matrix mesh
%
mmx = length(xg);
mmy = length(yg);
z_mesh=zeros(mmy,mmx);
% [y_mesh,x_mesh]=meshgrid(yg,xg);
[x_mesh,y_mesh]=meshgrid(xg,yg);
% xpts=reshape(x_mesh,mmx*mmy,1);
% ypts=reshape(y_mesh,mmx*mmy,1);
% xpts=reshape(x_mesh,1,mmx*mmy)';
% ypts=reshape(y_mesh,1,mmx*mmy)';

for i=1:mmy
    for j=1:mmx
        joff=j + (i-1)*mmx;
        z_mesh(i,j)=zg(joff);
%         y_mesh(i,j)=ypts(joff);
%         x_mesh(i,j)=xpts(joff);
    end
end

end

