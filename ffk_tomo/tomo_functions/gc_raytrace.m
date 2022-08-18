function [ tt ] = gc_raytrace( slow, rays, xnode, ynode )
% Calculate travel times along great circle paths
%
% slow: slowness model in vector form (s / km)
%

% Construct ray kernel
G = ray_kernel_build(rays, xnode, ynode);
% Calculate travel times along ray path
tt = G * slow;

end

