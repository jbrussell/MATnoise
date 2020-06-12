function [ cmap ] = tomo_cmap( n )
% Seismic tomography colormap
%
% n: number of elements in colormap
%
% JBR 2/2019
%

% % Positions
% pos01 =  linspace(0,1,7);
% % Colors
% clrs = [0 0 0; 224 0 0; 242 222 0; 255 255 255; 120 204 145; 0 157 219; 88 23 137]/255;

% Positions
pos01 = [0    0.1667    0.3333    0.5000    0.6667    0.75 0.95 1.0000];
% Colors
clrs = [1 1 1; 224 0 0; 242 222 0; 255 255 255; 120 204 145; 0 157 219; 88 23 137; 70 17 109]/255;

cmap = customcolormap(pos01, clrs, n);
cmap = flipud(cmap);

if 0
    figure(200);
    colorbar;
    colormap(mycolormap);
    axis off;
end


end

