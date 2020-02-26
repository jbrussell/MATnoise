function [ cmap ] = radon_cmap( n )
% Colormap for radon panel
%
% n: number of elements in colormap
%
% JBR 2/2020
%


% % Positions
% % pos01 = [0    0.1667    0.3333    0.5000    0.6667    0.75 0.9 1.0000];
% pos01 = linspace(0,1,8);
% % Colors
% clrs = [255 255 255; 224 0 0; 242 222 0; 120 204 145; 0 157 219; 88 23 137; 70 17 109; 0 0 0]/255;

% Positions
% pos01 = [0    0.1667    0.3333    0.5000    0.6667    0.75 0.9 1.0000];
% Colors
clrs = [255 255 255; 224 50 50; 242 214 0; 120 204 145; 0 157 219; 88 23 137; 0 0 0; 0 0 0]/255;
pos01 = linspace(0,1,size(clrs,1));

% pos01 = linspace(0,1,6);
% % Colors
% clrs = [255 255 255; 120 204 145; 0 157 219; 88 23 137; 70 17 109; 0 0 0]/255;


cmap = customcolormap(pos01, clrs, n);

if 0
    figure(200);
    colorbar;
    colormap(mycolormap);
    axis off;
end


end

