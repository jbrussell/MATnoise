function F = break_constraints(Fin, bound, xnode, ynode, N)
% Break flatness (first derivative) or smoothness (second derivative) kernels 
% over a boundary. Boundary is a single contour in map view.
%
% longitude : bound.lon
% latitude  : bound.lat
%
isf1 = 0; % Plot figure of the breaking matrix?
isf99 = 0; % Plot figure of the constraint equations?

Nx = length(xnode);
Ny = length(ynode);

[Xnode,Ynode] = meshgrid(xnode,ynode);
% [Ynode,Xnode] = meshgrid(ynode,xnode);

% Make a matrix of ones and zeros that will be used to delete equations
% from the matrix of prior constraints
mat_break = ones(Nx,Ny);
for ii = 1:length(bound.lat)
    if bound.lat(ii) > max(xnode) || bound.lat(ii) < min(xnode) || bound.lon(ii) > max(ynode) || bound.lon(ii) < min(ynode)
        continue
    end
    xlat = bound.lat(ii);
    ylon = bound.lon(ii);
    [~,Ixs] = min(min(abs(xlat - Xnode)));
    Ix(ii) = Ixs(1);    
    [~,Iys] = min(abs(ylon - Ynode));
    Iy(ii) = Iys(1); 
    mat_break(Ix(ii),Iy(ii)) = 0;
end
if isf1
    figure(1); clf;
    imagesc(ynode,xnode,mat_break); hold on;
    plot(bound.lon,bound.lat,'-w','linewidth',2);
    yticks(xnode);
    xticks(ynode);
    grid on
end

% Do the deleting, but make sure that every grid point still has at least
% one constraint equation attached to it.
F = Fin;
Fgrid = zeros(Nx,Ny);
Iequation = zeros(Nx,Ny);
for ieq = 1:size(Fin,1)
    breakflag = 0;
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            if mat_break(i,j)==0 && F(ieq,n)~=0 && Iequation(i,j)>0
                breakflag = breakflag + 1;
                F(ieq,n) = F(ieq,n) * mat_break(i,j);
            elseif mat_break(i,j)==0 && breakflag>0 && Iequation(i,j)==0
                F(ieq,n) = F(ieq,n) * mat_break(i,j);
            else
                F(ieq,n) = F(ieq,n); 
                Iequation(i,j) = Iequation(i,j)+abs(F(ieq,n));
            end
            Fgrid(i,j)= F(ieq,n);
        end
    end
    if breakflag > 0
        F(ieq,:) = 0;
        Fgrid(:,:) = 0;
    end
    if isf99
        figure(99); clf;
        subplot(1,2,1);
        imagesc(ynode,xnode,Fgrid); hold on;
        plot(bound.lon,bound.lat,'-w','linewidth',2);
        subplot(1,2,2);
        imagesc(ynode,xnode,Iequation); hold on;
        plot(bound.lon,bound.lat,'-w','linewidth',2);
        pause;
    end
end