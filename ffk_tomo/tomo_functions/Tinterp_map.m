function [ map ] = Tinterp_map( map, T_interp )
% Interpolate phase velocity maps to desired frequencies

[Nrow, Ncol] = size(map.lat);
phv_interp = zeros(Nrow,Ncol,length(T_interp));
for irow = 1:Nrow
    for icol = 1:Ncol
        T = map.T_vec(:);
        phv = squeeze(map.phv(irow,icol,:));
        phv_interp(irow,icol,:) = interp1(T,phv,T_interp);
    end
end
map.T_vec = T_interp;
map.phv = phv_interp;

end

