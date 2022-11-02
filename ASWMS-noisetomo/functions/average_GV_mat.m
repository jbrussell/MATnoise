function avgphv = average_GV_mat(GV_mat,raydense_mat,parameters)

periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
min_csgoodratio = parameters.min_csgoodratio;
min_phv_tol = parameters.min_phv_tol;
max_phv_tol = parameters.max_phv_tol;
is_raydense_weight = parameters.is_raydense_weight;
err_std_tol = parameters.err_std_tol;
min_event_num = parameters.min_event_num;
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);

for ip=1:length(periods)
	avgphv(ip).sumV = zeros(Nx,Ny);
	avgphv(ip).sumV_cor = zeros(Nx,Ny);
	avgphv(ip).sumweight = zeros(Nx,Ny);
	avgphv(ip).eventnum = zeros(Nx,Ny);
end

for ip = 1:length(periods)
	for ie = 1:size(GV_mat,3)
		raydense = raydense_mat(:,:,ie,ip);
		GV = GV_mat(:,:,ie,ip);
		ind = find(~isnan(GV));
		if is_raydense_weight
			avgphv(ip).sumV(ind) = avgphv(ip).sumV(ind) + 1./GV(ind).*raydense(ind);
		else
			avgphv(ip).sumV(ind) = avgphv(ip).sumV(ind) + 1./GV(ind);
		end
		avgphv(ip).sumweight(ind) = avgphv(ip).sumweight(ind) + raydense(ind);
		avgphv(ip).eventnum(ind) = avgphv(ip).eventnum(ind)+1;
	end
end

for ip=1:length(periods)
	if is_raydense_weight
		avgphv(ip).GV = avgphv(ip).sumV ./ avgphv(ip).sumweight;
	else
		avgphv(ip).GV = avgphv(ip).sumV ./ avgphv(ip).eventnum;
	end
	ind = find(avgphv(ip).eventnum < min_event_num);
	avgphv(ip).GV(ind) = NaN;
end	

for ip=1:length(periods)
	avgphv(ip).GV = 1./avgphv(ip).GV;
end
