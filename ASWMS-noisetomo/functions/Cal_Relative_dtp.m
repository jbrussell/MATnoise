function travel_time = Cal_relative_dtp(eventcs)
%% function to calculate the relative travel time based on one reference station 
%  with largest connection
%  written by Ge Jin, jinwar@gmail.com

setup_parameters

periods = parameters.periods;

% read in bad station list, if existed
if exist('badsta.lst')
	badstnms = textread('badsta.lst','%s');
	disp('Found Bad stations:')
	disp(badstnms)
end

if exist('badstnms','var')
	badstaids = find(ismember(eventcs.stnms,badstnms));
else
	badstaids = [];
end

for ip = 1:length(periods)
% 	ip
	clear coefmat dt
	sta_connect_num = zeros(length(eventcs.stlas),1);
	coefmat = zeros(1,length(eventcs.stlas));
	goodcsn = 0;
	for ics = 1:length(eventcs.CS)
		if eventcs.CS(ics).isgood(ip) > 0 && sum(ismember([eventcs.CS(ics).sta1 eventcs.CS(ics).sta2],badstaids)) == 0
			goodcsn = goodcsn + 1;
			sta1 = eventcs.CS(ics).sta1;
			sta2 = eventcs.CS(ics).sta2;
			coefmat(goodcsn,sta1) = 1;
			coefmat(goodcsn,sta2) = -1;
			dt(goodcsn) = eventcs.CS(ics).dtp(ip);
		end
	end
	if goodcsn < 1
		tp = zeros(length(eventcs.stlas),1);
		tp(:) = NaN;
		travel_time(ip).tp = tp;
		travel_time(ip).sta_connect_num = sta_connect_num;
		continue;
	end
	dt = dt(:);
	local_connect = diag(coefmat'*coefmat);
	sta_connect_num(find(local_connect==0)) = 1;
	while ~isempty(find(sta_connect_num == 0))
		refsta = find(sta_connect_num == 0,1);
		A = [coefmat; zeros(1,size(coefmat,2))];
		A(end,refsta) = 1;
		b = [dt;0];
		x0 = (A'*A + eye(size(coefmat,2))*1e-6)\(A'*b);
		b = [dt;1];
		x1 = (A'*A + eye(size(coefmat,2))*1e-6)\(A'*b);
		netsize = length(find((x1-x0)>0.9));
		sta_connect_num(find((x1-x0)>0.9)) = netsize;
	end
	[max_netsize beststaid] = max(sta_connect_num);
	A = [coefmat; zeros(1,size(coefmat,2))];
	A(end,beststaid) = 1;
	b = [dt;0];
	x0 = (A'*A + eye(size(coefmat,2))*1e-6)\(A'*b);
	b = [dt;1];
	x1 = (A'*A + eye(size(coefmat,2))*1e-6)\(A'*b);
	tp = x0;
	tp(find((x1-x0)<0.9)) = NaN;
	travel_time(ip).tp = tp;
	travel_time(ip).sta_connect_num = sta_connect_num;
end
