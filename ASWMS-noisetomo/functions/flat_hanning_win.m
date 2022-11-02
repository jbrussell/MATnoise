function windata = flat_hanning_win(taxis,data,winbgt,winendt,tapertime)
% function to apply a box-car window with a hanning taper end.

windata = data(:);

outind = find(taxis < winbgt | taxis > winendt);
windata(outind) = 0;

inind = find(taxis >= winbgt & taxis <= winendt);

taperlength = floor(tapertime./(taxis(2)-taxis(1)));

taper = hanning(taperlength*2);
windata(inind(1:taperlength)) = windata(inind(1:taperlength)).*taper(1:taperlength);
windata(inind(end-taperlength+1:end)) = windata(inind(end-taperlength+1:end)).*taper(taperlength+1:end);

end
