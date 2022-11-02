function windata = flat_hanning_win(taxis,data,wincenter,winlength)
% function to apply a box-car window with a hanning taper end.

windata = data(:);

winbgt = wincenter - winlength/2;
winendt = wincenter + winlength/2;

outind = find(taxis < winbgt | taxis > winendt);
windata(outind) = 0;

inind = find(taxis >= winbgt & taxis <= winendt);

taperlength = length(inind);

taper = hanning(taperlength);
windata(inind) = windata(inind).*taper;

end
