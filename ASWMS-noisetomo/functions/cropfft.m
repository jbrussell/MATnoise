function outfft = cropfft(infft,newN)
% function to crop fft due to the change of sample rate
% The time length of input and output fft should be the same. (same frequency interval)
if mod(length(infft),2) ~= mod(newN,2)
	disp('input and output have to be both even');
	return 
end
if length(infft) < newN
	disp('input length has to be longer than output length');
	return
end
if mod(length(infft),2) ==1
	disp('input and output have to be both even');
	return
end

if mod(newN,2) == 0
	outfft(1:newN/2+1) = infft(1:newN/2+1);
	outfft(newN/2+2:newN) = infft(end/2+2:end/2+newN/2);
end

end
