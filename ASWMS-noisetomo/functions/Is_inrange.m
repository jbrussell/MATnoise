function I = Is_inrange(stla,stlo,parameters)

	lalim = parameters.lalim;
	lolim = parameters.lolim;

	I = 1;
	if stla < min(lalim) || stla > max(lalim)
		I=0;
	end
	if stlo < min(lolim) || stlo > max(lolim)
		I=0;
	end

end
	
