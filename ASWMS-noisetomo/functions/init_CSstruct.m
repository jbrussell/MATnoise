function CS = init_CSstruct()

	setup_parameters;
	setup_ErrorCode;
	periods = parameters.periods;

	CS.sta1 = 0;
	CS.sta2 = 0;
	CS.win_cent_t = 0;
	CS.ddist = 0;
	CS.fitpara = zeros(5,length(periods));
	CS.fiterr = zeros(1,length(periods));
	CS.dtp = zeros(1,length(periods));
	CS.dtg = zeros(1,length(periods));
	CS.amp = zeros(1,length(periods));
	CS.w = zeros(1,length(periods));
	CS.sigma = zeros(1,length(periods));
	CS.exitflag = ones(1,length(periods))*ErrorCode.init_CS_struct;

end
