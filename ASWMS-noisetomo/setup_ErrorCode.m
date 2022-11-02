% script to setup error code structure

ErrorCode.low_cohere     = -1;
ErrorCode.high_tp_err    = -2;
ErrorCode.cs_fit_error   = -3;
ErrorCode.sta_fit_error  = -4;
ErrorCode.sta_outofrange = -5;
ErrorCode.sta_lackdata   = -6;
ErrorCode.sta_outofepidist   = -7;
ErrorCode.xcor_outofwin   = -8;
ErrorCode.init_CS_struct   = -10;
ErrorCode.min_stadist_wavelength = -11; % JBR
ErrorCode.max_stadist_wavelength = -12; % JBR
ErrorCode.near_node = -13; % JBR excitation A/A_max less than threshold, meaning event is near a node in the radiation pattern
% ErrorCode.A1_A0_max = -13; % JBR exceed max overtone interference allowed
% ErrorCode.corrfactor_thresh = -14; % JBR
