% Setup_parameters for ambient noise processing
%

addpath('../functions_beamforming/');
addpath('../../functions/');

parameters.ccfpath = './ccf_synth/';
[stalist, stalat, stalon, staz] = textread(['../../stations_oldorca.txt'],'%s %f %f %f\n');
parameters.stalist = stalist;
parameters.stalat = stalat;
parameters.stalon = stalon;
parameters.staz = staz;
parameters.nsta = length(parameters.stalist);

