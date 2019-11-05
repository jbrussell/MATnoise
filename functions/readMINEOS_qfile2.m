function [phV,grV,phVq] = readMINEOS_qfile_josh(qfile,periods,mode)
% [phV,grV] = readMINEOS_qfile(qfile,swperiods)
%  
%  Function to read MINEOS qfile file (qfile) with the fundamental mode
%  phase (and group) velocities listed by period, and then interpolate to
%  get velocities at the desired periods (swperiods).

fid = fopen(qfile,'r');
line=fgetl(fid); % read line with # of lines in q model
dum=str2num(line); %#ok<ST2NM>
nQline=dum(1); % will skip this many lines
C = textscan(fid,'%d %d %f %f %f %f %f %f %f %f','Headerlines',nQline); % n,l,omega,Q,?,c,U
fclose(fid);

% find mode of interest
I_mode = find(C{1}==mode);

l = C{2}(I_mode); % mode degree
freq_all = C{3}(I_mode)/2/pi; % all the frequencies
Q = C{4}(I_mode); % Q at each frequency
phV_all = C{6}(I_mode); % phase velocity at each frequency
grV_all = C{7}(I_mode); % group velocity at each frequency
phV_all_q = C{8}(I_mode);
freq_all_q = 1./C{9}(I_mode);
T = C{10}(I_mode);

% desired frequencies
freq_want = 1./periods;

% interpolate to get velocities
phV = interp1(freq_all,phV_all,freq_want);
grV = interp1(freq_all,grV_all,freq_want);

phVq = interp1(freq_all_q,phV_all_q,freq_want);


end

