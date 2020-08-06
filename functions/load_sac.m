function [S1t,S1raw,S1,tstart1] = load_sac(data1c)
% Load sac data and ensure it is referenced to the same time

S1 = readsac(data1c);
S1raw = S1.DATA1;
S1t = [0:S1.NPTS-1]'*S1.DELTA;
tstart1 = datetime(S1.NZYEAR,1,S1.NZJDAY,S1.NZHOUR,S1.NZMIN,S1.NZSEC,S1.NZMSEC);
            
end

