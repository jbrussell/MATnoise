function [t1, t2, d] = readsac(filenames)
% S = READSAC(FILENAMES); reads the given SAC data files. FILENAMES contains
% the names of SAC data files to be read. It may consist of a string for a
% single file or a cell array of strings for multiple file names. Wildcard '*' 
% is accepted. If FILENAMES is omitted, all SAC data files in the current 
% directory are read. The data sets and header variables of the SAC data
% files are returned in structure S.
%
% [INDEP, DEP] = READSAC(FILENAMES); reads the SAC data files and returns
% the independent and dependent variables in INDEP and DEP. INDEP and DEP are
% either cell arrays if more than one file is specified and data lengths are
% different for different files, or double vectors (matrices) if only one file is
% specified or if data have the same length for multiple input files.
% [INDEP, DEP1, DEP2] = READSAC(FILENAMES); reads spectral type of SAC data files and
% [INDEP1, INDEP2, DEP] = READSAC(FILENAMES); reads xyz (spectrogram) type 
% of SAC data files.
%
% *** As of this version, SAC does not define general xyz type files. ***
%
% Structure S is defined to have following fields corresponding to SAC
% header variables and the data set except the first one. Fields that are 
% set to undefined in the SAC data files are set to NaN or '', depending on
% the type of the variable, in the S structure. 
%
%    FILENAME   - File name of SAC data file; Path to the file, if exists,
%                 is removed.
%    DELTA      - Data sampling interval
%    DEPMIN     - Minimum value of dependent variable
%    DEPMAX     - Maximum value of dependent variable
%    SCALE      - Multiplying factor for dependent variable
%    ODELTA     - Observed increment if different from nominal
%    B          - Beginning value of independent variable     
%    E          - Ending value of independent variable     
%    O          - Event origin time                        
%    A          - First arrival time                       
%    INTERNAL1  - First internal variable                  
%    T0         - First user-defined time picks or markers 
%    T1         - Second user-defined time picks or markers 
%    T2         - Third user-defined time picks or markers 
%    T3         - Fourth user-defined time picks or markers 
%    T4         - Fifth user-defined time picks or markers 
%    T5         - Sixth user-defined time picks or markers 
%    T6         - Seventh user-defined time picks or markers 
%    T7         - Eighth user-defined time picks or markers 
%    T8         - Ninth user-defined time picks or markers 
%    T9         - Tenth user-defined time picks or markers 
%    F          - Final or end of event time               
%    RESP0      - First instrument response parameter      
%    RESP1      - Second instrument response parameter      
%    RESP2      - Third instrument response parameter      
%    RESP3      - Fourth instrument response parameter      
%    RESP4      - Fifth instrument response parameter      
%    RESP5      - Sixth instrument response parameter      
%    RESP6      - Seventh instrument response parameter      
%    RESP7      - Eighth instrument response parameter      
%    RESP8      - Ninth instrument response parameter      
%    RESP9      - Tenth instrument response parameter      
%    STLA       - Station latitude                         
%    STLO       - Station longitude                        
%    STEL       - Station elevation                        
%    STDP       - Station depth                            
%    EVLA       - Event latitude                           
%    EVLO       - Event longitude                          
%    EVEL       - Event elevation                          
%    EVDP       - Event depth                              
%    MAG        - Event magnitude                          
%    USER0      - First user-defined variable                   
%    USER1      - Second user-defined variable                   
%    USER2      - Third user-defined variable                   
%    USER3      - Fourth user-defined variable                   
%    USER4      - Fifth user-defined variable               
%    USER5      - Sixth user-defined variable                   
%    USER6      - Seventh user-defined variable                   
%    USER7      - Eighth user-defined variable                   
%    USER8      - Ninth user-defined variable                   
%    USER9      - Tenth user-defined variable                   
%    DIST       - Station-to-event distance (km)
%    AZ         - Event-to-station azimuth (degree)
%    BAZ        - Station-to-event azimuth (degree)
%    GCARC      - Station-to-event great-circle arc length (degree)
%    INTERNAL2  - Second internal variable                              
%    INTERNAL3  - Third internal variable                              
%    DEPMEN     - Mean value of dependent variable                     
%    CMPAZ      - Component azimuth                                    
%    CMPINC     - Component incident angle                             
%    XMINIMUM   - Minimum value of X (spectral file only)              
%    XMAXIMUM   - Maximum value of X (spectral file only)              
%    YMINIMUM   - Minimum value of Y (spectral file only)              
%    YMAXIMUM   - Maximum value of Y (spectral file only)              
%    UNUSED1    - First unused variable                                
%    UNUSED2    - Second unused variable                                
%    UNUSED3    - Third unused variable                                
%    UNUSED4    - Fourth unused variable                                
%    UNUSED5    - Fifth unused variable                                
%    UNUSED6    - Sixth unused variable                                
%    UNUSED7    - Seventh unused variable                                
%    NZYEAR     - GMT year corresponding to reference time               
%    NZJDAY     - GMT julian day corresponding to reference time               
%    NZHOUR     - GMT hour corresponding to reference time               
%    NZMIN      - GMT minute corresponding to reference time               
%    NZSEC      - GMT second corresponding to reference time               
%    NZMSEC     - GMT milisecond corresponding to reference time               
%    NVHDR      - Header version
%    NORID      - Origin ID     
%    NEVID      - Event ID      
%    NPTS       - Number of data points
%    INTERNAL4  - Fourth internal variable
%    NWFID      - Waveform ID
%    NXSIZE     - Spectral length (spectral file only)
%    NYSIZE     - Spectral width (spectral file only)
%    UNUSED8    - Eighth unused variable             
%    IFTYPE     - Type of file                       
%    IDEP       - Type of dependent variable            
%    IZTYPE     - Reference-time equivalence         
%    UNUSED9    - Ninth unused variable              
%    IINST      - Type of recording instrument       
%    ISTREG     - Station geographic region          
%    IEVREG     - Event geographic region          
%    IEVTYP     - Type of event                    
%    IQUAL      - Quality of data                  
%    ISYNTH     - Synthetic data flag              
%    IMAGTYP    - Magnitude type                   
%    IMAGSRC    - Magnitude source                 
%    UNUSED10   - Tenth unused variable              
%    UNUSED11   - Eleventh unused variable              
%    UNUSED12   - Twelveth unused variable              
%    UNUSED13   - Thirteenth unused variable              
%    UNUSED14   - Fourteenth unused variable              
%    UNUSED15   - Fifteenth unused variable              
%    UNUSED16   - Sixteenth unused variable              
%    UNUSED17   - Seventeenth unused variable              
%    LEVEN      - True if data is evenly spaced    
%    LPSPOL     - True if station polarity follows left-hand rule
%    LOVROK     - True if it is ok to overwrite this file on disk
%    LCALDA     - True if DIST, AZ, BAZ and GCARC are to be calculated from 
%                 station and event coordinates
%    UNUSED18   - Eighteenth unused variable              
%    KSTNM      - Station name
%    KEVNM      - Event name
%    KHOLE      - Hole ID for nuclear test
%    KO         - Event origin time ID    
%    KA         - First arrival time ID    
%    KT0        - First user-defined time pick ID    
%    KT1        - Second user-defined time pick ID    
%    KT2        - Third user-defined time pick ID    
%    KT3        - Fourth user-defined time pick ID    
%    KT4        - Fifth user-defined time pick ID    
%    KT5        - Sixth user-defined time pick ID    
%    KT6        - Seventh user-defined time pick ID    
%    KT7        - Eighth user-defined time pick ID    
%    KT8        - Ninth user-defined time pick ID    
%    KT9        - Tenth user-defined time pick ID    
%    KF         - Final or end event time ID    
%    KUSER0     - First user-defined text string
%    KUSER1     - Second user-defined text string
%    KUSER2     - Third user-defined text string
%    KCMPNM     - Component name
%    KNETWK     - Network name
%    KDATRD     - Date data were read onto computer
%    KINST      - Generic name of recording instrument
%    DATA1      - First data block
%
% If the SAC data file is for spectral or unevenly-spaced xy data, following
% fields may exist.
%
%    DATA2      - Second data block
%
% Values of enumerated I-type header variables are converted from integers to
% corresponding alphanumeric names (all upper case).

% *************************************************************************
% This program was developed originally from the program readsac.m of the 
% package MatSeis 1.6 by Sandia National Laboratories.
%
% Mark Harris, mharris@sandia.gov
% Copyright (c) 1996-2001 Sandia National Laboratories. All rights reserved.
% *************************************************************************
%
% Copyright, 19, The Board of Governors of the Los Alamos National Security, LLC.
% This software was produced under a U. S. Government contract (DE-AC52-06NA25396)
% by Los Alamos National Laboratory, which is operated by the Los Alamos
% National Security, LLC for the U. S. Department of Energy. The U. S. Government
% is licensed to use, reproduce, and distribute this software. Permission is granted
% to the public to copy and use this software without charge, provided that 
% this Notice and any statement of authorship are reproduced on all copies.
% Neither the Government nor the LANS makes any warranty, express or implied, 
% or assumes any liability or responsibility for the use of this software.
%
%        Xiaoning Yang, xyang@lanl.gov   2002, 2008
%

% process input argument
if nargin < 1
   filenames = dir;
   filenames = filenames(~[filenames.isdir]);
   filenames = {filenames.name};
else
   % check input argument
   error(nargchk(1, 1, nargin), 'string')
   if ~iscell(filenames) && (~ischar(filenames) || size(filenames, 1) ~= 1)
      error('Input must be a string or a cell array of strings !!!')
   elseif ischar(filenames)
      filenames = cellstr(filenames);
   end
   
   % treat wildcard
   filename_temp = {''};
   for i = 1:length(filenames)
      if ismember('*', filenames{i})
         [pathstr, name] = fileparts(filenames{i});   % separate path from name
         filenamei = dir(deblank(filenames{i}));
         if isempty(pathstr)
            filename_temp = [filename_temp; {filenamei.name}'];
         else  % add path to filenames
            filenamei = char({filenamei.name});
            pathstr = [pathstr, filesep];
            pathstr = repmat(pathstr, size(filenamei, 1), 1);
            filenamei = [pathstr filenamei];
            filename_temp = [filename_temp; cellstr(filenamei)];
         end
      else
         filename_temp = [filename_temp, filenames(i)];
      end
   end
   filenames = filename_temp(2:end);
end
nfiles = length(filenames);

% initialize output structure
s = sacstruct(nfiles);

% used to flag successfully read SAC files
count = false(nfiles, 1);

% read all files
for f = 1:nfiles

   % open SAC file
   fid = fopen(filenames{f}, 'r');

   if fid ~= -1
      % read SAC file header
      [A count1] = fread(fid, [70 1], 'float32');
      [B count2] = fread(fid, [40 1], 'int32');
      if ~ismember(B(7), 1:6) % non-native byte order
         fclose(fid);
         [comp, mxsize, endian] = computer;
         if strcmp(endian, 'B')
            fid = fopen(filenames{f}, 'r', 'l');
         else
            fid = fopen(filenames{f}, 'r', 'b');
         end
         [A count1] = fread(fid, [70 1], 'float32');
         [B count2] = fread(fid, [40 1], 'int32');
      end
      [C count3] = fread(fid, [1 192], 'char');
      C = char(C);

      % process valid SAC data file
      if count1 == 70 && count2 == 40 && count3 == 192 && ...
         B(10) >= 1 && (B(16) == 1 || B(16) == 2 || B(16) == 3 || ...
         B(16) == 4 || B(16) == 51) && (B(36) == 0 || B(36) == 1)

         % set undefined values
         A(A == -12345.0) = NaN;
         B(B == -12345) = NaN;
         C = cellstr(reshape(C, 8, 24)');
         C(strncmp('-12345', C, 6)) = {' '};

         % fill output structure
         [temp, name, ext] = fileparts(filenames{f});
         s(f).FILENAME = [name, ext];
         s(f).DELTA = A(1);
         s(f).DEPMIN = A(2);
         s(f).DEPMAX = A(3);
         s(f).SCALE = A(4);
         s(f).ODELTA = A(5);
         s(f).B = A(6);
         s(f).E = A(7);
         s(f).O = A(8);
         s(f).A = A(9);
         s(f).INTERNAL1 = A(10);
         s(f).T0 = A(11);
         s(f).T1 = A(12);
         s(f).T2 = A(13);
         s(f).T3 = A(14);
         s(f).T4 = A(15);
         s(f).T5 = A(16);
         s(f).T6 = A(17);
         s(f).T7 = A(18);
         s(f).T8 = A(19);
         s(f).T9 = A(20);
         s(f).F = A(21);
         s(f).RESP0 = A(22);
         s(f).RESP1 = A(23);
         s(f).RESP2 = A(24);
         s(f).RESP3 = A(25);
         s(f).RESP4 = A(26);
         s(f).RESP5 = A(27);
         s(f).RESP6 = A(28);
         s(f).RESP7 = A(29);
         s(f).RESP8 = A(30);
         s(f).RESP9 = A(31);
         s(f).STLA = A(32);
         s(f).STLO = A(33);
         s(f).STEL = A(34);
         s(f).STDP = A(35);
         s(f).EVLA = A(36);
         s(f).EVLO = A(37);
         s(f).EVEL = A(38);
         s(f).EVDP = A(39);
         s(f).MAG = A(40);
         s(f).USER0 = A(41);
         s(f).USER1 = A(42);
         s(f).USER2 = A(43);
         s(f).USER3 = A(44);
         s(f).USER4 = A(45);
         s(f).USER5 = A(46);
         s(f).USER6 = A(47);
         s(f).USER7 = A(48);
         s(f).USER8 = A(49);
         s(f).USER9 = A(50);
         s(f).DIST = A(51);
         s(f).AZ = A(52);
         s(f).BAZ = A(53);
         s(f).GCARC = A(54);
         s(f).INTERNAL2 = A(55);
         s(f).INTERNAL3 = A(56);
         s(f).DEPMEN = A(57);
         s(f).CMPAZ = A(58);
         s(f).CMPINC = A(59);
         s(f).XMINIMUM = A(60);
         s(f).XMAXIMUM = A(61);
         s(f).YMINIMUM = A(62);
         s(f).YMAXIMUM = A(63);
         s(f).UNUSED1 = A(64);
         s(f).UNUSED2 = A(65);
         s(f).UNUSED3 = A(66);
         s(f).UNUSED4 = A(67);
         s(f).UNUSED5 = A(68);
         s(f).UNUSED6 = A(69);
         s(f).UNUSED7 = A(70);
         s(f).NZYEAR = B(1);
         s(f).NZJDAY = B(2);
         s(f).NZHOUR = B(3);
         s(f).NZMIN = B(4);
         s(f).NZSEC = B(5);
         s(f).NZMSEC = B(6);
         s(f).NVHDR = B(7);
         s(f).NORID = B(8);
         s(f).NEVID = B(9);
         s(f).NPTS = B(10);
         s(f).INTERNAL4 = B(11);
         s(f).NWFID = B(12);
         s(f).NXSIZE = B(13);
         s(f).NYSIZE = B(14);
         s(f).UNUSED8 = B(15);
         switch B(16)
            case 1, filetype = 'ITIME';
            case 2, filetype = 'IRLIM';
            case 3, filetype = 'IAMPH';
            case 4, filetype = 'IXY';
            case 51, filetype = 'IXYZ';
         end
         s(f).IFTYPE = filetype;
         switch B(17)
            case 5, sigtype = 'IUNKN';
            case 6, sigtype = 'IDISP';
            case 7, sigtype = 'IVEL';
            case 8, sigtype = 'IACC';
            case 50, sigtype = 'IVOLTS';
            otherwise, sigtype = '';
         end
         s(f).IDEP = sigtype;
         switch B(18)
            case 5, refequiv = 'IUNKN';
            case 9, refequiv = 'IB';
            case 10, refequiv = 'IDAY';
            case 11, refequiv = 'IO';
            case 12, refequiv = 'IA';
            case 13, refequiv = 'IT0';
            case 14, refequiv = 'IT1';
            case 15, refequiv = 'IT2';
            case 16, refequiv = 'IT3';
            case 17, refequiv = 'IT4';
            case 18, refequiv = 'IT5';
            case 19, refequiv = 'IT6';
            case 20, refequiv = 'IT7';
            case 21, refequiv = 'IT8';
            case 22, refequiv = 'IT9';
            otherwise, refequiv = '';
         end
         s(f).IZTYPE = refequiv;
         s(f).UNUSED9 = B(19);
         s(f).IINST = B(20);
         s(f).ISTREG = B(21);
         s(f).IEVREG = B(22);
         switch B(23)
            case 5, evtype = 'IUNKN';
            case 37, evtype = 'INUCL';
            case 38, evtype = 'IPREN';
            case 39, evtype = 'IPOSTN';
            case 40, evtype = 'IQUAKE';
            case 41, evtype = 'IPREQ';
            case 42, evtype = 'IPOSTQ';
            case 43, evtype = 'ICHEM';
            case 44, evtype = 'IOTHER';
            case 70, evtype = 'IQB';
            case 71, evtype = 'IQB1';
            case 72, evtype = 'IQB2';
            case 73, evtype = 'IQBX';
            case 74, evtype = 'IQMT';
            case 75, evtype = 'IEQ';
            case 76, evtype = 'IEQ1';
            case 77, evtype = 'IEQ2';
            case 78, evtype = 'IME';
            case 79, evtype = 'IEX';
            case 80, evtype = 'INU';
            case 81, evtype = 'INC';
            case 82, evtype = 'IO_';
            case 83, evtype = 'IL';
            case 84, evtype = 'IR';
            case 85, evtype = 'IT';
            case 86, evtype = 'IU';
            otherwise, evtype = '';
         end
         s(f).IEVTYP = evtype;
         switch B(24)
            case 44, quality = 'IOTHER';
            case 45, quality = 'IGOOD';
            case 46, quality = 'IGLCH';
            case 47, quality = 'IDROP';
            case 48, quality = 'ILOWSN';
            otherwise, quality = '';
         end
         s(f).IQUAL = quality;
         switch B(25)
            case 49, synth = 'IRLDTA';
            otherwise, synth = '';
         end
         s(f).ISYNTH = synth;
         switch B(26)
            case 52, magtype = 'IMB';
            case 53, magtype = 'IMS';
            case 54, magtype = 'IML';
            case 55, magtype = 'IMW';
            case 56, magtype = 'IMD';
            case 57, magtype = 'IMX';
            otherwise, magtype = '';
         end
         s(f).IMAGTYP = magtype;
         switch B(27)
            case 58, magsrc = 'INEIC';
            case 59, magsrc = 'IPDE';
            case 60, magsrc = 'IISC';
            case 61, magsrc = 'IREB';
            case 62, magsrc = 'IUSGS';
            case 63, magsrc = 'IBRK';
            case 64, magsrc = 'ICALTECH';
            case 65, magsrc = 'ILLNL';
            case 66, magsrc = 'IEVLOC';
            case 67, magsrc = 'IJSOP';
            case 68, magsrc = 'IUSER';
            case 69, magsrc = 'IUNKNOWN';
            otherwise, magsrc = '';
         end
         s(f).IMAGSRC = magsrc;
         s(f).UNUSED10 = B(28);
         s(f).UNUSED11 = B(29);
         s(f).UNUSED12 = B(30);
         s(f).UNUSED13 = B(31);
         s(f).UNUSED14 = B(32);
         s(f).UNUSED15 = B(33);
         s(f).UNUSED16 = B(34);
         s(f).UNUSED17 = B(35);
         s(f).LEVEN = logical(B(36));
         if ~isnan(B(37))
            s(f).LPSPOL = logical(B(37));
         else
            s(f).LPSPOL = B(37);
         end
         if ~isnan(B(38))
            s(f).LOVROK = logical(B(38));
         else
            s(f).LOVROK = B(38);
         end
         if ~isnan(B(39))
            s(f).LCALDA = logical(B(39));
         else
            s(f).LCALDA = B(39);
         end
         if ~isnan(B(40))
            s(f).UNUSED18 = logical(B(40));
         else
            s(f).UNUSED18 = B(40);
         end
         s(f).KSTNM = C{1};
         s(f).KEVNM = [C{2:3}];
         s(f).KHOLE = C{4};
         s(f).KO = C{5};
         s(f).KA = C{6};
         s(f).KT0 = C{7};
         s(f).KT1 = C{8};
         s(f).KT2 = C{9};
         s(f).KT3 = C{10};
         s(f).KT4 = C{11};
         s(f).KT5 = C{12};
         s(f).KT6 = C{13};
         s(f).KT7 = C{14};
         s(f).KT8 = C{15};
         s(f).KT9 = C{16};
         s(f).KF = C{17};
         s(f).KUSER0 = C{18};
         s(f).KUSER1 = C{19};
         s(f).KUSER2 = C{20};
         s(f).KCMPNM = C{21};
         s(f).KNETWK = C{22};
         s(f).KDATRD = C{23};
         s(f).KINST = C{24};

         % read data blocks
         clear count2
         [D, count1] = fread(fid, [s(f).NPTS, 1], 'float32');
         s(f).DATA1 = D;
         if ~s(f).LEVEN || strcmp(s(f).IFTYPE, 'IAMPH') || ...
             strcmp(s(f).IFTYPE, 'IRLIM')
            [D, count2] = fread(fid, [s(f).NPTS, 1], 'float32');
            s(f).DATA2 = D;
         end

         % flag successfully read SAC file
         if count1 == s(f).NPTS
            if exist('count2','var')
               if count2 == s(f).NPTS
                   count(f) = true;
               end
            else
               count(f) = true;
            end
         end

      end

      % close SAC file.
      fclose(fid);
   end
end      
s = s(count);

% assign output
if nargout == 1
   t1 = s;
elseif nargout == 2
   [t1, t2] = getsacdata(s);
else
   [t1, t2, d] = getsacdata(s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = sacstruct(n)
% S = SACSTRUCT(N); returns an [Nx1] SAC structure with all fields
% initiated to their undefined (default) values.. N defaults to 1. 
%

%	Xiaoning Yang 2008 (modified from subfunction init_sac of readsac.m
% in the MATSEIS package)

% check arguments.
if (nargin < 1)
   n = 1;
end

% initialize output structure.
sacfields = {'FILENAME'; 'DELTA'; 'DEPMIN'; 'DEPMAX'; 'SCALE'; 'ODELTA'; ...
   'B'; 'E'; 'O'; 'A'; 'INTERNAL1'; 'T0'; 'T1'; 'T2'; 'T3'; 'T4'; 'T5'; ...
   'T6'; 'T7'; 'T8'; 'T9'; 'F'; 'RESP0'; 'RESP1'; 'RESP2'; 'RESP3'; ...
   'RESP4'; 'RESP5'; 'RESP6'; 'RESP7'; 'RESP8'; 'RESP9'; 'STLA'; 'STLO'; ...
   'STEL'; 'STDP'; 'EVLA'; 'EVLO'; 'EVEL'; 'EVDP'; 'MAG'; 'USER0'; ...
   'USER1'; 'USER2'; 'USER3'; 'USER4'; 'USER5'; 'USER6'; 'USER7'; 'USER8'; ...
   'USER9'; 'DIST'; 'AZ'; 'BAZ'; 'GCARC'; 'INTERNAL2'; 'INTERNAL3'; ...
   'DEPMEN'; 'CMPAZ'; 'CMPINC'; 'XMINIMUM'; 'XMAXIMUM'; 'YMINIMUM'; ...
   'YMAXIMUM'; 'UNUSED1'; 'UNUSED2'; 'UNUSED3'; 'UNUSED4'; 'UNUSED5'; ...
   'UNUSED6'; 'UNUSED7'; 'NZYEAR'; 'NZJDAY'; 'NZHOUR'; 'NZMIN'; 'NZSEC'; ...
   'NZMSEC'; 'NVHDR'; 'NORID'; 'NEVID'; 'NPTS'; 'INTERNAL4'; 'NWFID'; ...
   'NXSIZE'; 'NYSIZE'; 'UNUSED8'; 'IFTYPE'; 'IDEP'; 'IZTYPE'; 'UNUSED9'; ...
   'IINST'; 'ISTREG'; 'IEVREG'; 'IEVTYP'; 'IQUAL'; 'ISYNTH'; 'IMAGTYP'; ...
   'IMAGSRC'; 'UNUSED10'; 'UNUSED11'; 'UNUSED12'; 'UNUSED13'; 'UNUSED14'; ...
   'UNUSED15'; 'UNUSED16'; 'UNUSED17'; 'LEVEN'; 'LPSPOL'; 'LOVROK'; ...
   'LCALDA'; 'UNUSED18'; 'KSTNM'; 'KEVNM'; 'KHOLE'; 'KO'; 'KA'; 'KT0'; ...
   'KT1'; 'KT2'; 'KT3'; 'KT4'; 'KT5'; 'KT6'; 'KT7'; 'KT8'; 'KT9'; 'KF'; ...
   'KUSER0'; 'KUSER1'; 'KUSER2'; 'KCMPNM'; 'KNETWK'; 'KDATRD'; ...
   'KINST'; 'DATA1';};

cl = cell(size(sacfields,1),n);
cl(2:111, :) = {nan};
cl(112:end-1, :) = {' '};
s = cell2struct(cl, sacfields, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out1, out2, out3] = getsacdata(s)
% [INDEP, DEP] = GETSACDATA(S); extracts independent and dependent variables
% from SAC structure S and puts them in INDEP and DEP. INDEP and DEP can either
% be time and data or x and y. They are cell arrays if length(S) > 1 and 
% the lengths of the variables for different S elements are different, 
% or double vectors (matrices) if either length(S) == 1 or the lengths of
% the variables for different S elements are the same.
%
% [OUT1, OUT2, OUT3] = GETSACDATA(S); extracts either spectral data from spectral
% type of SAC files or xyz data from xyz (spectrogram) SAC files. In case 
% of spectral files, OUT1 is frequency vector, OUT2 is either real part or
% amplitude, and OUT3 is either imaginary part or phase. It is assumed, as
% the SAC default, that the complete complex spectrum of a signal, which 
% spans from zero to (NPTS-1)*df where df is frequence sampling rate, are
% present. OUT1 only contains frequency from zero to Fn, where Fn is the
% Nyquist frequency (usually NPTS/2*df.) For xyz data, OUT1 is x (time),
% OUT2 is y (frequency) and OUT3 is z (spectrogram).
%
% *** As of this version, SAC does not define general xyz type files. ***
%
% S is defined in m-file readsac.m and is usually the output from readsac.m.
% S must contain same type of data (time series, spectra, xy or xyz).

% check output arguments
if nargout < 2
   error(' Number of output arguments must be larger than one !!!')
end

% allocate output arguments
l = length(s);
if l == 0
   out1 = [];
   out2 = [];
   if nargout == 3
      out3 = [];
   end
   return
end
out1 = cell(l,1);
out2 = cell(l,1);
if nargout == 3
   out3 = cell(l,1);
end


% loop over s
ll = zeros(l, 1); % lengths of time series
for i = 1:l
   ll(i) = length(s(i).DATA1);
   if i == 1
      file_type = s(i).IFTYPE;
      
      % reference time for all traces
      if strcmp(file_type, 'ITIME')
         year0 = s(i).NZYEAR-1901;
         days0 = year0*365+round(year0/4)+s(i).NZJDAY;
         ref0 = days0*86400+s(i).NZHOUR*3600+s(i).NZMIN*60+s(i).NZSEC+s(i).NZMSEC/1000;
         if isnan(ref0)
            ref0 = 0;
         end
      end
        
   end

   % check file type
   if ~strcmp(s(i).IFTYPE, file_type)
      error(' File types in input structure must be same !!!')
   end

   if strcmp(file_type, 'IXYZ')   
      if nargout ~= 3
         error(' Number of output arguments for XYZ file must be three !!!')
      end
      out1(i) = {linspace(s(i).XMINIMUM, s(i).XMAXIMUM, s(i).NXSIZE)};
      out2(i) = {linspace(s(i).YMINIMUM, s(i).YMAXIMUM, s(i).NYSIZE)'};
      out3(i) = {reshape(s(i).DATA1, s(i).NXSIZE, s(i).NYSIZE)'};
   elseif strcmp(file_type, 'IXY') && ~s(i).LEVEN
      if nargout ~= 2
         error(' Number of output arguments for XY file must be two !!!')
      end
      out1(i) = {s(i).DATA2};   % see SAC Manual
      out2(i) = {s(i).DATA1};
   else     % spectrum and evenly-spaced time series         
      if isnan(s(i).B)
         B = 0;
      else
         B = s(i).B;
      end
      if strcmp(file_type, 'ITIME')

         if isnan(s(i).O)
            % reference time of this trace
            year = s(i).NZYEAR-1901;
            days = year*365+round(year/4)+s(i).NZJDAY;
            ref = days*86400+s(i).NZHOUR*3600+s(i).NZMIN*60+s(i).NZSEC+ ...
               s(i).NZMSEC/1000;
            if isnan(ref)
                ref = 0;
            end

            B = B+ref-ref0;
            if isnan(B)
               B = 0;
            end
            O = 0;
         else
            O = s(i).O;
         end
         out1(i) = {(0:s(i).NPTS-1)'*s(i).DELTA+B-O};
      elseif strcmp(file_type, 'IXY')
         out1(i) = {(0:s(i).NPTS-1)'*s(i).DELTA+B};
      else  % frequency vector for positive-frequency part of the spectrum
          if B ~= 0
              warning(' First frequency point is not zero Hz !!!')
          end
          npts = round(s(i).NPTS/2);
          if s(i).NPTS/2 ~= npts
              warning(' Spectrum has odd number of data points !!!')
          end
         out1(i) = {(0:npts)'*s(i).DELTA+B};
      end

      if (strcmp(file_type, 'ITIME') || strcmp(file_type, 'IXY')) && nargout ~= 2
         error([' Number of output arguments for time-series data ', ...
            'must be two !!!'])
      elseif (strcmp(file_type, 'IAMPH') || strcmp(file_type, 'IRLIM')) && nargout ~= 3
         error([' Number of output arguments for spectral data must be ', ...
            'three !!!'])
      else
         out2(i) = {s(i).DATA1};
         if strcmp(file_type, 'IAMPH') || strcmp(file_type, 'IRLIM')
            out3(i) = {s(i).DATA2};
         end
      end
   end
end

if l == 1
   out1 = out1{:};
   out2 = out2{:};
   if nargout == 3
      out3 = out3{:};
   end
elseif ~any(diff(ll))
   temp1 = zeros(ll(1), l);
   temp2 = temp1;
   if nargout == 3
      temp3 = temp1;
   end
   if strcmp(file_type, 'IRLIM') || strcmp(file_type, 'IAMPH')
       temp1 = zeros(round(ll(1)/2)+1, 1);
   end
   for i = 1:l
      temp1(:, i) = out1{i};
      temp2(:, i) = out2{i};
      if nargout == 3
         temp3(:, i) = out3{i};
      end
   end
   out1 = temp1;
   out2 = temp2;
   if nargout == 3
      out3 = temp3;
   end
end
