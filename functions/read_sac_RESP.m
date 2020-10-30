function [zz,pp,const,units] = read_sac_RESP(fnameRESP,units)
%  [zz,pp,constant,unit] = read_sac_RESP(respn, traces)
%
% Read poles, zeros, and sensitivity from IRIS RESP file
% (https://ds.iris.edu/ds/nodes/dmc/data/formats/resp/)
% By default, RESP file is in units of velocity (m/s). Can add an extra 
% 0+i0 zero to convert to units of displacement (m).
%
% INPUT
% fnameRESP: input RESP filename
% units: 'M' meters | 'M/S' meters per second
% 
% OUTPUT
% zz: zeros
% pp: poles
% const: (final sensitivity) * (A_0)
% units: 'M' or 'M/S'
%
% Josh Russell 20/2
% updated 20/10

    % Units for output poles and zeros
    % units = 'M'; % M (displacement) | M/S (velocity)
    
    % read file and Split on newlines
    str = fileread(fnameRESP);
    resp = regexp(str,'\r\n|\r|\n','split');
    resp = splitlines(resp);
    
    %initialize some flags and various variables
    zn = 0;
    pn = 0;
    zero_scan_flag = 1;
    pole_scan_flag = 1;
    A0_scan_flag = 1;
    z_re=[]; z_im=[]; z_ind=[]; iz=0;
    p_re=[]; p_im=[]; p_ind=[]; ip=0;
    sensitivity=[]; isensitivity=0;
    
    % Define RESP line codes
    A0_id = 'B053F07';
    zn_id = 'B053F09';
    pn_id = 'B053F14';
    z_id = 'B053F10-13';
    p_id = 'B053F15-18';
    sensitivity_id = 'B058F04';
    
    % loop over the entire file
    for irow = 1:length(resp)
        buff = resp{irow};
        if length(buff) < 10
            continue
        end
        
        % Get line id, killing trailing white spaces
        id = deblank(buff(1:10));
        
        % A0
        if strcmp(id,A0_id) && A0_scan_flag
            scan = textscan(buff,'%s A0 normalization factor: %f');
            A0 = scan{2};
            A0_scan_flag = 0;
            continue
        end
        
        % Zeros
        if zero_scan_flag
            if strcmp(id,zn_id)
                scan = textscan(buff,'%s Number of zeroes: %f');
                zn = scan{2};
            end
            if strcmp(id,z_id)
                iz = iz+1;
                scan = textscan(buff,'%s %d %f %f %f %f');
                z_ind(iz,:) = scan{2};
                z_re(iz,:) = scan{3};
                z_im(iz,:) = scan{4};
                if z_ind(iz) == zn-1
                    zero_scan_flag = 0;
                end
                continue
            end
        end
        
        % Poles
        if pole_scan_flag
            if strcmp(id,pn_id)
                scan = textscan(buff,'%s Number of poles: %f');
                pn = scan{2};
            end
            if strcmp(id,p_id)
                ip = ip+1;
                scan = textscan(buff,'%s %d %f %f %f %f');
                p_ind(ip,:) = scan{2};
                p_re(ip,:) = scan{3};
                p_im(ip,:) = scan{4};
                if p_ind(ip) == pn-1
                    pole_scan_flag = 0;
                end
                continue
            end
        end
        
        % Sensitivity
        if strcmp(id,sensitivity_id)
            isensitivity = isensitivity+1;
            scan = textscan(buff,'%s Sensitivity: %f');
            sensitivity(isensitivity,:) = scan{2};
            continue
        end
    end
    zz = z_re + 1i*z_im;
    pp = p_re + 1i*p_im;
    
    % Add an extra zero for units in displacement
    if strcmp(upper(units),'M')
        zz = [0 + 1i*0; zz];
    end
    
    % Check sensitivity
    final_sens = sensitivity(end);
    sens_tst = 1;
    for ii = 1:length(sensitivity)-1
        sens_tst = sens_tst*sensitivity(ii);
    end
%     if abs(sens_tst-final_sens)./final_sens*100>1e-1 && length(traces)==3
%         error('Instrument sensitivities don''t add up... check RESP');
%     end

    % Calculate constant gain
    const = final_sens * A0;
        
   
    return
