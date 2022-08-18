function [dtds_int, ier]=ffrq_2dkernel_lalo_fb(kmode,nfreqs,bw_frac,freq,c0,X1,X2, xg, yg,varargin)
%function [dtds, ier]=ffrq_2dkernel_lalo_fb(kmode,nfreqs,bw_frac,freq,c0,X1,X2, xg, yg,varargin)
%
%  Loop over ffrq_2dkernel using gaussian average to obtain "finite
%  bandwidth" estimate. This is less harsh than truncating the kernel at a
%  specified Fresnel zone.
%
% kmode: =0 Single frequency, all Fresnel zones
%      : >0 Single frequency calculated up to the kth Fresnel zone
%      : <0 Averge over bandwidth
%
% nfreqs: number of frequencies to include in the gaussian function
% bw_frac: bandwidth defined by fraction of central freq [0-1]
%
% Josh Russell 4/2020

type = 'analytical';
if ~isempty(varargin)
    s1 = varargin{1};
    s2 = varargin{2};
    type = 'empirical';
end

dtds_int = 0;

if kmode>=0 % single frequency truncated at kth Fresnel zone
    switch type
        case 'analytical'
            [dtds_int, ier]=ffrq_2dkernel_lalo(kmode,freq, c0,X1,X2, xg, yg);
        case 'empirical'
            [dtds_int, ier]=ffrq_2dkernel_lalo(kmode,freq, c0,X1,X2, xg, yg,s1,s2);
    end
    
% USING OUR QUICK & DIRTY METHOD. SPATIALLY RESTRICTS KERNELS MORE THAN LIN & RITZWOLLER (2010) eq. 14
else % Average over a bandwidth
    fvec = freq*linspace(1-bw_frac,1+bw_frac,nfreqs);
    fweight = gausswin(nfreqs)/sum(gausswin(nfreqs));
    switch type
        case 'analytical'
            for ii = 1:length(fvec)
                [dtds, ier]=ffrq_2dkernel_lalo(0,fvec(ii), c0,X1,X2, xg, yg);
                dtds_int = dtds_int + dtds.*fweight(ii);
            end
        case 'empirical'
            for ii = 1:length(fvec)
                [dtds, ier]=ffrq_2dkernel_lalo(0,fvec(ii), c0,X1,X2, xg, yg,s1,s2);
                dtds_int = dtds_int + dtds.*fweight(ii);
            end
    end
end

% % USING METHOD FROM LIN & RITZWOLLER 2010 (eq. 14)
% else % Average over a bandwidth
%     dtds_int = [];
%     fvec = freq*linspace(1-bw_frac,1+bw_frac,nfreqs);
% %     fweight = gausswin(nfreqs)/sum(gausswin(nfreqs));
%     gausfun = gausswin(nfreqs);
%     switch type
%         case 'analytical'
%             for ii = 1:length(fvec)
%                 [dtds, ier]=ffrq_2dkernel_lalo(0,fvec(ii), c0,X1,X2, xg, yg);
%                 dtds_int = [dtds_int, dtds*gausfun(ii).^2];
%             end
%         case 'empirical'
%             for ii = 1:length(fvec)
%                 [dtds, ier]=ffrq_2dkernel_lalo(0,fvec(ii), c0,X1,X2, xg, yg,s1,s2);
%                 dtds_int = [dtds_int, dtds*gausfun(ii).^2];
%             end
%     end
%     numer = trapz(fvec,dtds_int,2);
%     denom = trapz(fvec,gausfun.^2);
%     dtds_int = numer./denom;
% end

return
