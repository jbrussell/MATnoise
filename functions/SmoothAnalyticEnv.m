function [ sm_analytic_env ] = SmoothAnalyticEnv( x,data )
% Find analytic envelope of bessel function data. Fits a polynomial to the
% analytic envelope in order to get rid of spikes at the edges.
%
% Josh Russell
% github.com/jbrussell

analytic_env = abs(hilbert(data));
p = polyfit(x,analytic_env,2);
sm_analytic_env = polyval(p,x);

end

