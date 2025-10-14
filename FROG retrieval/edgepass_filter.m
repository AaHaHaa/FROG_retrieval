function filtered_FROGtrace = edgepass_filter(FROGtrace, x, cutonoff)
%EDGEPASS_SPECTRAL_FILTER Apply a edgefilter to a field
%
% Input:
%   FROGtrace - a (Nf, Nt) matrix
%   x - x or y axis data
%   cutonoff - a (1,2) array; the cutoff/cuton x of the edgefilter

cutonoff_slope = 3; % the slope at the cutoff/cuton point
gaussexpo = 1; % supergaussian exponent (~exp(-f^(2*gaussexpo)))

% Calculate the filter profile in frequency space
x_0 = (2*gaussexpo-1)/cutonoff_slope*nthroot(gaussexpo/(2*gaussexpo-1),2*gaussexpo)*exp(-(2*gaussexpo-1)/(2*gaussexpo));
center_x1 = cutonoff(1) - x_0*nthroot(log(2),2*gaussexpo); % highpass
center_x2 = cutonoff(2) + x_0*nthroot(log(2),2*gaussexpo); % lowpass
mult_factor1 = exp(-(x-center_x1).^(2*gaussexpo)/(2*x_0^(2*gaussexpo)));
mult_factor2 = exp(-(x-center_x2).^(2*gaussexpo)/(2*x_0^(2*gaussexpo)));
mult_factor1(x>=center_x1) = 1; % highpass
mult_factor2(x<=center_x2) = 1; % lowpass
mult_factor = mult_factor1.*mult_factor2;

% Apply the filter
filtered_FROGtrace = FROGtrace.*mult_factor;

end

