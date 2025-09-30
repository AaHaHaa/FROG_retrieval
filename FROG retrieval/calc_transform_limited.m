function [transform_limited_field,varargout] = calc_transform_limited( app,input_field,num_insert_points,t )
%DECHIRP_TRANSFORM_LIMITED It gives the dechirped result which is transform-limited.
%
% Input:
%   input_field: (Nt,...); the electric field in time domain (sqrt(W))
%   num_insert_points: (optional input argument)
%                      a scalar;
%                      the number of points inserted between two original data points
%                      (default: zero)
%
%   **Because the transform-limited pulse is really narrow, 
%     the original time spacing can be too large to correctly characterize the pulse.
%     This will be a problem with calculating the transform-limited duration.
%     Therefore, increasing the number of time grid points, i.e., reduing the time spacing, is highly recommended.
%
%   t: (optional input argument)
%      (Nt,1); time (ps)
%
% Output:
%   transform_limited_field
%   t_insert: the interpolated time grid points based on "num_insert_points"
%   transform_limited_FWHM: transform-limited pulse duration (fs)
%   pulse_FWHM: current pulse duration (fs)
% =========================================================================
% "num_insert_points" and "t" are optional input arguments, but "t" is
% required to calculate "t_insert" and "transform_limited_FWHM" as output.
% =========================================================================
% Usage:
%   transform_limited_field = calc_transform_limited(input_field);
%   transform_limited_field = calc_transform_limited(input_field,num_insert_points);
%   [transform_limited_field,t_insert,transform_limited_FWHM] = calc_transform_limited(input_field,num_insert_points,t);
%
% Note that this code employs the mathematical Fourier-transform convention.

calc_TL_duration = false;
switch nargin
    case 2
        num_insert_points = 0;
    case 4
        calc_TL_duration = true;
end

Nt = size(input_field,1);
insert_idx = [linspace(1,Nt,(num_insert_points+1)*(Nt-1)+1)'; Nt+(1:num_insert_points)'/(num_insert_points+1)];

% Interpolation
input_freq = fft(input_field);
%background_noise = feval(@(x) mean(abs(x(1:100)).^2)*2,fftshift(input_freq,1));
%input_freq(abs(input_freq).^2<background_noise) = 0; % remove background noise
if num_insert_points ~= 0
    input_freq = cat(1,input_freq(1:ceil(Nt/2),:),zeros(Nt*num_insert_points,1),input_freq(ceil(Nt/2)+1:end,:));
    input_field = reshape(ifft(input_freq),[Nt*(num_insert_points+1),1])*(num_insert_points+1);
else
    input_field = ifft(input_freq);
end

input_freq_TL = abs(fftshift(fft(ifftshift(input_field,1)),1));
transform_limited_field = fftshift(ifft(ifftshift(input_freq_TL,1)),1);

threshold = max(abs(input_field).^2)/1.5;
all_peaks = findpeaks(abs(transform_limited_field).^2,'MinPeakHeight',threshold,'WidthReference','halfheight');
if ~isempty(app) && length(all_peaks) ~= 1
    message = {'Please increase either the time or spectral window (most likely, the time window),',...
               'and re-run ''Resample'' and ''Retrieval'' to avoid aliasing in transform-limited computations.'};
    uialert(app.UIFigure,message,'Warning','Icon','warning');
end

if calc_TL_duration
    % Inserted time
    t_insert = interp1(t,insert_idx,'linear','extrap');
    
    % Current duration
    threshold = max(abs(input_field).^2)/1.01;
    [~,~,tmp_pulse_width,~] = findpeaks(abs(input_field).^2,t_insert*1e3,'MinPeakHeight',threshold,'WidthReference','halfheight');
    pulse_FWHM = tmp_pulse_width(1);
    
    % Transform-limited duration
    threshold = max(abs(transform_limited_field).^2)/1.0001;
    [~,~,transform_limited_FWHM,~] = findpeaks(abs(transform_limited_field).^2,t_insert*1e3,'MinPeakHeight',threshold,'WidthReference','halfheight');
    
    varargout = {t_insert,transform_limited_FWHM,pulse_FWHM};
else
    varargout = {};
end

end