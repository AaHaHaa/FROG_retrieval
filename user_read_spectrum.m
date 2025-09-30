function [wavelength,spectrum] = user_read_spectrum(filename)
%USER_READ_SPECTRUM This is the user-defined function to read the spectrum
%for the FROG retrieval code.
%
% There are cases when the spectrum file doesn't match the default of the
% FROG retrieval's GUI, so this function gives a user a freedom to design
% the function to read their own spectrum.
%
% Note that 
%   filename: a string specifying where the spectrum file is
%   wavelength: a column vector in "nm"
%   spectrum: a column vector in the "linear" scale.

% This is a sample file that reads the spectrum from a MATLAB's "mat" file.
% The data is saved in the "D" structure with D.wl the wavelength and
% D.spec the spectrum in "dB".
load(filename);

wavelength = D.wl'; % nm
spectrum = 10.^(D.spec'/10); % if D.spec is in log scale (dB); if it's linear scale, just do "spectrum = D.spec;"

end