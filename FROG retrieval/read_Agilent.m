function [wavelength,content] = read_Agilent( filename )
%READ_AGILENT It reads Agilent's spectrum csv files.

fid = fopen(filename);

% wavelength
wavelength_range = textscan(fid,'%s %f',2,'Delimiter',',','Headerlines',7,'EndOfLine','\r\n');
wavelength_start = wavelength_range{2}(1);
wavelength_end   = wavelength_range{2}(2);
wavelength_points = textscan(fid,'%s %f',1,'Delimiter',',','Headerlines',2,'EndOfLine','\r\n');
wavelength_points = wavelength_points{2};
wavelength = linspace(wavelength_start,wavelength_end,wavelength_points)';

% intensity
content = textscan(fid,'%f %f','Delimiter',',','Headerlines',16,'EndOfLine','\r\n');
content = content{2}; % this is the spectrum

fclose(fid);

end

