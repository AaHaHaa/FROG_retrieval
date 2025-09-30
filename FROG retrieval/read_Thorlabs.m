function [wavelength,content] = read_Thorlabs(filename)
%READ_THORLABS It reads Thorlabs's spectrum csv files.

fid = fopen(filename);

% Find the length of the data
read_string = {{''}};
while ~contains(read_string{1}{1},'#Length')
    read_string = textscan(fid,'%s',1);
end
wavelength_points = textscan(read_string{1}{1},'#Length %f','Delimiter',';');
wavelength_points = round(wavelength_points{1});

% The data is after the [Data] tag, so we read until find this tag, and
% then we read the data
while ~isequal(read_string{1}{1},'[Data]')
    read_string = textscan(fid,'%s',1);
end
data = textscan(fid,'%f %f',wavelength_points,'Delimiter',';','EndOfLine','\r\n');

wavelength = data{1};
content = data{2}; % this is the spectrum

fclose(fid);

end