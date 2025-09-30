function [wavelength,content] = read_ANDO(filename)
%READ_ANDO It reads ANDO's spectrum csv files.

fid = fopen(filename);

data = textscan(fid,'%f %f',1001,'Delimiter',',','Headerlines',3,'EndOfLine','\r\n');

wavelength = data{1};
content = data{2}; % this is the spectrum

fclose(fid);

end