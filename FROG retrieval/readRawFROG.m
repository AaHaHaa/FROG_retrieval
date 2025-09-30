function [delay, wavelength, FROGtrace] = readRawFROG(filename, device)
%READRAWFROG This reads the raw FROG trace.
%
%   filename: a string of filename
%   device: a string;
%           It currently supports reading raw files only from Mesaphotonics,
%           so only 'Mesaphotonics' is accepted.
%
% =========================================================================
% FROG raw files from Mesaphotonics includes
%   lines 1-4: basic information about this file, like
%       line 1: the number of lines for the wavelength data
%       line 2: the number of lines for the raw FROG trace data
%       line 3: 
%       line 4: 
%   lines 5-end: the raw FROG trace data
%
%   The raw FROG data needs to be processed as the following:
%       num_wavelength = line 1 number
%       num_delay = line 2 number / line 1 number
%       FROGtrace = reshape(FROGtrace, num_delay, num_wavelength)
%       FROGtrace = rot90(sqrt(FROGtrace),3)
%
% -------------------------------------------------------------------------
%
% Mesaphotonics VideoFROG software allows a user to export its retrieved
% result and "resampled" FROG trace, etc.
%
% Below are samples of its exported files:
%     FROGtrace 2_9_time_9_15_2575_0.txt
%     RetrievedFT 2_9_time_9_15_2575_0.txt
%     Pulse_Complex_Time_Domain 2_9_time_9_15_2575_0.txt
%     Pulse_Time_Domain 2_9_time_9_15_2575_0.txt
%     Pulse_Freq_Domain 2_9_time_9_15_2575_0.txt
%     Pulse_Stats 2_9_time_9_15_2575_0.txt
%     Gate 2_9_time_9_15_2575_0.txt
%
% FROGtrace: it includes the measured FROG trace, resampled for retrieval.
%            VideoFROG's export function doesn't export its raw FROG trace.
%            Raw FROG trace needs to be exported with another "save" button.
% RetrievedFT: it includes the retrieved FROG trace
% Pulse_Complex_Time_Domain: it includes the tabular data following [t, real(E), imag(E)]
% Pulse_Time_Domain: it includes the tabular data following [t, abs(Et).^2, angle(Et)]
% Pulse_Freq_Domain: it includes the tabular data following [f, lambda, abs(Ef).^2, angle(Ef)]
%                    "f" is the resampled frequency during retrieval, so it's centered at 0. However, lambda corresponds to the real wavelength.
% Pulse_Stats: it includes the numeric values, such as retrieved pulse duration
% Gate: (unknown)
%
% VideoFROG has resampled "t" the same as "delay" in the FROG trace.
% However, in principle, they are different. This VideoFROG's restriction
% self-limits its freedom to arbitrarily sample the field. For example,
% VideoFROG has a problem with retrieving a broadband chirped pulse that
% extends both in time and frequency. My FROG retrieval code releases this
% contraint.

switch device
    case 'Mesaphotonics (raw)'
        fid = fopen(filename,'r','ieee-be');
        
        % Basic information of this file
        specifications = fread(fid,4,'float32');
        num_wavelength = specifications(1);
        num_FROGtrace = specifications(2);
        num_delay = num_FROGtrace/num_wavelength;
        % Read data
        wavelength = fread(fid,num_wavelength,'float32');
        FROGtrace = fread(fid,num_FROGtrace,'float32');
        fclose(fid);
        
        % Remove outliers in FROGtrace, which is a problem after we install
        % the broadband OceanOptics spectrometer bought by Yi-Hao in 2023.
        % Use MATLAB's rmoutliers, which removes outliers based on local moving
        % medians. If the local window size is too large, the median will be
        % zero due to pulse existing in only a small temporal delay window.
        % This causes all signals, including outlier noise, to be removed.
        % Therefore, I find the best local window size by finding where the 
        % summation value stops to change significantly due to removal of most
        % outliers.
        min_window_size = 3;
        window_size_idx = 1:10;
        FROGtrace_sum = zeros(1,length(window_size_idx));
        for window_size_i = window_size_idx
            FROGtrace_sum(window_size_i) = sum(rmoutliers(FROGtrace,'movmedian',(window_size_i-1)+min_window_size));
        end
        FROGthreshold = mean(FROGtrace_sum(end-4:end))*2;
        idx = find(FROGtrace_sum>FROGthreshold,1,'last');
        if ~isempty(idx)
            window_size = (window_size_idx(idx)-1) + min_window_size;
            [~,rm_idx] = rmoutliers(FROGtrace,'movmedian',(window_size-1)*2+min_window_size);
            FROGtrace(rm_idx) = 0;
        end
        
        % FROG trace
        FROGtrace = reshape(FROGtrace,num_delay,num_wavelength);
        FROGtrace(FROGtrace<=0) = 0;
        FROGtrace = rot90( sqrt(FROGtrace), 3); % make it amplitude and correct the orientation
        % delay
        delay_step_size = specifications(3)*specifications(4)*1e-15; % s; delay step size
        delay = (-size(FROGtrace,2)/2*delay_step_size:delay_step_size:(size(FROGtrace,2)/2-1)*delay_step_size); % s
        delay_center = sum(delay.*sum(FROGtrace,1))/sum(sum(FROGtrace,1));
        delay = delay - delay_center;
    case 'Mesaphotonics (txt)'
        try
            % Load the folder
            if ispc
                sep_char = '\';
            else % unix
                sep_char = '/';
            end
            sep_pos = strfind(filename,sep_char);
            foldername = filename(1:sep_pos(end));

            % Mesaphotonics saves their data with different labels after
            % the whitespace, such as
            %   FROGtrace 2_9_time_9_15_2575_0.txt
            %   Gate 2_9_time_9_15_2575_0.txt
            %   Pulse_Complex_Time_Domain 2_9_time_9_15_2575_0.txt
            %   Pulse_Freq_Domain 2_9_time_9_15_2575_0.txt
            %   Pulse_Stats 2_9_time_9_15_2575_0.txt
            %   Pulse_Time_Domain 2_9_time_9_15_2575_0.txt
            %   RetrievedFT 2_9_time_9_15_2575_0.txt
            file_version = strsplit(filename,' ');
            file_version = file_version{end};

            dir_all = dir(foldername);
            for i = 1:length(dir_all)
                if ~dir_all(i).isdir
                    filetype = strsplit(dir_all(i).name,'.');
                    filetype = filetype{end};
                    if isequal(filetype,'txt') % read only 'txt'
                        filetype = strsplit(dir_all(i).name,' ');
                        if isequal(filetype{2},file_version) % read only the specified version
                            filetype = filetype{1};
                            if ~isequal(filetype,'Pulse_Stats')
                                data = load([foldername '/' dir_all(i).name]);
                            end
                            switch filetype
                                case 'FROGtrace'
                                    FROGtrace = flipud(data); % measured FROG trace
                                case 'Pulse_Complex_Time_Domain'
                                    delay = data(:,1)'*1e-15; % s (VideoFROG's output is in "fs")
                                case 'Pulse_Freq_Domain'
                                    wavelength = flipud(data(:,2))/2; % nm
                            end
                        end
                    end
                end
            end
            if ~exist('FROGtrace','var') || ~exist('delay','var') || ~exist('wavelength','var') % if there is any missing VideoFROG's file
                error('readRawFROG:FROGFileError','...'); % just to get into the "catch" block below
            end
        catch
            error('readRawFROG:FROGFileError',...
                  'The folder should contain all txt files generated by the VideoFROG.');
        end
end

end