function [ThrSlider_value,clean_FROGtrace] = find_optimal_threshold(FROGtrace,thr)
%FIND_OPTIMAL_THRESHOLD

if thr ~= 0
    total_pixels = numel(FROGtrace);
    
    i = 1;
    center_difference = zeros(991,1);
    noise_pixels = zeros(991,1);
    for ThrSlider_value = 1:0.1:100
        threshold = ThrSlider_value/100*thr*10;
        
        clean_FROGtrace = FROGtrace;
        below_threshold = FROGtrace<threshold;
        clean_FROGtrace(below_threshold) = 0;
        
        [~, ~, ~, x, y] = calcMFD(clean_FROGtrace, threshold);
        [~, indx] = min(abs(x));
        [~, indy] = min(abs(y));
        
        noise_pixels(i) = sum(below_threshold(:));
        if i > 1
            center_difference(i) = norm( [indx,indy] - [indx0,indy0] );
            
            pixel_diff = noise_pixels(i) - noise_pixels(i-1);
            if center_difference(i) < min(size(FROGtrace))/100 && ...
                    (noise_pixels(i) > total_pixels*0.2 && pixel_diff < total_pixels*0.01)
                break;
            end
        end
        indx0 = indx;
        indy0 = indy;
        
        i = i+1;
    end
    
    % Since only the central part is considered, I tend to clear all pixels
    % with small noise that is still slightly above the previously calculated
    % threshold. A new threshold is found here.
    [D4sigmaX, D4sigmaY, clean_FROGtrace, x, y] = calcMFD(clean_FROGtrace, threshold);
    noise_map = true(size(clean_FROGtrace));
    noise_map(abs(x)<D4sigmaX/2*1.5 & abs(y)<D4sigmaY/2*1.5) = false;
    noiseLevel_FROGtrace = clean_FROGtrace(noise_map);
    if ~isempty(noiseLevel_FROGtrace)
        threshold = max(noiseLevel_FROGtrace);
    end
    clean_FROGtrace(clean_FROGtrace<threshold) = 0;
    ThrSlider_value = max(1,min(100,ceil(threshold*100/thr/10)));
else % thr = 0 happens because the input FROGtrace has been cleaned already, which occurs when loading Mesaphotonics's exported ".txt" data, instead of its generated ".raw" data.
    ThrSlider_value = 100; % arbitrarily assign a value here
    clean_FROGtrace = FROGtrace; % FROG trace has been cleaned already
end

end

