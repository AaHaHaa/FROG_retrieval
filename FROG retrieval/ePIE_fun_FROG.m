function [bestObj, FROGtrace_error, Ir] = ePIE_fun_FROG(I0, delay, iterMax, freq, time, stopTol, guess, DispFlag, gpu_yes, UIFigure, use_measured_spectrum, ratio_of_measured_spectrum)
%EPIE_FUN_FROG ePIE function for SHG FROG
%It reconstructs a pulse function (in time) from a SHG FROG trace by use of
%the Ptychographic algorithm.
%
% Input:
%   I0      =   (Nf,Nt); Experimental / Simulated SHG FROG Trace
%   delay   =	(1 ,Nd); a vector of delays that corresponds to trace.
%   iterMax = Maximum number of iterations allowed (default = 1000).
%   freq    = (Nf,1); a vector of frequencies.
%   time    = (Nt,1); the coresponding vector in time domain.
%   stopTol = (Optioanl) Tolerence on the error (default: 1e-5).
%   guess   = initial guess/pulse for the retrieval algorithm
%   gpu_yes = whether to use GPU or not (default: false)
%
% Output:
%   bestObj = the objective function; the final reconstructed pulse field (in time)
%   error   = a vector of errors for each iteration
%   Ir      = the reconstructed FROG trace.
%
% Variables used in the code:
%   Obj = objective function; the reconstructed pulse field (in time)

dt = mean(diff(time))*1e12; % ps
freq_plot = freq*1e-12;

% ifftshift for fft or ifft
I0 = ifftshift(I0,1);
freq = ifftshift(freq,1);

I0 = ( I0+fliplr(I0) )/2; % it needs to be symmetric, especially numerically
I0 = I0/sqrt( sum(sum((I0).^2)) ); % normalization

if ~use_measured_spectrum
    I = ifftshift(imgaussfilt(fftshift(I0,1),5),1);
else
    I = I0;
end

 if DispFlag
     % Display the retrieved pulse at each iteration and also give a stop
     % button for stopping the retrieval process
     fig = figure('rend','painters','pos',[10 10 900 600]);      
     ann = annotation(gcf,'textbox',...
         [0.285 0.97 0.351 0.036],'String',...
         {'    ccc       '},'FitBoxToText','on', 'LineStyle','none',...
         'FontSize', 14);
     stopBtn = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Stop');
     stopBtn.Callback = @plotButtonPushed;
     stopBtnFlag = stopBtn.Value;
 else
     stopBtnFlag = 0;
 end
 
% Set maximum number of iterations
if ~exist('iterMax', 'var') || isempty(iterMax)
    iterMax = 1000;
end

% Set convergence limit
if ~exist('stopTol', 'var') || isempty(stopTol)
    stopTol = 0.05;
end

if ~exist('gpu_yes','var') || isempty(gpu_yes)
    gpu_yes = false;
end

[Nf,K] = size(I0);

Obj = guess;
Obj = Obj/max(abs(Obj)); % normalization
Ir = calc_FROGtrace(Obj, delay, freq); % its FROG trace

%gamma = 5e-6; % experimental soft-thresholding parameter
FROGtrace_error = zeros(iterMax+1, 1);
iter = 1;
FROGtrace_error(1) = sqrt( sum(sum((I0-Ir).^2)) );
bestObj = Obj;
Inoise = zeros(size(I0));%/15000;

if gpu_yes
    I0 = gpuArray(I0);
    I = gpuArray(I);
    delay = gpuArray(delay);
    freq = gpuArray(freq);
    Obj = gpuArray(Obj);
    Inoise = gpuArray(Inoise);
    Ir = gpuArray(Ir);
    FROGtrace_error = gpuArray(FROGtrace_error);
end
% Since we have measured the real spectrum, it's also taken into account 
% during the retrieval process.
if use_measured_spectrum
    guess = guess/max(abs(guess)); guess(abs(guess)<0.01) = 0;
    measured_spectrum = abs(fft(guess));
end

while iter <= iterMax
    s = randperm(K);
    
    %scal = 2+2*tanh((iter-iterMax/2)/20);
    scal = 9*exp(-4*log(2)*((iter - iterMax)/(iterMax*0.7))^2); % 4*log(2) = 2.77
    alpha = abs( 0.8+randn(1,1)/100 )/(scal+1);
    if mod(iter,10) == 0 && ~use_measured_spectrum
        % I = ifftshift(imgaussfilt(fftshift(I0,1),50/iter),1);
        I([ceil(size(I0,1)/2)+1:end,1:ceil(size(I0,1)/2)],:) = imgaussfilt(I0([ceil(size(I0,1)/2)+1:end,1:ceil(size(I0,1)/2)],:),50/iter);
    end
    
    % Take the measured spectrum into account during the retrieval iteration here.
    % It's considered every 10 iterations; however, its ratio becomes
    % smaller and smaller with increasing iterations.
    if mod(iter,10) == 0 && use_measured_spectrum
        ObjS = fft(Obj);
        Obj_measured = ifft(measured_spectrum*max(abs(ObjS))/max(measured_spectrum).*exp(1i*angle(ObjS)));
        Obj_measured(abs(Obj_measured)<max(abs(Obj_measured))/50) = 0;
        
        Obj = Obj*(1-ratio_of_measured_spectrum*10/iter) + Obj_measured*ratio_of_measured_spectrum*10/iter;
    end
    
    % Smooth the field:
    % The pulse is almost always smooth or well-defined when we use FROG
    % for pulse retrieval, so this artifitial smoothing is applied for a
    % faster convergence.
    %if exist('smooth','file') % "smooth" requires MATLAB "curve fitting toolbox"
    %    Obj = smooth(abs(Obj),ceil(length(Obj)/800)).*exp(1i*angle(Obj));
    %end
    
    % Update the objective field on the order of a randomly permuted delay
    for iterK =1:K
        delayed_Obj = ifft( fft(Obj).*exp(1i*2*pi*delay(s(iterK)).*freq) );
        SHG = Obj.*delayed_Obj;
        
        SHG_FROGtrace = fft(SHG)/Nf;
        SHG_FROGtrace = SHG_FROGtrace./abs(SHG_FROGtrace).*( I(:, s(iterK))-Inoise(:, s(iterK)) ); % replace the modulus of the reconstructed FROG trace by the measured one
        SHG_FROGtrace(isnan(SHG_FROGtrace)) = 0;
        % Experimental soft thresholding:
        % It's slow and seems not necessary for retrieval, so I took it off.
        %SHG_FROGtrace = (real(SHG_FROGtrace) - gamma*sign(real(SHG_FROGtrace))).*(abs(SHG_FROGtrace) >= gamma)+...
        %         1i*(imag(SHG_FROGtrace) - gamma*sign(imag(SHG_FROGtrace))).*(abs(SHG_FROGtrace) >= gamma);
        SHG_n = ifft(SHG_FROGtrace)*Nf;
        
        Inoise(:, s(iterK)) = Inoise(:, s(iterK)) + 1e-6*( SHG_FROGtrace - I(:, s(iterK)) );

        % Update the objective function based on the following weight function
        weight = conj(delayed_Obj)./max( (abs(delayed_Obj).^2) );
        Obj = Obj + alpha.*weight.*(SHG_n - SHG);
    end

    Ir = calc_FROGtrace(Obj, delay, freq);
    FROGtrace_error(iter+1) = sqrt( sum(sum(abs(I-Ir).^2)) );

    % Display the result
    %fprintf('Iter:%d   alpha=%d Error=%d\n',iter, alpha, FROGtrace_error(iter+1));
    if DispFlag
        if ~isvalid(fig)
            break
        end
        
        % Display the iteration and the error on the figure
        % When the iteration is larger than 50, the Gaussian smoothing is
        % turned off. The best retrieved pulse with the lowest error will 
        % be saved afterward.
        if iter >= 50
            minErroridx = 51;
        else
            minErroridx = 1;
        end
        ann.String = {sprintf('Iter:%d   alpha=%5.3f NMSE=%5.4f(%5.4f)\n',...
            iter, alpha, FROGtrace_error(iter+1), min(FROGtrace_error(minErroridx:iter+1)));};
        
        % I_Obj = abs(fftshift(Obj)).^2;
        I_Obj = abs(Obj([ceil(length(Obj)/2)+1:end,1:ceil(length(Obj)/2)])).^2;
        difference = I_Obj-max(I_Obj)/2;
        left = find(difference>0,1);
        right = find(difference>0,1,'last');
        if isempty(left) || isempty(right)
            pulse_FWHM = 0;
        else
            pulse_FWHM = (right-left)*dt;
        end
        
        
        if ~exist('ax1','var')
            subplot(2,2,1);
            ax1 = gca;
        end
        yyaxis(ax1,'left');
        % plot(ax11,time*1e12, I_Obj, 'LineWidth',2);
        plot(ax1,time*1e12, I_Obj, 'LineWidth',2);
        xlabel(ax1,'Time (ps)','FontSize',16);
        ylabel(ax1,'Power (a.u.)','FontSize',16);
        title(ax1,'Pulse');
        
        yyaxis(ax1,'right');
        time_angle = unwrap(angle(Obj([ceil(length(Obj)/2)+1:end,1:ceil(length(Obj)/2)])))/pi;
        time_angle(I_Obj<max(I_Obj)/200) = NaN;
        % plot(ax1,time*1e12, unwrap(angle(fftshift(Obj,1)))/pi, 'LineWidth',2);
        plot(ax1,time*1e12, time_angle, 'LineWidth',2);
        ylabel(ax1,'Phase (\pi)','FontSize',16);
        xlim(ax1,[min(time*1e12) max(time*1e12)]);
        title(ax1,sprintf('Approx. FWHM = %6.5f(ps)',pulse_FWHM));

        ObjS = fft(Obj);
        % ObjS = fftshift(ObjS,1)/max(ObjS);
        ObjS = ObjS([ceil(length(ObjS)/2)+1:end,1:ceil(length(ObjS)/2)])/max(ObjS);
        I_ObjS = abs(ObjS).^2;
        
        if ~exist('ax21','var')
            subplot(2,2,2)
            ax2 = gca;
        end
        yyaxis(ax2,'left');
        plot(ax2,freq_plot, I_ObjS, 'LineWidth',2);
        xlabel(ax2,'Frequency (THz)','FontSize',16);
        ylabel(ax2,'PSD (a.u.)','FontSize',16);
        title(ax2,'Pulse');
        
        yyaxis(ax2,'right');
        freq_angle = unwrap(angle(ObjS))/pi;
        freq_angle(I_ObjS<max(I_ObjS)/200) = NaN;
        plot(ax2,freq_plot, freq_angle, 'LineWidth',2);
        ylabel(ax2,'Phase (\pi)','FontSize',16);
        xlim(ax2,[min(freq_plot),max(freq_plot)]);

        if ~exist('ax3','var')
            subplot(2,2,3)
            ax3 = gca;
        end
        % imagesc(ax3,delay*1e12, freq_plot, abs(fftshift(I,1)) );
        imagesc(ax3,delay*1e12, freq_plot, abs(I([ceil(length(I)/2)+1:end,1:ceil(length(I)/2)],:)) );
        title(ax3,'Measured FT');
        xlabel(ax3,'Time (ps)','FontSize',16); ylabel(ax3,'Freq.(THz)','FontSize',16);
        colormap(ax3,'jet');

        if ~exist('ax4','var')
            subplot(2,2,4)
            ax4 = gca;
        end
        % imagesc(ax4,delay*1e12, freq_plot, abs(fftshift(Ir,1)) );
        imagesc(ax4,delay*1e12, freq_plot, abs(Ir([ceil(length(Ir)/2)+1:end,1:ceil(length(Ir)/2)],:)) );
        title(ax4,'Recovered FT');
        xlabel(ax4,'Time (ps)','FontSize',16); ylabel(ax4,'Freq.(THz)','FontSize',16);
        colormap(ax4,'jet');
        drawnow;
    else
        I_Obj = abs(Obj([ceil(length(Obj)/2)+1:end,1:ceil(length(Obj)/2)])).^2;
    end
    
    % Shift the pulse to the time center t=0
    [Obj,max_I] = pulse_centering(I_Obj,Obj,Nf);
    
    % After 50 iterations, the Gaussian smoothing of the original FROGtrace
    % will be turned off.
    if iter >= 50
        % Stop and show error if the time window is considered too small
        if I_Obj(1) > max_I/100 || I_ObjS(1) > max(I_ObjS)/100
            message = {'There''s aliasing! Please:','1. increase the time window, and resample again.','2. Increase the guess bandwidth that affects resampling, and resample again','3. reduce the guess pulse duration'};
            uialert(UIFigure,message,'Error','Icon','error');
            return
        end
        %{
        if error(iter+1) <= min(error(51:iter+1))
            bestObj = Obj;
        else
            Obj = (bestObj*error(iter+1)+Obj*min(error(51:iter+1)))/(min(error(51:iter+1))+error(iter+1));
        end
        %}

        % Stop if it's done and the field is successfully retrieved
        if FROGtrace_error(iter+1) < stopTol || stopBtnFlag
            break;
        end
    end
    
    iter = iter + 1;
    
    % Fix/Calibrate the delay vector based on the reconstructed field
    if mod(iter,25) == 0
        % [output, ~] = dftregistration(fft2(fftshift(Ir,1)),fft2(fftshift(I,1)),50);
        [output, ~] = dftregistration(fft2(Ir([ceil(size(Ir,1)/2)+1:end,1:ceil(size(Ir,1)/2)],:)),fft2(I([ceil(size(I,1)/2)+1:end,1:ceil(size(I,1)/2)],:)),50);
        delay = delay + mean(diff(delay))*output(end);
    end

end
bestObj = Obj;

% Remove background noise in the retrieved field
% This step is important for later computation of the transform-limited
% field to remove the sharp tiny temporal spike from the broadband
% spectral noise.
ObjS = fft(bestObj);
ObjS(abs(ObjS).^2<max(abs(ObjS).^2)/200) = 0;
bestObj = ifft(ObjS);

if gpu_yes
    bestObj = gather(bestObj);
    FROGtrace_error = gather(FROGtrace_error);
    Ir = gather(Ir);
end
bestObj = fftshift(bestObj,1);
Ir = fftshift(Ir,1);

% Stop button
% (It needs to be inside the main function.)
function plotButtonPushed(~,~) % (src,event); src and event are two necessary inputs for a callback
    stopBtnFlag = stopBtn.Value;
end

end

%% Helper functions
function I = calc_FROGtrace(pulse, delay, freq)
%CALC_FROGTRACE It calculates the FROG trace of a pulse.
% A FROG trace is the Fourier Transform of the product of a pulse and its
% delayed duplicate.

delayed_pulse = ifft( fft(pulse).*exp(1i*2*pi*delay.*freq) );
I = abs(fft(pulse.*delayed_pulse)/length(pulse));
%I = I/sqrt( sum(sum(I.^2)) ); % normalization

end

function [Obj,max_I] = pulse_centering(I_Obj,Obj,Nf)
%PULSE_CENTERING Shift the pulse to the time center t=0

[max_I,pulse_center] = max(I_Obj);
index_shift = mod(pulse_center-floor(Nf/2)-1,Nf);
Obj = Obj([index_shift+1:end,1:index_shift]);

end