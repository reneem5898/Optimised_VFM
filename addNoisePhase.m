function addNoisePhase(outDir, noiseDir, numNoise, percentages)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function adds Gaussian noise to a set of displacements and writes 
% the displacement fields (with noise) to a new text file
%
% Inputs: 1) modelDir - directory where to save resulting text files
%         2) noiseDir - where to save new noisey displacement files
%         3) numNoise - number of Monte-Carlo noise simulations
%         4) percentatges - list of percentages, varying amounts of noise
%   
% Output: text files - uComplex_Noise%d.txt (%d = counter)
%
% Noise calculated as:
% stdevNoise = p*stdevDisp;
% e.g. p = 0.15 (== 15%)
%
% Renee Miller
% Updated: 21 March 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open reference displacement field
dispFile = sprintf('%s/uComplex.txt', outDir);
U = load(dispFile);
U = U(:,1) + U(:,2)*1i;

% Interpolate complex displacements (U) at four phase offsets
noff = 4; % 4 phase offsets to interpolate displacements at
phase = zeros(length(U),noff);
for k = 1:noff
    phase(:,k) = abs(U).*cos(2*pi*(k-1)/noff+angle(U));
end

% phasePlot = zeros(length(U),200);
% for k = 1:200
%     phasePlot(:,k) = abs(U).*cos(2*pi*(k-1)/noff+angle(U));
% end

% Calculate standard deviation of absolute value of displacements
stdDisp = std(phase(:));

% Counter - for naming files
x = 0;

% Create waitbar
WH = waitbar(0, 'Creating noise copies...');

% For each value of p, generate standard deviation of Gaussian noise
for p = 1:length(percentages)

	% Standard deviation of noise to add to displacements
    stdNoise = percentages(p)*stdDisp;
    
    for i = 1:numNoise
        
        % Counter
        x = x + 1;
        
        % Update waitbar to give user an indication of time
        percentComplete = x/(length(percentages)*numNoise);
        waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))        
        
        % Create noise - random distribution centered at zero
        phaseNoise = normrnd(0, stdNoise, size(phase,1), size(phase,2));
        % Add noise to signal
        phaseNewNoise = phase + phaseNoise;
        
        % Take FFT of noisy phase offset data
        fourier = fft(phaseNewNoise,[],2);
        
        % 1st location of fft corresponds to 0th harmonic or DC, 2nd location in
        % fft corresponds to 1st harmonic, and so on for n harmonics. We need the first harmonic.
        fha = squeeze(fourier(:,2));
        fha = fha*2/noff;
        
        % Put into two column array
        uNoisey = [real(fha) imag(fha)];
        
        % Write noisy displacements to a file
        writeFile = sprintf('%s/uComplex_Noise%d.txt', noiseDir, x);
        dlmwrite(writeFile, uNoisey, 'delimiter', '\t', 'precision', '%1.8e');
        
    end
    
end

close(WH);