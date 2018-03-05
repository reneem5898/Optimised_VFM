function fibresNoise = addNoiseFibres(outDir, orientation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function adds Gaussian noise to a set of displacements and writes
% the displacement fields (with noise) to a new text file
%
% Inputs: 1) modelDir - directory where to save resulting text files
%         2) noiseDir - where to save new noisey displacement files
%         3) numNoise - number of Monte-Carlo noise simulations
%         4) percentatges - list of percentages, varying amounts of noise
%
% Output: text files - uComplex_noise_%d.txt (%d = counter)
%
% Noise calculated as:
% stdevNoise = p*stdevDisp;
% e.g. p = 0.15 (== 15%)
%
% Renee Miller
% Updated: 21 March 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate standard deviation of absolute value of displacements
stdFibres = std(orientation);

% Counter - for naming files
x = 0;

% Create waitbar
WH = waitbar(0, 'Creating noise copies...');

% Percentage of noise in fibre angles
perc = 0.12; % 12% noise -- results in a variation in fibre angles with 2*std = 5degree

% Standard deviation of noise = perc * stdFibres
stdNoise = perc*stdFibres;

% Initialise variable of fibres + noise
fibresNoise = orientation(:,1);

for d = 2:length(stdNoise) % Start at 2 since first column is element number
    
    % Counter
    x = x + 1;
    
    % Update waitbar to give user an indication of time
    percentComplete = x/length(stdNoise);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
    % Create noise - random distribution centered at zero
    Noise = normrnd(0, stdNoise(d), size(orientation,1), 1);
    
    % Add noise to signal
    fibresNoise = [fibresNoise, orientation(:,d) + Noise];
    
    % Write noisy displacements to a file
    writeFile = sprintf('%s/orientations-noise.txt', outDir);
    dlmwrite(writeFile, fibresNoise, 'delimiter', '\t', 'precision', '%1.8g');
    
end

close(WH);