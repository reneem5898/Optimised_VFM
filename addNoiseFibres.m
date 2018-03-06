function fibresNoise = addNoiseFibres(orDir, orientation, perc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function adds Gaussian noise to a set of fibre orientations and writes
% the orientations (with noise) to a new text file
%
% Inputs: 1) orDir - directory where to save resulting text files
%         2) orientation - material orientations without Gaussian noise
%         3) perc - percentage to multiply by standard deviation of orientations
%
% Output: fibresNoise - orientations with Gaussian noise applied
%
% Noise calculated as:
% stdevNoise = p*stdevDisp;
% e.g. p = 0.15 (== 15%)
%
% Renee Miller
% Updated: 6 March 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate standard deviation of absolute value of displacements
stdFibres = std(orientation);

% Standard deviation of noise = perc * stdFibres
stdNoise = perc*stdFibres;

% Initialise variable of fibres + noise
fibresNoise = orientation(:,1);

% Fibre noise directory - where to save orientations created
fibreOutDir = sprintf('%s/fibresNoise/', orDir);
if ~exist(fibreOutDir,'dir')
    mkdir(fibreOutDir);
end

disp('Adding Gaussian noise to fibre orientations...')
for d = 2:length(stdNoise) % Start at 2 since first column is element number
       
    % Create noise - random distribution centered at zero
    Noise = normrnd(0, stdNoise(d), size(orientation,1), 1);
    
    % Add noise to signal
    fibresNoise = [fibresNoise, orientation(:,d) + Noise];
    
    % Write noisy displacements to a file
    writeFile = sprintf('%s/orientations-noise.txt', fibreOutDir);
    dlmwrite(writeFile, fibresNoise, 'delimiter', '\t', 'precision', '%1.8g');
    
end