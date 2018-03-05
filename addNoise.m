function addNoise(outDir, noiseDir, numNoise, percentages)

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

% Open reference displacement field
dispFile = strcat(outDir, '/uComplex.txt');
U = load(dispFile);
%U = U(:,1) + U(:,2)*1i;

% Calculate standard deviation of absolute value of displacements
stdDisp = std(U);

% Counter - for naming files
x = 0;

% Create waitbar
WH = waitbar(0, 'Creating noise copies...');

% For each value of p, generate standard deviation of Gaussian noise
for p = 1:length(percentages)

    stdNoise = percentages(p)*stdDisp;
    
    for i = 1:numNoise
        
        % Counter
        x = x + 1;
        
        % Update waitbar to give user an indication of time
        percentComplete = x/(length(percentages)*numNoise);
        waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))        
        
        % Create noise - random distribution centered at zero
        realNoise = normrnd(0, stdNoise(1), size(U,1), 1);
        imagNoise = normrnd(0, stdNoise(2), size(U,1), 1);
         
        % Add noise to signal
        dispNewNoise = [U(:,1) + realNoise, U(:,2) + imagNoise];
               
        % Write noisy displacements to a file
        writeFile = sprintf('%s/uComplex_noise_%d.txt', noiseDir, x);
        dlmwrite(writeFile, dispNewNoise, 'delimiter', '\t', 'precision', '%1.8e');
        
    end
    
end

close(WH);