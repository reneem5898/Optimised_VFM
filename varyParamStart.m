function paramInit = varyParamStart(refParams)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function varies the starting parameter used as an initial guess in 
% the anisotropic VFM
%
% Input: reference parameters
% Output: initial parameter estimtes
%
% Renee Miller
% Date: 1 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise list of parameters
paramInit = length(refParams);

for i = 1:length(refParams) % Loop through parameters

    % Compute standard deviation = 20% of true value
    std = 0.20 * refParams(i); 
    
    % Get random starting value around true value +/- std
    paramInit(i) = normrnd(refParams(i), abs(std), 1); % Take care if stdev is negative: abs(std)
    
end