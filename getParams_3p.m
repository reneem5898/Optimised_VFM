function [G12, G13, damp, T] = getParams_3p(moduli, iter, maxIter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the elastic parameters - Young's moduli, shear moduli,
% Poisson's ratios and damping coefficients from the complex valued parameters
% of the elasticity matrix
%
% Inputs: 1) moduli - vector of complex moduli (C11, C33, C44, C66, C13)
%         2) iter - number of iterations it took to convergence
%         3) maxIter - maximum number of allowable iterations
%
% Outputs: parameters - G12, G13, T, damp
%
% Written by: Renee Miller
% Updated: 7 November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If analysis for subzone reached maximum number of iterations, save NaN's
if iter == maxIter
    G12 = NaN; G13 = NaN; T = NaN;
    damp = [NaN; NaN; NaN];
else
    
    % Parameters
    G12 = real(moduli(1));
    G13 = real(moduli(2));
    T = real(moduli(3));
    
    damp = imag(moduli)./real(moduli);
   
end