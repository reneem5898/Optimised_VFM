function [E, G12, G13, E1, E3, v12, v13, v31, damp] = getParams_5p(moduli, iter, maxIter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the elastic parameters - Young's moduli, shear moduli, 
% Poisson's ratios and damping coefficients from the complex valued parameters
% of the elasticity matrix
%
% Inputs: 1) moduli - vector of complex moduli (C11, C33, C44, C66, C13)
%         2) iter - number of iterations it took to convergence
%         3) maxIter - number of maximum iterations
% 
% Outputs: parameters - E (elasticity matrix), G12, G13, E1, E3, v12, v13, v31, damp
%
% Written by: Renee Miller
% Updated: 7 November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If analysis for subzone reached maximum number of iterations, save NaNs's
if iter == maxIter
    G12 = NaN; G13 = NaN; E1 = NaN; E3 = NaN;
    damp = [NaN; NaN; NaN; NaN];
    E = NaN(6,6);
    v12 = NaN; v13 = NaN; v31 = NaN;

else
    
    % Parameters of elasticity matrix
    C11 = moduli(1);
    C33 = moduli(2);
    C44 = moduli(3);
    C66 = moduli(4);
    C13 = moduli(5);
    
    % C12 is related to C11 and C44
    C12 = C11 - 2*C44;
     
    % Full complex elasticity matrix
    E = [C11    C12     C13     0   0   0; ...
         C12    C11     C13     0   0   0; ...
         C13    C13     C33     0   0   0; ...
         0      0       0       C44 0   0; ...
         0      0       0       0   C66 0; ...
         0      0       0       0   0   C66];
    
    % Compliance matrix = inverse of elasticity matrix
    C = inv(E);
    
    % Parameters E1, E3, G12, G13, and v12 are obtained from the compliance matrix
    E1 = 1/C(1,1);
    E3 = 1/C(3,3);
    G12 = 1/C(4,4);
    G13 = 1/C(6,6);
    v12 = -1*C(1,2)*E1;
    v31 = -1*C(1,3)*E3;
    v13 = -1*C(3,1)*E1;
    
    % Damping coefficients
    damp = [imag(E1)/real(E1); imag(E3)/real(E3); imag(G12)/real(G12); imag(G13)/real(G13)];
end
