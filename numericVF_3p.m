function [uVF1, uVF2, uVF3, eta] = numericVF_3p(U, nodesSubZone, elemSubZone, globalNodeNums, paramInit, orientation, surfaceNodes, elemType, GaussPoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates numeric virtual displacement field for 
% estimating three parameters which describe a transversely isotropic
% material
%
% The following conditions are specified:
% 1. C1: c_k 
% 2. C2: c_mu12
% 3: C3: c_mu13
% 4: C4: c_T
% 5. sum(uVFx) = 0, sum(uVFy) = 0, sum(uVFz) = 0 (Arb)
% 6. uVF(boundaries) = 0 (Acl)
% 7. Noise minimisation (H)
%
% Inputs: 1) u - MRE displacement field
%         2) nodesSubZone - list of node numbers and nodal coordinates (numNodes x 3)
%         3) elemsSubZone - (numElems x 8) - 8 node numbers that make up the
%         element
%         4) globalNodeNums - all nodes numbers in entire model
%         5) paramInit - initial parameter estimates (used in H)
%         6) orientation - material orientation
%         7) surfaceNodes - list of surface nodes (not required)
%         8) elemType - string denoting integration type to use: 'C3D8R', 'C3D8', 'C3D8F'
%         9) GaussPoints - number of integration points (per direction) to use
%            - input to 'C3D8F' integration type
%
% Outputs: 1) uVF1, uVF2, uVF3 - special and optimised virtual displacement fields
%          2) eta - sensitivity values
%
%
% Written by: Renee Miller
% Date: 28 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Construct Acl matrix - boundary nodes uVF = 0
% Find boundary nodes and set displacement condition to zero

disp('Constructing boundary node constraint matrix...');

% Get boundary nodes, unless given
if isempty(surfaceNodes)
    boundaryNodes = getBoundaryNodes(nodesSubZone, elemSubZone);
else
    boundaryNodes = surfaceNodes;
end

% Create boundary constraint matrices
[Acl, rhs_Acl] = createBoundaryConstraint(boundaryNodes, nodesSubZone);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: Construct Ak, Ag and H matrices

disp('Constructing Ak, Ag and Hg constraint matrices...');

% Uniform strain elements - C3D8R
if strcmp(elemType, 'C3D8R')
    [C1, C2, C3, C4, H] = C3D8R_3p(U, nodesSubZone, elemSubZone, globalNodeNums, paramInit, orientation);

% Selectively reduced integration type element - C3D8
elseif strcmp(elemType, 'C3D8')
    [C1, C2, C3, C4, H] = C3D8_3p(U, nodesSubZone, elemSubZone, globalNodeNums, paramInit, orientation);

% Fully integrated element 
else
    [C1, C2, C3, C4, H] = C3D8F_3p(U, nodesSubZone, elemSubZone, globalNodeNums, GaussPoints, paramInit, orientation);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3: Construct Arb - condition where the sum of all displacements in each
% direction: x, y and z, respectively are 0

disp('Constructing the constraint on the sum of virtual displacements in x, y and z...');

% Nodal DOFs
DOF = size(nodesSubZone,2) - 1;

% Construct repeated identity matrix: I I I I ... N times (N = number of
% nodes)
Arb = zeros(DOF,length(nodesSubZone)*DOF);
for i = 1:size(Arb,1)
    Arb(i,i:DOF:end) = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4: Compile constraint matrices and solve for virtual displacement fields

disp('Solving for the virtual displacement field...')

% Compile constraints
C = [C1; C2; C3; C4];
A = [Acl; Arb; C];
A = sparse(A);

% Assemble LHS
z1 = zeros(size(A,1),size(A,1));
% Need .' rather than ' since A is complex - want just transpose, not conjugate transpose
LHS_1 = [H A.'; A z1];
LHS_2 = [H A.'; A z1];
LHS_3 = [H A.'; A z1];

%   boundary   **** bulk wave noise cons *****     C1       C2        C3         C4       
Zg1 = [rhs_Acl; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 1.0 + 0*1i; 0 + 0*1i; 0 + 0*1i];
%   boundary   **** bulk wave noise cons *****     C1       C2        C3         C4       
Zg2 = [rhs_Acl; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 1.0 + 0*1i; 0 + 0*1i];
%   boundary   **** bulk wave noise cons *****     C1       C2        C3         C4       
Zg3 = [rhs_Acl; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 1.0 + 0*1i];

% Assemble RHS
z2 = complex(zeros(size(H,1),1)); % Minimise H
RHS_1 = [z2; Zg1];
RHS_2 = [z2; Zg2];
RHS_3 = [z2; Zg3];

% Solve for virtual displacements and lagrangian multipliers which satisfy conditions
x_1 = LHS_1\RHS_1;
x_2 = LHS_2\RHS_2;
x_3 = LHS_3\RHS_3;

% Return virtual displacement fields
uVF1 = x_1(1:size(H,1));
uVF2 = x_2(1:size(H,1));
uVF3 = x_3(1:size(H,1));

% Calculate eta values (sensitivities): conjugate transpose = uVF'
eta1 = uVF1' * H * uVF1;
eta2 = uVF2' * H * uVF2;
eta3 = uVF3' * H * uVF3;
eta = [eta1; eta2; eta3];

% The end