function [uVF, eta] = numericVF_Visco(uR, uI, nodesSubZone, elemSubZone, globalNodeNums, elemType, GaussPoints, surfaceNodes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates numeric virtual displacement field following the
% methods outlined in Connesson et al. for a harmonic displacement field
% INCLUDING VISCOELASTICITY
%
% The following conditions are specified:
% 1. fk = 0 (AkR and AkI)
% 2. fg = 1 (AgR and AgI)
% 3. sum(uVFx) = 0, sum(uVFy) = 0, sum(uVFz) = 0 (Arb)
% 4. uVF(boundaries) = 0 (Acl)
% 5. Noise minimisation (H)
%
% Inputs: 1) uR - real component of MRE displacement field
%         2) uI - imaginary component of MRE displacement field
%         3) nodes - list of node numbers and nodal coordinates (numNodes x 3)
%         4) elems - (numElems x 8) - 8 node numbers that make up the
%         element
%         5) DOF - number of degrees of freedom for each node
%
% Outputs: 1) uVF - special and optimised virtual displacement field
%       
%
% Written by: Renee Miller
% Date: 8 August 2016
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
    [AkR, AgR, AkI, AgI, H] = isoC3D8R_visco(uR, uI, nodesSubZone, elemSubZone, globalNodeNums);

% Selectively reduced integration type element - C3D8
elseif strcmp(elemType, 'C3D8')
    [AkR, AgR, AkI, AgI, H] = isoC3D8_visco(uR, uI, nodesSubZone, elemSubZone, globalNodeNums);

% Fully integrated element 
else
    [AkR, AgR, AkI, AgI, H] = isoC3D8F_visco(uR, uI, nodesSubZone, elemSubZone, globalNodeNums, GaussPoints);

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

%% Step 4: Compile constraint matrices and solve for virtual displacement field

disp('Solving for the virtual displacement field...')

A = [Acl; Arb; AkR; AkI; AgR; AgI];
A = sparse(A);

% Compile LHS
z1 = zeros(size(A,1),size(A,1));
LHS = [H A'; A z1];

% Compile RHS
%     Acl       %%%%%%%%% Arb %%%%%%%%%%%      AkR         AkI       AgR            AgI
Zg = [rhs_Acl; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 1.0 + 0.0*1i; 1.0 + 0.0*1i];
z2 = complex(zeros(size(H,1),1)); % Minimise H
RHS = [z2; Zg];

% Solve for virtual displacements and lagrangian multipliers which satisfy conditions
x = LHS\RHS;

% Return uVF (virtual displacement field)
uVF = x(1:size(H,1));

% Return eta (sensitivity)
eta = uVF.' * H * uVF;